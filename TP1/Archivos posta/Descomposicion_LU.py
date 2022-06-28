import numpy as np

# Resolución de un sistema de ecuaciones lineal con matriz LU.

#https://www.youtube.com/watch?v=34k1y31EZiY
def calcular_y(lower, b):
    cantidadDeCoeficientes = len(lower)

    y = np.zeros(cantidadDeCoeficientes)
    for i in range(cantidadDeCoeficientes):
        y[i] = b[i]
        for j in range(i):
            y[i] -= lower[i][j] * y[j]

    return y

def calcular_x(upper, y):
    cantidadDeCoeficientes = len(upper)

    x = np.zeros(cantidadDeCoeficientes)
    for i in range(cantidadDeCoeficientes-1, -1, -1):
        x[i] = y[i]
        for j in range(cantidadDeCoeficientes-1, i, -1):
            x[i] -= upper[i][j] * x[j]

        x[i] /= upper[i][i]

    return x

# Descomposicion LU

#https://www.youtube.com/watch?v=FpVeXhAQg9w
def factorizacionLU_DooLittle(mat):
    cantidadDeCoeficientes = len(mat)

    lower = np.zeros([cantidadDeCoeficientes, cantidadDeCoeficientes])
    upper = np.zeros([cantidadDeCoeficientes, cantidadDeCoeficientes])
    
    for i in range(0, cantidadDeCoeficientes):
        upper[0][i] = mat[0][i]
        for j in range(0, cantidadDeCoeficientes):
            if (i==j):
                lower[i][j] = 1
            if (i<j):
                factor = (mat[j][i]/mat[i][i])
                lower[j][i] = factor
                for k in range(0, cantidadDeCoeficientes):
                    mat[j][k] = mat[j][k] - factor * mat[i][k]
                    upper[j][k] = mat[j][k]

    return lower, upper

# Main Protocol

# Ax=B, separo en:
#   - Ly=x
#   - Ux=B
#   devuelvo x.
def solve(A,B):
    lower, upper = factorizacionLU_DooLittle(A)
    y = calcular_y(lower, B)
    return calcular_x(upper, y)

# Pruebas
#======================================================================

# Imprimir valores

def imprimirResolucionDeSistema(x):
    for i in range(len(x)):
        print("x" + str(i + 1) + " = ", x[i])
    print("x" " = ", x)

def imprimirDescomposicionLU(l, u):
    print("Matriz L")
    print(l)
    print("\n")

    print("Matriz U")
    print(u)
    print("\n")

    A = np.dot(l, u, out=None)
    print("Matriz Original")
    print(A)
    print("\n")

# testing

def pruebas():
    
    #mat = np.array([[7, 10, 4],
    #                [5, -2, 6],
    #                [3, 1, -1]])
    #b = np.array([-2, 38, 21])

    mat = np.array([[2, 4, 2, 6],
                    [4, 9, 6, 15],
                    [2, 6, 9, 18],
                    [6, 15, 18, 40]])
    b = np.array([9, 23, 22, 47])

    lower, upper = factorizacionLU_DooLittle(mat)
    imprimirDescomposicionLU(lower, upper)

    y = calcular_y(lower, b)
    x = calcular_x(upper, y)
    imprimirResolucionDeSistema(x)

    #y = np.linalg.solve(lower, b)
    #x = np.linalg.solve(upper, y)
    #print(x)



#pruebas()