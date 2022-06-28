from numpy import zeros

# Resuelve el sistema Ly=B para matriz L del metodo DooLittle.
# Devuelve el vector y.
def calcular_y(lower, B):
    cantidadDeCoeficientes = len(lower)

    y = zeros(cantidadDeCoeficientes)
    for i in range(cantidadDeCoeficientes):
        y[i] = B[i]
        for j in range(i):
            y[i] -= lower[i][j] * y[j]

    return y

# Resuelve el sistema Ux=y para matriz U del metodo DooLittle.
# Devuelve el vector x.
def calcular_x(upper, y):
    cantidadDeCoeficientes = len(upper)

    x = zeros(cantidadDeCoeficientes)
    for i in range(cantidadDeCoeficientes-1, -1, -1):
        x[i] = y[i]
        for j in range(cantidadDeCoeficientes-1, i, -1):
            x[i] -= upper[i][j] * x[j]

        x[i] /= upper[i][i]

    return x

# Descompone una matriz en dos matrices (LU) con la 
# forma del metodo DooLittle.
# Devuelve las matrices L y U.
# Fuente: https://www.youtube.com/watch?v=FpVeXhAQg9w.
def factorizacionLU_DooLittle(mat):
    cantidadDeCoeficientes = len(mat)

    lower = zeros([cantidadDeCoeficientes, cantidadDeCoeficientes])
    upper = zeros([cantidadDeCoeficientes, cantidadDeCoeficientes])
    
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
#======================================================================

# Resuelve un sistema LU.
def resolverSistemaLU(lower, upper, B):
    y = calcular_y(lower, B)
    return calcular_x(upper, y)

# Ax=B, separo en LUx=B:
#   - Ly=B
#   - Ux=y
#   resuelvo los sistemas y devuelvo x.
def resolverSistemaLineal(A,B):
    lower, upper = factorizacionLU_DooLittle(A)
    return resolverSistemaLU(lower, upper, B)
    








# Pruebas
#======================================================================
import numpy as np


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

    y = np.linalg.solve(lower, b)
    x = np.linalg.solve(upper, y)
    print(x)



#pruebas()