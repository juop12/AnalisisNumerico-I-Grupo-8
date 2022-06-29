from numpy import zeros


# Funciones auxiliares para resolver un sistema de ecuaciones
# lineal por el metodo de factorizacion LU DooLittle
#=====================================================================

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

# Descompone una matriz en dos matrices (L y U) con la 
# forma del metodo DooLittle.
# Devuelve las matrices L y U.
def factorizacionLU_DooLittle(mat):
    cantidadDeCoeficientes = len(mat)

    lower = zeros([cantidadDeCoeficientes, cantidadDeCoeficientes])
    upper = zeros([cantidadDeCoeficientes, cantidadDeCoeficientes])
    
    for i in range(cantidadDeCoeficientes):
        for j in range(i, cantidadDeCoeficientes):

            # Calculo valores de U.
            upper[i][j] = mat[i][j]
            for k in range(i):
                upper[i][j] -= lower[i][k] * upper[k][j]

            # Calculo valores de L.
            lower[j][i] = mat[j][i]
            for k in range(i):
                lower[j][i] -= lower[j][k] * upper[k][i]
            lower[j][i] /= upper[i][i]
    
    return lower, upper


# Main Protocol
#======================================================================

# Resuelve un sistema LU.
def resolverSistemaLU(lower, upper, B):
    y = calcular_y(lower, B)
    return calcular_x(upper, y)

# Ax=B, separa en LUx=B:
#   - Ly=B
#   - Ux=y
#   resuelve los sistemas y devuelve x.
def resolverSistemaLineal(A, B):
    lower, upper = factorizacionLU_DooLittle(A)
    return resolverSistemaLU(lower, upper, B)