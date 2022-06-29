from ResolucionDeSistemasLinealesLU_DooLittle import resolverSistemaLineal
import matplotlib.pyplot as plt
import tabulate as tab
import numpy as np


# Graficos
#======================================================================

# Imprime matrices A y B.
def imprimirMatricesAyB(A, B):
    print("Matriz A:")
    print(A, "\n")

    print("Matriz B:")
    print(str(B) + "\n")

# Imprime tabla de coeficientes de cada subfuncion de una Spline.
def mostrarTablaDeCoeficientes(a, b, c, d, cantidadDeNodos):
    datos = []
    titulos = ['S_j', 'a', 'b', 'c', 'd ']

    for i in range(cantidadDeNodos-1):
        datos.append(('S_'+str(i), a[i], b[i], c[i], d[i]))

    # Tabla ascii: tablefmt="github"
    print(tab.tabulate(datos, headers=titulos, floatfmt=".16f", tablefmt='fancy_grid'))
    print("\n", "="*100, "\n")

# Calculos para graficar Spline.

def obtenerValorSpline(a, b, c, d, x, x0):
    return a + b*(x-x0) + c*((x-x0)**2) + d*((x-x0)**3)

def obtenerPuntosCurvaSpline(a, b, c, d, x):
    resultado = []
    puntosEvaluados = []
    for i in range(len(x)-1):
        puntosAEvaluar = np.linspace(x[i], x[i+1], 20)
        resultadoActual = obtenerValorSpline(a[i], b[i], c[i], d[i], puntosAEvaluar, x[i])
        resultado = np.append(resultado, resultadoActual)
        puntosEvaluados = np.append(puntosEvaluados, puntosAEvaluar)

    return puntosEvaluados, resultado

def graficarSpline(a, b, c, d, x, y):
    puntosCurva, resultado = obtenerPuntosCurvaSpline(a, b, c, d, x)
    plt.plot(puntosCurva, resultado, 'k', lw=1)
    plt.scatter(x, y, color='red')


# Funciones auxiliares para calculo de Spline
#=====================================================================

# Calculo de matrices A y B.

def calcularMatrizA(h, cantidadDeNodos):
    A = np.zeros([cantidadDeNodos, cantidadDeNodos])

    # Primera fila
    A[0][0] = 2 * h[0]
    A[0][1] = h[0]

    # Tridiagonal
    for i in range(1, cantidadDeNodos - 1):
        A[i][i-1] = h[i-1]
        A[i][i] = 2*(h[i-1] + h[i])
        A[i][i+1] = h[i]

    # Ultima fila
    A[cantidadDeNodos-1][cantidadDeNodos-2] = h[cantidadDeNodos-2]
    A[cantidadDeNodos-1][cantidadDeNodos-1] = 2 * h[cantidadDeNodos-2]

    return A

def calcularMatrizB(a, h, cantidadDeNodos, pendienteInicial, pendienteFinal):
    # Primera posicion
    B = [3/h[0] * (a[1] - a[0]) - 3 * pendienteInicial]
    
    # Valores del centro
    for n in range(2, cantidadDeNodos):
        B.append(3/h[n-1] * (a[n] - a[n-1]) - 3/h[n-2] * (a[n-1] - a[n-2]))
    
    # Ultima posicion
    B.append(3 * pendienteFinal - 3/h[cantidadDeNodos-2] * (a[cantidadDeNodos-1] - a[cantidadDeNodos-2]))

    return B

# Calculo de coeficientes b y d.
def calcularCoeficientesRestantes(a, h, c, cantidadDeNodos):
    b = []
    for j in range(cantidadDeNodos - 1):
        b.append(1/h[j] * (a[j+1] - a[j]) - h[j]/3 * (2*c[j] + c[j+1]))

    d = []
    for j in range(cantidadDeNodos - 1):
        d.append(1/(3 * h[j]) * (c[j+1] - c[j]))

    return b, d


# Spline
#======================================================================

# Devuelve todos los coeficientes necesarios para armar una 
# spline ligada a partir de una tabla de valores.
def calcularCoeficientesDeSpline(x, a, pendienteInicial, pendienteFinal):
    cantidadDeNodos = len(x)
    h = [(x[n+1] - x[n]) for n in range(cantidadDeNodos - 1)]

    A = calcularMatrizA(h, cantidadDeNodos)
    B = calcularMatrizB(a, h, cantidadDeNodos, pendienteInicial, pendienteFinal)

    # Grafico de A y B.
    #imprimirMatricesAyB(A, B)

    # Resuelvo ecuacion
    c = resolverSistemaLineal(A, B)

    # Calculo coeficientes b, d
    b, d = calcularCoeficientesRestantes(a, h, c, cantidadDeNodos)

    return a, b, c, d

def spline(NombreDeCurva, x, y, pendienteInicial, pendienteFinal):
    print(NombreDeCurva)
    
    # Calculo los coeficientes de la Spline.
    a, b, c, d = calcularCoeficientesDeSpline(x, y, pendienteInicial, pendienteFinal)

    # Muestro la tabla de coeficientes de las funciones que conforman la Spline.
    mostrarTablaDeCoeficientes(a, b, c, d, len(x))
    
    # Creo el grafico de la Spline.
    graficarSpline(a, b, c, d, x, y)


# Main
#======================================================================

def main():

    # Tabla del enunciado hardcodeada.
    x1 = [  1,   2,   5,   6,   7,   8,  10,  13,  17]
    y1 = [3.0, 3.7, 3.9, 4.2, 5.7, 6.6, 7.1, 6.7, 4.5]

    x2 = [ 17,  20,  23,  24,  25,  27, 27.7]
    y2 = [4.5, 7.0, 6.1, 5.6, 5.8, 5.2,  4.1]

    x3 = [27.7,  28,  29,  30]
    y3 = [ 4.1, 4.3, 4.1, 3.0]

    spline("Curva 1", x1, y1, 1, -2/3)
    spline("Curva 2", x2, y2, 3, -4)
    spline("Curva 3", x3, y3, 1/3, -3/2)
    plt.vlines([1, 17, 27.7, 30], 1, 9, color="black", linestyle="dashed")
    plt.grid(True)
    plt.show()


main()