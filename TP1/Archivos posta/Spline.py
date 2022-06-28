import numpy as np
import matplotlib.pyplot as plt
import tabulate as tab
import Descomposicion_LU as lu
import math

# Calculo de matrices

def calcularMatrizA(h, cantidadDeNodos):
    A = [[0 for _ in range(cantidadDeNodos)] for _ in range(cantidadDeNodos)]

    A[0][0] = 2 * h[0]
    A[0][1] = h[0]

    for i in range(1, cantidadDeNodos - 1):
        A[i][i-1] = h[i-1]
        A[i][i] = 2*(h[i-1] + h[i])
        A[i][i+1] = h[i]

    A[cantidadDeNodos-1][cantidadDeNodos-2] = h[cantidadDeNodos-2]
    A[cantidadDeNodos-1][cantidadDeNodos-1] = 2 * h[cantidadDeNodos-2]

    return A

def calcularMatrizB(a, h, cantidadDeNodos, pendienteInicial, pendienteFinal):
    B = [3/h[0] * (a[1] - a[0]) - 3 * pendienteInicial]
    
    for n in range(2, cantidadDeNodos):
        B.append(3/h[n-1] * (a[n] - a[n-1]) - 3/h[n-2] * (a[n-1] - a[n-2]))
    
    B.append(3 * pendienteFinal - 3/h[cantidadDeNodos-2] * (a[cantidadDeNodos-1] - a[cantidadDeNodos-2]))

    return B

# Calculo de coeficientes b y d

def calcularCoeficientesRestantes(a, h, c, cantidadDeNodos):
    b = []
    for j in range(cantidadDeNodos - 1):
        b.append(1/h[j] * (a[j+1] - a[j]) - h[j]/3 * (2*c[j] + c[j+1]))

    d = []
    for j in range(cantidadDeNodos - 1):
        d.append(1/(3 * h[j]) * (c[j+1] - c[j]))

    return b, d

# Main protocol

def spline(x, a, pendienteInicial, pendienteFinal):
    cantidadDeNodos = len(x)
    h = [(x[n+1] - x[n]) for n in range(cantidadDeNodos - 1)]

    A = calcularMatrizA(h, cantidadDeNodos)
    B = calcularMatrizB(a, h, cantidadDeNodos, pendienteInicial, pendienteFinal)

    # Grafico de A y B.
    for x in range(cantidadDeNodos):
        print(A[x])

    print("\n" + str(B))
    print("\n")

    # Resuelvo ecuación
    c = lu.solve(A, B)

    # Calculo coeficientes b, d
    b, d = calcularCoeficientesRestantes(a, h, c, cantidadDeNodos)

    return a, b, c, d

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

# Imprimir tabla de coeficientes de una Spline

def mostrarTablaDeCoeficientes(a, b, c, d, cantidadDeNodos):
    datos = []
    titulos = ['S', 'a', 'b', 'c', 'd ']

    for i in range(cantidadDeNodos-1):
        datos.append(('S_'+str(i), a[i], b[i], c[i], d[i]))

    #ascii
    #print(tab.tabulate(datos, headers=titulos, floatfmt=".16f", tablefmt="github"))
    #fancy
    print(tab.tabulate(datos, headers=titulos, floatfmt=".16f", tablefmt='fancy_grid'))



def main():
    x1 = [1  ,   2,   5,   6,   7,   8,  10,  13,  17]
    y1 = [3.0, 3.7, 3.9, 4.2, 5.7, 6.6, 7.1, 6.7, 4.5]

    x2 = [17 ,  20,  23,  24,  25,  27, 27.7]
    y2 = [4.5, 7.0, 6.1, 5.6, 5.8, 5.2, 4.1]

    x3 = [27.7,  28,  29,  30]
    y3 = [4.1 , 4.3, 4.1, 3.0]

    print("Curva 1")
    a1, b1, c1, d1 = spline(x1, y1, 1, -2/3)
    mostrarTablaDeCoeficientes(a1, b1, c1, d1, len(x1))
    print("\n")
    print("="*100)
    print("\n")

    print("Curva 2")
    a2, b2, c2, d2 = spline(x2, y2, 3, -4)
    mostrarTablaDeCoeficientes(a2, b2, c2, d2, len(x2))
    print("\n")
    print("="*100)
    print("\n")
    
    print("Curva 3")
    a3, b3, c3, d3 = spline(x3, y3, 1/3, -3/2)
    mostrarTablaDeCoeficientes(a3, b3, c3, d3, len(x3))
    print("\n")
    print("="*100)
    print("\n")

    puntosCurva1, resultado1 = obtenerPuntosCurvaSpline(a1, b1, c1, d1, x1)
    puntosCurva2, resultado2 = obtenerPuntosCurvaSpline(a2, b2, c2, d2, x2)
    puntosCurva3, resultado3 = obtenerPuntosCurvaSpline(a3, b3, c3, d3, x3)

    plt.plot(puntosCurva1, resultado1, 'k', lw=1)
    plt.plot(puntosCurva2, resultado2, 'k', lw=1)
    plt.plot(puntosCurva3, resultado3, 'k', lw=1)
    plt.grid(True)
    plt.scatter(x1, y1, color='red')
    plt.scatter(x2, y2, color='red')
    plt.scatter(x3, y3, color='red')

    plt.show()


main()

def extra():
    euler = math.e
    x = [0, 1, 2, 3]
    y = [1, euler, euler**2, euler**3]

    print("Curva de prueba")
    a1, b1, c1, d1 = spline(x, y, 1, euler**3)

    puntosCurva, resultado = obtenerPuntosCurvaSpline(a1, b1, c1, d1, x)

    plt.plot(puntosCurva, resultado, 'k', lw=1)
    plt.grid(True)
    plt.scatter(x, y, color='red')

    plt.show()

#extra()