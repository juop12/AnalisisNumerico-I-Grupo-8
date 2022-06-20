import numpy as np
import matplotlib.pyplot as plt
import math as m
import tabulate as tab


def f(x):
    return (0.001 * x * ((x-1000)**2)) - 25000


def biseccion(funcion, cotaMin, cotaMax, tolerancia, maxIteraciones):
    if funcion(cotaMin) * funcion(cotaMax) >= 0:
        print("No hay cambio de signo en los extremos")
        return None

    historia = []
    p_anterior = cotaMin

    i = 1
    while i < maxIteraciones:
        p_candidata = (cotaMin + cotaMax) / 2
        historia.append((i, p_candidata))

        if np.abs(funcion(p_candidata) - funcion(p_anterior)) < tolerancia:
            return p_candidata, i, historia

        if funcion(cotaMin) * funcion(p_candidata) > 0:
            cotaMin = p_candidata
        else:
            cotaMax = p_candidata

        i += 1
        p_anterior = p_candidata

    print("No se alcanzo la tolerancia")
    return None


def calcularOrdenDeConvergencia(historia):
    ordenDeConvergencia = []
    tolerancia = 10**-16
    alfa = 0
    for indice in range(0, len(historia) - 1):
        if indice < 3:
            ordenDeConvergencia.append(alfa)
        else:
            error_nMas1 = abs(historia[indice][1] - historia[indice - 1][1])
            error_n = abs(historia[indice - 1][1] - historia[indice - 2][1])
            error_nMenos1 = abs(historia[indice - 2][1] - historia[indice - 3][1])
            """print("X n:", historia[indice][1], "Xn - 1:", historia[indice - 1][1], "\nXn-2:", historia[indice - 2][1], "Xn-3:", historia[indice - 3][1])
            print("\nXn - Xn-1 = Error n+1:", error_nMas1)
            print("Xn-1 - Xn-2 = Error n:", error_n)
            print("Xn-1 - Xn-2 = Error n-1:", error_nMenos1)"""
            if error_nMas1 > tolerancia or error_n > tolerancia or error_nMenos1 > tolerancia:
                logNumerador = abs(m.log(error_nMas1 / error_n))
                logDenominador = abs(m.log(error_n / error_nMenos1))
                """print("\nLog numerador = log(error n+1 / error n:", logNumerador)
                print("Log denominador = log(error n / errorr n-1:", logDenominador)"""
                alfa = logNumerador/logDenominador
                # print("\nAlfa:", alfa)
                ordenDeConvergencia.append(alfa)

            else:
                print("Errores muy chicos")
                ordenDeConvergencia.append(ordenDeConvergencia[-1])

    return ordenDeConvergencia


def calcularConstanteAsintotica(historia, alfa):
    constantesAsintoticas = []
    for indice in range(0, len(historia) - 1):
        if indice < 2:
            constantesAsintoticas.append(0)

        else:
            error_nMas1 = abs(historia[indice][1] - historia[indice - 1][1])
            error_n = abs(historia[indice - 1][1] - historia[indice - 2][1])
            constante = error_nMas1 / (error_n**alfa)
            constantesAsintoticas.append(constante)

    return constantesAsintoticas


def mostrarIteraciones(historia, ordenDeConvergencia, constanteAsintotica, tolerancia):
    """if historia:
        print("------------------------------------------------------------------------")
        print("Interaccion", "\t", "Raiz", "\t\t\t\t", "Convergencia", "\t", "Tolerancia: ", tolerancia)
        print("------------------------------------------------------------------------")

        if len(historia) <= 12:
            for x in range(len(historia)):
                print(historia[x][0], "\t\t\t\t", historia[x][1],  "\t\t\t\t\t\t", ordenDeConvergencia[x])
        else:
            for x in range(5):
                print(historia[x][0], "\t\t\t\t", historia[x][1],  "\t\t\t\t\t\t", ordenDeConvergencia[x])
            historiaInversa = historia[-5:]
            convergenciaInversa = ordenDeConvergencia[-5:]
            for x in range(5):
                print(historiaInversa[x][0], "\t\t\t\t", historiaInversa[x][1],  "\t\t\t\t\t\t", convergenciaInversa[x])"""

    if historia:
        datos = []
        titulos = ['Iteracion', 'Raiz', 'Convergencia', 'Constante Asintotica']
        if len(historia) <= 12:
            for x in range(len(historia)):
                datos.append((historia[x][0], historia[x][1], ordenDeConvergencia[x], constanteAsintotica[x]))

        else:
            for x in range(5):
                datos.append((historia[x][0], historia[x][1], ordenDeConvergencia[x], constanteAsintotica[x]))
            historiaInversa = historia[-5:]
            convergenciaInversa = ordenDeConvergencia[-5:]
            constanteInversa = constanteAsintotica[-5:]
            for x in range(5):
                datos.append((historiaInversa[x][0], historiaInversa[x][1], convergenciaInversa[x], constanteInversa[x]))

        print("\n\nTolerancia: ", tolerancia)
        print(tab.tabulate(datos, headers=titulos, floatfmt=".16f", tablefmt="github"))


def imprimirFuncion(desde, hasta, raiz):
    puntos = np.linspace(desde, hasta, 250)
    y = f(puntos)
    plt.plot(puntos, y, 'k', lw=1)
    plt.title('Utilidad Unitaria por kg')
    plt.grid(True)
    plt.xlabel('Kg producto')
    plt.ylabel('Utilidad unitaria')
    plt.scatter(raiz, f(raiz), color='red')
    plt.show()


def main():
    a = 827
    b = 1800
    tolerancia = 1e-5
    maxIteraciones = 100

    historia = biseccion(f, a, b, tolerancia, maxIteraciones)[2]
    mostrarIteraciones(historia, calcularOrdenDeConvergencia(historia), calcularConstanteAsintotica(historia, 1), tolerancia)

    tolerancia = 1e-13
    raiz, _, historia = biseccion(f, a, b, tolerancia, maxIteraciones)
    mostrarIteraciones(historia, calcularOrdenDeConvergencia(historia), calcularConstanteAsintotica(historia, 1), tolerancia)

    imprimirFuncion(600, b, raiz)


main()
