#import numpy as np
import matplotlib.pyplot as plt
import math as m
import tabulate as tab

#----------------------------------------------------------------------------------------------------
# Funciones de las que haremos uso

def funcionProduccion(x):
    #return (0.001 * x * ((x-1000)**2)) - 25000  
    return 0.001 * (x**3) - 2 * (x**2) + 1000 * x - 25000  

def funcionProduccionDerivada(x):
    return 0.003 * (x**2) - 4 * x + 1000

def funcionProduccionDerivadaSegunda(x):
    return 0.006 * x - 4

def funcionProduccionPuntoFijo(x):
    return 25 - 0.000001 * (x ** 3) + 0.002 * (x ** 2)

#----------------------------------------------------------------------------------------------------
# Metodos para la búsqueda de raices

def metodoDeBiseccion(funcion, semilla_1, semilla_2, tolerancia, maxIteraciones):
    """ Devuelve la raíz hallada a través del método de la bisección
        y un vector con todas las iteraciones de este."""
    iteraciones = [] 
    p_0 = semilla_1
    p_1 = semilla_2
    p_2 = (p_0 + p_1) / 2
    iteraciones.append(p_2)
    
    i = 1
    while abs(p_2 - p_1) > tolerancia and i < maxIteraciones:
        if funcion(p_0) * funcion(p_2) > 0:
            p_0 = p_2
        else:
            p_1 = p_2
        p_2 = (p_0 + p_1) / 2
        iteraciones.append(p_2)
        i += 1

    return iteraciones

def metodoDePuntoFijo(funcion, semilla, tolerancia, maxIteraciones):
    """ Devuelve la raíz hallada a través del método del punto fijo
        y un vector con todas las iteraciones de este."""
    iteraciones = []
    p_0 = semilla
    p_1 = funcion(p_0)
    iteraciones.append(p_1)
    i = 0

    while abs(p_1 - p_0) > tolerancia and i < maxIteraciones:
        p_0 = p_1
        p_1 = funcion(p_0)
        iteraciones.append(p_1)
        i += 1
    
    return iteraciones

def metodoDeNewtonRaphson(funcion, funcionDerivada, semilla, tolerancia, maxIteraciones):
    """ Devuelve la raíz hallada a través del método de Newton-Raphson
        y un vector con todas las iteraciones de este."""
    iteraciones = []
    p_0 = semilla
    p_1 = p_0 - funcion(p_0)/funcionDerivada(p_0)
    iteraciones.append(p_1)
    i = 0

    while abs(p_1 - p_0) > tolerancia and i < maxIteraciones:
        p_0 = p_1
        p_1 = p_0 - funcion(p_0)/funcionDerivada(p_0)
        iteraciones.append(p_1)
        i += 1
    
    return iteraciones

def metodoDeNewtonRaphsonModificado(funcion, funcionDerivada, funcionDerivadaSegunda, semilla, tolerancia, maxIteraciones):
    """ Devuelve la raíz hallada a través del método de Newton-Raphson modificado
        y un vector con todas las iteraciones de este."""
    iteraciones = []
    p_0 = semilla
    p_1 = p_0 - (funcion(p_0) * funcionDerivada(p_0))/((funcionDerivada(p_0))**2 - funcion(p_0) * funcionDerivadaSegunda(p_0))
    iteraciones.append(p_1)
    i = 0

    while abs(p_1 - p_0) > tolerancia and i < maxIteraciones:
        p_0 = p_1
        p_1 = p_0 - (funcion(p_0) * funcionDerivada(p_0))/((funcionDerivada(p_0))**2 - funcion(p_0) * funcionDerivadaSegunda(p_0))
        iteraciones.append(p_1)
        i += 1

    return iteraciones

def metodoDeLaSecante(funcion, semilla_1, semilla_2, tolerancia, maxIteraciones):
    """ Devuelve la raíz hallada a través del método de la secante
        secante y un vector con todas las iteraciones de este."""
    iteraciones = []
    p_0 = semilla_1
    p_1 = semilla_2
    p_2 = p_1 - ((funcion(p_1) * (p_1 - p_0)) / (funcion(p_1) - funcion(p_0)))
    iteraciones.append(p_2)
    i = 1

    while abs(p_2 - p_1) > tolerancia and i < maxIteraciones:
        p_0 = p_1
        p_1 = p_2
        p_2 = p_1 - ((funcion(p_1) * (p_1 - p_0)) / (funcion(p_1) - funcion(p_0)))
        iteraciones.append(p_2)

    return iteraciones

#-------------------------------------------------------------------------------------------
# Funciones de convergencia y constante asintótica

def calcularOrdenDeConvergencia(iteraciones):
    ordenDeConvergencia = []
    tolerancia = 10**-16
    alfa = 0
    for indice in range(len(iteraciones)):
        if indice < 3:
            ordenDeConvergencia.append(alfa)
        else:
            errorParMasReciente = abs(iteraciones[indice]     - iteraciones[indice - 1])
            errorParAnterior    = abs(iteraciones[indice - 1] - iteraciones[indice - 2])
            errorParMasAnterior = abs(iteraciones[indice - 2] - iteraciones[indice - 3])

            if errorParMasReciente > tolerancia and errorParAnterior > tolerancia and errorParMasAnterior > tolerancia and \
               errorParMasReciente != errorParAnterior and errorParAnterior != errorParMasAnterior:

                logNumerador   = abs(m.log(errorParMasReciente / errorParAnterior))
                logDenominador = abs(m.log(errorParAnterior    / errorParMasAnterior))
                alfa = logNumerador/logDenominador
                ordenDeConvergencia.append(alfa)
            else:
                ordenDeConvergencia.append(ordenDeConvergencia[-1])
                
    return ordenDeConvergencia

def calcularConstanteAsintotica(iteraciones, alfa):
    constantesAsintoticas = []
    tolerancia = 10**-12
    constante = 0
    for indice in range(len(iteraciones)):
        if indice < 2:
            constantesAsintoticas.append(constante)

        else:
            errorParActual   = abs(iteraciones[indice]     - iteraciones[indice - 1])
            errorParAnterior = abs(iteraciones[indice - 1] - iteraciones[indice - 2])
            if errorParActual > tolerancia and (errorParAnterior ** alfa) > tolerancia:
                constante = errorParActual / (errorParAnterior**alfa)
            constantesAsintoticas.append(constante)

    return constantesAsintoticas

#-------------------------------------------------------------------------------------------
# Función para imprimir los resultados

def imprimirOrdenDeConvergencia(nombreMetodo, tolerancia, ordenDeConvergencia):
    plt.plot(range(1, len(ordenDeConvergencia) + 1), ordenDeConvergencia)
    plt.grid(True)
    plt.title("Orden de Convergencia " + nombreMetodo + " con tolerancia " + str(tolerancia))
    plt.xlabel('Iteración')
    plt.ylabel('Orden de Convergencia')
    plt.ylim(0, 3)
    plt.show()

def imprimirConstanteAsintotica(nombreMetodo, tolerancia, constanteAsintotica):
    plt.plot(constanteAsintotica)
    plt.grid(True)
    plt.title("Constante Asintótica " + nombreMetodo + " con tolerancia " + str(tolerancia))
    plt.xlabel('Iteración')
    plt.ylabel('Constante Asintótica')
    plt.ylim(0, 2)
    plt.show()


def mostrarIteraciones(nombreMetodo, iteraciones, ordenDeConvergencia, constanteAsintotica, tolerancia):
    datos = []
    cantidadIteraciones = len(iteraciones)
    titulos = ['Iteracion', 'Raiz', 'Convergencia', 'Constante Asintotica']

    if cantidadIteraciones <= 12:
        for i in range(len(iteraciones)):
            datos.append((i + 1, iteraciones[i], ordenDeConvergencia[i], constanteAsintotica[i]))

    else:
        for i in range(5):
            datos.append((i + 1, iteraciones[i], ordenDeConvergencia[i], constanteAsintotica[i]))

        iteracionesInversa = iteraciones[-5:]
        convergenciaInversa = ordenDeConvergencia[-5:]
        constanteInversa = constanteAsintotica[-5:]

        for i in range(5):
            datos.append((cantidadIteraciones - 4 + i, iteracionesInversa[i], convergenciaInversa[i], constanteInversa[i]))

    print("\n\n" + nombreMetodo)
    print("Tolerancia: ", tolerancia)
    print(tab.tabulate(datos, headers=titulos, floatfmt=".16f", tablefmt="github"))
    #imprimirOrdenDeConvergencia(nombreMetodo, tolerancia, ordenDeConvergencia)
    #imprimirConstanteAsintotica(nombreMetodo, tolerancia, constanteAsintotica)

#-------------------------------------------------------------------------------------------
# Main

def main():
    tolerancia = [1e-5, 1e-13]
    maxIteraciones = 100

    for t in tolerancia:
        
        iteracionesBiseccion = metodoDeBiseccion(funcionProduccion, 827, 1800, t, maxIteraciones)
        convergenciaBiseccion = calcularOrdenDeConvergencia(iteracionesBiseccion)
        constanteAsintoticaBiseccion = calcularConstanteAsintotica(iteracionesBiseccion, convergenciaBiseccion[-1])

        mostrarIteraciones("Método Bisección", iteracionesBiseccion, convergenciaBiseccion, constanteAsintoticaBiseccion, t)


        iteracionesPuntoFijo = metodoDePuntoFijo(funcionProduccionPuntoFijo, 1000, t, maxIteraciones)
        convergenciaPuntoFijo = calcularOrdenDeConvergencia(iteracionesPuntoFijo)
        constanteAsintoticaPuntoFijo = calcularConstanteAsintotica(iteracionesPuntoFijo, convergenciaPuntoFijo[-1])

        mostrarIteraciones("Método Punto Fijo", iteracionesPuntoFijo, convergenciaPuntoFijo, constanteAsintoticaPuntoFijo, t)


        iteracionesNewtonRaphson = metodoDeNewtonRaphson(funcionProduccion, funcionProduccionDerivada, 1100, t, maxIteraciones)
        convergenciaNewtonRaphson = calcularOrdenDeConvergencia(iteracionesNewtonRaphson)
        constanteAsintoticaNewtonRaphson = calcularConstanteAsintotica(iteracionesNewtonRaphson, convergenciaNewtonRaphson[-1])
        
        mostrarIteraciones("Método Newton-Rapshon", iteracionesNewtonRaphson, convergenciaNewtonRaphson, constanteAsintoticaNewtonRaphson, t)


        iteracionesNewtonRaphsonModificado = metodoDeNewtonRaphsonModificado(funcionProduccion, funcionProduccionDerivada, funcionProduccionDerivadaSegunda, 1050, t, maxIteraciones)
        convergenciaNewtonRaphsonModificado = calcularOrdenDeConvergencia(iteracionesNewtonRaphsonModificado)
        constanteAsintoticaNewtonRaphsonModificado = calcularConstanteAsintotica(iteracionesNewtonRaphsonModificado, convergenciaNewtonRaphsonModificado
        [-1])

        mostrarIteraciones("Método Newton-Raphson Modificado", iteracionesNewtonRaphsonModificado, convergenciaNewtonRaphsonModificado, constanteAsintoticaNewtonRaphsonModificado, t)


        iteracionesSecante = metodoDeLaSecante(funcionProduccion, 1000, 1800, t, maxIteraciones)
        convergenciaSecante = calcularOrdenDeConvergencia(iteracionesSecante)
        constanteAsintoticaSecante = calcularConstanteAsintotica(iteracionesSecante, convergenciaSecante[-1])

        mostrarIteraciones("Método de la Secante", iteracionesSecante, convergenciaSecante, constanteAsintoticaSecante, t)


main()
