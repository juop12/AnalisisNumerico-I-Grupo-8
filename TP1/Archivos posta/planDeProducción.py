import numpy as np
import matplotlib.pyplot as plt
import math as m
import tabulate as tab


#----------------------------------------------------------------------------------------------------
# Funciones de las que haremos uso

def funcionProduccion(x): 
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
    p_comparacion = p_1
    iteraciones.append(p_2)
    
    i = 1
    while abs(p_2 - p_comparacion) > tolerancia and i < maxIteraciones:
        p_comparacion = p_2
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

    while abs(p_1 - p_0) > tolerancia and (not m.isclose(funcionDerivada(p_1), 0)) and i < maxIteraciones:
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
    """ Recibe las iteraciones de un método de busqueda de raíces 
        y devuelve una lista con los órdenes de convergencia"""
    ordenDeConvergencia = []
    tolerancia = 10**-10
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

def calcularConstanteAsintotica(iteraciones, convergencia):
    """ Recibe las iteraciones de un método de busqueda de raíces 
        y su convergencia. Devuelve una lista con las constantes 
        asintóticas."""
    constantesAsintoticas = []
    tolerancia = 10**-10
    constante = 0
    for indice in range(len(iteraciones)):
        if indice < 2:
            constantesAsintoticas.append(constante)

        else:
            errorParActual   = abs(iteraciones[indice]     - iteraciones[indice - 1])
            errorParAnterior = abs(iteraciones[indice - 1] - iteraciones[indice - 2])
            if errorParActual > tolerancia and (errorParAnterior ** convergencia) > tolerancia:
                constante = errorParActual / (errorParAnterior**convergencia)
            constantesAsintoticas.append(constante)

    return constantesAsintoticas

def obtenerListaErrorIteraciones(iteraciones):
    """ Recibe las iteraciones de un método de busqueda de raíces 
        y devuelve una lista con los errores de cada iteración."""
    errores = []
    errores.append(0)
    for i in range(1, len(iteraciones)):
        errores.append(abs(iteraciones[i] - iteraciones[i - 1]))
    return errores

#-------------------------------------------------------------------------------------------
# Función para imprimir los resultados

def mostrarComparacionIteraciones(iteracionesBiseccion, iteracionesPuntoFijo, iteracionesNewtonRaphson, iteracionesNewtonRaphsonModificado, iteracionesSecante, tolerancia):
    """ Imprime un gráfico con la evolución de la posible raíz a través de las 
        iteraciones de los diferentes métodos de busqueda de raíces"""
    plt.plot(range(1, len(iteracionesBiseccion) + 1), iteracionesBiseccion, label='Bisección', color='b')
    plt.plot(range(1, len(iteracionesPuntoFijo) + 1), iteracionesPuntoFijo, label='Punto fijo', color='y')
    plt.plot(range(1, len(iteracionesNewtonRaphson) + 1), iteracionesNewtonRaphson, label='Newton Raphson', color='g')
    plt.plot(range(1, len(iteracionesNewtonRaphsonModificado) + 1), iteracionesNewtonRaphsonModificado, label='Newton Raphson Modificado', color='m')
    plt.plot(range(1, len(iteracionesSecante) + 1), iteracionesSecante, label='Secante', color='r')
    plt.grid(True)
    plt.title('Comparación iteraciones de los distintos métodos iterativos con tolerancia' + str(tolerancia))
    plt.xlabel('Iteración')
    plt.ylabel('Aproximación raíz')
    plt.legend()
    plt.show()

def mostrarComparacionConvergencia(convergenciaBiseccion, convergenciaPuntoFijo, convergenciaNewtonRaphson, convergenciaNewtonRaphsonModificado, convergenciaSecante, tolerancia):
    """ Imprime un gráfico con la evolución de la convergencia a través de las 
        iteraciones de los diferentes métodos de busqueda de raíces"""
    plt.plot(range(1, len(convergenciaBiseccion) + 1), convergenciaBiseccion, label='Bisección', color='b')
    plt.plot(range(1, len(convergenciaPuntoFijo) + 1), convergenciaPuntoFijo, label='Punto fijo', color='y')
    plt.plot(range(1, len(convergenciaNewtonRaphson) + 1), convergenciaNewtonRaphson, label='Newton Raphson', color='g')
    plt.plot(range(1, len(convergenciaNewtonRaphsonModificado) + 1), convergenciaNewtonRaphsonModificado, label='Newton Raphson Modificado', color='m')
    plt.plot(range(1, len(convergenciaSecante) + 1), convergenciaSecante, label='Secante', color='r')
    plt.grid(True)
    plt.ylim(0, 3)
    plt.title('Comparación órdenes de convergencia de los distintos métodos iterativos con tolerancia' + str(tolerancia))
    plt.xlabel('Iteración')
    plt.ylabel('Orden de convergencia')
    plt.legend()
    plt.show()

def mostrarComparacionConstanteAsintotica(constanteAsintoticaBiseccion, constanteAsintoticaPuntoFijo, constanteAsintoticaNewtonRaphson, constanteAsintoticaNewtonRaphsonModificado, constanteAsintoticaSecante, tolerancia):
    """ Imprime un gráfico con la evolución de la constante asintótica a través de las 
        iteraciones de los diferentes métodos de busqueda de raíces"""
    plt.plot(range(1, len(constanteAsintoticaBiseccion) + 1), constanteAsintoticaBiseccion, label='Bisección', color='b')
    plt.plot(range(1, len(constanteAsintoticaPuntoFijo) + 1), constanteAsintoticaPuntoFijo, label='Punto fijo', color='y')
    plt.plot(range(1, len(constanteAsintoticaNewtonRaphson) + 1), constanteAsintoticaNewtonRaphson, label='Newton Raphson', color='g')
    plt.plot(range(1, len(constanteAsintoticaNewtonRaphsonModificado) + 1), constanteAsintoticaNewtonRaphsonModificado, label='Newton Raphson Modificado', color='m')
    plt.plot(range(1, len(constanteAsintoticaSecante) + 1), constanteAsintoticaSecante, label='Secante', color='r')
    plt.grid(True)
    plt.ylim(0, 1)
    plt.title('Comparación constantes asintóticas de los distintos métodos iterativos con tolerancia' + str(tolerancia))
    plt.xlabel('Iteración')
    plt.ylabel('Constante Asintótica')
    plt.legend()
    plt.show()

def mostrarComparacionErrores(erroresBiseccion, erroresPuntoFijo, erroresNewtonRaphson, erroresNewtonRaphsonModificado, erroresSecante, tolerancia):
    """ Imprime un gráfico con la evolución del error a través de las 
        iteraciones de los diferentes métodos de busqueda de raíces"""
    plt.plot(range(1, len(erroresBiseccion) + 1), erroresBiseccion, label='Bisección', color='b')
    plt.plot(range(1, len(erroresPuntoFijo) + 1), erroresPuntoFijo, label='Punto fijo', color='y')
    plt.plot(range(1, len(erroresNewtonRaphson) + 1), erroresNewtonRaphson, label='Newton Raphson', color='g')
    plt.plot(range(1, len(erroresNewtonRaphsonModificado) + 1), erroresNewtonRaphsonModificado, label='Newton Raphson Modificado', color='m')
    plt.plot(range(1, len(erroresSecante) + 1), erroresSecante, label='Secante', color='r')
    plt.grid(True)
    plt.title('Comparación log del error de cada iteración con tolerancia ' + str(tolerancia))
    plt.xlabel('Iteración')
    plt.ylabel('log(error)')
    plt.legend()
    plt.show()

def mostrarIteraciones(nombreMetodo, iteraciones, ordenDeConvergencia, constanteAsintotica, tolerancia):
    """ Imprime una tabla con la evolución de las posibles raíces, la convergencia y la
        constante asintótica de un método de búsqueda de raíces"""
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

#-------------------------------------------------------------------------------------------
# Main

def main():
    tolerancia = [1e-5, 1e-13]
    maxIteraciones = 100

    for t in tolerancia:
        iteracionesBiseccion = metodoDeBiseccion(funcionProduccion, 1000, 1200, t, maxIteraciones)
        erroresBiseccion = obtenerListaErrorIteraciones(iteracionesBiseccion)
        convergenciaBiseccion = calcularOrdenDeConvergencia(iteracionesBiseccion)
        constanteAsintoticaBiseccion = calcularConstanteAsintotica(iteracionesBiseccion, convergenciaBiseccion[-1])

        mostrarIteraciones("Método Bisección", iteracionesBiseccion, convergenciaBiseccion, constanteAsintoticaBiseccion, t)

        iteracionesPuntoFijo = metodoDePuntoFijo(funcionProduccionPuntoFijo, 1100, t, maxIteraciones)
        erroresPuntoFijo = obtenerListaErrorIteraciones(iteracionesPuntoFijo)
        convergenciaPuntoFijo = calcularOrdenDeConvergencia(iteracionesPuntoFijo)
        constanteAsintoticaPuntoFijo = calcularConstanteAsintotica(iteracionesPuntoFijo, convergenciaPuntoFijo[-1])

        mostrarIteraciones("Método Punto Fijo", iteracionesPuntoFijo, convergenciaPuntoFijo, constanteAsintoticaPuntoFijo, t)

        iteracionesNewtonRaphson = metodoDeNewtonRaphson(funcionProduccion, funcionProduccionDerivada, 1100, t, maxIteraciones)
        erroresNewtonRaphson = obtenerListaErrorIteraciones(iteracionesNewtonRaphson)
        convergenciaNewtonRaphson = calcularOrdenDeConvergencia(iteracionesNewtonRaphson)
        constanteAsintoticaNewtonRaphson = calcularConstanteAsintotica(iteracionesNewtonRaphson, convergenciaNewtonRaphson[-1])
        
        mostrarIteraciones("Método Newton-Rapshon", iteracionesNewtonRaphson, convergenciaNewtonRaphson, constanteAsintoticaNewtonRaphson, t)

        iteracionesNewtonRaphsonModificado = metodoDeNewtonRaphsonModificado(funcionProduccion, funcionProduccionDerivada, funcionProduccionDerivadaSegunda, 1050, t, maxIteraciones)
        erroresNewtonRaphsonModificado = obtenerListaErrorIteraciones(iteracionesNewtonRaphsonModificado)
        convergenciaNewtonRaphsonModificado = calcularOrdenDeConvergencia(iteracionesNewtonRaphsonModificado)
        constanteAsintoticaNewtonRaphsonModificado = calcularConstanteAsintotica(iteracionesNewtonRaphsonModificado, convergenciaNewtonRaphsonModificado[-1])

        mostrarIteraciones("Método Newton-Raphson Modificado", iteracionesNewtonRaphsonModificado, convergenciaNewtonRaphsonModificado, constanteAsintoticaNewtonRaphsonModificado, t)

        iteracionesSecante = metodoDeLaSecante(funcionProduccion, 1000, 1200, t, maxIteraciones)
        erroresSecante = obtenerListaErrorIteraciones(iteracionesSecante)
        convergenciaSecante = calcularOrdenDeConvergencia(iteracionesSecante)
        constanteAsintoticaSecante = calcularConstanteAsintotica(iteracionesSecante, convergenciaSecante[-1])

        mostrarIteraciones("Método de la Secante", iteracionesSecante, convergenciaSecante, constanteAsintoticaSecante, t)

        mostrarComparacionIteraciones(iteracionesBiseccion, iteracionesPuntoFijo, iteracionesNewtonRaphson, iteracionesNewtonRaphsonModificado, iteracionesSecante, t)
        mostrarComparacionErrores(erroresBiseccion, erroresPuntoFijo, erroresNewtonRaphson, erroresNewtonRaphsonModificado, erroresSecante, t)
        mostrarComparacionConvergencia(convergenciaBiseccion, convergenciaPuntoFijo, convergenciaNewtonRaphson, convergenciaNewtonRaphsonModificado, convergenciaSecante, t)
        mostrarComparacionConstanteAsintotica(constanteAsintoticaBiseccion, constanteAsintoticaPuntoFijo, constanteAsintoticaNewtonRaphson, constanteAsintoticaNewtonRaphsonModificado, constanteAsintoticaSecante, t)

main()
