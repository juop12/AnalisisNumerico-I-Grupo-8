from math import sin, cos
from pkgutil import ImpImporter

def función1(x):
    return ((x*x) / 4) - sin(x)

def función2(x):
    return x/2 - cos(x)

def función3(x):
    return 1/2 + sin(x)


def imprimir_tabla_de_raices(nombre_del_método, posibles_raices):
    """Recibe una lista de posibles raices y las imprime en forma de tabla"""

    print("Tabla de iteraciones de la busqueda de raíz a través del ", nombre_del_método, ":")

    numeros_de_iteracion = [str(numero) for numero in range(len(posibles_raices))]
    posibles_raices_strings = [str(raiz) for raiz in posibles_raices]

    lista_de_longitudes = [len(raiz) for raiz in posibles_raices_strings]
    ancho = max(lista_de_longitudes)

    print("|".join(["n".center(ancho+2), "p_n".center(ancho+2)]))
    print("-" * (ancho * 2 + 4))
    
    for iteración in zip(numeros_de_iteracion, posibles_raices_strings):
        línea = "|".join(elemento.center(ancho + 2) for elemento in iteración)
        print(línea)


def metodo_de_la_secante(semilla_1 , semilla_2, función, cota_de_error):
    """ Recibe dos semillas, una función y una cota de error.
        Devuelve la raíz hallada a través del método de la
        secante y un vector con todas las iteraciones de este."""
    posibles_raices = []
    posibles_raices.append(semilla_1)
    posibles_raices.append(semilla_2)

    p_0 = semilla_1
    p_1 = semilla_2

    p_2 = p_1 - ((función(p_1) * (p_1 - p_0)) / (función(p_1) - función(p_0)))
    posibles_raices.append(p_2)

    while (abs(p_2 - p_1) >= cota_de_error):
        p_0 = p_1
        p_1 = p_2
        p_2 = p_1 - ((función(p_1) * (p_1 - p_0)) / (función(p_1) - función(p_0)))
        posibles_raices.append(p_2)

    return p_2, posibles_raices

def metodo_de_Newton_Raphson(semilla, función, función_derivada, cota_de_error):
    """ Recibe una semilla, una función, su derivada y una cota de error.
        Devuelve la raíz hallada a través del método de Newton-Raphson
        y un vector con todas las iteraciones de este."""

    posibles_raices = []
    posibles_raices.append(semilla)

    p_0 = semilla
    p_1 = p_0 - función(p_0)/función_derivada(p_0)
    posibles_raices.append(p_1)

    while(abs(p_1 - p_0) >= cota_de_error):
        p_0 = p_1
        p_1 = p_0 - función(p_0)/función_derivada(p_0)
        posibles_raices.append(p_1)
    
    return p_1, posibles_raices

def metodo_de_Newton_Raphson_modificado(semilla, función, función_derivada, función_derivada_segunda, cota_de_error):
    posibles_raices = []
    posibles_raices.append(semilla)

    p_0 = semilla
    p_1 = p_0 - (función(p_0) * función_derivada(p_0))/((función_derivada(p_0))**2 - función(p_0) * función_derivada_segunda(p_0))
    posibles_raices.append(p_1)

    while(abs(p_1 - p_0) >= cota_de_error):
        p_0 = p_1
        p_1 = p_0 - (función(p_0) * función_derivada(p_0))/((función_derivada(p_0))**2 - función(p_0) * función_derivada_segunda(p_0))
        posibles_raices.append(p_1)

    return p_1, posibles_raices



def main():
    raiz, posibles_raices = metodo_de_la_secante(1.1, 3.5, función1, 0.00000000000000001)
    imprimir_tabla_de_raices("metodo de la secante",posibles_raices)

    raiz, posibles_raices = metodo_de_Newton_Raphson(2.6, función1, función2, 0.00000000000000001)
    imprimir_tabla_de_raices("metodo de Newton-Raphson", posibles_raices)

    raiz, posibles_raices = metodo_de_Newton_Raphson_modificado(2.6, función1, función2, función3, 0.000000000000000000001)
    imprimir_tabla_de_raices("metodo de Newton-Raphson", posibles_raices)

main()

