from math import sin
from pkgutil import ImpImporter

def función1(x):
    return ((x*x) / 4) - sin(x)

def imprimir_tabla_de_raices(posibles_raices):
	print("    n     |     p_n    \n--------------------")
	numero_de_iteracion = 0
	for raiz in posibles_raices:
		print("   ",numero_de_iteracion,"    |    ",raiz) 
		numero_de_iteracion += 1


def metodo_de_la_secante(semilla_1 , semilla_2, función, cota_de_error):
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

    return posibles_raices

posibles_raices = metodo_de_la_secante(1.1, 3.5, función1, 0.00000001)

imprimir_tabla_de_raices(posibles_raices)
        
