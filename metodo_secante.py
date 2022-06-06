from math import sin

def metodo_secante(p1, p2, función, cota_de_error, iteración):
    p3 = p2 - ((función(p2) * (p2 - p1)) / (función(p2) - función(p1)))
    iteración += 1
    if(abs(p3 - p2) >= cota_de_error):
        return metodo_secante(p2, p3, función, cota_de_error, iteración)
    else:
        return (p3, iteración)
    
def función1(x):
    return ((x*x) / 4) - sin(x)

p3, iteracion = metodo_secante(1.1, 3.5, función1, 0.001, 0)

print("La raíz es "+ str(p3) + " con " + str(iteracion) + " iteraciones")