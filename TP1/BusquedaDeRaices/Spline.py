import numpy as np


def spline(nodos):
    h = []
    for n in range(len(nodos) - 1):
        h.append(nodos[n+1][0] - nodos[n][0])

    print(h)
    print("\n")
    A = [[0 for _ in range(len(nodos) + 1)] for _ in range(len(nodos) + 1)]

    A[0][0] = 1
    for i in range(1, len(nodos)):
        A[i][i-1] = h[i-2]
        A[i][i] = 2*(h[i-2] + h[i-1])
        A[i][i+1] = h[i-1]
    A[len(nodos)][len(nodos)] = 1

    b = [0]
    for n in range(1, len(nodos)):
        b.append(3/h[n-1]*(nodos[n][1] - nodos[n - 1][1]) - 3/h[n-2]*(nodos[n-1][1] - nodos[n-2][1]))
    b += [0]

    for x in range(len(nodos)+1):
        print(A[x])

    print("\n" + str(b))

    coeficientes = np.linalg.solve(A, b)
    print("\n")
    print(coeficientes)

def main():
    nodosCurva1 = [(1, 3.0), (2, 3.7), (5, 3.9), (6, 4.2), (7, 5.7), (8, 6.6), (10, 7.1), (13, 6.7), (17, 4.5)]
    nodosCurva2 = [(17, 4.5), (20, 7.0), (23, 6.1), (24, 5.6), (25, 5.8), (27, 5.2), (27.7, 4.1)]
    nodosCurva3 = [(27.7, 4.1), (28, 4.3), (29, 4.1), (30, 3.0)]
    #derivadaExtremosCurva1 = [1, -2 / 3]
    #derivadaExtremosCurva2 = [3, -4]
    #derivadaExtremosCurva3 = [1 / 3, -3 / 2]
    print("Curva 1")
    spline(nodosCurva1)
    print("Curva 2")
    spline(nodosCurva2)
    print("Curva 3")
    spline(nodosCurva3)


main()
