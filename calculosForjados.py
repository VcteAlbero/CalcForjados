import numpy as np
import matplotlib.pyplot as plt

def clapeyron(luces, flag, carga_sup, carga_punt, x_punt, libre_izq, libre_der, cont):
    n_luces = len(luces)
    mapoyos = [0] * (n_luces+1)
    area = [0] * n_luces
    a = [0] * n_luces
    b = [0] * n_luces

    # FORJADO ISOSTÁTICO
    if cont == 0:
        if libre_izq == 1:
            if flag == 1:  # Carga superficial
                mapoyos[1] = carga_sup[0] * (luces[0] ** 2) / 2
            else:  # Carga puntual
                mapoyos[1] = carga_punt[0] * (luces[0] - x_punt[0])
        
        if libre_der == 1:
            if flag == 1:  # Carga superficial
                mapoyos[-2] = carga_sup[-1] * (luces[-1] ** 2) / 2
            else:  # Carga puntual
                mapoyos[-2] = carga_punt[-1] * x_punt[-1]
    
    # FORJADO CONTINUO
    else:
        long_matriz = n_luces - 1
        cont_m = 0
        # Momentos de voladizos si existen
        if libre_izq == 1:
            if flag == 1:  # Carga superficial
                mapoyos[1] = carga_sup[0] * (luces[0] ** 2) / 2
            else:  # Carga puntual
                mapoyos[1] = carga_punt[0] * (luces[0] - x_punt[0])
            long_matriz = long_matriz - 1
            cont_m = 1     
        if libre_der == 1:
            if flag == 1:  # Carga superficial
                mapoyos[-2] = carga_sup[-1] * (luces[-1] ** 2) / 2
            else:  # Carga puntual
                mapoyos[-2] = carga_punt[-1] * x_punt[-1]
            long_matriz = long_matriz - 1

        # Cálculo de áreas y centros de gravedad
        for i in range(n_luces):
            if flag == 1:  # Carga superficial
                area[i] = 2 / 3 * luces[i] * carga_sup[i] * (luces[i] ** 2) / 8
                a[i] = luces[i] / 2
                b[i] = luces[i] / 2
            else:  # Carga puntual
                area[i] = carga_punt[i] * x_punt[i] * (luces[i] - x_punt[i]) / 2
                a[i] = (luces[i] + x_punt[i]) / 3
                b[i] = luces[i] - a[i]
        if long_matriz > 0:
        # Matrices [A], [B] y resolución
            A = np.zeros((long_matriz, long_matriz))
            B = np.zeros((long_matriz, 1))

            for i in range(long_matriz):
                if i == 0:
                    A[i, i] = 2 * (luces[cont_m] + luces[cont_m + 1])
                    if i + 1 <= long_matriz-1:
                        A[i, i + 1] = luces[cont_m + 1]
                    B[i, 0] = -6 * area[cont_m] * a[cont_m] / luces[cont_m] \
                            - 6 * area[cont_m + 1] * b[cont_m + 1] / luces[cont_m + 1]+mapoyos[cont_m]*luces[cont_m]+mapoyos[cont_m+2]*luces[cont_m + 1]
                else:
                    A[i, i] = 2 * (luces[cont_m + i] + luces[cont_m + i + 1])
                    A[i, i - 1] = luces[cont_m + i]
                    if i < long_matriz-1:
                        A[i, i + 1] = luces[cont_m + i + 1]
                    B[i, 0] = -6 * area[cont_m + i] * a[cont_m + i] / luces[cont_m + i] \
                            - 6 * area[cont_m + i + 1] * b[cont_m + i + 1] / luces[cont_m + i + 1]+mapoyos[cont_m + i + 2]*luces[cont_m + i + 1]

            # Resolver sistema A·x = B
            moments = np.linalg.solve(A, B).flatten()

            # Asignar momentos calculados a los apoyos
            for i in range(long_matriz):
                mapoyos[cont_m + i + 1] = -moments[i]

    return mapoyos
def Ley_M(luces, flag, carga_sup, carga_punt, x_punt, libre_izq, libre_der, Mapoyos):
    """
    Traducción de la función Ley_M de VBA a Python.
    
    Args:
        luces (list[float]): Longitudes de las luces.
        flag (int): Indicador de tipo de carga (1: superficial, 0: puntual).
        carga_sup (list[float]): Cargas superficiales.
        carga_punt (list[float]): Cargas puntuales.
        x_punt (list[float]): Posiciones relativas de las cargas puntuales.
        libre_izq (int): Indicador de libertad en el extremo izquierdo (1: libre, 0: no libre).
        libre_der (int): Indicador de libertad en el extremo derecho (1: libre, 0: no libre).
        Mapoyos (list[float]): Valores de los apoyos.

    Returns:
        list[float]: Distribución de momentos.
    """
    # Calcular la longitud total
    T_long = sum(luces)
    Ley = [0] * int(T_long * 100)  # Leyes por centímetro

    dist = 0
    if flag == 1:  # CARGA SUPERFICIAL
        for i in range(len(luces)):
            if i == 0 and libre_izq == 1:  # Voladizo izquierdo
                for j in range(int(dist), int(dist + luces[i] * 100)):
                    Ley[j] = carga_sup[i] * (j / 100) ** 2 / 2
            elif i == len(luces) - 1 and libre_der == 1:  # Voladizo derecho
                for j in range(int(dist), int(dist + luces[i] * 100)):
                    Ley[j] = Mapoyos[i] - (Mapoyos[i] - carga_sup[i] * (luces[i] - ((j - dist) / 100)) ** 2 / 2)
            else:
                Vi = carga_sup[i] * luces[i] / 2 + (Mapoyos[i] - Mapoyos[i + 1]) / luces[i]
                for j in range(int(dist), int(dist + luces[i] * 100)):
                    Ley[j] = Mapoyos[i] - Vi * (j - dist) / 100 + carga_sup[i] * ((j - dist) / 100) ** 2 / 2
            dist += luces[i] * 100

    else:  # CARGA PUNTUAL
        for i in range(len(luces)):
            if i == 0 and libre_izq == 1:  # Voladizo izquierdo
                for j in range(int(dist), int(dist + x_punt[i] * 100)):
                    Ley[j] = 0
                for j in range(int(dist + x_punt[i] * 100), int(dist + luces[i] * 100)):
                    Ley[j] = carga_punt[i] * (j / 100 - (dist / 100 + x_punt[i]))
            elif i == len(luces) - 1 and libre_der == 1:  # Voladizo derecho
                for j in range(int(dist), int(dist + x_punt[i] * 100)):
                    Ley[j] = Mapoyos[i] - carga_punt[i] * ((j - dist) / 100)
                for j in range(int(dist + x_punt[i] * 100), int(dist + luces[i] * 100)):
                    Ley[j] = 0
            else:
                Vi = carga_punt[i] * (luces[i] - x_punt[i]) / luces[i] + (Mapoyos[i] - Mapoyos[i + 1]) / luces[i]
                for j in range(int(dist), int(dist + x_punt[i] * 100)):
                    Ley[j] = Mapoyos[i] - Vi * (j - dist) / 100
                for j in range(int(dist + x_punt[i] * 100), int(dist + luces[i] * 100)):
                    Ley[j] = Mapoyos[i] - Vi * (j - dist) / 100 + carga_punt[i] * ((j - dist) / 100 - x_punt[i])
            dist += luces[i] * 100

    return Ley
def Ley_V(luces, flag, carga_sup, carga_punt, x_punt, libre_izq, libre_der, Mapoyos):
    """
    Traducción de la función Ley_V de VBA a Python.

    Args:
        luces (list[float]): Longitudes de las luces.
        flag (int): Indicador de tipo de carga (1: superficial, 0: puntual).
        carga_sup (list[float]): Cargas superficiales.
        carga_punt (list[float]): Cargas puntuales.
        x_punt (list[float]): Posiciones relativas de las cargas puntuales.
        libre_izq (int): Indicador de libertad en el extremo izquierdo (1: libre, 0: no libre).
        libre_der (int): Indicador de libertad en el extremo derecho (1: libre, 0: no libre).
        Mapoyos (list[float]): Valores de los apoyos.

    Returns:
        list[float]: Distribución de cortantes.
    """
    # Calcular la longitud total
    T_long = sum(luces)
    Ley = [0] * int(T_long * 100)  # Leyes por centímetro

    dist = 0
    if flag == 1:  # CARGA SUPERFICIAL
        for i in range(len(luces)):
            if i == 0 and libre_izq == 1:  # Voladizo izquierdo
                for j in range(int(dist), int(dist + luces[i] * 100)):
                    Ley[j] = carga_sup[i] * (j / 100)
            elif i == len(luces) - 1 and libre_der == 1:  # Voladizo derecho
                for j in range(int(dist), int(dist + luces[i] * 100)):
                    Ley[j] = -1 * (carga_sup[i] * luces[i] - (carga_sup[i] * luces[i] - carga_sup[i] * (luces[i] - (j - dist) / 100)))
            else:
                Vi = -1 * (carga_sup[i] * luces[i] / 2 + (Mapoyos[i] - Mapoyos[i + 1]) / luces[i])
                for j in range(int(dist), int(dist + luces[i] * 100)):
                    Ley[j] = Vi + carga_sup[i] * ((j - dist) / 100)
            dist += luces[i] * 100

    else:  # CARGA PUNTUAL
        for i in range(len(luces)):
            if i == 0 and libre_izq == 1:  # Voladizo izquierdo
                for j in range(int(dist), int(dist + x_punt[i] * 100)):
                    Ley[j] = 0
                for j in range(int(dist + x_punt[i] * 100), int(dist + luces[i] * 100)):
                    Ley[j] = carga_punt[i]
            elif i == len(luces) - 1 and libre_der == 1:  # Voladizo derecho
                for j in range(int(dist), int(dist + x_punt[i] * 100)):
                    Ley[j] = -1 * carga_punt[i]
                for j in range(int(dist + x_punt[i] * 100), int(dist + luces[i] * 100)):
                    Ley[j] = 0
            else:
                Vi = carga_punt[i] * (luces[i] - x_punt[i]) / luces[i] + (Mapoyos[i] - Mapoyos[i + 1]) / luces[i]
                for j in range(int(dist), int(dist + x_punt[i] * 100)):
                    Ley[j] = -Vi
                for j in range(int(dist + x_punt[i] * 100), int(dist + luces[i] * 100)):
                    Ley[j] = -Vi + carga_punt[i]
            dist += luces[i] * 100

    return Ley
def deformada(luces, leyM, EI, libreizq, libreder):
    """
    Calcula la deformada de una viga a partir de los momentos flectores.
    
    Parámetros:
    - luces: lista de longitudes de los vanos (en metros).
    - leyM: lista de momentos flectores en la viga (en kNm).
    - EI: rigidez de la viga (E * I) en kNm^2.
    - libreizq: 1 si el extremo izquierdo es libre, 0 si no.
    - libreder: 1 si el extremo derecho es libre, 0 si no.
    
    Retorno:
    - def: lista de deformaciones interpoladas (en cm).
    """
    # Calcular la longitud total de la viga
    totalluz = sum(luces)
    num_puntos = int(totalluz * 100)  # Divisiones por cm
    deflect = np.zeros(num_puntos+1)
    
    # Paso 1: Cálculo de la deformada en cada tramo
    cont = 0
    for i, luz in enumerate(luces):
        R = 0
        for j in range(cont, cont + int(luz * 100)):
            R += leyM[j] * (cont + int(luz * 100) - j)
        R = -R / (luz * 100)
        
        for j in range(cont, cont + int(luz * 100), 20):
            deflect[j] = R * (j - cont)
            for k in range(cont, j):
                deflect[j] -= (-leyM[k]) * (j - k)
        
        cont += int(luz * 100)

    # Paso 2: Interpolación de los valores de la deformada
    cont = 0
    for i, luz in enumerate(luces):
        for j in range(cont, cont + int(luz * 100), 20):
            max_local = min(j + 20, cont + int(luz * 100))
            for k in range(j+1, max_local):
                deflect[k] = deflect[j] + (deflect[max_local] - deflect[j]) * (k - j) / (max_local - j)
        cont += int(luz * 100)
    
    # Paso 3: Extremos libres
    if libreder == 1:
        giro = np.arctan(deflect[int((totalluz - luces[-1]) * 100)] - deflect[int((totalluz - luces[-1]) * 100) - 1])
        for i in range(int((totalluz - luces[-1]) * 100), int(totalluz * 100)):
            deflect[i] = 0
            for j in range(int((totalluz - luces[-1]) * 100), i):
                deflect[i] -= (-leyM[j]) * (i - j)
            deflect[i] = np.tan(giro) * (i - int((totalluz - luces[-1]) * 100)) + deflect[i]
    
    if libreizq == 1:
        giro = np.arctan(deflect[int(luces[0] * 100)] - deflect[int(luces[0] * 100) + 1])
        for i in range(int(luces[0] * 100), -1, -1):
            deflect[i] = 0
            for j in range(int(luces[0] * 100), i - 1, -1):
                deflect[i] -= (-leyM[j]) * (i - j)
            deflect[i] = -np.tan(giro) * (i - int(luces[0] * 100)) + deflect[i]
    
    # Paso 4: Ajuste final de la deformada
    deflect = -1000 * (deflect/10**4) / (EI)  # Salida en mm

    return deflect.tolist()


# if __name__ == '__main__':

#     luces = [5,5, 5]
#     flag = 1
#     carga_sup = [2,2,2]
#     carga_punt = []
#     x_punt = []
#     libre_izq = 0
#     libre_der = 0
#     cont = 1
#     EI = 16800  # Ejemplo: Rigidez (kNm^2)

#     MApoyos = clapeyron(luces, flag, carga_sup, carga_punt, x_punt, libre_izq, libre_der, cont)
#     print("Momentos en los apoyos:", MApoyos)
#     # Calcular las distribuciones
#     ley_momentos = Ley_M(luces, flag, carga_sup, carga_punt, x_punt, libre_izq, libre_der, MApoyos)
#     ley_cortantes = Ley_V(luces, flag, carga_sup, carga_punt, x_punt, libre_izq, libre_der, MApoyos)
#     deformacion = deformada(luces, ley_momentos, EI, libre_izq, libre_der)

#     # Crear el eje horizontal: distancia acumulada en metros
#     distancia = [i / 100 for i in range(len(ley_momentos))]
#     distanciaD = [i / 100 for i in range(len(deformacion))]

#     # Figura 1: Ley de Momentos
#     fig1 = plt.figure(figsize=(10, 6))
#     plt.plot(distancia, ley_momentos, label="Ley de Momentos", color="blue", linewidth=2)
#     plt.title("Distribución de Momentos a lo largo de la estructura", fontsize=14)
#     plt.xlabel("Distancia (m)", fontsize=12)
#     plt.ylabel("Momento (kNm)", fontsize=12)
#     plt.axhline(0, color="black", linewidth=0.8, linestyle="--")  # Línea de referencia en 0
#     plt.grid(True, linestyle="--", alpha=0.6)
#     plt.legend(fontsize=12)

#     # Figura 2: Ley de Cortantes
#     fig2 = plt.figure(figsize=(10, 6))
#     plt.plot(distancia, ley_cortantes, label="Ley de Cortantes", color="red", linewidth=2)
#     plt.title("Distribución de Cortantes a lo largo de la estructura", fontsize=14)
#     plt.xlabel("Distancia (m)", fontsize=12)
#     plt.ylabel("Cortante (kN)", fontsize=12)
#     plt.axhline(0, color="black", linewidth=0.8, linestyle="--")  # Línea de referencia en 0
#     plt.grid(True, linestyle="--", alpha=0.6)
#     plt.legend(fontsize=12)

#     # Figura 3: Deformada
#     fig3 = plt.figure(figsize=(10, 6))
#     plt.plot(distanciaD, deformacion, label="Deformada", color="black", linewidth=2)
#     plt.title("Deformada a lo largo de la estructura", fontsize=14)
#     plt.xlabel("Distancia (m)", fontsize=12)
#     plt.ylabel("Deformación (mm)", fontsize=12)
#     plt.axhline(0, color="black", linewidth=0.8, linestyle="--")  # Línea de referencia en 0
#     plt.grid(True, linestyle="--", alpha=0.6)
#     plt.legend(fontsize=12)

#     # Mostrar ambas figuras
#     plt.show()
