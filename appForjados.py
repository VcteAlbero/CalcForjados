import streamlit as st
import calculosForjados
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Título de la aplicación
st.title("Cálculo de Forjados")

# Entrada de datos en sidebar
st.sidebar.header("Datos de Entrada")
nvanos = st.sidebar.number_input("Nº de vanos:", min_value=1, value=3)
#
checkbox_option1 = st.sidebar.checkbox("vol_izq")
libre_izq = 1 if checkbox_option1 else 0
checkbox_option2 = st.sidebar.checkbox("vol_der")
libre_der = 1 if checkbox_option2 else 0
#
radio_option = st.sidebar.radio("Tipo de carga:", ["SUP", "PUNT"])
flagCarga = 1 if radio_option == "SUP" else 0
radio_option = st.sidebar.radio("Tipo de cálculo:", ["CONT", "ISO"])
flagCONT = 1 if radio_option == "CONT" else 0
#
EI = st.sidebar.number_input("RIGIDEZ (kNm2):", min_value=1, value=16800)

# Generar una tabla editable en función de n vanos
st.write("### Introduce los datos en la tabla:")
if flagCarga ==1:
    data = pd.DataFrame({
        "Luz (m)": [5] * nvanos,
        "Carga (kN/m2)": [2]*nvanos,
    })
else:
        data = pd.DataFrame({
        "Luz (m)": [5] * nvanos,
        "Carga (kN)": [1]*nvanos,
        "x_dist (m)": [2.5]*nvanos,
    })

edited_data = st.data_editor(data, num_rows="dynamic")

#Almacenar datos
luces = np.array(edited_data["Luz (m)"])
if flagCarga ==1:
     carga_sup = np.array(edited_data["Carga (kN/m2)"])
     carga_punt = []
     x_punt = []
else:
     carga_sup = []
     carga_punt = np.array(edited_data["Carga (kN)"])
     x_punt = np.array(edited_data["x_dist (m)"])

#Calculos
MApoyos = calculosForjados.clapeyron(luces, flagCarga, carga_sup, carga_punt, x_punt, libre_izq, libre_der, flagCONT)
ley_momentos = calculosForjados.Ley_M(luces, flagCarga, carga_sup, carga_punt, x_punt, libre_izq, libre_der, MApoyos)
ley_cortantes = calculosForjados.Ley_V(luces, flagCarga, carga_sup, carga_punt, x_punt, libre_izq, libre_der, MApoyos)
deformacion = calculosForjados.deformada(luces, ley_momentos, EI, libre_izq, libre_der)

# Ley de momento flectores
st.write("### 1. Ley de momentos flectores:")
distancia = [i / 100 for i in range(len(ley_momentos))]
fig1 = plt.figure(figsize=(10, 3))
plt.plot(distancia, ley_momentos, label="Ley de Momentos", color="blue", linewidth=2)
plt.title("Distribución de Momentos a lo largo de la estructura", fontsize=14)
plt.xlabel("Distancia (m)", fontsize=12)
plt.ylabel("Momento (kNm)", fontsize=12)
plt.axhline(0, color="black", linewidth=0.8, linestyle="--")  # Línea de referencia en 0
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend(fontsize=12)
st.pyplot(fig1)

# Ley de cortantes
st.write("### 2. Ley de cortantes:")
fig2 = plt.figure(figsize=(10, 3))
plt.plot(distancia, ley_cortantes, label="Ley de Cortantes", color="red", linewidth=2)
plt.title("Distribución de Cortantes a lo largo de la estructura", fontsize=14)
plt.xlabel("Distancia (m)", fontsize=12)
plt.ylabel("Cortante (kN)", fontsize=12)
plt.axhline(0, color="black", linewidth=0.8, linestyle="--")  # Línea de referencia en 0
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend(fontsize=12)
st.pyplot(fig2)

# Deformada
st.write("### 3. Deformada:")
distanciaD = [i / 100 for i in range(len(deformacion))]
fig3 = plt.figure(figsize=(10, 1.5))
plt.plot(distanciaD, deformacion, label="Deformada", color="black", linewidth=2)
plt.title("Deformada a lo largo de la estructura", fontsize=14)
plt.xlabel("Distancia (m)", fontsize=12)
plt.ylabel("Deformación (mm)", fontsize=12)
plt.axhline(0, color="black", linewidth=0.8, linestyle="--")  # Línea de referencia en 0
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend(fontsize=12)
st.pyplot(fig3)

# Pie de página con información de contacto
st.markdown(
    """
    <hr style="border:1px solid gray;margin-top:50px;">
    <footer>
        <p style="text-align:center; font-size:12px; color:gray;">
        Esta es una aplicación en pruebas. Si tienes comentarios o preguntas, por favor contacta con: 
        <a href="mailto:valbero@uji.es">valbero@uji.es</a>.
        </p>
    </footer>
    """,
    unsafe_allow_html=True
)