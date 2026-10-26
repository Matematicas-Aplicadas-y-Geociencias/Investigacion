import logging
import numpy as np
import matlab.engine
import pandas as pd
import plotly.graph_objects as go
from sklearn.linear_model import RANSACRegressor
# from sklearn.ensemble import RandomForestRegressor
# from sklearn.impute import KNNImputer
# import filterpy
# from filterpy.kalman import KalmanFilter

# Configurar logging con formato estándar de timestamp, nivel y mensaje
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Iniciar el motor de Matlab
logger.info("Iniciando motor de MATLAB")
me = matlab.engine.start_matlab()
logger.info("MOTOR DE MATLAB SE HA INICIADO CORRECTAMENTE.")

# Cargar el archivo en el espacio de trabajo de matlab
me.eval("load('Datos_Presion.mat')", nargout=0)

# Cargar datos Conagua
# Presion atmosferica
matlab_pAtm = me.eval("data_Conagua.PAtm", nargout=1)
py_PAtm = np.array(matlab_pAtm).flatten()

# Tiempo
matlab_time = me.eval("posixtime(data_Conagua.time);", nargout=1)
py_time = pd.to_datetime(np.array(matlab_time).flatten(), unit="s")

data_Conagua = {"PAtm": py_PAtm, "time": py_time}

# Cargar datos Diverlink
# Presion atmosferica
matlab_pAtm = me.eval("data_Diverlink.EC09.PAtm", nargout=1)
py_PAtm = np.array(matlab_pAtm).flatten()

# Tiempo
matlab_time = me.eval("posixtime(data_Diverlink.EC09.time);", nargout=1)
py_time = pd.to_datetime(np.array(matlab_time).flatten(), unit="s")

data_DiverlinkpEC09 = {"PAtm": py_PAtm, "time": py_time}

# Cargar datos EM
# Presion atmosferica
matlab_pAtm_local_1 = me.eval("data_EM.SF01.PAtm.local_1", nargout=1)
py_PAtm_local_1 = np.array(matlab_pAtm_local_1).flatten()

matlab_pAtm_local_2 = me.eval("data_EM.SF01.PAtm.local_2", nargout=1)
py_PAtm_local_2 = np.array(matlab_pAtm_local_2).flatten()

matlab_pAtm_bq = me.eval("data_EM.SF01.PAtm.BigQuery", nargout=1)
py_PAtm_bq = np.array(matlab_pAtm_bq).flatten()

# Tiempo
matlab_time_l1 = me.eval("posixtime(data_EM.SF01.time.local_1);", nargout=1)
py_time_l1 = pd.to_datetime(np.array(matlab_time_l1).flatten(), unit="s")
matlab_time_l2 = me.eval("posixtime(data_EM.SF01.time.local_2);", nargout=1)
py_time_l2 = pd.to_datetime(np.array(matlab_time_l2).flatten(), unit="s")
matlab_time_bq = me.eval("posixtime(data_EM.SF01.time.BigQuery);", nargout=1)
py_time_bq = pd.to_datetime(np.array(matlab_time_bq).flatten(), unit="s")

data_EMpSF01_raw = {
    "PAtm": {
        "local_1": py_PAtm_local_1,
        "local_2": py_PAtm_local_2,
        "BigQuery": py_PAtm_bq,
    },
    "time": {
        "local_1": py_time_l1,
        "local_2": py_time_l2,
        "BigQuery": py_time_bq,
    },
}
data_EMpSF01 = {
    "PAtm": py_PAtm_bq,
    "time": py_time_bq,
}

me.quit()
logger.info("Finalizando motor de MATLAB")

dfConagua = pd.DataFrame(data_Conagua).set_index("time")
dfDLpEC09 = pd.DataFrame(data_DiverlinkpEC09).set_index("time")
dfEMpSF01 = pd.DataFrame(data_EMpSF01).set_index("time")

logger.info(f"Completando serie de tiempo de Conagua")
# Cargar datos RL
dfRL = pd.read_csv(
    "Rio_Lagartos.csv",
    header=None,
    skiprows=10,  # Saltar las 9 líneas de metadatos
    nrows=11676,
    names=[
        "Fecha_Local",
        "Fecha_UTC",
        "Dir_V1",
        "Dir_V2",
        "Vel_V1",
        "Vel_V2",
        "Temperatura",
        "Humedad",
        "PAtm",
        "Lluvia",
        "Radiacion",
    ],
)

# Convertir fechas de str a datetime
dfRL["Fecha_Local"] = pd.to_datetime(dfRL["Fecha_Local"], format="%d/%m/%Y %H:%M")

# Establecer la fecha como índice (usando Fecha_Local o Fecha_UTC según tu preferencia)
dfRL = dfRL.set_index("Fecha_Local").sort_index()

# Concatenar datos faltantes a dfConagua
dfConagua = pd.concat([dfConagua[["PAtm"]], dfRL[["PAtm"]]], axis=0)

logger.info(f"Serie de Tiempo Conagua Completa: ")
print(dfConagua.info())

logger.info("Graficando datos crudos.")
# Graficar
fig = go.Figure()
fig.add_trace(
    go.Scatter(
        x=dfConagua.index,
        y=dfConagua["PAtm"],
        name="Conagua",
        mode="markers",
        opacity=0.7,
        marker=dict(size=0.9, color="black"),
    )
)
fig.add_trace(
    go.Scatter(
        x=dfDLpEC09.index,
        y=dfDLpEC09["PAtm"],
        name="Diverlink.EC09",
        mode="markers",
        opacity=0.7,
        marker=dict(size=0.9, color="red"),
    )
)
fig.add_trace(
    go.Scatter(
        x=dfEMpSF01.index,
        y=dfEMpSF01["PAtm"],
        name="EM.SF01",
        mode="markers",
        opacity=0.7,
        marker=dict(size=0.9, color="blue"),
    )
)
fig.show()

"""
## Elegir un sensor de referencia (Gold Standard)
data_Conagua
"""
logger.info("Señal de referencia: Conagua")


"""
## Encontrar las zonas de "solapamiento".
"""

# 1. Resamplear a una frecuencia común (ej. cada 1 hora) para que coincidan los tiempos
freq = "10min"
ref = dfConagua["PAtm"].resample(freq).mean()
div = dfDLpEC09["PAtm"].resample(freq).mean()
em = dfEMpSF01["PAtm"].resample(freq).mean()


# 2. Unirlos en un solo DataFrame para ver dónde se cruzan
df_merged = pd.concat([ref, div, em], axis=1, keys=["Ref", "Div", "EM"], sort=True)

# Ver dónde hay datos de ambos sensores al mismo tiempo
overlap_ref_div = df_merged.dropna(subset=["Ref", "Div"])
overlap_ref_em = df_merged.dropna(subset=["Ref", "EM"])


"""
## Usar ML para calibrar / alinear los sensores
"""
logger.info("Usando ML para calibrar/alinear las señales.")
# Entrenar el modelo para convertir Diverlink -> Referencia (Conagua)
# X es el sensor a calibrar, y es el sensor de referencia
X_div = overlap_ref_div[["Div"]].rename(columns={"Div": "PAtm"}).dropna()
y_div = overlap_ref_div["Ref"].loc[X_div.index]


# Modelo ML: RANSAC (Ignora lecturas erráticas de los sensores)
ransac_div = RANSACRegressor(random_state=42)
ransac_div.fit(X_div, y_div)


# Hacer lo mismo para EM
X_em = overlap_ref_em[["EM"]].rename(columns={"EM": "PAtm"}).dropna()
y_em = overlap_ref_em["Ref"].loc[X_em.index]

ransac_em = RANSACRegressor(random_state=42)
ransac_em.fit(X_em, y_em)


# --- Aplicar la corrección a TODA la serie de tiempo ---
# Ahora transformamos los datos originales de Diverlink y EM
mask_valid_div = dfDLpEC09["PAtm"].notna()
dfDLpEC09.loc[mask_valid_div, "PAtm_Corregida"] = ransac_div.predict(
    dfDLpEC09.loc[mask_valid_div, ["PAtm"]]
)
mask_valid_em = dfEMpSF01["PAtm"].notna()
dfEMpSF01.loc[mask_valid_em, "PAtm_Corregida"] = ransac_em.predict(
    dfEMpSF01.loc[mask_valid_em, ["PAtm"]]
)


"""
## Unificar y rellenar huecos (imputación con ML)
"""
logger.info("Unificando y rellenando huecos.")

dfConaguaSD = dfConagua.groupby(dfConagua.index).mean()
dfDLpEC09SD = dfDLpEC09.groupby(dfDLpEC09.index).mean()
dfEMpSF01SD = dfEMpSF01.groupby(dfEMpSF01.index).mean()
dfSerieTiempoCompleta = (
    pd.concat(
        [
            dfConaguaSD[["PAtm"]].rename(columns={"PAtm": "Ref"}),
            dfDLpEC09SD[["PAtm_Corregida"]].rename(columns={"PAtm_Corregida": "Div"}),
            dfEMpSF01SD[["PAtm_Corregida"]].rename(columns={"PAtm_Corregida": "EM"}),
        ],
        axis=1,
        sort=True,
    )
    .resample(freq)
    .mean()
)


# 2. Crear la serie unificada (promedio de los sensores disponibles en cada hora)
dfSerieTiempoCompleta["PAtm_Unificada"] = dfSerieTiempoCompleta.mean(axis=1)


# 3. Rellenar huecos (Imputación)
# Opción A: KNN Imputer (ML clásico) - Busca horas similares en el pasado para rellenar
# imputer = KNNImputer(n_neighbors=144) # Busca las 24 horas más similares
# dfSerieTiempoCompleta['PAtm_Final'] = imputer.fit_transform(dfSerieTiempoCompleta[['PAtm_Unificada']])

# Opción B: Interpolación temporal (Más rápido y suele funcionar mejor para presión)
# dfSerieTiempoCompleta['PAtm_Final'] = dfSerieTiempoCompleta['PAtm_Unificada'].interpolate(method='time')
dfSerieTiempoCompleta["PAtm_Final"] = dfSerieTiempoCompleta[
    "PAtm_Unificada"
].interpolate(method="pchip")


"""
## Suavizar la transición
"""

# # Configuración básica de un Filtro de Kalman para suavizar
# kf = KalmanFilter(dim_x=2, dim_z=1)
# kf.x = np.array(
#     [dfSerieTiempoCompleta["PAtm_Final"].iloc[0], 0]
# )  # Estado inicial [presión, velocidad]
# kf.F = np.array([[1, 1], [0, 1]])  # Matriz de transición
# kf.H = np.array([[1, 0]])  # Matriz de observación
# kf.P *= 100  # Incertidumbre inicial
# kf.R = 5  # Ruido del sensor
# kf.Q = 0.1  # Ruido del proceso
#
# suavizado = []
# for z in dfSerieTiempoCompleta["PAtm_Final"]:
#     if not np.isnan(z):
#         kf.predict()
#         kf.update(z)
#         suavizado.append(kf.x[0])
#     else:
#         kf.predict()
#         suavizado.append(kf.x[0])
#
# dfSerieTiempoCompleta["PAtm_Kalman"] = suavizado


logger.info(f"Serie de Tiempo Completa: ")
print(dfSerieTiempoCompleta.info())
logger.info("Graficando serie completa.")

fig = go.Figure()
# fig.add_trace(
#     go.Scatter(
#         x=dfSerieTiempoCompleta.index,
#         y=dfSerieTiempoCompleta["PAtm_Kalman"],
#         mode="markers",
#         name="PAtm Kalman",
#         opacity=0.7,
#         marker=dict(size=2),
#     )
# )
fig.add_trace(
    go.Scatter(
        x=dfSerieTiempoCompleta.index,
        y=dfSerieTiempoCompleta["PAtm_Final"],
        mode="markers",
        name="PAtm Final",
        opacity=0.7,
        marker=dict(size=1, color="black"),
    )
)
fig.add_trace(
    go.Scatter(
        x=dfConagua.index,
        y=dfConagua["PAtm"],
        mode="markers",
        name="Conagua",
        opacity=0.7,
        marker=dict(size=0.9, color="red"),
    )
)
fig.add_trace(
    go.Scatter(
        x=dfEMpSF01.index,
        y=dfEMpSF01["PAtm"],
        mode="markers",
        name="EM.SF01",
        opacity=0.7,
        marker=dict(size=0.9, color="blue"),
    )
)
fig.add_trace(
    go.Scatter(
        x=dfDLpEC09.index,
        y=dfDLpEC09["PAtm"],
        mode="markers",
        name="Diverlink.EC09",
        opacity=0.7,
        marker=dict(size=0.9, color="green"),
    )
)
fig.show()

logger.info("Guardando datos para usar en Matlab.")
dfSerieTiempoCompleta.index.name = "Time"

dfSerieTiempoCompleta[["PAtm_Final"]].to_csv(
    "PresionAtmosfericaCompleta.csv", index=True
)

logger.info("✅ Archivo guardado con datetime nativo de MATLAB")
