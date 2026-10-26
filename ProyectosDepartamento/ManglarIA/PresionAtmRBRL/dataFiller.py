import marimo

__generated_with = "0.23.9"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    import scipy.io as sio
    import numpy as np
    import matlab.engine
    import pandas as pd
    import matplotlib.pyplot as plt
    import plotly.express as px
    import plotly.graph_objects as go
    from sklearn.linear_model import RANSACRegressor
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.impute import KNNImputer
    import filterpy
    from filterpy.kalman import KalmanFilter

    return KalmanFilter, RANSACRegressor, go, matlab, mo, np, pd, plt


@app.cell
def _(matlab):
    # Iniciar el motor de Matlab
    print("Iniciando MATLAB")
    me = matlab.engine.start_matlab()
    print("MATLAB SE HA INICIADO CORRECTAMENTE.")
    return (me,)


@app.cell
def _(me, np, pd):
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
    return data_Conagua, data_DiverlinkpEC09, data_EMpSF01, data_EMpSF01_raw


@app.cell
def _(data_Conagua, data_DiverlinkpEC09, data_EMpSF01_raw, go):
    # Graficar
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x = data_Conagua["time"],
            y = data_Conagua["PAtm"],
            name = "Conagua",
            mode = "markers"
        )
    )  

    fig.add_trace(
        go.Scatter(
            x = data_DiverlinkpEC09["time"],
            y = data_DiverlinkpEC09["PAtm"],
            name = "Diverlink.EC09",
            mode = "markers"
        )
    ) 

    fig.add_trace(
        go.Scatter(
            x = data_EMpSF01_raw["time"]["BigQuery"],
            y = data_EMpSF01_raw["PAtm"]["BigQuery"],
            name = "EM.SF01",
            mode = "markers"
        )
    )
    fig.show()
    return


@app.cell
def _(data_Conagua, data_DiverlinkpEC09, data_EMpSF01, pd):
    dfConagua = pd.DataFrame(data_Conagua).set_index("time")
    dfDLpEC09 = pd.DataFrame(data_DiverlinkpEC09).set_index("time")
    dfEMpSF01 = pd.DataFrame(data_EMpSF01).set_index("time")
    return dfConagua, dfDLpEC09, dfEMpSF01


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Elegir un sensor de referencia (Gold Standard)
    data_Conagua
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Encontrar las zonas de "solapamiento".
    """)
    return


@app.cell
def _(dfConagua, dfDLpEC09, dfEMpSF01):
    # 1. Resamplear a una frecuencia común (ej. cada 1 hora) para que coincidan los tiempos
    freq = "10min"
    ref = dfConagua["PAtm"].resample(freq).mean()
    div = dfDLpEC09["PAtm"].resample(freq).mean()
    em = dfEMpSF01["PAtm"].resample(freq).mean()
    return div, em, freq, ref


@app.cell
def _(div, em, pd, ref):
    # 2. Unirlos en un solo DataFrame para ver dónde se cruzan
    df_merged = pd.concat(
        [ref, div, em], axis=1, keys=["Ref", "Div", "EM"], sort=True
    )
    return (df_merged,)


@app.cell
def _(df_merged):
    # Ver dónde hay datos de ambos sensores al mismo tiempo
    overlap_ref_div = df_merged.dropna(subset=["Ref", "Div"])
    overlap_ref_em = df_merged.dropna(subset=["Ref", "EM"])
    return overlap_ref_div, overlap_ref_em


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Usar ML para calibrar / alinear los sensores
    """)
    return


@app.cell
def _(overlap_ref_div):
    # Entrenar el modelo para convertir Diverlink -> Referencia (Conagua)
    # X es el sensor a calibrar, y es el sensor de referencia
    X_div = overlap_ref_div[["Div"]].rename(columns={"Div": "PAtm"}).dropna()
    y_div = overlap_ref_div["Ref"].loc[X_div.index]
    return X_div, y_div


@app.cell
def _(RANSACRegressor, X_div, y_div):
    # Modelo ML: RANSAC (Ignora lecturas erráticas de los sensores)
    ransac_div = RANSACRegressor(random_state=42)
    ransac_div.fit(X_div, y_div)
    return (ransac_div,)


@app.cell
def _(RANSACRegressor, overlap_ref_em):
    # Hacer lo mismo para EM
    X_em = overlap_ref_em[["EM"]].rename(columns={"EM": "PAtm"}).dropna()
    y_em = overlap_ref_em["Ref"].loc[X_em.index]

    ransac_em = RANSACRegressor(random_state=42)
    ransac_em.fit(X_em, y_em)
    return (ransac_em,)


@app.cell
def _(dfDLpEC09, dfEMpSF01, ransac_div, ransac_em):
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
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Unificar y rellenar huecos (imputación con ML)
    """)
    return


@app.cell
def _(dfConagua, dfDLpEC09, dfEMpSF01, freq, pd):
    dfConaguaSD = dfConagua.groupby(dfConagua.index).mean()
    dfDLpEC09SD = dfDLpEC09.groupby(dfDLpEC09.index).mean()
    dfEMpSF01SD = dfEMpSF01.groupby(dfEMpSF01.index).mean()
    dfSerieTiempoCompleta = (
        pd.concat(
            [
                dfConaguaSD[["PAtm"]].rename(columns={"PAtm": "Ref"}),
                dfDLpEC09SD[["PAtm_Corregida"]].rename(
                    columns={"PAtm_Corregida": "Div"}
                ),
                dfEMpSF01SD[["PAtm_Corregida"]].rename(
                    columns={"PAtm_Corregida": "EM"}
                ),
            ],
            axis=1,
            sort=True,
        )
        .resample(freq)
        .mean()
    )
    return (dfSerieTiempoCompleta,)


@app.cell
def _(dfSerieTiempoCompleta):
    # 2. Crear la serie unificada (promedio de los sensores disponibles en cada hora)
    dfSerieTiempoCompleta["PAtm_Unificada"] = dfSerieTiempoCompleta.mean(axis=1)
    return


@app.cell
def _(dfSerieTiempoCompleta):
    # 3. Rellenar huecos (Imputación)
    # Opción A: KNN Imputer (ML clásico) - Busca horas similares en el pasado para rellenar
    # imputer = KNNImputer(n_neighbors=144) # Busca las 24 horas más similares
    # dfSerieTiempoCompleta['PAtm_Final'] = imputer.fit_transform(dfSerieTiempoCompleta[['PAtm_Unificada']])

    # Opción B: Interpolación temporal (Más rápido y suele funcionar mejor para presión)
    # dfSerieTiempoCompleta['PAtm_Final'] = dfSerieTiempoCompleta['PAtm_Unificada'].interpolate(method='time')
    dfSerieTiempoCompleta["PAtm_Final"] = dfSerieTiempoCompleta[
        "PAtm_Unificada"
    ].interpolate(method="pchip")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Suavizar la transición
    """)
    return


@app.cell
def _(KalmanFilter, dfSerieTiempoCompleta, np):
    # Configuración básica de un Filtro de Kalman para suavizar
    kf = KalmanFilter(dim_x=2, dim_z=1)
    kf.x = np.array(
        [dfSerieTiempoCompleta["PAtm_Final"].iloc[0], 0]
    )  # Estado inicial [presión, velocidad]
    kf.F = np.array([[1, 1], [0, 1]])  # Matriz de transición
    kf.H = np.array([[1, 0]])  # Matriz de observación
    kf.P *= 100  # Incertidumbre inicial
    kf.R = 5  # Ruido del sensor
    kf.Q = 0.1  # Ruido del proceso

    suavizado = []
    for z in dfSerieTiempoCompleta["PAtm_Final"]:
        if not np.isnan(z):
            kf.predict()
            kf.update(z)
            suavizado.append(kf.x[0])
        else:
            kf.predict()
            suavizado.append(kf.x[0])

    dfSerieTiempoCompleta["PAtm_Kalman"] = suavizado
    return


@app.cell
def _(dfSerieTiempoCompleta):
    dfSerieTiempoCompleta
    return


@app.cell
def _(dfEMpSF01, dfSerieTiempoCompleta, plt):
    plt.plot(dfSerieTiempoCompleta["PAtm_Kalman"], "ko", markersize=0.1)
    # plt.plot(dfSerieTiempoCompleta['PAtm_Final'], 'o', markersize=0.1)
    # plt.plot(dfConagua['PAtm'], 'o', markersize=0.1)
    plt.plot(dfEMpSF01["PAtm"], "o", markersize=0.1)
    # plt.plot(dfDLpEC09['PAtm'], 'o', markersize=0.1)
    return


if __name__ == "__main__":
    app.run()
