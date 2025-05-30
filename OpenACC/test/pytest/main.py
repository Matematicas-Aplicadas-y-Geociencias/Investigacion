import subprocess
import time
import re
import os
import pandas as pd


def cambiar_valor(variable: str, valor: int) -> bool:
    with open("../source/constantes.f90", "r") as archivo:
        contenido = archivo.read()

    patron = rf"({re.escape(variable)}\s*=\s*)\d+"

    def remplazo(match) -> str:
        return f"{match.group(1)}{valor}"

    contenido_modificado, numero_sustituciones = re.subn(patron, remplazo, contenido)

    if numero_sustituciones == 0:
        print(f"No se encontró el parámetro {variable} en el archivo.")
        return False

    with open("../source/constantes.f90", "w") as archivo:
        archivo.write(contenido_modificado)

    print(f"Parámetro {variable} modificado a {valor} exitosamente.")
    return True


def compilar_programa() -> None:
    print()
    print("=====Compilando el Programa=====")
    resultado = subprocess.run(
        ["make -C ../build"], shell=True, capture_output=True, text=True
    )
    print(f"{resultado.stdout}")
    print(f"{resultado.stderr}")


def medir_tiempo_ejecucion() -> float:
    start = time.perf_counter()
    resultado = subprocess.run(
        ["./../build/CALOR2DACC_ARR1D-gpu"], shell=True, capture_output=True, text=True
    )
    end = time.perf_counter()

    tiempo_ejecucion = end - start

    print("=====Tiempo de ejecución=====")
    print(f"Tiempo: {tiempo_ejecucion:.4f} segundos")

    return tiempo_ejecucion


def promedio_tiempos_ejecucion(tiempos_ejecucion: list[float]) -> tuple[float, float]:
    serie = pd.Series(tiempos_ejecucion)
    promedio = serie.mean()
    desviacion_estandar = serie.std()
    return promedio, desviacion_estandar


def datos_medicion(
    tiempos_ejecucion: list[float], parametros: str, numero_puntos: str
) -> None:
    tiempo_medido, desviacion_estandar = promedio_tiempos_ejecucion(tiempos_ejecucion)
    nombre_archivo = f"tiempos_medidos_{numero_puntos}.csv"

    if os.path.exists(nombre_archivo):
        df = pd.read_csv(nombre_archivo)
    else:
        df = pd.DataFrame(
            columns=["Parametros", "Tiempo Medido", "Desviacion Estandar"]
        )

    df.loc[len(df)] = [parametros, tiempo_medido, desviacion_estandar]

    df.to_csv(nombre_archivo, index=False)


def main() -> None:
    # Valores iniciales
    NUM_GANGSA = "NUM_GANGSA"
    VEC_LENGHTA = "VEC_LENGHTA"
    NUM_GANGSB = "NUM_GANGSB"
    VEC_LENGHTB = "VEC_LENGHTB"
    MI = "mi"
    NJ = "nj"

    val1 = 128
    val2 = 64

    for i in range(6):
        ch = 2 if i != 0 else 1

        val1 = val1 * ch
        val2 = val2 * ch

        cambiar_valor(MI, val1)
        cambiar_valor(NJ, val2)

        val3 = 16
        val4 = 32
        val5 = 32
        val6 = 32

        for j in range(5):
            ch = 2 if j != 0 else 1

            val3 = val3 * ch
            val4 = val4 // ch
            val5 = val5 * ch
            val6 = val6 // ch

            cambiar_valor(NUM_GANGSA, val3)
            cambiar_valor(VEC_LENGHTA, val4)
            cambiar_valor(NUM_GANGSB, val5)
            cambiar_valor(VEC_LENGHTB, val6)

            compilar_programa()

            tiempos_ejecucion = []
            for i in range(10):
                muestra = medir_tiempo_ejecucion()
                tiempos_ejecucion.append(muestra)

            parametros = f"[{NUM_GANGSA}={val3},{VEC_LENGHTA}={val4}],[{NUM_GANGSB}={val5},{VEC_LENGHTB}={val6}]"
            numero_puntos = f"{val1}x{val2}"
            datos_medicion(tiempos_ejecucion, parametros, numero_puntos)


if __name__ == "__main__":
    main()
