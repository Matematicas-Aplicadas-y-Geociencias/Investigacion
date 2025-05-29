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


def main() -> None:
    # Valores iniciales
    val1 = 16
    val2 = 32
    val3 = 32
    val4 = 32

    for i in range(6):
        if i != 0:
            ch = 2
        else:
            ch = 1

        var1 = "NUM_GANGSA"
        val1 = val1 * ch
        var2 = "VEC_LENGHTA"
        val2 = val2 // ch
        var3 = "NUM_GANGSB"
        val3 = val3 * ch
        var4 = "VEC_LENGHTB"
        val4 = val4 // ch

        cambiar_valor(var1, val1)
        cambiar_valor(var2, val2)
        cambiar_valor(var3, val3)
        cambiar_valor(var4, val4)

        compilar_programa()
        # tiempo_medido = medir_tiempo_ejecucion()

        tiempos_ejecucion = []
        for i in range(10):
            muestra = medir_tiempo_ejecucion()
            tiempos_ejecucion.append(muestra)

        tiempo_medido, desviacion_estandar = promedio_tiempos_ejecucion(
            tiempos_ejecucion
        )

        if os.path.exists("tiempos_medidos.csv"):
            df = pd.read_csv("tiempos_medidos.csv")
        else:
            df = pd.DataFrame(
                columns=["Parametros", "Tiempo Medido", "Desviacion Estandar"]
            )

        parametros = f"[{var1}={val1},{var2}={val2}],[{var3}={val3},{var4}={val4}]"
        df.loc[len(df)] = [parametros, tiempo_medido, desviacion_estandar]

        df.to_csv("tiempos_medidos.csv", index=False)


if __name__ == "__main__":
    main()
