import subprocess
import time
import re
import os
import pandas as pd

def cambiar_valor(variable: str, valor: int) -> bool:
    with open('../source/constantes.f90', 'r') as archivo:
        contenido = archivo.read()

    patron: str = rf"({re.escape(variable)}\s*=\s*)\d+"

    def remplazo(match) -> str:
        return f"{match.group(1)}{valor}"
    
    contenido_modificado, numero_sustituciones = re.subn(patron, remplazo, contenido)

    if numero_sustituciones == 0:
        print(f"No se encontró el parámetro {variable} en el archivo.")
        return False

    with open("../source/constantes.f90", 'w') as archivo:
        archivo.write(contenido_modificado)

    print(f"Parámetro {variable} modificado a {valor} exitosamente.")
    return True

def compilar_programa() -> None:
    print()
    print(f"=====Compilando el Programa=====")
    resultado = subprocess.run(
            ["make -C ../build"],
            shell=True, capture_output=True,
            text=True)
    print(f"{resultado.stdout}")
    print(f"{resultado.stderr}")

def medir_tiempo_ejecucion() -> float:
    start = time.perf_counter()
    resultado = subprocess.run(
            ["./../build/CALOR2DACC_ARR1D-gpu"],
            shell=True, capture_output=True,
            text=True)
    end = time.perf_counter()

    tiempo_ejecucion = end - start

    print(f"=====Tiempo de ejecución=====")
    print(f"Tiempo: {tiempo_ejecucion:.4f} segundos")

    return tiempo_ejecucion

def main() -> None:
    val1=16
    val2=32
    val3=32
    val4=32

    for i in range(6):
        if i == 0:
            ch = 1
        else:
            ch = 2

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
        tiempo_medido = medir_tiempo_ejecucion()
        
        if os.path.exists("tiempos_medidos.csv"):
            df = pd.read_csv("tiempos_medidos.csv")
        else:
            df = pd.DataFrame(columns=["Parametros", "Tiempo Medido"])

        parametros = f"[{var1}={val1},{var2}={val2}],[{var3}={val3},{var4}={val4}]"
        df.loc[len(df)] = [parametros, tiempo_medido]

        df.to_csv("tiempos_medidos.csv", index=False)

if __name__ == "__main__":
    main()
