A continuación se dan comandos para compilar y perfilar adecuadamente programas desarrollados con OpenMP

compilaciòn:
```bash
nvfortran -Minfo=all -mp=ompt *.f90 -o <nombre_ejecutable>
```

ejecución con perfilador
```bash
nsys profile -o report -t osrt,openmp --stats=true -b fp  --force-overwrite true <nombre_ejecutable>
```



