# README

Se usa **test** como nombre del ejecutable para que -git- no rastree el ejecutable.

### Directivas para compilar usando la libreria lapacke

```bash

$ gcc eq_heat_1d_lapacke_solver.c -o test -llapacke

```

### Directivas para compilar usando la libreria lapacke

```bash

$ nvc -acc=gpu -gpu=ccXX -o test eq_heat_1d_cusolver.c -lcusparse -lcusolver

```

### Otras directivas por si son necesarias

```bash

nvc -acc=gpu -gpu=ccXX -o <name_executable> <source_name> -lcudart -lcublas -lcusparse -lcusolver

```
