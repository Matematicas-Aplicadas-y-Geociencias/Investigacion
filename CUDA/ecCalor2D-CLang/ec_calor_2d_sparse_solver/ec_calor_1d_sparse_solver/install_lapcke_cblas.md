# Install LAPACK, LAPACKE y CBLAS

### Directivas para instalar LAPACK y BLAS tanto para FORTRAN como para C
```bash

sudo apt update 
sudo apt install liblapack-dev
sudo apt install liblapack-doc 
sudo apt install liblapack3 
sudo apt install libopenblas-base 
sudo apt install libopenblas-dev 
sudo apt install liblapacke
sudo apt install liblapacke-dev 
sudo apt install libblas-dev
sudo apt install libblas3

```

### Ejecutar programa con lapacke y cblas en C

```bash


$ gcc <source_name>.c -o <executable_name> -llapacke -lblas


```
