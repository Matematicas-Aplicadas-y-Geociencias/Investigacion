Sugiero revisar el parámetro ```Default Target``` de la GPU con el siguiente comando:
```bash
$ nvaccelinfo
```
En mi caso, el ```Defaul Target``` para mi tarjeta gráfica es
```bash
Default Target:		cc75
```

La directiva para compilar **OpenAcc** en lengua C es:
```bash
$ nvc -fast -acc=gpu -gpu=ccXX -Minfo=all <file_name>.c -o <executable_name>
```

La directiva para compilar **OpenAcc** en lengua Fortran es:
```bash
$ nvfortran -fast -acc=gpu -gpu=ccXX -Minfo=all <file_name>.f90 -o <executable_name>
```
