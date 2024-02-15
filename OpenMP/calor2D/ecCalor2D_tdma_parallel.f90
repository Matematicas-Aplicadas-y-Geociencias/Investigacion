PROGRAM Calor2D
!
use constantes
use omp_lib
implicit none !(type, external)
! external :: sgesv
integer :: ii, jj, kk
!
! Variables de tamaño para la placa
double precision :: a, b, deltax, deltay
!
! Parámetros físicos del problema
double precision :: cond_ter, temp_ini, temp_fin
double precision :: flux_aba, flux_arr
!
! Arreglos para incógnitas, término fuente y posición
double precision, dimension(mi)    :: sourcex, xx
double precision, dimension(nj)    :: sourcey, yy
double precision, dimension(mi)    :: resultx,tempx
double precision, dimension(nj)    :: resulty,tempy
double precision, dimension(mi,nj) :: temper,temp_ant
!
! Matriz a invertir
double precision, dimension(mi)    :: AI,AD
double precision, dimension(mi)    :: AC
double precision, dimension(nj)    :: BI,BD
double precision, dimension(nj)    :: BC
!
! Se crea la malla 2D
!
a = 10.d0
b = 5.d0
do ii = 1, mi
    xx(ii) = a*dfloat(ii-1)/dfloat(mi-1)
!     print*, ii, xx(ii)
end do
do jj = 1, nj
    yy(jj) = b*dfloat(jj-1)/dfloat(nj-1)
end do
deltax = 1.d0/dfloat(mi-1)
deltay = 1.d0/dfloat(nj-1)
!
! Se definen los parámetros físicos del problema
!
cond_ter = 97.5d-6
temp_ini = 308.d0
temp_fin = 298.d0
flux_aba = 0.d0
flux_arr = 0.d0
!
! Inicialización de los arreglos a utilizar
!
AI = 0.d0; AC = 0.d0; AD = 0.d0
BI = 0.d0; BC = 0.d0; BD = 0.d0
temper = (temp_fin+temp_ini)/2.d0
resultx = 0.0d0
resulty = 0.0d0
!
!
! Se muestra el máximo número permitido de hilos
    !max_proc = omp_get_max_threads()
    !print*, "El máximo número de hilos permitido es: ", max_proc
!
! Se obtiene el identificador de hilo
    !id = omp_get_thread_num()
    !print*, " Yo soy el hilo: ", id
!
!! Bucle iterativo
!
do kk = 1, 50
    !
    !
    !!! Abrimos la región paralela
    !
    !
    !$OMP PARALLEL DEFAULT(NONE)&
    !$OMP& SHARED(cond_ter, temp_ini, temp_fin, deltax, deltay)&
    !$OMP& SHARED(temper, flux_aba, flux_arr, kk)&
    !$OMP& PRIVATE(AI, AC, AD, BI, BC, BD, resultx, resulty)&
    !$OMP& PRIVATE(tempx, tempy)
    !
    !
    !!! PARTE 1
    !
    !
    ! Inicia el ciclo que recorre la coordenada y resolviendo
    ! problemas 1D en la dirección de x
    !
    !$OMP DO
    do jj = 2, nj-1
        !
        ! Se definen las condiciones de frontera
        !
        AC(1)      = 1.0d0
        AD(1)      = 0.0d0
        resultx(1)  = temp_ini
        AI(mi)     = 0.d0
        AC(mi)     = 1.0d0
        resultx(mi) = temp_fin
        !
        ! Ensamblado de la matriz tridiagonal
        ! y del vector de resultados
        !
        do ii=2, mi-1
            AI(ii)     =-1.0d0*cond_ter/(deltax*deltax)
            AC(ii)     = 2.0d0*cond_ter*(1.d0/(deltax*deltax)+1.d0/(deltay*deltay))
            AD(ii)     =-1.0d0*cond_ter/(deltax*deltax)
            resultx(ii) = cond_ter/(deltay*deltay)*temper(ii,jj+1)+cond_ter/(deltay*deltay)*temper(ii,jj-1)
        end do
        !
        ! Llamamos al TDMA
        call tri(AI,AC,AD,resultx,tempx,mi)
        do ii =1, mi
            temper(ii,jj) = tempx(ii)
        end do
    end do
    !$OMP END DO
    !
    !
    !!! PARTE 2
    !
    !
    ! Inicia el ciclo que recorre la coordenada x resolviendo
    ! problemas 1D en la dirección de y
    !
    !$OMP DO
    do ii = 2, mi-1
        !
        ! Se definen las condiciones de frontera
        !
        BC(1)      =-1.0d0/deltay
        BD(1)      = 1.0d0/deltay
        resulty(1) = flux_aba
        BI(nj)     =-1.d0/deltay
        BC(nj)     = 1.0d0/deltay
        resulty(nj) = flux_arr !(temp_fin+temp_ini)/2.d0
        !
        ! Ensamblado de la matriz tridiagonal
        ! y del vector de resultados
        !
        do jj=2, nj-1
            BI(jj)     =-1.0d0*cond_ter/(deltay*deltay)
            BC(jj)     = 2.0d0*cond_ter*(1.d0/(deltay*deltay)+1.d0/(deltax*deltax))
            BD(jj)     =-1.0d0*cond_ter/(deltay*deltay)
            resulty(jj) = cond_ter/(deltax*deltax)*temper(ii+1,jj)+cond_ter/(deltax*deltax)*temper(ii-1,jj)
        end do
        !
        ! Llamamos al TDMA
        !
        call tri(BI,BC,BD,resulty,tempy,nj)
        do jj =1, nj
            temper(ii,jj) = tempy(jj)
        end do
    end do
    !$OMP END DO
    !
    !
    !$OMP END PARALLEL
    !
    !!! Cerramos la región paralela
    !
    !print*, "Residuo: ", kk, maxval(temper-temp_ant)
    temp_ant = temper
    !
end do
    !
    !! Escritura de resultados
    !
!!$    do ii = 1, mi
!!$        do jj = 1, nj
!!$            write(101,*) xx(ii), yy(jj), temper(ii,jj)
!!$        end do
!!$            write(101,*)
!!$    end do
    !
    ! Obtenemos la diferencia entre la solución analítica y la solución numérica, así como el error porcentual
    !do ii =2, mi-1
    !    do jj = 2, nj-1
    !        !
    !        write(102,*) xx(ii), yy(jj), temper(ii,jj)-(-xx(ii)+308), (temper(ii,jj)-(-xx(ii)+308))/(-xx(ii)+308)
    !    end do
    !        write(102,*)
    !end do
    !
    !
    ! Obtenemos la diferencia entre la solución analítica y la solución numérica, así como el error porcentual
    ! do ii =2, mi-1
    !     do jj = 2, nj-1
    !         !
    !         write(102,*) xx(ii), yy(jj), temper(ii,jj)-(-xx(ii)+308), (temper(ii,jj)-(-xx(ii)+308))/(temp_fin-temp_ini)
    !     end do
    !         write(102,*)
    ! end do
    !
END PROGRAM Calor2D
