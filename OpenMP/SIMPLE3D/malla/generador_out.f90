PROGRAM MEDIDOR
USE constantes
IMPLICIT NONE
INTEGER :: i,j,k,i_o,i_1
!**********************************
! Variables del flujo e inc'ognitas
!y variables de la malla
REAL(kind=DBL), DIMENSION(mi+1,nj+1) :: temp,u,v
REAL(kind=DBL), DIMENSION(mi,nj+1)   :: uo
REAL(kind=DBL), DIMENSION(mi+1,nj)   :: vo
REAL(kind=DBL), DIMENSION(mi+1)      :: x
REAL(kind=DBL), DIMENSION(nj+1)      :: y
REAL(kind=DBL), DIMENSION(mi)        :: xu
REAL(kind=DBL), DIMENSION(nj)        :: yv
!************************************
!Variables para el archivo de entrada
CHARACTER(len=10):: entrada1
CHARACTER(len=22):: entrada2, entrada3
!*******************
entrada1='t08000.dat'
entrada2='out_no76m151_Rxxxu.dat'
entrada3='out_no76m151_Rxxxv.dat'
!*************************
OPEN(unit=14,file=entrada2)
DO j=1,nj+1
  DO i=1,mi
    READ(14,*) xu(i),y(j),uo(i,j)
  END DO
END DO
CLOSE(unit=14)
!*************************
OPEN(unit=15,file=entrada3)
DO j=1,nj
  DO i=1,mi+1
    READ(15,*) x(i),yv(j),vo(i,j)
  END DO
END DO
CLOSE(unit=15)
!*************************
OPEN(unit=13,file=entrada1)
DO j=1,nj+1
  DO i=1,mi+1
    READ(13,*) y(j),x(i),temp(i,j),v(i,j),u(i,j)
  END DO
END DO
CLOSE(unit=13)
u  =-u
uo = 0
vo = 0
!*************
DO i = 2, mi-1
  DO j = 2, nj
    uo(i,j)=2._DBL*u(i,j)-uo(i-1,j)
  END DO
END DO
!*************
DO i = 2, mi
  DO j = 2, nj-1
    vo(i,j)=2._DBL*v(i,j)-vo(i,j-1)
  END DO
END DO
!**********************************
!Formatos de escritura y de lectura
21 FORMAT (1D23.15)
25 FORMAT (5D23.15)
27 FORMAT (8D23.15)
END PROGRAM