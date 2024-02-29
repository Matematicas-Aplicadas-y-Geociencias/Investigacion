SUBROUTINE nusselt(xo,yo,d_xuo,d_yvo,temp_o,nusselt0_o,nusselt1_o,i_oo,i_1o)
USE constantes
IMPLICIT NONE
INTEGER :: i,j,k
INTEGER, INTENT(in) :: i_oo,i_1o
!*****************************************
!Variables de malla, nusselt y temperatura
REAL(kind=DBL), DIMENSION(mi+1,nj+1), INTENT(in) :: temp_o
REAL(kind=DBL), DIMENSION(mi+1), INTENT(in)      :: xo
REAL(kind=DBL), DIMENSION(nj+1), INTENT(in)      :: yo
REAL(kind=DBL), DIMENSION(mi-1), INTENT(in)      :: d_xuo
REAL(kind=DBL), DIMENSION(nj-1), INTENT(in)      :: d_yvo
REAL(kind=DBL), INTENT(out) :: nusselt0_o, nusselt1_o
!******************************************
!Variables de interpolaci'on para derivadas
REAL(kind=DBL) a,b,c,dx1,dx2,dy1,dy2
REAL(kind=DBL), DIMENSION(mi+1) :: derivada
!*******************************************
!Variables de interpolaci'on para integrales
REAL(kind=DBL) alpha,beta,gamma
REAL(kind=DBL) x_o,x_1,x_2,y_o,y_1,y_2
  !******************
  !C'alculos para y=0
  !********************************
  !C'alculo de la derivada respecto
  !a la normal con interpolaci'on 
  !cuadr'atica
  derivada = 0._DBL
  DO i = i_oo, i_1o
    j = 1
    dy1 = temp_o(i,j+1)-temp_o(i,j)
    dy2 = temp_o(i,j+2)-temp_o(i,j+1)
    dx1 = yo(j+1)-yo(j)
    dx2 = yo(j+2)-yo(j+1)
    a = (dy2/dx2-dy1/dx1)/(dx1+dx2)
    b = dy2/dx2-a*(yo(j+1)+yo(j+2))
    c = temp_o(i,j)-a*yo(j)*yo(j)-b*yo(j)
    derivada(i) = 2._DBL*a*yo(j)+b
  END DO
  !*****************************************
  !C'alculo de la integral con interpolacion
  !cuadrática
  nusselt0_o = 0._DBL
  DO i = i_oo, i_1o, 2
    x_o   = xo(i)
    x_1   = xo(i+1)
    x_2   = xo(i+2)
    y_o   = derivada(i)
    y_1   = derivada(i+1)
    y_2   = derivada(i+2)
    alpha = ((y_2-y_1)/(x_2-x_1)-(y_1-y_o)/(x_1-x_o))/(x_2-x_o)
    beta  = (y_2-y_1)/(x_2-x_1)-alpha*(x_1+x_2)
    gamma = y_o-alpha*x_o*x_o-beta*x_o
    nusselt0_o = nusselt0_o+alpha/3._DBL*(x_2*x_2*x_2-x_o*x_o*x_o)+beta/2._DBL*(x_2*x_2-x_o*x_o)+gamma*(x_2-x_o)
  END DO
  nusselt0_o=-nusselt0_o
  !******************
  !C'alculos para y=1
  !********************************
  !C'alculo de la derivada respecto 
  !a la normal con interpolaci'on 
  !cuadr'atica
  derivada = 0._DBL
  DO i = i_oo, i_1o
    j = nj-1
    dy1 = temp_o(i,j+1)-temp_o(i,j)
    dy2 = temp_o(i,j+2)-temp_o(i,j+1)
    dx1 = yo(j+1)-yo(j)
    dx2 = yo(j+2)-yo(j+1)
    a = (dy2/dx2-dy1/dx1)/(dx1+dx2)
    b = dy2/dx2-a*(yo(j+1)+yo(j+2))
    c = temp_o(i,j)-a*yo(j)*yo(j)-b*yo(j)
    derivada(i) = 2._DBL*a*yo(j+2)+b
  END DO
  !*****************************************
  !C'alculo de la integral con interpolacion
  !cuadrática
  nusselt1_o = 0._DBL
  DO i = i_oo, i_1o, 2
    x_o = xo(i)
    x_1 = xo(i+1)
    x_2 = xo(i+2)
    y_o = derivada(i)
    y_1 = derivada(i+1)
    y_2 = derivada(i+2)
    alpha = ((y_2-y_1)/(x_2-x_1)-(y_1-y_o)/(x_1-x_o))/(x_2-x_o)
    beta  = (y_2-y_1)/(x_2-x_1)-alpha*(x_1+x_2)
    gamma = y_o-alpha*x_o*x_o-beta*x_o
    nusselt1_o = nusselt1_o+alpha/3._DBL*(x_2*x_2*x_2-x_o*x_o*x_o)+beta/2._DBL*(x_2*x_2-x_o*x_o)+gamma*(x_2-x_o)
  END DO
  nusselt1_o=-nusselt1_o
END SUBROUTINE
