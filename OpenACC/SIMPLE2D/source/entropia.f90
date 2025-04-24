SUBROUTINE entropia_cvt(xo,yo,u_o,xuo,v_o,yvo,temp_o,entr_calor_o,entr_viscosa_o,entr_o,entr_int_o,temp_int_o,a_ent_o,lambda_ent_o)
USE constantes
IMPLICIT NONE
INTEGER :: i,j,k
!**********************************
! Variables del flujo e inc'ognitas
!y variables de la malla
REAL(kind=DBL), DIMENSION(mi+1,nj+1), INTENT(in)  :: temp_o
REAL(kind=DBL), DIMENSION(mi,nj+1),   INTENT(in)  :: u_o
REAL(kind=DBL), DIMENSION(mi+1,nj),   INTENT(in)  :: v_o
REAL(kind=DBL), DIMENSION(mi+1),      INTENT(in)  :: xo
REAL(kind=DBL), DIMENSION(nj+1),      INTENT(in)  :: yo
REAL(kind=DBL), DIMENSION(mi),        INTENT(in)  :: xuo
REAL(kind=DBL), DIMENSION(nj),        INTENT(in)  :: yvo
REAL(kind=DBL), DIMENSION(mi+1,nj+1), INTENT(out) :: entr_o,entr_calor_o,entr_viscosa_o
REAL(kind=DBL), INTENT(in)  :: a_ent_o,lambda_ent_o
REAL(kind=DBL), INTENT(out) :: entr_int_o,temp_int_o
!***************************
!Variables de interpolaci'on
REAL(kind=DBL), DIMENSION(mi+1,nj+1) :: du_dx,du_dy,dv_dx,dv_dy,dtemp_dx,dtemp_dy
REAL(kind=DBL), DIMENSION(mi,nj+1)   :: du_dyo
REAL(kind=DBL), DIMENSION(mi+1,nj)   :: dv_dxo
REAL(kind=DBL), DIMENSION(nj+1)      :: temp_intx,areax,entr_cal_intx,entr_vis_intx
REAL(kind=DBL) :: a,b,c,dy1,dy2,dx1,dx2,x_o,x_1,x_2,y_o,y_1,y_2,alpha,beta,gamma,area
REAL(kind=DBL) :: entr_cal_int,entr_vis_int
!**********************
entr_calor_o   = 0._DBL
entr_viscosa_o = 0._DBL
entr_o         = 0._DBL
!**************************
!C'alculo de derivadas de u
DO j = 1, nj+1
  i = 1
  dy1 = u_o(i+1,j)-u_o(i,j)
  dy2 = u_o(i+2,j)-u_o(i+1,j)
  dx1 = xuo(i+1)-xuo(i)
  dx2 = xuo(i+2)-xuo(i+1)
  a   = (dy2/dx2-dy1/dx1)/(dx1+dx2)
  b   = dy2/dx2-a*(xuo(i+1)+xuo(i+2))
  c   = u_o(i,j)-a*xuo(i)*xuo(i)-b*xuo(i)
  du_dx(i,j) = 2._DBL*a*xo(i)+b
  du_dx(2,j) = 2._DBL*a*xo(2)+b
  DO i = 2, mi-2
    !**************************
    !Interpolaciones necesarias
    !derivada en x
    dy1 = u_o(i+1,j)-u_o(i,j)
    dy2 = u_o(i+2,j)-u_o(i+1,j)
    dx1 = xuo(i+1)-xuo(i)
    dx2 = xuo(i+2)-xuo(i+1)
    a   = (dy2/dx2-dy1/dx1)/(dx1+dx2)
    b   = dy2/dx2-a*(xuo(i+1)+xuo(i+2))
    c   = u_o(i,j)-a*xuo(i)*xuo(i)-b*xuo(i)
    du_dx(i+1,j) = 2._DBL*a*xo(i+1)+b
  END DO
  du_dx(mi,j) = 2._DBL*a*xo(mi+1)+b
END DO
DO i = 1, mi
  j = 1
  dy1 = u_o(i,j+1)-u_o(i,j)
  dy2 = u_o(i,j+2)-u_o(i,j+1)
  dx1 = yo(j+1)-yo(j)
  dx2 = yo(j+2)-yo(j+1)
  a   = (dy2/dx2-dy1/dx1)/(dx1+dx2)
  b   = dy2/dx2-a*(yo(j+1)+yo(j+2))
  c   = u_o(i,j)-a*yo(j)*yo(j)-b*yo(j)
  du_dyo(i,j) = 2._DBL*a*yo(j)+b
  du_dyo(i,2) = 2._DBL*a*yo(2)+b
  DO j = 2, nj-1
    !**************************
    !Interpolaciones necesarias
    !derivada en y
    dy1 = u_o(i,j+1)-u_o(i,j)
    dy2 = u_o(i,j+2)-u_o(i,j+1)
    dx1 = yo(j+1)-yo(j)
    dx2 = yo(j+2)-yo(j+1)
    a   = (dy2/dx2-dy1/dx1)/(dx1+dx2)
    b   = dy2/dx2-a*(yo(j+1)+yo(j+2))
    c   = u_o(i,j)-a*yo(j)*yo(j)-b*yo(j)
    du_dyo(i,j+1) = 2._DBL*a*yo(j+1)+b
  END DO
  du_dyo(i,nj+1) = 2._DBL*a*yo(nj+1)+b
END DO
!*************
DO j = 1, nj+1
  du_dy(1,j) = du_dyo(1,j)
  DO i = 2, mi
    du_dy(i,j) = (du_dyo(i,j)+du_dyo(i-1,j))/2._DBL
  END DO
  du_dy(mi+1,j) = du_dyo(mi,j)
END DO
!**************************
!C'alculo de derivadas de v
DO j = 1, nj
  i = 1
  dy1 = v_o(i+1,j)-v_o(i,j)
  dy2 = v_o(i+2,j)-v_o(i+1,j)
  dx1 = xo(i+1)-xo(i)
  dx2 = xo(i+2)-xo(i+1)
  a   = (dy2/dx2-dy1/dx1)/(dx1+dx2)
  b   = dy2/dx2-a*(xo(i+1)+xo(i+2))
  c   = v_o(i,j)-a*xo(i)*xo(i)-b*xo(i)
  dv_dxo(i,j) = 2._DBL*a*xo(i)+b
  dv_dxo(2,j) = 2._DBL*a*xo(2)+b
  DO i = 2, mi-1
    !**************************
    !Interpolaciones necesarias
    !derivada en x
    dy1 = v_o(i+1,j)-v_o(i,j)
    dy2 = v_o(i+2,j)-v_o(i+1,j)
    dx1 = xo(i+1)-xo(i)
    dx2 = xo(i+2)-xo(i+1)
    a   = (dy2/dx2-dy1/dx1)/(dx1+dx2)
    b   = dy2/dx2-a*(xo(i+1)+xo(i+2))
    c   = v_o(i,j)-a*xo(i)*xo(i)-b*xo(i)
    dv_dxo(i+1,j) = 2._DBL*a*xo(i+1)+b
  END DO
  dv_dxo(mi+1,j) = 2._DBL*a*xo(mi+1)+b
END DO
!*************
DO i = 1, mi+1
  dv_dx(i,1) = dv_dxo(i,1)
  DO j = 2, nj
    dv_dx(i,j) = (dv_dxo(i,j)+dv_dxo(i,j-1))/2._DBL
  END DO
  dv_dx(i,nj+1) = dv_dxo(i,nj)
END DO
!*************
DO i = 1, mi+1
  j = 1
  dy1 = v_o(i,j+1)-v_o(i,j)
  dy2 = v_o(i,j+2)-v_o(i,j+1)
  dx1 = yvo(j+1)-yvo(j)
  dx2 = yvo(j+2)-yvo(j+1)
  a   = (dy2/dx2-dy1/dx1)/(dx1+dx2)
  b   = dy2/dx2-a*(yvo(j+1)+yvo(j+2))
  c   = v_o(i,j)-a*yvo(j)*yvo(j)-b*yvo(j)
  dv_dy(i,j) = 2._DBL*a*yo(j)+b
  dv_dy(i,2) = 2._DBL*a*yo(2)+b
  DO j = 2, nj-2
    !**************************
    !Interpolaciones necesarias
    !derivada en y
    dy1 = v_o(i,j+1)-v_o(i,j)
    dy2 = v_o(i,j+2)-v_o(i,j+1)
    dx1 = yvo(j+1)-yvo(j)
    dx2 = yvo(j+2)-yvo(j+1)
    a   = (dy2/dx2-dy1/dx1)/(dx1+dx2)
    b   = dy2/dx2-a*(yvo(j+1)+yvo(j+2))
    c   = v_o(i,j)-a*yvo(j)*yvo(j)-b*yvo(j)
    dv_dy(i,j+1) = 2._DBL*a*yo(j+1)+b
  END DO
  dv_dy(i,nj) = 2._DBL*a*yo(nj+1)+b
END DO
!*****************************
!C'alculo de derivadas de temp
DO j = 1, nj+1
  i = 1
    dy1 = temp_o(i+1,j)-temp_o(i,j)
    dy2 = temp_o(i+2,j)-temp_o(i+1,j)
    dx1 = xo(i+1)-xo(i)
    dx2 = xo(i+2)-xo(i+1)
    a   = (dy2/dx2-dy1/dx1)/(dx1+dx2)
    b   = dy2/dx2-a*(xo(i+1)+xo(i+2))
    c   = temp_o(i,j)-a*xo(i)*xo(i)-b*xo(i)
    dtemp_dx(i,j) = 2._DBL*a*xo(i)+b
    dtemp_dx(2,j) = 2._DBL*a*xo(2)+b
  DO i = 2, mi-1
    !**************************
    !Interpolaciones necesarias
    !derivada en x
    dy1 = temp_o(i+1,j)-temp_o(i,j)
    dy2 = temp_o(i+2,j)-temp_o(i+1,j)
    dx1 = xo(i+1)-xo(i)
    dx2 = xo(i+2)-xo(i+1)
    a   = (dy2/dx2-dy1/dx1)/(dx1+dx2)
    b   = dy2/dx2-a*(xo(i+1)+xo(i+2))
    c   = temp_o(i,j)-a*xo(i)*xo(i)-b*xo(i)
    dtemp_dx(i+1,j) = 2._DBL*a*xo(i+1)+b
  END DO
  dtemp_dx(mi+1,j) = 2._DBL*a*xo(mi+1)+b
END DO
!*************
DO i = 1, mi+1
  j = 1
    dy1 = temp_o(i,j+1)-temp_o(i,j)
    dy2 = temp_o(i,j+2)-temp_o(i,j+1)
    dx1 = yo(j+1)-yo(j)
    dx2 = yo(j+2)-yo(j+1)
    a   = (dy2/dx2-dy1/dx1)/(dx1+dx2)
    b   = dy2/dx2-a*(yo(j+1)+yo(j+2))
    c   = temp_o(i,j)-a*yo(j)*yo(j)-b*yo(j)
    dtemp_dy(i,j) = 2._DBL*a*yo(j)+b
    dtemp_dy(i,2) = 2._DBL*a*yo(2)+b
  DO j = 2, nj-1
    !**************************
    !Interpolaciones necesarias
    !derivada en y
    dy1 = temp_o(i,j+1)-temp_o(i,j)
    dy2 = temp_o(i,j+2)-temp_o(i,j+1)
    dx1 = yo(j+1)-yo(j)
    dx2 = yo(j+2)-yo(j+1)
    a   = (dy2/dx2-dy1/dx1)/(dx1+dx2)
    b   = dy2/dx2-a*(yo(j+1)+yo(j+2))
    c   = temp_o(i,j)-a*yo(j)*yo(j)-b*yo(j)
    dtemp_dy(i,j+1) = 2._DBL*a*yo(j+1)+b
  END DO
  dtemp_dy(i,nj+1) = 2._DBL*a*yo(nj+1)+b
END DO 
!***************************************
!c'alculo de prod. entropia flujo calor,
!prod entropia viscosa, entropia total.
DO i = 1, mi+1
  DO j = 1, nj+1
    entr_calor_o(i,j) = (1._DBL/((temp_o(i,j)+a_ent_o)*(temp_o(i,j)+a_ent_o)))*&
    &(dtemp_dx(i,j)*dtemp_dx(i,j)+dtemp_dy(i,j)*dtemp_dy(i,j))
    entr_viscosa_o(i,j) = (lambda_ent_o/(2._DBL*(temp_o(i,j)+a_ent_o)))*(4._DBL*(du_dx(i,j)*du_dx(i,j)+&
    &dv_dy(i,j)*dv_dy(i,j))+2._DBL*(du_dy(i,j)*du_dy(i,j)+dv_dx(i,j)*dv_dx(i,j)+&
    &2._DBL*du_dy(i,j)*dv_dx(i,j)))
    entr_o(i,j) = entr_calor_o(i,j) + entr_viscosa_o(i,j)
  END DO
END DO
!******************************************
!******************************************
!Se calculan las integrales de la entrop'ia
!y de la temperatura
!******************************************
!se calcula el area sobre la que se integra
areax = 0._DBL
area  = 0._DBL
DO j = 1, nj+1
  DO i = 1, mi-1, 2
    x_o = xo(i)
    x_1 = xo(i+1)
    x_2 = xo(i+2)
    y_o = 1.0_DBL
    y_1 = 1.0_DBL
    y_2 = 1.0_DBL
    alpha = ((y_2-y_1)/(x_2-x_1)-(y_1-y_o)/(x_1-x_o))/(x_2-x_o)
    beta  = (y_2-y_1)/(x_2-x_1)-alpha*(x_1+x_2)
    gamma = y_o-alpha*x_o*x_o-beta*x_o
    areax(j) = areax(j)+alpha/3._DBL*(x_2*x_2*x_2-x_o*x_o*x_o)+beta/2._DBL*(x_2*x_2-x_o*x_o)+gamma*(x_2-x_o)
  END DO
  IF(mod(mi+1,2)==0)THEN
    i = mi-1
    x_o = xo(i)
    x_1 = xo(i+1)
    x_2 = xo(i+2)
    y_o = 1.0_DBL
    y_1 = 1.0_DBL
    y_2 = 1.0_DBL
    alpha = ((y_2-y_1)/(x_2-x_1)-(y_1-y_o)/(x_1-x_o))/(x_2-x_o)
    beta  = (y_2-y_1)/(x_2-x_1)-alpha*(x_1+x_2)
    gamma = y_o-alpha*x_o*x_o-beta*x_o
    areax(j) = areax(j)+alpha/3._DBL*(x_2*x_2*x_2-x_1*x_1*x_1)+beta/2._DBL*(x_2*x_2-x_1*x_1)+gamma*(x_2-x_1)
  END IF
END DO
! !***************************
!integral area en direcci'on y
DO j = 1, nj-1, 2
  x_o = yo(j)
  x_1 = yo(j+1)
  x_2 = yo(j+2)
  y_o = areax(j)
  y_1 = areax(j+1)
  y_2 = areax(j+2)
  alpha = ((y_2-y_1)/(x_2-x_1)-(y_1-y_o)/(x_1-x_o))/(x_2-x_o)
  beta  = (y_2-y_1)/(x_2-x_1)-alpha*(x_1+x_2)
  gamma = y_o-alpha*x_o*x_o-beta*x_o
  area  = area+alpha/3._DBL*(x_2*x_2*x_2-x_o*x_o*x_o)+beta/2._DBL*(x_2*x_2-x_o*x_o)+gamma*(x_2-x_o)
END DO
IF(mod(nj+1,2)==0)THEN
  j = nj-1
  x_o = yo(j)
  x_1 = yo(j+1)
  x_2 = yo(j+2)
  y_o = areax(j)
  y_1 = areax(j+1)
  y_2 = areax(j+2)
  alpha = ((y_2-y_1)/(x_2-x_1)-(y_1-y_o)/(x_1-x_o))/(x_2-x_o)
  beta  = (y_2-y_1)/(x_2-x_1)-alpha*(x_1+x_2)
  gamma = y_o-alpha*x_o*x_o-beta*x_o
  area  = area+alpha/3._DBL*(x_2*x_2*x_2-x_1*x_1*x_1)+beta/2._DBL*(x_2*x_2-x_1*x_1)+gamma*(x_2-x_1)
END IF
!*************************************
!*************************************
!integral de la entropia por flujos de 
!calor en direcci'on x
entr_cal_int  = 0._DBL
entr_cal_intx = 0._DBL
DO j = 1, nj+1
  DO i = 1, mi-1, 2
    x_o = xo(i)
    x_1 = xo(i+1)
    x_2 = xo(i+2)
    y_o = entr_calor_o(i,j)
    y_1 = entr_calor_o(i+1,j)
    y_2 = entr_calor_o(i+2,j)
    alpha = ((y_2-y_1)/(x_2-x_1)-(y_1-y_o)/(x_1-x_o))/(x_2-x_o)
    beta  = (y_2-y_1)/(x_2-x_1)-alpha*(x_1+x_2)
    gamma = y_o-alpha*x_o*x_o-beta*x_o
    entr_cal_intx(j) = entr_cal_intx(j)+alpha/3._DBL*(x_2*x_2*x_2-x_o*x_o*x_o)+beta/2._DBL*(x_2*x_2-x_o*x_o)+gamma*(x_2-x_o)
  END DO
  IF(mod(mi+1,2)==0)THEN
    i = mi-1
    x_o = xo(i)
    x_1 = xo(i+1)
    x_2 = xo(i+2)
    y_o = entr_calor_o(i,j)
    y_1 = entr_calor_o(i+1,j)
    y_2 = entr_calor_o(i+2,j)
    alpha = ((y_2-y_1)/(x_2-x_1)-(y_1-y_o)/(x_1-x_o))/(x_2-x_o)
    beta  = (y_2-y_1)/(x_2-x_1)-alpha*(x_1+x_2)
    gamma = y_o-alpha*x_o*x_o-beta*x_o
    entr_cal_intx(j) = entr_cal_intx(j)+alpha/3._DBL*(x_2*x_2*x_2-x_1*x_1*x_1)+beta/2._DBL*(x_2*x_2-x_1*x_1)+gamma*(x_2-x_1)
  END IF
END DO
! !************************************
!integral entropia_c en la direcci'on y
DO j = 1, nj-1, 2
  x_o = yo(j)
  x_1 = yo(j+1)
  x_2 = yo(j+2)
  y_o = entr_cal_intx(j)
  y_1 = entr_cal_intx(j+1)
  y_2 = entr_cal_intx(j+2)
  alpha = ((y_2-y_1)/(x_2-x_1)-(y_1-y_o)/(x_1-x_o))/(x_2-x_o)
  beta  = (y_2-y_1)/(x_2-x_1)-alpha*(x_1+x_2)
  gamma = y_o-alpha*x_o*x_o-beta*x_o
  entr_cal_int = entr_cal_int+alpha/3._DBL*(x_2*x_2*x_2-x_o*x_o*x_o)+beta/2._DBL*(x_2*x_2-x_o*x_o)+gamma*(x_2-x_o)
END DO
IF(mod(nj+1,2)==0)THEN
  j = nj-1
  x_o = yo(j)
  x_1 = yo(j+1)
  x_2 = yo(j+2)
  y_o = entr_cal_intx(j)
  y_1 = entr_cal_intx(j+1)
  y_2 = entr_cal_intx(j+2)
  alpha = ((y_2-y_1)/(x_2-x_1)-(y_1-y_o)/(x_1-x_o))/(x_2-x_o)
  beta  = (y_2-y_1)/(x_2-x_1)-alpha*(x_1+x_2)
  gamma = y_o-alpha*x_o*x_o-beta*x_o
  entr_cal_int = entr_cal_int+alpha/3._DBL*(x_2*x_2*x_2-x_1*x_1*x_1)+beta/2._DBL*(x_2*x_2-x_1*x_1)+gamma*(x_2-x_1)
END IF
!***********************************
!***********************************
!integral de la entropia por efectos 
!viscosos en direcci'on x
entr_vis_int  = 0._DBL
entr_vis_intx = 0._DBL
DO j = 1, nj+1
  DO i = 1, mi-1, 2
    x_o = xo(i)
    x_1 = xo(i+1)
    x_2 = xo(i+2)
    y_o = entr_viscosa_o(i,j)
    y_1 = entr_viscosa_o(i+1,j)
    y_2 = entr_viscosa_o(i+2,j)
    alpha = ((y_2-y_1)/(x_2-x_1)-(y_1-y_o)/(x_1-x_o))/(x_2-x_o)
    beta  = (y_2-y_1)/(x_2-x_1)-alpha*(x_1+x_2)
    gamma = y_o-alpha*x_o*x_o-beta*x_o
    entr_vis_intx(j) = entr_vis_intx(j)+alpha/3._DBL*(x_2*x_2*x_2-x_o*x_o*x_o)+beta/2._DBL*(x_2*x_2-x_o*x_o)+gamma*(x_2-x_o)
  END DO
  IF(mod(mi+1,2)==0)THEN
    i = mi-1
    x_o = xo(i)
    x_1 = xo(i+1)
    x_2 = xo(i+2)
    y_o = entr_viscosa_o(i,j)
    y_1 = entr_viscosa_o(i+1,j)
    y_2 = entr_viscosa_o(i+2,j)
    alpha = ((y_2-y_1)/(x_2-x_1)-(y_1-y_o)/(x_1-x_o))/(x_2-x_o)
    beta  = (y_2-y_1)/(x_2-x_1)-alpha*(x_1+x_2)
    gamma = y_o-alpha*x_o*x_o-beta*x_o
    entr_vis_intx(j) = entr_vis_intx(j)+alpha/3._DBL*(x_2*x_2*x_2-x_1*x_1*x_1)+beta/2._DBL*(x_2*x_2-x_1*x_1)+gamma*(x_2-x_1)
  END IF
END DO
! !************************************
!integral entropia_v en la direcci'on y
DO j = 1, nj-1, 2
  x_o = yo(j)
  x_1 = yo(j+1)
  x_2 = yo(j+2)
  y_o = entr_vis_intx(j)
  y_1 = entr_vis_intx(j+1)
  y_2 = entr_vis_intx(j+2)
  alpha = ((y_2-y_1)/(x_2-x_1)-(y_1-y_o)/(x_1-x_o))/(x_2-x_o)
  beta  = (y_2-y_1)/(x_2-x_1)-alpha*(x_1+x_2)
  gamma = y_o-alpha*x_o*x_o-beta*x_o
  entr_vis_int = entr_vis_int+alpha/3._DBL*(x_2*x_2*x_2-x_o*x_o*x_o)+beta/2._DBL*(x_2*x_2-x_o*x_o)+gamma*(x_2-x_o)
END DO
IF(mod(nj+1,2)==0)THEN
  j = nj-1
  x_o = yo(j)
  x_1 = yo(j+1)
  x_2 = yo(j+2)
  y_o = entr_vis_intx(j)
  y_1 = entr_vis_intx(j+1)
  y_2 = entr_vis_intx(j+2)
  alpha = ((y_2-y_1)/(x_2-x_1)-(y_1-y_o)/(x_1-x_o))/(x_2-x_o)
  beta  = (y_2-y_1)/(x_2-x_1)-alpha*(x_1+x_2)
  gamma = y_o-alpha*x_o*x_o-beta*x_o
  entr_vis_int = entr_vis_int+alpha/3._DBL*(x_2*x_2*x_2-x_1*x_1*x_1)+beta/2._DBL*(x_2*x_2-x_1*x_1)+gamma*(x_2-x_1)
END IF  
entr_int_o = entr_cal_int + entr_vis_int
!*****************************************
!*****************************************
!integral de la temperatura en direccion x
temp_int_o = 0._DBL
temp_intx  = 0._DBL
DO j = 1, nj+1
  DO i = 1, mi-1, 2
    x_o = xo(i)
    x_1 = xo(i+1)
    x_2 = xo(i+2)
    y_o = temp_o(i,j)
    y_1 = temp_o(i+1,j)
    y_2 = temp_o(i+2,j)
    alpha = ((y_2-y_1)/(x_2-x_1)-(y_1-y_o)/(x_1-x_o))/(x_2-x_o)
    beta  = (y_2-y_1)/(x_2-x_1)-alpha*(x_1+x_2)
    gamma = y_o-alpha*x_o*x_o-beta*x_o
    temp_intx(j) = temp_intx(j)+alpha/3._DBL*(x_2*x_2*x_2-x_o*x_o*x_o)+beta/2._DBL*(x_2*x_2-x_o*x_o)+gamma*(x_2-x_o)
  END DO
  IF(mod(mi+1,2)==0)THEN
    i=mi-1
    x_o = xo(i)
    x_1 = xo(i+1)
    x_2 = xo(i+2)
    y_o = temp_o(i,j)
    y_1 = temp_o(i+1,j)
    y_2 = temp_o(i+2,j)
    alpha = ((y_2-y_1)/(x_2-x_1)-(y_1-y_o)/(x_1-x_o))/(x_2-x_o)
    beta  = (y_2-y_1)/(x_2-x_1)-alpha*(x_1+x_2)
    gamma = y_o-alpha*x_o*x_o-beta*x_o
    temp_intx(j) = temp_intx(j)+alpha/3._DBL*(x_2*x_2*x_2-x_1*x_1*x_1)+beta/2._DBL*(x_2*x_2-x_1*x_1)+gamma*(x_2-x_1)
  END IF
END DO
!*****************************************
!integral de la temperatura en direccion y
DO j = 1, nj-1, 2
  x_o = yo(j)
  x_1 = yo(j+1)
  x_2 = yo(j+2)
  y_o = temp_intx(j)
  y_1 = temp_intx(j+1)
  y_2 = temp_intx(j+2)
  alpha = ((y_2-y_1)/(x_2-x_1)-(y_1-y_o)/(x_1-x_o))/(x_2-x_o)
  beta  = (y_2-y_1)/(x_2-x_1)-alpha*(x_1+x_2)
  gamma = y_o-alpha*x_o*x_o-beta*x_o
  temp_int_o = temp_int_o+alpha/3._DBL*(x_2*x_2*x_2-x_o*x_o*x_o)+beta/2._DBL*(x_2*x_2-x_o*x_o)+gamma*(x_2-x_o)
END DO
IF(mod(nj+1,2)==0)THEN
  j = nj-1
  x_o = yo(j)
  x_1 = yo(j+1)
  x_2 = yo(j+2)
  y_o = temp_intx(j)
  y_1 = temp_intx(j+1)
  y_2 = temp_intx(j+2)
  alpha = ((y_2-y_1)/(x_2-x_1)-(y_1-y_o)/(x_1-x_o))/(x_2-x_o)
  beta  = (y_2-y_1)/(x_2-x_1)-alpha*(x_1+x_2)
  gamma = y_o-alpha*x_o*x_o-beta*x_o
  temp_int_o = temp_int_o+alpha/3._DBL*(x_2*x_2*x_2-x_1*x_1*x_1)+beta/2._DBL*(x_2*x_2-x_1*x_1)+gamma*(x_2-x_1)
END IF
temp_int_o = temp_int_o / area
entr_int_o = entr_int_o / area
!*************
END SUBROUTINE