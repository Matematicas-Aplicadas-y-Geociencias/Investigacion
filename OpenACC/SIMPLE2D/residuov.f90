SUBROUTINE residuov(xo,yvo,d_yvo,d2_yvo,d_xuo,u_o,v_o,v_anto,pres_o,Resv_o,gamma_uo,dt_o)
USE constantes
IMPLICIT NONE
INTEGER :: i,j,k
!**********************************
! Variables del flujo e inc'ognitas
REAL(kind=DBL), DIMENSION(mi+1,nj), INTENT(inout) :: Resv_o
REAL(kind=DBL), DIMENSION(mi+1,nj), INTENT(in)    :: v_o,v_anto
REAL(kind=DBL), DIMENSION(mi,nj+1), INTENT(in)    :: u_o
REAL(kind=DBL), DIMENSION(mi+1,nj+1), INTENT(in)  :: pres_o
REAL(kind=DBL), DIMENSION(mi+1,nj) :: fv
!**********************************
!Variables de la tridiagonal
REAL(kind=DBL) :: AC,AI,AD,AS,AN
!*********************
!Variables de interpolaci'on
REAL(kind=DBL) :: ui,ud,vn,vs,di,dd,ds,dn,delta_x,delta_y
!*********************
!Variables de la malla,volumen de control,incremento de tiempo y num Richardson
REAL(kind=DBL), DIMENSION(mi+1), INTENT(in) :: xo
REAL(kind=DBL), DIMENSION(nj),   INTENT(in) :: yvo
REAL(kind=DBL), DIMENSION(mi-1), INTENT(in) :: d_xuo
REAL(kind=DBL), DIMENSION(nj-1), INTENT(in) :: d_yvo,d2_yvo
REAL(kind=DBL), INTENT(in) :: dt_o,gamma_uo
!********************
!C'alculo del residuo
Resv_o=0._DBL
DO j=2,nj-1
  DO i=2,mi
    !**************************
    !Interpolaciones necesarias
    ud=u_o(i,j)+d_yvo(j-1)/d2_yvo(j)*(u_o(i,j+1)-u_o(i,j))
    ui=u_o(i-1,j)+d_yvo(j-1)/d2_yvo(j)*(u_o(i-1,j+1)-u_o(i-1,j))
    vn=(v_o(i,j+1)+v_o(i,j))/2._DBL
    vs=(v_o(i,j)+v_o(i,j-1))/2._DBL
    di=xo(i)-xo(i-1)
    dd=xo(i+1)-xo(i)
    ds=d_yvo(j-1)
    dn=d_yvo(j)
    delta_x=d_xuo(i-1)
    delta_y=d2_yvo(j)/2._DBL
    !************************
    AI=gamma_uo*delta_y/di*DMAX1(0._DBL,(1._DBL-0.1_DBL*dabs(ui*di/gamma_uo))**5._DBL)+&
    &DMAX1(0._DBL,ui*delta_y)
    AD=gamma_uo*delta_y/dd*DMAX1(0._DBL,(1._DBL-0.1_DBL*dabs(ud*dd/gamma_uo))**5._DBL)+&
    &DMAX1(0._DBL,-ud*delta_y)
    AS=gamma_uo*delta_x/ds*DMAX1(0._DBL,(1._DBL-0.1_DBL*dabs(vs*ds/gamma_uo))**5._DBL)+&
    &DMAX1(0._DBL,vs*delta_x)
    AN=gamma_uo*delta_x/dn*DMAX1(0._DBL,(1._DBL-0.1_DBL*dabs(vn*dn/gamma_uo))**5._DBL)+&
    &DMAX1(0._DBL,-vn*delta_x)
    AC=AI+AD+AS+AN+delta_x*delta_y/dt_o
    Resv_o(i,j)=AI*v_o(i-1,j)+AD*v_o(i+1,j)+AS*v_o(i,j-1)+AN*v_o(i,j+1)+delta_x*delta_y*v_anto(i,j)/dt_o+&
    &(pres_o(i,j)-pres_o(i,j+1))*delta_x-AC*v_o(i,j)
    !**************************
  END DO
END DO
END SUBROUTINE
