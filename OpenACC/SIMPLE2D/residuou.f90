SUBROUTINE residuou(res_fluido_uo,xuo,yo,fey,d_xuo,d2_xuo,d_yvo,u_o,u_anto,v_o,temp_o,pres_o,Resu_o,gamma_uo,Ri_o,dt_o)
USE malla
IMPLICIT NONE
!$acc routine
INTEGER :: i,j,k
!**********************************
! Variables del flujo e inc'ognitas
REAL(kind=DBL), DIMENSION(mi,nj+1), INTENT(out)   :: Resu_o
REAL(kind=DBL), DIMENSION(mi,nj+1), INTENT(in)    :: u_o,u_anto,gamma_uo,Ri_o
REAL(kind=DBL), DIMENSION(mi+1,nj), INTENT(in)    :: v_o
REAL(kind=DBL), DIMENSION(mi+1,nj+1), INTENT(in)  :: temp_o,pres_o
!**********************************
!Variables de la tridiagonal
REAL(kind=DBL) :: AC,AI,AD,AS,AN
!*********************
!Variables de interpolaci'on
REAL(kind=DBL) :: ui,ud,vn,vs,di,dd,ds,dn,gamma_i,gamma_d,gamma_s,gamma_n
REAL(kind=DBL) :: delta_x,delta_y,temp_int
!*********************
!Variables de la malla,volumen de control,incremento de tiempo y num Richardson
REAL(kind=DBL), DIMENSION(mi),   INTENT(in) :: xuo
REAL(kind=DBL), DIMENSION(nj+1), INTENT(in) :: yo,fey
REAL(kind=DBL), DIMENSION(mi-1), INTENT(in) :: d_xuo,d2_xuo
REAL(kind=DBL), DIMENSION(nj-1), INTENT(in) :: d_yvo
REAL(kind=DBL), INTENT(in) :: dt_o
LOGICAL, INTENT(out)       :: res_fluido_uo
!********************
!C'alculo de residuos
Resu_o=0._DBL
res_fluido_uo = .TRUE.
!$OMP  PARALLEL DO DEFAULT(NONE)&
!$OMP& PRIVATE(AI,AC,AD,AS,AN,ui,ud,vs,vn,di,dd,ds,dn,gamma_i,gamma_d,gamma_s,gamma_n,delta_x,delta_y,temp_int)&
!$OMP& SHARED(Resu_o,u_o,u_anto,pres_o,gamma_uo,v_o,d_xuo,d2_xuo,yo,fey,d_yvo,temp_o,Ri_o,dt_o)
DO j=2,nj
  DO i=2,mi-1
    !**************************
    !Interpolaciones necesarias
    ui=(u_o(i-1,j)+u_o(i,j))/2._DBL
    ud=(u_o(i,j)+u_o(i+1,j))/2._DBL
    vs=v_o(i,j-1)+d_xuo(i-1)/d2_xuo(i)*(v_o(i+1,j-1)-v_o(i,j-1))
    vn=v_o(i,j)+d_xuo(i-1)/d2_xuo(i)*(v_o(i+1,j)-v_o(i,j))
    di=d_xuo(i-1)
    dd=d_xuo(i)
    ds=yo(j)-yo(j-1)
    dn=yo(j+1)-yo(j)
    gamma_i = 2._DBL * gamma_uo(i-1,j) * gamma_uo(i,j) / (gamma_uo(i-1,j) + gamma_uo(i,j))
    gamma_d = 2._DBL * gamma_uo(i+1,j) * gamma_uo(i,j) / (gamma_uo(i+1,j) + gamma_uo(i,j))
    gamma_s = 1._DBL/((1._DBL-fey(j)) / gamma_uo(i,j-1) + fey(j) / gamma_uo(i,j))
    gamma_n = 1._DBL/((1._DBL-fey(j+1)) / gamma_uo(i,j) + fey(j+1) / gamma_uo(i,j+1))
    delta_x=d2_xuo(i)/2._DBL
    delta_y=d_yvo(j-1)
    temp_int=(1._DBL - di/(2._DBL*delta_x)) * temp_o(i,j) + di/(2._DBL*delta_x) * temp_o(i+1,j)
!     (temp_o(i,j)+temp_o(i+1,j))/2._DBL
    !************************
    AI=gamma_i*delta_y/di*DMAX1(0._DBL,(1._DBL-0.1_DBL*DABS(ui*di/gamma_i))**5._DBL)+&
    &DMAX1(0._DBL,ui*delta_y)
    AD=gamma_d*delta_y/dd*DMAX1(0._DBL,(1._DBL-0.1_DBL*DABS(ud*dd/gamma_d))**5._DBL)+&
    &DMAX1(0._DBL,-ud*delta_y)
    AS=gamma_s*delta_x/ds*DMAX1(0._DBL,(1._DBL-0.1_DBL*DABS(vs*ds/gamma_s))**5._DBL)+&
    &DMAX1(0._DBL, vs*delta_x)
    AN=gamma_n*delta_x/dn*DMAX1(0._DBL,(1._DBL-0.1_DBL*DABS(vn*dn/gamma_n))**5._DBL)+&
    &DMAX1(0._DBL,-vn*delta_x)
    AC=AI+AD+AS+AN+delta_x*delta_y/dt_o
    Resu_o(i,j)=AI*u_o(i-1,j)+AD*u_o(i+1,j)+AS*u_o(i,j-1)+AN*u_o(i,j+1)-AC*u_o(i,j)+&
    &delta_x*delta_y*u_anto(i,j)/dt_o-delta_x*delta_y*Ri_o(i,j)*temp_int+(pres_o(i,j)-&
    &pres_o(i+1,j))*delta_y
    !**************************
  END DO
END DO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
DO i = nsolid+1, mi-nsolid
  DO j = nsolid+1, nj+1-nsolid
    IF(DABS(Resu_o(i,j))>1.e-5) THEN 
      res_fluido_uo = .FALSE.
    END IF
  END DO
END DO
!$OMP END PARALLEL DO

! !$OMP PARALLEL DO
! DO i = 1, mi
!   DO j = 1, nsolid
!     IF(DABS(u_o(i,j)) > 1.e-20_DBL .OR. DABS(u_o(i,nj+2-j)) > 1.e-20_DBL .OR. DABS(Resu_o(i,j)) > 1.e-7) THEN 
!       res_fluido_uo = .FALSE.
!     END IF
!   END DO
! END DO
! !$OMP END PARALLEL DO
! !$OMP PARALLEL DO
! DO i = 1, nsolid 
!   DO j = 1, nj+1 
!     IF(DABS(u_o(i,j)) > 1.e-20_DBL .OR. DABS(u_o(mi+1-i,j) ) > 1.e-20_DBL .OR. DABS(Resu_o(i,j)) > 1.e-7) THEN 
!       res_fluido_uo = .FALSE.
!     END IF
!   END DO
! END DO
! !$OMP END PARALLEL DO
END SUBROUTINE
