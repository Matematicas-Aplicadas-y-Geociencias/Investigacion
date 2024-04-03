SUBROUTINE temperatura(xo,yo,zo,fexu,feyv,fezw,d_xuo,d_yvo,d_zwo,u_o,v_o,w_o,temp_o,temp_anto,gamma_to,dt_o,dtemp_o,i_oo,i_1o,rel_to,j_oo,j_1o)
USE constantes
! USE mkl95_LAPACK
IMPLICIT NONE
INTEGER :: i,j,k,info
INTEGER, INTENT(in) :: i_oo,i_1o,j_oo,j_1o
!*********************************
! Variables del flujo e inc'ognita
REAL(kind=DBL), DIMENSION(mi+1,nj+1,lk+1), INTENT(inout) :: temp_o
REAL(kind=DBL), DIMENSION(mi+1,nj+1,lk+1), INTENT(out)   :: dtemp_o
REAL(kind=DBL), DIMENSION(mi+1,nj+1,lk+1), INTENT(in)    :: temp_anto,gamma_to
REAL(kind=DBL), DIMENSION(mi,nj+1,lk+1),   INTENT(in)    :: u_o
REAL(kind=DBL), DIMENSION(mi+1,nj,lk+1),   INTENT(in)    :: v_o
REAL(kind=DBL), DIMENSION(mi+1,nj+1,lk),   INTENT(in)    :: w_o
REAL(kind=DBL), DIMENSION(mi+1,nj+1,lk+1)                :: ftemp
!***************************
!Variables de la tridiagonal
REAL(kind=DBL), DIMENSION(mi+1) :: APi,Rxi
REAL(kind=DBL), DIMENSION(mi)   :: AW,AE
REAL(kind=DBL), DIMENSION(nj+1) :: APj,Ryj
REAL(kind=DBL), DIMENSION(nj)   :: AS,AN
REAL(kind=DBL), DIMENSION(lk+1) :: APk,Rzk
REAL(kind=DBL), DIMENSION(lk)   :: AT,AB
!***************************
!Variables de interpolaci'on
REAL(kind=DBL) :: uw,ue,vs,vn,wb,wt
REAL(kind=DBL) :: dw,de,ds,dn,db,dtop
REAL(kind=DBL) :: gamma_w,gamma_e,gamma_s,gamma_n,gamma_b,gamma_top
REAL(kind=DBL) :: delta_x,delta_y,delta_z
!***************************************************************
!Variables de la malla, volumen de control, incremento de tiempo
REAL(kind=DBL), DIMENSION(mi+1), INTENT(in) :: xo,fexu
REAL(kind=DBL), DIMENSION(nj+1), INTENT(in) :: yo,feyv
REAL(kind=DBL), DIMENSION(lk+1), INTENT(in) :: zo,fezw
REAL(kind=DBL), DIMENSION(mi-1), INTENT(in) :: d_xuo
REAL(kind=DBL), DIMENSION(nj-1), INTENT(in) :: d_yvo
REAL(kind=DBL), DIMENSION(lk-1), INTENT(in) :: d_zwo
REAL(kind=DBL), INTENT(in)                  :: dt_o,rel_to
!********************
!Variables auxiliares
REAL(kind=DBL) alpha,beta,kappa,sigma
!****************************
!auxiliar para calcular dtemp
ftemp = temp_o
!*******************
!TDMA en direccion x
!$OMP  PARALLEL DO DEFAULT(NONE)&
!$OMP& PRIVATE(AW,APi,AE,Rxi,uw,ue,vs,vn,wb,wt,dw,de,ds,dn,db,dtop,gamma_w,gamma_e,gamma_s,gamma_n,gamma_b,gamma_top,delta_x,delta_y,delta_z,alpha,beta,kappa,sigma)&
!$OMP& SHARED(u_o,v_o,w_o,temp_o,temp_anto,d_xuo,d_yvo,d_zwo,xo,yo,zo,fexu,feyv,fezw,dt_o,gamma_to,rel_to,i_oo,i_1o,j_oo,j_1o)
DO k = 2, lk
  DO j = 2, nj
    !*************************
    !Llenado de la matriz en x
    DO i = 2, mi
      !**************************
      !Interpolaciones necesarias
      ue = u_o(i,j,k)
      uw = u_o(i-1,j,k)
      vn = v_o(i,j,k)
      vs = v_o(i,j-1,k)
      wt = w_o(i,j,k)
      wb = w_o(i,j,k-1)
      de = xo(i+1)-xo(i)
      dw = xo(i)-xo(i-1)
      dn = yo(j+1)-yo(j)
      ds = yo(j)-yo(j-1)
      dtop = zo(k+1)-zo(k)
      db = zo(k)-zo(k-1)
      gamma_e = 1._DBL/((1._DBL-fexu(i+1)) / gamma_to(i,j,k) + fexu(i+1) / gamma_to(i+1,j,k))
      gamma_w = 1._DBL/((1._DBL-fexu(i)) / gamma_to(i-1,j,k) + fexu(i) / gamma_to(i,j,k))
      gamma_n = 1._DBL/((1._DBL-feyv(j+1)) / gamma_to(i,j,k) + feyv(j+1) / gamma_to(i,j+1,k))
      gamma_s = 1._DBL/((1._DBL-feyv(j)) / gamma_to(i,j-1,k) + feyv(j) / gamma_to(i,j,k))
      gamma_top = 1._DBL/((1._DBL-fezw(k+1)) / gamma_to(i,j,k) + fezw(k+1) / gamma_to(i,j,k+1))
      gamma_b   = 1._DBL/((1._DBL-fezw(k)) / gamma_to(i,j,k-1) + fezw(k) / gamma_to(i,j,k))
      delta_x = d_xuo(i-1)
      delta_y = d_yvo(j-1)
      delta_z = d_zwo(k-1)
      !*******************
      AE(i)   =-(gamma_e*delta_y*delta_z/de*DMAX1(cero,(1._DBL-0.1_DBL*DABS(ue*de/gamma_e))**5._DBL)+&
      &DMAX1(cero,-ue*delta_y*delta_z))
      AW(i-1) =-(gamma_w*delta_y*delta_z/dw*DMAX1(cero,(1._DBL-0.1_DBL*DABS(uw*dw/gamma_w))**5._DBL)+&
      &DMAX1(cero,uw*delta_y*delta_z))
      alpha   =  gamma_n*delta_x*delta_z/dn*DMAX1(cero,(1._DBL-0.1_DBL*DABS(vn*dn/gamma_n))**5._DBL)+&
      &DMAX1(cero,-vn*delta_x*delta_z)
      beta    =  gamma_s*delta_x*delta_z/ds*DMAX1(cero,(1._DBL-0.1_DBL*DABS(vs*ds/gamma_s))**5._DBL)+&
      &DMAX1(cero,vs*delta_x*delta_z)
      kappa   =  gamma_top*delta_x*delta_y/dtop*DMAX1(cero,(1._DBL-0.1_DBL*DABS(wt*dtop/gamma_top))**5._DBL)+&
      &DMAX1(cero,-wt*delta_x*delta_y)
      sigma   =  gamma_b*delta_x*delta_y/db*DMAX1(cero,(1._DBL-0.1_DBL*DABS(wb*db/gamma_b))**5._DBL)+&
      &DMAX1(cero,wb*delta_x*delta_y)
      APi(i)  =(-AE(i)-AW(i-1)+alpha+beta+kappa+sigma+delta_x*delta_y*delta_z/dt_o) / rel_to
      Rxi(i)  = alpha*temp_o(i,j+1,k)+beta*temp_o(i,j-1,k)+kappa*temp_o(i,j,k+1)+sigma*temp_o(i,j,k-1)+delta_x*delta_y*delta_z*temp_anto(i,j,k)/dt_o+&
      &APi(i)*(1._DBL-rel_to)*temp_o(i,j,k)
      !************************************
    END DO
    !***********************
    !Condiciones de frontera
    IF (k>=i_oo .AND. k<=i_1o .and. j>=j_oo .AND. j<=j_1o) THEN
      APi(1) = 1._DBL
      AE(1)  = cero
      Rxi(1) = 1._DBL
    ELSE
      APi(1) =-1._DBL
      AE(1)  = 1._DBL
      Rxi(1) = temp_o(1,j,k)*delta_x !cero
    END IF
    !*****************
    IF (k>=i_oo .and. k<=i_1o .and. j>=j_oo .AND. j<=j_1o) THEN
      APi(mi+1) = 1._DBL
      AW(mi)    = 0._DBL !-1._DBL ! cero
      Rxi(mi+1) = 1.0_DBL !-temp_o(mi+1,j,k)*100.d0*delta_x !cero   ! 1._DBL
    ELSE
      APi(mi+1) = 1._DBL
      AW(mi)    =-1._DBL
      Rxi(mi+1) =-temp_o(mi+1,j,k)*delta_x !cero
    END IF 
    !***************************
    CALL tri(AW,APi,AE,Rxi,mi+1)
    DO i = 1, mi+1
      temp_o(i,j,k) = Rxi(i)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO
!********************
!TDMA en direccion  y
!$OMP  PARALLEL DO DEFAULT(NONE)&
!$OMP& PRIVATE(AS,APj,AN,Ryj,uw,ue,vs,vn,wb,wt,dw,de,ds,dn,db,dtop,gamma_w,gamma_e,gamma_s,gamma_n,gamma_b,gamma_top,delta_x,delta_y,delta_z,alpha,beta,kappa,sigma)&
!$OMP& SHARED(u_o,v_o,w_o,temp_o,temp_anto,gamma_to,d_xuo,d_yvo,d_zwo,xo,yo,zo,fexu,feyv,fezw,dt_o,i_oo,i_1o,rel_to)
DO k = 2, lk
  DO i = 2, mi
    !*************************
    !Llenado de la matriz en x
    DO j = 2, nj
      !**************************
      !Interpolaciones necesarias
      ue = u_o(i,j,k)
      uw = u_o(i-1,j,k)
      vn = v_o(i,j,k)
      vs = v_o(i,j-1,k)
      wt = w_o(i,j,k)
      wb = w_o(i,j,k-1)
      de = xo(i+1)-xo(i)
      dw = xo(i)-xo(i-1)
      dn = yo(j+1)-yo(j)
      ds = yo(j)-yo(j-1)
      dtop = zo(k+1)-zo(k)
      db = zo(k)-zo(k-1)
      gamma_e = 1._DBL/((1._DBL-fexu(i+1)) / gamma_to(i,j,k) + fexu(i+1) / gamma_to(i+1,j,k))
      gamma_w = 1._DBL/((1._DBL-fexu(i)) / gamma_to(i-1,j,k) + fexu(i) / gamma_to(i,j,k))
      gamma_n = 1._DBL/((1._DBL-feyv(j+1)) / gamma_to(i,j,k) + feyv(j+1) / gamma_to(i,j+1,k))
      gamma_s = 1._DBL/((1._DBL-feyv(j)) / gamma_to(i,j-1,k) + feyv(j) / gamma_to(i,j,k))
      gamma_top = 1._DBL/((1._DBL-fezw(k+1)) / gamma_to(i,j,k) + fezw(k+1) / gamma_to(i,j,k+1))
      gamma_b   = 1._DBL/((1._DBL-fezw(k)) / gamma_to(i,j,k-1) + fezw(k) / gamma_to(i,j,k))
      delta_x = d_xuo(i-1)
      delta_y = d_yvo(j-1)
      delta_z = d_zwo(k-1)
      !*******************
      alpha   =  gamma_e*delta_y*delta_z/de*DMAX1(cero,(1._DBL-0.1_DBL*DABS(ue*de/gamma_e))**5._DBL)+&
      &DMAX1(cero,-ue*delta_y*delta_z)
      beta    =  gamma_w*delta_y*delta_z/dw*DMAX1(cero,(1._DBL-0.1_DBL*DABS(uw*dw/gamma_w))**5._DBL)+&
      &DMAX1(cero,uw*delta_y*delta_z)
      AN(j)   =-(gamma_n*delta_x*delta_z/dn*DMAX1(cero,(1._DBL-0.1_DBL*DABS(vn*dn/gamma_n))**5._DBL)+&
      &DMAX1(cero,-vn*delta_x*delta_z))
      AS(j-1) =-(gamma_s*delta_x*delta_z/ds*DMAX1(cero,(1._DBL-0.1_DBL*DABS(vs*ds/gamma_s))**5._DBL)+&
      &DMAX1(cero,vs*delta_x*delta_z))
      kappa   =  gamma_top*delta_x*delta_y/dtop*DMAX1(cero,(1._DBL-0.1_DBL*DABS(wt*dtop/gamma_top))**5._DBL)+&
      &DMAX1(cero,-wt*delta_x*delta_y)
      sigma   =  gamma_b*delta_x*delta_y/db*DMAX1(cero,(1._DBL-0.1_DBL*DABS(wb*db/gamma_b))**5._DBL)+&
      &DMAX1(cero,wb*delta_x*delta_y)
      APj(j)  = (alpha+beta-AN(j)-AS(j-1)+kappa+sigma+delta_x*delta_y*delta_z/dt_o) / rel_to
      Ryj(j)  = alpha*temp_o(i+1,j,k)+beta*temp_o(i-1,j,k)+kappa*temp_o(i,j,k+1)+sigma*temp_o(i,j,k-1)+delta_x*delta_y*delta_z*temp_anto(i,j,k)/dt_o+&
      &APj(j)*(1._DBL-rel_to)*temp_o(i,j,k)
      !************************************
    END DO
    !***********************
    !Condiciones de frontera
    ! IF (j>=j_oo .AND. j<=j_1o) THEN
    !   APj(1) = 1._DBL
    !   AN(1)  = cero
    !   Ryj(1) = 1._DBL
    !   ! ------------
    !   APj(nj+1) = 1._DBL
    !   AN(nj)    = cero
    !   Ryj(nj+1) = 1._DBL
    ! ELSE
    !   APj(1) =-1._DBL
    !   AN(1)  = 1._DBL
    !   Ryj(1) = temp_o(i,1,k)*delta_y !cero
    !   ! ------------
    !   APj(nj+1) = 1._DBL
    !   AS(nj)    =-1._DBL
    !   Ryj(nj+1) = temp_o(i,nj+1,k)*delta_y !cero
    ! END IF
    ! !*****************
    ! ! IF (k>=i_oo .and. k<=i_1o) THEN
    ! !   APi(mi+1) = 1._DBL
    ! !   AW(mi)    = 0._DBL !-1._DBL ! cero
    ! !   Rxi(mi+1) = 1.0_DBL !-temp_o(mi+1,j,k)*100.d0*delta_x !cero   ! 1._DBL
    ! ! ELSE
    ! !   APi(mi+1) = 1._DBL
    ! !   AW(mi)    =-1._DBL
    ! !   Rxi(mi+1) =-temp_o(mi+1,j,k)*100.d0*delta_x !cero
    ! ! END IF
    !***********************
    !Condiciones de frontera
    APj(1) =-1._DBL
    AN(1)  = 1._DBL
    Ryj(1) = temp_o(i,1,k)*1.d0*delta_y !cero
    !*****************
    AS(nj)    =-1._DBL
    APj(nj+1) = 1._DBL
    Ryj(nj+1) =-temp_o(i,nj+1,k)*1.d0*delta_y !cero
    !***************************
    CALL tri(AS,APj,AN,Ryj,nj+1)
    DO j = 1, nj+1
      temp_o(i,j,k) = Ryj(j)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO
!********************
!TDMA en direccion  z
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(AB,APk,AT,Rzk,uw,ue,vs,vn,wb,wt,dw,de,ds,dn,db,dtop,gamma_w,gamma_e,gamma_s,gamma_n,gamma_b,gamma_top,delta_x,delta_y,delta_z,alpha,beta,kappa,sigma) SHARED(u_o,v_o,w_o,temp_o,temp_anto,gamma_to,d_xuo,d_yvo,d_zwo,xo,yo,zo,fexu,feyv,fezw,dt_o,i_oo,i_1o,rel_to)
DO j = 2, nj
  DO i = 2, mi
    !*************************
    !Llenado de la matriz en x
    DO k = 2, lk
      !**************************
      !Interpolaciones necesarias
      ue = u_o(i,j,k)
      uw = u_o(i-1,j,k)
      vn = v_o(i,j,k)
      vs = v_o(i,j-1,k)
      wt = w_o(i,j,k)
      wb = w_o(i,j,k-1)
      de = xo(i+1)-xo(i)
      dw = xo(i)-xo(i-1)
      dn = yo(j+1)-yo(j)
      ds = yo(j)-yo(j-1)
      dtop = zo(k+1)-zo(k)
      db = zo(k)-zo(k-1)
      gamma_e = 1._DBL/((1._DBL-fexu(i+1)) / gamma_to(i,j,k) + fexu(i+1) / gamma_to(i+1,j,k))
      gamma_w = 1._DBL/((1._DBL-fexu(i)) / gamma_to(i-1,j,k) + fexu(i) / gamma_to(i,j,k))
      gamma_n = 1._DBL/((1._DBL-feyv(j+1)) / gamma_to(i,j,k) + feyv(j+1) / gamma_to(i,j+1,k))
      gamma_s = 1._DBL/((1._DBL-feyv(j)) / gamma_to(i,j-1,k) + feyv(j) / gamma_to(i,j,k))
      gamma_top = 1._DBL/((1._DBL-fezw(k+1)) / gamma_to(i,j,k) + fezw(k+1) / gamma_to(i,j,k+1))
      gamma_b   = 1._DBL/((1._DBL-fezw(k)) / gamma_to(i,j,k-1) + fezw(k) / gamma_to(i,j,k))
      delta_x = d_xuo(i-1)
      delta_y = d_yvo(j-1)
      delta_z = d_zwo(k-1)
      !*******************
      alpha   =  gamma_e*delta_y*delta_z/de*DMAX1(cero,(1._DBL-0.1_DBL*DABS(ue*de/gamma_e))**5._DBL)+&
      &DMAX1(cero,-ue*delta_y*delta_z)
      beta    =  gamma_w*delta_y*delta_z/dw*DMAX1(cero,(1._DBL-0.1_DBL*DABS(uw*dw/gamma_w))**5._DBL)+&
      &DMAX1(cero,uw*delta_y*delta_z)
      kappa   =  gamma_n*delta_x*delta_z/dn*DMAX1(cero,(1._DBL-0.1_DBL*DABS(vn*dn/gamma_n))**5._DBL)+&
      &DMAX1(cero,-vn*delta_x*delta_z)
      sigma   =  gamma_s*delta_x*delta_z/ds*DMAX1(cero,(1._DBL-0.1_DBL*DABS(vs*ds/gamma_s))**5._DBL)+&
      &DMAX1(cero,vs*delta_x*delta_z)
      AT(k)   =-(gamma_top*delta_x*delta_y/dtop*DMAX1(cero,(1._DBL-0.1_DBL*DABS(wt*dtop/gamma_top))**5._DBL)+&
      &DMAX1(cero,-wt*delta_x*delta_y))
      AB(k-1) =-(gamma_b*delta_x*delta_y/db*DMAX1(cero,(1._DBL-0.1_DBL*DABS(wb*db/gamma_b))**5._DBL)+&
      &DMAX1(cero,wb*delta_x*delta_y))
      APk(k)  = (alpha+beta+kappa+sigma-AT(k)-AB(k-1)+delta_x*delta_y*delta_z/dt_o) / rel_to
      Rzk(k)  = alpha*temp_o(i+1,j,k)+beta*temp_o(i-1,j,k)+kappa*temp_o(i,j+1,k)+sigma*temp_o(i,j-1,k)+delta_x*delta_y*delta_z*temp_anto(i,j,k)/dt_o+&
      &APk(k)*(1._DBL-rel_to)*temp_o(i,j,k)
      !************************************
    END DO
    !***********************
    !Condiciones de frontera
    APk(1) =-1._DBL
    AT(1)  = 1._DBL
    Rzk(1) = cero
    !*****************
    AB(lk)    = cero !-1._DBL
    APk(lk+1) = 1._DBL
    Rzk(lk+1) = cero
    !***************************
    CALL tri(AB,APk,AT,Rzk,lk+1)
    DO k = 1, lk+1
      temp_o(i,j,k) = Rzk(k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO
!********************
! WHERE(temp_o /= cero)
!   dtemp_o = (ftemp-temp_o)/temp_o
! ELSEWHERE
dtemp_o = ftemp-temp_o
! ENDWHERE
END SUBROUTINE
