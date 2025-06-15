PROGRAM SIMPLE
USE constantes, only : itermax, cero
! USE mkl95_LAPACK
!
! Variables de la malla, volumen de control y factores de interpolaci\'on
!
use malla, only : mi, nj, lk, DBL
use malla, only : mic, njc, lkc
use malla, only : xu, yv, zw, xp, yp, zp
use malla, only : deltaxp, deltayp, deltazp
use malla, only : deltaxu, deltaxv, deltaxw
use malla, only : deltayu, deltayv, deltayw
use malla, only : deltazu, deltazv, deltazw
use malla, only : fexp, feyp, fezp
use malla, only : fexu, feyv, fezw
use malla, only : form21, form26
use malla, only : lectura_mallas_escalonadas
use malla, only : indexu, indexp
use malla, only : indeyu
use malla, only : indeyv, indeyp
use malla, only : indezw, indezp
!
use ec_continuidad, only : pres, corr_pres
use ec_continuidad, only : dcorr_pres, fcorr_pres
use ec_continuidad, only : b_o
use ec_continuidad, only : ensambla_corr_pres_x
use ec_continuidad, only : ensambla_corr_pres_y
use ec_continuidad, only : ensambla_corr_pres_z
!
use ec_momento, only : u, u_ant, du, au, Resu, fu
use ec_momento, only : v, v_ant, dv, av, Resv, fv
use ec_momento, only : w, w_ant, dw, aw, Resw, fw
use ec_momento, only : fuente_con_u, fuente_lin_u
use ec_momento, only : fuente_con_v, fuente_lin_v
use ec_momento, only : fuente_con_w, fuente_lin_w
use ec_momento, only : gamma_momen, Ri
use ec_momento, only : ensambla_velu_x
use ec_momento, only : ensambla_velu_y
use ec_momento, only : ensambla_velu_z
!
use solucionador, only : tridiagonal
!
IMPLICIT NONE
INCLUDE 'omp_lib.h'
!
! Variables de modularizaci\'on
!
!
! Coeficientes para las matrices 
!
real(kind=DBL), dimension(mi+1,nj+1,lk+1) :: AA, BB, CC, RR
real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)) :: a1, b1, c1, r1
!
! --------------------------------------------------------------------------------
!
integer :: ii, jj, kk
INTEGER :: i,j,k,tt,kl,l,itera_total,itera,itera_inicial,i_o,i_1,j_o,j_1,paq_itera
INTEGER :: millar,centena,decena,unidad,decima,id,nthreads
!*******************************************
! Variables del flujo,entropia,nusselt e inc'ognitas,residuos,relajaci'on y convergencia
REAL(kind=DBL), DIMENSION(mi,nj+1,lk+1)   :: gamma_u
REAL(kind=DBL), DIMENSION(mi+1,nj,lk+1)   :: gamma_v
REAL(kind=DBL), DIMENSION(mi+1,nj+1,lk)   :: gamma_w
REAL(kind=DBL), DIMENSION(mi+1,nj+1,lk+1) :: temp,temp_ant,dtemp,Restemp,gamma_t
!,pres,corr_pres,dcorr_pres
REAL(kind=DBL), DIMENSION(mi+1,nj+1,lk+1) :: entropia_calor,entropia_viscosa,entropia,uf,vf,wf
REAL(kind=DBL) :: temp_med,nusselt0,nusselt1,entropia_int,temp_int,gamma_s
REAL(kind=DBL) :: conv_u,conv_p,conv_t,conv_resi,conv_paso,rel_pres,rel_vel,rel_tem
!********************************************
!Variables de la malla, volumen de control, incremento de tiempo y nums Reynolds, Peclet
!Rayleigh, Richardson, Prandtl, valores de las constantes de difusividad
! REAL(kind=DBL), DIMENSION(mi)   :: xu
! REAL(kind=DBL), DIMENSION(nj)   :: yv
! REAL(kind=DBL), DIMENSION(lk)   :: zw
! REAL(kind=DBL), DIMENSION(mi+1) :: x,fexu
! REAL(kind=DBL), DIMENSION(nj+1) :: y,feyv
! REAL(kind=DBL), DIMENSION(lk+1) :: z,fezw
REAL(kind=DBL), DIMENSION(mi-1) :: d_xu,d2_xu
REAL(kind=DBL), DIMENSION(nj-1) :: d_yv,d2_yv
REAL(kind=DBL), DIMENSION(lk-1) :: d_zw,d2_zw
REAL(kind=DBL)   :: ao
REAL(kind=DBL)   :: tiempo,tiempo_inicial,dt,Pe,Re,Ra,Pr,Rin
REAL(kind=DBL)   :: a_ent,lambda_ent
CHARACTER(len=1) :: dec,un,de,ce,m
! CHARACTER(len=3) :: mic,njc,lkc,
character(len=3) :: Rac
CHARACTER(len=32):: entrada_u,entrada_v,entrada_w,entrada_tp,entrada_xyz
character(len=40):: archivo=repeat(' ',40)
LOGICAL          :: res_fluido_u
!****************************************
!declaraci´on de variable DBL
REAL(kind=DBL) :: var2=cero
!*******************************************
!Se muestra cu'antos procesadores hay en uso
!$OMP PARALLEL private(id)
id = omp_get_thread_num()
WRITE(*,*) 'Este es el thread no. ', id
!$OMP BARRIER
IF(id == 0)THEN
  nthreads = omp_get_num_threads()
  WRITE(*,*) 'Se usan ', nthreads, ' threads'
END IF
!$OMP END PARALLEL
!*************************************
! Par'ametros para Convecci'on natural
OPEN(unit=10,file='parametros.dat')
  READ (10,*) Ra          ! n'umero de Rayleigh
  READ (10,*) Rac         ! caracter de Ra
  READ (10,*) Pr          ! n'umero de Prandtl
  READ (10,*) dt          ! incremento de tiempo
  READ (10,*) paq_itera   ! paquete de iteraciones
  READ (10,*) Rin         ! n'umero de Richardson
  READ (10,*) rel_pres    ! relajaci'on de la presi'on
  READ (10,*) rel_vel     ! relajaci'on de la velocidad
  READ (10,*) rel_tem     ! relajaci'on de la temperatura
  READ (10,*) conv_u      ! convergencia de la velocidad
  READ (10,*) conv_t      ! convergencia de la temperatura
  READ (10,*) conv_p      ! convergencia de la presi'on
  READ (10,*) conv_resi   ! convergencia del residuo
  READ (10,*) conv_paso   ! convergencia del paso de tiempo
  READ (10,*) entrada_u   ! archivo de entrada para u
  READ (10,*) entrada_v   ! archivo de entrada para v
  READ (10,*) entrada_w   ! archivo de entrada para w
  READ (10,*) entrada_tp  ! archivo de entrada para t y p
CLOSE(unit=10)
IF(mi < 100)THEN
  WRITE(mic,170) int(mi)
  mic = '0'//mic
ELSE
  WRITE(mic,160) int(mi)
ENDIF
IF(nj < 100)THEN
  WRITE(njc,170) int(nj);170 format(I2)
  njc = '0'//njc
ELSE
  WRITE(njc,160) int(nj);160 format(I3)
ENDIF
IF(lk < 100)THEN
  WRITE(lkc,170) int(lk)
  lkc = '0'//lkc
ELSE
  WRITE(lkc,160) int(lk)
ENDIF
gamma_s = 600._DBL*sqrt(1._DBL/(Pr*Ra))
gamma_t = 1._DBL/(Pr*Ra) !sqrt(1._DBL/(Pr*Ra))
gamma_u = 1._DBL/Ra      !sqrt(Pr/Ra)
gamma_momen = 1._DBL/Ra
fuente_con_u = 0._DBL
fuente_lin_u = 0._DBL
gamma_v = 1._DBL/Ra      !sqrt(Pr/Ra)
gamma_w = 1._DBL/Ra      !sqrt(Pr/Ra)
Ri      = Rin
!***************************
OPEN(unit=11,file=entrada_u)
READ(11,*)i_o,i_1,itera_inicial,ao
DO k = 1, lk+1
  DO j = 1, nj+1
    DO i = 1, mi
      READ(11,24) xu(i),yp(j),zp(k),u_ant(i,j,k)
    END DO
  END DO
END DO
CLOSE(unit=11)
!***************************
OPEN(unit=12,file=entrada_v)
READ(12,*)j_o,j_1,itera_inicial,ao
DO k =1, lk+1
  DO j = 1, nj
    DO i = 1, mi+1
      READ(12,24) xp(i),yv(j),zp(k),v_ant(i,j,k)
    END DO
  END DO
END DO
CLOSE(unit=12)
!***************************
OPEN(unit=13,file=entrada_w)
READ(13,*)i_o,i_1,itera_inicial,ao
DO k =1, lk
  DO j = 1, nj+1
    DO i = 1, mi+1
      READ(13,24) xp(i),yp(j),zw(k),w_ant(i,j,k)
    END DO
  END DO
END DO
CLOSE(unit=13)
!****************************
OPEN(unit=14,file=entrada_tp)
READ(14,*)i_o,i_1,itera_inicial,ao
DO k = 1, lk+1
  DO j = 1, nj+1
    DO i = 1, mi+1
      READ(14,25) xp(i),yp(j),zp(k),temp_ant(i,j,k),pres(i,j,k)
    END DO
  END DO
END DO
CLOSE(unit=14)
! temp_ant(1,1)    = 1._DBL
! temp_ant(mi+1,1) = 1._DBL
!******************************
!
!        incrementos
!
d_xu(1) = xu(2)-xu(1)
DO i = 2, mi-1
  d_xu(i)  = xu(i+1)-xu(i)
  d2_xu(i) = xu(i+1)-xu(i-1)
END DO
d_yv(1) = yv(2)-yv(1)
DO j = 2, nj-1
  d_yv(j)  = yv(j+1)-yv(j)
  d2_yv(j) = yv(j+1)-yv(j-1)
END DO
d_zw(1) = zw(2)-zw(1)
DO k = 2, lk-1
  d_zw(k)  = zw(k+1)-zw(k)
  d2_zw(k) = zw(k+1)-zw(k-1)
END DO
!
!------------------------------------------------
!
! Lectura de mallas escalonadas
!
write(*,*) "Inicia lectura de mallas"
call lectura_mallas_escalonadas(entrada_u,entrada_v,entrada_w,&
     &entrada_tp,entrada_xyz,&
     &u_ant,v_ant,w_ant,pres,temp_ant,&
     &xp,yp,zp,xu,yv,zw,&
     &deltaxp,deltayp,deltazp,&
     &deltaxu,deltayu,deltazu,&
     &deltaxv,deltayv,deltazv,&
     &deltaxw,deltayw,deltazw,&
     &fexp,feyp,fezp,fexu,feyv,fezw,&
     &ao,itera_inicial)
write(*,*) "Finaliza lectura de mallas"
write(*,*) "--------------------------"
!*****************
!valores iniciales
! DO k= 1, lk+1
!   pres(:,:,k) = 12._DBL/Ra*z(k)
! END DO
tiempo_inicial = itera_inicial*dt
! u_ant = cero
! v_ant = cero
! w_ant =-1.0_DBL
! temp_ant = cero !temp_ant
corr_pres = cero
au    = 1._DBL
av    = 1._DBL
aw    = 1._DBL
! pres  = cero
b_o   = cero
itera = 0
entropia = cero
!************************************************
!escribe las caracter´isticas de las variable DBL
WRITE(*,100) 'Doble',KIND(var2),PRECISION(var2),RANGE(var2)
100 FORMAT(1X,A,': kind= ',I2,', Precision= ',I2,' Rango= ',I3)
WRITE(*,*)' '
!escribe informaci'on de los parametros usados
WRITE(*,101) Ra,Pr,rel_pres,rel_vel
WRITE(*,102) itera_inicial,mi,nj
WRITE(*,106) lambda_ent,a_ent
WRITE(*,*)' '
101 FORMAT(1X,'Ra=',F12.3,', Pr=',F8.3', rel_pres=',F8.3', rel_vel=',F8.3)
102 FORMAT(1X,'Iteracion inicial=',I7,', mi=',I3,', nj=',I3)
106 FORMAT(1X,'No. de Eckert=',F13.10,', a_ent=',F15.3)
!*********************************************************
DO l=1,itermax/paq_itera   !inicio del repetidor principal
   DO kl=1,paq_itera          !inicio del paquete iteraciones
      ALGORITMO_SIMPLE: DO       !inicio del algoritmo SIMPLE
         DO tt= 1, 100
            !
            !----------------------------------------------------------
            !----------------------------------------------------------
            !
            !      Se resuelve la ecuaci'on de momento
            !
            !----------------------------------------------------------
            !----------------------------------------------------------
            !
            !$acc parallel loop gang collapse(2) !async(stream1)
            inicializacion_fu: do kk = 1, lk+1
               do jj = 1, nj+1
                  do ii = 1, mi
                     fu(ii,jj,kk) = u(ii,jj,kk)
                  end do
               end do
            end do inicializacion_fu
            !
            !----------------------------------------
            !
            ! Se ensamblan las matrices tridiagonales
            ! en la direcci'on de y
            !
            !$acc parallel loop gang !async(stream2)
            ensa_velu_dir_y: do kk = 2, lk
               do ii = 2, mi-1
                  do jj = 2, nj
                     call ensambla_velu_y(&
                          &deltaxu,&
                          &deltayu,&
                          &deltazu,&
                          &deltaxp,&
                          &deltayv,&
                          &deltazw,&
                          &fexp,&
                          &feyp,&
                          &fezp,&
                          &fexu,&
                          &gamma_momen,&
                          &u,&
                          &u_ant,&
                          &v,&
                          &w,&
                          &temp,&
                          &pres,&
                          &fuente_con_u,&
                          &fuente_lin_u,&
                          &Ri,&
                          &dt,&
                          &rel_vel,&
                          &a1,b1,c1,r1,&
                          &jj,ii,kk&
                          &)
                  end do
               end do
            end do ensa_velu_dir_y
            !
            !-------------------------
            !
            ! Condiciones de frontera direcci\'on y
            !
            !$acc parallel loop vector !async(stream1)
            cond_fron_u_direc_y: do kk = 2, lk
               do ii = 2, mi-1
                  !***********************
                  !Condiciones de frontera
                  a1(indeyu(1,ii,kk))    = 0.0_DBL
                  b1(indeyu(1,ii,kk))    = 1.0_DBL
                  c1(indeyu(1,ii,kk))    = 0.0_DBL
                  r1(indeyu(1,ii,kk))    = 0.0_DBL
                  au(ii,1,kk)            = 1.0e40_DBL
                  !
                  a1(indeyu(nj+1,ii,kk)) = 0.0_DBL
                  b1(indeyu(nj+1,ii,kk)) = 1.0_DBL
                  c1(indeyu(nj+1,ii,kk)) = 0.0_DBL
                  r1(indeyu(nj+1,ii,kk)) = 0.0_DBL
                  au(ii,nj+1,kk)         = 1.0e40_DBL
                  !
               end do
            end do cond_fron_u_direc_y
            !
            !----------------------------------------------
            !
            ! Soluci\'on de la ec. de momento para u
            ! en direcci\'on y
            !
            !$acc parallel loop gang async(stream1) wait(stream2)
            sol_u_dir_y: do kk = 2, lk
               do ii = 2, mi-1
                  !
                  call tridiagonal(&
                       &a1(indeyu(1,ii,kk):indeyu(nj+1,ii,kk)),&
                       &b1(indeyu(1,ii,kk):indeyu(nj+1,ii,kk)),&
                       &c1(indeyu(1,ii,kk):indeyu(nj+1,ii,kk)),&
                       &r1(indeyu(1,ii,kk):indeyu(nj+1,ii,kk)),&
                       &nj+1)
                  !
               end do
            end do sol_u_dir_y
            !
            !----------------------------------
            !
            ! Actualizaci\'on de la velocidad u
            ! en direcci\'on y
            !
            !$acc parallel loop gang collapse(2) !async(stream2) wait(stream1)
            do kk = 2, lk
               do ii = 2, mi-1
                  do jj = 1, nj+1
                     u(ii,jj,kk) = r1(indeyu(jj,ii,kk))
                  end do
               end do
            end do
            !
            !
            !----------------------------------------
            !
            ! Se ensamblan las matrices tridiagonales
            ! en la direcci'on de z
            !
            !$acc parallel loop gang !async(stream2)
            ensa_velu_dir_z: do ii = 2, mi-1
               do jj = 2, nj
                  do kk = 2, lk
                     call ensambla_velu_z(&
                          &deltaxu,&
                          &deltayu,&
                          &deltazu,&
                          &deltaxp,&
                          &deltayv,&
                          &deltazw,&
                          &fexp,&
                          &feyp,&
                          &fezp,&
                          &fexu,&
                          &gamma_momen,&
                          &u,&
                          &u_ant,&
                          &v,&
                          &w,&
                          &temp,&
                          &pres,&
                          &fuente_con_u,&
                          &fuente_lin_u,&
                          &Ri,&
                          &dt,&
                          &rel_vel,&
                          &a1,b1,c1,r1,&
                          &kk,jj,ii&
                          &)
                  end do
               end do
            end do ensa_velu_dir_z
            !
            !-------------------------
            !
            ! Condiciones de frontera direcci\'on z
            !
            !$acc parallel loop vector !async(stream1)
            cond_fron_u_direc_z: do ii = 2, mi-1
               do jj = 2, nj
                  !***********************
                  !Condiciones de frontera
                  a1(indezp(1,jj,ii))    = 0.0_DBL
                  b1(indezp(1,jj,ii))    = 1.0_DBL
                  c1(indezp(1,jj,ii))    = 0.0_DBL
                  r1(indezp(1,jj,ii))    = 1.0_DBL
                  au(ii,jj,1)            = 1.0e40_DBL
                  !
                  a1(indezp(lk+1,jj,ii)) = 0.0_DBL
                  b1(indezp(lk+1,jj,ii)) = 1.0_DBL
                  c1(indezp(lk+1,jj,ii)) = 0.0_DBL
                  r1(indezp(lk+1,jj,ii)) = 0.0_DBL
                  au(ii,jj,lk+1)         = 1.0e40_DBL
                  !
               end do
            end do cond_fron_u_direc_z
            !
            !----------------------------------------------
            !
            ! Soluci\'on de la ec. de momento para u
            ! en direcci\'on z
            !
            !$acc parallel loop gang async(stream1) wait(stream2)
            sol_u_dir_z: do ii = 2, mi-1
               do jj = 2, nj
                  !
                  call tridiagonal(&
                       &a1(indezp(1,jj,ii):indezp(lk+1,jj,ii)),&
                       &b1(indezp(1,jj,ii):indezp(lk+1,jj,ii)),&
                       &c1(indezp(1,jj,ii):indezp(lk+1,jj,ii)),&
                       &r1(indezp(1,jj,ii):indezp(lk+1,jj,ii)),&
                       &lk+1)
                  !
               end do
            end do sol_u_dir_z
            !
            !----------------------------------
            !
            ! Actualizaci\'on de la velocidad u
            ! en direcci\'on z
            !
            !$acc parallel loop gang collapse(2) !async(stream2) wait(stream1)
            do ii = 2, mi-1
               do jj = 2, nj
                  do kk = 1, lk+1
                     u(ii,jj,kk) = r1(indezp(kk,jj,ii))
                  end do
               end do
            end do
            !
            !----------------------------------------
            !
            ! Se ensamblan las matrices tridiagonales
            ! en la direcci'on de x
            !
            !$acc parallel loop gang !async(stream2)
            ensa_velu_dir_x: do kk = 2, lk
               do jj = 2, nj
                  do ii = 2, mi-1
                     call ensambla_velu_x(&
                          &deltaxu,&
                          &deltayu,&
                          &deltazu,&
                          &deltaxp,&
                          &deltayv,&
                          &deltazw,&
                          &fexp,&
                          &feyp,&
                          &fezp,&
                          &fexu,&
                          &gamma_momen,&
                          &u,&
                          &u_ant,&
                          &v,&
                          &w,&
                          &temp,&
                          &pres,&
                          &fuente_con_u,&
                          &fuente_lin_u,&
                          &Ri,&
                          &dt,&
                          &rel_vel,&
                          &a1,b1,c1,r1,&
                          &au,&
                          &ii,jj,kk&
                          &)
                  end do
               end do
            end do ensa_velu_dir_x
            !
            !-------------------------
            !
            ! Condiciones de frontera direcci\'on x
            !
            !$acc parallel loop vector !async(stream1)
            cond_fron_u_direc_x: do kk = 2, lk
               do jj = 2, nj
                  !***********************
                  !Condiciones de frontera
                  a1(indexu(1,jj,kk))    = 0.0_DBL
                  b1(indexu(1,jj,kk))    = 1.0_DBL
                  c1(indexu(1,jj,kk))    = 0.0_DBL
                  r1(indexu(1,jj,kk))    = 0.0_DBL
                  au(1,jj,kk)            = 1.0e40_DBL
                  !
                  a1(indexu(mi,jj,kk))   = 0.0_DBL
                  b1(indexu(mi,jj,kk))   = 1.0_DBL
                  c1(indexu(mi,jj,kk))   = 0.0_DBL
                  r1(indexu(mi,jj,kk))   = 0.0_DBL
                  au(mi,jj,kk)           = 1.0e40_DBL
                  !
               end do
            end do cond_fron_u_direc_x
            !
            !----------------------------------------------
            !
            ! Soluci\'on de la ec. de momento para u
            ! en direcci\'on x
            !
            !$acc parallel loop gang async(stream1) wait(stream2)
            sol_u_dir_x: do kk = 2, lk
               do jj = 2, nj
                  !
                  call tridiagonal(&
                       &a1(indexu(1,jj,kk):indexu(mi,jj,kk)),&
                       &b1(indexu(1,jj,kk):indexu(mi,jj,kk)),&
                       &c1(indexu(1,jj,kk):indexu(mi,jj,kk)),&
                       &r1(indexu(1,jj,kk):indexu(mi,jj,kk)),&
                       &mi)
                  !
               end do
            end do sol_u_dir_x
            !
            !----------------------------------
            !
            ! Actualizaci\'on de la velocidad u
            ! en direcci\'on x
            !
            !$acc parallel loop gang collapse(2) !async(stream2) wait(stream1)
            do kk = 2, lk
               do jj = 2, nj
                  do ii = 1, mi
                     u(ii,jj,kk) = r1(indexu(ii,jj,kk))
                  end do
               end do
            end do
            !
            !
            calcula_fu: do kk = 1, lk+1
               do jj = 1, nj+1
                  do ii = 1, mi
                     fu(ii,jj,kk) = fu(ii,jj,kk)-u(ii,jj,kk)
                  end do
               end do
            end do calcula_fu
            !
            ! CALL vel_u(xu,yp,zp,feyv,fezw,d_xu,d2_xu,d_yv,d_zw,u,u_ant,&
            !      &v,w,pres,gamma_u,dt,du,au,rel_vel)
            CALL vel_v(xp,yv,zp,fexu,fezw,d_yv,d2_yv,d_xu,d_zw,u,v,v_ant,&
                 &w,pres,gamma_v,dt,dv,av,rel_vel)
            CALL vel_w(xp,yp,zw,fexu,feyv,d_zw,d2_zw,d_xu,d_yv,w,w_ant,u,v&
                 &,pres,temp,gamma_w,Ri,dt,dw,aw,rel_vel)
            !****************************************
            !Criterio de convergencia de la velocidad
            IF(MAXVAL(DABS(fu))<conv_u.and.MAXVAL(DABS(dv))<conv_u.and.&
                 &MAXVAL(DABS(dw))<conv_u)EXIT
            WRITE(*,*) 'velocidad ',itera, MAXVAL(DABS(fu))
         END DO
    !
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    !
    !       Se calcula la correcci'on de la presi'on
    !
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    !
    !$acc parallel loop gang collapse(2) !async(stream2)
    inicializa_corrector_presion: do kk=1, lk+1
       do jj = 1, nj+1
          do ii = 1, mi+1
             corr_pres(ii,jj,kk) = 0.0_DBL
             fcorr_pres(ii,jj,kk)= 0.0_DBL
          end do
       end do
    end do inicializa_corrector_presion
    !
    correccion_presion: do tt = 1, 10
       !
       !$acc parallel loop gang collapse(2) ! async(stream1) wait(stream2)
       inicializa_fcorr_press: do kk=2, lk
          do jj=2, nj
             do ii = 2, mi
                fcorr_pres(ii,jj,kk) = corr_pres(ii,jj,kk)
             end do
          end do
       end do inicializa_fcorr_press
       !
       !------------------------------------------
       !
       ! Se ensambla la ecuaci\'on de correcci\'on
       ! de la presi\'on en direcci\'on x
       !
       !$acc parallel loop gang !async(stream2)
       ensa_corr_dir_x: do kk = 2, lk
          do jj = 2, nj
             !$acc loop vector
             do ii = 2, mi
                call ensambla_corr_pres_x(&
                     &deltaxp,&
                     &deltayp,&
                     &deltazp,&
                     &deltaxu,&
                     &deltayv,&
                     &deltazw,&
                     &u,v,w,&
                     &corr_pres,&
                     &rel_pres,&
                     &a1,b1,c1,r1,&
                     &au,av,aw,&
                     &ii,jj,kk)
             end do
          end do
       end do ensa_corr_dir_x
       !
       !---------------------------------------
       !
       ! Condiciones de frontera direcci\'on x
       !
       !$acc parallel loop vector !async(stream1)
       cond_fron_direc_x: do kk = 2, lk
          do jj = 2, nj
             !***********************
             !Condiciones de frontera
             a1(indexp(1,jj,kk))    = 0.0_DBL
             b1(indexp(1,jj,kk))    = 1.0_DBL
             c1(indexp(1,jj,kk))    = 0.0_DBL
             r1(indexp(1,jj,kk))    = 0.0_DBL
             !
             a1(indexp(mi+1,jj,kk)) = 0.0_DBL
             b1(indexp(mi+1,jj,kk)) = 1.0_DBL
             c1(indexp(mi+1,jj,kk)) = 0.0_DBL
             r1(indexp(mi+1,jj,kk)) = 0.0_DBL
             !
          end do
       end do cond_fron_direc_x
       !
       !----------------------------------------------
       !
       ! Soluci\'on de la correcci\'on de la presi\'on
       !
       !$acc parallel loop gang async(stream1) wait(stream2)
       sol_corr_dir_x: do kk = 2, lk
          do jj = 2,nj
             !
             call tridiagonal(&
                  &a1(indexp(1,jj,kk):indexp(mi+1,jj,kk)),&
                  &b1(indexp(1,jj,kk):indexp(mi+1,jj,kk)),&
                  &c1(indexp(1,jj,kk):indexp(mi+1,jj,kk)),&
                  &r1(indexp(1,jj,kk):indexp(mi+1,jj,kk)),&
                  &mi+1)
             ! corr_pres(1,jj,kk)    = r1(indexp(1,jj,kk))
             ! corr_pres(mi+1,jj,kk) = r1(indexp(mi+1,jj,kk))
             !
          end do
       end do sol_corr_dir_x
       !
       ! Actualizaci\'on del corrector de la presi\'on
       ! en direcci\'on x
       !
       do kk = 2, lk
          do jj = 2, nj
             do ii =1, mi+1
                corr_pres(ii,jj,kk) = r1(indexp(ii,jj,kk))
             end do
          end do
       end do
       !
       !------------------------------------------
       !
       ! Se ensambla la ecuaci\'on de correcci\'on
       ! de la presi\'on en direcci\'on x
       !
       !$acc parallel loop gang !async(stream2)
       ensa_corr_dir_y: do kk = 2, lk
          do ii = 2, mi
             !$acc loop vector
             do jj = 2, nj
                call ensambla_corr_pres_y(&
                     &deltaxp,&
                     &deltayp,&
                     &deltazp,&
                     &deltaxu,&
                     &deltayv,&
                     &deltazw,&
                     &u,v,w,&
                     &corr_pres,&
                     &rel_pres,&
                     &a1,b1,c1,r1,&
                     &au,av,aw,&
                     &jj,ii,kk)
             end do
          end do
       end do ensa_corr_dir_y
       !
       !--------------------------------------------
       !--------------------------------------------
       !
       ! Condiciones de frontera direcci\'on y
       !
       !$acc parallel loop vector !async(stream1)
       cond_fron_direc_y: do kk = 2, lk
          do ii = 2, mi
             !***********************
             !Condiciones de frontera
             a1(indeyp(1,ii,kk))    = 0.0_DBL
             b1(indeyp(1,ii,kk))    = 1.0_DBL
             c1(indeyp(1,ii,kk))    = 0.0_DBL
             !
             a1(indeyp(nj+1,ii,kk)) = 0.0_DBL
             b1(indeyp(nj+1,ii,kk)) = 1.0_DBL
             c1(indeyp(nj+1,ii,kk)) = 0.0_DBL
             !
          end do
       end do cond_fron_direc_y
       !
       !----------------------------------------------
       !
       ! Soluci\'on de la correcci\'on de la presi\'on
       ! en direcci\'on y
       !
       !$acc parallel loop gang async(stream1) wait(stream2)
       sol_corr_dir_y: do kk = 2, lk
          do ii = 2, mi
             !
             call tridiagonal(&
                  &a1(indeyp(1,ii,kk):indeyp(nj+1,ii,kk)),&
                  &b1(indeyp(1,ii,kk):indeyp(nj+1,ii,kk)),&
                  &c1(indeyp(1,ii,kk):indeyp(nj+1,ii,kk)),&
                  &r1(indeyp(1,ii,kk):indeyp(nj+1,ii,kk)),&
                  &nj+1)
             ! corr_pres(ii,1,kk)    = r1(indeyp(1,ii,kk))
             ! corr_pres(ii,nj+1,kk) = r1(indeyp(nj+1,ii,kk))
             !
          end do
       end do sol_corr_dir_y
       !
       ! Actualizaci\'on del corrector de la presi\'on
       ! en direcci\'on y
       !
       do kk = 2, lk
          do ii = 2, mi
             do jj =2, nj
                corr_pres(ii,jj,kk) = r1(indeyp(jj,ii,kk))
             end do
          end do
       end do
       !
       !------------------------------------------
       !
       ! Se ensambla la ecuaci\'on de correcci\'on
       ! de la presi\'on en direcci\'on x
       !
       !$acc parallel loop gang !async(stream2)
       ensa_corr_dir_z: do ii = 2, mi
          do jj = 2, nj
             !$acc loop vector
             do kk = 2, lk
                call ensambla_corr_pres_z(&
                     &deltaxp,&
                     &deltayp,&
                     &deltazp,&
                     &deltaxu,&
                     &deltayv,&
                     &deltazw,&
                     &u,v,w,&
                     &corr_pres,&
                     &b_o,&
                     &rel_pres,&
                     &a1,b1,c1,r1,&
                     &au,av,aw,&
                     &kk,jj,ii)
             end do
          end do
       end do ensa_corr_dir_z
       !
       !------------------------------------------------------------------------
       !------------------------------------------------------------------------
       !
       ! Condiciones de frontera direcci\'on z
       !
       !$acc parallel loop vector !async(stream1)
       cond_fron_direc_z: do ii = 2, mi
          do jj = 2, nj
             !***********************
             !Condiciones de frontera
             a1(indezp(1,jj,ii))    = 0.0_DBL
             b1(indezp(1,jj,ii))    = 1.0_DBL
             c1(indezp(1,jj,ii))    = 0.0_DBL
             !
             a1(indezp(lk+1,jj,ii)) = 0.0_DBL
             b1(indezp(lk+1,jj,ii)) = 1.0_DBL
             c1(indezp(lk+1,jj,ii)) = 0.0_DBL
             !
          end do
       end do cond_fron_direc_z
       !
       !----------------------------------------------
       !
       ! Soluci\'on de la correcci\'on de la presi\'on
       ! en direcci\'on y
       !
       !$acc parallel loop gang async(stream1) wait(stream2)
       sol_corr_dir_z: do ii = 2, mi
          do jj = 2, nj
             !
             call tridiagonal(&
                  &a1(indezp(1,jj,ii):indezp(lk+1,jj,ii)),&
                  &b1(indezp(1,jj,ii):indezp(lk+1,jj,ii)),&
                  &c1(indezp(1,jj,ii):indezp(lk+1,jj,ii)),&
                  &r1(indezp(1,jj,ii):indezp(lk+1,jj,ii)),&
                  &lk+1)
             ! corr_pres(ii,jj,1)    = r1(indezp(1,jj,ii))
             ! corr_pres(ii,jj,lk+1) = r1(indezp(lk+1,jj,ii))
             !
          end do
       end do sol_corr_dir_z
       !
       ! Actualizaci\'on del corrector de la presi\'on
       ! en direcci\'on y
       !
       do ii = 2, mi
          do jj = 2, nj
             do kk =1, lk+1
                corr_pres(ii,jj,kk) = r1(indezp(kk,jj,ii))
             end do
          end do
       end do
       !
       calcula_fcorr_press: do kk=2, lk
          do jj=2, nj
             do ii = 2, mi
                fcorr_pres(ii,jj,kk) = fcorr_pres(ii,jj,kk)-corr_pres(ii,jj,kk)
             end do
          end do
       end do calcula_fcorr_press
       !-----------------------------------------------------------------------------
       !-----------------------------------------------------------------------------
       !
       !****************************************************
       !critero de convergencia del corrector de la presi'on
       IF(MAXVAL(DABS(fcorr_pres))<conv_p)EXIT
       WRITE(*,*) 'corrector presion ', MAXVAL(DABS(fcorr_pres))
    END DO correccion_presion
    !*********************
    !se corrige la presion
    !$OMP PARALLEL DO
    DO k = 2, lk
      DO j = 2, nj
        DO i = 2, mi
          pres(i,j,k) = pres(i,j,k)+corr_pres(i,j,k)
        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO
    !*****************************
    !se actualizan las velocidades
    !$OMP PARALLEL DO
    DO k = 2, lk
      DO j = 2, nj
        DO i = 2, mi-1
          u(i,j,k) = u(i,j,k)+d_yv(j-1)*d_zw(k-1)*(corr_pres(i,j,k)-corr_pres(i+1,j,k))/au(i,j,k)
        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO
    DO k = 2, lk
      DO j = 2, nj-1
        DO i = 2, mi
          v(i,j,k) = v(i,j,k)+d_xu(i-1)*d_zw(k-1)*(corr_pres(i,j,k)-corr_pres(i,j+1,k))/av(i,j,k)
        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO
    DO k = 2, lk-1
      DO j = 2, nj
        DO i = 2, mi
          w(i,j,k) = w(i,j,k)+d_xu(i-1)*d_yv(j-1)*(corr_pres(i,j,k)-corr_pres(i,j,k+1))/aw(i,j,k)
        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO
    !*************************
    !se calcula la temperatura
    DO tt = 1, 3
      CALL temperatura(xp,yp,zp,fexu,feyv,fezw,d_xu,d_yv,d_zw,u,v,w,temp,temp_ant,gamma_t,dt,dtemp,i_o,i_1,rel_tem,j_o,j_1)
      !************************************
      !Criterio de convergencia temperatura
      IF(MAXVAL(DABS(dtemp))<conv_t)EXIT
      ! WRITE(*,*) 'temp', MAXVAL(DABS(dtemp))
    END DO
    !*******************************************
    !Criterio de convergencia del paso de tiempo
!     CALL residuou(res_fluido_u,xu,y,feyv,d_xu,d2_xu,d_yv,u,u_ant,v,temp,pres,Resu,gamma_u,Ri,dt)
    WRITE(*,*) 'tiempo ',itera,MAXVAL(DABS(b_o))
    IF( MAXVAL(DABS(b_o))<conv_paso )EXIT
    !*************************************************
  END DO ALGORITMO_SIMPLE  !final del algoritmo SIMPLE
  itera = itera + 1
  IF( mod(itera,100)==0 )WRITE(*,*) 'tiempa ',itera,res_fluido_u,MAXVAL(DABS(Resu)),MAXVAL(DABS(b_o))
!   IF(mod(itera,10)==0)THEN
!     CALL entropia_cvt(x,y,u,xu,v,yv,temp,entropia_calor,entropia_viscosa,entropia,entropia_int,temp_int,a_ent,lambda_ent)
!     CALL nusselt(x,y,d_xu,d_yv,temp,nusselt0,nusselt1,i_o,i_1)
!     temp_med = (temp((i_o+i_1)/2,nj/2+1)+temp((i_o+i_1)/2+1,nj/2+1))/2._DBL
!     OPEN(unit = 5,file = 'nuss_sim_n'//njc//'m'//mic//'_R'//Rac//'.dat',access = 'append')
!     WRITE(5,26) tiempo_inicial+itera*dt,nusselt0,nusselt1,temp_med,temp_int,entropia_int
!     CLOSE(unit = 5)
!   ENDIF
  !*********************************
  tiempo   = tiempo_inicial+itera*dt
  temp_ant = temp
  u_ant    = u
  v_ant    = v
  w_ant    = w
!   !************************************
!   !*** Formato de escritura Tecplot ***
!   OPEN(unit=1,file='tprueba.plt')
!   WRITE(1,*)'TITLE     = "Tempe" '
!   WRITE(1,*)'VARIABLES = "x"'
!   WRITE(1,*)'"y"'
!   WRITE(1,*)'"z"'
!   WRITE(1,*)'"Temperature"'        ! variable
!   WRITE(1,*)'"Temperature_a"'        ! variable
!   !WRITE(1,*)'"XD"'
!   !WRITE(1,*)'"DTEMP"'
!   WRITE(1,*)'ZONE T= "Matrix"'
!   WRITE(1,*)'I=',mi+1,' J=',nj+1,' K=',lk+1,' F=POINT'
!   DO kl = 1, lk+1
!     DO j = 1, nj+1
!       DO i = 1, mi+1
! 	WRITE(1,25) x(i),y(j),z(kl),temp(i,j,kl),temp_ant(i,j,kl)
!       END DO
!     END DO
!   END DO
!   CLOSE(unit=1)
END DO !*************termina el paquete de iteraciones
!*****************************************************
!*****************************************************
itera_total = itera_inicial+itera
millar      = itera_total/(1000*paq_itera)
centena     = (itera_total-millar*1000*paq_itera)/(100*paq_itera)
decena      = (itera_total-millar*1000*paq_itera-centena*100*paq_itera)/(10*paq_itera)
unidad      = (itera_total-millar*1000*paq_itera-centena*100*paq_itera-decena*10*paq_itera)/(paq_itera)
decima      = (itera_total-millar*1000*paq_itera-centena*100*paq_itera-decena*10*paq_itera-unidad*paq_itera)/(paq_itera/10)
WRITE(dec,16)decima;16 format(I1)
WRITE(un,16) unidad
WRITE(de,16) decena
WRITE(ce,16) centena
WRITE(m,16)  millar
DO k = 1, lk+1
  DO j = 1, nj+1
    DO i = 2, mi
      uf(i,j,k) = (u(i,j,k)+u(i-1,j,k))/2._DBL
    END DO
  END DO
END DO
DO k = 1, lk+1
  DO i = 1, mi+1
    DO j = 2, nj
      vf(i,j,k) = (v(i,j,k)+v(i,j-1,k))/2._DBL
    END DO
  END DO
END DO
DO j = 1, nj+1
  DO i = 1, mi+1
    DO k = 2, lk
      wf(i,j,k) = (w(i,j,k)+w(i,j,k-1))/2._DBL
    END DO
  END DO
END DO
DO k = 1, lk+1
  DO j = 1, nj+1
    uf(1,j,k)    = u(1,j,k)
    uf(mi+1,j,k) = u(mi,j,k)
  END DO
END DO
DO k = 1, lk+1
  DO i = 1, mi+1
    vf(i,1,k)    = v(i,1,k)
    vf(i,nj+1,k) = v(i,nj,k)
  END DO
END DO
DO j = 1, nj+1
  DO i = 1, mi+1
    wf(i,j,1)    = w(i,j,1)
    wf(i,j,lk+1) = w(i,j,lk)
  END DO
END DO
! !************************************
WRITE(*,*) 'itera_total=',itera_total
WRITE(*,103) nusselt0,nusselt1
WRITE(*,104) MAXVAL(DABS(b_o)),MAXVAL(DABS(Resu))
WRITE(*,105) MAXVAL(DABS(Restemp)),MAXVAL(DABS(Resv))
WRITE(*,*)' '
103 FORMAT(1X,'N_Izq=',D23.15,', N_Der=',D23.15)
104 FORMAT(1X,'b_o  =',D23.15,', Res_u=',D23.15)
105 FORMAT(1X,'Res_T=',D23.15,', Res_v=',D23.15)
! !********************************
! !*** Formato de escritura dat ***
!archivo de escritura de la malla u
OPEN(unit=1,file='out_n'//njc//'m'//mic//'k'//lkc//'_R'//Rac//'u.dat')
WRITE(1,*) i_o,i_1,itera_total,ao
DO k=1,lk+1
  DO j=1,nj+1
    DO i=1,mi
      WRITE(1,24) xu(i),yp(j),zp(k),u(i,j,k)
    END DO
  END DO
END DO
CLOSE(unit=1)
!**********************************
!archivo de escritura de la malla v
OPEN(unit=2,file='out_n'//njc//'m'//mic//'k'//lkc//'_R'//Rac//'v.dat')
WRITE(2,*) j_o,j_1,itera_total,ao
DO k=1,lk+1
  DO j=1,nj
    DO i=1,mi+1
      WRITE(2,24) xp(i),yv(j),zp(k),v(i,j,k)
    END DO
  END DO
END DO
CLOSE(unit=2)
!**********************************
!archivo de escritura de la malla w
OPEN(unit=3,file='out_n'//njc//'m'//mic//'k'//lkc//'_R'//Rac//'w.dat')
WRITE(3,*) i_o,i_1,itera_total,ao
DO k=1,lk
  DO j=1,nj+1
    DO i=1,mi+1
      WRITE(3,24) xp(i),yp(j),zw(k),w(i,j,k)
    END DO
  END DO
END DO
CLOSE(unit=3)
!**********************************
!archivo de escritura de la malla p
OPEN(unit=4,file='out_n'//njc//'m'//mic//'k'//lkc//'_R'//Rac//'p.dat')
WRITE(4,*) i_o,i_1,itera_total,ao
DO k=1,lk+1
  DO j=1,nj+1
    DO i=1,mi+1
      WRITE(4,25) xp(i),yp(j),zp(k),temp(i,j,k),pres(i,j,k)
    END DO
  END DO
END DO
CLOSE(unit=4)
!********************************
!*** Formato de escritura VTK ***
archivo = 'n'//njc//'m'//mic//'k'//lkc//'R'//Rac//'/t_'//m//ce//de//un//dec//'.vtk'
CALL postprocess_vtk(xp,yp,zp,uf,vf,wf,pres,temp,b_o,trim(archivo))
! !************************************
! !*** Formato de escritura Tecplot ***
! OPEN(unit=1,file='n'//njc//'m'//mic//'R'//Rac//'/t'//m//ce//de//un//dec//'.plt')
! WRITE(1,*)'TITLE     = "Tempe" '
! WRITE(1,*)'VARIABLES = "x"'
! WRITE(1,*)'"y"'
! WRITE(1,*)'"z"'
! WRITE(1,*)'"Temperature"'        ! variable
! WRITE(1,*)'"U"'
! WRITE(1,*)'"V"'
! WRITE(1,*)'"W"'
! WRITE(1,*)'"Pressure"'
! WRITE(1,*)'"Entropy"'
! WRITE(1,*)'"Residual"'
! !WRITE(1,*)'"XD"'
! !WRITE(1,*)'"DTEMP"'
! WRITE(1,*)'ZONE T= "Matrix"'
! WRITE(1,*)'I=',mi+1,' J=',nj+1,' K=',lk+1,' F=POINT'
! DO k = 1, lk+1
!   DO j = 1, nj+1
!     DO i = 1, mi+1
!       WRITE(1,27) x(i),y(j),z(k),temp(i,j,k),uf(i,j,k),vf(i,j,k),wf(i,j,k),pres(i,j,k),entropia(i,j,k),b_o(i,j,k)
!     END DO
!   END DO
! END DO
! CLOSE(unit=1)
!************************************
!*** Formato de escritura Tecplot ***
! OPEN(unit=1,file='n'//njc//'m'//mic//'R'//Rac//'/res'//m//ce//de//un//dec//'.plt')
! WRITE(1,*)'TITLE     = "vel" '
! WRITE(1,*)'VARIABLES = "x"'
! WRITE(1,*)'"y"'
! WRITE(1,*)'"z"'
! WRITE(1,*)'"w"'
! !WRITE(1,*)'"XD"'
! !WRITE(1,*)'"DTEMP"'
! WRITE(1,*)'ZONE T= "Matrix"'
! WRITE(1,*)'I=',mi+1,' J=',nj+1,' K=',lk,' F=POINT'
! DO k =1, lk
!   DO j = 1, nj+1
!     DO i = 1, mi+1
!       WRITE(1,24) x(i),y(j),zw(k),w(i,j,k)
!     END DO
!   END DO
! END DO
! CLOSE(unit=1)
!*******************************
!formatos de escritura y lectura
24 FORMAT (4D23.15)
25 FORMAT (5D23.15)
26 FORMAT (6D23.15)
27 FORMAT (10e15.7)
END DO!***********final del repetidor principal
!**********************************************
END PROGRAM SIMPLE
