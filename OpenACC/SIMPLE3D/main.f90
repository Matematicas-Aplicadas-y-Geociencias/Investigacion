PROGRAM SIMPLE
!
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
use malla, only : indexu, indeyu, indexp
use malla, only : indexv, indeyv, indezv
use malla, only : indeyp
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
use ec_momento, only : ensambla_velv_x
use ec_momento, only : ensambla_velv_y
use ec_momento, only : ensambla_velv_z
use ec_momento, only : ensambla_velw_x
use ec_momento, only : ensambla_velw_y
use ec_momento, only : ensambla_velw_z
use ec_momento, only : residuo_u
!
use ec_energia, only : temp, temp_ant, ftemp
use ec_energia, only : fuente_con_temp, fuente_lin_temp
use ec_energia, only : gamma_ener
use ec_energia, only : rel_ener
use ec_energia, only : ensambla_energia_x
use ec_energia, only : ensambla_energia_y
use ec_energia, only : ensambla_energia_z
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
real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)) :: a1, b1, c1, r1
!
real(kind=DBL) :: residuo, maxbo, error, erro1
!
! --------------------------------------------------------------------------------
!
integer :: ii, jj, kk
integer :: itermax, ecuamax
integer :: simpmax, iter_simp, tt
INTEGER :: i,j,k,kl,l,itera_total,itera,itera_inicial,i_o,i_1,j_o,j_1,paq_itera
INTEGER :: millar,centena,decena,unidad,decima,id,nthreads
!*******************************************
! Variables del flujo,entropia,nusselt e inc'ognitas,residuos,relajaci'on y convergencia
REAL(kind=DBL), DIMENSION(mi,nj+1,lk+1)   :: gamma_u
REAL(kind=DBL), DIMENSION(mi+1,nj,lk+1)   :: gamma_v
REAL(kind=DBL), DIMENSION(mi+1,nj+1,lk)   :: gamma_w
REAL(kind=DBL), DIMENSION(mi+1,nj+1,lk+1) :: dtemp,Restemp,gamma_t
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
CHARACTER(len=64):: entrada_u,entrada_v,entrada_w,entrada_tp,entrada_xyz
character(len=64):: archivo=repeat(' ',64)
LOGICAL          :: res_fluido_u
!****************************************
!declaraci´on de variable DBL
REAL(kind=DBL) :: var2=0.0_DBL
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
  read (10,*) itermax     ! iteraciones m'aximas
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
  read (10,*) simpmax     ! iteraciones m'aximas de SIMPLE
  read (10,*) ecuamax     ! iteraciones m'aximas de las ecuaciones
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
gamma_ener  = 1._DBL/(Pr*Ra)
fuente_con_u = 0._DBL
fuente_lin_u = 0._DBL
fuente_con_v = 0._DBL
fuente_lin_v = 0._DBL
fuente_con_w = 0._DBL
fuente_lin_w = 0._DBL
fuente_con_temp = 0._DBL
fuente_lin_temp = 0._DBL
rel_ener = rel_tem
gamma_v = 1._DBL/Ra      !sqrt(Pr/Ra)
gamma_w = 1._DBL/Ra      !sqrt(Pr/Ra)
Ri      = Rin
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
tiempo = tiempo_inicial
! u_ant = 0.0_DBL
! v_ant = 0.0_DBL
! w_ant =-1.0_DBL
! temp_ant = 0.0_DBL !temp_ant
corr_pres = 0.0_DBL
au    = 1._DBL
av    = 1._DBL
aw    = 1._DBL
! pres  = 0.0_DBL
b_o   = 0._DBL
itera = 1
itera_total = itera_total+itera
entropia = 0.0_DBL
! simpmax = 200
! itermax = 6000
!************************************************
!escribe las caracter´isticas de las variable DBL
WRITE(*,100) 'Doble',KIND(var2),PRECISION(var2),RANGE(var2)
100 FORMAT(1X,A,': kind= ',I2,', Precision= ',I2,' Rango= ',I3)
WRITE(*,*)' '
!escribe informaci'on de los parametros usados
WRITE(*,101) Ra,Pr,rel_pres,rel_vel
WRITE(*,102) itera_inicial,mi,nj,lk
WRITE(*,106) dt,itermax
WRITE(*,*)' '
101 FORMAT(1X,'Ra=',F12.3,', Pr=',F8.3', rel_pres=',F8.3', rel_vel=',F8.3)
102 FORMAT(1X,'Iteracion inicial=',I7,', mi=',I4,', nj=',I4,', lk=',I4)
106 FORMAT(1X,'Paso de tiempo dt=',F8.6,', iter final itermax=',I8)
!*********************************************************
DO l=1,itermax/paq_itera   !inicio del repetidor principal
   DO kl=1,paq_itera          !inicio del paquete iteraciones
      ALGORITMO_SIMPLE: DO  iter_simp = 1, simpmax     !inicio del algoritmo SIMPLE
         DO tt= 1, ecuamax
            !
            !----------------------------------------------------------
            !----------------------------------------------------------
            !
            !      Se resuelve la ecuaci'on de momento
            !
            !----------------------------------------------------------
            !----------------------------------------------------------
            !
            !--------------------------
            !--------------------------
            !
            !     Ecuaci'on para u
            !
            !--------------------------
            !--------------------------
            !
            !$acc parallel loop gang collapse(2) !async(stream1)
            !$OMP PARALLEL DO COLLAPSE(3)
            inicializacion_fu: do kk = 1, lk+1
               do jj = 1, nj+1
                  do ii = 1, mi
                     fu(ii,jj,kk) = u(ii,jj,kk)
                  end do
               end do
            end do inicializacion_fu
            !$OMP END PARALLEL DO
            !
            !----------------------------------------
            !
            ! Se ensamblan las matrices tridiagonales
            ! en la direcci'on de y para u
            !
            !$acc parallel loop gang !async(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
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
            !$OMP END PARALLEL DO
            !
            !-------------------------
            !
            ! Condiciones de frontera direcci\'on y
            !
            !$acc parallel loop vector !async(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED)
            cond_fron_u_direc_y: do kk = 2, lk
               do ii = 2, mi-1
                  !***********************
                  !Condiciones de frontera
                  a1(indeyu(1,ii,kk))    = 0.0_DBL
                  b1(indeyu(1,ii,kk))    = 1.0_DBL !-1.0_DBL !
                  c1(indeyu(1,ii,kk))    = 0.0_DBL ! 1.0_DBL
                  r1(indeyu(1,ii,kk))    = 0.0_DBL
                  au(ii,1,kk)            = 1.0e40_DBL
                  !
                  a1(indeyu(nj+1,ii,kk)) = 0.0_DBL !-1.0_DBL
                  b1(indeyu(nj+1,ii,kk)) = 1.0_DBL
                  c1(indeyu(nj+1,ii,kk)) = 0.0_DBL
                  r1(indeyu(nj+1,ii,kk)) = 0.0_DBL
                  au(ii,nj+1,kk)         = 1.0e40_DBL
                  !
               end do
            end do cond_fron_u_direc_y
            !$OMP END PARALLEL DO
            !
            !----------------------------------------------
            !
            ! Soluci\'on de la ec. de momento para u
            ! en direcci\'on y
            !
            !$acc parallel loop gang async(stream1) wait(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
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
            !$OMP END PARALLEL DO
            !
            !----------------------------------
            !
            ! Actualizaci\'on de la velocidad u
            ! en direcci\'on y
            !
            !$acc parallel loop gang collapse(2) !async(stream2) wait(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(3)
            do kk = 2, lk
               do ii = 2, mi-1
                  do jj = 1, nj+1
                     u(ii,jj,kk) = r1(indeyu(jj,ii,kk))
                  end do
               end do
            end do
            !$OMP END PARALLEL DO
            !
            !----------------------------------------
            !
            ! Se ensamblan las matrices tridiagonales
            ! en la direcci'on de z para u
            !
            !$acc parallel loop gang !async(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
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
            !$OMP END PARALLEL DO
            !
            !-------------------------
            !
            ! Condiciones de frontera direcci\'on z
            !
            !$acc parallel loop vector !async(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED)
            cond_fron_u_direc_z: do ii = 2, mi-1
               do jj = 2, nj
                  !***********************
                  !Condiciones de frontera
                  a1(indezp(1,jj,ii))    = 0.0_DBL
                  b1(indezp(1,jj,ii))    = 1.0_DBL
                  c1(indezp(1,jj,ii))    = 0.0_DBL
                  r1(indezp(1,jj,ii))    = 0.0_DBL*&
                       &dtanh((tiempo+dt)/0.1)
                  au(ii,jj,1)            = 1.0e40_DBL
                  !
                  a1(indezp(lk+1,jj,ii)) = 0.0_DBL
                  b1(indezp(lk+1,jj,ii)) = 1.0_DBL
                  c1(indezp(lk+1,jj,ii)) = 0.0_DBL
                  r1(indezp(lk+1,jj,ii)) = 1.0_DBL*&
                       &dtanh((tiempo+dt)/0.3)
                  au(ii,jj,lk+1)         = 1.0e40_DBL
                  !
               end do
            end do cond_fron_u_direc_z
            !$OMP END PARALLEL DO
            !
            !----------------------------------------------
            !
            ! Soluci\'on de la ec. de momento para u
            ! en direcci\'on z
            !
            !$acc parallel loop gang async(stream1) wait(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
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
            !$OMP END PARALLEL DO
            !
            !----------------------------------
            !
            ! Actualizaci\'on de la velocidad u
            ! en direcci\'on z
            !
            !$acc parallel loop gang collapse(2) !async(stream2) wait(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(3)
            do ii = 2, mi-1
               do jj = 2, nj
                  do kk = 1, lk+1
                     u(ii,jj,kk) = r1(indezp(kk,jj,ii))
                  end do
               end do
            end do
            !$OMP END PARALLEL DO
            !
            !----------------------------------------
            !
            ! Se ensamblan las matrices tridiagonales
            ! en la direcci'on de x para u
            !
            !$acc parallel loop gang !async(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
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
            !$OMP END PARALLEL DO
            !
            !-------------------------
            !
            ! Condiciones de frontera direcci\'on x
            !
            !$acc parallel loop vector !async(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED)
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
            !OMP END PARALLEL DO
            !
            !----------------------------------------------
            !
            ! Soluci\'on de la ec. de momento para u
            ! en direcci\'on x
            !
            !$acc parallel loop gang async(stream1) wait(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
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
            !$OMP END PARALLEL DO
            !
            !----------------------------------
            !
            ! Actualizaci\'on de la velocidad u
            ! en direcci\'on x
            !
            !$acc parallel loop gang collapse(2) !async(stream2) wait(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(3)
            do kk = 2, lk
               do jj = 2, nj
                  do ii = 1, mi
                     u(ii,jj,kk) = r1(indexu(ii,jj,kk))
                  end do
               end do
            end do
            !$OMP END PARALLEL DO
            !
            error = 0.0_DBL
            !$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:error)
            calcula_fu: do kk = 2, lk
               do jj = 2, nj
                  do ii = 2, mi-1
                     error = error + (fu(ii,jj,kk)-u(ii,jj,kk))*&
                          &(fu(ii,jj,kk)-u(ii,jj,kk))
                  end do
               end do
            end do calcula_fu
            !$OMP END PARALLEL DO
            error = dsqrt(error)
            !
            !--------------------------
            !--------------------------
            !
            !   Ecuaci'on para v
            !
            !--------------------------
            !--------------------------
            !
            !$acc parallel loop gang collapse(2) !async(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(3)
            inicializacion_fv: do kk = 1, lk+1
               do jj = 1, nj
                  do ii = 1, mi+1
                     fv(ii,jj,kk) = v(ii,jj,kk)
                  end do
               end do
            end do inicializacion_fv
            !OMP END PARALLEL DO
            !
            !----------------------------------------
            !
            ! Se ensamblan las matrices tridiagonales
            ! en la direcci'on de x para v
            !
            !$acc parallel loop gang !async(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
            ensa_velv_dir_z: do ii = 2, mi
               do jj = 2, nj-1
                  do kk = 2, lk
                     call ensambla_velv_z(&
                          &deltaxv,&
                          &deltayv,&
                          &deltazv,&
                          &deltaxu,&
                          &deltayp,&
                          &deltazw,&
                          &fexp,&
                          &feyp,&
                          &fezp,&
                          &feyv,&
                          &gamma_momen,&
                          &u,&
                          &v,&
                          &v_ant,&
                          &w,&
                          &temp,&
                          &pres,&
                          &fuente_con_v,&
                          &fuente_lin_v,&
                          &Ri,&
                          &dt,&
                          &rel_vel,&
                          &a1,b1,c1,r1,&
                          &kk,jj,ii&
                          &)
                  end do
               end do
            end do ensa_velv_dir_z
            !$OMP END PARALLEL DO
            !
            !-----------------------------------------
            !
            ! Condiciones de frontera v direcci\'on y
            !
            !$acc parallel loop vector !async(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED)
            cond_fron_v_direc_z: do ii = 2, mi
               do jj = 2, nj-1
                  !***********************
                  !Condiciones de frontera
                  a1(indezv(1,jj,ii))    = 0.0_DBL
                  b1(indezv(1,jj,ii))    = 1.0_DBL
                  c1(indezv(1,jj,ii))    = 0.0_DBL
                  r1(indezv(1,jj,ii))    = 0.0_DBL
                  av(ii,jj,1)            = 1.0e40_DBL
                  !
                  a1(indezv(lk+1,jj,ii)) = 0.0_DBL
                  b1(indezv(lk+1,jj,ii)) = 1.0_DBL
                  c1(indezv(lk+1,jj,ii)) = 0.0_DBL
                  r1(indezv(lk+1,jj,ii)) = 0.0_DBL
                  av(ii,jj,lk+1)         = 1.0e40_DBL
                  !
               end do
            end do cond_fron_v_direc_z
            !$OMP END PARALLEL DO
            !
            !----------------------------------------------
            !
            ! Soluci\'on de la ec. de momento para v
            ! en direcci\'on y
            !
            !$acc parallel loop gang async(stream1) wait(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
            sol_v_dir_z: do ii = 2, mi
               do jj = 2, nj-1
                  !
                  call tridiagonal(&
                       &a1(indezv(1,jj,ii):indezv(lk+1,jj,ii)),&
                       &b1(indezv(1,jj,ii):indezv(lk+1,jj,ii)),&
                       &c1(indezv(1,jj,ii):indezv(lk+1,jj,ii)),&
                       &r1(indezv(1,jj,ii):indezv(lk+1,jj,ii)),&
                       &lk+1)
                  !
               end do
            end do sol_v_dir_z
            !$OMP END PARALLEL DO
            !
            !----------------------------------
            !
            ! Actualizaci\'on de la velocidad v
            ! en direcci\'on y
            !
            !$acc parallel loop gang collapse(2) !async(stream2) wait(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(3)
            do ii = 2, mi
               do jj = 2, nj-1
                  do kk = 1, lk+1
                     v(ii,jj,kk) = r1(indezv(kk,jj,ii))
                  end do
               end do
            end do
            !$OMP END PARALLEL DO
            !
            !----------------------------------------
            !
            ! Se ensamblan las matrices tridiagonales
            ! en la direcci'on de x para v
            !
            !$acc parallel loop gang !async(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
            ensa_velv_dir_y: do kk = 2, lk
               do ii = 2, mi
                  do jj = 2, nj-1
                     call ensambla_velv_y(&
                          &deltaxv,&
                          &deltayv,&
                          &deltazv,&
                          &deltaxu,&
                          &deltayp,&
                          &deltazw,&
                          &fexp,&
                          &feyp,&
                          &fezp,&
                          &feyv,&
                          &gamma_momen,&
                          &u,&
                          &v,&
                          &v_ant,&
                          &w,&
                          &temp,&
                          &pres,&
                          &fuente_con_v,&
                          &fuente_lin_v,&
                          &Ri,&
                          &dt,&
                          &rel_vel,&
                          &a1,b1,c1,r1,&
                          &jj,ii,kk&
                          &)
                  end do
               end do
            end do ensa_velv_dir_y
            !$OMP END PARALLEL DO
            !
            !-----------------------------------------
            !
            ! Condiciones de frontera v direcci\'on y
            !
            !$acc parallel loop vector !async(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED)
            cond_fron_v_direc_y: do kk = 2, lk
               do ii = 2, mi
                  !***********************
                  !Condiciones de frontera
                  a1(indeyv(1,ii,kk))    = 0.0_DBL
                  b1(indeyv(1,ii,kk))    = 1.0_DBL
                  c1(indeyv(1,ii,kk))    = 0.0_DBL
                  r1(indeyv(1,ii,kk))    = 0.0_DBL
                  av(ii,1,kk)            = 1.0e40_DBL
                  !
                  a1(indeyv(nj,ii,kk))   = 0.0_DBL
                  b1(indeyv(nj,ii,kk))   = 1.0_DBL
                  c1(indeyv(nj,ii,kk))   = 0.0_DBL
                  r1(indeyv(nj,ii,kk))   = 0.0_DBL
                  av(ii,nj,kk)           = 1.0e40_DBL
                  !
               end do
            end do cond_fron_v_direc_y
            !$OMP END PARALLEL DO
            !
            !----------------------------------------------
            !
            ! Soluci\'on de la ec. de momento para v
            ! en direcci\'on y
            !
            !$acc parallel loop gang async(stream1) wait(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
            sol_v_dir_y: do kk = 2, lk
               do ii = 2, mi
                  !
                  call tridiagonal(&
                       &a1(indeyv(1,ii,kk):indeyv(nj,ii,kk)),&
                       &b1(indeyv(1,ii,kk):indeyv(nj,ii,kk)),&
                       &c1(indeyv(1,ii,kk):indeyv(nj,ii,kk)),&
                       &r1(indeyv(1,ii,kk):indeyv(nj,ii,kk)),&
                       &nj)
                  !
               end do
            end do sol_v_dir_y
            !$OMP END PARALLEL DO
            !
            !----------------------------------
            !
            ! Actualizaci\'on de la velocidad v
            ! en direcci\'on y
            !
            !$acc parallel loop gang collapse(2) !async(stream2) wait(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(3)
            do kk = 2, lk
               do ii = 2, mi
                  do jj = 1, nj
                     v(ii,jj,kk) = r1(indeyv(jj,ii,kk))
                  end do
               end do
            end do
            !$OMP END PARALLEL DO
            !
            !----------------------------------------
            !
            ! Se ensamblan las matrices tridiagonales
            ! en la direcci'on de x para v
            !
            !$acc parallel loop gang !async(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
            ensa_velv_dir_x: do kk = 2, lk
               do jj = 2, nj-1
                  do ii = 2, mi
                     call ensambla_velv_x(&
                          &deltaxv,&
                          &deltayv,&
                          &deltazv,&
                          &deltaxu,&
                          &deltayp,&
                          &deltazw,&
                          &fexp,&
                          &feyp,&
                          &fezp,&
                          &feyv,&
                          &gamma_momen,&
                          &u,&
                          &v,&
                          &v_ant,&
                          &w,&
                          &temp,&
                          &pres,&
                          &fuente_con_v,&
                          &fuente_lin_v,&
                          &Ri,&
                          &dt,&
                          &rel_vel,&
                          &a1,b1,c1,r1,&
                          &av,&
                          &ii,jj,kk&
                          &)
                  end do
               end do
            end do ensa_velv_dir_x
            !$OMP END PARALLEL DO
            !
            !-----------------------------------------
            !
            ! Condiciones de frontera v direcci\'on x
            !
            !$acc parallel loop vector !async(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED)
            cond_fron_v_direc_x: do kk = 2, lk
               do jj = 2, nj-1
                  !***********************
                  !Condiciones de frontera
                  a1(indexv(1,jj,kk))      = 0.0_DBL
                  b1(indexv(1,jj,kk))      = 1.0_DBL
                  c1(indexv(1,jj,kk))      = 0.0_DBL
                  r1(indexv(1,jj,kk))      = 0.0_DBL
                  av(1,jj,kk)              = 1.0e40_DBL
                  !
                  a1(indexv(mi+1,jj,kk))   = 0.0_DBL
                  b1(indexv(mi+1,jj,kk))   = 1.0_DBL
                  c1(indexv(mi+1,jj,kk))   = 0.0_DBL
                  r1(indexv(mi+1,jj,kk))   = 0.0_DBL
                  av(mi+1,jj,kk)           = 1.0e40_DBL
                  !
               end do
            end do cond_fron_v_direc_x
            !$OMP END PARALLEL DO
            !
            !----------------------------------------------
            !
            ! Soluci\'on de la ec. de momento para v
            ! en direcci\'on x
            !
            !$acc parallel loop gang async(stream1) wait(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
            sol_v_dir_x: do kk = 2, lk
               do jj = 2, nj-1
                  !
                  call tridiagonal(&
                       &a1(indexv(1,jj,kk):indexv(mi+1,jj,kk)),&
                       &b1(indexv(1,jj,kk):indexv(mi+1,jj,kk)),&
                       &c1(indexv(1,jj,kk):indexv(mi+1,jj,kk)),&
                       &r1(indexv(1,jj,kk):indexv(mi+1,jj,kk)),&
                       &mi+1)
                  !
               end do
            end do sol_v_dir_x
            !$OMP END PARALLEL DO
            !
            !----------------------------------
            !
            ! Actualizaci\'on de la velocidad v
            ! en direcci\'on x
            !
            !$acc parallel loop gang collapse(2) !async(stream2) wait(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(3)
            do kk = 2, lk
               do jj = 2, nj-1
                  do ii = 1, mi+1
                     v(ii,jj,kk) = r1(indexv(ii,jj,kk))
                  end do
               end do
            end do
            !$OMP END PARALLEL DO
            !
            ! Se suma el error al error de la ecuaci'on de u
            !
            erro1 = 0.0_DBL
            !$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:error)
            calcula_fv: do kk = 2, lk
               do jj = 2, nj-1
                  do ii = 2, mi
                     erro1 = erro1 + (fv(ii,jj,kk)-v(ii,jj,kk))*&
                          &(fv(ii,jj,kk)-v(ii,jj,kk))
                  end do
               end do
            end do calcula_fv
            !$OMP END PARALLEL DO
            erro1=dsqrt(erro1)+error
            !
            !--------------------------
            !--------------------------
            !
            !     Ecuaci'on para w
            !
            !--------------------------
            !--------------------------
            !
            !$acc parallel loop gang collapse(2) !async(stream1)
            !$OMP PARALLEL DO COLLAPSE(3)
            inicializacion_fw: do kk = 1, lk
               do jj = 1, nj+1
                  do ii = 1, mi+1
                     fw(ii,jj,kk) = w(ii,jj,kk)
                  end do
               end do
            end do inicializacion_fw
            !$OMP END PARALLEL DO
            !
            !----------------------------------------
            !
            ! Se ensamblan las matrices tridiagonales
            ! en la direcci'on de y para w
            !
            !$acc parallel loop gang !async(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
            ensa_velw_dir_z: do ii = 2, mi
               do jj = 2, nj
                  do kk = 2, lk-1
                     call ensambla_velw_z(&
                          &deltaxw,&
                          &deltayw,&
                          &deltazw,&
                          &deltaxu,&
                          &deltayv,&
                          &deltazp,&
                          &fexp,&
                          &feyp,&
                          &fezp,&
                          &fezw,&
                          &gamma_momen,&
                          &u,&
                          &v,&
                          &w,&
                          &w_ant,&
                          &temp,&
                          &pres,&
                          &fuente_con_w,&
                          &fuente_lin_w,&
                          &Ri,&
                          &dt,&
                          &rel_vel,&
                          &a1,b1,c1,r1,&
                          &kk,jj,ii&
                          &)
                  end do
               end do
            end do ensa_velw_dir_z
            !$OMP END PARALLEL DO
            !
            !-------------------------
            !
            ! Condiciones de frontera direcci\'on x
            !
            !$acc parallel loop vector !async(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED)
            cond_fron_w_direc_z: do ii = 2, mi
               do jj = 2, nj
                  !***********************
                  !Condiciones de frontera
                  a1(indezw(1,jj,ii))  = 0.0_DBL
                  b1(indezw(1,jj,ii))  = 1.0_DBL
                  c1(indezw(1,jj,ii))  = 0.0_DBL
                  r1(indezw(1,jj,ii))  = 0.0_DBL
                  aw(ii,jj,1)          = 1.0e40_DBL
                  !
                  a1(indezw(lk,jj,ii)) = 0.0_DBL
                  b1(indezw(lk,jj,ii)) = 1.0_DBL
                  c1(indezw(lk,jj,ii)) = 0.0_DBL
                  r1(indezw(lk,jj,ii)) = 0.0_DBL
                  aw(ii,jj,lk)         = 1.0e40_DBL
                  !
               end do
            end do cond_fron_w_direc_z
            !$OMP END PARALLEL DO
            !
            !----------------------------------------------
            !
            ! Soluci\'on de la ec. de momento para u
            ! en direcci\'on y
            !
            !$acc parallel loop gang async(stream1) wait(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
            sol_w_dir_z: do ii = 2, mi
               do jj = 2, nj
                  !
                  call tridiagonal(&
                       &a1(indezw(1,jj,ii):indezw(lk,jj,ii)),&
                       &b1(indezw(1,jj,ii):indezw(lk,jj,ii)),&
                       &c1(indezw(1,jj,ii):indezw(lk,jj,ii)),&
                       &r1(indezw(1,jj,ii):indezw(lk,jj,ii)),&
                       &lk)
                  !
               end do
            end do sol_w_dir_z
            !$OMP END PARALLEL DO
            !
            !----------------------------------
            !
            ! Actualizaci\'on de la velocidad w
            ! en direcci\'on y
            !
            !$acc parallel loop gang collapse(2) !async(stream2) wait(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(3)
            do ii = 2, mi
               do jj = 2, nj
                  do kk = 1, lk
                     w(ii,jj,kk) = r1(indezw(kk,jj,ii))
                  end do
               end do
            end do
            !$OMP END PARALLEL DO
            !
            !----------------------------------------
            !
            ! Se ensamblan las matrices tridiagonales
            ! en la direcci'on de y para w
            !
            !$acc parallel loop gang !async(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
            ensa_velw_dir_y: do kk = 2, lk-1
               do ii = 2, mi
                  do jj = 2, nj
                     call ensambla_velw_y(&
                          &deltaxw,&
                          &deltayw,&
                          &deltazw,&
                          &deltaxu,&
                          &deltayv,&
                          &deltazp,&
                          &fexp,&
                          &feyp,&
                          &fezp,&
                          &fezw,&
                          &gamma_momen,&
                          &u,&
                          &v,&
                          &w,&
                          &w_ant,&
                          &temp,&
                          &pres,&
                          &fuente_con_w,&
                          &fuente_lin_w,&
                          &Ri,&
                          &dt,&
                          &rel_vel,&
                          &a1,b1,c1,r1,&
                          &jj,ii,kk&
                          &)
                  end do
               end do
            end do ensa_velw_dir_y
            !$OMP END PARALLEL DO
            !
            !-------------------------
            !
            ! Condiciones de frontera direcci\'on x
            !
            !$acc parallel loop vector !async(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED)
            cond_fron_w_direc_y: do kk = 2, lk-1
               do ii = 2, mi
                  !***********************
                  !Condiciones de frontera
                  a1(indeyp(1,ii,kk))    = 0.0_DBL
                  b1(indeyp(1,ii,kk))    = 1.0_DBL !-1.0_DBL
                  c1(indeyp(1,ii,kk))    = 0.0_DBL ! 1.0_DBL
                  r1(indeyp(1,ii,kk))    = 0.0_DBL
                  aw(ii,1,kk)            = 1.0e40_DBL
                  !
                  a1(indeyp(nj+1,ii,kk)) = 0.0_DBL !-1.0_DBL
                  b1(indeyp(nj+1,ii,kk)) = 1.0_DBL
                  c1(indeyp(nj+1,ii,kk)) = 0.0_DBL
                  r1(indeyp(nj+1,ii,kk)) = 0.0_DBL
                  aw(ii,nj+1,kk)         = 1.0e40_DBL
                  !
               end do
            end do cond_fron_w_direc_y
            !$OMP END PARALLEL DO
            !
            !----------------------------------------------
            !
            ! Soluci\'on de la ec. de momento para u
            ! en direcci\'on y
            !
            !$acc parallel loop gang async(stream1) wait(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
            sol_w_dir_y: do kk = 2, lk-1
               do ii = 2, mi
                  !
                  call tridiagonal(&
                       &a1(indeyp(1,ii,kk):indeyp(nj+1,ii,kk)),&
                       &b1(indeyp(1,ii,kk):indeyp(nj+1,ii,kk)),&
                       &c1(indeyp(1,ii,kk):indeyp(nj+1,ii,kk)),&
                       &r1(indeyp(1,ii,kk):indeyp(nj+1,ii,kk)),&
                       &nj+1)
                  !
               end do
            end do sol_w_dir_y
            !$OMP END PARALLEL DO
            !
            !----------------------------------
            !
            ! Actualizaci\'on de la velocidad w
            ! en direcci\'on y
            !
            !$acc parallel loop gang collapse(2) !async(stream2) wait(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(3)
            do kk = 2, lk-1
               do ii = 2, mi
                  do jj = 1, nj+1
                     w(ii,jj,kk) = r1(indeyp(jj,ii,kk))
                  end do
               end do
            end do
            !$OMP END PARALLEL DO
            !
            !
            !----------------------------------------
            !
            ! Se ensamblan las matrices tridiagonales
            ! en la direcci'on de x para w
            !
            !$acc parallel loop gang !async(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
            ensa_velw_dir_x: do kk = 2, lk-1
               do jj = 2, nj
                  do ii = 2, mi
                     call ensambla_velw_x(&
                          &deltaxw,&
                          &deltayw,&
                          &deltazw,&
                          &deltaxu,&
                          &deltayv,&
                          &deltazp,&
                          &fexp,&
                          &feyp,&
                          &fezp,&
                          &fezw,&
                          &gamma_momen,&
                          &u,&
                          &v,&
                          &w,&
                          &w_ant,&
                          &temp,&
                          &pres,&
                          &fuente_con_w,&
                          &fuente_lin_w,&
                          &Ri,&
                          &dt,&
                          &rel_vel,&
                          &a1,b1,c1,r1,&
                          &aw,&
                          &ii,jj,kk&
                          &)
                  end do
               end do
            end do ensa_velw_dir_x
            !$OMP END PARALLEL DO
            !
            !-------------------------
            !
            ! Condiciones de frontera direcci\'on x
            !
            !$acc parallel loop vector !async(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED)
            cond_fron_w_direc_x: do kk = 2, lk-1
               do jj = 2, nj
                  !***********************
                  !Condiciones de frontera
                  a1(indexp(1,jj,kk))    = 0.0_DBL
                  b1(indexp(1,jj,kk))    = 1.0_DBL
                  c1(indexp(1,jj,kk))    = 0.0_DBL
                  r1(indexp(1,jj,kk))    = 0.0_DBL
                  aw(1,jj,kk)            = 1.0e40_DBL
                  !
                  a1(indexp(mi+1,jj,kk)) = 0.0_DBL
                  b1(indexp(mi+1,jj,kk)) = 1.0_DBL
                  c1(indexp(mi+1,jj,kk)) = 0.0_DBL
                  r1(indexp(mi+1,jj,kk)) = 0.0_DBL
                  aw(mi+1,jj,kk)         = 1.0e40_DBL
                  !
               end do
            end do cond_fron_w_direc_x
            !$OMP END PARALLEL DO
            !
            !----------------------------------------------
            !
            ! Soluci\'on de la ec. de momento para u
            ! en direcci\'on y
            !
            !$acc parallel loop gang async(stream1) wait(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
            sol_w_dir_x: do kk = 2, lk-1
               do jj = 2, nj
                  !
                  call tridiagonal(&
                       &a1(indexp(1,jj,kk):indexp(mi+1,jj,kk)),&
                       &b1(indexp(1,jj,kk):indexp(mi+1,jj,kk)),&
                       &c1(indexp(1,jj,kk):indexp(mi+1,jj,kk)),&
                       &r1(indexp(1,jj,kk):indexp(mi+1,jj,kk)),&
                       &mi+1)
                  !
               end do
            end do sol_w_dir_x
            !$OMP END PARALLEL DO
            !
            !----------------------------------
            !
            ! Actualizaci\'on de la velocidad w
            ! en direcci\'on x
            !
            !$acc parallel loop gang collapse(2) !async(stream2) wait(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(3)
            do kk = 2, lk-1
               do jj = 2, nj
                  do ii = 1, mi+1
                     w(ii,jj,kk) = r1(indexp(ii,jj,kk))
                  end do
               end do
            end do
            !$OMP END PARALLEL DO
            !
            ! Se suma el error al error de la ecuaci'on de u y de v
            !
            error = 0.0_DBL
            !$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION (+:error)
            calcula_fw: do kk = 2, lk-1
               do jj = 2, nj
                  do ii = 2, mi
                     error = error + (fw(ii,jj,kk)-w(ii,jj,kk))*&
                          &(fw(ii,jj,kk)-w(ii,jj,kk))
                  end do
               end do
            end do calcula_fw
            !$OMP END PARALLEL DO
            error = dsqrt(error) + erro1
            !****************************************
            ! Criterio de convergencia de la velocidad
            ! WRITE(*,*) 'velocidad ',itera, error
            IF( error<conv_u )EXIT
         END DO
         !
         !-------------------------------------------------------------------------------
         !-------------------------------------------------------------------------------
         !
         !                   Se calcula la correcci'on de la presi'on
         !
         !-------------------------------------------------------------------------------
         !-------------------------------------------------------------------------------
         !
         !$acc parallel loop gang collapse(2) !async(stream2)
         !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(3)    
         inicializa_corrector_presion: do kk=1, lk+1
            do jj = 1, nj+1
               do ii = 1, mi+1
                  corr_pres(ii,jj,kk) = 0.0_DBL
                  fcorr_pres(ii,jj,kk)= 0.0_DBL
               end do
            end do
         end do inicializa_corrector_presion
         !$OMP END PARALLEL DO
         !
         correccion_presion: do tt = 1, ecuamax
            !
            !$acc parallel loop gang collapse(2) ! async(stream1) wait(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(3)
            inicializa_fcorr_press: do kk=2, lk
               do jj=2, nj
                  do ii = 2, mi
                     fcorr_pres(ii,jj,kk) = corr_pres(ii,jj,kk)
                  end do
               end do
            end do inicializa_fcorr_press
            !$OMP END PARALLEL DO
            !
            !------------------------------------------
            !
            ! Se ensambla la ecuaci\'on de correcci\'on
            ! de la presi\'on en direcci\'on x
            !
            !$acc parallel loop gang !async(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
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
            !OMP END PARALLEL DO
            !
            !---------------------------------------
            !
            ! Condiciones de frontera direcci\'on x
            !
            !$acc parallel loop vector !async(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED)
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
            !$OMP END PARALLEL DO
            !
            !----------------------------------------------
            !
            ! Soluci\'on de la correcci\'on de la presi\'on
            !
            !$acc parallel loop gang async(stream1) wait(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
            sol_corr_dir_x: do kk = 2, lk
               do jj = 2,nj
                  !
                  call tridiagonal(&
                       &a1(indexp(1,jj,kk):indexp(mi+1,jj,kk)),&
                       &b1(indexp(1,jj,kk):indexp(mi+1,jj,kk)),&
                       &c1(indexp(1,jj,kk):indexp(mi+1,jj,kk)),&
                       &r1(indexp(1,jj,kk):indexp(mi+1,jj,kk)),&
                       &mi+1)
                  !
               end do
            end do sol_corr_dir_x
            !$OMP END PARALLEL DO
            !
            ! Actualizaci\'on del corrector de la presi\'on
            ! en direcci\'on x
            !
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(3)
            do kk = 2, lk
               do jj = 2, nj
                  do ii =1, mi+1
                     corr_pres(ii,jj,kk) = r1(indexp(ii,jj,kk))
                  end do
               end do
            end do
            !$OMP END PARALLEL DO
            !
            !------------------------------------------
            !
            ! Se ensambla la ecuaci\'on de correcci\'on
            ! de la presi\'on en direcci\'on y
            !
            !$acc parallel loop gang !async(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
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
            !$OMP END PARALLEL DO
            !
            !--------------------------------------------
            !--------------------------------------------
            !
            ! Condiciones de frontera direcci\'on y
            !
            !$acc parallel loop vector !async(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED)
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
            !$OMP END PARALLEL DO
            !
            !----------------------------------------------
            !
            ! Soluci\'on de la correcci\'on de la presi\'on
            ! en direcci\'on y
            !
            !$acc parallel loop gang async(stream1) wait(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
            sol_corr_dir_y: do kk = 2, lk
               do ii = 2, mi
                  !
                  call tridiagonal(&
                       &a1(indeyp(1,ii,kk):indeyp(nj+1,ii,kk)),&
                       &b1(indeyp(1,ii,kk):indeyp(nj+1,ii,kk)),&
                       &c1(indeyp(1,ii,kk):indeyp(nj+1,ii,kk)),&
                       &r1(indeyp(1,ii,kk):indeyp(nj+1,ii,kk)),&
                       &nj+1)
                  !
               end do
            end do sol_corr_dir_y
            !$OMP END PARALLEL DO
            !
            ! Actualizaci\'on del corrector de la presi\'on
            ! en direcci\'on y
            !
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(3)
            do kk = 2, lk
               do ii = 2, mi
                  do jj =2, nj
                     corr_pres(ii,jj,kk) = r1(indeyp(jj,ii,kk))
                  end do
               end do
            end do
            !$OMP END PARALLEL DO
            !
            !------------------------------------------
            !
            ! Se ensambla la ecuaci\'on de correcci\'on
            ! de la presi\'on en direcci\'on z
            !
            !$acc parallel loop gang !async(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
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
            !$OMP END PARALLEL DO
            !
            !------------------------------------------------------------------------
            !------------------------------------------------------------------------
            !
            ! Condiciones de frontera direcci\'on z
            !
            !$acc parallel loop vector !async(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED)
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
            !$OMP END PARALLEL DO
            !
            !----------------------------------------------
            !
            ! Soluci\'on de la correcci\'on de la presi\'on
            ! en direcci\'on y
            !
            !$acc parallel loop gang async(stream1) wait(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
            sol_corr_dir_z: do ii = 2, mi
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
            end do sol_corr_dir_z
            !$OMP END PARALLEL DO
            !
            ! Actualizaci\'on del corrector de la presi\'on
            ! en direcci\'on y
            !
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(3)
            do ii = 2, mi
               do jj = 2, nj
                  do kk =1, lk+1
                     corr_pres(ii,jj,kk) = r1(indezp(kk,jj,ii))
                  end do
               end do
            end do
            !$OMP END PARALLEL DO
            !
            error=0._DBL
            !$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:error)
            calcula_fcorr_press: do kk=2, lk
               do jj=2, nj
                  do ii = 2, mi
                     error = (fcorr_pres(ii,jj,kk)-corr_pres(ii,jj,kk))*&
                          &(fcorr_pres(ii,jj,kk)-corr_pres(ii,jj,kk))
                  end do
               end do
            end do calcula_fcorr_press
            !$OMP END PARALLEL DO
            error =dsqrt(error)
            !
            !-----------------------------------------------------------------------------
            !-----------------------------------------------------------------------------
            !
            !****************************************************
            !critero de convergencia del corrector de la presi'on
            ! WRITE(*,*) 'corrector presion ', error
            IF(error<conv_p)EXIT
         END DO correccion_presion
         !*********************
         !se corrige la presion
         !$OMP PARALLEL DO COLLAPSE(3)
         DO k = 2, lk
            DO j = 2, nj
               DO i = 2, mi
                  pres(i,j,k) = pres(i,j,k)+0.6_DBL*corr_pres(i,j,k)
               END DO
            END DO
         END DO
         !$OMP END PARALLEL DO
         !*****************************
         !se actualizan las velocidades
         !$OMP PARALLEL DO COLLAPSE(3)
         DO kk = 2, lk
            DO jj = 2, nj
               DO ii = 2, mi-1
                  u(ii,jj,kk) = u(ii,jj,kk)+deltayu(jj)*deltazu(kk)*&
                       &(corr_pres(ii,jj,kk)-corr_pres(ii+1,jj,kk))/au(ii,jj,kk)
               END DO
            END DO
         END DO
         !$OMP END PARALLEL DO
         !$OMP PARALLEL DO COLLAPSE(3)
         DO kk = 2, lk
            DO jj = 2, nj-1
               DO ii = 2, mi
                  v(ii,jj,kk) = v(ii,jj,kk)+deltaxv(ii)*deltazv(kk)*&
                       &(corr_pres(ii,jj,kk)-corr_pres(ii,jj+1,kk))/av(ii,jj,kk)
               END DO
            END DO
         END DO
         !$OMP END PARALLEL DO
         !$OMP PARALLEL DO COLLAPSE(3)
         DO kk = 2, lk-1
            DO jj = 2, nj
               DO ii = 2, mi
                  w(ii,jj,kk) = w(ii,jj,kk)+deltaxw(ii)*deltayw(jj)*&
                       &(corr_pres(ii,jj,kk)-corr_pres(ii,jj,kk+1))/aw(ii,jj,kk)
               END DO
            END DO
         END DO
         !$OMP END PARALLEL DO
         !*************************
         !
         !------------------------------------------------------------------------------
         !------------------------------------------------------------------------------
         !
         !                   Se resuelve la ec. de la energ\'ia
         !
         !------------------------------------------------------------------------------
         !------------------------------------------------------------------------------
         !
         DO tt = 1, ecuamax
            !$acc parallel loop gang collapse(2) !async(stream1)
            !$OMP PARALLEL DO COLLAPSE(3)
            inicializacion_ftemp: do kk = 1, lk+1
               do jj = 1, nj+1
                  do ii = 1, mi+1
                     ftemp(ii,jj,kk) = temp(ii,jj,kk)
                  end do
               end do
            end do inicializacion_ftemp
            !$OMP END PARALLEL DO
            !
!             CALL temperatura(xp,yp,zp,fexu,feyv,fezw,d_xu,d_yv,d_zw,&
!                  &u,v,w,temp,temp_ant,gamma_t,dt,dtemp,i_o,i_1,rel_tem,j_o,j_1)
            !
            !----------------------------------------
            !
            ! Se ensamblan las matrices tridiagonales
            ! en la direcci'on de x para la energía
            !
            !$acc parallel loop gang !async(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
            ensa_ener_dir_x: do kk = 2, lk
               do jj = 2, nj
                  do ii = 2, mi
                     call ensambla_energia_x(&
                          &deltaxp,&
                          &deltayp,&
                          &deltazp,&
                          &deltaxu,&
                          &deltayv,&
                          &deltazw,&
                          &fexp,&
                          &feyp,&
                          &fezp,&
                          &gamma_ener,&
                          &u,&
                          &v,&
                          &w,&
                          &temp,&
                          &temp_ant,&
                          &fuente_con_temp,&
                          &fuente_lin_temp,&
                          &dt,&
                          &rel_ener,&
                          &a1,b1,c1,r1,&
                          &ii,jj,kk&
                          &)
                  end do
               end do
            end do ensa_ener_dir_x
            !$OMP END PARALLEL DO
            !
            !-------------------------
            !
            ! Condiciones de frontera direcci\'on x
            !
            !$acc parallel loop vector !async(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED)
            cond_fron_ener_direc_x: do kk = 2, lk
               do jj = 2, nj
                  !***********************
                  !Condiciones de frontera
                  a1(indexp(1,jj,kk))    = 0.0_DBL
                  b1(indexp(1,jj,kk))    = 1.0_DBL
                  c1(indexp(1,jj,kk))    = 0.0_DBL
                  r1(indexp(1,jj,kk))    = 1.0_DBL
                  !
                  a1(indexp(mi+1,jj,kk)) = 0.0_DBL
                  b1(indexp(mi+1,jj,kk)) = 1.0_DBL
                  c1(indexp(mi+1,jj,kk)) = 0.0_DBL
                  r1(indexp(mi+1,jj,kk)) = 0.0_DBL
                  !
               end do
            end do cond_fron_ener_direc_x
            !$OMP END PARALLEL DO
            !
            !----------------------------------------------
            !
            ! Soluci\'on de la ec. de la energ\'ia
            ! en direcci\'on x
            !
            !$acc parallel loop gang async(stream1) wait(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
            sol_ener_dir_x: do kk = 2, lk
               do jj = 2, nj
                  !
                  call tridiagonal(&
                       &a1(indexp(1,jj,kk):indexp(mi+1,jj,kk)),&
                       &b1(indexp(1,jj,kk):indexp(mi+1,jj,kk)),&
                       &c1(indexp(1,jj,kk):indexp(mi+1,jj,kk)),&
                       &r1(indexp(1,jj,kk):indexp(mi+1,jj,kk)),&
                       &mi+1)
                  !
               end do
            end do sol_ener_dir_x
            !$OMP END PARALLEL DO
            !
            !----------------------------------
            !
            ! Actualizaci\'on de la temperatura
            ! en direcci\'on x
            !
            !$acc parallel loop gang collapse(2) !async(stream2) wait(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(3)
            do kk = 2, lk
               do jj = 2, nj
                  do ii = 1, mi+1
                     temp(ii,jj,kk) = r1(indexp(ii,jj,kk))
                  end do
               end do
            end do
            !$OMP END PARALLEL DO
            !
            !----------------------------------------
            !
            ! Se ensamblan las matrices tridiagonales
            ! en la direcci'on de y para la energía
            !
            !$acc parallel loop gang !async(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
            ensa_ener_dir_y: do kk = 2, lk
               do ii = 2, mi
                  do jj = 2, nj
                     call ensambla_energia_y(&
                          &deltaxp,&
                          &deltayp,&
                          &deltazp,&
                          &deltaxu,&
                          &deltayv,&
                          &deltazw,&
                          &fexp,&
                          &feyp,&
                          &fezp,&
                          &gamma_ener,&
                          &u,&
                          &v,&
                          &w,&
                          &temp,&
                          &temp_ant,&
                          &fuente_con_temp,&
                          &fuente_lin_temp,&
                          &dt,&
                          &rel_ener,&
                          &a1,b1,c1,r1,&
                          &jj,ii,kk&
                          &)
                  end do
               end do
            end do ensa_ener_dir_y
            !$OMP END PARALLEL DO
            !
            !-------------------------
            !
            ! Condiciones de frontera direcci\'on y
            !
            !$acc parallel loop vector !async(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED)
            cond_fron_ener_direc_y: do kk = 2, lk
               do ii = 2, mi
                  !***********************
                  !Condiciones de frontera
                  a1(indeyp(1,ii,kk))    = 0.0_DBL
                  b1(indeyp(1,ii,kk))    =-1.0_DBL
                  c1(indeyp(1,ii,kk))    = 1.0_DBL
                  r1(indeyp(1,ii,kk))    = 0.0_DBL
                  !
                  a1(indeyp(nj+1,ii,kk)) =-1.0_DBL
                  b1(indeyp(nj+1,ii,kk)) = 1.0_DBL
                  c1(indeyp(nj+1,ii,kk)) = 0.0_DBL
                  r1(indeyp(nj+1,ii,kk)) = 0.0_DBL
                  !
               end do
            end do cond_fron_ener_direc_y
            !$OMP END PARALLEL DO
            !
            !----------------------------------------------
            !
            ! Soluci\'on de la ec. de la energ\'ia
            ! en direcci\'on y
            !
            !$acc parallel loop gang async(stream1) wait(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
            sol_ener_dir_y: do kk = 2, lk
               do ii = 2, mi
                  !
                  call tridiagonal(&
                       &a1(indeyp(1,ii,kk):indeyp(nj+1,ii,kk)),&
                       &b1(indeyp(1,ii,kk):indeyp(nj+1,ii,kk)),&
                       &c1(indeyp(1,ii,kk):indeyp(nj+1,ii,kk)),&
                       &r1(indeyp(1,ii,kk):indeyp(nj+1,ii,kk)),&
                       &nj+1)
                  !
               end do
            end do sol_ener_dir_y
            !$OMP END PARALLEL DO
            !
            !----------------------------------
            !
            ! Actualizaci\'on de la temperatura
            ! en direcci\'on z
            !
            !$acc parallel loop gang collapse(2) !async(stream2) wait(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(3)
            do kk = 2, lk
               do ii = 2, mi
                  do jj = 1, nj+1
                     temp(ii,jj,kk) = r1(indeyp(jj,ii,kk))
                  end do
               end do
            end do
            !$OMP END PARALLEL DO
            !
            !----------------------------------------
            !
            ! Se ensamblan las matrices tridiagonales
            ! en la direcci'on de z para la energía
            !
            !$acc parallel loop gang !async(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
            ensa_ener_dir_z: do ii = 2, mi
               do jj = 2, nj
                  do kk = 2, lk
                     call ensambla_energia_z(&
                          &deltaxp,&
                          &deltayp,&
                          &deltazp,&
                          &deltaxu,&
                          &deltayv,&
                          &deltazw,&
                          &fexp,&
                          &feyp,&
                          &fezp,&
                          &gamma_ener,&
                          &u,&
                          &v,&
                          &w,&
                          &temp,&
                          &temp_ant,&
                          &fuente_con_temp,&
                          &fuente_lin_temp,&
                          &dt,&
                          &rel_ener,&
                          &a1,b1,c1,r1,&
                          &kk,jj,ii&
                          &)
                  end do
               end do
            end do ensa_ener_dir_z
            !$OMP END PARALLEL DO
            !
            !-------------------------
            !
            ! Condiciones de frontera direcci\'on z
            !
            !$acc parallel loop vector !async(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED)
            cond_fron_ener_direc_z: do ii = 2, mi
               do jj = 2, nj
                  !***********************
                  !Condiciones de frontera
                  a1(indezp(1,jj,ii))    = 0.0_DBL
                  b1(indezp(1,jj,ii))    = 1.0_DBL !-1.0_DBL
                  c1(indezp(1,jj,ii))    = 0.0_DBL !1.0_DBL
                  r1(indezp(1,jj,ii))    = 0.0_DBL
                  !
                  a1(indezp(lk+1,jj,ii)) =-1.0_DBL !0.0_DBL
                  b1(indezp(lk+1,jj,ii)) = 1.0_DBL
                  c1(indezp(lk+1,jj,ii)) = 0.0_DBL
                  r1(indezp(lk+1,jj,ii)) = 0.0_DBL
                  !
               end do
            end do cond_fron_ener_direc_z
            !$OMP END PARALLEL DO
            !
            !----------------------------------------------
            !
            ! Soluci\'on de la ec. de la energ\'ia
            ! en direcci\'on z
            !
            !$acc parallel loop gang async(stream1) wait(stream2)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
            sol_ener_dir_z: do ii = 2, mi
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
            end do sol_ener_dir_z
            !$OMP END PARALLEL DO
            !
            !----------------------------------
            !
            ! Actualizaci\'on de la temperatura
            ! en direcci\'on z
            !
            !$acc parallel loop gang collapse(2) !async(stream2) wait(stream1)
            !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(3)
            do ii = 2, mi
               do jj = 2, nj
                  do kk = 1, lk+1
                     temp(ii,jj,kk) = r1(indezp(kk,jj,ii))
                  end do
               end do
            end do
            !$OMP END PARALLEL DO
            !
            error = 0.0_DBL
            !$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:error)
            calcula_ftemp: do kk = 1, lk+1
               do jj = 1, nj+1
                  do ii = 1, mi+1
                     error = error + (ftemp(ii,jj,kk)-temp(ii,jj,kk))*&
                          &(ftemp(ii,jj,kk)-temp(ii,jj,kk))
                  end do
               end do
            end do calcula_ftemp
            !$OMP END PARALLEL DO
            error = dsqrt(error)
            !
            !*************************************
            ! Criterio de convergencia temperatura
            ! WRITE(*,*) 'temp', error
            IF(error<conv_t)EXIT
            !
         END DO
         !
         !
         !---------------------------------------------
         !---------------------------------------------
         !
         ! Criterios de convergencia del paso de tiempo
         !
         !---------------------------------------------
         !---------------------------------------------
         !
         !$acc parallel !async(stream1)
         !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)
         bucle_direccion_z: do kk = 2, lk
            !
            bucle_direccion_y: do jj = 2, nj
               !
               !$acc loop vector
               bucle_direccion_x: do ii = 2, mi-1
                  call residuo_u(&
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
                       &ii,jj,kk&
                       &)
               end do bucle_direccion_x
               !
            end do bucle_direccion_y
            !
         end do bucle_direccion_z
         !$OMP END PARALLEL DO
         !$acc end parallel
         !
         !-------------------------
         !
         ! residuo del algoritmo
         !
         maxbo   = 0.0_DBL
         !$OMP PARALLEL DO REDUCTION(+:maxbo)
         calculo_maxbo: do ii = 2, mi
            do jj = 2, nj
               do kk = 2, lk
                  maxbo = maxbo + b_o(kk,jj,ii)*b_o(kk,jj,ii)
               end do
            end do
         end do calculo_maxbo
         !$OMP END PARALLEL DO
         !
         maxbo =dsqrt(maxbo)
         !
         residuo = 0.0_DBL
         !$acc parallel loop reduction(+:residuo) !async(stream1)
         !$OMP PARALLEL DO REDUCTION(+:residuo)
         calculo_residuou: do kk = 2, lk
            do jj = 2, nj
               do ii = 2, mi-1
                  residuo = residuo + r1(indexu(ii,jj,kk))*r1(indexu(ii,jj,kk))
               end do
            end do
         end do calculo_residuou
         !$OMP END PARALLEL DO
         !
         residuo =dsqrt(residuo)
         !
         write(102,*) 'SIMPLE', iter_simp, maxbo, residuo
         IF( maxbo<conv_paso .and. residuo < conv_resi)EXIT
         !*************************************************
      END DO ALGORITMO_SIMPLE  !final del algoritmo SIMPLE
      !
      ! Mensaje de convergencia
      !
      WRITE(*,*) 'tiempo ',itera,iter_simp,maxbo,residuo
      itera = itera + 1
      !
      !*********************************
      tiempo   = tiempo_inicial+itera*dt
      !
      !-------------------------------------------------------
      !
      ! Se actualizan los arreglos del paso de tiempo anterior
      !
      !-------------------------------------------------------
      !
      !$OMP PARALLEL DO
      actualiza_temp: do kk = 2, lk
         do jj = 2, nj
            do ii = 2, mi
               temp_ant(ii,jj,kk) = temp(ii,jj,kk)
            end do
         end do
      end do actualiza_temp
      !$OMP END PARALLEL DO
      !
      !$OMP PARALLEL DO
      actualiza_u: do kk = 2, lk
         do jj = 2, nj
            do ii = 2, mi-1
               u_ant(ii,jj,kk) = u(ii,jj,kk)
            end do
         end do
      end do actualiza_u
      !$OMP END PARALLEL DO
      !
      !$OMP PARALLEL DO
      actualiza_v: do kk = 2, lk
         do jj = 2, nj-1
            do ii = 2, mi
               v_ant(ii,jj,kk) = v(ii,jj,kk)
            end do
         end do
      end do actualiza_v
      !$OMP END PARALLEL DO
      !
      !$OMP PARALLEL DO
      actualiza_w: do kk = 2, lk-1
         do jj = 2, nj
            do ii = 2, mi
               w_ant(ii,jj,kk) = w(ii,jj,kk)
            end do
         end do
      end do actualiza_w
      !$OMP END PARALLEL DO
      !
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
CALL postprocess_vtk(xp,yp,zp,uf,vf,wf,pres,temp,b_o,archivo)
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
