PROGRAM SIMPLE2D
  use malla
  use ensamblaje
  use solucionador
  use residuos
  use postproceso
  ! USE mkl95_LAPACK
  IMPLICIT NONE
  INCLUDE 'omp_lib.h'
  INTEGER :: i,j,k,l,itera_total,itera,itera_inicial,i_1,paq_itera
  integer :: iter_ecuaci, iter_ecuaci_max
  integer :: iter_simple, iter_simple_max
  INTEGER :: millar,centena,decena,unidad,decima,id,nthreads
  integer :: ii,jj, iter, auxiliar
  integer :: stream1 = 1, stream2 = 2
  ! ------------------------------------------------------------
  !
  ! Variables para los archivos de la entrada de datos
  !
  CHARACTER(len=22) :: entrada_u,entrada_v,entrada_tp
  !*******************************************
  ! Variables del flujo,entropia,nusselt e inc'ognitas,residuos,relajaci'on y convergencia
  REAL(kind=DBL), DIMENSION(mi,nj+1)   ::  u,u_ant,du,au,Resu,gamma_u,Ri,au_aux,fu
  REAL(kind=DBL), DIMENSION(mi+1,nj)   ::  v,v_ant,dv,av,Resv,gamma_v,fv
  REAL(kind=DBL), DIMENSION(mi+1,nj+1) ::  temp,temp_ant,dtemp,Restemp,pres,corr_pres,dcorr_pres
  REAL(kind=DBL), DIMENSION(mi+1,nj+1) ::  fcorr_pres,ftemp
  REAL(kind=DBL), DIMENSION(mi+1,nj+1) ::  entropia_calor,entropia_viscosa,entropia,uf,vf,b_o,gamma_t
  REAL(kind=DBL) :: temp_med,nusselt0,nusselt1,entropia_int,temp_int,gamma_s,residuo,error,maxbo
  REAL(kind=DBL) :: conv_u,conv_p,conv_t,conv_resi,conv_paso,rel_pres,rel_vel,rel_ener
  !
  ! Variables de la malla, volumen de control y factores de interpolaci\'on
  !
  REAL(kind=DBL), DIMENSION(mi)   ::  xu
  REAL(kind=DBL), DIMENSION(nj)   ::  yv
  REAL(kind=DBL), DIMENSION(mi+1) ::  xp
  REAL(kind=DBL), DIMENSION(nj+1) ::  yp
  real(kind=DBL), dimension(mi)   ::  deltaxp
  real(kind=DBL), dimension(nj)   ::  deltayp
  real(kind=DBL), dimension(mi)   ::  deltaxu
  real(kind=DBL), dimension(nj)   ::  deltayu
  real(kind=DBL), dimension(mi)   ::  deltaxv   
  real(kind=DBL), dimension(nj)   ::  deltayv
  REAL(kind=DBL), DIMENSION(mi)   ::  fexp
  REAL(kind=DBL), DIMENSION(nj)   ::  feyp
  REAL(kind=DBL), DIMENSION(mi)   ::  fexu
  REAL(kind=DBL), DIMENSION(nj)   ::  feyv
  !!!!!!!!!!!!!!!!!!1
  REAL(kind=DBL), DIMENSION(mi-1) :: d_xu,d2_xu
  REAL(kind=DBL), DIMENSION(nj-1) :: d_yv,d2_yv
  !!!!!!!!!!!!!!!!!!!
  !
  ! Coeficientes de difusi\'on para las ecuaciones de energia
  ! y balance de momento
  !
  real(kind=DBL), dimension(mi+1,nj+1)  ::  gamma_momen, Riy
  real(kind=DBL), dimension(mi+1,nj+1)  ::  gamma_energ
  !
  ! Tamanio del dominio y ubicaci\'on de las fuentes de calor
  !
  REAL(kind=DBL)   :: ao
  integer          :: placa_min, placa_max
  !
  ! Coeficientes para las matrices 
  !
  real(kind=DBL), dimension(mi+1,nj+1) :: AI, AC, AD, Rx
  real(kind=DBL), dimension(nj+1,mi+1) :: BS, BC, BN, Ry
  !
  !
  !
  REAL(kind=DBL)   :: tiempo,tiempo_inicial,dt,Ra,Pr,Ri_1
  REAL(kind=DBL)   :: a_ent,lambda_ent
  CHARACTER(len=1) :: dec,un,de,ce,m
  CHARACTER(len=3) :: njc,mic,Rec
  CHARACTER(len=5) :: sample

  CHARACTER(len=46):: archivo=repeat(' ',46)
  LOGICAL          :: res_fluido_u
  !****************************************
  !Variables de caracterizaci'on del fluido
  REAL(kind=DBL) :: temp_ref,visc_cin,dif_term,cond_ter,cons_gra,coef_exp,long_ref,dens_ref
  !****************************
  !declaraci´on de variable DBL
  REAL(kind=DBL) :: var2=0.0_DBL
  !*******************************************
  !Se muestra cu'antos procesadores hay en uso
  !$OMP PARALLEL private(id)
!!$  id = omp_get_thread_num()
!!$  WRITE(*,*) 'Este es el thread no. ', id
!!$  !$OMP BARRIER
!!$  IF(id == 0)THEN
!!$     nthreads = omp_get_num_threads()
!!$     WRITE(*,*) 'Se usan ', nthreads, ' threads'
!!$  END IF
  !$OMP END PARALLEL
  !*************************************
  ! Par'ametros para Convecci'on mixta
  OPEN(unit=10,file='parametros.dat')
  READ (10,*) Ra          ! n'umero de Reynolds
  READ (10,*) Rec         ! caracter de Re
  READ (10,*) Pr          ! n'umero de Prandtl
  READ (10,*) Ri_1        ! n'umero de Richardson
  READ (10,*) dt          ! incremento de tiempo
  READ (10,*) paq_itera   ! paquete de iteraciones
  READ (10,*) rel_pres    ! relajaci'on de la presi'on
  READ (10,*) rel_vel     ! relajaci'on de la velocidad
  READ (10,*) rel_ener    ! relajaci'on de la temperatura
  READ (10,*) conv_u      ! convergencia de la velocidad
  READ (10,*) conv_t      ! convergencia de la temperatura
  READ (10,*) conv_p      ! convergencia de la presi'on
  READ (10,*) conv_resi   ! convergencia del residuo
  READ (10,*) conv_paso   ! convergencia del paso de tiempo
  read (10,*) iter_ecuaci_max ! iter m\'aximas para las ecuaciones
  read (10,*) iter_simple_max ! iter m\'aximas algoritmo SIMPLE
  READ (10,*) entrada_u   ! archivo de entrada para u
  READ (10,*) entrada_v   ! archivo de entrada para v
  READ (10,*) entrada_tp  ! archivo de entrada para t y p
  CLOSE(unit=10)
  IF(nj < 100)THEN
     WRITE(njc,170) int(nj);170 format(I2)
     njc = '0'//njc
  ELSE
     WRITE(njc,160) int(nj);160 format(I3)
  ENDIF
  IF(mi < 100)THEN
     WRITE(mic,170) int(mi)
     mic = '0'//mic
  ELSE
     WRITE(mic,160) int(mi)
  ENDIF
  gamma_momen = 1.0_DBL/(Ra)
  !gamma_s = 10._DBL*(1._DBL/(Re*Pr))
  gamma_energ = 1.0_DBL/(Ra*Pr)
  
  gamma_t = 1._DBL/(Ra*Pr) !sqrt(1._DBL/(Pr*Ra))
  gamma_u = 1._DBL/(Ra)    !sqrt(Pr/Ra)
  gamma_v = 1._DBL/(Ra)    !sqrt(Pr/Ra)
  Ri      = Ri_1
  Riy     = 0._DBL
  !
  ! Lectura de las mallas escalonadas e inicializaci\'on de arreglos
  !
  call lectura_mallas_escalonadas(entrada_u,entrada_v,entrada_tp,&
       &u_ant,v_ant,pres,temp_ant,&
       &xp,yp,xu,yv,&
       &deltaxp,deltayp,&
       &deltaxu,deltayu,&
       &deltaxv,deltayv,&
       &fexp,feyp,fexu,feyv,&
       &ao,placa_min,placa_max,itera_inicial)
  !------------------------------
  !
  ! incrementos (codigo temporal)
  !
  d_xu(1) = xu(2)-xu(1)
  DO i = 2, mi-1
     d_xu(i)  = xu(i+1)-xu(i)        ! deltaxp(ii+1)
     d2_xu(i) = xu(i+1)-xu(i-1)
  END DO
  d_yv(1) = yv(2)-yv(1)
  DO j = 2, nj-1
     d_yv(j)  = yv(j+1)-yv(j)        !deltayp(jj+1)
     d2_yv(j) = yv(j+1)-yv(j-1)
  END DO
  !*****************
  !valores iniciales
  tiempo_inicial = itera_inicial*dt
  u_ant = 1.0_DBL
  ! v_ant = 0.0_DBL
  ! temp_ant = 0.0_DBL
  u     = u_ant
  v     = v_ant
  uf    = 0.0_DBL
  vf    = 0.0_DBL
  temp  = temp_ant
  Resu      = 0.0_DBL
  corr_pres = 0.0_DBL * dfloat(mi+1-i)/dfloat(mi)
  au    = 1._DBL
  av    = 1._DBL
  b_o   = 0.0_DBL
  itera = 0
  iter_ecuaci = 0
  iter_simple = 0
  !************************************************
  !escribe las caracter´isticas de las variable DBL
  WRITE(*,100) 'Doble',KIND(var2),PRECISION(var2),RANGE(var2)
100 FORMAT(1X,A,': kind= ',I2,', Precision= ',I2,' Rango= ',I3)
  WRITE(*,*)' '
  !escribe informaci'on de los parametros usados
  WRITE(*,101) Ra,Pr,Ri_1,rel_pres,rel_vel
  WRITE(*,102) itera_inicial,mi,nj
  WRITE(*,106) lambda_ent,a_ent
  WRITE(*,*)' '
101 FORMAT(1X,'Re=',F12.3,', Pr=',F8.3', Ri=',F8.3', rel_pres=',F8.3', rel_vel=',F8.3)
102 FORMAT(1X,'Iteracion inicial=',I7,', mi=',I3,', nj=',I3)
106 FORMAT(1X,'No. de Eckert=',F13.10,', a_ent=',F15.3)
  !--------------------------------------------
  !
  ! Inicio del repetidor principal
  !
  do l=1,itermax/paq_itera
     !
     ! Inicio del paquete de iteraciones
     !
     !------------------------------------------
     !
     ! Apertura de la regi\'on de datos paralela
     !
     !$acc data copy(&
     !$acc &         u(1:mi,1:nj+1),v(1:mi+1,1:nj),                            &
     !$acc &         u_ant(1:mi,1:nj+1),v_ant(1:mi+1,1:nj),                    &
     !$acc &         temp_ant(1:mi+1,1:nj+1),                                  &
     !$acc &         pres(1:mi+1,1:nj+1),temp(1:mi+1,1:nj+1),                  &
     !$acc &         corr_pres(1:mi+1,1:nj+1),Resu(1:mi,1:nj+1),               &
     !$acc &         au(1:mi,1:nj+1),av(1:mi+1,1:nj),b_o(1:mi+1,1:nj+1),       &
     !$acc &         gamma_momen(1:mi+1,1:nj+1),gamma_energ(1:mi+1,1:nj+1),    &
     !$acc &         tiempo_inicial,nusselt0,nusselt1                          &
     !$acc &         )                                                         &
     !$acc & copyin(&
     !$acc &        deltaxp(1:mi),deltayp(1:nj),deltaxu(1:mi),deltayu(1:nj),   &
     !$acc &        deltaxv(1:mi),deltayv(1:nj),                               &
     !$acc &        fexp(1:mi),feyp(1:nj),fexu(1:mi),feyv(1:nj),               &
     !$acc &        Ri(1:mi,1:nj+1),Riy(1:mi+1,1:nj+1),dt,                     &
     !$acc &        rel_vel,conv_u,conv_p,                                     &
     !$acc &        rel_ener,conv_t,placa_min,placa_max                        &
     !$acc &        )&
     !$acc & create(AI(1:mi+1,1:nj+1),AC(1:mi+1,1:nj+1),                       &
     !$acc &        AD(1:mi+1,1:nj+1),Rx(1:mi+1,1:nj+1),                       &
     !$acc &        BS(1:nj+1,1:mi+1),BC(1:nj+1,1:mi+1),                       &
     !$acc &        BN(1:nj+1,1:mi+1),Ry(1:nj+1,1:mi+1),                       &
     !$acc &        fu(1:mi,1:nj+1),du(1:mi,1:nj+1),                           &
     !$acc &        fv(1:mi+1,1:nj),dv(1:mi+1,1:nj),                           &
     !$acc &        fcorr_pres(1:mi+1,1:nj+1),dcorr_pres(1:mi+1,1:nj+1),       &
     !$acc &        ftemp(1:mi+1,1:nj+1),dtemp(1:mi+1,1:nj+1)                  &
     !$acc &)
     !
     do k=1,paq_itera
        !
        ! Inicio del algoritmo SIMPLE
        !
        ALGORITMO_SIMPLE: do
           !--------------------------------
           !
           ! residuo del algoritmo
           residuo = 0.0_DBL
           !
           ecuacion_momento: do
              !
              ! error de la ecuacion de momento
              error = 0.0_DBL
              !
              !$acc parallel loop gang
              inicializacion_fu: do jj=2, nj
                 !
                 !$acc loop vector
                 do ii = 2, mi
                    fu(ii,jj) = u(ii,jj)
                    fv(ii,jj) = v(ii,jj)
                 end do
              end do inicializacion_fu
              !
              !$acc parallel
              call ensambla_velu(deltaxu,deltayu,deltaxp,&
                   &deltayv,fexp,feyp,fexu,gamma_momen,&
                   &u,u_ant,v,&
                   &temp,pres,Ri,dt,rel_vel,&
                   &AI,AC,AD,Rx,BS,BC,BN,Ry,au&
                   &)
              !$acc end parallel
              !
              ! Soluci\'on del sistema de ecuaciones
              !
              !$acc parallel loop gang async(stream1)
              !
              ! $acc loop gang
              solucion_momento_ux: do jj = 2, nj
                 
                 call tridiagonal(AI(1:mi,jj),AC(1:mi,jj),AD(1:mi,jj),Rx(1:mi,jj),mi)
                 
              end do solucion_momento_ux
              !
              !$acc parallel loop gang  async(stream2)            
              solucion_momento_uy: do ii = 2, mi-1
                 
                 call tridiagonal(BS(1:nj+1,ii),BC(1:nj+1,ii),BN(1:nj+1,ii),Ry(1:nj+1,ii),nj+1)

              end do solucion_momento_uy
              !
              ! $acc end parallel
              !
              !$acc wait
              !
              ! Actualizaci\'on de la velocidad u
              !
              !$acc parallel
              !
              !$acc loop vector
              do jj= 2, nj
                 u(1, jj) = Rx(1,jj)
                 u(mi,jj) = Rx(mi,jj)
              end do
              !
              !$acc loop vector
              do ii = 2, mi-1
                 u(ii,1)    = Ry(1,ii)
                 u(ii,nj+1) = Ry(nj+1,ii)
              end do
              !
              !$acc loop gang
              do jj = 2, nj
                 !
                 !$acc loop vector
                 do ii = 2, mi-1
                    u(ii,jj) = 0.5_DBL*Rx(ii,jj)+0.5_DBL*Ry(jj,ii)
                 end do
              end do
              !
              !$acc end parallel 
              !
              !$acc parallel
              call ensambla_velv(deltaxv,deltayv,deltaxu,&
                   &deltayp,fexp,feyp,feyv,gamma_momen,&
                   &v,v_ant,u,&
                   &temp,pres,Riy,dt,rel_vel,&
                   &AI,AC,AD,Rx,BS,BC,BN,Ry,av&
                   &)
              !$acc end parallel
              !
              ! Soluci\'on de las ecs. de momento
              !
              !$acc parallel loop gang async(stream1)
              solucion_momento_vx: do jj = 2, nj-1

                 call tridiagonal(AI(1:mi+1,jj),AC(1:mi+1,jj),AD(1:mi+1,jj),Rx(1:mi+1,jj),mi+1)

              end do solucion_momento_vx
              !
              !$acc parallel loop gang async(stream2)
              solucion_momento_vy: do ii = 2, mi

                 call tridiagonal(BS(1:nj,ii),BC(1:nj,ii),BN(1:nj,ii),Ry(1:nj,ii),nj)

              end do solucion_momento_vy
              !$acc wait
              !
              ! actualizaci\'on de la velocidad v
              !
              !$acc parallel
              !
              !$acc loop vector
              do jj= 2, nj-1
                 v(1,jj)    = Rx(1,jj)
                 v(mi+1,jj) = Rx(mi+1,jj)
              end do
              !
              !$acc loop vector
              do ii = 2, mi
                 v(ii,1)    = Ry(1,ii)
                 v(ii,nj) = Ry(nj,ii)
              end do
              !
              !$acc loop gang
              do jj = 2, nj-1
                 !
                 !$acc loop vector
                 do ii = 1, mi
                    v(ii,jj) = 0.5_DBL*Rx(ii,jj)+0.5_DBL*Ry(jj,ii)
                 end do
              end do
              !
              !$acc end parallel             
              !
              !$acc parallel loop reduction(max:error) async(stream1)
              calculo_diferencias_dv: do jj=2, nj
                 !
                 !$acc loop
                 do ii = 2, mi
                    error = max(error, dabs(u(ii,jj)-fu(ii,jj)))
                    error = max(error, dabs(v(ii,jj)-fv(ii,jj)))
                    ! du(ii,jj) = u(ii,jj)-fu(ii,jj)
                    ! dv(ii,jj) = v(ii,jj)-fv(ii,jj)
                 end do
              end do calculo_diferencias_dv
              !
              !Criterio de convergencia de la velocidad
              !
              if ( error < conv_u ) then
                 iter_ecuaci = 0
                 exit
              else if (iter_ecuaci > iter_ecuaci_max) then
                 iter_ecuaci = 0
                 write(*,*) "ADVER. MOMEN: convergencia no alcanzada ", error
                 exit
              else
                 iter_ecuaci = iter_ecuaci+1
                 ! write(*,*) 'velocidad ',k,MAXVAL(DABS(du)),MAXVAL(DABS(dv))  
              end if
              !            
           end do ecuacion_momento
           !-----------------------------------------
           !-----------------------------------------
           !
           ! Se calcula la correcci'on de la presi'on
           !
           !$acc parallel loop gang
           inicializa_corrector_presion: do jj = 1, nj+1
              !
              !$acc loop vector
              do ii = 1, mi+1
                 corr_pres(ii,jj) = 0.0_DBL
                 fcorr_pres(ii,jj)= 0.0_DBL
              end do
           end do inicializa_corrector_presion
           !
           correccion_presion: do
              !
              error = 0.0_DBL
              maxbo = 0.0_DBL
              !
              !$acc parallel loop gang
              inicializa_fcorr_press: do jj=2, nj
                 !
                 !$acc loop vector
                 do ii = 2, mi
                    fcorr_pres(ii,jj) = corr_pres(ii,jj)
                 end do
              end do inicializa_fcorr_press
              !
              ! $acc parallel loop gang
              ! barrido_direccion_yp: do jj=1, nj
              !
              !$acc parallel
              call ensambla_corr_pres(deltaxp,deltayp,&
                   &deltaxu,deltayv,&
                   &u,v,b_o,&
                   &corr_pres,0.85_DBL,&
                   &AI,AC,AD,Rx,BS,BC,BN,Ry,au,av&
                   &)
              !$acc end parallel
              !
              !-------------------------------------------------------
              !
              ! Soluci\'on de la correcci\'on de la presi\'on
              !
              !$acc parallel loop gang async(stream1)
              solucion_presion_x: do jj = 2, nj

                 call tridiagonal(AI(1:mi+1,jj),AC(1:mi+1,jj),AD(1:mi+1,jj),Rx(1:mi+1,jj),mi+1)
  
              end do solucion_presion_x
              !
              !$acc parallel loop gang async(stream2)
              solucion_presion_y: do ii = 2, mi

                 call tridiagonal(BS(1:nj+1,ii),BC(1:nj+1,ii),BN(1:nj+1,ii),Ry(1:nj+1,ii),nj+1)
  
              end do solucion_presion_y
              !$acc wait
              !
              ! Actualizaci\'on del corrector de la presi\'on
              !
              !$acc parallel
              !
              !$acc loop vector
              do jj= 2, nj
                 corr_pres(1,jj)    = Rx(1,jj)
                 corr_pres(mi+1,jj) = Rx(mi+1,jj)
              end do
              !
              !$acc loop vector
              do ii = 2, mi
                 corr_pres(ii,1)    = Ry(1,ii)
                 corr_pres(ii,nj+1) = Ry(nj+1,ii)
              end do
              !
              !$acc loop gang
              do jj = 2, nj
                 !
                 !$acc loop vector
                 do ii = 1, mi
                    corr_pres(ii,jj) = 0.5_DBL*Rx(ii,jj)+0.5_DBL*Ry(jj,ii)
                 end do
              end do
              !
              !$acc end parallel             
              !              
              !
              ! C\'alculo de diferencias y criterio de convergencia
              !
              !$acc parallel loop gang reduction(max:error,maxbo)
              calculo_dif_corr_pres: do jj=1, nj+1
                 !
                 !$acc loop vector
                 do ii=1, mi+1
                    !dcorr_pres(ii,jj) = (corr_pres(ii,jj)-fcorr_pres(ii,jj))
                    error = max(error,dabs(corr_pres(ii,jj)-fcorr_pres(ii,jj)))
                    maxbo = max(maxbo,dabs(b_o(ii,jj)))
                 end do
              end do calculo_dif_corr_pres
              !-----------------------------------------------------
              !
              ! Critero de convergencia del corrector de la presi'on
              !
              if( error<conv_p )then
                 iter_ecuaci = 0
                 exit
              else if (iter_ecuaci > iter_ecuaci_max) then
                 iter_ecuaci = 0
                 write(*,*) "ADVER. PRES: convergencia no alcanzada ", &
                      error
                 exit
              else
                 iter_ecuaci = iter_ecuaci+1
                 ! write(*,*) 'corrector presion ', MAXVAL(DABS(dcorr_pres)), MAXVAL(DABS(b_o))
              end if
           end do correccion_presion
           !--------------------------------------------
           !
           ! Se corrige la presion
           !
           !$acc parallel
           !
           !$acc loop gang
           do jj = 2, nj
              !$acc loop vector
              do ii = 2, mi
                 pres(ii,jj) = pres(ii,jj) + corr_pres(ii,jj)
              end do
           end do
           !
           !---------------------------------
           !
           ! Se corrigen las velocidades
           !
           !$acc loop gang
           do jj = 2, nj-1
              !
              !$acc loop vector
              do ii = 2, mi-1
                 u(ii,jj) = u(ii,jj)+deltayu(jj)*(corr_pres(ii,jj)-corr_pres(ii+1,jj))/au(ii,jj)
                 v(ii,jj) = v(ii,jj)+deltaxv(ii)*(corr_pres(ii,jj)-corr_pres(ii,jj+1))/av(ii,jj)
              end do
              v(mi,jj) = v(mi,jj)+deltaxv(mi)*(corr_pres(mi,jj)-corr_pres(mi,jj+1))/av(mi,jj)
           end do
           !
           !$acc loop vector
           do ii = 2, mi-1
              u(ii,nj) = u(ii,nj)+deltayu(nj)*(corr_pres(ii,nj)-corr_pres(ii+1,nj))/au(ii,nj)
           end do
           !
           !$acc end parallel
           !
           !--------------------------------------
           !
           ! Se resuelve la ecuaci\'on de balance de energ\'ia
           !
           solucion_energia: do
              !
              error = 0.0_DBL
              !
              !$acc parallel loop gang
              inicializacion_ftemp: do jj=2, nj
                 !
                 !$acc loop vector
                 do ii = 2, mi
                    ftemp(ii,jj) = temp(ii,jj)
                 end do
              end do inicializacion_ftemp
              !
              !$acc parallel
              call ensambla_energia(deltaxp,deltayp,&
                   &deltaxu,deltayv,fexu,feyv,gamma_energ,&
                   &u,v,&
                   &temp,temp_ant,dt,&
                   &rel_ener,placa_min,placa_max,&
                   &AI,AC,AD,Rx,BS,BC,BN,Ry)
              !$acc end parallel
              !
              !
              !$acc parallel loop gang async(stream1)
              solucion_energia_x: do jj = 2, nj
                 
                 call tridiagonal(AI(1:mi+1,jj),AC(1:mi+1,jj),AD(1:mi+1,jj),Rx(1:mi+1,jj),mi+1)

              end do solucion_energia_x
              !
              !$acc parallel loop gang async(stream2)
              solucion_energia_y: do ii = 2, mi

                 call tridiagonal(BS(1:nj+1,ii),BC(1:nj+1,ii),BN(1:nj+1,ii),Ry(1:nj+1,ii),nj+1)

              end do solucion_energia_y
              !$acc wait
              !
              !$acc parallel
              !
              !$acc loop vector
              do jj= 2, nj
                 temp(1,jj)    = Rx(1,jj)
                 temp(mi+1,jj) = Rx(mi+1,jj)
              end do
              !
              !$acc loop vector
              do ii = 2, mi
                 temp(ii,1)    = Ry(1,ii)
                 temp(ii,nj+1) = Ry(nj+1,ii)
              end do
              !
              !$acc loop gang
              do jj = 2, nj
                 !
                 !$acc loop vector
                 do ii = 1, mi
                    temp(ii,jj) = 0.5_DBL*Rx(ii,jj)+0.5_DBL*Ry(jj,ii)
                 end do
              end do
              !
              !$acc end parallel             
              !
              !
              !$acc parallel loop reduction(max:error)
              calculo_diferencias_dtemp: do jj = 2, nj
                 !
                 !$acc loop
                 do ii = 2, mi
                    error = max(error, dabs(temp(ii,jj)-ftemp(ii,jj)))
                    ! dtemp(ii,jj) = temp(ii,jj)-ftemp(ii,jj)
                 end do
              end do calculo_diferencias_dtemp
              !
              ! ------------------------------------------
              !
              ! Criterio de convergencia temperatura
              !
              if( error < conv_t )then
                 iter_ecuaci = 0
                 exit
              else if( iter_ecuaci > iter_ecuaci_max )then
                 iter_ecuaci = 0
                 write(*,*) "Adver. ENERG: convergencia no alcanzada ", &
                      error
                 exit
              else
                 iter_ecuaci = iter_ecuaci + 1
                 ! write(*,*) 'temp', maxval(dabs(dtemp))
              end if
              !
           end do solucion_energia
           !
           !--------------------------------------------
           !
           ! Criterio de convergencia del paso de tiempo
           !
           ! call residuou(res_fluido_u,xu,yp,feyp,d_xu,d2_xu,d_yv,u,&
           !      u_ant,v,temp,pres,Resu,gamma_u,Ri,dt)
           !$acc parallel
           call residuo_u(deltaxu,deltayu,deltaxp,&
                &deltayv,fexp,feyp,fexu,gamma_momen,&
                &u,u_ant,v,&
                &temp,pres,Ri,dt,rel_vel,&
                &Resu&
                &)
           !$acc end parallel
           !
           !$acc parallel loop reduction(max:error)
           calculo_maximo_residuou: do jj=2, nj
              !
              !$acc loop
              do ii = 2, mi-1
                 residuo = max(residuo, dabs(Resu(ii,jj)))
              end do
              !
           end do calculo_maximo_residuou
           !
           if ( residuo < conv_resi .and. maxbo<conv_paso)then
              iter_simple = 0
              exit
           else if ( iter_simple > iter_simple_max ) then
              iter_simple = 0
              write(*,*) "Advertencia SIMPLE: convergencia no alcanzada, ", residuo,maxbo   
              exit
           else
              iter_simple = iter_simple + 1
              ! write(*,*) 'tiempo ',itera,res_fluido_u,MAXVAL(DABS(Resu)),MAXVAL(DABS(b_o)),&
                   ! &MAXVAL(DABS(pres))
           end if
           !
        end do ALGORITMO_SIMPLE  !final del algoritmo SIMPLE
        !
        ! Escritura de mensajes y postprocesos
        !
        itera = itera + 1
        if( mod(itera,500)==0 )write(*,*) 'ITERACION: ',itera,maxval(dabs(b_o)),&
             &maxval(dabs(Resu))

        if( mod(itera,100)==0 )then
           !     CALL entropia_cvt(x,y,u,xu,v,yv,temp,entropia_calor,entropia_viscosa,entropia,&
           !&entropia_int,temp_int,a_ent,lambda_ent)
           ! call nusselt(xp,yp,d_xu,d_yv,temp,nusselt0,nusselt1,placa_min,placa_max)
           !$acc parallel
           call nusselt_promedio_y(&
                &xp,yp,deltaxp,deltayp,&
                &temp,nusselt0,nusselt1,&
                &placa_min,placa_max&
                &)
           !$acc end parallel
           temp_med = (temp((placa_min+placa_max)/2,nj/2+1)+&
                &temp((placa_min+placa_max)/2+1,nj/2+1))/2._DBL
           OPEN(unit = 5,file = 'nuss_sim_n'//njc//'m'//mic//'_R'//Rec//'.dat',access = 'append')
           WRITE(5,form26) tiempo_inicial+itera*dt,nusselt0,&
                &-nusselt1,temp_med,temp_int,entropia_int
           CLOSE(unit = 5)
        ENDIF
        !*********************************
        tiempo   = tiempo_inicial+itera*dt
        !
        !$acc parallel loop gang
        do jj=1, nj+1
           !
           !$acc loop vector
           do ii=1, mi+1
              temp_ant(ii,jj) = temp(ii,jj)
           end do
        end do
        !
        !$acc parallel loop gang
        do jj=1, nj+1
           !
           !$acc loop vector
           do ii=1, mi
              u_ant(ii,jj) = u(ii,jj)
           end do
        end do
        !
        !$acc parallel loop gang
        do jj=1, nj
           !
           !$acc loop vector
           do ii=1, mi+1
              v_ant(ii,jj) = v(ii,jj)
           end do
        end do
     END DO
     !--------------------------------------------
     !
     ! Se cierra la regi\'on paralela de datos
     !
     !$acc end data
     !
     !*************       termina el paquete de iteraciones
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
     DO j = 2, nj
        DO i = 2, mi
           uf(i,j) = (u(i,j)+u(i-1,j))/2._DBL
        END DO
     END DO
     DO i = 2, mi
        DO j = 2, nj
           vf(i,j) = (v(i,j)+v(i,j-1))/2._DBL
        END DO
     END DO
     DO j = 1, nj+1
        uf(1,j)    = u(1,j)
        uf(mi+1,j) = u(mi,j)
     END DO
     DO i = 1, mi+1
        vf(i,1)    = v(i,1)
        vf(i,nj+1) = v(i,nj)
     END DO
     !************************************
     WRITE(*,*) 'itera_total=',itera_total
     WRITE(*,103) nusselt0,nusselt1
     WRITE(*,104) MAXVAL(DABS(b_o)),MAXVAL(DABS(Resu))
     WRITE(*,105) MAXVAL(DABS(Restemp)),MAXVAL(DABS(Resv))
     WRITE(*,*)' '
103  FORMAT(1X,'N_Izq=',D23.15,', N_Der=',D23.15)
104  FORMAT(1X,'b_o  =',D23.15,', Res_u=',D23.15)
105  FORMAT(1X,'Res_T=',D23.15,', Res_v=',D23.15)
     !********************************
     !*** Formato de escritura dat ***
     OPEN(unit=2,file='out_n'//njc//'m'//mic//'_R'//Rec//'u.dat')
     WRITE(2,*) placa_min,placa_max,itera_total,ao
     DO j = 1, nj+1
        DO i = 1, mi
           WRITE(2,form24) xu(i),yp(j),u(i,j)
        END DO
     END DO
     CLOSE(unit=2)
     OPEN(unit=3,file='out_n'//njc//'m'//mic//'_R'//Rec//'v.dat')
     WRITE(3,*) placa_min,placa_max,itera_total,ao
     DO j = 1, nj
        DO i = 1, mi+1
           WRITE(3,form24) xp(i),yv(j),v(i,j)
        END DO
     END DO
     CLOSE(unit=3)
     OPEN(unit=4,file='out_n'//njc//'m'//mic//'_R'//Rec//'p.dat')
     WRITE(4,*) placa_min,placa_max,itera_total,ao
     DO j = 1, nj+1
        DO i = 1, mi+1
           WRITE(4,form25) xp(i),yp(j),temp(i,j),pres(i,j)
        END DO
     END DO
     CLOSE(unit=4)
     ! *************************************
     ! *** Formato escritura VTK ***********
     ! sample  = m//ce//de//un//dec
     archivo = 'n'//njc//'m'//mic//'R'//Rec//'/t_'//m//ce//de//un//dec//'.vtk'
     CALL postprocess_vtk(xp,yp,uf,vf,pres,temp,b_o,archivo)
     !
     ! !************************************
     ! !*** Formato de escritura Tecplot ***
     ! OPEN(unit=1,file='n'//njc//'m'//mic//'R'//Rec//'/t'//m//ce//de//un//dec//'.plt')
     ! WRITE(1,*)'TITLE     = "Tempe" '
     ! WRITE(1,*)'VARIABLES = "x"'
     ! WRITE(1,*)'"y"'
     ! WRITE(1,*)'"Temperature"'        ! variable
     ! WRITE(1,*)'"U"'
     ! WRITE(1,*)'"V"'
     ! WRITE(1,*)'"Pressure"'
     ! WRITE(1,*)'"corr_pres"'
     ! WRITE(1,*)'"dcorr_pres"'
     ! WRITE(1,*)'"Entropy"'
     ! WRITE(1,*)'"Residual"'
     ! !WRITE(1,*)'"XD"'
     ! !WRITE(1,*)'"DTEMP"'
     ! WRITE(1,*)'ZONE T= "Matrix"'
     ! WRITE(1,*)'I=',mi+1,' J=',nj+1,' K=1,F=POINT'
     ! DO j = 1, nj+1
     !   DO i = 1, mi+1
     !     WRITE(1,27) y(j),ao-x(i),temp(i,j),vf(i,j),-uf(i,j),pres(i,j),corr_pres(i,j),dcorr_pres(i,j),b_o(i,j),b_o(i,j)
     !   END DO
     ! END DO
     ! CLOSE(unit=1)
     !************************************
     !*** Formato de escritura Tecplot ***
     ! OPEN(unit=1,file='n'//njc//'m'//mic//'R'//Rec//'/res'//m//ce//de//un//dec//'.plt')
     ! WRITE(1,*)'TITLE     = "Resu" '
     ! WRITE(1,*)'VARIABLES = "x"'
     ! WRITE(1,*)'"y"'
     ! WRITE(1,*)'"U"'
     ! WRITE(1,*)'"Residual"'
     ! !WRITE(1,*)'"XD"'
     ! !WRITE(1,*)'"DTEMP"'
     ! WRITE(1,*)'ZONE T= "Matrix"'
     ! WRITE(1,*)'I=',mi,' J=',nj+1,' K=1,F=POINT'
     ! DO j = 1, nj+1
     !   DO i = 1, mi
     !     WRITE(1,27) y(j),ao-xu(i),-u(i,j),resu(i,j)
     !   END DO
     ! END DO
     ! CLOSE(unit=1)
     !*******************************
     !formatos de escritura y lectura
  end do !***********final del repetidor principal

end program SIMPLE2D
