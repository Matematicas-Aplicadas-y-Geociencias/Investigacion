!
! Programa que implementa el algoritmo SIMPLE para resolver las ecuaciones
! de Navier-Stokes y la energ\'ia
!
PROGRAM SIMPLE2D
  !
  ! Variables de la malla, volumen de control y factores de interpolaci\'on
  !
  use malla, only : mi, nj, DBL
  use malla, only : mic, njc, zkc
  use malla, only : xu, yv, xp, yp
  use malla, only : deltaxp, deltayp, deltaxu
  use malla, only : deltayu, deltaxv, deltayv
  use malla, only : fexp, feyp, fexu, feyv
  use malla, only : form24, form25, form26, form27
  use malla, only : lectura_mallas_escalonadas
  !
  ! Componentes de velocidad, presi\'on
  ! residuos de las ecuaciones de momento y correcci\'on de la presi\'on
  ! coeficientes de difusi\'on y criterios e convergencia
  !
  use ec_momento, only : u, u_ant, du, au, Resu, Ri, Riy, fu
  use ec_momento, only : v, v_ant, dv, av, Resv, fv
  use ec_momento, only : pres, corr_pres, dcorr_pres, fcorr_pres
  use ec_momento, only : uf, vf, b_o
  use ec_momento, only : maxbo, conv_u, conv_p, conv_paso, rel_pres, rel_vel
  use ec_momento, only : gamma_momen, Riy
  !
  ! Rutinas de ensamblaje de la ec. de momento, correcci\'on de presi\'on y residuo
  !
  use ec_momento, only : ensambla_velu, ensambla_velv, ensambla_corr_pres
  use ec_momento, only : residuo_u
  !
  ! Variables para la ecuaci\'on de la energ\'ia
  ! temperatura, coeficiente de difusi\'on y criterios de convergencia
  !
  use ec_energia, only : temp, temp_ant, dtemp, Restemp
  use ec_energia, only : ftemp, gamma_energ, conv_t,rel_ener
  ! 
  ! Rutina de ensamblaje de la ec. de energ\'ia
  !
  use ec_energia, only : ensambla_energia
  !
  ! Rutinas de soluci\'on de ecuaciones
  !
  use solucionador, only : tridiagonal
  !
  ! Subrutinas de postproceso
  ! 
  use postproceso, only  : nusselt_promedio_y, postproceso_vtk
  use postproceso, only  : postproceso_bin, entero_caracter
  !
  ! Variables auxiliares para bucles y n\'umero de iteraciones
  !
  IMPLICIT NONE
  INCLUDE 'omp_lib.h'
  INTEGER :: itera_total,itera,itera_inicial,i_1,paq_itera,itermax
  integer :: iter_ecuaci, iter_ecuaci_max
  integer :: iter_simple, iter_simple_max
  INTEGER :: millar,centena,decena,unidad,decima,id,nthreads
  integer :: ii, jj, kk, ll, iter, auxiliar
  integer :: stream1 = 1, stream2 = 2, stream3 = 3
  ! ------------------------------------------------------------
  !
  ! Variables para los archivos de la entrada de datos
  !
  CHARACTER(len=22) :: entrada_u,entrada_v,entrada_tp
  !*******************************************
  !
  REAL(kind=DBL), DIMENSION(mi+1,nj+1) :: entropia_calor,entropia_viscosa,entropia,gamma_t
  REAL(kind=DBL) :: temp_med,nusselt0,nusselt1,entropia_int,temp_int,gamma_s,residuo,error
  REAL(kind=DBL) :: conv_resi
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
  ! CHARACTER(len=3) :: njc,mic,Rec
  CHARACTER(len=6) :: Rec=repeat(' ',6)
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
      !
    ! Auxiliares de interpolaci\'on
    !
    real(kind=DBL) :: ui, ud, vs, vn
    real(kind=DBL) :: di, dd, ds, dn
    real(kind=DBL) :: gammai, gammad
    real(kind=DBL) :: gammas, gamman
    real(kind=DBL) :: deltax, deltay
    ! real(kind=DBL) :: temp_int
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
  read (10,*) itermax     ! iteraciones totales de la ejecuci\'on
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
  ! IF(nj < 100)THEN
  !    WRITE(njc,170) int(nj);170 format(I2)
  !    njc = '0'//njc
  ! ELSE
  !    WRITE(njc,160) int(nj);160 format(I3)
  ! ENDIF
  ! IF(mi < 100)THEN
  !    WRITE(mic,170) int(mi)
  !    mic = '0'//mic
  ! ELSE
  !    WRITE(mic,160) int(mi)
  ! ENDIF
  mic = entero_caracter(mi)
  njc = entero_caracter(nj)
  Rec = entero_caracter(ceiling(Ra))
  gamma_momen = 1.0_DBL/(Ra)
  !gamma_s = 10._DBL*(1._DBL/(Re*Pr))
  gamma_energ = 1.0_DBL/(Ra*Pr)
  
  ! gamma_t = 1._DBL/(Ra*Pr) !sqrt(1._DBL/(Pr*Ra))
  ! gamma_u = 1._DBL/(Ra)    !sqrt(Pr/Ra)
  ! gamma_v = 1._DBL/(Ra)    !sqrt(Pr/Ra)
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
  ! !*****************
  !valores iniciales
  tiempo_inicial = itera_inicial*dt
  ! u_ant = 1.0_DBL
  ! v_ant = 0.0_DBL
  ! temp_ant = 0.0_DBL
  u     = u_ant
  v     = v_ant
  uf    = 0.0_DBL
  vf    = 0.0_DBL
  temp  = temp_ant
  Resu      = 0.0_DBL
  corr_pres = 0.0_DBL !* dfloat(mi+1-ii)/dfloat(mi)
  au    = 1._DBL
  av    = 1._DBL
  b_o   = 0.0_DBL
  itera = 0
  iter_ecuaci = 0
  iter_simple = 0
  maxbo   = 0.0_DBL
  error   = 0.0_DBL
  residuo = 0.0_DBL
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
  do ll = 1, itermax/paq_itera
     !
     ! Inicio del paquete de iteraciones
     !
     !------------------------------------------
     !
     ! Apertura de la regi\'on de datos paralela
     !
     !$acc data copy(&
     !$acc &         u(1:mi,1:nj+1),v(1:mi+1,1:nj),                            &
     !$acc &         pres(1:mi+1,1:nj+1),temp(1:mi+1,1:nj+1),                  &
     !$acc &         corr_pres(1:mi+1,1:nj+1),                                 &
     !$acc &         u_ant(1:mi,1:nj+1),v_ant(1:mi+1,1:nj),                    &
     !$acc &         temp_ant(1:mi+1,1:nj+1)                                   &
     !$acc &         )&     
     !$acc & copyin(&
     !$acc &        Resu(1:mi,1:nj+1),                                         &
     !$acc &        tiempo_inicial,                                            &
     !$acc &        au(1:mi,1:nj+1),av(1:mi+1,1:nj),b_o(1:mi+1,1:nj+1),        &
     !$acc &        gamma_momen(1:mi+1,1:nj+1),gamma_energ(1:mi+1,1:nj+1),     &
     !$acc &        deltaxp(1:mi),deltayp(1:nj),deltaxu(1:mi),deltayu(1:nj),   &
     !$acc &        deltaxv(1:mi),deltayv(1:nj),                               &
     !$acc &        fexp(1:mi),feyp(1:nj),fexu(1:mi),feyv(1:nj),               &
     !$acc &        Ri(1:mi,1:nj+1),Riy(1:mi+1,1:nj+1),dt,                     &
     !$acc &        rel_vel,conv_u,conv_p,                 &
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
     do kk=1,paq_itera
        !
        ! Inicio del algoritmo SIMPLE
        !
        ALGORITMO_SIMPLE: do
           !
           !-------------------------------------------
           !-------------------------------------------
           !
           ! Soluci\'on de la ecuaci\'on de momento
           !
           !-------------------------------------------
           !-------------------------------------------           
           ecuacion_momento: do
              !
              !$acc parallel loop gang collapse(2) !async(stream1)
              inicializacion_fu: do jj=1, nj
                 do ii = 1, mi
                    fu(ii,jj) = u(ii,jj)
                    fv(ii,jj) = v(ii,jj)
                 end do
              end do inicializacion_fu
              !
              !------------------------------------------
              !
              ! Se ensambla la ecuaci\'on de momento en u
              !
              ! !$acc parallel !async(stream2)
              !
              ! Llenado de la matriz
              !
              !$acc parallel loop gang !async(stream2)
              bucle_direccion_y: do jj = 2, nj
                 !
                 ! Llenado de la matriz
                 !
                 !$acc loop vector
                 bucle_direccion_x: do ii = 2, mi-1
                    call ensambla_velu(deltaxu,deltayu,deltaxp,&
                         &deltayv,fexp,feyp,fexu,gamma_momen,&
                         &u,u_ant,v,&
                         &temp,pres,Ri,dt,rel_vel,&
                         &AI,AC,AD,Rx,BS,BC,BN,Ry,au,&
                         &ii,jj&
                         &)
                    ! $acc end parallel
                 end do bucle_direccion_x
              end do bucle_direccion_y
              ! !***********************
              ! !Condiciones de frontera
              ! AC(mi,jj) = 1._DBL
              ! AI(mi,jj) =-1._DBL !
              ! Rx(mi,jj) = 0.0_DBL
              ! au(mi,jj) = 1.e40_DBL !ACi(mi)
              !
              ! Condiciones de frontera para la direcci\'on y
              !
              !$acc parallel loop vector !async(stream1)
              bucle_direccionx1: do jj = 2, nj
                 !***********************
                 !Condiciones de frontera
                 !***********************
                 AC(1,jj) = 1._DBL
                 AD(1,jj) = 0._DBL
                 Rx(1,jj) = 1._DBL ! *tanh((itera_total-1)*dt/3.0_DBL) 
                 au(1,jj) = 1.e40_DBL !ACi(1)
                 AC(mi,jj) = 1._DBL
                 AI(mi,jj) =-1._DBL !
                 Rx(mi,jj) = 0.0_DBL
                 au(mi,jj) = 1.e40_DBL !ACi(mi)      
              end do bucle_direccionx1
              !
              ! Condiciones de frontera para la direcci\'on y
              !
              !$acc parallel loop vector !async(stream1)
              bucle_direccionx: do ii = 2, mi-1
                 !***********************
                 !Condiciones de frontera
                 BC(1,ii)     = 1._DBL
                 BN(1,ii)     = 0.0_DBL
                 Ry(1,ii)     = 0.0_DBL
                 au(ii,1)     = 1.e40_DBL !ACj(1)
                 BC(nj+1,ii)  = 1._DBL
                 BS(nj+1,ii)  = 0.0_DBL
                 Ry(nj+1,ii)  = 0.0_DBL
                 au(ii,nj+1)  = 1.e40_DBL !ACj(nj+1)       
              end do bucle_direccionx
              !
              !-------------------------------------
              !
              ! Soluci\'on del sistema de ecuaciones
              !
              !$acc parallel loop gang async(stream1) !wait(stream2)
               solucion_momento_ux: do jj = 2, nj
                 
                 call tridiagonal(AI(1:mi,jj),AC(1:mi,jj),AD(1:mi,jj),Rx(1:mi,jj),mi)
                 u(1, jj) = Rx(1,jj)
                 u(mi,jj) = Rx(mi,jj)
                 
              end do solucion_momento_ux
              !
              !$acc parallel loop gang async(stream2)            
              solucion_momento_uy: do ii = 2, mi-1
                 
                 call tridiagonal(BS(1:nj+1,ii),BC(1:nj+1,ii),BN(1:nj+1,ii),Ry(1:nj+1,ii),nj+1)
                 u(ii,1)    = Ry(1,ii)
                 u(ii,nj+1) = Ry(nj+1,ii)

              end do solucion_momento_uy
              !$acc wait
              !
              !----------------------------------
              !
              ! Actualizaci\'on de la velocidad u
              !
              !$acc parallel loop gang collapse(2) !async(stream2) wait(stream1)
              do jj = 2, nj
                 do ii = 2, mi-1
                    u(ii,jj) = 0.5_DBL*Rx(ii,jj)+0.5_DBL*Ry(jj,ii)
                 end do
              end do
              !$acc wait
              !
              !---------------------------
              !
              ! Se ensambla la velocidad v
              !
              !$acc parallel vector_length(64) !async(stream1)
              !
              ! Condiciones de frontera para la direcci\'on y
              !
              !$acc loop vector
              bucle_direccionxv: do ii = 2, mi
                 !***********************
                 !Condiciones de frontera
                 BC(1,ii)     = 1._DBL
                 BN(1,ii)     = 0.0_DBL
                 Ry(1,ii)     = 0.0_DBL
                 av(ii,1)     = 1.e40_DBL !ACj(1)
                 BC(nj,ii)    = 1._DBL
                 BS(nj,ii)    = 0.0_DBL
                 Ry(nj,ii)    = 0.0_DBL
                 av(ii,nj)    = 1.e40_DBL !ACj(nj+1)       
              end do bucle_direccionxv
              !$acc loop vector
              do jj = 2, nj-1
                 !Condiciones de frontera
                 AC(1,jj) = 1._DBL
                 AD(1,jj) = 0._DBL
                 Rx(1,jj) = 0._DBL 
                 av(1,jj) = 1.e40_DBL !ACi(1)
                 !***********************
                 !Condiciones de frontera
                 AC(mi+1,jj) = 1.0_DBL
                 AI(mi+1,jj) = 0.0_DBL !
                 Rx(mi+1,jj) = 0.0_DBL
                 av(mi+1,jj) = 1.e40_DBL !ACi(mi)
              end do
              !$acc end parallel
              !
              !$acc parallel loop gang !async(stream2)
              do jj = 2, nj-1
                 !$acc loop vector
                 do ii = 2, mi
                    call ensambla_velv(deltaxv,deltayv,deltaxu,&
                         &deltayp,fexp,feyp,feyv,gamma_momen,&
                         &v,v_ant,u,&
                         &temp,pres,Riy,dt,rel_vel,&
                         &AI,AC,AD,Rx,BS,BC,BN,Ry,av,&
                         &ii,jj&
                         &)
                 end do
              end do
              ! $acc end parallel
              !
              !------------------------------------
              !
              ! Soluci\'on de las ecs. de momento v
              !
              !$acc parallel loop gang async(stream1)  !wait(stream2)
              solucion_momento_vx: do jj = 2, nj-1

                 call tridiagonal(AI(1:mi+1,jj),AC(1:mi+1,jj),AD(1:mi+1,jj),Rx(1:mi+1,jj),mi+1)
                 v(1,jj)    = Rx(1,jj)
                 v(mi+1,jj) = Rx(mi+1,jj)

              end do solucion_momento_vx
              !
              !$acc parallel loop gang async(stream2)
              solucion_momento_vy: do ii = 2, mi

                 call tridiagonal(BS(1:nj,ii),BC(1:nj,ii),BN(1:nj,ii),Ry(1:nj,ii),nj)
                 v(ii,1)    = Ry(1,ii)
                 v(ii,nj)   = Ry(nj,ii)

              end do solucion_momento_vy
              !$acc wait
              !
              !----------------------------------
              !
              ! Actualizaci\'on de la velocidad v
              !
              !$acc parallel loop gang collapse(2) !async(stream1) wait(stream2)
              do jj = 2, nj-1
                 do ii = 2, mi
                    v(ii,jj) = 0.5_DBL*Rx(ii,jj)+0.5_DBL*Ry(jj,ii)
                 end do
              end do
              !$acc wait
              !
              ! error de la ecuacion de momento
              !
              error = 0.0_DBL
              !$acc parallel loop gang reduction(max:error) ! async(stream1)
              calculo_diferencias_dv: do jj=2, nj-1
                 !
                 !$acc loop vector
                  do ii = 2, mi
                    error = max(error, abs(v(ii,jj)-fv(ii,jj)))
                  end do
              end do calculo_diferencias_dv
              !
              ! Criterio de convergencia de la velocidad
              !
              if ( error < conv_u ) then
                 iter_ecuaci = 0
                 write(101,*) 'velocidad ', error
                 exit
              else if (iter_ecuaci > iter_ecuaci_max) then
                 iter_ecuaci = 0
                 write(*,*) "ADVER. MOMEN: convergencia no alcanzada ", error
                 exit
              else
                 iter_ecuaci = iter_ecuaci+1
                 write(101,*) 'velocidad ', error
              end if
              !            
           end do ecuacion_momento
           !$acc wait
           !-----------------------------------------
           !-----------------------------------------
           !
           ! Se calcula la correcci'on de la presi'on
           !
           !-----------------------------------------
           !-----------------------------------------
           !
           !$acc parallel loop gang collapse(2) !async(stream2)
           inicializa_corrector_presion: do jj = 1, nj+1
              do ii = 1, mi+1
                 corr_pres(ii,jj) = 0.0_DBL
                 fcorr_pres(ii,jj)= 0.0_DBL
              end do
           end do inicializa_corrector_presion
           !
           correccion_presion: do
              !
              !$acc parallel loop gang collapse(2) ! async(stream1) wait(stream2)
              inicializa_fcorr_press: do jj=2, nj
                 do ii = 2, mi
                    fcorr_pres(ii,jj) = corr_pres(ii,jj)
                 end do
              end do inicializa_fcorr_press
              !
              !-------------------------
              !
              ! Condiciones de frontera
              !
              !$acc parallel loop vector !async(stream1)
              bucle_direccionxe: do ii = 2, mi
                 !***********************
                 !Condiciones de frontera
                 BC(1,ii)     = 1._DBL
                 BN(1,ii)     = 0.0_DBL
                 Ry(1,ii)     = 0.0_DBL
                 !
                 BC(nj+1,ii)  = 1._DBL
                 BS(nj+1,ii)  = 0.0_DBL
                 Ry(nj+1,ii)  = 0.0_DBL
              end do bucle_direccionxe
              !
              !$acc parallel loop vector !async(stream1)
              do jj = 2, nj
                 !------------------------
                 ! Condiciones de frontera
                 AC(1,jj) = 1._DBL
                 AD(1,jj) = 0._DBL
                 Rx(1,jj) = 0._DBL
                 !
                 AI(mi+1,jj) = 0.0_DBL
                 AC(mi+1,jj) = 1.0_DBL
                 Rx(mi+1,jj) = 0.0_DBL
              end do
              !------------------------------------------
              !
              ! Se ensambla la ecuaci\'on de la presi\'on
              !
              !$acc parallel loop gang !async(stream2)
              do jj = 2, nj
                 !$acc loop vector
                 do ii = 2, mi
                    call ensambla_corr_pres(deltaxp,deltayp,&
                         &deltaxu,deltayv,&
                         &u,v,b_o,&
                         &corr_pres,0.85_DBL,&
                         &AI,AC,AD,Rx,BS,BC,BN,Ry,au,av,&
                         &ii,jj)
                 end do
              end do
              !
              !----------------------------------------------
              !
              ! Soluci\'on de la correcci\'on de la presi\'on
              !
              !$acc parallel loop gang async(stream1) wait(stream2)
              solucion_presion_x: do jj = 2, nj

                 call tridiagonal(AI(1:mi+1,jj),AC(1:mi+1,jj),AD(1:mi+1,jj),Rx(1:mi+1,jj),mi+1)
                 corr_pres(1,jj)    = Rx(1,jj)
                 corr_pres(mi+1,jj) = Rx(mi+1,jj)
                 
              end do solucion_presion_x
              !
              !$acc parallel loop gang async(stream2)
              solucion_presion_y: do ii = 2, mi

                 call tridiagonal(BS(1:nj+1,ii),BC(1:nj+1,ii),BN(1:nj+1,ii),Ry(1:nj+1,ii),nj+1)
                 corr_pres(ii,1)    = Ry(1,ii)
                 corr_pres(ii,nj+1) = Ry(nj+1,ii)
  
              end do solucion_presion_y
              !$acc wait
              !
              ! Actualizaci\'on del corrector de la presi\'on
              !
              !$acc parallel loop gang collapse(2) !async(stream1) wait(stream2)
              do jj = 2, nj
                 do ii = 2, mi
                    corr_pres(ii,jj) = 0.5_DBL*Rx(ii,jj)+0.5_DBL*Ry(jj,ii)
                 end do
              end do
              !$acc wait
              !
              ! C\'alculo de diferencias y criterio de convergencia
              !
              error = 0.0_DBL
              maxbo = 0.0_DBL
              !
              !$acc parallel loop gang reduction(max:error) !async(stream1)
              calculo_dif_corr_pres: do jj=2, nj
                 do ii=2, mi
                    error = max(error,abs(corr_pres(ii,jj)-fcorr_pres(ii,jj)))
                    ! maxbo = max(maxbo,abs(b_o(ii,jj)))
                 end do
              end do calculo_dif_corr_pres
              !
              !$acc parallel loop gang reduction(max:maxbo) !async(stream2)
              calculo_dif_maxbo: do jj=2, nj
                 do ii=2, mi
                    ! error = max(error,abs(corr_pres(ii,jj)-fcorr_pres(ii,jj)))
                    maxbo = max(maxbo,abs(b_o(ii,jj)))
                 end do
              end do calculo_dif_maxbo
              !-----------------------------------------------------
              !
              ! Critero de convergencia del corrector de la presi'on
              !
              ! $acc wait
              if( error<conv_p )then
                 ! write(*,*) "PRES: convergencia ", error, " con ", iter_ecuaci," iteraciones"
                 iter_ecuaci = 0
                 exit
              else if (iter_ecuaci > iter_ecuaci_max) then
                 iter_ecuaci = 0
                 ! write(*,*) "ADVER. PRES: convergencia no alcanzada ", &
                 !      error
                 exit
              else
                 iter_ecuaci = iter_ecuaci+1
                 ! write(*,*) 'corrector presion ', MAXVAL(ABS(dcorr_pres)), MAXVAL(ABS(b_o))
              end if
           end do correccion_presion
           !$acc wait
           !--------------------------------------------
           !
           ! Se corrige la presion
           !
           !$acc parallel loop gang collapse(2) !async(stream2)
           do jj = 2, nj
              do ii = 2, mi
                 pres(ii,jj) = pres(ii,jj) + corr_pres(ii,jj)
              end do
           end do
           !
           !---------------------------------
           !
           ! Se corrigen las velocidades
           !
           !$acc parallel loop gang collapse(2) !async(stream1)
           do jj = 2, nj-1
              do ii = 2, mi-1
                 u(ii,jj) = u(ii,jj)+deltayu(jj)*(corr_pres(ii,jj)-corr_pres(ii+1,jj))/au(ii,jj)
                 v(ii,jj) = v(ii,jj)+deltaxv(ii)*(corr_pres(ii,jj)-corr_pres(ii,jj+1))/av(ii,jj)
              end do
           end do
           !
           !$acc parallel loop vector !async(stream1)
           do ii = 2, mi-1
              u(ii,nj) = u(ii,nj)+deltayu(nj)*(corr_pres(ii,nj)-corr_pres(ii+1,nj))/au(ii,nj)
           end do
           !
           !$acc parallel loop vector !async(stream2)
           do jj = 2, nj-1
              v(mi,jj) = v(mi,jj)+deltaxv(mi)*(corr_pres(mi,jj)-corr_pres(mi,jj+1))/av(mi,jj)
           end do
           !$acc wait
           !
           !--------------------------------------------------
           !--------------------------------------------------
           !
           ! Se resuelve la ecuaci\'on de balance de energ\'ia
           !
           !--------------------------------------------------
           !--------------------------------------------------
           solucion_energia: do
              !
              !$acc parallel loop gang collapse(2) !async(stream2)
              inicializacion_ftemp: do jj=2, nj
                 do ii = 2, mi
                    ftemp(ii,jj) = temp(ii,jj)
                 end do
              end do inicializacion_ftemp
              !
              !-------------------------
              !
              ! Condiciones de frontera
              !
              !$acc parallel loop vector !async(stream1)
              bucle_direccionxp: do ii = 2, mi
                 !***********************
                 !Condiciones de frontera
                 if (ii>=placa_min .and. ii<=placa_max) then
                    BC(1,ii) = 1.0_DBL
                    BN(1,ii) = 0.0_DBL
                    Ry(1,ii) = 1.0_DBL
                 else
                    BC(1,ii) =-1.0_DBL
                    BN(1,ii) = 1.0_DBL
                    Ry(1,ii) = 0.0_DBL
                 end if
                 !-----------------
                 if (ii>=placa_min .and. ii<=placa_max) then
                    BC(nj+1,ii) = 1.0_DBL
                    BN(nj+1,ii) = 0.0_DBL
                    Ry(nj+1,ii) = 1.0_DBL
                 else
                    BC(nj+1,ii)  = 1.0_DBL
                    BS(nj+1,ii)  =-1.0_DBL
                    Ry(nj+1,ii)  = 0.0_DBL
                 end if
                 ! ------------------
              end do bucle_direccionxp
              !
              !$acc parallel loop vector !async(stream1)
              do jj = 2, nj
                 !------------------------
                 ! Condiciones de frontera
                 AC(1,jj) =-1._DBL
                 AD(1,jj) = 1._DBL
                 Rx(1,jj) = 0._DBL
                 !
                 AI(mi+1,jj) =-1._DBL 
                 AC(mi+1,jj) = 1._DBL
                 Rx(mi+1,jj) = 0.0_DBL
                 !
              end do
              !------------------------------------------
              !
              ! Se ensambla la ecuaci\'on de la energ\'ia
              !
              !$acc parallel loop gang !async(stream1)
              do jj = 2, nj
                 !$acc loop vector
                 do ii = 2, mi
                    call ensambla_energia(deltaxp,deltayp,&
                         &deltaxu,deltayv,fexu,feyv,gamma_energ,&
                         &u,v,&
                         &temp,temp_ant,dt,&
                         &rel_ener,placa_min,placa_max,&
                         &AI,AC,AD,Rx,BS,BC,BN,Ry,&
                         &ii,jj)
                 end do
              end do
              !
              !---------------------------------------------
              !
              ! Soluci\'on de la ecuaci\'on de la energ\'ia
              !
              !$acc parallel loop gang async(stream2) wait(stream1)
              solucion_energia_x: do jj = 2, nj
                 
                 call tridiagonal(AI(1:mi+1,jj),AC(1:mi+1,jj),AD(1:mi+1,jj),Rx(1:mi+1,jj),mi+1)
                 temp(1,jj)    = Rx(1,jj)
                 temp(mi+1,jj) = Rx(mi+1,jj)

              end do solucion_energia_x
              !
              !$acc parallel loop gang async(stream1) 
              solucion_energia_y: do ii = 2, mi

                 call tridiagonal(BS(1:nj+1,ii),BC(1:nj+1,ii),BN(1:nj+1,ii),Ry(1:nj+1,ii),nj+1)
                 temp(ii,1)    = Ry(1,ii)
                 temp(ii,nj+1) = Ry(nj+1,ii)

              end do solucion_energia_y
              !$acc wait
              !
              !$acc parallel loop gang collapse(2) !async(stream1) wait(stream2)
              do jj = 2, nj
                 do ii = 2, mi
                    temp(ii,jj) = 0.5_DBL*Rx(ii,jj)+0.5_DBL*Ry(jj,ii)
                 end do
              end do
              !$acc wait
              !
              ! error de la ecuaci\'on de la energ\'ia
              !
              error = 0.0_DBL
              !
              !$acc parallel loop reduction(max:error) !async(stream1) wait(stream2)
              calculo_diferencias_dtemp: do jj = 2, nj
                 do ii = 2, mi
                    error = max(error, abs(temp(ii,jj)-ftemp(ii,jj)))
                 end do
              end do calculo_diferencias_dtemp
              !
              !------------------------------------------
              !
              ! Criterio de convergencia de la energ\'ia
              !
              !$acc wait
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
                 ! write(*,*) 'temp', maxval(abs(dtemp))
              end if
              !
           end do solucion_energia
           !$acc wait
           !
           !--------------------------------------------
           !--------------------------------------------
           !
           ! Criterio de convergencia del paso de tiempo
           !
           !--------------------------------------------
           !--------------------------------------------
           !
           !$acc parallel !async(stream1)
           call residuo_u(deltaxu,deltayu,deltaxp,&
                &deltayv,fexp,feyp,fexu,gamma_momen,&
                &u,u_ant,v,&
                &temp,pres,Ri,dt,rel_vel,&
                &Resu&
                &)
           !$acc end parallel
           !
           !--------------------------------
           !
           ! residuo del algoritmo
           residuo = 0.0_DBL
           !$acc parallel loop reduction(max:residuo) !async(stream1)
           calculo_maximo_residuou: do jj=2, nj
              do ii = 2, mi-1
                 residuo = max(residuo, abs(Resu(ii,jj)))
              end do
           end do calculo_maximo_residuou
           !
           !$acc wait
           if ( maxbo<conv_paso .and. residuo<conv_resi)then
              iter_simple = 0
              write(102,*) 'SIMPLE', maxbo, residuo
              exit
           else if ( iter_simple > iter_simple_max ) then
              iter_simple = 0
              ! write(*,*) "Advertencia SIMPLE: convergencia no alcanzada, ", residuo,maxbo   
              exit
           else
              iter_simple = iter_simple + 1
              write(102,*) 'SIMPLE', maxbo, residuo
              ! write(*,*) 'tiempo ',itera,res_fluido_u,MAXVAL(ABS(Resu)),MAXVAL(ABS(b_o)),&
                   ! &MAXVAL(ABS(pres))
           end if
           !
        end do ALGORITMO_SIMPLE  !final del algoritmo SIMPLE
        !
        !-------------------------------------
        !
        ! Escritura de mensajes y postprocesos
        !
        itera = itera + 1
        if( mod(itera,500)==0 )write(*,*) 'ITERACION: ',itera,residuo,maxbo

        if( mod(itera,1000)==0 )then
           !     CALL entropia_cvt(x,y,u,xu,v,yv,temp,entropia_calor,entropia_viscosa,entropia,&
           !&entropia_int,temp_int,a_ent,lambda_ent)
           !
           !$acc update self(temp(1:mi+1,1:nj+1)) ! !async(stream2)
           ! $acc parallel !async(stream2)
           call nusselt_promedio_y(&
                &xp,yp,deltaxp,deltayp,&
                &temp,nusselt0,nusselt1,&
                &placa_min,placa_max&
                &)
           ! $acc end parallel
           ! $acc wait
           !
           temp_med = (temp((placa_min+placa_max)/2,nj/2+1)+&
                &temp((placa_min+placa_max)/2+1,nj/2+1))/2._DBL
           OPEN(unit = 5,file='nuss_sim_n'//trim(njc)//'m'//trim(mic)//'_R'&
                &//trim(Rec)//'.dat',access = 'append')
           WRITE(5,form26) tiempo_inicial+itera*dt,nusselt0,&
                &-nusselt1,temp_med,temp_int,entropia_int
           CLOSE(unit = 5)
        end if
        !*********************************
        tiempo   = tiempo_inicial+itera*dt
        !
        !$acc parallel loop gang collapse(2) !async(stream1)
        do jj=1, nj+1
           do ii=1, mi+1
              temp_ant(ii,jj) = temp(ii,jj)
           end do
        end do
        !
        !$acc parallel loop gang collapse(2) !async(stream1)
        do jj=1, nj+1
           do ii=1, mi
              u_ant(ii,jj) = u(ii,jj)
           end do
        end do
        !
        !$acc parallel loop gang collapse(2) !async(stream2)
        do jj=1, nj
           do ii=1, mi+1
              v_ant(ii,jj) = v(ii,jj)
           end do
        end do
        !
     end do
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
     unidad      = (itera_total-millar*1000*paq_itera-centena*100*paq_itera-decena*10*paq_itera)&
          &/(paq_itera)
     decima      = (itera_total-millar*1000*paq_itera-centena*100*paq_itera-decena*10*paq_itera&
          &-unidad*paq_itera)/(paq_itera/10)
     WRITE(dec,16)decima;16 format(I1)
     WRITE(un,16) unidad
     WRITE(de,16) decena
     WRITE(ce,16) centena
     WRITE(m,16)  millar
     DO jj = 2, nj
        DO ii = 2, mi
           uf(ii,jj) = (u(ii,jj)+u(ii-1,jj))/2._DBL
        END DO
     END DO
     DO ii = 2, mi
        DO jj = 2, nj
           vf(ii,jj) = (v(ii,jj)+v(ii,jj-1))/2._DBL
        END DO
     END DO
     DO jj = 1, nj+1
        uf(1,jj)    = u(1,jj)
        uf(mi+1,jj) = u(mi,jj)
     END DO
     DO ii = 1, mi+1
        vf(ii,1)    = v(ii,1)
        vf(ii,nj+1) = v(ii,nj)
     END DO
     !************************************
     WRITE(*,*) 'itera_total=',itera_total
     WRITE(*,103) nusselt0,nusselt1
     WRITE(*,104) maxbo,residuo
     WRITE(*,105) MAXVAL(ABS(Restemp)),MAXVAL(ABS(Resv))
     WRITE(*,*)' '
103  FORMAT(1X,'N_Izq=',D23.15,', N_Der=',D23.15)
104  FORMAT(1X,'b_o  =',D23.15,', Res_u=',D23.15)
105  FORMAT(1X,'Res_T=',D23.15,', Res_v=',D23.15)
     !********************************
     !*** Formato de escritura dat ***
     ! OPEN(unit=2,file='out_n'//njc//'m'//mic//'_R'//Rec//'u.dat')
     ! WRITE(2,*) placa_min,placa_max,itera_total,ao
     ! DO j = 1, nj+1
     !    DO i = 1, mi
     !       WRITE(2,form24) xu(i),yp(j),u(i,j)
     !    END DO
     ! END DO
     ! CLOSE(unit=2)
     ! OPEN(unit=3,file='out_n'//njc//'m'//mic//'_R'//Rec//'v.dat')
     ! WRITE(3,*) placa_min,placa_max,itera_total,ao
     ! DO j = 1, nj
     !    DO i = 1, mi+1
     !       WRITE(3,form24) xp(i),yv(j),v(i,j)
     !    END DO
     ! END DO
     ! CLOSE(unit=3)
     ! OPEN(unit=4,file='out_n'//njc//'m'//mic//'_R'//Rec//'p.dat')
     ! WRITE(4,*) placa_min,placa_max,itera_total,ao
     ! DO j = 1, nj+1
     !    DO i = 1, mi+1
     !       WRITE(4,form25) xp(i),yp(j),temp(i,j),pres(i,j)
     !    END DO
     ! END DO
     ! CLOSE(unit=4)
     ! *************************************
     ! *** Formato escritura VTK ***********
     ! sample  = m//ce//de//un//dec
     archivo = 'n'//trim(njc)//'m'//trim(mic)//'R'//trim(Rec)//'/t_'//m//ce//de//un//dec//'.vtk'
     call postproceso_bin(xu,yv,xp,yp,u,v,pres,temp,b_o,Ra)
     CALL postproceso_vtk(xp,yp,uf,vf,pres,temp,b_o,archivo)
     !
     !
  end do !*********** final del repetidor principal
  !
end program SIMPLE2D
