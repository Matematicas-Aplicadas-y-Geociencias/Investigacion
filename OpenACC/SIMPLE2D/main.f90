PROGRAM SIMPLE
USE constantes
! USE mkl95_LAPACK
IMPLICIT NONE
INCLUDE 'omp_lib.h'
INTEGER :: i,j,k,l,itera_total,itera,itera_inicial,i_o,i_1,paq_itera
INTEGER :: millar,centena,decena,unidad,decima,id,nthreads
!*******************************************
! Variables del flujo,entropia,nusselt e inc'ognitas,residuos,relajaci'on y convergencia
REAL(kind=DBL), DIMENSION(mi,nj+1)   :: u,u_ant,du,au,Resu,gamma_u,Ri
REAL(kind=DBL), DIMENSION(mi+1,nj)   :: v,v_ant,dv,av,Resv,gamma_v
REAL(kind=DBL), DIMENSION(mi+1,nj+1) :: temp,temp_ant,dtemp,Restemp,pres,corr_pres,dcorr_pres
REAL(kind=DBL), DIMENSION(mi+1,nj+1) :: entropia_calor,entropia_viscosa,entropia,uf,vf,b_o,gamma_t
REAL(kind=DBL) :: temp_med,nusselt0,nusselt1,entropia_int,temp_int,gamma_s
REAL(kind=DBL) :: conv_u,conv_p,conv_t,conv_resi,conv_paso,rel_pres,rel_vel,rel_tem
!********************************************
!Variables de la malla, volumen de control, incremento de tiempo y nums Reynolds, Peclet
!Rayleigh, Richardson, Prandtl, valores de las constantes de difusividad
REAL(kind=DBL), DIMENSION(mi)   :: xu
REAL(kind=DBL), DIMENSION(nj)   :: yv
REAL(kind=DBL), DIMENSION(mi+1) :: x,fexu
REAL(kind=DBL), DIMENSION(nj+1) :: y,feyv
REAL(kind=DBL), DIMENSION(mi-1) :: d_xu,d2_xu
REAL(kind=DBL), DIMENSION(nj-1) :: d_yv,d2_yv
REAL(kind=DBL)   :: ao
REAL(kind=DBL)   :: tiempo,tiempo_inicial,dt,Ra,Pr,Ri_1
REAL(kind=DBL)   :: a_ent,lambda_ent
CHARACTER(len=1) :: dec,un,de,ce,m
CHARACTER(len=3) :: njc,mic,Rec
CHARACTER(len=5) :: sample
CHARACTER(len=22):: entrada_u,entrada_v,entrada_tp
CHARACTER(len=46):: archivo=repeat(' ',46)
LOGICAL          :: res_fluido_u
!****************************************
!Variables de caracterizaci'on del fluido
REAL(kind=DBL) :: temp_ref,visc_cin,dif_term,cond_ter,cons_gra,coef_exp,long_ref,dens_ref
!****************************
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
  READ (10,*) rel_tem     ! relajaci'on de la temperatura
  READ (10,*) conv_u      ! convergencia de la velocidad
  READ (10,*) conv_t      ! convergencia de la temperatura
  READ (10,*) conv_p      ! convergencia de la presi'on
  READ (10,*) conv_resi   ! convergencia del residuo
  READ (10,*) conv_paso   ! convergencia del paso de tiempo
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
!gamma_s = 10._DBL*(1._DBL/(Re*Pr))
gamma_t = 1._DBL/(Ra*Pr) !sqrt(1._DBL/(Pr*Ra))
gamma_u = 1._DBL/(Ra)    !sqrt(Pr/Ra)
gamma_v = 1._DBL/(Ra)    !sqrt(Pr/Ra)
Ri      = Ri_1
!***************************
OPEN(unit=11,file=entrada_u)
READ(11,*)i_o,i_1,itera_inicial,ao
DO j = 1, nj+1
  DO i = 1, mi
    READ(11,24) xu(i),y(j),u_ant(i,j)
  END DO
END DO
CLOSE(unit=11)
!***************************
Open(unit=12,file=entrada_v)
READ(12,*)i_o,i_1,itera_inicial,ao
DO j = 1, nj
  DO i = 1, mi+1
    READ(12,24) x(i),yv(j),v_ant(i,j)
  END DO
END DO
CLOSE(unit=12)
!****************************
OPEN(unit=13,file=entrada_tp)
READ(13,*)i_o,i_1,itera_inicial,ao
DO j = 1, nj+1
  DO i = 1, mi+1
    READ(13,25) x(i),y(j),temp_ant(i,j),pres(i,j)
  END DO
END DO
CLOSE(unit=13)
!***********
!incrementos
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
fexu(1)    = 0_DBL
DO i = 2, mi+1
  fexu(i)  = (x(i) - xu(i-1)) / (x(i) - x(i-1))
END DO
feyv(1)    = 0_DBL
DO j = 2, nj+1
  feyv(j)  = (y(j) - yv(j-1)) / (y(j) - y(j-1))
END DO
!*****************
!valores iniciales
tiempo_inicial = itera_inicial*dt
! u_ant = cero
! v_ant = cero
! temp_ant = cero
u     = u_ant
v     = v_ant
uf    = cero
vf    = cero
temp  = temp_ant
corr_pres = cero * dfloat(mi+1-i)/dfloat(mi)
au    = 1._DBL
av    = 1._DBL
b_o   = cero
itera = 0
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
!*********************************************************
DO l=1,itermax/paq_itera   !inicio del repetidor principal
DO k=1,paq_itera           !inicio del paquete iteraciones
ALGORITMO_SIMPLE: DO       !inicio del algoritmo SIMPLE
    DO
      CALL vel_u(y,feyv,d_xu,d2_xu,d_yv,u,u_ant,v,temp,pres,gamma_u,Ri,dt,du,au,rel_vel)
      CALL vel_v(x,fexu,d_yv,d2_yv,d_xu,u,v,v_ant,pres,gamma_v,dt,dv,av,rel_vel)
      !****************************************
      !Criterio de convergencia de la velocidad
      IF(MAXVAL(DABS(du))<conv_u.and.MAXVAL(DABS(dv))<conv_u)EXIT
!       WRITE(*,*) 'velocidad ',MAXVAL(DABS(du)),MAXVAL(DABS(dv))
    END DO
    !****************************************
    !se calcula la correcci'on de la presi'on
    corr_pres = cero
    DO
      CALL corrector_presion(corr_pres,d_xu,d_yv,u,v,b_o,au,av,dcorr_pres)
      !****************************************************
      !critero de convergencia del corrector de la presi'on
      IF(MAXVAL(DABS(dcorr_pres))<conv_p)EXIT
!       WRITE(*,*) 'corrector presion ', MAXVAL(DABS(dcorr_pres)), MAXVAL(DABS(b_o))!, MAXVAL(DABS(corr_pres))
    END DO
    corr_pres = rel_pres * corr_pres
    !*********************
    !se corrige la presion
    !$OMP PARALLEL DO
    DO i = 2, mi
      DO j = 2, nj
        pres(i,j) = pres(i,j) + corr_pres(i,j)
      END DO
    END DO
    !$OMP END PARALLEL DO
    !*****************************
    !se actualizan las velocidades
    !$OMP PARALLEL DO
    DO i = 2, mi-1
      DO j = 2, nj-1
        u(i,j) = u(i,j)+d_yv(j-1)*(corr_pres(i,j)-corr_pres(i+1,j))/au(i,j)
        v(i,j) = v(i,j)+d_xu(i-1)*(corr_pres(i,j)-corr_pres(i,j+1))/av(i,j)
      END DO
    END DO
    !$OMP END PARALLEL DO
    DO i = 2, mi-1
      u(i,nj) = u(i,nj)+d_yv(nj-1)*(corr_pres(i,nj)-corr_pres(i+1,nj))/au(i,nj)
    END DO
    DO j = 2, nj-1
      v(mi,j) = v(mi,j)+d_xu(mi-1)*(corr_pres(mi,j)-corr_pres(mi,j+1))/av(mi,j)
    END DO
!     corr_pres = corr_pres/rel_pres
    !*************************
    !se calcula la temperatura
    DO
      CALL temperatura(x,y,fexu,feyv,d_xu,d_yv,u,v,temp,temp_ant,gamma_t,dt,dtemp,i_o,i_1,rel_tem)
      !************************************
      !Criterio de convergencia temperatura
      IF(MAXVAL(DABS(dtemp))<conv_t)EXIT
!       WRITE(*,*) 'temp', MAXVAL(DABS(dtemp))
    END DO
    !*******************************************
    !Criterio de convergencia del paso de tiempo
    CALL residuou(res_fluido_u,xu,y,feyv,d_xu,d2_xu,d_yv,u,u_ant,v,temp,pres,Resu,gamma_u,Ri,dt)
    IF(res_fluido_u .EQV. .TRUE. .AND. MAXVAL(DABS(Resu))<conv_resi .AND. MAXVAL(DABS(b_o))<conv_paso)EXIT
!     IF( MAXVAL(DABS(Resu))<conv_resi .AND. MAXVAL(DABS(b_o))<conv_paso)EXIT
    ! WRITE(*,*) 'tiempo ',itera,res_fluido_u,MAXVAL(DABS(Resu)),MAXVAL(DABS(b_o)),MAXVAL(DABS(pres))
  END DO ALGORITMO_SIMPLE  !final del algoritmo SIMPLE
  itera = itera + 1
  IF(mod(itera,100)==0) WRITE(*,*) 'tiempa ',itera,res_fluido_u,MAXVAL(DABS(Resu)),MAXVAL(DABS(b_o)),MAXVAL(DABS(pres))

  IF(mod(itera,100)==0)THEN
!     CALL entropia_cvt(x,y,u,xu,v,yv,temp,entropia_calor,entropia_viscosa,entropia,entropia_int,temp_int,a_ent,lambda_ent)
    CALL nusselt(x,y,d_xu,d_yv,temp,nusselt0,nusselt1,i_o,i_1)
    temp_med = (temp((i_o+i_1)/2,nj/2+1)+temp((i_o+i_1)/2+1,nj/2+1))/2._DBL
    OPEN(unit = 5,file = 'nuss_sim_n'//njc//'m'//mic//'_R'//Rec//'.dat',access = 'append')
    WRITE(5,26) tiempo_inicial+itera*dt,nusselt0,-nusselt1,temp_med,temp_int,entropia_int
    CLOSE(unit = 5)
  ENDIF
  !*********************************
  tiempo   = tiempo_inicial+itera*dt
  temp_ant = temp
  u_ant    = u
  v_ant    = v
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
103 FORMAT(1X,'N_Izq=',D23.15,', N_Der=',D23.15)
104 FORMAT(1X,'b_o  =',D23.15,', Res_u=',D23.15)
105 FORMAT(1X,'Res_T=',D23.15,', Res_v=',D23.15)
!********************************
!*** Formato de escritura dat ***
OPEN(unit=2,file='out_n'//njc//'m'//mic//'_R'//Rec//'u.dat')
WRITE(2,*) i_o,i_1,itera_total,ao
DO j = 1, nj+1
  DO i = 1, mi
    WRITE(2,24) xu(i),y(j),u(i,j)
  END DO
END DO
CLOSE(unit=2)
OPEN(unit=3,file='out_n'//njc//'m'//mic//'_R'//Rec//'v.dat')
WRITE(3,*) i_o,i_1,itera_total,ao
DO j = 1, nj
  DO i = 1, mi+1
    WRITE(3,24) x(i),yv(j),v(i,j)
  END DO
END DO
CLOSE(unit=3)
OPEN(unit=4,file='out_n'//njc//'m'//mic//'_R'//Rec//'p.dat')
WRITE(4,*) i_o,i_1,itera_total,ao
DO j = 1, nj+1
  DO i = 1, mi+1
    WRITE(4,25) x(i),y(j),temp(i,j),pres(i,j)
  END DO
END DO
CLOSE(unit=4)
! *************************************
! *** Formato escritura VTK ***********
! sample  = m//ce//de//un//dec
archivo = 'n'//njc//'m'//mic//'R'//Rec//'/t_'//m//ce//de//un//dec//'.vtk'
CALL postprocess_vtk(x,y,uf,vf,pres,temp,b_o,archivo)
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
24 FORMAT (3D23.15)
25 FORMAT (4D23.15)
26 FORMAT (6D23.15)
27 FORMAT (10e15.7)
END DO!***********final del repetidor principal
!**********************************************
END PROGRAM SIMPLE
