!
!
! M\'odulo de herramientas de malla para el c\'odigo SIMPLE2D
!
! Este m\'odulo contiene herramientas que se utilizan para leer las mallas, hacer interpolaciones,
! definir arreglos de incrementos del dominio computacional para mallas escalonadas
!
!
MODULE malla
  !
  implicit none
  !
  ! Dimensiones de la malla
  !
  INTEGER(4), PARAMETER :: mi=1001, nj=251, nsolid = 7
  !
  ! Tipo de variables reales, definen precici'on y rango
  !
  INTEGER(4), PARAMETER :: SGL=SELECTED_REAL_KIND(P=6,R=37)
  INTEGER(4), PARAMETER :: DBL=SELECTED_REAL_KIND(P=15,R=307)
  !
  ! formatos de escritura
  !
  character(len=9), parameter :: form24="(3D23.15)",form25="(4D23.15)"
  character(len=9), parameter :: form26="(6D23.15)",form27="(10e15.7)"
  !
  ! viscosidad para frontera inmersa
  !
  REAL(kind=DBL),   PARAMETER :: cero=0._DBL, visc_solido = 1.e40_DBL
  !
  ! caracteres para los nombres de archivos de salida
  !
  character(20) :: mic, njc, zkc
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
  !
contains
  !
  !-----------------------------------------------------------------------------
  !
  ! Esta subrutina abre los archivos para leer las mallas escalonadas
  ! para u, v, p y t. Lee los arreglos que definen las coordenadas de los nodos,
  ! define los incrementos que dependen \'unicamente de cantidades
  ! geom\'etricas y los devuelve como variables de salida. Tambi\'en lee los
  ! valores iniciales de las variables f\'isicas u,v,p,T.
  !
  subroutine lectura_mallas_escalonadas(entrada_uo,entrada_vo,entrada_tpo,&
       &u_anto,v_anto,presso,temp_anto,&
       &xpo,ypo,xuo,yvo,&
       &deltaxpo,deltaypo,&
       &deltaxuo,deltayuo,&
       &deltaxvo,deltayvo,&
       &fexpo,feypo,fexuo,feyvo,&
       &ao,i_0,i_1,iter_inio)
    ! ------------------------------------------------------------
    !
    ! Variables para los archivos de la entrada de datos
    !
    CHARACTER(len=28), intent(in):: entrada_uo,entrada_vo,entrada_tpo
    !
    ! Variables para las velocidades, temperatura y presion iniciales
    !
    real(kind=DBL), dimension(mi,nj+1),   intent(out) :: u_anto
    real(kind=DBL), dimension(mi+1,nj),   intent(out) :: v_anto
    real(kind=DBL), dimension(mi+1,nj+1), intent(out) :: presso
    real(kind=DBL), dimension(mi+1,nj+1), intent(out) :: temp_anto
    !
    ! Variables de la malla, volumen de control e interpolaciones
    !
    real(kind=DBL), dimension(mi+1),      intent(out) :: xpo
    real(kind=DBL), dimension(nj+1),      intent(out) :: ypo
    REAL(kind=DBL), DIMENSION(mi),        intent(out) :: xuo
    REAL(kind=DBL), DIMENSION(nj),        intent(out) :: yvo
    real(kind=DBL), dimension(mi),        intent(out) :: deltaxpo
    real(kind=DBL), dimension(nj),        intent(out) :: deltaypo
    real(kind=DBL), dimension(mi),        intent(out) :: deltaxuo
    real(kind=DBL), dimension(nj),        intent(out) :: deltayuo
    real(kind=DBL), dimension(mi),        intent(out) :: deltaxvo   
    real(kind=DBL), dimension(nj),        intent(out) :: deltayvo
    REAL(kind=DBL), DIMENSION(mi),        intent(out) :: fexpo
    REAL(kind=DBL), DIMENSION(nj),        intent(out) :: feypo
    REAL(kind=DBL), DIMENSION(mi-1),      intent(out) :: fexuo
    REAL(kind=DBL), DIMENSION(nj-1),      intent(out) :: feyvo
    !
    ! Dimensiones del dominio, ubicaci\'on de las placas calientes, iteracion inicial
    !
    real(kind=DBL) :: ao
    integer        :: i_0, i_1
    integer        :: iter_inio
    !
    ! Variables auxiliares
    !
    integer :: ii, jj
    !***********************************************
    !
    ! Inicializaci\'on de los arreglos geom\'etricos
    !
    xpo      =  0.0_DBL
    ypo      =  0.0_DBL
    xuo      =  0.0_DBL
    yvo      =  0.0_DBL
    deltaxpo =  0.0_DBL
    deltaypo =  0.0_DBL
    deltaxuo =  0.0_DBL
    deltayuo =  0.0_DBL
    deltaxvo =  0.0_DBL
    deltayvo =  0.0_DBL
    fexpo    =  0.0_DBL
    feypo    =  0.0_DBL
    fexuo    =  0.0_DBL
    feyvo    =  0.0_DBL
    !
    ! Inicializaci\'on de los arreglos f\'isicos
    !
    u_anto    =  0.0_DBL
    v_anto    =  0.0_DBL
    temp_anto =  0.0_DBL
    presso    =  0.0_DBL
    !
    ! Apertura de los archivos y lectura de los arreglos
    !
    OPEN(unit=11,file=entrada_uo)
    READ(11,*)i_0,i_1,iter_inio,ao
    DO jj = 1, nj+1
       DO ii = 1, mi
          READ(11,form24) xuo(ii),ypo(jj),u_anto(ii,jj)
       END DO
    END DO
    CLOSE(unit=11)
    !***************************
    Open(unit=12,file=entrada_vo)
    READ(12,*)i_0,i_1,iter_inio,ao
    DO jj = 1, nj
       DO ii = 1, mi+1
          READ(12,form24) xpo(ii),yvo(jj),v_anto(ii,jj)
       END DO
    END DO
    CLOSE(unit=12)
    !****************************
    OPEN(unit=13,file=entrada_tpo)
    READ(13,*)i_0,i_1,iter_inio,ao
    DO jj = 1, nj+1
       DO ii = 1, mi+1
          READ(13,form25) xpo(ii),ypo(jj),temp_anto(ii,jj),presso(ii,jj)
       END DO
    END DO
    CLOSE(unit=13)
    !
    ! C\'alculo de los tamanios de los vol\'umenes de control
    ! malla de la presi\'on y la temperatura
    !
    do ii = 2, mi
       deltaxpo(ii) = xuo(ii)-xuo(ii-1)
    end do
    do jj = 2, nj
       deltaypo(jj) = yvo(jj)-yvo(jj-1)
    end do
    !
    ! C\'alculo de los tamanios de los vol\'umenes de control
    ! malla para la velocidad u(mi,nj+1)
    !
    do ii = 1, mi
       deltaxuo(ii) = xpo(ii+1)-xpo(ii)
    end do
    do jj = 2, nj
       deltayuo(jj) = yvo(jj)-yvo(jj-1)
    end do
    !
    ! C\'alculo de los tamanios de los vol\'umenes de control
    ! malla para la velocidad v(mi+1,nj)
    !
    do ii = 2, mi
       deltaxvo(ii) = xuo(ii)-xuo(ii-1)
    end do
    do jj = 1, nj
       deltayvo(jj) = ypo(jj+1)-ypo(jj)
    end do
    !
    ! C\'alculo de los arreglos de interpolaci\'on para la malla
    ! de la presi\'on y la temperatura
    !
    do ii = 1, mi
       fexpo(ii) = ( xpo(ii+1)-xuo(ii) ) / ( xpo(ii+1)-xpo(ii)  )
    end do
    do jj = 1, nj
       feypo(jj) = ( ypo(jj+1)-yvo(jj) ) / ( ypo(jj+1)-ypo(jj)  )
    end do
    !
    ! C\'alculo de los arreglos de interpolaci\'on para la malla
    ! de la velocidad u(mi,nj+1). Solo se necesita un arreglo
    ! en direcci\'on ii
    !
    do ii = 1, mi-1
       fexuo(ii) = ( xuo(ii+1)-xpo(ii+1) ) / ( xuo(ii+1)-xuo(ii)  )
    end do
    !
    ! C\'alculo de los arreglos de interpolaci\'on para la malla
    ! de la velocidad v(mi+1,nj). Solo se necesita un arreglo
    ! en direcci\'on jj
    !
    do jj = 1, nj-1
       feyvo(jj) = ( yvo(jj+1)-ypo(jj+1) ) / ( yvo(jj+1)-yvo(jj)  )
    end do
    
  end subroutine lectura_mallas_escalonadas

end MODULE malla
