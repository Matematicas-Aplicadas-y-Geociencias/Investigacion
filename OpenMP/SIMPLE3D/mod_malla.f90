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
  INTEGER(4), PARAMETER :: mi=100, nj=100, lk=100, nsolid = 7
  !
  ! Tipo de variables reales, definen precici'on y rango
  !
  INTEGER(4), PARAMETER :: SGL=SELECTED_REAL_KIND(P=6,R=37)
  INTEGER(4), PARAMETER :: DBL=SELECTED_REAL_KIND(P=15,R=307)
  !
  ! formatos de escritura
  !
  character(len=9), parameter :: form21="(1D23.15)",form22="(2D23.15)"
  character(len=9), parameter :: form23="(3D23.15)",form24="(4D23.15)"
  character(len=9), parameter :: form25="(5D23.15)",form26="(6D23.15)"
  character(len=9), parameter :: form27="(10e15.7)"
  !
  ! viscosidad para frontera inmersa
  !
  REAL(kind=DBL),   PARAMETER :: cero=0._DBL, visc_solido = 1.e40_DBL
  !
  ! caracteres para los nombres de archivos de salida
  !
  character(3) :: mic, njc, lkc
  !
  ! Variables de la malla, volumen de control y factores de interpolaci\'on
  !
  REAL(kind=DBL), DIMENSION(mi)   ::  xu
  REAL(kind=DBL), DIMENSION(nj)   ::  yv
  real(kind=DBL), dimension(lk)   ::  zw
  REAL(kind=DBL), DIMENSION(mi+1) ::  xp
  REAL(kind=DBL), DIMENSION(nj+1) ::  yp
  real(kind=DBL), dimension(lk+1) ::  zp
  real(kind=DBL), dimension(mi)   ::  deltaxp
  real(kind=DBL), dimension(nj)   ::  deltayp
  real(kind=DBL), dimension(lk)   ::  deltazp
  real(kind=DBL), dimension(mi)   ::  deltaxu
  real(kind=DBL), dimension(nj)   ::  deltayu
  real(kind=DBL), dimension(lk)   ::  deltazu
  real(kind=DBL), dimension(mi)   ::  deltaxv   
  real(kind=DBL), dimension(nj)   ::  deltayv
  real(kind=DBL), dimension(lk)   ::  deltazv
  real(kind=DBL), dimension(mi)   ::  deltaxw
  real(kind=DBL), dimension(nj)   ::  deltayw
  real(kind=DBL), dimension(lk)   ::  deltazw
  REAL(kind=DBL), DIMENSION(mi)   ::  fexp
  REAL(kind=DBL), DIMENSION(nj)   ::  feyp
  real(kind=DBL), dimension(lk)   ::  fezp
  REAL(kind=DBL), DIMENSION(mi-1) ::  fexu
  REAL(kind=DBL), DIMENSION(nj-1) ::  feyv
  real(kind=DBL), dimension(lk-1) ::  fezw
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
  subroutine lectura_mallas_escalonadas(entrada_uo,entrada_vo,entrada_wo,&
       &entrada_tpo,entrada_xyzo,&
       &u_anto,v_anto,w_anto,presso,temp_anto,&
       &xpo,ypo,zpo,xuo,yvo,zwo,&
       &deltaxpo,deltaypo,deltazpo,&
       &deltaxuo,deltayuo,deltazuo,&
       &deltaxvo,deltayvo,deltazvo,&
       &deltaxwo,deltaywo,deltazwo,&
       &fexpo,feypo,fezpo,fexuo,feyvo,fezwo,&
       &ao,iter_inio)
    ! ------------------------------------------------------------
    !
    ! Variables para los archivos de la entrada de datos
    !
    character(len=64), intent(in) :: entrada_uo,entrada_vo
    character(len=64), intent(in) :: entrada_wo,entrada_tpo
    character(len=64), intent(in) :: entrada_xyzo
    !
    ! Variables para las velocidades, temperatura y presion iniciales
    !
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(out) :: u_anto
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(out) :: v_anto
    real(kind=DBL), dimension(mi+1,nj+1,lk),   intent(out) :: w_anto
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(out) :: presso
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(out) :: temp_anto
    !
    ! Variables de la malla, volumen de control e interpolaciones
    !
    real(kind=DBL), dimension(mi+1),      intent(out) :: xpo
    real(kind=DBL), dimension(nj+1),      intent(out) :: ypo
    real(kind=DBL), dimension(lk+1),      intent(out) :: zpo
    REAL(kind=DBL), DIMENSION(mi),        intent(out) :: xuo
    REAL(kind=DBL), DIMENSION(nj),        intent(out) :: yvo
    real(kind=DBL), dimension(lk),        intent(out) :: zwo
    real(kind=DBL), dimension(mi),        intent(out) :: deltaxpo
    real(kind=DBL), dimension(nj),        intent(out) :: deltaypo
    real(kind=DBL), dimension(lk),        intent(out) :: deltazpo
    real(kind=DBL), dimension(mi),        intent(out) :: deltaxuo
    real(kind=DBL), dimension(nj),        intent(out) :: deltayuo
    real(kind=DBL), dimension(lk),        intent(out) :: deltazuo
    real(kind=DBL), dimension(mi),        intent(out) :: deltaxvo   
    real(kind=DBL), dimension(nj),        intent(out) :: deltayvo
    real(kind=DBL), dimension(lk),        intent(out) :: deltazvo
    real(kind=DBL), dimension(mi),        intent(out) :: deltaxwo   
    real(kind=DBL), dimension(nj),        intent(out) :: deltaywo
    real(kind=DBL), dimension(lk),        intent(out) :: deltazwo
    REAL(kind=DBL), DIMENSION(mi),        intent(out) :: fexpo
    REAL(kind=DBL), DIMENSION(nj),        intent(out) :: feypo
    real(kind=DBL), dimension(lk),        intent(out) :: fezpo
    REAL(kind=DBL), DIMENSION(mi-1),      intent(out) :: fexuo
    REAL(kind=DBL), DIMENSION(nj-1),      intent(out) :: feyvo
    real(kind=DBL), dimension(lk-1),      intent(out) :: fezwo
    !
    ! Dimensiones del dominio, ubicaci\'on de las placas calientes, iteracion inicial
    !
    real(kind=DBL) :: ao
    integer        :: i_0, i_1
    integer        :: iter_inio
    !
    ! Variables auxiliares
    !
    integer :: ii, jj, kk
    !***********************************************
    !
    ! Inicializaci\'on de los arreglos geom\'etricos
    !
    xpo      =  0.0_DBL
    ypo      =  0.0_DBL
    zpo      =  0.0_DBL
    xuo      =  0.0_DBL
    yvo      =  0.0_DBL
    zwo      =  0.0_DBL
    deltaxpo =  0.0_DBL
    deltaypo =  0.0_DBL
    deltazpo =  0.0_DBL
    deltaxuo =  0.0_DBL
    deltayuo =  0.0_DBL
    deltazuo =  0.0_DBL
    deltaxvo =  0.0_DBL
    deltayvo =  0.0_DBL
    deltazvo =  0.0_DBL
    deltaxwo =  0.0_DBL
    deltaywo =  0.0_DBL
    deltazwo =  0.0_DBL    
    fexpo    =  0.0_DBL
    feypo    =  0.0_DBL
    fezpo    =  0.0_DBL
    fexuo    =  0.0_DBL
    feyvo    =  0.0_DBL
    fezwo    =  0.0_DBL
    !
    ! Inicializaci\'on de los arreglos f\'isicos
    !
    u_anto    =  0.0_DBL
    v_anto    =  0.0_DBL
    w_anto    =  0.0_DBL
    temp_anto =  0.0_DBL
    presso    =  0.0_DBL
    !
    ! Apertura de los archivos y lectura de los arreglos
    !
    open(unit=11,file=entrada_uo)
    read(11,*)i_0,i_1,iter_inio,ao
    do kk = 1, lk+1
       do jj = 1, nj+1
          do ii = 1, mi
             read(11,form24) xuo(ii), ypo(jj), zpo(kk), u_anto(ii,jj,kk)
          end do
       end do
    end do
    close(unit=11)
    !***************************
    Open(unit=12,file=entrada_vo)
    READ(12,*)i_0,i_1,iter_inio,ao
    do kk = 1, lk+1
       do jj = 1, nj
          do ii = 1, mi+1
             read(12,form24) xpo(ii), yvo(jj), zpo(kk), v_anto(ii,jj,kk)
          end do
       end do
    end do
    close(unit=12)
    !****************************
    open(unit=13,file=entrada_wo)
    READ(13,*)i_0,i_1,iter_inio,ao
    do kk = 1, lk
       do jj = 1, nj+1
          do ii = 1, mi+1
             read(13,form24) xpo(ii), ypo(jj), zwo(kk),w_anto(ii,jj,kk)
          end do
       end do
    end do
    close(unit=13)
    !****************************
    open(unit=14,file=entrada_tpo)
    READ(14,*)i_0,i_1,iter_inio,ao
    do kk = 1, lk+1
       do jj = 1, nj+1
          do ii = 1, mi+1
             read(14,form25) xpo(ii),ypo(jj),zpo(kk),temp_anto(ii,jj,kk),presso(ii,jj,kk)
          end do
       end do
    end do
    close(unit=14)
    ! !**********************************
    ! !
    ! ! Lectura de los arreglos de malla
    ! !
    ! open(unit=15,file=entrada_xyzo)
    ! do ii = 1, mi
    !    read(15,form21) xuo
    ! end do
    ! do ii = 1, mi+1
    !    read(15,form21) xpo
    ! end do
    ! do jj = 1, nj
    !    read(15,form21) yvo
    ! end do
    ! do jj = 1, nj+1
    !    read(15,form21) ypo
    ! end do
    ! do kk = 1, lk
    !    read(15,form21) zwo
    ! end do
    ! do ii = 1, lk+1
    !    read(15,form21) zpo
    ! end do
    ! close(unit=15)
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
    do kk = 2, lk
       deltazpo(kk) = zwo(kk)-zwo(kk-1)
    end do
    !
    ! C\'alculo de los tamanios de los vol\'umenes de control
    ! malla para la velocidad u(mi,nj+1,lk+1)
    !
    do ii = 1, mi
       deltaxuo(ii) = xpo(ii+1)-xpo(ii)
    end do
    do jj = 2, nj
       deltayuo(jj) = yvo(jj)-yvo(jj-1)
    end do
    do kk = 2, lk
       deltazuo(kk) = zwo(kk)-zwo(kk-1)
    end do
    !
    ! C\'alculo de los tamanios de los vol\'umenes de control
    ! malla para la velocidad v(mi+1,nj,lk+1)
    !
    do ii = 2, mi
       deltaxvo(ii) = xuo(ii)-xuo(ii-1)
    end do
    do jj = 1, nj
       deltayvo(jj) = ypo(jj+1)-ypo(jj)
    end do
    do kk = 2, lk
       deltazvo(kk) = zwo(kk)-zwo(kk-1)
    end do
    !
    ! C\'alculo de los tamanios de los vol\'umenes de control
    ! malla para la velocidad w(mi+1,nj+1,lk)
    !
    do ii = 2, mi
       deltaxwo(ii) = xuo(ii)-xuo(ii-1)
    end do
    do jj = 2, nj
       deltaywo(jj) = yvo(jj)-yvo(jj-1)
    end do
    do kk = 1, lk
       deltazwo(kk) = zpo(kk+1)-zpo(kk)
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
    do kk = 1, lk
       fezpo(kk) = ( zpo(kk+1)-zwo(kk) ) / ( zpo(kk+1)-zpo(kk)  )
    end do
    !
    ! C\'alculo de los arreglos de interpolaci\'on para la malla
    ! de la velocidad u(mi,nj+1,lk+1). Solo se necesita un arreglo
    ! en direcci\'on ii
    !
    do ii = 1, mi-1
       fexuo(ii) = ( xuo(ii+1)-xpo(ii+1) ) / ( xuo(ii+1)-xuo(ii)  )
    end do
    !
    ! C\'alculo de los arreglos de interpolaci\'on para la malla
    ! de la velocidad v(mi+1,nj,lk+1). Solo se necesita un arreglo
    ! en direcci\'on jj
    !
    do jj = 1, nj-1
       feyvo(jj) = ( yvo(jj+1)-ypo(jj+1) ) / ( yvo(jj+1)-yvo(jj)  )
    end do
    !
    ! C\'alculo de los arreglos de interpolaci\'on para la malla
    ! de la velocidad w(mi+1,nj,lk). Solo se necesita un arreglo
    ! en direcci\'on kk
    !
    do kk = 1, lk-1
       fezwo(kk) = ( zwo(kk+1)-zpo(kk+1) ) / ( zwo(kk+1)-zwo(kk)  )
    end do
    
  end subroutine lectura_mallas_escalonadas
  !
  !-----------------------------------------------------------------------------
  !
  ! Estas funciones devuelve \'indices 1D para arreglos que representan
  ! cantidades tridimensionales con el objetivo de garantizar que los datos
  ! se encuentran en elementos consecutivos de memoria:
  !
  !     index(ii,jj,kk) = (dimx*dimy)*(kk-1)+dimx*(jj-1)+ii
  !     indey(jj,ii,kk) = (dimx*dimy)*(kk-1)+dimy*(ii-1)+jj
  !     indez(kk,jj,ii) = (dimy*dimz)*(ii-1)+dimz*(jj-1)+kk
  !
  function indexu(ii,jj,kk)
    !
    implicit none
    !
    integer             :: indexu
    integer, intent(in) :: ii,jj,kk
    !
    indexu = (mi*(nj+1))*(kk-1)+mi*(jj-1)+ii
    !
  end function indexu
  !
  function indeyu(jj,ii,kk)
    !
    implicit none
    !
    integer             :: indeyu
    integer, intent(in) :: ii,jj,kk
    !
    indeyu = (mi*(nj+1))*(kk-1)+(nj+1)*(ii-1)+jj
    !
  end function indeyu
  !
  function indexp(ii,jj,kk)
    !
    implicit none
    !
    integer             :: indexp
    integer, intent(in) :: ii,jj,kk
    !
    indexp = ((mi+1)*(nj+1))*(kk-1)+(mi+1)*(jj-1)+ii
    !
  end function indexp
  !
  function indexv(ii,jj,kk)
    !
    implicit none
    !
    integer             :: indexv
    integer, intent(in) :: ii,jj,kk
    !
    indexv = ((mi+1)*(nj))*(kk-1)+(mi+1)*(jj-1)+ii
    !
  end function indexv
  !
  function indeyv(jj,ii,kk)
    !
    implicit none
    !
    integer             :: indeyv
    integer, intent(in) :: ii,jj,kk
    !
    indeyv = (mi+1)*nj*(kk-1)+(nj)*(ii-1)+jj
    !
  end function indeyv
  !
  function indezv(kk,jj,ii)
    !
    implicit none
    !
    integer             :: indezv
    integer, intent(in) :: ii,jj,kk
    !
    indezv = (nj*(lk+1))*(ii-1)+lk*(jj-1)+kk
    !
  end function indezv
  !
  function indeyp(jj,ii,kk)
    !
    implicit none
    !
    integer             :: indeyp
    integer, intent(in) :: ii,jj,kk
    !
    indeyp = ((mi+1)*(nj+1))*(kk-1)+(nj+1)*(ii-1)+jj
    !
  end function indeyp
   !
  function indezw(kk,jj,ii)
    !
    implicit none
    !
    integer             :: indezw
    integer, intent(in) :: ii,jj,kk
    !
    indezw = (lk*(nj+1))*(ii-1)+lk*(jj-1)+kk
    !
  end function indezw
  !
  function indezp(kk,jj,ii)
    !
    implicit none
    !
    integer             :: indezp
    integer, intent(in) :: ii,jj,kk
    !
    indezp = ((lk+1)*(nj+1))*(ii-1)+(lk+1)*(jj-1)+kk
    !
  end function indezp
  !
end MODULE malla
