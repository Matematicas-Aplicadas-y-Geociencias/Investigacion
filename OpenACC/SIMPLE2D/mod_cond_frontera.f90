!
!
! M\'odulo de herramientas para la lectura de condiciones de frontera del programa SIMPLE 2D
!
! Este m\'odulo contiene herramientas que se utilizan para leer las condiciones de frontera
! en un dominio rectangular
!
MODULE cond_frontera
  !
  use malla, only : mi, nj, DBL
  use malla, only : 
  !
  implicit none
  !
  ! Arreglos  para imponer las condiciones en las matrices
  !
  real(kind=DBL), dimension(nj+1,2) :: ai_u,ac_u,ad_u
  real(kind=DBL), dimension(mi,2)   :: bs_u,bc_u,bn_u
  !
  ! Arreglos para la condici\'on de frontera de u
  !
  real(kind=DBL), dimension(2,nj+1) :: cond_front_ux
  real(kind=DBL), dimension(2,mi)   :: cond_front_uy
  !
  ! Arreglos  para imponer las condiciones en las matrices
  !
  real(kind=DBL), dimension(nj,2)   :: ai_v,ac_v,ad_v
  real(kind=DBL), dimension(mi+1,2) :: bs_v,bc_v,bn_v
  !
  ! Arreglos para la condici\'on de frontera de v
  !
  real(kind=DBL), dimension(2,nj)   :: cond_front_vx
  real(kind=DBL), dimension(2,mi+1) :: cond_front_vy
  !
  ! Arreglos  para imponer las condiciones en las matrices
  !
  real(kind=DBL), dimension(nj+1,2) :: ai_t,ac_t,ad_t
  real(kind=DBL), dimension(mi+1,2) :: bs_t,bc_t,bn_t
  !
  ! Arreglos para la condici\'on de frontera de temp
  !
  real(kind=DBL), dimension(2,nj+1) :: cond_front_tx
  real(kind=DBL), dimension(2,mi+1) :: cond_front_ty
  !
  ! Nombres de los archivos para condiciones de frontera
  !
  character(len=18) :: entrada_front_u
  character(len=18) :: entrada_front_v
  character(len=18) :: entrada_front_t
  !
contains
  !
  !-----------------------------------------------------------------------------
  !
  ! Esta subrutina abre los archivos para leer las condiciones de frontera
  ! y crea los arreglos necesarios para imponerlas durante la construcci\'on
  ! de las matrices
  !
  subroutine lectura_cond_frontera(entrada_front_uuo,&
       &xx,yy,&
       &mm,nn,&
       &cond_front_uux,&
       &cond_front_uuy,&
       &ai,ac,ad,&
       &bs,bc,bn &
       &)
    !
    implicit none
    !
    character(len=18), intent(in)                :: entrada_front_uuo
    integer, intent(in)                          :: mm,nn
    !
    real(kind=DBL), dimension(mm),   intent(in)  :: xx
    real(kind=DBL), dimension(nn),   intent(in)  :: yy
    real(kind=DBL), dimension(nn,2), intent(out) :: cond_front_uux
    real(kind=DBL), dimension(mm,2), intent(out) :: cond_front_uuy
    real(kind=DBL), dimension(nn,2), intent(out) :: ai,ac,ad
    real(kind=DBL), dimension(mm,2), intent(out) :: bs,bc,bn
    !
    character(len=1) :: lado
    character(len=4) :: variable
    character(len=4) :: tipo_condicion
    !
    real(kind=DBL)   :: x0, x1, valor
    !
    integer          :: divisiones
    integer          :: indice
    integer          :: ii, jj, kk
    !
    !***************************************************************
    !
    ! Inicializaci\'on de los arreglos para condiciones de frontera
    !
    cond_front_uux = 0.0_DBL
    cond_front_uuy = 0.0_DBL
    !
    ! Inicializaci\'on de las variables para recorrer segmentos
    !
    jj = 1
    indice = 1 ! Distingue entre los lados a y c, o b y d
    !
    ! Apertura de los archivos y lectura de los arreglos
    !
    open(unit=111,file=entrada_front_uuo)
    !
    ! Se esperan 4 lados en el rect\'angulo
    !
    do kk = 1, 4
       !
       read(111,*) variable, lado, divisiones
       !
       write(*,*) "SIMPLE2D: Condiciones de frontera para ", variable, &
            &" con ", divisiones," divisiones en lado ", lado
       !
       lado: if( lado == 'a' .or. lado == 'c' )then
          !
          if( lado == 'c' ) indice = 2
          !
          divisiones_x: do ii = 1, divisiones
             !
             read(111,*) x0, x1, tipo_condicion, valor
             !
             write(*,*) "SIMPLE 2D: condici\'on tipo ", tipo_condicion, &
                  &" entre ", x0, " y ", x1," en lado ", lado
             !
             tipo_condicion_x: select case( tipo_condicion )
                !
             case( 'diri' )
                !
                do while( x0 <= xx(jj) .and. xx(jj) <= x1 )
                   ai(jj,indice) = 0.0_DBL
                   ac(jj,indice) = 1.0_DBL
                   ad(jj,indice) = 0.0_DBL
                   cond_front_uux(jj,indice) = valor
                   jj = jj + 1
                end do
                !
             case ( 'neum' )
                !
                if( lado == 'a' )then
                   !
                   do while( x0 <= xx(jj) .and. xx(jj) <= x1 )
                      ai(jj,indice) = 0.0_DBL
                      ac(jj,indice) =-1.0_DBL
                      ad(jj,indice) = 1.0_DBL
                      cond_front_uux(jj,indice) = valor
                      jj = jj + 1
                   end do
                   !
                else
                   !
                   do while( x0 <= xx(jj) .and. xx(jj) <= x1 )
                      ai(jj,indice) =-1.0_DBL
                      ac(jj,indice) = 1.0_DBL
                      ad(jj,indice) = 0.0_DBL
                      cond_front_uux(jj,indice) = valor
                      jj = jj + 1
                   end do
                   !
                end if
                !
             end select tipo_condicion_x
             !
             ! Se devuelve el valor de 1 al indice para
             ! la condici\'on en la otra direcci\'on
             !
             indice = 1
             !
          end do divisiones_x
          !
       else if( lado == 'b' .or. lado == 'd' )then
          !
          if( lado == 'd' ) indice = 2
          !
          divisiones_y: do ii = 1, divisiones
             !
             read(111,*) x0, x1, tipo_condicion, valor
             write(*,*) "SIMPLE 2D: condici\'on tipo ", tipo_condicion, &
                  &" entre ", x0, " y ", x1," en lado ", lado
             !
             tipo_condicion_y: select case( tipo_condicion )
                !
             case( 'diri' )
                !
                do while( x0 <= xx(jj) .and. xx(jj) <= x1 )
                   !
                   bs(jj,indice) = 0.0_DBL
                   bc(jj,indice) = 1.0_DBL
                   bn(jj,indice) = 0.0_DBL
                   cond_front_uuy(jj,indice) = valor
                   jj = jj + 1
                   !
                end do
                !
             case ( 'neum' )
                !
                if( lado == 'b' )then
                   !
                   do while( x0 <= xx(jj) .and. xx(jj) <= x1 )
                      bs(jj,indice) = 0.0_DBL
                      bc(jj,indice) =-1.0_DBL
                      bn(jj,indice) = 1.0_DBL
                      cond_front_uuy(jj,indice) = valor
                      jj = jj + 1
                   end do
                   !
                else
                   !
                   do while( x0 <= xx(jj) .and. xx(jj) <= x1 )
                      bs(jj,indice) =-1.0_DBL
                      bc(jj,indice) = 1.0_DBL
                      bn(jj,indice) = 0.0_DBL
                      cond_front_uuy(jj,indice) = valor
                      jj = jj + 1
                   end do
                   !                  
                end if
                !
             end select tipo_condicion_y
             !
             ! Se devuelve el valor de 1 al indice para
             ! la condici\'on en la otra direcci\'on
             !
             indice = 1
             !
          end do divisiones_y
          !
       end if lado
       !
    end do
    !
    close(unit=111)
    !
  end subroutine lectura_cond_frontera
  !
end MODULE cond_frontera
