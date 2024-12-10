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
  !
  implicit none
  !
  ! Arreglos para la condici\'on de frontera de u
  !
  real(kind=DBL), dimension(2,nj+1) :: cond_front_ux
  real(kind=DBL), dimension(2,mi)   :: cond_front_uy
  !
  ! Arreglos para la condici\'on de frontera de v
  !
  real(kind=DBL), dimension(2,nj)   :: cond_front_vx
  real(kind=DBL), dimension(2,mi+1) :: cond_front_vy
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
       &cond_front_uux, cond_front_uuy,nn&
       &)
    !
    implicit none
    !
    character(len=18), intent(in)                 :: entrada_front_uuo
    !
    real(kind=DBL), dimension(2,mm), intent(in)   :: cond_front_uux
    real(kind=DBL), dimension(2,nn), intent(in)   :: cond_front_uuy
    !
    character(len=1) :: segmento
    !
    real(kind=DBL)   :: x0, x1, y0, y1
    !
    !***************************************************************
    !
    ! Inicializaci\'on de los arreglos para condiciones de frontera
    !
    cond_front_uux = 0.0_DBL
    cond_front_uuy = 0.0_DBL
    !
    !
    ! Apertura de los archivos y lectura de los arreglos
    !
    open(unit=111,file=entrada_front_uuo)
    read(111,*) segmento, x0, x1
    do jj = 1, nj+1
       do ii = 1, mi
          read(111,*) xuo(ii),ypo(jj),u_anto(ii,jj)
       end do
    end do
    close(unit=111)
    !
  end subroutine lectura_cond_frontera
  !
end MODULE cond_frontera
