!
! M\'odulo de frontera  inmersa
!
! Este m\'odulo contiene las subrutinas que permiten incluir s\'olidos en el dominio del fluido 
! Se utiliza la media geom\'etrica de los coeficientes de difusi\'on en las caras de los
! vol\'umenes de control como se describe en Patankar 1987.
!
! Autor: J.C. Cajas
!
module frontera_inmersa
  !
  use malla, only : mi, nj, DBL     ! dimensiones de la malla y tipo de doble
  use malla, only : xp, yp          ! coordenadas en la malla de la presi\'on
  use ec_momento, only : u, u_ant   ! provisional para probar frontera inmersa
  !
  implicit none
  !
  ! Interfaces para poder usar funciones como argumentos en subrutinas
  !
  abstract interface
     !
     function funcion1d(x) result(y)
       !
       import DBL
       implicit none
       !
       real(kind=DBL), intent(in) :: x
       real(kind=DBL)             :: y
     end function funcion1d
     !
  end interface
  !  
contains
  !
  !*******************************************************************
  !
  ! definir_cuerpo
  !
  ! Subrutina que define la regi\'on de la malla que est\'a dentro de
  ! un cuerpo s\'olido o no. Asigna los valores correspondientes a las
  ! constantes de difusi\'on gamma.
  !
  !*******************************************************************
  !
  subroutine definir_cuerpo(gamma_momeno, gamma_enero, opcion)
    !
    ! import xp, yp, DBL
    ! import en_region_tipo1
    ! import recta_m1xb1, recta_mxb2
    implicit none
    !
    real(kind=DBL), dimension(mi+1,nj+1), intent(inout) :: gamma_momeno, gamma_enero
    character(5), intent(in)                            :: opcion
    !
    integer         :: ii, jj
    real(kind=DBL)  :: xv, yv, hh
    real(kind=DBL)  :: m1, m2 ! Pend. tri\'angulo
    ! procedure(fun    :: f1, f2
    !
    ! Mensaje
    !
    print*, " Se define una frontera inmersa con la opcion: ",opcion
    if( opcion == 'trian' )then
       !
       ! Tri\'angulo con v\'ertice en p(xv,yv) y altura h
       ! definido usando las funciones recta_mxb1 y recta_mxb2
       !
       yv = 0.2_DBL
       hh = 0.2_DBL
       do jj = 1, nj+1
          if( yv-hh <= yp(jj) .and. yp(jj) <= yv )then
             do ii = 1, mi+1
                ! xx = yp(jj) Se utiliza como una regi\'on de tipo II
                ! yy = xp(ii)
                if( en_region_tipo1(yp(jj),xp(ii),yv-hh,yv,recta_mxb1,recta_mxb2) )then
                   ! print*, "DEBUG: dentro", xp(ii), yp(jj)
                   gamma_momeno(ii,jj) = 10.0e6_DBL
                   u(ii,jj)     = 0.0_DBL
                   u_ant(ii,jj) = 0.0_DBL
                end if
             end do
          end if
       end do
       !
    end if
       !
  end subroutine definir_cuerpo
  !
  !
  !*******************************************************************
  !
  ! en_region_tipo1
  !
  ! Funci\'on para decidir si un punto est\'a en una regi\'on en el
  ! plano entre dos valores fijos a <= x <= b y dos funciones
  ! f1(x) <= y <= f2(x)
  !
  !*******************************************************************
  !
  function en_region_tipo1(xx,yy,aa,bb,f1,f2)
    !
    ! import DBL
    ! implicit none
    !
    logical                    :: en_region_tipo1 ! true para dentro de la region
    real(kind=DBL), intent(in) :: xx, yy          ! coordenadas del punto de prueba
    real(kind=DBL), intent(in) :: aa, bb          ! extremos para x
    !
    procedure(funcion1d)       :: f1, f2          ! extremos para y
    !
    en_region_tipo1=.false.
    !
    if( aa <= xx .and. xx <= bb )then
       !
       if( f1(xx) <= yy .and. yy<= f2(xx) )then
          !
          en_region_tipo1 = .true.
          !
       end if
       !
    end if
    !    
  end function en_region_tipo1
  !*******************************************************************
  !
  ! funciones para definir cuerpos inmersos
  !
  !*******************************************************************
  !
  ! Recta de pendiente mm y ordenada al origen bb
  !
  real(kind=DBL) function recta_mxb1(xx)
    !
    ! import DBL
    ! implicit none
    !
    real(kind=DBL), intent(in) :: xx
    real(kind=DBL)             :: mm, bb
    !
    mm = 1.0_DBL/4.0_DBL
    bb =-( 0.2_DBL - mm * 6.0_DBL ) / mm
    recta_mxb1 = xx/mm + bb
    !
  end function recta_mxb1
  !
  ! Recta de pendiente mm y ordenada al origen bb
  !
  real(kind=DBL) function recta_mxb2(xx)
    !
    ! import DBL
    implicit none
    !
    real(kind=DBL), intent(in) :: xx
    real(kind=DBL)             :: mm, bb
    !
    mm =-1.0_DBL/4.0_DBL
    bb =-( 0.2_DBL - mm * 6.0_DBL ) / mm
    recta_mxb2 = xx/mm + bb
    !
  end function recta_mxb2
  !  
end module frontera_inmersa
