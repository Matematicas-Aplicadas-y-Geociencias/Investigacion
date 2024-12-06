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
  use malla, only : mi, nj, DBL       ! dimensiones de la malla y tipo de doble
  use malla, only : xp, yp            ! coordenadas en la malla de la presi\'on
  !
  ! Variables para fijar las fronteras inmersas
  !
  use ec_momento, only : u, u_ant, v, v_ant
  use ec_momento, only : au, av  
  use ec_momento, only : fuente_con_u, fuente_lin_u
  use ec_momento, only : fuente_con_v, fuente_lin_v
  !
  ! Variables para fijar temperatura de las fronteras inmersas
  !
  use ec_energia, only : fuente_con_t, fuente_lin_t
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
    !
    if( opcion == 'trian' )then
       !
       ! Tri\'angulo con v\'ertice en p(xv,yv) y altura hh
       ! definido usando las funciones recta_mxb1 y recta_mxb2
       !
       yv = 0.1_DBL
       hh = 0.1_DBL
       do jj = 1, nj+1
          if( yv-hh <= yp(jj) .and. yp(jj) <= yv )then
             do ii = 1, mi+1
                ! xx = yp(jj) Se utiliza como una regi\'on de tipo II
                ! yy = xp(ii)
                if( en_region_tipo1(yp(jj),xp(ii),yv-hh,yv,recta_mxb1,recta_mxb2) )then
                   ! print*, "DEBUG: dentro", xp(ii), yp(jj)
                   !
                   fuente_lin_u(ii,jj) =-10.0e50_DBL
                   fuente_con_u(ii,jj) = 10.0e50_DBL*10.0e-12_DBL
                   !
                   fuente_lin_v(ii,jj) =-10.0e50_DBL
                   fuente_con_v(ii,jj) = 10.0e50_DBL*10.0e-12_DBL 
                   !
                end if
             end do
          end if
       end do
       !
    else if( opcion == 'cuadr' )then
       !
       ! Cuadrado con centro en xv, yv y lado hh
       !
       yv = 0.05_DBL
       xv = 6.0_DBL
       hh = 0.1_DBL
       !
       do jj = 1, nj+1
          if( yv-hh/2.0_DBL <= yp(jj) .and. yp(jj) <= yv+hh/2.0_DBL )then
             do ii = 1, mi+1
                ! xx = yp(jj) Se utiliza como una regi\'on de tipo II
                ! yy = xp(ii)
                if( xv-hh/2.0_DBL <= xp(ii) .and. xp(ii) <= xv+hh/2.0_DBL )then
                   ! gamma_momeno(ii,jj) = 10.0e6_DBL
                   !
                   ! u(ii,jj)            = 10.0_DBL
                   ! v(ii,jj)            = 5.0_DBL
                   !                   
                   ! u_ant(ii,jj)        = 10.0_DBL
                   ! v_ant(ii,jj)        = 5.0_DBL
                   !
                   fuente_lin_u(ii,jj) =-10.0e50_DBL
                   fuente_con_u(ii,jj) = 10.0e50_DBL*10.0e-12_DBL
                   !
                   fuente_lin_v(ii,jj) =-10.0e50_DBL
                   fuente_con_v(ii,jj) = 10.0e50_DBL*10.0e-12_DBL 
                   !
                   !print*, "DEBUG: dentro", ii,jj, u(ii,jj)
                end if
                !
             end do
             !
          end if
          !
       end do
       !
    else if( opcion == 'chime' )then
       !
       ! varillas inmersas para chimenea solar
       !
       yv = 0.5_DBL
       xv = 8.5_DBL
       hh = 0.3_DBL
       !
       do jj = 1, nj+1
          if( yv-hh/2.0_DBL <= yp(jj) .and. yp(jj) <= yv+hh/2.0_DBL )then
             do ii = 1, mi+1
                ! xx = yp(jj) Se utiliza como una regi\'on de tipo II
                ! yy = xp(ii)
                if( xv-hh/2.0_DBL <= xp(ii) .and. xp(ii) <= xv+hh/2.0_DBL )then
                   ! gamma_momeno(ii,jj) = 10.0e6_DBL
                   !
                   ! u(ii,jj)            = 10.0_DBL
                   ! v(ii,jj)            = 5.0_DBL
                   !                   
                   ! u_ant(ii,jj)        = 10.0_DBL
                   ! v_ant(ii,jj)        = 5.0_DBL
                   !
                   fuente_lin_u(ii,jj) =-10.0e50_DBL
                   fuente_con_u(ii,jj) = 10.0e50_DBL*10.0e-12_DBL
                   au(ii,jj)           = 10.0e40_DBL
                   !
                   fuente_lin_v(ii,jj) =-10.0e50_DBL
                   fuente_con_v(ii,jj) = 10.0e50_DBL*10.0e-12_DBL
                   av(ii,jj)           = 10.0e40_DBL
                   !
                   fuente_lin_t(ii,jj) =-10e50_DBL
                   fuente_con_t(ii,jj) = 10e50_DBL*1.0_DBL
                   !print*, "DEBUG: dentro", ii,jj, u(ii,jj)
                end if
                !
             end do
             !
          end if
          !
       end do
       !       
    end if ! Selecci\'on del caso del cuerpo inmerso
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
    mm = 1.0_DBL/2.0_DBL
    bb =-( 0.1_DBL - mm * 6.0_DBL ) / mm
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
    mm =-1.0_DBL/2.0_DBL
    bb =-( 0.1_DBL - mm * 6.0_DBL ) / mm
    recta_mxb2 = xx/mm + bb
    !
  end function recta_mxb2
  !  
end module frontera_inmersa
