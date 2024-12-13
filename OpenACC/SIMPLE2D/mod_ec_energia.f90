!
!
! M\'odulo para la ecuaci\'on de energ\'ia
!
! Este m\'odulo contiene:
! - Las variables de la ecuaci\'on de la energ\'ia
! - Las subrutinas que calculan los coeficientes de la ecuaci\'on discretizada
!
module ec_energia
  !
  use malla, only : mi, nj, DBL
  use malla, only : xp, yp
  !
  use cond_frontera, only : tipo_cond_front
  use cond_frontera, only : inicializa_cond_front
  use cond_frontera, only : lectura_cond_frontera
  !
  implicit none
  !
  ! Variables para la ecuaci\'on de la energ\'ia
  ! temperatura, coeficiente de difusi\'on y criterios de convergencia
  ! 
  REAL(kind=DBL), DIMENSION(mi+1,nj+1) :: temp, temp_ant, dtemp, Restemp
  REAL(kind=DBL), DIMENSION(mi+1,nj+1) :: ftemp, gamma_energ
  REAL(kind=DBL)  :: conv_t,rel_ener
  !
  ! Variables para los t\'erminos fuente de la ec de energ\'ia.
  ! El t\'ermino lineal fuente_lin debe ser negativo para favorecer
  ! convergencia y estabilidad
  real(kind=DBL), dimension(mi+1,nj+1)   :: fuente_con_t, fuente_lin_t
  !
  ! Estructuras para guardar la informaci\'on de las condiciones de frontera
  !
  type( tipo_cond_front ) :: cond_front_ta, cond_front_tb, cond_front_tc, cond_front_td
  !
contains
  !
  !*******************************************************************
  !
  ! ini_frontera_t
  !
  ! Subrutina que inicializa los arreglos para las condiciones
  ! de frontera de t, lee las condiciones de los archivos
  ! de entrada cond_fronterat.dat
  !
  !*******************************************************************  
  subroutine ini_frontera_t()
    !
    implicit none
    !
    ! arreglos de u
    !
    call inicializa_cond_front(cond_front_ta)
    call inicializa_cond_front(cond_front_tb)
    call inicializa_cond_front(cond_front_tc)
    call inicializa_cond_front(cond_front_td)
    call lectura_cond_frontera('cond_fronterat.dat',&
         & cond_front_ta, &
         & cond_front_tb, &
         & cond_front_tc, &
         & cond_front_td, &
         & xp, yp,        &
         & mi+1, nj+1     &
         & )
  end subroutine ini_frontera_t
  !
  !*******************************************************************
  !
  ! ensambla_energia
  !
  ! Subrutina que calcula los coeficientes de la matriz tridiagonal
  ! para la ecuacion de la energ\'ia
  !
  !*******************************************************************
  subroutine ensambla_energia(deltaxpo,deltaypo,&
       &deltaxuo,deltayvo,fexpo,feypo,gamma_ener_o,&
       &fuente_con_to,fuente_lin_to,&
       &u_o,v_o,&
       &temp_o,temp_anto,dt_o,&
       &rel_ener,placa_mino,placa_maxo,&
       &AI_o,AC_o,AD_o,Rx_o,BS_o,BC_o,BN_o,Ry_o,&
       &ii,jj&
       &)
    implicit none
    !$acc routine
    !
    ! Tama\~no del volumen de control
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxpo
    real(kind=DBL), dimension(nj), intent(in) :: deltaypo
    !
    ! Distancia entre nodos contiguos de la malla de p en direcci\'on horizontal
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxuo
    !
    ! Distancia entre nodos contiguos de la malla de p en direcci\'on vertical
    !
    real(kind=DBL), dimension(nj), intent(in) :: deltayvo
    !
    ! Coeficientes de interpolaci\'on
    !
    real(kind=DBL), DIMENSION(mi), intent(in) :: fexpo    
    real(kind=DBL), DIMENSION(nj), intent(in) :: feypo
    !
    ! Coeficiente de difusi\'on
    !
    real(kind=DBL), dimension(mi+1,nj+1), intent(in) ::  gamma_ener_o
    !
    ! Velocidad, presi\'on, temperatura
    !
    real(kind=DBL), dimension(mi,nj+1),   intent(in) :: u_o
    real(kind=DBL), dimension(mi+1,nj),   intent(in) :: v_o
    real(kind=DBL), dimension(mi+1,nj+1), intent(in) :: temp_o, temp_anto
    !
    ! T\'erminos fuente
    !
    real(kind=DBL), dimension(mi+1,nj+1), intent(in) :: fuente_con_to
    real(kind=DBL), dimension(mi+1,nj+1), intent(in) :: fuente_lin_to 
    !
    ! paso de tiempo y coeficiente de relajaci\'on
    !
    real(kind=DBL), intent(in) :: dt_o
    real(kind=DBL), intent(in) :: rel_ener
    !
    ! l\'imites de las placas calientes
    !
    integer,        intent(in) :: placa_mino, placa_maxo
    integer,        intent(in) :: ii,jj
    !
    ! Coeficientes de las matrices
    !
    ! ** Estos coeficientes est\'an sobredimensionados para reducir el uso de memoria
    ! en la gpu, los arreglos que se reciben en esta subrutina se usan para las ecs.
    ! de momento en, energ\'ia y la correcci\'on de la presi\'on **
    !
    real(kind=DBL), dimension(mi+1,nj+1), intent(out) :: AI_o, AC_o, AD_o, Rx_o
    real(kind=DBL), dimension(nj+1,mi+1), intent(out) :: BS_o, BC_o, BN_o, Ry_o
    !
    ! Variables auxiliares
    !
    ! integer :: ii, jj
    !
    ! Auxiliares de interpolaci\'on
    !
    real(kind=DBL) :: ui, ud, vs, vn
    real(kind=DBL) :: di, dd, ds, dn
    real(kind=DBL) :: gammai, gammad
    real(kind=DBL) :: gammas, gamman
    !
    ! C\'alculo de los coeficientes
    !
    ! $acc loop gang
    !bucle_direccion_y: do jj = 2, nj
    !------------------------
    ! Condiciones de frontera
    ! AC_o(1,jj) =-1._DBL
    ! AD_o(1,jj) = 1._DBL
    ! Rx_o(1,jj) = 0._DBL
    !
    ! Llenado de la matriz
    !
    ! $acc loop vector
    !bucle_direccion_x: do ii = 2, mi
    !
    ! Interpolaciones necesarias
    !
    ! u
    !
    ud = u_o(ii,jj)
    ui = u_o(ii-1,jj)
    !
    ! v
    !
    vn = v_o(ii,jj)
    vs = v_o(ii,jj-1)
    !
    ! distancias entre nodos contiguos
    !
    di = deltaxuo(ii-1)
    dd = deltaxuo(ii)
    ds = deltayvo(jj-1)
    dn = deltayvo(jj)
    !
    ! Coeficientes de difusi\'on
    !
    gammad = ( gamma_ener_o(ii+1,jj) * gamma_ener_o(ii,jj) ) / &
         &( gamma_ener_o(ii+1,jj) * (1._DBL-fexpo(ii))+gamma_ener_o(ii,jj)*fexpo(ii) )
    gammai = ( gamma_ener_o(ii,jj) * gamma_ener_o(ii-1,jj) ) / &
         &( gamma_ener_o(ii,jj) * (1._DBL-fexpo(ii-1))+gamma_ener_o(ii-1,jj)*fexpo(ii-1) )
    gamman = ( gamma_ener_o(ii,jj+1) * gamma_ener_o(ii,jj) ) / &
         &( gamma_ener_o(ii,jj+1) * (1._DBL-feypo(jj))+gamma_ener_o(ii,jj)*feypo(jj) )
    gammas = ( gamma_ener_o(ii,jj) * gamma_ener_o(ii,jj-1) ) / &
         &( gamma_ener_o(ii,jj) * (1._DBL-feypo(jj-1))+gamma_ener_o(ii,jj-1)*feypo(jj-1) )
    !
    ! Tama\~no de los vol\'umenes de control para la energia
    !
    ! delta_x = deltaxpo(ii)
    ! delta_y = deltaypo(jj)
    !
    ! -------------------------
    !
    ! Coeficientes de la matriz
    !
    AI_o(ii,jj) =-(gammai*deltaypo(jj)/di*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ui*di/gammai))**5)+&
         &DMAX1(0.0_DBL,ui*deltaypo(jj)))
    AD_o(ii,jj) =-(gammad*deltaypo(jj)/dd*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ud*dd/gammad))**5)+&
         &DMAX1(0.0_DBL,-ud*deltaypo(jj)))
    BS_o(jj,ii) =-(gammas*deltaxpo(ii)/ds*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vs*ds/gammas))**5)+&
         &DMAX1(0.0_DBL, vs*deltaxpo(ii)))
    BN_o(jj,ii) =-(gamman*deltaxpo(ii)/dn*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vn*dn/gamman))**5)+&
         &DMAX1(0.0_DBL,-vn*deltaxpo(ii)))
    AC_o(ii,jj) = ( -AI_o(ii,jj) - AD_o(ii,jj) - BS_o(jj,ii) - BN_o(jj,ii)-&
         &deltaxpo(ii)*deltaypo(jj)*fuente_lin_to(ii,jj) +&
         &deltaxpo(ii)*deltaypo(jj)/dt_o ) / rel_ener
    Rx_o(ii,jj) =-BS_o(jj,ii)*temp_o(ii,jj-1) - BN_o(jj,ii)*temp_o(ii,jj+1)+&
         &deltaxpo(ii)*deltaypo(jj)*fuente_con_to(ii,jj) +&
         &deltaxpo(ii)*deltaypo(jj)*temp_anto(ii,jj)/dt_o+&
         &AC_o(ii,jj)*(1._DBL-rel_ener)*temp_o(ii,jj)
    BC_o(jj,ii) = AC_o(ii,jj)
    Ry_o(jj,ii) =-AI_o(ii,jj)*temp_o(ii-1,jj) - AD_o(ii,jj)*temp_o(ii+1,jj)+&
         &deltaxpo(ii)*deltaypo(jj)*fuente_con_to(ii,jj) +&
         &deltaxpo(ii)*deltaypo(jj)*temp_anto(ii,jj)/dt_o+&
         &BC_o(jj,ii)*(1._DBL-rel_ener)*temp_o(ii,jj)
    ! end do bucle_direccion_x
    !------------------------
    ! Condiciones de frontera
    ! AI_o(mi+1,jj) =-1._DBL 
    ! AC_o(mi+1,jj) = 1._DBL
    ! Rx_o(mi+1,jj) = 0.0_DBL
    ! ---------------------

  end subroutine ensambla_energia
  !
end module ec_energia
