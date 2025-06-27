!
!
! M\'odulo para la ecuaci\'on de la energ\'ia
!
! Este m\'odulo contiene:
!
! - Las variables de la ecuaci\'on de la energ\'ia
! - Las subrutinas que calculan los coeficientes de las ecuaciones discretizadas
! - Las subrutinas que sirven para leer e imponer condiciones de frontera
!
module ec_energia
  !
  use malla, only : mi, nj, lk, DBL
  use malla, only : indexp, indeyp, indezp
  !
  ! use cond_frontera, only : tipo_cond_front
  ! !
  ! use cond_frontera, only : inicializa_cond_front
  ! use cond_frontera, only : lectura_cond_frontera
  !
  implicit none
  !
  ! Variables para la ecuaci\'on de la energ\'ia
  !
  real(kind=DBL), dimension(mi+1,nj+1,lk+1) :: temp       ! Temperatura
  real(kind=DBL), dimension(mi+1,nj+1,lk+1) :: temp_ant   ! Temperatura del paso anterior
  real(kind=DBL), dimension(mi+1,nj+1,lk+1) :: ftemp      ! Diferencia de temperatura
  !
  real(kind=DBL), dimension(mi+1,nj+1,lk+1) :: gamma_ener ! Coeficiente de difusi\'on
  !
  real(kind=DBL) :: conv_t    ! Criterio de convergencia
  real(kind=DBL) :: rel_ener  ! Par\'ametro de relajaci\'on
  !
  ! Estructuras para guardar la informaci\'on de las condiciones de frontera
  !
  ! type( tipo_cond_front ) :: cond_front_ua, cond_front_ub, cond_front_uc, cond_front_ud
  ! type( tipo_cond_front ) :: cond_front_va, cond_front_vb, cond_front_vc, cond_front_vd
  !
contains
  !
  !*******************************************************************
  !
  ! ensambla_ec_energ_x
  !
  ! Subrutina que calcula los coeficientes de la matriz tridiagonal
  ! para la ec. de la energ\'ia en direcci\'on x 
  !
  !*******************************************************************
  subroutine ensambla_energ_x(&
       &deltaxpo,&
       &deltaypo,&
       &deltazpo,&
       &deltaxuo,&
       &deltayvo,&
       &deltazwo,&
       &fexpo,&
       &feypo,&
       &fezpo,&
       &temper_o,&
       &u_o,v_o,w_o,&
       &rel_energo,&
       &AI_o,AC_o,AD_o,Rx_o,&
       &ii,jj,kk)
    implicit none
    !$acc routine
    !
    ! Tama\~no del volumen de control
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxpo
    real(kind=DBL), dimension(nj), intent(in) :: deltaypo
    real(kind=DBL), dimension(lk), intent(in) :: deltazpo
    !
    ! Distancia entre nodos contiguos de la malla de p en direcci\'on horizontal x
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxuo
    !
    ! Distancia entre nodos contiguos de la malla de p en direcci\'on horizontal y
    !
    real(kind=DBL), dimension(nj), intent(in) :: deltayvo
    !
    !
    ! Distancia entre nodos contiguos de la malla de p en direcci\'on vertical z
    !
    real(kind=DBL), dimension(nj), intent(in) :: deltazwo
    !
    ! Coeficientes de interpolaci\'on
    !
    real(kind=DBL), DIMENSION(mi), intent(in) :: fexpo
    real(kind=DBL), DIMENSION(nj), intent(in) :: feypo
    real(kind=DBL), DIMENSION(nj), intent(in) :: fezpo
    !
    ! Velocidad, presi\'on, t\'ermino fuente b
    !
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in)  :: u_o
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in)  :: v_o
    real(kind=DBL), dimension(mi+1,nj+1,lk),   intent(in)  :: w_o
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in)  :: temper_o
    !
    ! coeficiente de relajaci\'on
    !
    real(kind=DBL), intent(in) :: rel_energo
    !
    ! Coeficientes de las matrices
    !
    ! ** Estos coeficientes est\'an sobredimensionados para reducir el uso de memoria
    ! en la gpu, los arreglos que se reciben en esta subrutina se usan para las ecs.
    ! de momento en, energ\'ia y la correcci\'on de la presi\'on **
    !
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AI_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AC_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AD_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: RX_o
    !
    ! \'Indice para recorrer las direcciones x, y. z
    !
    integer,   intent(in) :: ii, jj, kk
    !
    ! Variables auxiliares
    !
    ! Auxiliares de interpolaci\'on
    !
    real(kind=DBL) :: ui, ud, vs, vn, wt, wb
    real(kind=DBL) :: di, dd, ds, dn, db, da
    real(kind=DBL) :: gammai, gammad
    real(kind=DBL) :: gammas, gamman
    real(kind=DBL) :: gammab, gammat
    real(kind=DBL) :: alpha,beta,gamma,delta
    !
    ! Interpolaciones necesarias
    !
    ! u
    !
    ud = u_o(ii,jj,kk)
    ui = u_o(ii-1,jj,kk)
    !
    ! v
    !
    vn = v_o(ii,jj,kk)
    vs = v_o(ii,jj-1,kk)
    !
    ! w
    !
    wt = w_o(ii,jj,kk)
    wb = w_o(ii,jj,kk-1)
    !
    ! distancias entre nodos contiguos
    !
    di = deltaxuo(ii-1)
    dd = deltaxuo(ii)
    ds = deltayvo(jj-1)
    dn = deltayvo(jj)
    db = deltazwo(kk-1)
    da = deltazwo(kk)
    !
    ! Coeficientes de difusi\'on
    !
    gammad = ( gamma_ener(ii+1,jj,kk) * gamma_ener(ii,jj,kk) ) / &
         &( gamma_ener(ii+1,jj,kk) * (1._DBL-fexpo(ii))+gamma_ener(ii,jj,kk)*fexpo(ii) )
    gammai = ( gamma_ener(ii,jj,kk) * gamma_ener(ii-1,jj,kk) ) / &
         &( gamma_ener(ii,jj,kk) * (1._DBL-fexpo(ii-1))+gamma_ener(ii-1,jj,kk)*fexpo(ii-1) )
    gamman = ( gamma_ener(ii,jj+1,kk) * gamma_ener(ii,jj,kk) ) / &
         &( gamma_ener(ii,jj+1,kk) * (1._DBL-feypo(jj))+gamma_ener(ii,jj,kk)*feypo(jj) )
    gammas = ( gamma_ener(ii,jj,kk) * gamma_ener(ii,jj-1,kk) ) / &
         &( gamma_ener(ii,jj,kk) * (1._DBL-feypo(jj-1))+gamma_ener(ii,jj-1,kk)*feypo(jj-1) )
    gammat = ( gamma_ener(ii,jj,kk+1) * gamma_ener(ii,jj,kk) ) / &
         &( gamma_ener(ii,jj,kk+1) * (1._DBL-fezpo(kk))+gamma_ener(ii,jj,kk)*fezpo(kk) )
    gammat = ( gamma_ener(ii,jj,kk) * gamma_ener(ii,jj,kk-1) ) / &
         &( gamma_ener(ii,jj,kk) * (1._DBL-fezpo(kk-1))+gamma_ener(ii,jj,kk-1)*fezpo(kk-1) )
    !
    ! Tama\~no de los vol\'umenes de control para la velocidad u
    !
    ! delta_x = deltaxpo(ii)
    ! delta_y = deltaypo(jj)
    ! delta_z = deltazpo(kk)
    !
    ! -------------------------
    !
    ! Coeficientes de la matriz
    !
    AI_o(indexp(ii,jj,kk)) =-gammai * deltaypo(jj)*deltazpo(kk)/di
    AD_o(indexp(ii,jj,kk)) =-gammad * deltaypo(jj)*deltazpo(kk)/dd
    alpha                  =-gammas * deltaxpo(ii)*deltazpo(kk)/ds
    beta                   =-gamman * deltaxpo(ii)*deltazpo(kk)/dn
    gamma                  =-gammab * deltaxpo(ii)*deltaypo(jj)/db
    delta                  =-gammat * deltaxpo(ii)*deltaypo(jj)/da
    !
    AC_o(indexp(ii,jj,kk)) = ( -AI_o(indexp(ii,jj,kk)) - AD_o(indexp(ii,jj,kk))-&
         &alpha - beta - gamma - delta ) / rel_energo
    !
    Rx_o(indexp(ii,jj,kk)) =-alpha-&
         &beta +&
         &gamma+&
         &delta+&
         &(1._DBL-rel_energo)*AC_o(indexp(ii,jj,kk))

  end subroutine ensambla_energ_x
  !
end module ec_energia
