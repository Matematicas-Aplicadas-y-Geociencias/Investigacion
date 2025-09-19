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
  ! Variables para los t\'erminos fuente de la ec de la energ\'ia
  ! El t\'ermino lineal fuente_lin debe ser negativo para favorecer
  ! convergencia y estabilidad
  !
  real(kind=DBL), dimension(mi+1,nj+1,lk+1)   :: fuente_con_temp, fuente_lin_temp
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
  subroutine ensambla_energia_x(&
       &deltaxpo,&
       &deltaypo,&
       &deltazpo,&
       &deltaxuo,&
       &deltayvo,&
       &deltazwo,&
       &fexpo,&
       &feypo,&
       &fezpo,&
       &gamma_enero,&
       &u_o,&
       &v_o,&
       &w_o,&
       &temper_o,&
       &temper_anto,&
       &fuente_con_tempo,&
       &fuente_lin_tempo,&
       &dt_o,&
       &rel_energo,&
       &AI_o,AC_o,AD_o,Rx_o,&
       &ii,jj,kk&
       &)
    !
    implicit none
    !
    !$omp declare target
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
    ! Distancia entre nodos contiguos de la malla de p en direcci\'on vertical z
    !
    real(kind=DBL), dimension(lk), intent(in) :: deltazwo
    !
    ! Coeficientes de interpolaci\'on
    !
    real(kind=DBL), DIMENSION(mi), intent(in) :: fexpo
    real(kind=DBL), DIMENSION(nj), intent(in) :: feypo
    real(kind=DBL), DIMENSION(lk), intent(in) :: fezpo
    !
    ! Coeficiente de difusi\'on para la energ\'ia
    !
    real(kind=DBL), dimension(mi+1,nj+1,lk+1) :: gamma_enero
    !
    ! Velocidad, presi\'on, t\'ermino fuente b
    !
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in)  :: u_o
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in)  :: v_o
    real(kind=DBL), dimension(mi+1,nj+1,lk),   intent(in)  :: w_o
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in)  :: temper_o
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in)  :: temper_anto
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in)  :: fuente_con_tempo
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in)  :: fuente_lin_tempo
    !
    ! coeficiente de relajaci\'on
    !
    real(kind=DBL), intent(in) :: rel_energo
    !
    ! Incremento de tiempo
    !
    real(kind=DBL), intent(in) :: dt_o
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
    real(kind=DBL) :: deltax, deltay, deltaz
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
    gammad = ( gamma_enero(ii+1,jj,kk) * gamma_enero(ii,jj,kk) ) / &
         &( gamma_enero(ii+1,jj,kk)*(1._DBL-fexpo(ii))+&
         &gamma_enero(ii,jj,kk)*fexpo(ii) )
    gammai = ( gamma_enero(ii,jj,kk) * gamma_enero(ii-1,jj,kk) ) / &
         &( gamma_enero(ii,jj,kk)*(1._DBL-fexpo(ii-1))+&
         &gamma_enero(ii-1,jj,kk)*fexpo(ii-1) )
    gamman = ( gamma_enero(ii,jj+1,kk) * gamma_enero(ii,jj,kk) ) / &
         &( gamma_enero(ii,jj+1,kk) * (1._DBL-feypo(jj))+&
         &gamma_enero(ii,jj,kk)*feypo(jj) )
    gammas = ( gamma_enero(ii,jj,kk) * gamma_enero(ii,jj-1,kk) ) / &
         &( gamma_enero(ii,jj,kk) * (1._DBL-feypo(jj-1))+&
         &gamma_enero(ii,jj-1,kk)*feypo(jj-1) )
    gammat = ( gamma_enero(ii,jj,kk+1) * gamma_enero(ii,jj,kk) ) / &
         &( gamma_enero(ii,jj,kk+1) * (1._DBL-fezpo(kk))+&
         &gamma_enero(ii,jj,kk)*fezpo(kk) )
    gammab = ( gamma_enero(ii,jj,kk) * gamma_enero(ii,jj,kk-1) ) / &
         &( gamma_enero(ii,jj,kk) * (1._DBL-fezpo(kk-1))+&
         &gamma_enero(ii,jj,kk-1)*fezpo(kk-1) )
    !
    ! Tama\~no de los vol\'umenes de control para la temperatura
    !
    deltax = deltaxpo(ii)
    deltay = deltaypo(jj)
    deltaz = deltazpo(kk)
    !
    ! -------------------------
    !
    ! Coeficientes de la matriz
    !
    AI_o(indexp(ii,jj,kk)) =-(gammai*deltay*deltaz/di*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ui*di/gammai))**5)+&
         &DMAX1(0.0_DBL, ui*deltay*deltaz))
    !
    AD_o(indexp(ii,jj,kk)) =-(gammad*deltay*deltaz/dd*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ud*dd/gammad))**5)+&
         &DMAX1(0.0_DBL,-ud*deltay*deltaz))
    !
    alpha                  =-(gammas*deltax*deltaz/ds*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vs*ds/gammas))**5)+&
         &DMAX1(0.0_DBL, vs*deltax*deltaz))
    !
    beta                   =-(gamman*deltax*deltaz/dn*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vn*dn/gamman))**5)+&
         &DMAX1(0.0_DBL,-vn*deltax*deltaz))
    !
    gamma                  =-(gammab*deltax*deltay/db*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wb*db/gammab))**5)+&
         &DMAX1(0.0_DBL, wb*deltax*deltay))
    !
    delta                  =-(gammat*deltax*deltay/da*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wt*da/gammat))**5)+&
         &DMAX1(0.0_DBL,-wt*deltax*deltay))
    !
    AC_o(indexp(ii,jj,kk)) = ( -AI_o(indexp(ii,jj,kk)) - AD_o(indexp(ii,jj,kk))-&
         &alpha - beta - gamma - delta +&
         &deltax*deltay*deltaz*fuente_lin_tempo(ii,jj,kk)+&
         &deltax*deltay*deltaz/dt_o) / rel_energo
    !
    Rx_o(indexp(ii,jj,kk)) =-alpha*temper_o(ii,jj-1,kk)-&
         &beta *temper_o(ii,jj+1,kk) -&
         &gamma*temper_o(ii,jj,kk-1) -&
         &delta*temper_o(ii,jj,kk+1) +&
         &deltax*deltay*deltaz*fuente_con_tempo(ii,jj,kk)+&
         &deltax*deltay*deltaz*temper_anto(ii,jj,kk)/dt_o+&
         &AC_o(indexp(ii,jj,kk))*(1._DBL-rel_energo)*temper_o(ii,jj,kk)
    !
  end subroutine ensambla_energia_x
  !
  !*******************************************************************
  !
  ! ensambla_ec_energ_y
  !
  ! Subrutina que calcula los coeficientes de la matriz tridiagonal
  ! para la ec. de la energ\'ia en direcci\'on y
  !
  !*******************************************************************
  subroutine ensambla_energia_y(&
       &deltaxpo,&
       &deltaypo,&
       &deltazpo,&
       &deltaxuo,&
       &deltayvo,&
       &deltazwo,&
       &fexpo,&
       &feypo,&
       &fezpo,&
       &gamma_enero,&
       &u_o,&
       &v_o,&
       &w_o,&
       &temper_o,&
       &temper_anto,&
       &fuente_con_tempo,&
       &fuente_lin_tempo,&
       &dt_o,&
       &rel_energo,&
       &AI_o,AC_o,AD_o,Rx_o,&
       &jj,ii,kk&
       &)
    !
    implicit none
    !
    !$omp declare target
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
    ! Distancia entre nodos contiguos de la malla de p en direcci\'on vertical z
    !
    real(kind=DBL), dimension(lk), intent(in) :: deltazwo
    !
    ! Coeficientes de interpolaci\'on
    !
    real(kind=DBL), DIMENSION(mi), intent(in) :: fexpo
    real(kind=DBL), DIMENSION(nj), intent(in) :: feypo
    real(kind=DBL), DIMENSION(lk), intent(in) :: fezpo
    !
    ! Coeficiente de difusi\'on para la energ\'ia
    !
    real(kind=DBL), dimension(mi+1,nj+1,lk+1) :: gamma_enero
    !
    ! Velocidad, presi\'on, t\'ermino fuente b
    !
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in)  :: u_o
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in)  :: v_o
    real(kind=DBL), dimension(mi+1,nj+1,lk),   intent(in)  :: w_o
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in)  :: temper_o
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in)  :: temper_anto
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in)  :: fuente_con_tempo
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in)  :: fuente_lin_tempo
    !
    ! coeficiente de relajaci\'on
    !
    real(kind=DBL), intent(in) :: rel_energo
    !
    ! Incremento de tiempo
    !
    real(kind=DBL), intent(in) :: dt_o
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
    integer,   intent(in) :: jj, ii, kk
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
    real(kind=DBL) :: deltax, deltay, deltaz
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
    gammad = ( gamma_enero(ii+1,jj,kk) * gamma_enero(ii,jj,kk) ) / &
         &( gamma_enero(ii+1,jj,kk)*(1._DBL-fexpo(ii))+&
         &gamma_enero(ii,jj,kk)*fexpo(ii) )
    gammai = ( gamma_enero(ii,jj,kk) * gamma_enero(ii-1,jj,kk) ) / &
         &( gamma_enero(ii,jj,kk)*(1._DBL-fexpo(ii-1))+&
         &gamma_enero(ii-1,jj,kk)*fexpo(ii-1) )
    gamman = ( gamma_enero(ii,jj+1,kk) * gamma_enero(ii,jj,kk) ) / &
         &( gamma_enero(ii,jj+1,kk) * (1._DBL-feypo(jj))+&
         &gamma_enero(ii,jj,kk)*feypo(jj) )
    gammas = ( gamma_enero(ii,jj,kk) * gamma_enero(ii,jj-1,kk) ) / &
         &( gamma_enero(ii,jj,kk) * (1._DBL-feypo(jj-1))+&
         &gamma_enero(ii,jj-1,kk)*feypo(jj-1) )
    gammat = ( gamma_enero(ii,jj,kk+1) * gamma_enero(ii,jj,kk) ) / &
         &( gamma_enero(ii,jj,kk+1) * (1._DBL-fezpo(kk))+&
         &gamma_enero(ii,jj,kk)*fezpo(kk) )
    gammab = ( gamma_enero(ii,jj,kk) * gamma_enero(ii,jj,kk-1) ) / &
         &( gamma_enero(ii,jj,kk) * (1._DBL-fezpo(kk-1))+&
         &gamma_enero(ii,jj,kk-1)*fezpo(kk-1) )
    !
    ! Tama\~no de los vol\'umenes de control para la temperatura
    !
    deltax = deltaxpo(ii)
    deltay = deltaypo(jj)
    deltaz = deltazpo(kk)
    !
    ! -------------------------
    !
    ! Coeficientes de la matriz
    !
    alpha                  =-(gammai*deltay*deltaz/di*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ui*di/gammai))**5)+&
         &DMAX1(0.0_DBL, ui*deltay*deltaz))
    !
    beta                   =-(gammad*deltay*deltaz/dd*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ud*dd/gammad))**5)+&
         &DMAX1(0.0_DBL,-ud*deltay*deltaz))
    !
    AI_o(indeyp(jj,ii,kk)) =-(gammas*deltax*deltaz/ds*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vs*ds/gammas))**5)+&
         &DMAX1(0.0_DBL, vs*deltax*deltaz))
    !
    AD_o(indeyp(jj,ii,kk)) =-(gamman*deltax*deltaz/dn*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vn*dn/gamman))**5)+&
         &DMAX1(0.0_DBL,-vn*deltax*deltaz))
    !
    gamma                  =-(gammab*deltax*deltay/db*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wb*db/gammab))**5)+&
         &DMAX1(0.0_DBL, wb*deltax*deltay))
    !
    delta                  =-(gammat*deltax*deltay/da*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wt*da/gammat))**5)+&
         &DMAX1(0.0_DBL,-wt*deltax*deltay))
    !
    AC_o(indeyp(jj,ii,kk)) = ( -AI_o(indeyp(jj,ii,kk)) - AD_o(indeyp(jj,ii,kk))-&
         &alpha - beta - gamma - delta +&
         &deltax*deltay*deltaz*fuente_lin_tempo(ii,jj,kk)+&
         &deltax*deltay*deltaz/dt_o) / rel_energo
    !
    Rx_o(indeyp(jj,ii,kk)) =-alpha*temper_o(ii-1,jj,kk)-&
         &beta *temper_o(ii+1,jj,kk) -&
         &gamma*temper_o(ii,jj,kk-1) -&
         &delta*temper_o(ii,jj,kk+1) +&
         &deltax*deltay*deltaz*fuente_con_tempo(ii,jj,kk)+&
         &deltax*deltay*deltaz*temper_anto(ii,jj,kk)/dt_o+&
         &AC_o(indeyp(jj,ii,kk))*(1._DBL-rel_energo)*temper_o(ii,jj,kk)
    !
  end subroutine ensambla_energia_y
  !
  !*******************************************************************
  !
  ! ensambla_ec_energ_z
  !
  ! Subrutina que calcula los coeficientes de la matriz tridiagonal
  ! para la ec. de la energ\'ia en direcci\'on z
  !
  !*******************************************************************
  subroutine ensambla_energia_z(&
       &deltaxpo,&
       &deltaypo,&
       &deltazpo,&
       &deltaxuo,&
       &deltayvo,&
       &deltazwo,&
       &fexpo,&
       &feypo,&
       &fezpo,&
       &gamma_enero,&
       &u_o,&
       &v_o,&
       &w_o,&
       &temper_o,&
       &temper_anto,&
       &fuente_con_tempo,&
       &fuente_lin_tempo,&
       &dt_o,&
       &rel_energo,&
       &AI_o,AC_o,AD_o,Rx_o,&
       &kk,jj,ii&
       &)
    !
    implicit none
    !
    !$omp declare target
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
    ! Distancia entre nodos contiguos de la malla de p en direcci\'on vertical z
    !
    real(kind=DBL), dimension(lk), intent(in) :: deltazwo
    !
    ! Coeficientes de interpolaci\'on
    !
    real(kind=DBL), DIMENSION(mi), intent(in) :: fexpo
    real(kind=DBL), DIMENSION(nj), intent(in) :: feypo
    real(kind=DBL), DIMENSION(lk), intent(in) :: fezpo
    !
    ! Coeficiente de difusi\'on para la energ\'ia
    !
    real(kind=DBL), dimension(mi+1,nj+1,lk+1) :: gamma_enero
    !
    ! Velocidad, presi\'on, t\'ermino fuente b
    !
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in)  :: u_o
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in)  :: v_o
    real(kind=DBL), dimension(mi+1,nj+1,lk),   intent(in)  :: w_o
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in)  :: temper_o
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in)  :: temper_anto
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in)  :: fuente_con_tempo
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in)  :: fuente_lin_tempo
    !
    ! coeficiente de relajaci\'on
    !
    real(kind=DBL), intent(in) :: rel_energo
    !
    ! Incremento de tiempo
    !
    real(kind=DBL), intent(in) :: dt_o
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
    integer,   intent(in) :: jj, ii, kk
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
    real(kind=DBL) :: deltax, deltay, deltaz
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
    gammad = ( gamma_enero(ii+1,jj,kk) * gamma_enero(ii,jj,kk) ) / &
         &( gamma_enero(ii+1,jj,kk)*(1._DBL-fexpo(ii))+&
         &gamma_enero(ii,jj,kk)*fexpo(ii) )
    gammai = ( gamma_enero(ii,jj,kk) * gamma_enero(ii-1,jj,kk) ) / &
         &( gamma_enero(ii,jj,kk)*(1._DBL-fexpo(ii-1))+&
         &gamma_enero(ii-1,jj,kk)*fexpo(ii-1) )
    gamman = ( gamma_enero(ii,jj+1,kk) * gamma_enero(ii,jj,kk) ) / &
         &( gamma_enero(ii,jj+1,kk) * (1._DBL-feypo(jj))+&
         &gamma_enero(ii,jj,kk)*feypo(jj) )
    gammas = ( gamma_enero(ii,jj,kk) * gamma_enero(ii,jj-1,kk) ) / &
         &( gamma_enero(ii,jj,kk) * (1._DBL-feypo(jj-1))+&
         &gamma_enero(ii,jj-1,kk)*feypo(jj-1) )
    gammat = ( gamma_enero(ii,jj,kk+1) * gamma_enero(ii,jj,kk) ) / &
         &( gamma_enero(ii,jj,kk+1) * (1._DBL-fezpo(kk))+&
         &gamma_enero(ii,jj,kk)*fezpo(kk) )
    gammab = ( gamma_enero(ii,jj,kk) * gamma_enero(ii,jj,kk-1) ) / &
         &( gamma_enero(ii,jj,kk) * (1._DBL-fezpo(kk-1))+&
         &gamma_enero(ii,jj,kk-1)*fezpo(kk-1) )
    !
    ! Tama\~no de los vol\'umenes de control para la temperatura
    !
    deltax = deltaxpo(ii)
    deltay = deltaypo(jj)
    deltaz = deltazpo(kk)
    !
    ! -------------------------
    !
    ! Coeficientes de la matriz
    !
    alpha                  =-(gammai*deltay*deltaz/di*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ui*di/gammai))**5)+&
         &DMAX1(0.0_DBL, ui*deltay*deltaz))
    !
    beta                   =-(gammad*deltay*deltaz/dd*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ud*dd/gammad))**5)+&
         &DMAX1(0.0_DBL,-ud*deltay*deltaz))
    !
    gamma                   =-(gammas*deltax*deltaz/ds*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vs*ds/gammas))**5)+&
         &DMAX1(0.0_DBL, vs*deltax*deltaz))
    !
    delta                   =-(gamman*deltax*deltaz/dn*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vn*dn/gamman))**5)+&
         &DMAX1(0.0_DBL,-vn*deltax*deltaz))
    !
    AI_o(indezp(kk,jj,ii)) =-(gammab*deltax*deltay/db*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wb*db/gammab))**5)+&
         &DMAX1(0.0_DBL, wb*deltax*deltay))
    !
    AD_o(indezp(kk,jj,ii)) =-(gammat*deltax*deltay/da*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wt*da/gammat))**5)+&
         &DMAX1(0.0_DBL,-wt*deltax*deltay))
    !
    AC_o(indezp(kk,jj,ii)) = ( -AI_o(indezp(kk,jj,ii)) - AD_o(indezp(kk,jj,ii))-&
         &alpha - beta - gamma - delta +&
         &deltax*deltay*deltaz*fuente_lin_tempo(ii,jj,kk)+&
         &deltax*deltay*deltaz/dt_o) / rel_energo
    !
    Rx_o(indezp(kk,jj,ii)) =-alpha*temper_o(ii-1,jj,kk)-&
         &beta *temper_o(ii+1,jj,kk) -&
         &gamma*temper_o(ii,jj-1,kk) -&
         &delta*temper_o(ii,jj+1,kk) +&
         &deltax*deltay*deltaz*fuente_con_tempo(ii,jj,kk)+&
         &deltax*deltay*deltaz*temper_anto(ii,jj,kk)/dt_o+&
         &AC_o(indezp(kk,jj,ii))*(1._DBL-rel_energo)*temper_o(ii,jj,kk)
    !
  end subroutine ensambla_energia_z
  !
end module ec_energia
