!
!
! M\'odulo para la ecuaci\'on de momento
!
! Este m\'odulo contiene:
!
! - Las variables de la ecuaci\'on de momento
! - Las subrutinas que calculan los coeficientes de las ecuaciones discretizadas
! - Las subrutinas que sirven para leer e imponer condiciones de frontera
!
module ec_momento
  !
  use malla, only : mi, nj, lk, DBL
  use malla, only : indexu, indeyu, indezp  ! indices para u
  use malla, only : indexv, indeyv, indezv  ! indices para v
  use malla, only : indexp, indeyp, indezw  ! indices para w
  !
  ! use cond_frontera, only : tipo_cond_front
  !
  ! use cond_frontera, only : inicializa_cond_front
  ! use cond_frontera, only : lectura_cond_frontera
  !
  implicit none
  !
  ! Componentes de velocidad, presi\'on
  ! residuos de las ecuaciones de momento y correcci\'on de la presi\'on
  ! coeficientes de difusi\'on y criterios e convergencia
  !
  real(kind=DBL), dimension(mi,nj+1,lk+1)   :: u, u_ant, du, au, Resu, fu
  real(kind=DBL), dimension(mi+1,nj,lk+1)   :: v, v_ant, dv, av, Resv, fv
  real(kind=DBL), dimension(mi+1,nj+1,lk)   :: w, w_ant, dw, aw, Resw, fw
  real(kind=DBL), dimension(mi+1,nj+1,lk+1) :: gamma_momen, Ri
  real(kind=DBL) :: conv_u, rel_vel
  !
  ! Variables para los t\'erminos fuente de la ec de momento
  ! El t\'ermino lineal fuente_lin debe ser negativo para favorecer
  ! convergencia y estabilidad
  !
  real(kind=DBL), dimension(mi,nj+1,lk+1)   :: fuente_con_u, fuente_lin_u
  real(kind=DBL), dimension(mi+1,nj,lk+1)   :: fuente_con_v, fuente_lin_v
  real(kind=DBL), dimension(mi+1,nj+1,lk)   :: fuente_con_w, fuente_lin_w
  !
  ! Variables de postproceso
  !
  real(kind=DBL), dimension(mi+1,nj+1,lk+1) :: uf, vf, wf
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
  ! condicion_inicial_uv
  !
  ! Subrutina que inicializa los arreglos para u, v y p de acuerdo a
  ! las conndiciones iniciales
  !
  !*******************************************************************  
  subroutine condicion_inicial_flujo(cond_inicial)
    !
    implicit none
    ! $acc routine
    !
    integer :: ii, jj, kk
    character(len=8), intent(in) :: cond_inicial
    !
    flujo_inicial: select case( trim(cond_inicial) )
       !
    case('nulo')
       !
       !$acc loop gang
       do kk = 1, lk+1
          do jj = 1, nj+1
             do ii = 1, mi
                u_ant(ii,jj,kk) = 0.0_DBL
             end do
          end do
       end do
       !$acc loop gang
       do kk = 1, lk+1
          do jj = 1, nj
             do ii = 1, mi+1
                v_ant(ii,jj,kk) = 0.0_DBL
             end do
          end do
       end do
       !$acc loop gang
       do kk = 1, lk
          do jj = 1, nj+1
             do ii = 1, mi+1
                w_ant(ii,jj,kk) = 0.0_DBL
             end do
          end do
       end do       
       !
    case('archivo')
       !
       return
       !
    end select flujo_inicial
  end subroutine condicion_inicial_flujo
  !
  !*******************************************************************
  !
  ! ini_frontera_uv
  !
  ! Subrutina que inicializa los arreglos para las condiciones
  ! de frontera de u y de v, y lee las condiciones de los archivos
  ! de entrada cond_fronterau.dat y cond_fronterav.dat
  !
  !*******************************************************************  
  ! subroutine ini_frontera_uv()
  !   !
  !   implicit none
  !   !
  !   ! arreglos de u
  !   !
  !   call inicializa_cond_front(cond_front_ua)
  !   call inicializa_cond_front(cond_front_ub)
  !   call inicializa_cond_front(cond_front_uc)
  !   call inicializa_cond_front(cond_front_ud)
  !   call lectura_cond_frontera('cond_fronterau.dat',&
  !        & cond_front_ua, &
  !        & cond_front_ub, &
  !        & cond_front_uc, &
  !        & cond_front_ud, &
  !        & xu, yp,        &
  !        & mi, nj+1       &
  !        & )
  !   !
  !   ! arreglos de v
  !   !
  !   call inicializa_cond_front(cond_front_va)
  !   call inicializa_cond_front(cond_front_vb)
  !   call inicializa_cond_front(cond_front_vc)
  !   call inicializa_cond_front(cond_front_vd)
  !   !
  !   call lectura_cond_frontera('cond_fronterav.dat',&
  !        & cond_front_va, &
  !        & cond_front_vb, &
  !        & cond_front_vc, &
  !        & cond_front_vd, &
  !        & xp, yv,        &
  !        & mi+1, nj       &
  !        & )
  !   !
  ! end subroutine ini_frontera_uv
  !
  !-------------------------------------------------------------------  
  !*******************************************************************
  !
  ! ensambla_velw_z
  !
  ! Subrutina que calcula los coeficientes de la matriz tridiagonal
  ! para la velocidad w en la direcci\'on z
  !
  !*******************************************************************
  !-------------------------------------------------------------------
  subroutine ensambla_velw_z(&
       &deltaxwo,&
       &deltaywo,&
       &deltazwo,&
       &deltaxuo,&
       &deltayvo,&
       &deltazpo,&
       &fexpo,&
       &feypo,&
       &fezpo,&
       &fezwo,&
       &gamma_momento,&
       &u_o,&
       &v_o,&
       &w_o,&
       &w_anto,&
       &temp_o,&
       &pres_o,&
       &fuente_con_wo,&
       &fuente_lin_wo,&
       &Ri_o,&
       &dt_o,&
       &rel_vo,&
       &AI_o,AC_o,AD_o,Rx_o,&
       &kk,jj,ii&
       &)
    implicit none
    !$acc routine
    !
    ! Tama\~no del volumen de control
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxwo
    real(kind=DBL), dimension(nj), intent(in) :: deltaywo
    real(kind=DBL), dimension(lk), intent(in) :: deltazwo
    !
    ! Distancia entre nodos contiguos de la malla de w en direcci\'on horizontal x
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxuo
    !
    ! Distancia entre nodos contiguos de la malla de w en direcci\'on horizontal y
    !
    real(kind=DBL), dimension(nj), intent(in) :: deltayvo
    !
    ! Distancia entre nodos contiguos de la malla de w en direcci\'on vertical z
    !
    real(kind=DBL), dimension(lk), intent(in) :: deltazpo
    !
    ! Coeficientes para interpolaci\'on
    !
    real(kind=DBL), DIMENSION(mi),   intent(in) :: fexpo
    real(kind=DBL), DIMENSION(nj),   intent(in) :: feypo
    real(kind=DBL), dimension(lk),   intent(in) :: fezpo
    real(kind=DBL), DIMENSION(lk-1), intent(in) :: fezwo
    !
    ! Coeficiente de difusi\'on
    !
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in) :: gamma_momento
    !
    ! Velocidad, presi\'on y temperatura
    !
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in) :: u_o
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in) :: v_o
    real(kind=DBL), dimension(mi+1,nj+1,lk),   intent(in) :: w_o, w_anto
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in) :: temp_o, pres_o
    !
    ! T\'erminos fuente
    !
    real(kind=DBL), dimension(mi+1,nj+1,lk),   intent(in) :: fuente_con_wo
    real(kind=DBL), dimension(mi+1,nj+1,lk),   intent(in) :: fuente_lin_wo    
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in) :: Ri_o
    
    !
    ! Incremento de tiempo y coeficiente de relajaci\'on
    !
    real(kind=DBL), intent(in) :: dt_o, rel_vo
    !
    ! \'Indices para recorrer las direcciones x, y, z
    !
    integer, intent(in)        :: kk, jj, ii
    !
    ! Coeficientes de las matrices
    !
    ! ** Estos coeficientes est\'an sobredimensionados para reducir el uso de memoria
    ! en la gpu, los arreglos que se reciben en esta subrutina se usan para las ecs.
    ! de momento, energ\'ia y la correcci\'on de la presi\'on **
    !
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AI_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AC_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AD_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: RX_o
    !
    ! Variables auxiliares
    !
    integer :: info
    !
    ! Auxiliares de interpolaci\'on y coeficientes
    !
    real(kind=DBL) :: ui, ud, vs, vn, wb, wt
    real(kind=DBL) :: di, dd, ds, dn, db, dt
    real(kind=DBL) :: gammai, gammad
    real(kind=DBL) :: gammas, gamman
    real(kind=DBL) :: gammab, gammat
    real(kind=DBL) :: alpha, beta, gamma, delta
    real(kind=DBL) :: deltax, deltay, deltaz
    real(kind=DBL) :: temp_int
    !
    ! Interpolaciones necesarias
    !
    ! u
    !
    ud = fezpo(kk)*u_o(ii,jj,kk+1)  +(1.0_DBL-fezpo(kk))*u_o(ii,jj,kk)
    ui = fezpo(kk)*u_o(ii-1,jj,kk+1)+(1.0_DBL-fezpo(kk))*u_o(ii-1,jj,kk)
    !
    ! v
    !
    vn = fezpo(kk)*v_o(ii,jj,kk+1)  +(1.0_DBL-fezpo(kk))*v_o(ii,jj,kk)
    vs = fezpo(kk)*v_o(ii,jj-1,kk+1)+(1.0_DBL-fezpo(kk))*v_o(ii,jj-1,kk)
    !
    ! w
    !
    wt = fezwo(kk)  *w_o(ii,jj,kk+1)+(1.0_DBL-fezwo(kk))  *w_o(ii,jj,kk)
    wb = fezwo(kk-1)*w_o(ii,jj,kk)  +(1.0_DBL-fezwo(kk-1))*w_o(ii,jj,kk-1)
    !
    ! gamma_n 
    !
    ! ** se utilizan las constantes gammai y gammad como auxiliares para
    ! calcular gamman, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammad = ( gamma_momento(ii,jj+1,kk+1) * gamma_momento(ii,jj,kk+1) ) / &
         &(gamma_momento(ii,jj+1,kk+1)*(1._DBL-feypo(jj))+&
         &gamma_momento(ii,jj,kk+1)*feypo(jj) )
    gammai = ( gamma_momento(ii,jj+1,kk) * gamma_momento(ii,jj,kk) ) / &
         &( gamma_momento(ii,jj+1,kk) * (1._DBL-feypo(jj))+&
         &gamma_momento(ii,jj,kk)*feypo(jj) )
    !
    gamman = gammai*gammad / (gammad * (1._DBL-fezpo(kk)) + gammai * fezpo(kk))
    !
    ! gamma_s 
    !
    ! ** se utilizan las constantes gammai y gammad como auxiliares para
    ! calcular gamman, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammad = ( gamma_momento(ii,jj,kk+1) * gamma_momento(ii,jj-1,kk+1) ) / &
         &(gamma_momento(ii,jj,kk+1)*(1._DBL-feypo(jj-1))+&
         &gamma_momento(ii,jj-1,kk+1)*feypo(jj-1))
    gammai = ( gamma_momento(ii,jj,kk) * gamma_momento(ii,jj-1,kk) ) / &
         &(gamma_momento(ii,jj,kk)*(1._DBL-feypo(jj-1))+&
         &gamma_momento(ii,jj-1,kk)*feypo(jj-1))
    !
    gammas = gammai*gammad / (gammad * (1._DBL-fezpo(kk)) + gammai * fezpo(kk))
    !
    ! gamma_d 
    !
    ! ** se utilizan las constantes gammat y gammab como auxiliares para
    ! calcular gammad, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammat = ( gamma_momento(ii+1,jj,kk+1) * gamma_momento(ii,jj,kk+1) ) / &
         &(gamma_momento(ii+1,jj,kk+1)*(1._DBL-fexpo(ii))+&
         &gamma_momento(ii,jj,kk+1)*fexpo(ii))
    gammab = ( gamma_momento(ii+1,jj,kk) * gamma_momento(ii,jj,kk) ) / &
         &(gamma_momento(ii+1,jj,kk)*(1._DBL-fexpo(ii))+&
         &gamma_momento(ii,jj,kk)*fexpo(ii))
    !
    gammad = gammat*gammab / (gammat * (1._DBL-fezpo(kk)) + gammab * fezpo(kk))
    !
    ! gamma_i
    !
    ! ** se utilizan las constantes gammat y gammab como auxiliares para
    ! calcular gammai, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammat = ( gamma_momento(ii,jj,kk+1) * gamma_momento(ii-1,jj,kk+1) ) / &
         &( gamma_momento(ii,jj,kk+1)*(1._DBL-fexpo(ii-1))+&
         &gamma_momento(ii-1,jj,kk+1)*fexpo(ii-1) )
    gammab = ( gamma_momento(ii,jj,kk) * gamma_momento(ii-1,jj,kk) ) / &
         &( gamma_momento(ii,jj,kk) * (1._DBL-fexpo(ii-1))+&
         &gamma_momento(ii-1,jj,kk)*fexpo(ii-1) )
    !
    gammai = gammat*gammab / (gammat * (1._DBL-fezpo(kk)) + gammab * fezpo(kk))
    !
    ! gamma_t
    !
    gammat = gamma_momento(ii,jj,kk+1)
    !
    ! gamma_b
    !
    gammab = gamma_momento(ii,jj,kk)
    !
    ! distancias entre nodos contiguos
    !
    di = deltaxuo(ii-1)
    dd = deltaxuo(ii)
    ds = deltayvo(jj-1)
    dn = deltayvo(jj)
    db = deltazpo(kk)
    dt = deltazpo(kk+1)
    !
    ! Tama\~no de los vol\'umenes de control para la velocidad u
    !
    deltax = deltaxwo(ii)
    deltay = deltaywo(jj)
    deltaz = deltazwo(kk)
    !
    ! Interpolaci\'on para la temperatura
    !
    temp_int = fezpo(kk)*temp_o(ii,jj,kk+1) + (1.0_DBL-fezpo(kk))*temp_o(ii,jj,kk)
    !
    ! *************************
    !
    ! Coeficientes de la matriz
    !
    alpha =-(gammai*deltay*deltaz/di*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ui*di/gammai))**5)+&
         &DMAX1(0.0_DBL, ui*deltay*deltaz))
    !
    beta  =-(gammad*deltay*deltaz/dd*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ud*dd/gammad))**5)+&
         &DMAX1(0.0_DBL,-ud*deltay*deltaz))
    !
    gamma =-(gammas*deltax*deltaz/ds*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vs*ds/gammas))**5)+&
         &DMAX1(0.0_DBL, vs*deltax*deltaz))
    !
    delta =-(gamman*deltax*deltaz/dn*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vn*dn/gamman))**5)+&
         &DMAX1(0.0_DBL,-vn*deltax*deltaz))
    !
    AI_o(indezw(kk,jj,ii)) =-(gammab*deltax*deltay/db*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wb*db/gammab))**5)+&
         &DMAX1(0.0_DBL, wb*deltax*deltay))
    !
    AD_o(indezw(kk,jj,ii)) =-(gammat*deltax*deltay/dt*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wt*dt/gammat))**5)+&
         &DMAX1(0.0_DBL,-wt*deltax*deltay)) 
    !
    AC_o(indezw(kk,jj,ii)) = ( -AI_o(indezw(kk,jj,ii)) - AD_o(indezw(kk,jj,ii)) &
         &- alpha - beta - gamma - delta - &
         &deltax*deltay*deltaz*fuente_lin_wo(ii,jj,kk)+&
         &deltax*deltay*deltaz/dt_o ) / rel_vo
    !
    Rx_o(indezw(kk,jj,ii)) =-alpha*w_o(ii-1,jj,kk) -&
         &beta  * w_o(ii+1,jj,kk) - &
         &gamma * w_o(ii,jj-1,kk) - &
         &delta * w_o(ii,jj+1,kk) - &
         &deltax*deltay*deltaz*Ri_o(ii,jj,kk)*temp_int+&
         &deltax*deltay*deltaz*fuente_con_wo(ii,jj,kk)+&
         &deltax*deltay*deltaz*w_anto(ii,jj,kk)/dt_o+&
         &(pres_o(ii,jj,kk)-pres_o(ii,jj,kk+1))*deltax*deltay+&
         &AC_o(indezw(kk,jj,ii))*(1._DBL-rel_vo)*w_o(ii,jj,kk)
    !
  end subroutine ensambla_velw_z
  !
  !-------------------------------------------------------------------  
  !*******************************************************************
  !
  ! ensambla_velw_y
  !
  ! Subrutina que calcula los coeficientes de la matriz tridiagonal
  ! para la velocidad w en la direcci\'on y
  !
  !*******************************************************************
  !-------------------------------------------------------------------
  subroutine ensambla_velw_y(&
       &deltaxwo,&
       &deltaywo,&
       &deltazwo,&
       &deltaxuo,&
       &deltayvo,&
       &deltazpo,&
       &fexpo,&
       &feypo,&
       &fezpo,&
       &fezwo,&
       &gamma_momento,&
       &u_o,&
       &v_o,&
       &w_o,&
       &w_anto,&
       &temp_o,&
       &pres_o,&
       &fuente_con_wo,&
       &fuente_lin_wo,&
       &Ri_o,&
       &dt_o,&
       &rel_vo,&
       &AI_o,AC_o,AD_o,Rx_o,&
       &jj,ii,kk&
       &)
    implicit none
    !$acc routine
    !
    ! Tama\~no del volumen de control
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxwo
    real(kind=DBL), dimension(nj), intent(in) :: deltaywo
    real(kind=DBL), dimension(lk), intent(in) :: deltazwo
    !
    ! Distancia entre nodos contiguos de la malla de w en direcci\'on horizontal x
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxuo
    !
    ! Distancia entre nodos contiguos de la malla de w en direcci\'on horizontal y
    !
    real(kind=DBL), dimension(nj), intent(in) :: deltayvo
    !
    ! Distancia entre nodos contiguos de la malla de w en direcci\'on vertical z
    !
    real(kind=DBL), dimension(lk), intent(in) :: deltazpo
    !
    ! Coeficientes para interpolaci\'on
    !
    real(kind=DBL), DIMENSION(mi),   intent(in) :: fexpo
    real(kind=DBL), DIMENSION(nj),   intent(in) :: feypo
    real(kind=DBL), dimension(lk),   intent(in) :: fezpo
    real(kind=DBL), DIMENSION(lk-1), intent(in) :: fezwo
    !
    ! Coeficiente de difusi\'on
    !
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in) :: gamma_momento
    !
    ! Velocidad, presi\'on y temperatura
    !
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in) :: u_o
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in) :: v_o
    real(kind=DBL), dimension(mi+1,nj+1,lk),   intent(in) :: w_o, w_anto
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in) :: temp_o, pres_o
    !
    ! T\'erminos fuente
    !
    real(kind=DBL), dimension(mi+1,nj+1,lk),   intent(in) :: fuente_con_wo
    real(kind=DBL), dimension(mi+1,nj+1,lk),   intent(in) :: fuente_lin_wo    
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in) :: Ri_o
    
    !
    ! Incremento de tiempo y coeficiente de relajaci\'on
    !
    real(kind=DBL), intent(in) :: dt_o, rel_vo
    !
    ! \'Indices para recorrer las direcciones x, y, z
    !
    integer, intent(in)        :: jj, ii, kk
    !
    ! Coeficientes de las matrices
    !
    ! ** Estos coeficientes est\'an sobredimensionados para reducir el uso de memoria
    ! en la gpu, los arreglos que se reciben en esta subrutina se usan para las ecs.
    ! de momento, energ\'ia y la correcci\'on de la presi\'on **
    !
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AI_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AC_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AD_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: RX_o
    !
    ! Variables auxiliares
    !
    integer :: info
    !
    ! Auxiliares de interpolaci\'on y coeficientes
    !
    real(kind=DBL) :: ui, ud, vs, vn, wb, wt
    real(kind=DBL) :: di, dd, ds, dn, db, dt
    real(kind=DBL) :: gammai, gammad
    real(kind=DBL) :: gammas, gamman
    real(kind=DBL) :: gammab, gammat
    real(kind=DBL) :: alpha, beta, gamma, delta
    real(kind=DBL) :: deltax, deltay, deltaz
    real(kind=DBL) :: temp_int
    !
    ! Interpolaciones necesarias
    !
    ! u
    !
    ud = fezpo(kk)*u_o(ii,jj,kk+1)  +(1.0_DBL-fezpo(kk))*u_o(ii,jj,kk)
    ui = fezpo(kk)*u_o(ii-1,jj,kk+1)+(1.0_DBL-fezpo(kk))*u_o(ii-1,jj,kk)
    !
    ! v
    !
    vn = fezpo(kk)*v_o(ii,jj,kk+1)  +(1.0_DBL-fezpo(kk))*v_o(ii,jj,kk)
    vs = fezpo(kk)*v_o(ii,jj-1,kk+1)+(1.0_DBL-fezpo(kk))*v_o(ii,jj-1,kk)
    !
    ! w
    !
    wt = fezwo(kk)  *w_o(ii,jj,kk+1)+(1.0_DBL-fezwo(kk))  *w_o(ii,jj,kk)
    wb = fezwo(kk-1)*w_o(ii,jj,kk)  +(1.0_DBL-fezwo(kk-1))*w_o(ii,jj,kk-1)
    !
    ! gamma_n 
    !
    ! ** se utilizan las constantes gammai y gammad como auxiliares para
    ! calcular gamman, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammad = ( gamma_momento(ii,jj+1,kk+1) * gamma_momento(ii,jj,kk+1) ) / &
         &(gamma_momento(ii,jj+1,kk+1)*(1._DBL-feypo(jj))+&
         &gamma_momento(ii,jj,kk+1)*feypo(jj) )
    gammai = ( gamma_momento(ii,jj+1,kk) * gamma_momento(ii,jj,kk) ) / &
         &( gamma_momento(ii,jj+1,kk) * (1._DBL-feypo(jj))+&
         &gamma_momento(ii,jj,kk)*feypo(jj) )
    !
    gamman = gammai*gammad / (gammad * (1._DBL-fezpo(kk)) + gammai * fezpo(kk))
    !
    ! gamma_s 
    !
    ! ** se utilizan las constantes gammai y gammad como auxiliares para
    ! calcular gamman, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammad = ( gamma_momento(ii,jj,kk+1) * gamma_momento(ii,jj-1,kk+1) ) / &
         &(gamma_momento(ii,jj,kk+1)*(1._DBL-feypo(jj-1))+&
         &gamma_momento(ii,jj-1,kk+1)*feypo(jj-1))
    gammai = ( gamma_momento(ii,jj,kk) * gamma_momento(ii,jj-1,kk) ) / &
         &(gamma_momento(ii,jj,kk)*(1._DBL-feypo(jj-1))+&
         &gamma_momento(ii,jj-1,kk)*feypo(jj-1))
    !
    gammas = gammai*gammad / (gammad * (1._DBL-fezpo(kk)) + gammai * fezpo(kk))
    !
    ! gamma_d 
    !
    ! ** se utilizan las constantes gammat y gammab como auxiliares para
    ! calcular gammad, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammat = ( gamma_momento(ii+1,jj,kk+1) * gamma_momento(ii,jj,kk+1) ) / &
         &(gamma_momento(ii+1,jj,kk+1)*(1._DBL-fexpo(ii))+&
         &gamma_momento(ii,jj,kk+1)*fexpo(ii))
    gammab = ( gamma_momento(ii+1,jj,kk) * gamma_momento(ii,jj,kk) ) / &
         &(gamma_momento(ii+1,jj,kk)*(1._DBL-fexpo(ii))+&
         &gamma_momento(ii,jj,kk)*fexpo(ii))
    !
    gammad = gammat*gammab / (gammat * (1._DBL-fezpo(kk)) + gammab * fezpo(kk))
    !
    ! gamma_i
    !
    ! ** se utilizan las constantes gammat y gammab como auxiliares para
    ! calcular gammai, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammat = ( gamma_momento(ii,jj,kk+1) * gamma_momento(ii-1,jj,kk+1) ) / &
         &( gamma_momento(ii,jj,kk+1)*(1._DBL-fexpo(ii-1))+&
         &gamma_momento(ii-1,jj,kk+1)*fexpo(ii-1) )
    gammab = ( gamma_momento(ii,jj,kk) * gamma_momento(ii-1,jj,kk) ) / &
         &( gamma_momento(ii,jj,kk) * (1._DBL-fexpo(ii-1))+&
         &gamma_momento(ii-1,jj,kk)*fexpo(ii-1) )
    !
    gammai = gammat*gammab / (gammat * (1._DBL-fezpo(kk)) + gammab * fezpo(kk))
    !
    ! gamma_t
    !
    gammat = gamma_momento(ii,jj,kk+1)
    !
    ! gamma_b
    !
    gammab = gamma_momento(ii,jj,kk)
    !
    ! distancias entre nodos contiguos
    !
    di = deltaxuo(ii-1)
    dd = deltaxuo(ii)
    ds = deltayvo(jj-1)
    dn = deltayvo(jj)
    db = deltazpo(kk)
    dt = deltazpo(kk+1)
    !
    ! Tama\~no de los vol\'umenes de control para la velocidad u
    !
    deltax = deltaxwo(ii)
    deltay = deltaywo(jj)
    deltaz = deltazwo(kk)
    !
    ! Interpolaci\'on para la temperatura
    !
    temp_int = fezpo(kk)*temp_o(ii,jj,kk+1) + (1.0_DBL-fezpo(kk))*temp_o(ii,jj,kk)
    !
    ! *************************
    !
    ! Coeficientes de la matriz
    !
    alpha =-(gammai*deltay*deltaz/di*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ui*di/gammai))**5)+&
         &DMAX1(0.0_DBL, ui*deltay*deltaz))
    !
    beta  =-(gammad*deltay*deltaz/dd*&
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
    gamma =-(gammab*deltax*deltay/db*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wb*db/gammab))**5)+&
         &DMAX1(0.0_DBL, wb*deltax*deltay))
    !
    delta =-(gammat*deltax*deltay/dt*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wt*dt/gammat))**5)+&
         &DMAX1(0.0_DBL,-wt*deltax*deltay)) 
    !
    AC_o(indeyp(jj,ii,kk)) = ( -AI_o(indeyp(jj,ii,kk)) - AD_o(indeyp(jj,ii,kk)) &
         &- alpha - beta - gamma - delta - &
         &deltax*deltay*deltaz*fuente_lin_wo(ii,jj,kk)+&
         &deltax*deltay*deltaz/dt_o ) / rel_vo
    !
    Rx_o(indeyp(jj,ii,kk)) =-alpha*w_o(ii-1,jj,kk) -&
         &beta  * w_o(ii+1,jj,kk) - &
         &gamma * w_o(ii,jj,kk-1) - &
         &delta * w_o(ii,jj,kk+1) - &
         &deltax*deltay*deltaz*Ri_o(ii,jj,kk)*temp_int+&
         &deltax*deltay*deltaz*fuente_con_wo(ii,jj,kk)+&
         &deltax*deltay*deltaz*w_anto(ii,jj,kk)/dt_o+&
         &(pres_o(ii,jj,kk)-pres_o(ii,jj,kk+1))*deltax*deltay+&
         &AC_o(indeyp(jj,ii,kk))*(1._DBL-rel_vo)*w_o(ii,jj,kk)
    !
  end subroutine ensambla_velw_y
  !
  !-------------------------------------------------------------------  
  !*******************************************************************
  !
  ! ensambla_velw_x
  !
  ! Subrutina que calcula los coeficientes de la matriz tridiagonal
  ! para la velocidad w en la direcci\'on x
  !
  !*******************************************************************
  !-------------------------------------------------------------------
  subroutine ensambla_velw_x(&
       &deltaxwo,&
       &deltaywo,&
       &deltazwo,&
       &deltaxuo,&
       &deltayvo,&
       &deltazpo,&
       &fexpo,&
       &feypo,&
       &fezpo,&
       &fezwo,&
       &gamma_momento,&
       &u_o,&
       &v_o,&
       &w_o,&
       &w_anto,&
       &temp_o,&
       &pres_o,&
       &fuente_con_wo,&
       &fuente_lin_wo,&
       &Ri_o,&
       &dt_o,&
       &rel_vo,&
       &AI_o,AC_o,AD_o,Rx_o,&
       &aw_o,&
       &ii,jj,kk&
       &)
    implicit none
    !$acc routine
    !
    ! Tama\~no del volumen de control
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxwo
    real(kind=DBL), dimension(nj), intent(in) :: deltaywo
    real(kind=DBL), dimension(lk), intent(in) :: deltazwo
    !
    ! Distancia entre nodos contiguos de la malla de w en direcci\'on horizontal x
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxuo
    !
    ! Distancia entre nodos contiguos de la malla de w en direcci\'on horizontal y
    !
    real(kind=DBL), dimension(nj), intent(in) :: deltayvo
    !
    ! Distancia entre nodos contiguos de la malla de w en direcci\'on vertical z
    !
    real(kind=DBL), dimension(lk), intent(in) :: deltazpo
    !
    ! Coeficientes para interpolaci\'on
    !
    real(kind=DBL), DIMENSION(mi),   intent(in) :: fexpo
    real(kind=DBL), DIMENSION(nj),   intent(in) :: feypo
    real(kind=DBL), dimension(lk),   intent(in) :: fezpo
    real(kind=DBL), DIMENSION(lk-1), intent(in) :: fezwo
    !
    ! Coeficiente de difusi\'on
    !
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in) :: gamma_momento
    !
    ! Velocidad, presi\'on y temperatura
    !
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in) :: u_o
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in) :: v_o
    real(kind=DBL), dimension(mi+1,nj+1,lk),   intent(in) :: w_o, w_anto
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in) :: temp_o, pres_o
    !
    ! T\'erminos fuente
    !
    real(kind=DBL), dimension(mi+1,nj+1,lk),   intent(in) :: fuente_con_wo
    real(kind=DBL), dimension(mi+1,nj+1,lk),   intent(in) :: fuente_lin_wo    
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in) :: Ri_o
    
    !
    ! Incremento de tiempo y coeficiente de relajaci\'on
    !
    real(kind=DBL), intent(in) :: dt_o, rel_vo
    !
    ! \'Indices para recorrer las direcciones x, y, z
    !
    integer, intent(in)        :: ii, jj, kk
    !
    ! Coeficientes de las matrices
    !
    ! ** Estos coeficientes est\'an sobredimensionados para reducir el uso de memoria
    ! en la gpu, los arreglos que se reciben en esta subrutina se usan para las ecs.
    ! de momento, energ\'ia y la correcci\'on de la presi\'on **
    !
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AI_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AC_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AD_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: RX_o
    !
    real(kind=DBL), dimension(mi+1,nj+1,lk),         intent(out) :: aw_o
    !
    ! Variables auxiliares
    !
    integer :: info
    !
    ! Auxiliares de interpolaci\'on y coeficientes
    !
    real(kind=DBL) :: ui, ud, vs, vn, wb, wt
    real(kind=DBL) :: di, dd, ds, dn, db, dt
    real(kind=DBL) :: gammai, gammad
    real(kind=DBL) :: gammas, gamman
    real(kind=DBL) :: gammab, gammat
    real(kind=DBL) :: alpha, beta, gamma, delta
    real(kind=DBL) :: deltax, deltay, deltaz
    real(kind=DBL) :: temp_int
    !
    ! Interpolaciones necesarias
    !
    ! u
    !
    ud = fezpo(kk)*u_o(ii,jj,kk+1)  +(1.0_DBL-fezpo(kk))*u_o(ii,jj,kk)
    ui = fezpo(kk)*u_o(ii-1,jj,kk+1)+(1.0_DBL-fezpo(kk))*u_o(ii-1,jj,kk)
    !
    ! v
    !
    vn = fezpo(kk)*v_o(ii,jj,kk+1)  +(1.0_DBL-fezpo(kk))*v_o(ii,jj,kk)
    vs = fezpo(kk)*v_o(ii,jj-1,kk+1)+(1.0_DBL-fezpo(kk))*v_o(ii,jj-1,kk)
    !
    ! w
    !
    wt = fezwo(kk)  *w_o(ii,jj,kk+1)+(1.0_DBL-fezwo(kk))  *w_o(ii,jj,kk)
    wb = fezwo(kk-1)*w_o(ii,jj,kk)  +(1.0_DBL-fezwo(kk-1))*w_o(ii,jj,kk-1)
    !
    ! gamma_n 
    !
    ! ** se utilizan las constantes gammai y gammad como auxiliares para
    ! calcular gamman, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammad = ( gamma_momento(ii,jj+1,kk+1) * gamma_momento(ii,jj,kk+1) ) / &
         &(gamma_momento(ii,jj+1,kk+1)*(1._DBL-feypo(jj))+&
         &gamma_momento(ii,jj,kk+1)*feypo(jj) )
    gammai = ( gamma_momento(ii,jj+1,kk) * gamma_momento(ii,jj,kk) ) / &
         &( gamma_momento(ii,jj+1,kk) * (1._DBL-feypo(jj))+&
         &gamma_momento(ii,jj,kk)*feypo(jj) )
    !
    gamman = gammai*gammad / (gammad * (1._DBL-fezpo(kk)) + gammai * fezpo(kk))
    !
    ! gamma_s 
    !
    ! ** se utilizan las constantes gammai y gammad como auxiliares para
    ! calcular gamman, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammad = ( gamma_momento(ii,jj,kk+1) * gamma_momento(ii,jj-1,kk+1) ) / &
         &(gamma_momento(ii,jj,kk+1)*(1._DBL-feypo(jj-1))+&
         &gamma_momento(ii,jj-1,kk+1)*feypo(jj-1))
    gammai = ( gamma_momento(ii,jj,kk) * gamma_momento(ii,jj-1,kk) ) / &
         &(gamma_momento(ii,jj,kk)*(1._DBL-feypo(jj-1))+&
         &gamma_momento(ii,jj-1,kk)*feypo(jj-1))
    !
    gammas = gammai*gammad / (gammad * (1._DBL-fezpo(kk)) + gammai * fezpo(kk))
    !
    ! gamma_d 
    !
    ! ** se utilizan las constantes gammat y gammab como auxiliares para
    ! calcular gammad, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammat = ( gamma_momento(ii+1,jj,kk+1) * gamma_momento(ii,jj,kk+1) ) / &
         &(gamma_momento(ii+1,jj,kk+1)*(1._DBL-fexpo(ii))+&
         &gamma_momento(ii,jj,kk+1)*fexpo(ii))
    gammab = ( gamma_momento(ii+1,jj,kk) * gamma_momento(ii,jj,kk) ) / &
         &(gamma_momento(ii+1,jj,kk)*(1._DBL-fexpo(ii))+&
         &gamma_momento(ii,jj,kk)*fexpo(ii))
    !
    gammad = gammat*gammab / (gammat * (1._DBL-fezpo(kk)) + gammab * fezpo(kk))
    !
    ! gamma_i
    !
    ! ** se utilizan las constantes gammat y gammab como auxiliares para
    ! calcular gammai, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammat = ( gamma_momento(ii,jj,kk+1) * gamma_momento(ii-1,jj,kk+1) ) / &
         &( gamma_momento(ii,jj,kk+1)*(1._DBL-fexpo(ii-1))+&
         &gamma_momento(ii-1,jj,kk+1)*fexpo(ii-1) )
    gammab = ( gamma_momento(ii,jj,kk) * gamma_momento(ii-1,jj,kk) ) / &
         &( gamma_momento(ii,jj,kk) * (1._DBL-fexpo(ii-1))+&
         &gamma_momento(ii-1,jj,kk)*fexpo(ii-1) )
    !
    gammai = gammat*gammab / (gammat * (1._DBL-fezpo(kk)) + gammab * fezpo(kk))
    !
    ! gamma_t
    !
    gammat = gamma_momento(ii,jj,kk+1)
    !
    ! gamma_b
    !
    gammab = gamma_momento(ii,jj,kk)
    !
    ! distancias entre nodos contiguos
    !
    di = deltaxuo(ii-1)
    dd = deltaxuo(ii)
    ds = deltayvo(jj-1)
    dn = deltayvo(jj)
    db = deltazpo(kk)
    dt = deltazpo(kk+1)
    !
    ! Tama\~no de los vol\'umenes de control para la velocidad u
    !
    deltax = deltaxwo(ii)
    deltay = deltaywo(jj)
    deltaz = deltazwo(kk)
    !
    ! Interpolaci\'on para la temperatura
    !
    temp_int = fezpo(kk)*temp_o(ii,jj,kk+1) + (1.0_DBL-fezpo(kk))*temp_o(ii,jj,kk)
    !
    ! *************************
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
    alpha =-(gammas*deltax*deltaz/ds*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vs*ds/gammas))**5)+&
         &DMAX1(0.0_DBL, vs*deltax*deltaz))
    !
    beta  =-(gamman*deltax*deltaz/dn*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vn*dn/gamman))**5)+&
         &DMAX1(0.0_DBL,-vn*deltax*deltaz))
    !
    gamma =-(gammab*deltax*deltay/db*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wb*db/gammab))**5)+&
         &DMAX1(0.0_DBL, wb*deltax*deltay))
    !
    delta =-(gammat*deltax*deltay/dt*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wt*dt/gammat))**5)+&
         &DMAX1(0.0_DBL,-wt*deltax*deltay)) 
    !
    AC_o(indexp(ii,jj,kk)) = ( -AI_o(indexp(ii,jj,kk)) - AD_o(indexp(ii,jj,kk)) &
         &- alpha - beta - gamma - delta - &
         &deltax*deltay*deltaz*fuente_lin_wo(ii,jj,kk)+&
         &deltax*deltay*deltaz/dt_o ) / rel_vo
    !
    Rx_o(indexp(ii,jj,kk)) =-alpha*w_o(ii,jj-1,kk) -&
         &beta  * w_o(ii,jj+1,kk) - &
         &gamma * w_o(ii,jj,kk-1) - &
         &delta * w_o(ii,jj,kk+1) - &
         &deltax*deltay*deltaz*Ri_o(ii,jj,kk)*temp_int+&
         &deltax*deltay*deltaz*fuente_con_wo(ii,jj,kk)+&
         &deltax*deltay*deltaz*w_anto(ii,jj,kk)/dt_o+&
         &(pres_o(ii,jj,kk)-pres_o(ii,jj,kk+1))*deltax*deltay+&
         &AC_o(indexp(ii,jj,kk))*(1._DBL-rel_vo)*w_o(ii,jj,kk)
    !
    aw_o(ii,jj,kk) = AC_o(indexp(ii,jj,kk)) * rel_vo
    !
  end subroutine ensambla_velw_x
  !
  !-------------------------------------------------------------------  
  !*******************************************************************
  !
  ! ensambla_velv_x
  !
  ! Subrutina que calcula los coeficientes de la matriz tridiagonal
  ! para la velocidad u en la direcci\'on x
  !
  !*******************************************************************
  !-------------------------------------------------------------------
  subroutine ensambla_velv_z(&
       &deltaxvo,&
       &deltayvo,&
       &deltazvo,&
       &deltaxuo,&
       &deltaypo,&
       &deltazwo,&
       &fexpo,&
       &feypo,&
       &fezpo,&
       &feyvo,&
       &gamma_momento,&
       &u_o,&
       &v_o,&
       &v_anto,&
       &w_o,&
       &temp_o,&
       &pres_o,&
       &fuente_con_vo,&
       &fuente_lin_vo,&
       &Ri_o,&
       &dt_o,&
       &rel_vo,&
       &AI_o,AC_o,AD_o,Rx_o,&
       &kk,jj,ii&
       &)
    implicit none
    !$acc routine
    !
    ! Tama\~no del volumen de control
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxvo
    real(kind=DBL), dimension(nj), intent(in) :: deltayvo
    real(kind=DBL), dimension(lk), intent(in) :: deltazvo
    !
    ! Distancia entre nodos contiguos de la malla de v en direcci\'on horizontal x
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxuo
    !
    ! Distancia entre nodos contiguos de la malla de v en direcci\'on horizontal y
    !
    real(kind=DBL), dimension(nj), intent(in) :: deltaypo
    !
    ! Distancia entre nodos contiguos de la malla de u en direcci\'on vertical z
    !
    real(kind=DBL), dimension(lk), intent(in) :: deltazwo
    !
    ! Coeficientes para interpolaci\'on
    !
    real(kind=DBL), DIMENSION(mi),   intent(in) :: fexpo
    real(kind=DBL), DIMENSION(nj),   intent(in) :: feypo
    real(kind=DBL), dimension(lk),   intent(in) :: fezpo
    real(kind=DBL), DIMENSION(nj-1), intent(in) :: feyvo
    !
    ! Coeficiente de difusi\'on
    !
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in) :: gamma_momento
    !
    ! Velocidad, presi\'on y temperatura
    !
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in) :: u_o
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in) :: v_o, v_anto
    real(kind=DBL), dimension(mi+1,nj+1,lk),   intent(in) :: w_o
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in) :: temp_o, pres_o
    !
    ! T\'erminos fuente
    !
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in) :: fuente_con_vo
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in) :: fuente_lin_vo    
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in) :: Ri_o
    
    !
    ! Incremento de tiempo y coeficiente de relajaci\'on
    !
    real(kind=DBL), intent(in) :: dt_o, rel_vo
    !
    ! \'Indices para recorrer las direcciones x, y, z
    !
    integer, intent(in)        :: ii, jj, kk
    !
    ! Coeficientes de las matrices
    !
    ! ** Estos coeficientes est\'an sobredimensionados para reducir el uso de memoria
    ! en la gpu, los arreglos que se reciben en esta subrutina se usan para las ecs.
    ! de momento, energ\'ia y la correcci\'on de la presi\'on **
    !
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AI_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AC_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AD_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: RX_o
    !
    ! Variables auxiliares
    !
    integer :: info
    !
    ! Auxiliares de interpolaci\'on y coeficientes
    !
    real(kind=DBL) :: ui, ud, vs, vn, wb, wt
    real(kind=DBL) :: di, dd, ds, dn, db, dt
    real(kind=DBL) :: gammai, gammad
    real(kind=DBL) :: gammas, gamman
    real(kind=DBL) :: gammab, gammat
    real(kind=DBL) :: alpha, beta, gamma, delta
    real(kind=DBL) :: deltax, deltay, deltaz
    real(kind=DBL) :: temp_int
    !
    ! Interpolaciones necesarias
    !
    ! u
    !
    ud = feypo(jj)*u_o(ii,jj+1,kk)  +(1.0_DBL-feypo(jj)) * u_o(ii,jj,kk)
    ui = feypo(jj)*u_o(ii-1,jj+1,kk)+(1.0_DBL-feypo(jj)) * u_o(ii-1,jj,kk)
    !
    ! v
    !
    vn = feyvo(jj)  *v_o(ii,jj+1,kk)+(1.0_DBL-feyvo(jj))   * v_o(ii,jj,kk)
    vs = feyvo(jj-1)*v_o(ii,jj,kk)  +(1.0_DBL-feyvo(jj-1)) * v_o(ii,jj-1,kk)
    !
    ! w
    !
    wt = feypo(jj)*w_o(ii,jj+1,kk)  +(1.0_DBL-feypo(jj)) * w_o(ii,jj,kk)
    wb = feypo(jj)*w_o(ii,jj+1,kk-1)+(1.0_DBL-feypo(jj)) * w_o(ii,jj,kk-1)
    !
    ! gamma_d
    !
    ! ** se utilizan las constantes gamman y gammas como auxiliares para
    ! calcular gammad, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gamman = ( gamma_momento(ii+1,jj+1,kk) * gamma_momento(ii,jj+1,kk) ) / &
         &(gamma_momento(ii+1,jj+1,kk)*(1._DBL-fexpo(ii))+&
         &gamma_momento(ii,jj+1,kk)*fexpo(ii) )
    gammas = ( gamma_momento(ii+1,jj,kk) * gamma_momento(ii,jj,kk) ) / &
         &( gamma_momento(ii+1,jj,kk) * (1._DBL-fexpo(ii))+&
         &gamma_momento(ii,jj,kk)*fexpo(jj) )
    !
    gammad = gammas*gamman / (gammas * (1._DBL-feypo(jj)) + gamman * feypo(jj))
    !
    ! gamma_i
    !
    ! ** se utilizan las constantes gamman y gammas como auxiliares para
    ! calcular gammai, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gamman = ( gamma_momento(ii,jj+1,kk) * gamma_momento(ii-1,jj+1,kk) ) / &
         &(gamma_momento(ii,jj+1,kk)*(1._DBL-fexpo(ii-1))+&
         &gamma_momento(ii-1,jj+1,kk)*fexpo(ii-1))
    gammas = ( gamma_momento(ii,jj,kk) * gamma_momento(ii-1,jj,kk) ) / &
         &(gamma_momento(ii,jj,kk)*(1._DBL-fexpo(ii-1))+&
         &gamma_momento(ii-1,jj,kk)*fexpo(ii-1))
    !
    gammai = gamman*gammas / (gammas * (1._DBL-feypo(jj)) + gamman * feypo(jj))
    !
    ! gamma_t
    !
    ! ** se utilizan las constantes gamman y gammas como auxiliares para
    ! calcular gammat, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gamman = ( gamma_momento(ii,jj+1,kk+1) * gamma_momento(ii,jj+1,kk) ) / &
         &(gamma_momento(ii,jj+1,kk+1)*(1._DBL-fezpo(kk))+&
         &gamma_momento(ii,jj+1,kk)*fezpo(kk) )
    !
    gammas = ( gamma_momento(ii,jj,kk+1) * gamma_momento(ii,jj,kk) ) / &
         &( gamma_momento(ii,jj,kk+1) * (1._DBL-fezpo(kk))+&
         &gamma_momento(ii,jj,kk)*fezpo(kk) ) 
    !
    gammat = gamman*gammas / (gammas * (1._DBL-feypo(jj)) + gamman * feypo(jj))
    !
    ! gamma_b 
    !
    ! ** se utilizan las constantes gamman y gammas como auxiliares para
    ! calcular gammab, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gamman = ( gamma_momento(ii,jj+1,kk) * gamma_momento(ii,jj+1,kk-1) ) / &
         &(gamma_momento(ii,jj+1,kk)*(1._DBL-fezpo(kk-1))+&
         &gamma_momento(ii,jj+1,kk-1)*fezpo(kk-1))
    !
    gammas = ( gamma_momento(ii,jj,kk) * gamma_momento(ii,jj,kk-1) ) / &
         &(gamma_momento(ii,jj,kk)*(1._DBL-fezpo(kk-1))+&
         &gamma_momento(ii,jj,kk-1)*fezpo(kk-1))
    !
    gammab =  gamman*gammas / (gammas * (1._DBL-feypo(jj)) + gamman * feypo(jj))
    !
    ! gamma_n
    !
    gamman = gamma_momento(ii,jj+1,kk)
    !
    ! gamma_s
    !
    gammai = gamma_momento(ii,jj,kk)
    !
    ! distancias entre nodos contiguos
    !
    di = deltaxuo(ii-1)
    dd = deltaxuo(ii)
    ds = deltaypo(jj)
    dn = deltaypo(jj+1)
    db = deltazwo(kk-1)
    dt = deltazwo(kk)
    !
    ! Tama\~no de los vol\'umenes de control para la velocidad u
    !
    deltax = deltaxvo(ii)
    deltay = deltayvo(jj)
    deltaz = deltazvo(kk)
    !
    ! Interpolaci\'on para la temperatura
    !
    temp_int = feypo(jj)*temp_o(ii,jj+1,kk) + (1.0_DBL-feypo(jj))*temp_o(ii,jj,kk)
    !
    ! *************************
    !
    ! Coeficientes de la matriz
    !
    alpha =-(gammai*deltay*deltaz/di*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ui*di/gammai))**5)+&
         &DMAX1(0.0_DBL, ui*deltay*deltaz))
    !
    beta =-(gammad*deltay*deltaz/dd*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ud*dd/gammad))**5)+&
         &DMAX1(0.0_DBL,-ud*deltay*deltaz))
    !
    gamma =-(gammas*deltax*deltaz/ds*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vs*ds/gammas))**5)+&
         &DMAX1(0.0_DBL, vs*deltax*deltaz))
    !
    delta =-(gamman*deltax*deltaz/dn*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vn*dn/gamman))**5)+&
         &DMAX1(0.0_DBL,-vn*deltax*deltaz))
    !
    AI_o(indezv(kk,jj,ii)) =-(gammab*deltax*deltay/db*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wb*db/gammab))**5)+&
         &DMAX1(0.0_DBL, wb*deltax*deltay))
    !
    AD_o(indezv(kk,jj,ii)) =-(gammat*deltax*deltay/dt*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wt*dt/gammat))**5)+&
         &DMAX1(0.0_DBL,-wt*deltax*deltay)) 
    !
    AC_o(indezv(kk,jj,ii)) = ( -AI_o(indezv(kk,jj,ii)) - AD_o(indezv(kk,jj,ii)) &
         &- alpha - beta - gamma - delta - &
         &deltax*deltay*deltaz*fuente_lin_vo(ii,jj,kk)+&
         &deltax*deltay*deltaz/dt_o ) / rel_vo
    !
    Rx_o(indezv(kk,jj,ii)) =-alpha*v_o(ii-1,jj,kk) -&
         &beta  * v_o(ii+1,jj,kk) - &
         &gamma * v_o(ii,jj-1,kk) - &
         &delta * v_o(ii,jj+1,kk) - &
         &deltax*deltay*deltaz*Ri_o(ii,jj,kk)*temp_int+&
         &deltax*deltay*deltaz*fuente_con_vo(ii,jj,kk)+&
         &deltax*deltay*deltaz*v_anto(ii,jj,kk)/dt_o+&
         &(pres_o(ii,jj,kk)-pres_o(ii,jj+1,kk))*deltax*deltaz+&
         &AC_o(indezv(kk,jj,ii))*(1._DBL-rel_vo)*v_o(ii,jj,kk)
    !
  end subroutine ensambla_velv_z
  !
  !
  !-------------------------------------------------------------------  
  !*******************************************************************
  !
  ! ensambla_velv_x
  !
  ! Subrutina que calcula los coeficientes de la matriz tridiagonal
  ! para la velocidad u en la direcci\'on x
  !
  !*******************************************************************
  !-------------------------------------------------------------------
  subroutine ensambla_velv_y(&
       &deltaxvo,&
       &deltayvo,&
       &deltazvo,&
       &deltaxuo,&
       &deltaypo,&
       &deltazwo,&
       &fexpo,&
       &feypo,&
       &fezpo,&
       &feyvo,&
       &gamma_momento,&
       &u_o,&
       &v_o,&
       &v_anto,&
       &w_o,&
       &temp_o,&
       &pres_o,&
       &fuente_con_vo,&
       &fuente_lin_vo,&
       &Ri_o,&
       &dt_o,&
       &rel_vo,&
       &AI_o,AC_o,AD_o,Rx_o,&
       &jj,ii,kk&
       &)
    implicit none
    !$acc routine
    !
    ! Tama\~no del volumen de control
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxvo
    real(kind=DBL), dimension(nj), intent(in) :: deltayvo
    real(kind=DBL), dimension(lk), intent(in) :: deltazvo
    !
    ! Distancia entre nodos contiguos de la malla de v en direcci\'on horizontal x
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxuo
    !
    ! Distancia entre nodos contiguos de la malla de v en direcci\'on horizontal y
    !
    real(kind=DBL), dimension(nj), intent(in) :: deltaypo
    !
    ! Distancia entre nodos contiguos de la malla de u en direcci\'on vertical z
    !
    real(kind=DBL), dimension(lk), intent(in) :: deltazwo
    !
    ! Coeficientes para interpolaci\'on
    !
    real(kind=DBL), DIMENSION(mi),   intent(in) :: fexpo
    real(kind=DBL), DIMENSION(nj),   intent(in) :: feypo
    real(kind=DBL), dimension(lk),   intent(in) :: fezpo
    real(kind=DBL), DIMENSION(nj-1), intent(in) :: feyvo
    !
    ! Coeficiente de difusi\'on
    !
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in) :: gamma_momento
    !
    ! Velocidad, presi\'on y temperatura
    !
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in) :: u_o
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in) :: v_o, v_anto
    real(kind=DBL), dimension(mi+1,nj+1,lk),   intent(in) :: w_o
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in) :: temp_o, pres_o
    !
    ! T\'erminos fuente
    !
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in) :: fuente_con_vo
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in) :: fuente_lin_vo    
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in) :: Ri_o
    
    !
    ! Incremento de tiempo y coeficiente de relajaci\'on
    !
    real(kind=DBL), intent(in) :: dt_o, rel_vo
    !
    ! \'Indices para recorrer las direcciones x, y, z
    !
    integer, intent(in)        :: ii, jj, kk
    !
    ! Coeficientes de las matrices
    !
    ! ** Estos coeficientes est\'an sobredimensionados para reducir el uso de memoria
    ! en la gpu, los arreglos que se reciben en esta subrutina se usan para las ecs.
    ! de momento, energ\'ia y la correcci\'on de la presi\'on **
    !
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AI_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AC_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AD_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: RX_o
    !
    ! Variables auxiliares
    !
    integer :: info
    !
    ! Auxiliares de interpolaci\'on y coeficientes
    !
    real(kind=DBL) :: ui, ud, vs, vn, wb, wt
    real(kind=DBL) :: di, dd, ds, dn, db, dt
    real(kind=DBL) :: gammai, gammad
    real(kind=DBL) :: gammas, gamman
    real(kind=DBL) :: gammab, gammat
    real(kind=DBL) :: alpha, beta, gamma, delta
    real(kind=DBL) :: deltax, deltay, deltaz
    real(kind=DBL) :: temp_int
    !
    ! Interpolaciones necesarias
    !
    ! u
    !
    ud = feypo(jj)*u_o(ii,jj+1,kk)  +(1.0_DBL-feypo(jj)) * u_o(ii,jj,kk)
    ui = feypo(jj)*u_o(ii-1,jj+1,kk)+(1.0_DBL-feypo(jj)) * u_o(ii-1,jj,kk)
    !
    ! v
    !
    vn = feyvo(jj)  *v_o(ii,jj+1,kk)+(1.0_DBL-feyvo(jj))   * v_o(ii,jj,kk)
    vs = feyvo(jj-1)*v_o(ii,jj,kk)  +(1.0_DBL-feyvo(jj-1)) * v_o(ii,jj-1,kk)
    !
    ! w
    !
    wt = feypo(jj)*w_o(ii,jj+1,kk)  +(1.0_DBL-feypo(jj)) * w_o(ii,jj,kk)
    wb = feypo(jj)*w_o(ii,jj+1,kk-1)+(1.0_DBL-feypo(jj)) * w_o(ii,jj,kk-1)
    !
    ! gamma_d
    !
    ! ** se utilizan las constantes gamman y gammas como auxiliares para
    ! calcular gammad, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gamman = ( gamma_momento(ii+1,jj+1,kk) * gamma_momento(ii,jj+1,kk) ) / &
         &(gamma_momento(ii+1,jj+1,kk)*(1._DBL-fexpo(ii))+&
         &gamma_momento(ii,jj+1,kk)*fexpo(ii) )
    gammas = ( gamma_momento(ii+1,jj,kk) * gamma_momento(ii,jj,kk) ) / &
         &( gamma_momento(ii+1,jj,kk) * (1._DBL-fexpo(ii))+&
         &gamma_momento(ii,jj,kk)*fexpo(jj) )
    !
    gammad = gammas*gamman / (gammas * (1._DBL-feypo(jj)) + gamman * feypo(jj))
    !
    ! gamma_i
    !
    ! ** se utilizan las constantes gamman y gammas como auxiliares para
    ! calcular gammai, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gamman = ( gamma_momento(ii,jj+1,kk) * gamma_momento(ii-1,jj+1,kk) ) / &
         &(gamma_momento(ii,jj+1,kk)*(1._DBL-fexpo(ii-1))+&
         &gamma_momento(ii-1,jj+1,kk)*fexpo(ii-1))
    gammas = ( gamma_momento(ii,jj,kk) * gamma_momento(ii-1,jj,kk) ) / &
         &(gamma_momento(ii,jj,kk)*(1._DBL-fexpo(ii-1))+&
         &gamma_momento(ii-1,jj,kk)*fexpo(ii-1))
    !
    gammai = gamman*gammas / (gammas * (1._DBL-feypo(jj)) + gamman * feypo(jj))
    !
    ! gamma_t
    !
    ! ** se utilizan las constantes gamman y gammas como auxiliares para
    ! calcular gammat, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gamman = ( gamma_momento(ii,jj+1,kk+1) * gamma_momento(ii,jj+1,kk) ) / &
         &(gamma_momento(ii,jj+1,kk+1)*(1._DBL-fezpo(kk))+&
         &gamma_momento(ii,jj+1,kk)*fezpo(kk) )
    !
    gammas = ( gamma_momento(ii,jj,kk+1) * gamma_momento(ii,jj,kk) ) / &
         &( gamma_momento(ii,jj,kk+1) * (1._DBL-fezpo(kk))+&
         &gamma_momento(ii,jj,kk)*fezpo(kk) ) 
    !
    gammat = gamman*gammas / (gammas * (1._DBL-feypo(jj)) + gamman * feypo(jj))
    !
    ! gamma_b 
    !
    ! ** se utilizan las constantes gamman y gammas como auxiliares para
    ! calcular gammab, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gamman = ( gamma_momento(ii,jj+1,kk) * gamma_momento(ii,jj+1,kk-1) ) / &
         &(gamma_momento(ii,jj+1,kk)*(1._DBL-fezpo(kk-1))+&
         &gamma_momento(ii,jj+1,kk-1)*fezpo(kk-1))
    !
    gammas = ( gamma_momento(ii,jj,kk) * gamma_momento(ii,jj,kk-1) ) / &
         &(gamma_momento(ii,jj,kk)*(1._DBL-fezpo(kk-1))+&
         &gamma_momento(ii,jj,kk-1)*fezpo(kk-1))
    !
    gammab =  gamman*gammas / (gammas * (1._DBL-feypo(jj)) + gamman * feypo(jj))
    !
    ! gamma_n
    !
    gamman = gamma_momento(ii,jj+1,kk)
    !
    ! gamma_s
    !
    gammai = gamma_momento(ii,jj,kk)
    !
    ! distancias entre nodos contiguos
    !
    di = deltaxuo(ii-1)
    dd = deltaxuo(ii)
    ds = deltaypo(jj)
    dn = deltaypo(jj+1)
    db = deltazwo(kk-1)
    dt = deltazwo(kk)
    !
    ! Tama\~no de los vol\'umenes de control para la velocidad u
    !
    deltax = deltaxvo(ii)
    deltay = deltayvo(jj)
    deltaz = deltazvo(kk)
    !
    ! Interpolaci\'on para la temperatura
    !
    temp_int = feypo(jj)*temp_o(ii,jj+1,kk) + (1.0_DBL-feypo(jj))*temp_o(ii,jj,kk)
    !
    ! *************************
    !
    ! Coeficientes de la matriz
    !
    alpha =-(gammai*deltay*deltaz/di*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ui*di/gammai))**5)+&
         &DMAX1(0.0_DBL, ui*deltay*deltaz))
    !
    beta =-(gammad*deltay*deltaz/dd*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ud*dd/gammad))**5)+&
         &DMAX1(0.0_DBL,-ud*deltay*deltaz))
    !
    AI_o(indeyv(jj,ii,kk)) =-(gammas*deltax*deltaz/ds*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vs*ds/gammas))**5)+&
         &DMAX1(0.0_DBL, vs*deltax*deltaz))
    !
    AD_o(indeyv(jj,ii,kk))  =-(gamman*deltax*deltaz/dn*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vn*dn/gamman))**5)+&
         &DMAX1(0.0_DBL,-vn*deltax*deltaz))
    !
    gamma =-(gammab*deltax*deltay/db*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wb*db/gammab))**5)+&
         &DMAX1(0.0_DBL, wb*deltax*deltay))
    !
    delta =-(gammat*deltax*deltay/dt*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wt*dt/gammat))**5)+&
         &DMAX1(0.0_DBL,-wt*deltax*deltay)) 
    !
    AC_o(indeyv(jj,ii,kk)) = ( -AI_o(indeyv(jj,ii,kk)) - AD_o(indeyv(jj,ii,kk)) &
         &- alpha - beta - gamma - delta - &
         &deltax*deltay*deltaz*fuente_lin_vo(ii,jj,kk)+&
         &deltax*deltay*deltaz/dt_o ) / rel_vo
    !
    Rx_o(indeyv(jj,ii,kk)) =-alpha*v_o(ii-1,jj,kk) -&
         &beta  * v_o(ii+1,jj,kk) - &
         &gamma * v_o(ii,jj,kk-1) - &
         &delta * v_o(ii,jj,kk+1) - &
         &deltax*deltay*deltaz*Ri_o(ii,jj,kk)*temp_int+&
         &deltax*deltay*deltaz*fuente_con_vo(ii,jj,kk)+&
         &deltax*deltay*deltaz*v_anto(ii,jj,kk)/dt_o+&
         &(pres_o(ii,jj,kk)-pres_o(ii,jj+1,kk))*deltax*deltaz+&
         &AC_o(indeyv(jj,ii,kk))*(1._DBL-rel_vo)*v_o(ii,jj,kk)
    !
  end subroutine ensambla_velv_y
  !
  !
  !-------------------------------------------------------------------  
  !*******************************************************************
  !
  ! ensambla_velu_x
  !
  ! Subrutina que calcula los coeficientes de la matriz tridiagonal
  ! para la velocidad u en la direcci\'on x
  !
  !*******************************************************************
  !-------------------------------------------------------------------
  subroutine ensambla_velv_x(&
       &deltaxvo,&
       &deltayvo,&
       &deltazvo,&
       &deltaxuo,&
       &deltaypo,&
       &deltazwo,&
       &fexpo,&
       &feypo,&
       &fezpo,&
       &feyvo,&
       &gamma_momento,&
       &u_o,&
       &v_o,&
       &v_anto,&
       &w_o,&
       &temp_o,&
       &pres_o,&
       &fuente_con_vo,&
       &fuente_lin_vo,&
       &Ri_o,&
       &dt_o,&
       &rel_vo,&
       &AI_o,AC_o,AD_o,Rx_o,&
       &av_o,&
       &ii,jj,kk&
       &)
    implicit none
    !$acc routine
    !
    ! Tama\~no del volumen de control
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxvo
    real(kind=DBL), dimension(nj), intent(in) :: deltayvo
    real(kind=DBL), dimension(lk), intent(in) :: deltazvo
    !
    ! Distancia entre nodos contiguos de la malla de v en direcci\'on horizontal x
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxuo
    !
    ! Distancia entre nodos contiguos de la malla de v en direcci\'on horizontal y
    !
    real(kind=DBL), dimension(nj), intent(in) :: deltaypo
    !
    ! Distancia entre nodos contiguos de la malla de u en direcci\'on vertical z
    !
    real(kind=DBL), dimension(lk), intent(in) :: deltazwo
    !
    ! Coeficientes para interpolaci\'on
    !
    real(kind=DBL), DIMENSION(mi),   intent(in) :: fexpo
    real(kind=DBL), DIMENSION(nj),   intent(in) :: feypo
    real(kind=DBL), dimension(lk),   intent(in) :: fezpo
    real(kind=DBL), DIMENSION(nj-1), intent(in) :: feyvo
    !
    ! Coeficiente de difusi\'on
    !
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in) :: gamma_momento
    !
    ! Velocidad, presi\'on y temperatura
    !
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in) :: u_o
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in) :: v_o, v_anto
    real(kind=DBL), dimension(mi+1,nj+1,lk),   intent(in) :: w_o
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in) :: temp_o, pres_o
    !
    ! T\'erminos fuente
    !
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in) :: fuente_con_vo
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in) :: fuente_lin_vo    
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in) :: Ri_o
    
    !
    ! Incremento de tiempo y coeficiente de relajaci\'on
    !
    real(kind=DBL), intent(in) :: dt_o, rel_vo
    !
    ! \'Indices para recorrer las direcciones x, y, z
    !
    integer, intent(in)        :: ii, jj, kk
    !
    ! Coeficientes de las matrices
    !
    ! ** Estos coeficientes est\'an sobredimensionados para reducir el uso de memoria
    ! en la gpu, los arreglos que se reciben en esta subrutina se usan para las ecs.
    ! de momento, energ\'ia y la correcci\'on de la presi\'on **
    !
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AI_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AC_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AD_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: RX_o
    !
    real(kind=DBL), dimension(mi+1,nj,lk+1),         intent(out) :: av_o
    !
    ! Variables auxiliares
    !
    integer :: info
    !
    ! Auxiliares de interpolaci\'on y coeficientes
    !
    real(kind=DBL) :: ui, ud, vs, vn, wb, wt
    real(kind=DBL) :: di, dd, ds, dn, db, dt
    real(kind=DBL) :: gammai, gammad
    real(kind=DBL) :: gammas, gamman
    real(kind=DBL) :: gammab, gammat
    real(kind=DBL) :: alpha, beta, gamma, delta
    real(kind=DBL) :: deltax, deltay, deltaz
    real(kind=DBL) :: temp_int
    !
    ! Interpolaciones necesarias
    !
    ! u
    !
    ud = feypo(jj)*u_o(ii,jj+1,kk)  +(1.0_DBL-feypo(jj)) * u_o(ii,jj,kk)
    ui = feypo(jj)*u_o(ii-1,jj+1,kk)+(1.0_DBL-feypo(jj)) * u_o(ii-1,jj,kk)
    !
    ! v
    !
    vn = feyvo(jj)  *v_o(ii,jj+1,kk)+(1.0_DBL-feyvo(jj))   * v_o(ii,jj,kk)
    vs = feyvo(jj-1)*v_o(ii,jj,kk)  +(1.0_DBL-feyvo(jj-1)) * v_o(ii,jj-1,kk)
    !
    ! w
    !
    wt = feypo(jj)*w_o(ii,jj+1,kk)  +(1.0_DBL-feypo(jj)) * w_o(ii,jj,kk)
    wb = feypo(jj)*w_o(ii,jj+1,kk-1)+(1.0_DBL-feypo(jj)) * w_o(ii,jj,kk-1)
    !
    ! gamma_d
    !
    ! ** se utilizan las constantes gamman y gammas como auxiliares para
    ! calcular gammad, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gamman = ( gamma_momento(ii+1,jj+1,kk) * gamma_momento(ii,jj+1,kk) ) / &
         &(gamma_momento(ii+1,jj+1,kk)*(1._DBL-fexpo(ii))+&
         &gamma_momento(ii,jj+1,kk)*fexpo(ii) )
    gammas = ( gamma_momento(ii+1,jj,kk) * gamma_momento(ii,jj,kk) ) / &
         &( gamma_momento(ii+1,jj,kk) * (1._DBL-fexpo(ii))+&
         &gamma_momento(ii,jj,kk)*fexpo(jj) )
    !
    gammad = gammas*gamman / (gammas * (1._DBL-feypo(jj)) + gamman * feypo(jj))
    !
    ! gamma_i
    !
    ! ** se utilizan las constantes gamman y gammas como auxiliares para
    ! calcular gammai, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gamman = ( gamma_momento(ii,jj+1,kk) * gamma_momento(ii-1,jj+1,kk) ) / &
         &(gamma_momento(ii,jj+1,kk)*(1._DBL-fexpo(ii-1))+&
         &gamma_momento(ii-1,jj+1,kk)*fexpo(ii-1))
    gammas = ( gamma_momento(ii,jj,kk) * gamma_momento(ii-1,jj,kk) ) / &
         &(gamma_momento(ii,jj,kk)*(1._DBL-fexpo(ii-1))+&
         &gamma_momento(ii-1,jj,kk)*fexpo(ii-1))
    !
    gammai = gamman*gammas / (gammas * (1._DBL-feypo(jj)) + gamman * feypo(jj))
    !
    ! gamma_t
    !
    ! ** se utilizan las constantes gamman y gammas como auxiliares para
    ! calcular gammat, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gamman = ( gamma_momento(ii,jj+1,kk+1) * gamma_momento(ii,jj+1,kk) ) / &
         &(gamma_momento(ii,jj+1,kk+1)*(1._DBL-fezpo(kk))+&
         &gamma_momento(ii,jj+1,kk)*fezpo(kk) )
    !
    gammas = ( gamma_momento(ii,jj,kk+1) * gamma_momento(ii,jj,kk) ) / &
         &( gamma_momento(ii,jj,kk+1) * (1._DBL-fezpo(kk))+&
         &gamma_momento(ii,jj,kk)*fezpo(kk) ) 
    !
    gammat = gamman*gammas / (gammas * (1._DBL-feypo(jj)) + gamman * feypo(jj))
    !
    ! gamma_b 
    !
    ! ** se utilizan las constantes gamman y gammas como auxiliares para
    ! calcular gammab, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gamman = ( gamma_momento(ii,jj+1,kk) * gamma_momento(ii,jj+1,kk-1) ) / &
         &(gamma_momento(ii,jj+1,kk)*(1._DBL-fezpo(kk-1))+&
         &gamma_momento(ii,jj+1,kk-1)*fezpo(kk-1))
    !
    gammas = ( gamma_momento(ii,jj,kk) * gamma_momento(ii,jj,kk-1) ) / &
         &(gamma_momento(ii,jj,kk)*(1._DBL-fezpo(kk-1))+&
         &gamma_momento(ii,jj,kk-1)*fezpo(kk-1))
    !
    gammab =  gamman*gammas / (gammas * (1._DBL-feypo(jj)) + gamman * feypo(jj))
    !
    ! gamma_n
    !
    gamman = gamma_momento(ii,jj+1,kk)
    !
    ! gamma_s
    !
    gammai = gamma_momento(ii,jj,kk)
    !
    ! distancias entre nodos contiguos
    !
    di = deltaxuo(ii-1)
    dd = deltaxuo(ii)
    ds = deltaypo(jj)
    dn = deltaypo(jj+1)
    db = deltazwo(kk-1)
    dt = deltazwo(kk)
    !
    ! Tama\~no de los vol\'umenes de control para la velocidad u
    !
    deltax = deltaxvo(ii)
    deltay = deltayvo(jj)
    deltaz = deltazvo(kk)
    !
    ! Interpolaci\'on para la temperatura
    !
    temp_int = feypo(jj)*temp_o(ii,jj+1,kk) + (1.0_DBL-feypo(jj))*temp_o(ii,jj,kk)
    !
    ! *************************
    !
    ! Coeficientes de la matriz
    !
    AI_o(indexv(ii,jj,kk)) =-(gammai*deltay*deltaz/di*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ui*di/gammai))**5)+&
         &DMAX1(0.0_DBL, ui*deltay*deltaz))
    !
    AD_o(indexv(ii,jj,kk)) =-(gammad*deltay*deltaz/dd*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ud*dd/gammad))**5)+&
         &DMAX1(0.0_DBL,-ud*deltay*deltaz))
    !
    alpha =-(gammas*deltax*deltaz/ds*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vs*ds/gammas))**5)+&
         &DMAX1(0.0_DBL, vs*deltax*deltaz))
    !
    beta  =-(gamman*deltax*deltaz/dn*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vn*dn/gamman))**5)+&
         &DMAX1(0.0_DBL,-vn*deltax*deltaz))
    !
    gamma =-(gammab*deltax*deltay/db*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wb*db/gammab))**5)+&
         &DMAX1(0.0_DBL, wb*deltax*deltay))
    !
    delta =-(gammat*deltax*deltay/dt*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wt*dt/gammat))**5)+&
         &DMAX1(0.0_DBL,-wt*deltax*deltay)) 
    !
    AC_o(indexv(ii,jj,kk)) = ( -AI_o(indexv(ii,jj,kk)) - AD_o(indexv(ii,jj,kk)) &
         &- alpha - beta - gamma - delta - &
         &deltax*deltay*deltaz*fuente_lin_vo(ii,jj,kk)+&
         &deltax*deltay*deltaz/dt_o ) / rel_vo
    !
    Rx_o(indexv(ii,jj,kk)) =-alpha*v_o(ii,jj-1,kk) -&
         &beta  * v_o(ii,jj+1,kk) - &
         &gamma * v_o(ii,jj,kk-1) - &
         &delta * v_o(ii,jj,kk+1) - &
         &deltax*deltay*deltaz*Ri_o(ii,jj,kk)*temp_int+&
         &deltax*deltay*deltaz*fuente_con_vo(ii,jj,kk)+&
         &deltax*deltay*deltaz*v_anto(ii,jj,kk)/dt_o+&
         &(pres_o(ii,jj,kk)-pres_o(ii,jj+1,kk))*deltax*deltaz+&
         &AC_o(indexv(ii,jj,kk))*(1._DBL-rel_vo)*v_o(ii,jj,kk)
    !
    av_o(ii,jj,kk) = AC_o(indexv(ii,jj,kk)) * rel_vo
    !
  end subroutine ensambla_velv_x
  !
  !
  !-------------------------------------------------------------------  
  !*******************************************************************
  !
  ! ensambla_velu_x
  !
  ! Subrutina que calcula los coeficientes de la matriz tridiagonal
  ! para la velocidad u en la direcci\'on x
  !
  !*******************************************************************
  !-------------------------------------------------------------------
  subroutine ensambla_velu_x(&
       &deltaxuo,&
       &deltayuo,&
       &deltazuo,&
       &deltaxpo,&
       &deltayvo,&
       &deltazwo,&
       &fexpo,&
       &feypo,&
       &fezpo,&
       &fexuo,&
       &gamma_momento,&
       &u_o,&
       &u_anto,&
       &v_o,&
       &w_o,&
       &temp_o,&
       &pres_o,&
       &fuente_con_uo,&
       &fuente_lin_uo,&
       &Ri_o,&
       &dt_o,&
       &rel_vo,&
       &AI_o,AC_o,AD_o,Rx_o,&
       &au_o,&
       &ii,jj,kk&
       &)
    implicit none
    !$acc routine
    !
    ! Tama\~no del volumen de control
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxuo
    real(kind=DBL), dimension(nj), intent(in) :: deltayuo
    real(kind=DBL), dimension(lk), intent(in) :: deltazuo
    !
    ! Distancia entre nodos contiguos de la malla de u en direcci\'on horizontal x
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxpo
    !
    ! Distancia entre nodos contiguos de la malla de u en direcci\'on horizontal y
    !
    real(kind=DBL), dimension(nj), intent(in) :: deltayvo
    !
    ! Distancia entre nodos contiguos de la malla de u en direcci\'on vertical z
    !
    real(kind=DBL), dimension(lk), intent(in) :: deltazwo
    !
    ! Coeficientes para interpolaci\'on
    !
    real(kind=DBL), DIMENSION(mi),   intent(in) :: fexpo
    real(kind=DBL), DIMENSION(nj),   intent(in) :: feypo
    real(kind=DBL), dimension(lk),   intent(in) :: fezpo
    real(kind=DBL), DIMENSION(mi-1), intent(in) :: fexuo
    !
    ! Coeficiente de difusi\'on
    !
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in) :: gamma_momento
    !
    ! Velocidad, presi\'on y temperatura
    !
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in) :: u_o, u_anto
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in) :: v_o
    real(kind=DBL), dimension(mi+1,nj+1,lk),   intent(in) :: w_o
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in) :: temp_o, pres_o
    !
    ! T\'erminos fuente
    !
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in) :: fuente_con_uo
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in) :: fuente_lin_uo    
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in) :: Ri_o
    
    !
    ! Incremento de tiempo y coeficiente de relajaci\'on
    !
    real(kind=DBL), intent(in) :: dt_o, rel_vo
    !
    ! \'Indices para recorrer las direcciones x, y, z
    !
    integer, intent(in)        :: ii, jj, kk
    !
    ! Coeficientes de las matrices
    !
    ! ** Estos coeficientes est\'an sobredimensionados para reducir el uso de memoria
    ! en la gpu, los arreglos que se reciben en esta subrutina se usan para las ecs.
    ! de momento, energ\'ia y la correcci\'on de la presi\'on **
    !
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AI_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AC_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AD_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: RX_o
    !
    real(kind=DBL), dimension(mi,nj+1,lk+1),         intent(out) :: au_o
    !
    ! Variables auxiliares
    !
    integer :: info
    !
    ! Auxiliares de interpolaci\'on y coeficientes
    !
    real(kind=DBL) :: ui, ud, vs, vn, wb, wt
    real(kind=DBL) :: di, dd, ds, dn, db, dt
    real(kind=DBL) :: gammai, gammad
    real(kind=DBL) :: gammas, gamman
    real(kind=DBL) :: gammab, gammat
    real(kind=DBL) :: alpha, beta, gamma, delta
    real(kind=DBL) :: deltax, deltay, deltaz
    real(kind=DBL) :: temp_int
    !
    ! Interpolaciones necesarias
    !
    ! u
    !
    ud = fexuo(ii)  *u_o(ii+1,jj,kk)+(1.0_DBL-fexuo(ii))  *u_o(ii,jj,kk)
    ui = fexuo(ii-1)*u_o(ii,jj,kk)  +(1.0_DBL-fexuo(ii-1))*u_o(ii-1,jj,kk)
    !
    ! v
    !
    vn = fexpo(ii)*v_o(ii+1,jj,kk)  +(1.0_DBL-fexpo(ii))  *v_o(ii,jj,kk)
    vs = fexpo(ii)*v_o(ii+1,jj-1,kk)+(1.0_DBL-fexpo(ii))  *v_o(ii,jj-1,kk)
    !
    ! w
    !
    wt = fexpo(ii)*w_o(ii+1,jj,kk)  +(1.0_DBL-fexpo(ii))  *w_o(ii,jj,kk)
    wb = fexpo(ii)*w_o(ii+1,jj,kk-1)+(1.0_DBL-fexpo(ii))  *w_o(ii,jj,kk-1)
    !
    ! gamma_n 
    !
    ! ** se utilizan las constantes gammai y gammad como auxiliares para
    ! calcular gamman, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammad = ( gamma_momento(ii+1,jj+1,kk) * gamma_momento(ii+1,jj,kk) ) / &
         &(gamma_momento(ii+1,jj+1,kk)*(1._DBL-feypo(jj))+&
         &gamma_momento(ii+1,jj,kk)*feypo(jj) )
    gammai = ( gamma_momento(ii,jj+1,kk) * gamma_momento(ii,jj,kk) ) / &
         &( gamma_momento(ii,jj+1,kk) * (1._DBL-feypo(jj))+&
         &gamma_momento(ii,jj,kk)*feypo(jj) )
    !
    gamman = gammai*gammad / (gammad * (1._DBL-fexpo(ii)) + gammai * fexpo(ii))
    !
    ! gamma_s 
    !
    ! ** se utilizan las constantes gammai y gammad como auxiliares para
    ! calcular gamman, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammad = ( gamma_momento(ii+1,jj,kk) * gamma_momento(ii+1,jj-1,kk) ) / &
         &(gamma_momento(ii+1,jj,kk)*(1._DBL-feypo(jj-1))+&
         &gamma_momento(ii+1,jj-1,kk)*feypo(jj-1))
    gammai = ( gamma_momento(ii,jj,kk) * gamma_momento(ii,jj-1,kk) ) / &
         &(gamma_momento(ii,jj,kk)*(1._DBL-feypo(jj-1))+&
         &gamma_momento(ii,jj-1,kk)*feypo(jj-1))
    !
    gammas = gammai*gammad / (gammad * (1._DBL-fexpo(ii)) + gammai * fexpo(ii))
    !
    ! gamma_t
    !
    ! ** se utilizan las constantes gammai y gammad como auxiliares para
    ! calcular gammat, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammad = ( gamma_momento(ii+1,jj,kk+1) * gamma_momento(ii+1,jj,kk) ) / &
         &( gamma_momento(ii+1,jj,kk+1)*(1._DBL-fezpo(kk))+&
         &gamma_momento(ii+1,jj,kk)*fezpo(kk) )
    gammai = ( gamma_momento(ii,jj,kk+1) * gamma_momento(ii,jj,kk) ) / &
         &( gamma_momento(ii,jj,kk+1) * (1._DBL-fezpo(kk))+&
         &gamma_momento(ii,jj,kk)*fezpo(kk) )
    !
    gammat = gammai*gammad / (gammad * (1._DBL-fexpo(ii)) + gammai * fexpo(ii))
    !
    ! gamma_b 
    !
    ! ** se utilizan las constantes gammai y gammad como auxiliares para
    ! calcular gamman, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammad = ( gamma_momento(ii+1,jj,kk) * gamma_momento(ii+1,jj,kk-1) ) / &
         &(gamma_momento(ii+1,jj,kk)*(1._DBL-fezpo(kk-1))+&
         &gamma_momento(ii+1,jj,kk-1)*fezpo(kk-1))
    gammai = ( gamma_momento(ii,jj,kk+1) * gamma_momento(ii,jj,kk) ) / &
         &(gamma_momento(ii,jj,kk)*(1._DBL-fezpo(kk-1))+&
         &gamma_momento(ii,jj,kk-1)*fezpo(kk-1))
    !
    gammab = gammai*gammad / (gammad * (1._DBL-fexpo(ii)) + gammai * fexpo(ii))
    !
    ! gamma_i
    !
    gammai = gamma_momento(ii,jj,kk)
    !
    ! gamma_d
    !
    gammad = gamma_momento(ii+1,jj,kk)
    !
    ! distancias entre nodos contiguos
    !
    di = deltaxpo(ii)
    dd = deltaxpo(ii+1)
    ds = deltayvo(jj-1)
    dn = deltayvo(jj)
    db = deltazwo(kk-1)
    dt = deltazwo(kk)
    !
    ! Tama\~no de los vol\'umenes de control para la velocidad u
    !
    deltax = deltaxuo(ii)
    deltay = deltayuo(jj)
    deltaz = deltazuo(kk)
    !
    ! Interpolaci\'on para la temperatura
    !
    temp_int = fexpo(ii)*temp_o(ii+1,jj,kk) + (1.0_DBL-fexpo(ii))*temp_o(ii,jj,kk)
    !
    ! *************************
    !
    ! Coeficientes de la matriz
    !
    AI_o(indexu(ii,jj,kk)) =-(gammai*deltay*deltaz/di*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ui*di/gammai))**5)+&
         &DMAX1(0.0_DBL, ui*deltay*deltaz))
    !
    AD_o(indexu(ii,jj,kk)) =-(gammad*deltay*deltaz/dd*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ud*dd/gammad))**5)+&
         &DMAX1(0.0_DBL,-ud*deltay*deltaz))
    !
    alpha =-(gammas*deltax*deltaz/ds*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vs*ds/gammas))**5)+&
         &DMAX1(0.0_DBL, vs*deltax*deltaz))
    !
    beta  =-(gamman*deltax*deltaz/dn*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vn*dn/gamman))**5)+&
         &DMAX1(0.0_DBL,-vn*deltax*deltaz))
    !
    gamma =-(gammab*deltax*deltay/db*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wb*db/gammab))**5)+&
         &DMAX1(0.0_DBL, wb*deltax*deltay))
    !
    delta =-(gammat*deltax*deltay/dt*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wt*dt/gammat))**5)+&
         &DMAX1(0.0_DBL,-wt*deltax*deltay)) 
    !
    AC_o(indexu(ii,jj,kk)) = ( -AI_o(indexu(ii,jj,kk)) - AD_o(indexu(ii,jj,kk)) &
         &- alpha - beta - gamma - delta - &
         &deltax*deltay*deltaz*fuente_lin_uo(ii,jj,kk)+&
         &deltax*deltay*deltaz/dt_o ) / rel_vo
    !
    Rx_o(indexu(ii,jj,kk)) =-alpha*u_o(ii,jj-1,kk) -&
         &beta  * u_o(ii,jj+1,kk) - &
         &gamma * u_o(ii,jj,kk-1) - &
         &delta * u_o(ii,jj,kk+1) - &
         &deltax*deltay*deltaz*Ri_o(ii,jj,kk)*temp_int+&
         &deltax*deltay*deltaz*fuente_con_uo(ii,jj,kk)+&
         &deltax*deltay*deltaz*u_anto(ii,jj,kk)/dt_o+&
         &(pres_o(ii,jj,kk)-pres_o(ii+1,jj,kk))*deltay*deltaz+&
         &AC_o(indexu(ii,jj,kk))*(1._DBL-rel_vo)*u_o(ii,jj,kk)
    !
    au_o(ii,jj,kk) = AC_o(indexu(ii,jj,kk)) * rel_vo
    !
  end subroutine ensambla_velu_x
  !
  !-------------------------------------------------------------------
  !*******************************************************************
  !
  ! ensambla_velu_y
  !
  ! Subrutina que calcula los coeficientes de la matriz tridiagonal
  ! para la velocidad u en la direcci\'on y
  !
  !*******************************************************************
  !-------------------------------------------------------------------
  !
  subroutine ensambla_velu_y(&
       &deltaxuo,&
       &deltayuo,&
       &deltazuo,&
       &deltaxpo,&
       &deltayvo,&
       &deltazwo,&
       &fexpo,&
       &feypo,&
       &fezpo,&
       &fexuo,&
       &gamma_momento,&
       &u_o,&
       &u_anto,&
       &v_o,&
       &w_o,&
       &temp_o,&
       &pres_o,&
       &fuente_con_uo,&
       &fuente_lin_uo,&
       &Ri_o,&
       &dt_o,&
       &rel_vo,&
       &AI_o,AC_o,AD_o,Rx_o,&
       &jj,ii,kk&
       &)
    implicit none
    !$acc routine
    !
    ! Tama\~no del volumen de control
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxuo
    real(kind=DBL), dimension(nj), intent(in) :: deltayuo
    real(kind=DBL), dimension(lk), intent(in) :: deltazuo
    !
    ! Distancia entre nodos contiguos de la malla de u en direcci\'on horizontal x
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxpo
    !
    ! Distancia entre nodos contiguos de la malla de u en direcci\'on horizontal y
    !
    real(kind=DBL), dimension(nj), intent(in) :: deltayvo
    !
    ! Distancia entre nodos contiguos de la malla de u en direcci\'on vertical z
    !
    real(kind=DBL), dimension(lk), intent(in) :: deltazwo
    !
    ! Coeficientes para interpolaci\'on
    !
    real(kind=DBL), DIMENSION(mi),   intent(in) :: fexpo
    real(kind=DBL), DIMENSION(nj),   intent(in) :: feypo
    real(kind=DBL), dimension(lk),   intent(in) :: fezpo
    real(kind=DBL), DIMENSION(mi-1), intent(in) :: fexuo
    !
    ! Coeficiente de difusi\'on
    !
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in) :: gamma_momento
    !
    ! Velocidad, presi\'on y temperatura
    !
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in) :: u_o, u_anto
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in) :: v_o
    real(kind=DBL), dimension(mi+1,nj+1,lk),   intent(in) :: w_o
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in) :: temp_o, pres_o
    !
    ! T\'erminos fuente
    !
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in) :: fuente_con_uo
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in) :: fuente_lin_uo    
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in) :: Ri_o
    
    !
    ! Incremento de tiempo y coeficiente de relajaci\'on
    !
    real(kind=DBL), intent(in) :: dt_o, rel_vo
    !
    ! \'Indices para recorrer las direcciones x, y, z
    !
    integer, intent(in)        :: ii, jj, kk
    !
    ! Coeficientes de las matrices
    !
    ! ** Estos coeficientes est\'an sobredimensionados para reducir el uso de memoria
    ! en la gpu, los arreglos que se reciben en esta subrutina se usan para las ecs.
    ! de momento, energ\'ia y la correcci\'on de la presi\'on **
    !
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AI_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AC_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AD_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: RX_o
    !
    ! real(kind=DBL), dimension(mi,nj+1,lk+1),         intent(out) :: au_o
    !
    ! Variables auxiliares
    !
    integer :: info
    !
    ! Auxiliares de interpolaci\'on y coeficientes
    !
    real(kind=DBL) :: ui, ud, vs, vn, wb, wt
    real(kind=DBL) :: di, dd, ds, dn, db, dt
    real(kind=DBL) :: gammai, gammad
    real(kind=DBL) :: gammas, gamman
    real(kind=DBL) :: gammab, gammat
    real(kind=DBL) :: alpha, beta, gamma, delta
    real(kind=DBL) :: deltax, deltay, deltaz
    real(kind=DBL) :: temp_int
    !
    ! Interpolaciones necesarias
    !
    ! u
    !
    ud = fexuo(ii)  *u_o(ii+1,jj,kk)+(1.0_DBL-fexuo(ii))  *u_o(ii,jj,kk)
    ui = fexuo(ii-1)*u_o(ii,jj,kk)  +(1.0_DBL-fexuo(ii-1))*u_o(ii-1,jj,kk)
    !
    ! v
    !
    vn = fexpo(ii)*v_o(ii+1,jj,kk)  +(1.0_DBL-fexpo(ii))  *v_o(ii,jj,kk)
    vs = fexpo(ii)*v_o(ii+1,jj-1,kk)+(1.0_DBL-fexpo(ii))  *v_o(ii,jj-1,kk)
    !
    ! w
    !
    wt = fexpo(ii)*w_o(ii+1,jj,kk)  +(1.0_DBL-fexpo(ii))  *w_o(ii,jj,kk)
    wb = fexpo(ii)*w_o(ii+1,jj,kk-1)+(1.0_DBL-fexpo(ii))  *w_o(ii,jj,kk-1)
    !
    ! gamma_n 
    !
    ! ** se utilizan las constantes gammai y gammad como auxiliares para
    ! calcular gamman, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammad = ( gamma_momento(ii+1,jj+1,kk) * gamma_momento(ii+1,jj,kk) ) / &
         &(gamma_momento(ii+1,jj+1,kk)*(1._DBL-feypo(jj))+&
         &gamma_momento(ii+1,jj,kk)*feypo(jj) )
    gammai = ( gamma_momento(ii,jj+1,kk) * gamma_momento(ii,jj,kk) ) / &
         &( gamma_momento(ii,jj+1,kk) * (1._DBL-feypo(jj))+&
         &gamma_momento(ii,jj,kk)*feypo(jj) )
    !
    gamman = gammai*gammad / (gammad * (1._DBL-fexpo(ii)) + gammai * fexpo(ii))
    !
    ! gamma_s 
    !
    ! ** se utilizan las constantes gammai y gammad como auxiliares para
    ! calcular gamman, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammad = ( gamma_momento(ii+1,jj,kk) * gamma_momento(ii+1,jj-1,kk) ) / &
         &(gamma_momento(ii+1,jj,kk)*(1._DBL-feypo(jj-1))+&
         &gamma_momento(ii+1,jj-1,kk)*feypo(jj-1))
    gammai = ( gamma_momento(ii,jj,kk) * gamma_momento(ii,jj-1,kk) ) / &
         &(gamma_momento(ii,jj,kk)*(1._DBL-feypo(jj-1))+&
         &gamma_momento(ii,jj-1,kk)*feypo(jj-1))
    !
    gammas = gammai*gammad / (gammad * (1._DBL-fexpo(ii)) + gammai * fexpo(ii))
    !
    ! gamma_t
    !
    ! ** se utilizan las constantes gammai y gammad como auxiliares para
    ! calcular gammat, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammad = ( gamma_momento(ii+1,jj,kk+1) * gamma_momento(ii+1,jj,kk) ) / &
         &( gamma_momento(ii+1,jj,kk+1)*(1._DBL-fezpo(kk))+&
         &gamma_momento(ii+1,jj,kk)*fezpo(kk) )
    gammai = ( gamma_momento(ii,jj,kk+1) * gamma_momento(ii,jj,kk) ) / &
         &( gamma_momento(ii,jj,kk+1) * (1._DBL-fezpo(kk))+&
         &gamma_momento(ii,jj,kk)*fezpo(kk) )
    !
    gammat = gammai*gammad / (gammad * (1._DBL-fexpo(ii)) + gammai * fexpo(ii))
    !
    ! gamma_b 
    !
    ! ** se utilizan las constantes gammai y gammad como auxiliares para
    ! calcular gamman, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammad = ( gamma_momento(ii+1,jj,kk) * gamma_momento(ii+1,jj,kk-1) ) / &
         &(gamma_momento(ii+1,jj,kk)*(1._DBL-fezpo(kk-1))+&
         &gamma_momento(ii+1,jj,kk-1)*fezpo(kk-1))
    gammai = ( gamma_momento(ii,jj,kk+1) * gamma_momento(ii,jj,kk) ) / &
         &(gamma_momento(ii,jj,kk)*(1._DBL-fezpo(kk-1))+&
         &gamma_momento(ii,jj,kk-1)*fezpo(kk-1))
    !
    gammab = gammai*gammad / (gammad * (1._DBL-fexpo(ii)) + gammai * fexpo(ii))
    !
    ! gamma_i
    !
    gammai = gamma_momento(ii,jj,kk)
    !
    ! gamma_d
    !
    gammad = gamma_momento(ii+1,jj,kk)
    !
    ! distancias entre nodos contiguos
    !
    di = deltaxpo(ii)
    dd = deltaxpo(ii+1)
    ds = deltayvo(jj-1)
    dn = deltayvo(jj)
    db = deltazwo(kk-1)
    dt = deltazwo(kk)
    !
    ! Tama\~no de los vol\'umenes de control para la velocidad u
    !
    deltax = deltaxuo(ii)
    deltay = deltayuo(jj)
    deltaz = deltazuo(kk)
    !
    ! Interpolaci\'on para la temperatura
    !
    temp_int = fexpo(ii)*temp_o(ii+1,jj,kk) + (1.0_DBL-fexpo(ii))*temp_o(ii,jj,kk)
    !
    ! *************************
    !
    ! Coeficientes de la matriz
    !
    alpha =-(gammai*deltay*deltaz/di*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ui*di/gammai))**5)+&
         &DMAX1(0.0_DBL, ui*deltay*deltaz))
    !
    beta =-(gammad*deltay*deltaz/dd*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ud*dd/gammad))**5)+&
         &DMAX1(0.0_DBL,-ud*deltay*deltaz))
    !
    AI_o(indeyu(jj,ii,kk)) =-(gammas*deltax*deltaz/ds*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vs*ds/gammas))**5)+&
         &DMAX1(0.0_DBL, vs*deltax*deltaz))
    !
    AD_o(indeyu(jj,ii,kk)) =-(gamman*deltax*deltaz/dn*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vn*dn/gamman))**5)+&
         &DMAX1(0.0_DBL,-vn*deltax*deltaz))
    !
    gamma =-(gammab*deltax*deltay/db*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wb*db/gammab))**5)+&
         &DMAX1(0.0_DBL, wb*deltax*deltay))
    !
    delta =-(gammat*deltax*deltay/dt*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wt*dt/gammat))**5)+&
         &DMAX1(0.0_DBL,-wt*deltax*deltay)) 
    !
    AC_o(indeyu(jj,ii,kk)) = ( -AI_o(indeyu(jj,ii,kk)) - AD_o(indeyu(jj,ii,kk)) &
         &- alpha - beta - gamma - delta - &
         &deltax*deltay*deltaz*fuente_lin_uo(ii,jj,kk)+&
         &deltax*deltay*deltaz/dt_o ) / rel_vo
    !
    Rx_o(indeyu(jj,ii,kk)) =-alpha*u_o(ii-1,jj,kk) -&
         &beta  * u_o(ii+1,jj,kk) - &
         &gamma * u_o(ii,jj,kk-1) - &
         &delta * u_o(ii,jj,kk+1) - &
         &deltax*deltay*deltaz*Ri_o(ii,jj,kk)*temp_int+&
         &deltax*deltay*deltaz*fuente_con_uo(ii,jj,kk)+&
         &deltax*deltay*deltaz*u_anto(ii,jj,kk)/dt_o+&
         &(pres_o(ii,jj,kk)-pres_o(ii+1,jj,kk))*deltay*deltaz+&
         &AC_o(indeyu(jj,ii,kk))*(1._DBL-rel_vo)*u_o(ii,jj,kk)
    !
    !au_o(ii,jj,kk) = AC_o(indexu(ii,jj,kk)) * rel_vo
    !
  end subroutine ensambla_velu_y
  !
  !-------------------------------------------------------------------
  !*******************************************************************
  !
  ! ensambla_velu_z
  !
  ! Subrutina que calcula los coeficientes de la matriz tridiagonal
  ! para la velocidad u en la direcci\'on z
  !
  ! En esta subrutina se utiliza el indice indezp para la malla de u
  !
  !*******************************************************************
  !-------------------------------------------------------------------
  !
  subroutine ensambla_velu_z(&
       &deltaxuo,&
       &deltayuo,&
       &deltazuo,&
       &deltaxpo,&
       &deltayvo,&
       &deltazwo,&
       &fexpo,&
       &feypo,&
       &fezpo,&
       &fexuo,&
       &gamma_momento,&
       &u_o,&
       &u_anto,&
       &v_o,&
       &w_o,&
       &temp_o,&
       &pres_o,&
       &fuente_con_uo,&
       &fuente_lin_uo,&
       &Ri_o,&
       &dt_o,&
       &rel_vo,&
       &AI_o,AC_o,AD_o,Rx_o,&
       &kk,jj,ii&
       &)
    implicit none
    !$acc routine
    !
    ! Tama\~no del volumen de control
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxuo
    real(kind=DBL), dimension(nj), intent(in) :: deltayuo
    real(kind=DBL), dimension(lk), intent(in) :: deltazuo
    !
    ! Distancia entre nodos contiguos de la malla de u en direcci\'on horizontal x
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxpo
    !
    ! Distancia entre nodos contiguos de la malla de u en direcci\'on horizontal y
    !
    real(kind=DBL), dimension(nj), intent(in) :: deltayvo
    !
    ! Distancia entre nodos contiguos de la malla de u en direcci\'on vertical z
    !
    real(kind=DBL), dimension(lk), intent(in) :: deltazwo
    !
    ! Coeficientes para interpolaci\'on
    !
    real(kind=DBL), DIMENSION(mi),   intent(in) :: fexpo
    real(kind=DBL), DIMENSION(nj),   intent(in) :: feypo
    real(kind=DBL), dimension(lk),   intent(in) :: fezpo
    real(kind=DBL), DIMENSION(mi-1), intent(in) :: fexuo
    !
    ! Coeficiente de difusi\'on
    !
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in) :: gamma_momento
    !
    ! Velocidad, presi\'on y temperatura
    !
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in) :: u_o, u_anto
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in) :: v_o
    real(kind=DBL), dimension(mi+1,nj+1,lk),   intent(in) :: w_o
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in) :: temp_o, pres_o
    !
    ! T\'erminos fuente
    !
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in) :: fuente_con_uo
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in) :: fuente_lin_uo    
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in) :: Ri_o
    
    !
    ! Incremento de tiempo y coeficiente de relajaci\'on
    !
    real(kind=DBL), intent(in) :: dt_o, rel_vo
    !
    ! \'Indices para recorrer las direcciones x, y, z
    !
    integer, intent(in)        :: ii, jj, kk
    !
    ! Coeficientes de las matrices
    !
    ! ** Estos coeficientes est\'an sobredimensionados para reducir el uso de memoria
    ! en la gpu, los arreglos que se reciben en esta subrutina se usan para las ecs.
    ! de momento, energ\'ia y la correcci\'on de la presi\'on **
    !
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AI_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AC_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AD_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: RX_o
    !
    ! Variables auxiliares
    !
    integer :: info
    !
    ! Auxiliares de interpolaci\'on y coeficientes
    !
    real(kind=DBL) :: ui, ud, vs, vn, wb, wt
    real(kind=DBL) :: di, dd, ds, dn, db, dt
    real(kind=DBL) :: gammai, gammad
    real(kind=DBL) :: gammas, gamman
    real(kind=DBL) :: gammab, gammat
    real(kind=DBL) :: alpha, beta, gamma, delta
    real(kind=DBL) :: deltax, deltay, deltaz
    real(kind=DBL) :: temp_int
    !
    ! Interpolaciones necesarias
    !
    ! u
    !
    ud = fexuo(ii)  *u_o(ii+1,jj,kk)+(1.0_DBL-fexuo(ii))  *u_o(ii,jj,kk)
    ui = fexuo(ii-1)*u_o(ii,jj,kk)  +(1.0_DBL-fexuo(ii-1))*u_o(ii-1,jj,kk)
    !
    ! v
    !
    vn = fexpo(ii)*v_o(ii+1,jj,kk)  +(1.0_DBL-fexpo(ii))  *v_o(ii,jj,kk)
    vs = fexpo(ii)*v_o(ii+1,jj-1,kk)+(1.0_DBL-fexpo(ii))  *v_o(ii,jj-1,kk)
    !
    ! w
    !
    wt = fexpo(ii)*w_o(ii+1,jj,kk)  +(1.0_DBL-fexpo(ii))  *w_o(ii,jj,kk)
    wb = fexpo(ii)*w_o(ii+1,jj,kk-1)+(1.0_DBL-fexpo(ii))  *w_o(ii,jj,kk-1)
    !
    ! gamma_n 
    !
    ! ** se utilizan las constantes gammai y gammad como auxiliares para
    ! calcular gamman, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammad = ( gamma_momento(ii+1,jj+1,kk) * gamma_momento(ii+1,jj,kk) ) / &
         &(gamma_momento(ii+1,jj+1,kk)*(1._DBL-feypo(jj))+&
         &gamma_momento(ii+1,jj,kk)*feypo(jj) )
    gammai = ( gamma_momento(ii,jj+1,kk) * gamma_momento(ii,jj,kk) ) / &
         &( gamma_momento(ii,jj+1,kk) * (1._DBL-feypo(jj))+&
         &gamma_momento(ii,jj,kk)*feypo(jj) )
    !
    gamman = gammai*gammad / (gammad * (1._DBL-fexpo(ii)) + gammai * fexpo(ii))
    !
    ! gamma_s 
    !
    ! ** se utilizan las constantes gammai y gammad como auxiliares para
    ! calcular gamman, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammad = ( gamma_momento(ii+1,jj,kk) * gamma_momento(ii+1,jj-1,kk) ) / &
         &(gamma_momento(ii+1,jj,kk)*(1._DBL-feypo(jj-1))+&
         &gamma_momento(ii+1,jj-1,kk)*feypo(jj-1))
    gammai = ( gamma_momento(ii,jj,kk) * gamma_momento(ii,jj-1,kk) ) / &
         &(gamma_momento(ii,jj,kk)*(1._DBL-feypo(jj-1))+&
         &gamma_momento(ii,jj-1,kk)*feypo(jj-1))
    !
    gammas = gammai*gammad / (gammad * (1._DBL-fexpo(ii)) + gammai * fexpo(ii))
    !
    ! gamma_t
    !
    ! ** se utilizan las constantes gammai y gammad como auxiliares para
    ! calcular gammat, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammad = ( gamma_momento(ii+1,jj,kk+1) * gamma_momento(ii+1,jj,kk) ) / &
         &( gamma_momento(ii+1,jj,kk+1)*(1._DBL-fezpo(kk))+&
         &gamma_momento(ii+1,jj,kk)*fezpo(kk) )
    gammai = ( gamma_momento(ii,jj,kk+1) * gamma_momento(ii,jj,kk) ) / &
         &( gamma_momento(ii,jj,kk+1) * (1._DBL-fezpo(kk))+&
         &gamma_momento(ii,jj,kk)*fezpo(kk) )
    !
    gammat = gammai*gammad / (gammad * (1._DBL-fexpo(ii)) + gammai * fexpo(ii))
    !
    ! gamma_b 
    !
    ! ** se utilizan las constantes gammai y gammad como auxiliares para
    ! calcular gamman, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammad = ( gamma_momento(ii+1,jj,kk) * gamma_momento(ii+1,jj,kk-1) ) / &
         &(gamma_momento(ii+1,jj,kk)*(1._DBL-fezpo(kk-1))+&
         &gamma_momento(ii+1,jj,kk-1)*fezpo(kk-1))
    gammai = ( gamma_momento(ii,jj,kk+1) * gamma_momento(ii,jj,kk) ) / &
         &(gamma_momento(ii,jj,kk)*(1._DBL-fezpo(kk-1))+&
         &gamma_momento(ii,jj,kk-1)*fezpo(kk-1))
    !
    gammab = gammai*gammad / (gammad * (1._DBL-fexpo(ii)) + gammai * fexpo(ii))
    !
    ! gamma_i
    !
    gammai = gamma_momento(ii,jj,kk)
    !
    ! gamma_d
    !
    gammad = gamma_momento(ii+1,jj,kk)
    !
    ! distancias entre nodos contiguos
    !
    di = deltaxpo(ii)
    dd = deltaxpo(ii+1)
    ds = deltayvo(jj-1)
    dn = deltayvo(jj)
    db = deltazwo(kk-1)
    dt = deltazwo(kk)
    !
    ! Tama\~no de los vol\'umenes de control para la velocidad u
    !
    deltax = deltaxuo(ii)
    deltay = deltayuo(jj)
    deltaz = deltazuo(kk)
    !
    ! Interpolaci\'on para la temperatura
    !
    temp_int = fexpo(ii)*temp_o(ii+1,jj,kk) + (1.0_DBL-fexpo(ii))*temp_o(ii,jj,kk)
    !
    ! *************************
    !
    ! Coeficientes de la matriz
    !
    alpha =-(gammai*deltay*deltaz/di*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ui*di/gammai))**5)+&
         &DMAX1(0.0_DBL, ui*deltay*deltaz))
    !
    beta  =-(gammad*deltay*deltaz/dd*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ud*dd/gammad))**5)+&
         &DMAX1(0.0_DBL,-ud*deltay*deltaz))
    !
    gamma =-(gammas*deltax*deltaz/ds*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vs*ds/gammas))**5)+&
         &DMAX1(0.0_DBL, vs*deltax*deltaz))
    !
    delta =-(gamman*deltax*deltaz/dn*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vn*dn/gamman))**5)+&
         &DMAX1(0.0_DBL,-vn*deltax*deltaz))
    !
    AI_o(indezp(kk,jj,ii)) =-(gammab*deltax*deltay/db*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wb*db/gammab))**5)+&
         &DMAX1(0.0_DBL, wb*deltax*deltay))
    !
    AD_o(indezp(kk,jj,ii)) =-(gammat*deltax*deltay/dt*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wt*dt/gammat))**5)+&
         &DMAX1(0.0_DBL,-wt*deltax*deltay)) 
    !
    AC_o(indezp(kk,jj,ii)) = ( -AI_o(indezp(kk,jj,ii)) - AD_o(indezp(kk,jj,ii)) &
         &- alpha - beta - gamma - delta - &
         &deltax*deltay*deltaz*fuente_lin_uo(ii,jj,kk)+&
         &deltax*deltay*deltaz/dt_o ) / rel_vo
    !
    Rx_o(indezp(kk,jj,ii)) =-alpha*u_o(ii-1,jj,kk) -&
         &beta  * u_o(ii+1,jj,kk) - &
         &gamma * u_o(ii,jj-1,kk) - &
         &delta * u_o(ii,jj+1,kk) - &
         &deltax*deltay*deltaz*Ri_o(ii,jj,kk)*temp_int+&
         &deltax*deltay*deltaz*fuente_con_uo(ii,jj,kk)+&
         &deltax*deltay*deltaz*u_anto(ii,jj,kk)/dt_o+&
         &(pres_o(ii,jj,kk)-pres_o(ii+1,jj,kk))*deltay*deltaz+&
         &AC_o(indezp(kk,jj,ii))*(1._DBL-rel_vo)*u_o(ii,jj,kk)
    !
    !au_o(ii,jj,kk) = AC_o(indexu(ii,jj,kk)) * rel_vo
    !
  end subroutine ensambla_velu_z
  !
  !*******************************************************************
  !
  ! residuo_u
  !
  ! Subrutina que calcula el residuo de la ecuaci\'on de momento en u
  !
  !*******************************************************************
  subroutine residuo_u(&
       &deltaxuo,&
       &deltayuo,&
       &deltazuo,&
       &deltaxpo,&
       &deltayvo,&
       &deltazwo,&
       &fexpo,&
       &feypo,&
       &fezpo,&
       &fexuo,&
       &gamma_momento,&
       &u_o,&
       &u_anto,&
       &v_o,&
       &w_o,&
       &temp_o,&
       &pres_o,&
       &fuente_con_uo,&
       &fuente_lin_uo,&
       &Ri_o,&
       &dt_o,&
       &rel_vo,&
       &AI_o,AC_o,AD_o,Rx_o,&
       &ii,jj,kk&
       &)
    implicit none
    !$acc routine gang
        !
    ! Tama\~no del volumen de control
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxuo
    real(kind=DBL), dimension(nj), intent(in) :: deltayuo
    real(kind=DBL), dimension(lk), intent(in) :: deltazuo
    !
    ! Distancia entre nodos contiguos de la malla de u en direcci\'on horizontal x
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxpo
    !
    ! Distancia entre nodos contiguos de la malla de u en direcci\'on horizontal y
    !
    real(kind=DBL), dimension(nj), intent(in) :: deltayvo
    !
    ! Distancia entre nodos contiguos de la malla de u en direcci\'on vertical z
    !
    real(kind=DBL), dimension(lk), intent(in) :: deltazwo
    !
    ! Coeficientes para interpolaci\'on
    !
    real(kind=DBL), DIMENSION(mi),   intent(in) :: fexpo
    real(kind=DBL), DIMENSION(nj),   intent(in) :: feypo
    real(kind=DBL), dimension(lk),   intent(in) :: fezpo
    real(kind=DBL), DIMENSION(mi-1), intent(in) :: fexuo
    !
    ! Coeficiente de difusi\'on
    !
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in) :: gamma_momento
    !
    ! Velocidad, presi\'on y temperatura
    !
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in) :: u_o, u_anto
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in) :: v_o
    real(kind=DBL), dimension(mi+1,nj+1,lk),   intent(in) :: w_o
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in) :: temp_o, pres_o
    !
    ! T\'erminos fuente
    !
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in) :: fuente_con_uo
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in) :: fuente_lin_uo    
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in) :: Ri_o
    
    !
    ! Incremento de tiempo y coeficiente de relajaci\'on
    !
    real(kind=DBL), intent(in) :: dt_o, rel_vo
    !
    ! Variables para bucles
    !
    integer, intent(in) :: ii,jj,kk
    !
    ! Coeficientes de las matrices
    !
    ! ** Estos coeficientes est\'an sobredimensionados para reducir el uso de memoria
    ! en la gpu, los arreglos que se reciben en esta subrutina se usan para las ecs.
    ! de momento, energ\'ia y la correcci\'on de la presi\'on **
    !
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AI_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AC_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AD_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: RX_o
    !
    ! Auxiliares de interpolaci\'on y coeficientes
    !
    real(kind=DBL) :: ui, ud, vs, vn, wb, wt
    real(kind=DBL) :: di, dd, ds, dn, db, dt
    real(kind=DBL) :: gammai, gammad
    real(kind=DBL) :: gammas, gamman
    real(kind=DBL) :: gammab, gammat
    real(kind=DBL) :: alpha, beta, gamma, delta
    real(kind=DBL) :: deltax, deltay, deltaz
    real(kind=DBL) :: temp_int
    !
    !
    ! Interpolaciones necesarias
    !
    ! u
    !
    ud = fexuo(ii)  *u_o(ii+1,jj,kk)+(1.0_DBL-fexuo(ii))  *u_o(ii,jj,kk)
    ui = fexuo(ii-1)*u_o(ii,jj,kk)  +(1.0_DBL-fexuo(ii-1))*u_o(ii-1,jj,kk)
    !
    ! v
    !
    vn = fexpo(ii)*v_o(ii+1,jj,kk)  +(1.0_DBL-fexpo(ii))  *v_o(ii,jj,kk)
    vs = fexpo(ii)*v_o(ii+1,jj-1,kk)+(1.0_DBL-fexpo(ii))  *v_o(ii,jj-1,kk)
    !
    ! w
    !
    wt = fexpo(ii)*w_o(ii+1,jj,kk)  +(1.0_DBL-fexpo(ii))  *w_o(ii,jj,kk)
    wb = fexpo(ii)*w_o(ii+1,jj,kk-1)+(1.0_DBL-fexpo(ii))  *w_o(ii,jj,kk-1)
    !
    ! gamma_n 
    !
    ! ** se utilizan las constantes gammai y gammad como auxiliares para
    ! calcular gamman, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammad = ( gamma_momento(ii+1,jj+1,kk) * gamma_momento(ii+1,jj,kk) ) / &
         &(gamma_momento(ii+1,jj+1,kk)*(1._DBL-feypo(jj))+&
         &gamma_momento(ii+1,jj,kk)*feypo(jj) )
    gammai = ( gamma_momento(ii,jj+1,kk) * gamma_momento(ii,jj,kk) ) / &
         &( gamma_momento(ii,jj+1,kk) * (1._DBL-feypo(jj))+&
         &gamma_momento(ii,jj,kk)*feypo(jj) )
    !
    gamman = gammai*gammad / (gammad * (1._DBL-fexpo(ii)) + gammai * fexpo(ii))
    !
    ! gamma_s 
    !
    ! ** se utilizan las constantes gammai y gammad como auxiliares para
    ! calcular gamman, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammad = ( gamma_momento(ii+1,jj,kk) * gamma_momento(ii+1,jj-1,kk) ) / &
         &(gamma_momento(ii+1,jj,kk)*(1._DBL-feypo(jj-1))+&
         &gamma_momento(ii+1,jj-1,kk)*feypo(jj-1))
    gammai = ( gamma_momento(ii,jj,kk) * gamma_momento(ii,jj-1,kk) ) / &
         &(gamma_momento(ii,jj,kk)*(1._DBL-feypo(jj-1))+&
         &gamma_momento(ii,jj-1,kk)*feypo(jj-1))
    !
    gammas = gammai*gammad / (gammad * (1._DBL-fexpo(ii)) + gammai * fexpo(ii))
    !
    ! gamma_t
    !
    ! ** se utilizan las constantes gammai y gammad como auxiliares para
    ! calcular gammat, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammad = ( gamma_momento(ii+1,jj,kk+1) * gamma_momento(ii+1,jj,kk) ) / &
         &( gamma_momento(ii+1,jj,kk+1)*(1._DBL-fezpo(kk))+&
         &gamma_momento(ii+1,jj,kk)*fezpo(kk) )
    gammai = ( gamma_momento(ii,jj,kk+1) * gamma_momento(ii,jj,kk) ) / &
         &( gamma_momento(ii,jj,kk+1) * (1._DBL-fezpo(kk))+&
         &gamma_momento(ii,jj,kk)*fezpo(kk) )
    !
    gammat = gammai*gammad / (gammad * (1._DBL-fexpo(ii)) + gammai * fexpo(ii))
    !
    ! gamma_b 
    !
    ! ** se utilizan las constantes gammai y gammad como auxiliares para
    ! calcular gamman, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammad = ( gamma_momento(ii+1,jj,kk) * gamma_momento(ii+1,jj,kk-1) ) / &
         &(gamma_momento(ii+1,jj,kk)*(1._DBL-fezpo(kk-1))+&
         &gamma_momento(ii+1,jj,kk-1)*fezpo(kk-1))
    gammai = ( gamma_momento(ii,jj,kk+1) * gamma_momento(ii,jj,kk) ) / &
         &(gamma_momento(ii,jj,kk)*(1._DBL-fezpo(kk-1))+&
         &gamma_momento(ii,jj,kk-1)*fezpo(kk-1))
    !
    gammab = gammai*gammad / (gammad * (1._DBL-fexpo(ii)) + gammai * fexpo(ii))
    !
    ! gamma_i
    !
    gammai = gamma_momento(ii,jj,kk)
    !
    ! gamma_d
    !
    gammad = gamma_momento(ii+1,jj,kk)
    !
    ! distancias entre nodos contiguos
    !
    di = deltaxpo(ii)
    dd = deltaxpo(ii+1)
    ds = deltayvo(jj-1)
    dn = deltayvo(jj)
    db = deltazwo(kk-1)
    dt = deltazwo(kk)
    !
    ! Tama\~no de los vol\'umenes de control para la velocidad u
    !
    deltax = deltaxuo(ii)
    deltay = deltayuo(jj)
    deltaz = deltazuo(kk)
    !
    ! Interpolaci\'on para la temperatura
    !
    temp_int = fexpo(ii)*temp_o(ii+1,jj,kk) + (1.0_DBL-fexpo(ii))*temp_o(ii,jj,kk)
    !
    ! *************************
    !
    ! Coeficientes de la matriz
    !
    AI_o(indexu(ii,jj,kk)) =-(gammai*deltay*deltaz/di*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ui*di/gammai))**5)+&
         &DMAX1(0.0_DBL, ui*deltay*deltaz))
    !
    AD_o(indexu(ii,jj,kk)) =-(gammad*deltay*deltaz/dd*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ud*dd/gammad))**5)+&
         &DMAX1(0.0_DBL,-ud*deltay*deltaz))
    !
    alpha =-(gammas*deltax*deltaz/ds*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vs*ds/gammas))**5)+&
         &DMAX1(0.0_DBL, vs*deltax*deltaz))
    !
    beta  =-(gamman*deltax*deltaz/dn*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vn*dn/gamman))**5)+&
         &DMAX1(0.0_DBL,-vn*deltax*deltaz))
    !
    gamma =-(gammab*deltax*deltay/db*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wb*db/gammab))**5)+&
         &DMAX1(0.0_DBL, wb*deltax*deltay))
    !
    delta =-(gammat*deltax*deltay/dt*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(wt*dt/gammat))**5)+&
         &DMAX1(0.0_DBL,-wt*deltax*deltay)) 
    !
    AC_o(indexu(ii,jj,kk)) = -AI_o(indexu(ii,jj,kk)) - AD_o(indexu(ii,jj,kk)) &
         &- alpha - beta - gamma - delta - &
         &deltax*deltay*deltaz*fuente_lin_uo(ii,jj,kk)+&
         &deltax*deltay*deltaz/dt_o 
    !
    Rx_o(indexu(ii,jj,kk)) =-AI_o(indexu(ii,jj,kk)) * u_o(ii-1,jj,kk) - &
         & AC_o(indexu(ii,jj,kk)) * u_o(ii,jj,kk)   - &
         & AD_o(indexu(ii,jj,kk)) * u_o(ii+1,jj,kk) - &
         & alpha * u_o(ii,jj-1,kk) - &
         & beta  * u_o(ii,jj+1,kk) - &
         & gamma * u_o(ii,jj,kk-1) - &
         & delta * u_o(ii,jj,kk+1) - &
         & deltax*deltay*deltaz*Ri_o(ii,jj,kk)*temp_int + &
         & deltax*deltay*deltaz*fuente_con_uo(ii,jj,kk) + &
         & deltax*deltay*deltaz*u_anto(ii,jj,kk)/dt_o   + &
         & (pres_o(ii,jj,kk)-pres_o(ii+1,jj,kk))*deltay*deltaz
    !
  end subroutine residuo_u
  !
end module ec_momento
