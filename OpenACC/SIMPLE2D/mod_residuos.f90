!
! M\'odulo de residuos
!
! Este m\'odulo contiene las subrutinas que calculan residuos 
! de ecuaciones de movimiento
!
!
module residuos
  use malla, only :  mi, nj, DBL
  implicit none
  
contains
  !*******************************************************************
  !
  ! residuo_u
  !
  ! Subrutina que calcula el residuo de la ecuaci\'on de momento en u
  !
  !*******************************************************************
  subroutine residuo_u(deltaxuo,deltayuo,deltaxpo,&
       &deltayvo,fexpo,feypo,fexuo,gamma_momento,&
       &u_o,u_anto,v_o,&
       &temp_o,pres_o,Ri_o,dt_o,rel_vo,&
       &Rx_o&
       &)
    implicit none
    !$acc routine gang
    !
    ! Tama\~no del volumen de control
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxuo
    real(kind=DBL), dimension(nj), intent(in) :: deltayuo
    !
    ! Distancia entre nodos contiguos de la malla de u en direcci\'on horizontal
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxpo
    !
    ! Distancia entre nodos contiguos de la malla de u en direcci\'on vertical
    !
    real(kind=DBL), dimension(nj), intent(in) :: deltayvo
    !
    ! Coeficientes para interpolaci\'on
    !
    real(kind=DBL), DIMENSION(mi),   intent(in) :: fexpo
    real(kind=DBL), DIMENSION(nj),   intent(in) :: feypo
    real(kind=DBL), DIMENSION(mi-1), intent(in) :: fexuo
    !
    ! Coeficiente de difusi\'on
    !
    real(kind=DBL), dimension(mi+1,nj+1), intent(in) :: gamma_momento
    !
    ! Velocidad, presi\'on y temperatura
    !
    real(kind=DBL), dimension(mi,nj+1),   intent(in) :: u_o,    u_anto
    real(kind=DBL), dimension(mi+1,nj),   intent(in) :: v_o
    real(kind=DBL), dimension(mi+1,nj+1), intent(in) :: temp_o, pres_o
    !
    ! T\'erminos fuente
    !
    real(kind=DBL), dimension(mi,nj+1),   intent(in) :: Ri_o
    !
    ! Incremento de tiempo y coeficiente de relajaci\'on
    !
    real(kind=DBL), intent(in) :: dt_o, rel_vo
    !
    ! Residuo de la ecuaci\'on
    !
    real(kind=DBL), dimension(mi,nj+1), intent(out)   :: Rx_o
    !
    ! Variables auxiliares
    !
    integer :: ii, jj, info
    !
    ! Auxiliares de interpolaci\'on
    !
    real(kind=DBL) :: ui, ud, vs, vn
    real(kind=DBL) :: di, dd, ds, dn
    real(kind=DBL) :: gammai, gammad
    real(kind=DBL) :: gammas, gamman
    real(kind=DBL) :: temp_int
    !
    ! Coeficientes auxiliares
    !
    real(kind=DBL) :: AI_o, AC_o, AD_o
    real(kind=DBL) :: BS_o, BN_o
    !
    !$acc loop gang
    bucle_direccion_y: do jj = 2, nj
       !
       ! Llenado de la matriz
       !
       !$acc loop vector
       bucle_direccion_x: do ii = 2, mi-1
          !
          ! Interpolaciones necesarias
          !
          ! u
          !
          ud = fexuo(ii)  *u_o(ii+1,jj)+(1.0_DBL-fexuo(ii))  *u_o(ii,jj)
          ui = fexuo(ii-1)*u_o(ii,jj)  +(1.0_DBL-fexuo(ii-1))*u_o(ii-1,jj)
          !
          ! v
          !
          vn = fexpo(ii)*v_o(ii+1,jj)  +(1.0_DBL-fexpo(ii))  *v_o(ii,jj)
          vs = fexpo(ii)*v_o(ii+1,jj-1)+(1.0_DBL-fexpo(ii))  *v_o(ii,jj-1)
          !
          ! gamma_n 
          !
          ! ** se utilizan las constantes gammai y gammad como auxiliares para
          ! calcular gamman, despu\'es se utilizan para el coeficiente gamma que
          ! corresponde **
          gammad = ( gamma_momento(ii+1,jj+1) * gamma_momento(ii+1,jj) ) / &
               &( gamma_momento(ii+1,jj+1) * (1._DBL-feypo(jj))+gamma_momento(ii+1,jj)*feypo(jj) )
          gammai = ( gamma_momento(ii,jj+1) * gamma_momento(ii,jj) ) / &
               &( gamma_momento(ii,jj+1) * (1._DBL-feypo(jj))+gamma_momento(ii,jj)*feypo(jj) )
          !
          gamman = gammai*gammad / (gammad * (1._DBL-fexpo(ii)) + gammai * fexpo(ii))
          !
          ! gamma_s 
          !
          ! ** se utilizan las constantes gammai y gammad como auxiliares para
          ! calcular gamman, despu\'es se utilizan para el coeficiente gamma que
          ! corresponde **
          gammad = ( gamma_momento(ii+1,jj) * gamma_momento(ii+1,jj-1) ) / &
               &(gamma_momento(ii+1,jj)*(1._DBL-feypo(jj-1))+gamma_momento(ii+1,jj-1)*feypo(jj-1))
          gammai = ( gamma_momento(ii,jj+1) * gamma_momento(ii,jj) ) / &
               &(gamma_momento(ii,jj)*(1._DBL-feypo(jj-1))+gamma_momento(ii,jj-1)*feypo(jj-1))
          !
          gammas = gammai*gammad / (gammad * (1._DBL-fexpo(ii)) + gammai * fexpo(ii))         
          !
          ! gamma_i
          !
          gammai = gamma_momento(ii,jj)
          !
          ! gamma_d
          !
          gammai = gamma_momento(ii+1,jj)
          !
          ! distancias entre nodos contiguos
          !
          di = deltaxpo(ii)
          dd = deltaxpo(ii+1)
          ds = deltayvo(jj-1)
          dn = deltayvo(jj)
          !
          ! Tama\~no de los vol\'umenes de control para la velocidad u
          !
          ! delta_x = deltaxuo(ii)
          ! delta_y = deltayuo(jj)
          !
          ! Interpolaci\'on para la temperatura
          !
          temp_int = fexpo(ii)*temp_o(ii+1,jj) + (1.0_DBL-fexpo(ii))*temp_o(ii,jj)
          !
          ! temp_int = (1._DBL-di/(2._DBL*delta_x))*temp_o(i,j)+di/(2._DBL*delta_x)*temp_o(i+1,j)
          ! !     (temp_o(i,j)+temp_o(i+1,j))/2._DBL
          !
          ! *************************
          !
          ! Coeficientes de la matriz
          !
          AI_o = (gammai*deltayuo(jj)/di*&
               &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ui*di/gammai))**5)+&
               &DMAX1(0.0_DBL,ui*deltayuo(jj)))
          AD_o = (gammad*deltayuo(jj)/dd*&
               &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ud*dd/gammad))**5)+&
               &DMAX1(0.0_DBL,-ud*deltayuo(jj)))
          BS_o = (gammas*deltaxuo(ii)/ds*&
               &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vs*ds/gammas))**5)+&
               &DMAX1(0.0_DBL, vs*deltaxuo(ii)))
          BN_o = (gamman*deltaxuo(ii)/dn*&
               &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vn*dn/gamman))**5)+&
               &DMAX1(0.0_DBL,-vn*deltaxuo(ii)))
          AC_o = AI_o + AD_o + BS_o + BN_o+&
               &deltaxuo(ii)*deltayuo(jj)/dt_o
          Rx_o(ii,jj) = AI_o*u_o(ii-1,jj)-AC_o*u_o(ii,jj)+AD_o*u_o(ii+1,jj)+&
               &BS_o*u_o(ii,jj-1)+BN_o*u_o(ii,jj+1)-&
               &deltaxuo(ii)*deltayuo(jj)*Ri_o(ii,jj)*temp_int+&
               &deltaxuo(ii)*deltayuo(jj)*u_anto(ii,jj)/dt_o+&
               &(pres_o(ii,jj)-pres_o(ii+1,jj))*deltayuo(jj)
       end do bucle_direccion_x
    end do bucle_direccion_y
    !
   end subroutine residuo_u

end module residuos
