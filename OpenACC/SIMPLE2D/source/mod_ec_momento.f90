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
  use malla,      only : mi, nj, DBL
  use malla,      only : xu, yv, xp, yp
  !
  use cond_frontera, only : tipo_cond_front
  !
  use cond_frontera, only : inicializa_cond_front
  use cond_frontera, only : lectura_cond_frontera
  !
  implicit none
  !
  ! Componentes de velocidad, presi\'on
  ! residuos de las ecuaciones de momento y correcci\'on de la presi\'on
  ! coeficientes de difusi\'on y criterios e convergencia
  !
  REAL(kind=DBL), DIMENSION(mi,nj+1)   :: u, u_ant, du, au, Resu, fu, Ri
  REAL(kind=DBL), DIMENSION(mi+1,nj)   :: v, v_ant, dv, av, Resv, fv
  !
  ! Variables para los t\'erminos fuente de la ec de momento
  ! El t\'ermino lineal fuente_lin debe ser negativo para favorecer
  ! convergencia y estabilidad
  !
  real(kind=DBL), dimension(mi,nj+1)   :: fuente_con_u, fuente_lin_u
  real(kind=DBL), dimension(mi+1,nj)   :: fuente_con_v, fuente_lin_v
  !
  ! Variables para la presi\'on
  !
  REAL(kind=DBL), DIMENSION(mi+1,nj+1) :: pres, corr_pres, dcorr_pres, fcorr_pres
  REAL(kind=DBL), DIMENSION(mi+1,nj+1) :: uf, vf, b_o
  real(kind=DBL), dimension(mi+1,nj+1) :: gamma_momen, Riy
  REAL(kind=DBL) :: maxbo, conv_u, conv_p, conv_paso, rel_pres, rel_vel
  !
  ! Estructuras para guardar la informaci\'on de las condiciones de frontera
  !
  type( tipo_cond_front ) :: cond_front_ua, cond_front_ub, cond_front_uc, cond_front_ud
  type( tipo_cond_front ) :: cond_front_va, cond_front_vb, cond_front_vc, cond_front_vd
  !
contains
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
  subroutine ini_frontera_uv()
    !
    implicit none
    !
    ! arreglos de u
    !
    call inicializa_cond_front(cond_front_ua)
    call inicializa_cond_front(cond_front_ub)
    call inicializa_cond_front(cond_front_uc)
    call inicializa_cond_front(cond_front_ud)
    call lectura_cond_frontera('cond_fronterau.dat',&
         & cond_front_ua, &
         & cond_front_ub, &
         & cond_front_uc, &
         & cond_front_ud, &
         & xu, yp,        &
         & mi, nj+1       &
         & )
    !
    ! arreglos de v
    !
    call inicializa_cond_front(cond_front_va)
    call inicializa_cond_front(cond_front_vb)
    call inicializa_cond_front(cond_front_vc)
    call inicializa_cond_front(cond_front_vd)
    !
    call lectura_cond_frontera('cond_fronterav.dat',&
         & cond_front_va, &
         & cond_front_vb, &
         & cond_front_vc, &
         & cond_front_vd, &
         & xp, yv,        &
         & mi+1, nj       &
         & )
    !
  end subroutine ini_frontera_uv
  !
  !*******************************************************************
  !
  ! ensambla_corr_pres
  !
  ! Subrutina que calcula los coeficientes de la matriz tridiagonal
  ! para la correcci\'on de la presi\'on
  !
  !*******************************************************************
  subroutine ensambla_corr_pres(deltaxpo,deltaypo,&
       &deltaxuo,deltayvo,&
       &u_o,v_o,b_o,&
       &corr_pres_o,rel_vo,&
       &AI_o,AC_o,AD_o,Rx_o,BS_o,BC_o,BN_o,Ry_o,au_o,av_o,&
       &ii,jj)
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
    ! Velocidad, presi\'on, t\'ermino fuente b
    !
    real(kind=DBL), dimension(mi,nj+1),   intent(in)  :: u_o
    real(kind=DBL), dimension(mi+1,nj),   intent(in)  :: v_o
    real(kind=DBL), dimension(mi+1,nj+1), intent(in)  :: corr_pres_o
    real(kind=DBL), dimension(mi+1,nj+1), intent(out) :: b_o
    !
    ! coeficiente de relajaci\'on
    !
    real(kind=DBL), intent(in) :: rel_vo
    !
    ! Coeficientes de las matrices
    !
    ! ** Estos coeficientes est\'an sobredimensionados para reducir el uso de memoria
    ! en la gpu, los arreglos que se reciben en esta subrutina se usan para las ecs.
    ! de momento en, energ\'ia y la correcci\'on de la presi\'on **
    !
    real(kind=DBL), dimension(mi+1,nj+1), intent(out) :: AI_o, AC_o, AD_o, Rx_o
    real(kind=DBL), dimension(nj+1,mi+1), intent(out) :: BS_o, BC_o, BN_o, Ry_o
    real(kind=DBL), dimension(mi,nj+1),   intent(in)  :: au_o
    real(kind=DBL), dimension(mi+1,nj),   intent(in)  :: av_o
    !
    ! \'Indice para recorrer la direcci\'on y
    !
    integer, intent(in)        :: ii, jj
    !
    ! Variables auxiliares
    !
    ! integer :: ii, jj
    !
    ! Auxiliares de interpolaci\'on
    !
    real(kind=DBL) :: ui, ud, vs, vn
    real(kind=DBL) :: ai, ad, as, an
    ! real(kind=DBL) :: di, dd, ds, dn
    !
    ! C\'alculo de los coeficientes
    !
    ! $acc loop gang
    ! bucle_direccion_y: do jj = 2, nj
    !------------------------
    ! Condiciones de frontera
    ! AC_o(1,jj) = 1._DBL
    ! AD_o(1,jj) = 0._DBL
    ! Rx_o(1,jj) = 0._DBL
    !
    ! Llenado de la matriz
    !
    ! $acc loop vector
    ! bucle_direccion_x: do ii = 2, mi
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
    ! coeficientes de la ecuaci\'on de momento
    !
    ai = au_o(ii-1,jj)
    ad = au_o(ii,jj)
    as = av_o(ii,jj-1)
    an = av_o(ii,jj)
    !
    ! Tama\~no de los vol\'umenes de control para la velocidad u
    !
    ! delta_x = deltaxpo(ii)
    ! delta_y = deltaypo(jj)
    !
    ! -------------------------
    !
    ! Coeficientes de la matriz
    !
    AI_o(ii,jj) =-deltaypo(jj)*deltaypo(jj)/ai
    AD_o(ii,jj) =-deltaypo(jj)*deltaypo(jj)/ad
    BS_o(jj,ii) =-deltaxpo(ii)*deltaxpo(ii)/as
    BN_o(jj,ii) =-deltaxpo(ii)*deltaxpo(ii)/an
    AC_o(ii,jj) = ( -AI_o(ii,jj) - AD_o(ii,jj)-&
         &BS_o(jj,ii) - BN_o(jj,ii) ) / rel_vo
    BC_o(jj,ii) = AC_o(ii,jj)
    !
    b_o(ii,jj)  = (ui-ud)*deltaypo(jj)+(vs-vn)*deltaxpo(ii)
    !
    Rx_o(ii,jj) =-BS_o(jj,ii)*corr_pres_o(ii,jj-1)-&
         &BN_o(jj,ii)*corr_pres_o(ii,jj+1)+&
         b_o(ii,jj) + (1._DBL-rel_vo)*AC_o(ii,jj)*corr_pres_o(ii,jj)
    Ry_o(jj,ii) =-AI_o(ii,jj)*corr_pres_o(ii-1,jj)-&
         &AD_o(ii,jj)*corr_pres_o(ii+1,jj)+&
         b_o(ii,jj) + (1._DBL-rel_vo)*BC_o(jj,ii)*corr_pres_o(ii,jj)

    ! end do bucle_direccion_x
    ! !
    ! ! Condicion frontera
    ! !
    ! AI_o(mi+1,jj) = 0.0_DBL
    ! AC_o(mi+1,jj) = 1.0_DBL
    ! Rx_o(mi+1,jj) = 0.0_DBL

    ! end do bucle_direccion_y
    !
    ! Condiciones de frontera para la direcci\'on y
    !
    ! bucle_direccionx: do ii = 2, mi
    !    !***********************
    !    !Condiciones de frontera
    !    BC_o(1,ii)     = 1._DBL
    !    BN_o(1,ii)     = 0.0_DBL
    !    Ry_o(1,ii)     = 0.0_DBL
    !    !
    !    BC_o(nj+1,ii)  = 1._DBL
    !    BS_o(nj+1,ii)  = 0.0_DBL
    !    Ry_o(nj+1,ii)  = 0.0_DBL
    ! end do bucle_direccionx
  end subroutine ensambla_corr_pres
  !  
  !*******************************************************************
  !
  ! ensambla_velu
  !
  ! Subrutina que calcula los coeficientes de la matriz tridiagonal
  ! para la velocidad u
  !
  !*******************************************************************
  subroutine ensambla_velu(deltaxuo,deltayuo,deltaxpo,&
       &deltayvo,fexpo,feypo,fexuo,gamma_momento,&
       &fuente_con_uo,fuente_lin_uo,&
       &u_o,u_anto,v_o,&
       &temp_o,pres_o,Ri_o,dt_o,rel_vo,&
       &AI_o,AC_o,AD_o,Rx_o,BS_o,BC_o,BN_o,Ry_o,au_o,&
       &ii,jj&
       &)
    implicit none
    !$acc routine
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
    real(kind=DBL), dimension(mi,nj+1),   intent(in) :: fuente_con_uo
    real(kind=DBL), dimension(mi,nj+1),   intent(in) :: fuente_lin_uo    
    real(kind=DBL), dimension(mi,nj+1),   intent(in) :: Ri_o
    
    !
    ! Incremento de tiempo y coeficiente de relajaci\'on
    !
    real(kind=DBL), intent(in) :: dt_o, rel_vo
    ! !
    ! ! \'Indice para recorrer la direcci\'on y
    ! !
    integer, intent(in)        :: ii, jj
    !
    ! Coeficientes de las matrices
    !
    ! ** Estos coeficientes est\'an sobredimensionados para reducir el uso de memoria
    ! en la gpu, los arreglos que se reciben en esta subrutina se usan para las ecs.
    ! de momento, energ\'ia y la correcci\'on de la presi\'on **
    !
    real(kind=DBL), dimension(mi+1,nj+1), intent(out) :: AI_o, AC_o, AD_o, Rx_o
    real(kind=DBL), dimension(nj+1,mi+1), intent(out) :: BS_o, BC_o, BN_o, Ry_o
    real(kind=DBL), dimension(mi,nj+1),   intent(out) :: au_o
    !
    ! Variables auxiliares
    !
    integer :: kk, info
    !
    ! Auxiliares de interpolaci\'on
    !
    real(kind=DBL) :: ui, ud, vs, vn
    real(kind=DBL) :: di, dd, ds, dn
    real(kind=DBL) :: gammai, gammad
    real(kind=DBL) :: gammas, gamman
    real(kind=DBL) :: deltax, deltay
    real(kind=DBL) :: temp_int
    !
    ! Interpolaciones necesarias
    !
    ! u
    !
    ud = fexuo(ii)  *u_o(ii+1,jj)+(1.0_DBL-fexuo(ii))  *u_o(ii,jj)
    ui = fexuo(ii-1)*u_o(ii,jj)  +(1.0_DBL-fexuo(ii-1))*u_o(ii-1,jj)
    ! ui = (u_o(ii-1,jj)+u_o(ii,jj))/2._DBL
    ! ud = (u_o(ii,jj)+u_o(ii+1,jj))/2._DBL
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
    deltax = deltaxuo(ii)
    deltay = deltayuo(jj)
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
    AI_o(ii,jj) =-(gammai*deltay/di*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ui*di/gammai))**5)+&
         &DMAX1(0.0_DBL,ui*deltay))
    
    AD_o(ii,jj) =-(gammad*deltay/dd*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ud*dd/gammad))**5)+&
         &DMAX1(0.0_DBL,-ud*deltay))

    BS_o(jj,ii) =-(gammas*deltax/ds*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vs*ds/gammas))**5)+&
         &DMAX1(0.0_DBL, vs*deltax))

    BN_o(jj,ii) =-(gamman*deltax/dn*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vn*dn/gamman))**5)+&
         &DMAX1(0.0_DBL,-vn*deltax))

    AC_o(ii,jj) = ( -AI_o(ii,jj) - AD_o(ii,jj) - BS_o(jj,ii) - BN_o(jj,ii) -&
         &deltax*deltay*fuente_lin_uo(ii,jj)+&
         &deltax*deltay/dt_o) / rel_vo

    Rx_o(ii,jj) =-BS_o(jj,ii)*u_o(ii,jj-1) - BN_o(jj,ii)*u_o(ii,jj+1)-&
         &deltax*deltay*Ri_o(ii,jj)*temp_int+&
         &deltax*deltay*fuente_con_uo(ii,jj)+&
         &deltax*deltay*u_anto(ii,jj)/dt_o+&
         &(pres_o(ii,jj)-pres_o(ii+1,jj))*deltay+&
         &AC_o(ii,jj)*(1._DBL-rel_vo)*u_o(ii,jj)

    au_o(ii,jj) = AC_o(ii,jj) * rel_vo

    BC_o(jj,ii) = AC_o(ii,jj)

    Ry_o(jj,ii) =-AI_o(ii,jj)*u_o(ii-1,jj)-AD_o(ii,jj)*u_o(ii+1,jj)-&
         &deltax*deltay*Ri_o(ii,jj)*temp_int+&
         &deltax*deltay*fuente_con_uo(ii,jj)+&
         &deltax*deltay*u_anto(ii,jj)/dt_o+&
         &(pres_o(ii,jj)-pres_o(ii+1,jj))*deltay+&
         &BC_o(jj,ii)*(1._DBL-rel_vo)*u_o(ii,jj)
    !
  end subroutine ensambla_velu


  !*******************************************************************
  !*******************************************************************
  !
  ! ensambla_velv
  !
  ! Subrutina que calcula los coeficientes de la matriz tridiagonal
  ! para la velocidad v 
  !
  !*******************************************************************
  !*******************************************************************
  subroutine ensambla_velv(deltaxvo,deltayvo,deltaxuo,&
       &deltaypo,fexpo,feypo,feyvo,gamma_momento,&
       &fuente_con_vo,fuente_lin_vo,&
       &v_o,v_anto,u_o,&
       &temp_o,pres_o,Ri_o,dt_o,rel_vo,&
       &AI_o,AC_o,AD_o,Rx_o,BS_o,BC_o,BN_o,Ry_o,av_o,&
       &ii,jj&
       &)
    !$acc routine
    !
    implicit none
    !
    ! Tama\~no del volumen de control
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxvo
    real(kind=DBL), dimension(nj), intent(in) :: deltayvo
    !
    ! Distancia entre nodos contiguos de la malla de v en direcci\'on horizontal
    !
    real(kind=DBL), dimension(mi), intent(in) :: deltaxuo
    !
    ! Distancia entre nodos contiguos de la malla de v en direcci\'on vertical
    !
    real(kind=DBL), dimension(nj), intent(in) :: deltaypo
    !
    ! Coeficientes para interpolaci\'on
    !
    real(kind=DBL), DIMENSION(mi),   intent(in) :: fexpo
    real(kind=DBL), DIMENSION(nj),   intent(in) :: feypo
    real(kind=DBL), DIMENSION(nj-1), intent(in) :: feyvo
    !
    ! Coeficiente de difusi\'on
    !
    real(kind=DBL), dimension(mi+1,nj+1), intent(in) :: gamma_momento
    !
    ! Velocidad, presi\'on y temperatura
    !
    real(kind=DBL), dimension(mi,nj+1),   intent(in) :: u_o
    real(kind=DBL), dimension(mi+1,nj),   intent(in) :: v_o,    v_anto
    real(kind=DBL), dimension(mi+1,nj+1), intent(in) :: temp_o, pres_o
    !
    ! T\'erminos fuente
    !
    real(kind=DBL), dimension(mi+1,nj),   intent(in) :: fuente_con_vo
    real(kind=DBL), dimension(mi+1,nj),   intent(in) :: fuente_lin_vo 
    real(kind=DBL), dimension(mi+1,nj+1), intent(in) :: Ri_o
    !
    ! Incremento de tiempo y coeficiente de relajaci\'on
    !
    real(kind=DBL), intent(in) :: dt_o, rel_vo
    !
    ! Coeficientes de las matrices
    !
    ! ** Estos coeficientes est\'an sobredimensionados para reducir el uso de memoria
    ! en la gpu, los arreglos que se reciben en esta subrutina se usan para las ecs.
    ! de momento, energ\'ia y la correcci\'on de la presi\'on **
    !
    real(kind=DBL), dimension(mi+1,nj+1), intent(out) :: AI_o, AC_o, AD_o, Rx_o
    real(kind=DBL), dimension(nj+1,mi+1), intent(out) :: BS_o, BC_o, BN_o, Ry_o
    real(kind=DBL), dimension(mi+1,nj),   intent(out) :: av_o
    ! !
    ! ! \'Indice para recorrer la direcci\'on y
    ! !
    integer, intent(in) :: ii, jj
    !
    ! Variables auxiliares
    !
    integer :: kk, info
    !
    ! Auxiliares de interpolaci\'on
    !
    real(kind=DBL) :: ui, ud, vs, vn
    real(kind=DBL) :: di, dd, ds, dn
    real(kind=DBL) :: gammai, gammad
    real(kind=DBL) :: gammas, gamman
    real(kind=DBL) :: temp_int
    !
    !
    ! Interpolaciones necesarias
    !
    ! u
    !
    ud = feypo(jj)*u_o(ii,jj+1)  +(1.0_DBL-feypo(jj))*u_o(ii,jj)
    ui = feypo(jj)*u_o(ii-1,jj+1)+(1.0_DBL-feypo(jj))*u_o(ii-1,jj)
    ! ui = (u_o(ii-1,jj)+u_o(ii,jj))/2._DBL
    ! ud = (u_o(ii,jj)+u_o(ii+1,jj))/2._DBL
    !
    ! v
    !
    vn = feyvo(jj)  *v_o(ii,jj+1)+(1.0_DBL-feyvo(jj))  *v_o(ii,jj)
    vs = feyvo(jj-1)*v_o(ii,jj)  +(1.0_DBL-feyvo(jj-1))*v_o(ii,jj-1)
    !
    ! gamma_d
    !
    ! ** se utilizan las constantes gamman y gammas como auxiliares para
    ! calcular gammai, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammas = ( gamma_momento(ii+1,jj) * gamma_momento(ii,jj) ) / &
         &( gamma_momento(ii+1,jj) * (1._DBL-fexpo(ii))+gamma_momento(ii,jj)*fexpo(ii) )
    gamman = ( gamma_momento(ii+1,jj+1) * gamma_momento(ii,jj+1) ) / &
         &( gamma_momento(ii+1,jj+1) * (1._DBL-fexpo(ii))+gamma_momento(ii,jj+1)*fexpo(ii) )
    gammad = ( gammas * gamman ) / &
         &( gammas*(1._DBL-feypo(jj)) + gamman*feypo(jj) )
    !
    ! gamma_i
    !
    ! ** se utilizan las constantes gamman y gammas como auxiliares para
    ! calcular gammai, despu\'es se utilizan para el coeficiente gamma que
    ! corresponde **
    gammas = ( gamma_momento(ii,jj) * gamma_momento(ii-1,jj) ) / &
         &( gamma_momento(ii,jj) * (1._DBL-fexpo(ii-1))+gamma_momento(ii-1,jj)*fexpo(ii-1) )
    gamman = ( gamma_momento(ii,jj+1) * gamma_momento(ii-1,jj+1) ) / &
         &(gamma_momento(ii,jj+1)*(1._DBL-fexpo(ii-1))+gamma_momento(ii-1,jj+1)*fexpo(ii-1) )
    gammai = ( gammas * gamman ) / &
         &( gammas*(1._DBL-feypo(jj)) + gamman*feypo(jj) )
    !
    ! gamma_n 
    !
    gamman = gamma_momento(ii,jj+1)
    !
    ! gamma_s 
    !
    gammas = gamma_momento(ii,jj)       
    !
    ! distancias entre nodos contiguos
    !
    di = deltaxuo(ii-1)
    dd = deltaxuo(ii)
    ds = deltaypo(jj)
    dn = deltaypo(jj+1)
    !
    ! Tama\~no de los vol\'umenes de control para la velocidad v
    !
    ! delta_x = deltaxvo(ii)
    ! delta_y = deltayvo(jj)
    !
    ! Interpolaci\'on para la temperatura
    !
    temp_int = feypo(jj)*temp_o(ii,jj+1) + (1.0_DBL-feypo(jj))*temp_o(ii,jj)
    !
    ! *************************
    !
    ! Coeficientes de la matriz
    !
    AI_o(ii,jj) =-(gammai*deltayvo(jj)/di*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ui*di/gammai))**5)+&
         &DMAX1(0.0_DBL,ui*deltayvo(jj)))

    AD_o(ii,jj) =-(gammad*deltayvo(jj)/dd*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ud*dd/gammad))**5)+&
         &DMAX1(0.0_DBL,-ud*deltayvo(jj)))

    BS_o(jj,ii) =-(gammas*deltaxvo(ii)/ds*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vs*ds/gammas))**5)+&
         &DMAX1(0.0_DBL, vs*deltaxvo(ii)))

    BN_o(jj,ii) =-(gamman*deltaxvo(ii)/dn*&
         &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vn*dn/gamman))**5)+&
         &DMAX1(0.0_DBL,-vn*deltaxvo(ii)))

    AC_o(ii,jj) = ( -AI_o(ii,jj) - AD_o(ii,jj) - BS_o(jj,ii) - BN_o(jj,ii)-&
         &deltaxvo(ii)*deltayvo(jj)*fuente_lin_vo(ii,jj)+&
         &deltaxvo(ii)*deltayvo(jj)/dt_o) / rel_vo

    Rx_o(ii,jj) =-BS_o(jj,ii)*v_o(ii,jj-1) - BN_o(jj,ii)*v_o(ii,jj+1)-&
         &deltaxvo(ii)*deltayvo(jj)*Ri_o(ii,jj)*temp_int+&
         &deltaxvo(ii)*deltayvo(jj)*fuente_con_vo(ii,jj)+&
         &deltaxvo(ii)*deltayvo(jj)*v_anto(ii,jj)/dt_o+&
         &(pres_o(ii,jj)-pres_o(ii,jj+1))*deltaxvo(ii)+&
         &AC_o(ii,jj)*(1._DBL-rel_vo)*v_o(ii,jj)

    av_o(ii,jj) = AC_o(ii,jj) * rel_vo

    BC_o(jj,ii) = AC_o(ii,jj)

    Ry_o(jj,ii) =-AI_o(ii,jj)*v_o(ii-1,jj)-AD_o(ii,jj)*v_o(ii+1,jj)-&
         &deltaxvo(ii)*deltayvo(jj)*Ri_o(ii,jj)*temp_int+&
         &deltaxvo(ii)*deltayvo(jj)*fuente_con_vo(ii,jj)+&
         &deltaxvo(ii)*deltayvo(jj)*v_anto(ii,jj)/dt_o+&
         &(pres_o(ii,jj)-pres_o(ii,jj+1))*deltaxvo(ii)+&
         &BC_o(jj,ii)*(1._DBL-rel_vo)*v_o(ii,jj)
    ! end do bucle_direccion_x
    !
    ! end do bucle_direccion_y
    !
  end subroutine ensambla_velv
  !
  !*******************************************************************
  !
  ! residuo_u
  !
  ! Subrutina que calcula el residuo de la ecuaci\'on de momento en u
  !
  !*******************************************************************
  subroutine residuo_u(deltaxuo,deltayuo,deltaxpo,&
       &deltayvo,fexpo,feypo,fexuo,gamma_momento,&
       &fuente_con_uo,fuente_lin_uo,&
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
    real(kind=DBL), dimension(mi,nj+1),   intent(in) :: fuente_con_uo
    real(kind=DBL), dimension(mi,nj+1),   intent(in) :: fuente_lin_uo
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
          AC_o = AI_o + AD_o + BS_o + BN_o-&
               &deltaxuo(ii)*deltayuo(jj)*fuente_lin_uo(ii,jj)+&
               &deltaxuo(ii)*deltayuo(jj)/dt_o
          Rx_o(ii,jj) = AI_o*u_o(ii-1,jj)-AC_o*u_o(ii,jj)+AD_o*u_o(ii+1,jj)+&
               &BS_o*u_o(ii,jj-1)+BN_o*u_o(ii,jj+1)-&
               &deltaxuo(ii)*deltayuo(jj)*Ri_o(ii,jj)*temp_int+&
               &deltaxuo(ii)*deltayuo(jj)*fuente_con_uo(ii,jj)+&
               &deltaxuo(ii)*deltayuo(jj)*u_anto(ii,jj)/dt_o+&
               &(pres_o(ii,jj)-pres_o(ii+1,jj))*deltayuo(jj)
       end do bucle_direccion_x
       !
    end do bucle_direccion_y
    !
   end subroutine residuo_u
   !
end module ec_momento
