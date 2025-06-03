!
!
! M\'odulo para la ecuaci\'on de continuidad
!
! Este m\'odulo contiene:
!
! - Las variables de la ecuaci\'on de continuidad
! - Las subrutinas que calculan los coeficientes de las ecuaciones discretizadas
! - Las subrutinas que sirven para leer e imponer condiciones de frontera
!
module ec_continuidad
  !
  use malla,         only : mi, nj, lk, DBL
  use malla,         only : xu, yv, zw, xp, yp, zp
  !
  use cond_frontera, only : tipo_cond_front
  !
  use cond_frontera, only : inicializa_cond_front
  use cond_frontera, only : lectura_cond_frontera
  !
  implicit none
  !
  ! !
  ! ! Componentes de velocidad, presi\'on
  ! ! residuos de las ecuaciones de momento y correcci\'on de la presi\'on
  ! ! coeficientes de difusi\'on y criterios e convergencia
  ! !
  ! real(kind=DBL), dimension(mi,nj+1,lk+1)   :: u, u_ant, du, au, Resu, fu
  ! real(kind=DBL), dimension(mi+1,nj,lk+1)   :: v, v_ant, dv, av, Resv, fv
  ! real(kind=DBL), dimension(mi+1,nj+1,lk)   :: w, w_ant, dw, aw, Resw, fw
  ! real(kind=DBL), dimension(mi+1,nj+1,lk+1) :: gamma_momen, Ri
  ! !
  ! ! Variables para los t\'erminos fuente de la ec de momento
  ! ! El t\'ermino lineal fuente_lin debe ser negativo para favorecer
  ! ! convergencia y estabilidad
  ! !
  ! real(kind=DBL), dimension(mi,nj+1,lk+1)   :: fuente_con_u, fuente_lin_u
  ! real(kind=DBL), dimension(mi+1,nj,lk+1)   :: fuente_con_v, fuente_lin_v
  ! real(kind=DBL), dimension(mi+1,nj+1,lk)   :: fuente_con_w, fuente_lin_w
  !
  ! Variables para la presi\'on
  !
  real(kind=DBL), dimension(mi+1,nj+1,lk+1) :: pres, corr_pres
  real(kind=DBL), dimension(mi+1,nj+1,lk+1) :: dcorr_pres, fcorr_pres
  real(kind=DBL), dimension(mi+1,nj+1,lk+1) :: b_o
  !
  ! Variables de convergencia y relajaci\'on
  !
  real(kind=DBL) :: maxbo, conv_p, rel_pres
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
  ! ensambla_corr_pres
  !
  ! Subrutina que calcula los coeficientes de la matriz tridiagonal
  ! para la correcci\'on de la presi\'on 
  !
  !*******************************************************************
  subroutine ensambla_corr_pres_x(&
       &deltaxpo,&
       &deltaypo,&
       &deltazpo,&
       &deltaxuo,&
       &deltayvo,&
       &deltazwo,&
       &u_o,v_o,w_o,&
       &corr_pres_o,&
       &rel_vo,&
       &AI_o,AC_o,AD_o,Rx_o,&
       &au_o,av_o,aw_o,&
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
    ! Velocidad, presi\'on, t\'ermino fuente b
    !
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in)  :: u_o
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in)  :: v_o
    real(kind=DBL), dimension(mi+1,nj+1,lk),   intent(in)  :: w_o
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in)  :: corr_pres_o
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
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(out) :: AI_o, AC_o, AD_o, Rx_o
    ! real(kind=DBL), dimension(nj+1,mi+1), intent(out) :: BS_o, BC_o, BN_o, Ry_o
    real(kind=DBL), dimension(mi,nj+1,lk+1), intent(in) :: au_o
    real(kind=DBL), dimension(mi+1,nj,lk+1), intent(in) :: av_o
    real(kind=DBL), dimension(mi+1,nj+1,lk), intent(in) :: aw_o
    !
    ! \'Indice para recorrer las direcciones x, y. z
    !
    integer,   intent(in) :: ii, jj, kk
    !
    ! Variables auxiliares
    !
    ! integer :: ii, jj
    !
    ! Auxiliares de interpolaci\'on
    !
    real(kind=DBL) :: ui, ud, vs, vn, wt, wb
    real(kind=DBL) :: ai, ad, as, an, at, ab
    real(kind=DBL) :: alpha,beta,gamma,delta
    real(kind=DBL) :: b_o
    ! real(kind=DBL) :: di, dd, ds, dn
    !
    ! C\'alculo de los coeficientes
    !
    ! $acc loop gang
    ! bucle_direccion_y: do jj = 2, nj
    !------------------------
    ! Condiciones de frontera
    AC_o(1,jj,kk) = 1._DBL
    AD_o(1,jj,kk) = 0._DBL
    Rx_o(1,jj,kk) = 0._DBL
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
    ! coeficientes de la ecuaci\'on de momento
    !
    ai = au_o(ii-1,jj,kk)
    ad = au_o(ii,jj,kk)
    as = av_o(ii,jj-1,kk)
    an = av_o(ii,jj,kk)
    ab = aw_o(ii,jj,kk-1)
    at = aw_o(ii,jj,kk)
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
    AI_o(ii,jj,kk) =-deltaypo(jj)*deltaypo(jj)*deltazpo(kk)*deltazpo(kk)/ai
    AD_o(ii,jj,kk) =-deltaypo(jj)*deltaypo(jj)*deltazpo(kk)*deltazpo(kk)/ad
    alpha          =-deltaxpo(ii)*deltaxpo(ii)*deltazpo(kk)*deltazpo(kk)/as
    beta           =-deltaxpo(ii)*deltaxpo(ii)*deltazpo(kk)*deltazpo(kk)/an
    gamma          =-deltaxpo(ii)*deltaxpo(ii)*deltaypo(jj)*deltaypo(jj)/ab
    delta          =-deltaxpo(ii)*deltaxpo(ii)*deltaypo(jj)*deltaypo(jj)/at
    AC_o(ii,jj,kk) = ( -AI_o(ii,jj,kk) - AD_o(ii,jj,kk)-&
         &alpha - beta - gamma - delta ) / rel_vo
    !
    b_o = (ui-ud)*deltaypo(jj)*deltazpo(kk)+&
         &(vs-vn)*deltaxpo(ii)*deltazpo(kk)+&
         &(wb-wt)*deltaxpo(ii)*deltaypo(jj)
    !
    Rx_o(ii,jj,kk) =-alpha*corr_pres_o(ii,jj-1,kk)-&
         &beta *corr_pres_o(ii,jj+1,kk)+&
         &gamma*corr_pres_o(ii,jj,kk-1)+&
         &delta*corr_pres_o(ii,jj,kk+1)+&
         &b_o +&
         &(1._DBL-rel_vo)*AC_o(ii,jj,kk)*corr_pres_o(ii,jj,kk)
    ! end do bucle_direccion_x
    !
    ! Condicion frontera
    !
    AI_o(mi+1,jj,kk) = 0.0_DBL
    AC_o(mi+1,jj,kk) = 1.0_DBL
    Rx_o(mi+1,jj,kk) = 0.0_DBL
    !
  end subroutine ensambla_corr_pres_x
  !
  !*******************************************************************
  !
  ! ensambla_corr_pres
  !
  ! Subrutina que calcula los coeficientes de la matriz tridiagonal
  ! para la correcci\'on de la presi\'on en una l\'inea en direcci\'on
  ! y
  !
  !*******************************************************************
  subroutine ensambla_corr_pres_y(&
       &deltaxpo,&
       &deltaypo,&
       &deltazpo,&
       &deltaxuo,&
       &deltayvo,&
       &deltazwo,&
       &u_o,v_o,w_o,&
       &corr_pres_o,&
       &rel_vo,&
       &AI_o,AC_o,AD_o,Rx_o,&
       &au_o,av_o,aw_o,&
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
    ! Velocidad, presi\'on, t\'ermino fuente b
    !
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in)  :: u_o
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in)  :: v_o
    real(kind=DBL), dimension(mi+1,nj+1,lk),   intent(in)  :: w_o
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in)  :: corr_pres_o
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
    real(kind=DBL), dimension(nj+1,mi+1,lk+1), intent(out) :: AI_o, AC_o, AD_o, Rx_o
    ! real(kind=DBL), dimension(nj+1,mi+1), intent(out) :: BS_o, BC_o, BN_o, Ry_o
    real(kind=DBL), dimension(mi,nj+1,lk+1), intent(in) :: au_o
    real(kind=DBL), dimension(mi+1,nj,lk+1), intent(in) :: av_o
    real(kind=DBL), dimension(mi+1,nj+1,lk), intent(in) :: aw_o
    !
    ! \'Indice para recorrer las direcciones x, y. z
    !
    integer,   intent(in) :: ii, jj, kk
    !
    ! Variables auxiliares
    !
    ! integer :: ii, jj
    !
    ! Auxiliares de interpolaci\'on
    !
    real(kind=DBL) :: ui, ud, vs, vn, wt, wb
    real(kind=DBL) :: ai, ad, as, an, at, ab
    real(kind=DBL) :: alpha,beta,gamma,delta
    real(kind=DBL) :: b_o
    ! real(kind=DBL) :: di, dd, ds, dn
    !
    ! C\'alculo de los coeficientes
    !
    ! $acc loop gang
    ! bucle_direccion_y: do jj = 2, nj
    !------------------------
    ! Condiciones de frontera
    AC_o(1,ii,kk) = 1._DBL
    AD_o(1,ii,kk) = 0._DBL
    Rx_o(1,ii,kk) = 0._DBL
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
    ! coeficientes de la ecuaci\'on de momento
    !
    ai = au_o(ii-1,jj,kk)
    ad = au_o(ii,jj,kk)
    as = av_o(ii,jj-1,kk)
    an = av_o(ii,jj,kk)
    ab = aw_o(ii,jj,kk-1)
    at = aw_o(ii,jj,kk)
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
    alpha          =-deltaypo(jj)*deltaypo(jj)*deltazpo(kk)*deltazpo(kk)/ai
    beta           =-deltaypo(jj)*deltaypo(jj)*deltazpo(kk)*deltazpo(kk)/ad
    AI_o(jj,ii,kk) =-deltaxpo(ii)*deltaxpo(ii)*deltazpo(kk)*deltazpo(kk)/as
    AD_o(jj,ii,kk) =-deltaxpo(ii)*deltaxpo(ii)*deltazpo(kk)*deltazpo(kk)/an
    gamma          =-deltaxpo(ii)*deltaxpo(ii)*deltaypo(jj)*deltaypo(jj)/ab
    delta          =-deltaxpo(ii)*deltaxpo(ii)*deltaypo(jj)*deltaypo(jj)/at
    AC_o(jj,ii,kk) = ( -AI_o(jj,ii,kk) - AD_o(jj,ii,kk)-&
         &alpha - beta - gamma - delta ) / rel_vo
    !
    b_o = (ui-ud)*deltaypo(jj)*deltazpo(kk)+&
         &(vs-vn)*deltaxpo(ii)*deltazpo(kk)+&
         &(wb-wt)*deltaxpo(ii)*deltaypo(jj)
    !
    Rx_o(jj,ii,kk) =-alpha*corr_pres_o(ii,jj-1,kk)-&
         &beta *corr_pres_o(ii,jj+1,kk)+&
         &gamma*corr_pres_o(ii,jj,kk-1)+&
         &delta*corr_pres_o(ii,jj,kk+1)+&
         &b_o+&
         &(1._DBL-rel_vo)*AC_o(jj,ii,kk)*corr_pres_o(ii,jj,kk)
    ! end do bucle_direccion_x
    !
    ! Condicion frontera
    !
    AI_o(nj+1,ii,kk) = 0.0_DBL
    AC_o(nj+1,ii,kk) = 1.0_DBL
    Rx_o(nj+1,ii,kk) = 0.0_DBL
    !
  end subroutine ensambla_corr_pres_y
  !
  !*******************************************************************
  !
  ! ensambla_corr_pres
  !
  ! Subrutina que calcula los coeficientes de la matriz tridiagonal
  ! para la correcci\'on de la presi\'on en una l\'inea en direcci\'on
  ! y
  !
  !*******************************************************************
  subroutine ensambla_corr_pres_z(&
       &deltaxpo,&
       &deltaypo,&
       &deltazpo,&
       &deltaxuo,&
       &deltayvo,&
       &deltazwo,&
       &u_o,v_o,w_o,&
       &corr_pres_o,&
       &b_o,&
       &rel_vo,&
       &AI_o,AC_o,AD_o,Rx_o,&
       &au_o,av_o,aw_o,&
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
    ! Velocidad, presi\'on, t\'ermino fuente b
    !
    real(kind=DBL), dimension(mi,nj+1,lk+1),   intent(in)  :: u_o
    real(kind=DBL), dimension(mi+1,nj,lk+1),   intent(in)  :: v_o
    real(kind=DBL), dimension(mi+1,nj+1,lk),   intent(in)  :: w_o
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(in)  :: corr_pres_o
    real(kind=DBL), dimension(mi+1,nj+1,lk+1), intent(out) :: b_o
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
    real(kind=DBL), dimension(lk+1,mi+1,nj+1), intent(out) :: AI_o, AC_o, AD_o, Rx_o
    ! real(kind=DBL), dimension(nj+1,mi+1), intent(out) :: BS_o, BC_o, BN_o, Ry_o
    real(kind=DBL), dimension(mi,nj+1,lk+1), intent(in) :: au_o
    real(kind=DBL), dimension(mi+1,nj,lk+1), intent(in) :: av_o
    real(kind=DBL), dimension(mi+1,nj+1,lk), intent(in) :: aw_o
    !
    ! \'Indice para recorrer las direcciones x, y. z
    !
    integer,   intent(in) :: ii, jj, kk
    !
    ! Variables auxiliares
    !
    ! integer :: ii, jj
    !
    ! Auxiliares de interpolaci\'on
    !
    real(kind=DBL) :: ui, ud, vs, vn, wt, wb
    real(kind=DBL) :: ai, ad, as, an, at, ab
    real(kind=DBL) :: alpha,beta,gamma,delta
    ! real(kind=DBL) :: di, dd, ds, dn
    !
    ! C\'alculo de los coeficientes
    !
    ! $acc loop gang
    ! bucle_direccion_y: do jj = 2, nj
    !------------------------
    ! Condiciones de frontera
    AC_o(1,ii,jj) = 1._DBL
    AD_o(1,ii,jj) = 0._DBL
    Rx_o(1,ii,jj) = 0._DBL
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
    ! coeficientes de la ecuaci\'on de momento
    !
    ai = au_o(ii-1,jj,kk)
    ad = au_o(ii,jj,kk)
    as = av_o(ii,jj-1,kk)
    an = av_o(ii,jj,kk)
    ab = aw_o(ii,jj,kk-1)
    at = aw_o(ii,jj,kk)
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
    alpha          =-deltaypo(jj)*deltaypo(jj)*deltazpo(kk)*deltazpo(kk)/ai
    beta           =-deltaypo(jj)*deltaypo(jj)*deltazpo(kk)*deltazpo(kk)/ad
    gamma          =-deltaxpo(ii)*deltaxpo(ii)*deltazpo(kk)*deltazpo(kk)/as
    delta          =-deltaxpo(ii)*deltaxpo(ii)*deltazpo(kk)*deltazpo(kk)/an
    AI_o(kk,ii,jj) =-deltaxpo(ii)*deltaxpo(ii)*deltaypo(jj)*deltaypo(jj)/ab
    AD_o(kk,ii,jj) =-deltaxpo(ii)*deltaxpo(ii)*deltaypo(jj)*deltaypo(jj)/at
    AC_o(kk,ii,jj) = ( -AI_o(kk,ii,jj) - AD_o(kk,ii,jj)-&
         &alpha - beta - gamma - delta ) / rel_vo
    !
    b_o(kk,ii,jj) = (ui-ud)*deltaypo(jj)*deltazpo(kk)+&
         &(vs-vn)*deltaxpo(ii)*deltazpo(kk)+&
         &(wb-wt)*deltaxpo(ii)*deltaypo(jj)
    !
    Rx_o(kk,ii,jj) =-alpha*corr_pres_o(ii,jj-1,kk)-&
         &beta *corr_pres_o(ii,jj+1,kk)+&
         &gamma*corr_pres_o(ii,jj,kk-1)+&
         &delta*corr_pres_o(ii,jj,kk+1)+&
         &b_o+&
         &(1._DBL-rel_vo)*AC_o(kk,ii,jj)*corr_pres_o(ii,jj,kk)
    ! end do bucle_direccion_x
    !
    ! Condicion frontera
    !
    AI_o(lk+1,ii,jj) = 0.0_DBL
    AC_o(lk+1,ii,jj) = 1.0_DBL
    Rx_o(lk+1,ii,jj) = 0.0_DBL
    !
  end subroutine ensambla_corr_pres_z
  !
end module ec_continuidad
