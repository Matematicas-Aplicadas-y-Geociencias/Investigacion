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
  ! Variables para la presi\'on
  !
  real(kind=DBL), dimension(mi+1,nj+1,lk+1) :: pres, corr_pres
  real(kind=DBL), dimension(mi+1,nj+1,lk+1) :: dcorr_pres, fcorr_pres
  real(kind=DBL), dimension(lk+1,nj+1,mi+1) :: b_o
  !
  ! Variables de convergencia y relajaci\'on
  !
  real(kind=DBL) :: maxbo, conv_p, rel_pres
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
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AI_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AC_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AD_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: RX_o
    !
    real(kind=DBL), dimension(mi,nj+1,lk+1), intent(in)          :: au_o
    real(kind=DBL), dimension(mi+1,nj,lk+1), intent(in)          :: av_o
    real(kind=DBL), dimension(mi+1,nj+1,lk), intent(in)          :: aw_o
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
    real(kind=DBL) :: ai, ad, as, an, at, ab
    real(kind=DBL) :: alpha,beta,gamma,delta
    real(kind=DBL) :: bo
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
    AI_o(indexp(ii,jj,kk)) =-deltaypo(jj)*deltaypo(jj)*deltazpo(kk)*deltazpo(kk)/ai
    AD_o(indexp(ii,jj,kk)) =-deltaypo(jj)*deltaypo(jj)*deltazpo(kk)*deltazpo(kk)/ad
    alpha                  =-deltaxpo(ii)*deltaxpo(ii)*deltazpo(kk)*deltazpo(kk)/as
    beta                   =-deltaxpo(ii)*deltaxpo(ii)*deltazpo(kk)*deltazpo(kk)/an
    gamma                  =-deltaxpo(ii)*deltaxpo(ii)*deltaypo(jj)*deltaypo(jj)/ab
    delta                  =-deltaxpo(ii)*deltaxpo(ii)*deltaypo(jj)*deltaypo(jj)/at
    !
    AC_o(indexp(ii,jj,kk)) = ( -AI_o(indexp(ii,jj,kk)) - AD_o(indexp(ii,jj,kk))-&
         &alpha - beta - gamma - delta ) / rel_vo
    !
    bo = (ui-ud)*deltaypo(jj)*deltazpo(kk)+&
         &(vs-vn)*deltaxpo(ii)*deltazpo(kk)+&
         &(wb-wt)*deltaxpo(ii)*deltaypo(jj)
    !
    Rx_o(indexp(ii,jj,kk)) =-alpha*corr_pres_o(ii,jj-1,kk)-&
         &beta *corr_pres_o(ii,jj+1,kk)+&
         &gamma*corr_pres_o(ii,jj,kk-1)+&
         &delta*corr_pres_o(ii,jj,kk+1)+&
         &bo +&
         &(1._DBL-rel_vo)*AC_o(indexp(ii,jj,kk))*corr_pres_o(ii,jj,kk)
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
       &jj,ii,kk)
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
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AI_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AC_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AD_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: RX_o
    !
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
    real(kind=DBL) :: bo
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
    alpha                  =-deltaypo(jj)*deltaypo(jj)*deltazpo(kk)*deltazpo(kk)/ai
    beta                   =-deltaypo(jj)*deltaypo(jj)*deltazpo(kk)*deltazpo(kk)/ad
    AI_o(indeyp(jj,ii,kk)) =-deltaxpo(ii)*deltaxpo(ii)*deltazpo(kk)*deltazpo(kk)/as
    AD_o(indeyp(jj,ii,kk)) =-deltaxpo(ii)*deltaxpo(ii)*deltazpo(kk)*deltazpo(kk)/an
    gamma                  =-deltaxpo(ii)*deltaxpo(ii)*deltaypo(jj)*deltaypo(jj)/ab
    delta                  =-deltaxpo(ii)*deltaxpo(ii)*deltaypo(jj)*deltaypo(jj)/at
    AC_o(indeyp(jj,ii,kk)) = ( -AI_o(indeyp(jj,ii,kk)) - AD_o(indeyp(jj,ii,kk))-&
         &alpha - beta - gamma - delta ) / rel_vo
    !
    bo = (ui-ud)*deltaypo(jj)*deltazpo(kk)+&
         &(vs-vn)*deltaxpo(ii)*deltazpo(kk)+&
         &(wb-wt)*deltaxpo(ii)*deltaypo(jj)
    !
    Rx_o(indeyp(jj,ii,kk)) =-alpha*corr_pres_o(ii-1,jj,kk)-&
         &beta *corr_pres_o(ii+1,jj,kk)+&
         &gamma*corr_pres_o(ii,jj,kk-1)+&
         &delta*corr_pres_o(ii,jj,kk+1)+&
         &bo+&
         &(1._DBL-rel_vo)*AC_o(indeyp(jj,ii,kk))*corr_pres_o(ii,jj,kk)
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
       &b_oo,&
       &rel_vo,&
       &AI_o,AC_o,AD_o,Rx_o,&
       &au_o,av_o,aw_o,&
       &kk,jj,ii)
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
    real(kind=DBL), dimension(lk+1,mi+1,nj+1), intent(out) :: b_oo
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
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AI_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AC_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: AD_o
    real(kind=DBL), dimension((mi+1)*(nj+1)*(lk+1)), intent(out) :: RX_o
    !    
    real(kind=DBL), dimension(mi,nj+1,lk+1), intent(in) :: au_o
    real(kind=DBL), dimension(mi+1,nj,lk+1), intent(in) :: av_o
    real(kind=DBL), dimension(mi+1,nj+1,lk), intent(in) :: aw_o
    !
    ! \'Indice para recorrer las direcciones x, y. z
    !
    integer,   intent(in) :: ii, jj, kk
    !
    ! Auxiliares de interpolaci\'on
    !
    real(kind=DBL) :: ui, ud, vs, vn, wt, wb
    real(kind=DBL) :: ai, ad, as, an, at, ab
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
    AI_o(indezp(kk,jj,ii)) =-deltaxpo(ii)*deltaxpo(ii)*deltaypo(jj)*deltaypo(jj)/ab
    AD_o(indezp(kk,jj,ii)) =-deltaxpo(ii)*deltaxpo(ii)*deltaypo(jj)*deltaypo(jj)/at
    AC_o(indezp(kk,jj,ii)) = ( -AI_o(indezp(kk,jj,ii)) - AD_o(indezp(kk,jj,ii))-&
         &alpha - beta - gamma - delta ) / rel_vo
    !
    b_oo(kk,jj,ii) = (ui-ud)*deltaypo(jj)*deltazpo(kk)+&
         &(vs-vn)*deltaxpo(ii)*deltazpo(kk)+&
         &(wb-wt)*deltaxpo(ii)*deltaypo(jj)
    !
    Rx_o(indezp(kk,jj,ii)) =-alpha*corr_pres_o(ii-1,jj,kk)-&
         &beta *corr_pres_o(ii+1,jj,kk)-&
         &gamma*corr_pres_o(ii,jj-1,kk)-&
         &delta*corr_pres_o(ii,jj+1,kk)+&
         &b_oo(kk,jj,ii)+&
         &(1._DBL-rel_vo)*AC_o(indezp(kk,jj,ii))*corr_pres_o(ii,jj,kk)
    !
  end subroutine ensambla_corr_pres_z
  !
end module ec_continuidad
