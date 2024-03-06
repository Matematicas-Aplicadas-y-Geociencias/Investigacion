!
!
! M\'odulo de ensamblaje
!
! Este m\'odulo contiene las subrutinas que calculan los coeficientes de las matrices
! tridiagonales a invertir en el barrido l\'inea por l\'inea
!
!
MODULE ensamblaje
  use malla
  implicit none
  
contains
  !*******************************************************************
  !
  ! ensambla_velu
  !
  ! Subrutina que calcula los coeficientes de la matriz tridiagonal
  ! para la velocidad u barriendo la direcci\'on x
  !
  !*******************************************************************
  subroutine ensambla_velu(deltaxuo,deltayuo,deltaxpo,&
       &deltayvo,fexpo,feypo,fexuo,gamma_momento,&
       &u_o,u_anto,v_o,&
       &temp_o,pres_o,Ri_o,dt_o,rel_vo,&
       &AI_o,AC_o,AD_o,Rx_o,BS_o,BC_o,BN_o,Ry_o,au_o)
    implicit none
    !
    ! Tama\~no del volumen de control
    !
    real(kind=DBL), dimension(mi),   intent(in) :: deltaxuo
    real(kind=DBL), dimension(nj),   intent(in) :: deltayuo
    !
    ! Distancia entre nodos contiguos de la malla de u en direcci\'on horizontal
    !
    real(kind=DBL), dimension(mi),   intent(in) :: deltaxpo
    !
    ! Distancia entre nodos contiguos de la malla de u en direcci\'on vertical
    !
    real(kind=DBL), dimension(nj),   intent(in) :: deltayvo
    !
    ! Coeficientes para interpolaci\'on
    !
    real(kind=DBL), DIMENSION(mi), intent(in)   ::  fexpo
    real(kind=DBL), DIMENSION(nj), intent(in)   ::  feypo
    real(kind=DBL), DIMENSION(mi-1), intent(in) ::  fexuo
    !
    ! Coeficiente de difusi\'on
    !
    real(kind=DBL), dimension(mi+1,nj+1), intent(in) ::  gamma_momento
    !
    ! Velocidad, presi\'on y temperatura
    !
    real(kind=DBL), dimension(mi,nj+1),   intent(in)    :: u_o,    u_anto
    real(kind=DBL), dimension(mi+1,nj),   intent(in)    :: v_o
    real(kind=DBL), dimension(mi+1,nj+1), intent(in)    :: temp_o, pres_o
    !
    ! T\'erminos fuente
    !
    real(kind=DBL), dimension(mi,nj+1),   intent(in)    :: Ri_o
    !
    ! Incremento de tiempo y coeficiente de relajaci\'on
    !
    real(kind=DBL), intent(in) :: dt_o, rel_vo
    !
    ! Coeficientes de las matrices
    !
    ! ** Estos coeficientes est\'an sobredimensionados para reducir el uso de memoria
    ! en la gpu, los arreglos que se reciben en esta subrutina se usan para las ecs.
    ! de momento en y, energ\'ia y la correcci\'on de la presi\'on **
    !
    real(kind=DBL), dimension(mi+1,nj+1), intent(out)   :: AI_o, AC_o, AD_o, Rx_o
    real(kind=DBL), dimension(nj+1,mi+1), intent(out)   :: BS_o, BC_o, BN_o, Ry_o
    real(kind=DBL), dimension(mi,nj+1),   intent(out)   :: au_o
    !
    ! Variables auxiliares
    !
    integer :: ii, jj, kk, info
    !
    ! Auxiliares de interpolaci\'on
    !
    real(kind=DBL) :: ui, ud, vs, vn
    real(kind=DBL) :: di, dd, ds, dn
    real(kind=DBL) :: gammai, gammad
    real(kind=DBL) :: gammas, gamman
    real(kind=DBL) :: temp_int
    !
    ! ********************
    ! TDMA en direccion x
    ! $OMP  PARALLEL DO DEFAULT(NONE)&
    ! $OMP& PRIVATE(AI,ACi,AD,Rxi,ui,ud,vs,vn,di,dd,ds,dn,gamma_i,gamma_d,gamma_s,gamma_n,&
    ! $OMP& delta_x,delta_y,temp_int,alpha,beta,info)&
    ! $OMP& SHARED(u_o,fu,u_anto,pres_o,gamma_uo,au_o,v_o,d_xuo,d2_xuo,yo,fey,d_yvo,temp_o,&
    ! $OMP& rel_vo,Ri_o,dt_o)
    bucle_direccion_y: do jj = 2, nj
       !***********************
       !Condiciones de frontera
       AC_o(1,jj) = 1._DBL
       AD_o(1,jj) = 0._DBL
       Rx_o(1,jj) = 1._DBL !-6.d0*dfloat(j-1)/dfloat(nj)*(dfloat(j-1)/dfloat(nj)-1.d0)! 1._DBL !
       au_o(1,jj) = 1.e40_DBL !ACi(1)
       !
       ! Llenado de la matriz
       !
       bucle_direccion_x: do ii = 2, mi-1
          !
          ! Interpolaciones necesarias
          !
          ! u
          !
          ! ud = fexuo(ii)  *u_o(ii+1,jj)+(1.0_DBL-fexuo(ii))  *u_o(ii,jj)
          ! ui = fexuo(ii-1)*u_o(ii,jj)  +(1.0_DBL-fexuo(ii-1))*u_o(ii-1,jj)
          ui = (u_o(ii-1,jj)+u_o(ii,jj))/2._DBL
          ud = (u_o(ii,jj)+u_o(ii+1,jj))/2._DBL
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
          gammad = ( gamma_momento(ii,jj+1) * gamma_momento(ii,jj) ) / &
               &( gamma_momento(ii,jj+1) * (1._DBL-feypo(jj))+gamma_momento(ii,jj)*feypo(jj) )
          gammai = ( gamma_momento(ii-1,jj+1) * gamma_momento(ii-1,jj) ) / &
               &( gamma_momento(ii-1,jj+1) * (1._DBL-feypo(jj))+gamma_momento(ii-1,jj)*feypo(jj) )
          !
          gamman = gammai*gammad / (gammad * (1._DBL-fexpo(ii)) + gammai * fexpo(ii))
          !
          ! gamma_s 
          !
          ! ** se utilizan las constantes gammai y gammad como auxiliares para
          ! calcular gamman, despu\'es se utilizan para el coeficiente gamma que
          ! corresponde **
          gammad = ( gamma_momento(ii,jj) * gamma_momento(ii,jj-1) ) / &
               &( gamma_momento(ii,jj) * (1._DBL-feypo(jj-1))+gamma_momento(ii,jj-1)*feypo(jj-1) )
          gammai = ( gamma_momento(ii-1,jj+1) * gamma_momento(ii-1,jj) ) / &
               &(gamma_momento(ii-1,jj)*(1._DBL-feypo(jj-1))+gamma_momento(ii-1,jj-1)*feypo(jj-1))
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
          AI_o(ii,jj) =-(gammai*deltayuo(jj)/di*&
               &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ui*di/gammai))**5)+&
               &DMAX1(0.0_DBL,ui*deltayuo(jj)))
          AD_o(ii,jj) =-(gammad*deltayuo(jj)/dd*&
               &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ud*dd/gammad))**5)+&
               &DMAX1(0.0_DBL,-ud*deltayuo(jj)))
          BS_o(jj,ii) =-gammas*deltaxuo(ii)/ds*&
               &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vs*ds/gammas))**5)+&
               &DMAX1(0.0_DBL, vs*deltaxuo(ii))
          BN_o(jj,ii) =-gamman*deltaxuo(ii)/dn*&
               &DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vn*dn/gamman))**5)+&
               &DMAX1(0.0_DBL,-vn*deltaxuo(ii))
          AC_o(ii,jj) = ( -AI_o(ii,jj) - AD_o(ii,jj) - BS_o(jj,ii) - BN_o(jj,ii)+&
               &deltaxuo(ii)*deltayuo(jj)/dt_o) / rel_vo
          Rx_o(ii,jj) =-BS_o(jj,ii)*u_o(ii,jj-1) - BN_o(jj,ii)*u_o(ii,jj+1)-&
               &deltaxuo(ii)*deltayuo(jj)*Ri_o(ii,jj)*temp_int+&
               &deltaxuo(ii)*deltayuo(jj)*u_anto(ii,jj)/dt_o+&
               &(pres_o(ii,jj)-pres_o(ii+1,jj))*deltayuo(jj)+&
               &AC_o(ii,jj)*(1._DBL-rel_vo)*u_o(ii,jj)
          au_o(ii,jj) = AC_o(ii,jj) * rel_vo
          BC_o(jj,ii) = AC_o(ii,jj)
          Ry_o(jj,ii) =-AI_o(ii,jj)*u_o(ii-1,jj)-AD_o(ii,jj)*u_o(ii+1,jj)-&
               &deltaxuo(ii)*deltayuo(jj)*Ri_o(ii,jj)*temp_int+&
               &deltaxuo(ii)*deltayuo(jj)*u_anto(ii,jj)/dt_o+&
               &(pres_o(ii,jj)-pres_o(ii+1,jj))*deltayuo(jj)+&
               &BC_o(jj,ii)*(1._DBL-rel_vo)*u_o(ii,jj)
       end do bucle_direccion_x
       !***********************
       !Condiciones de frontera
       AC_o(mi,jj) = 1._DBL
       AI_o(mi,jj) =-1._DBL !
       Rx_o(mi,jj) = 0.0_DBL
       au_o(mi,jj) = 1.e40_DBL !ACi(mi)
    end do bucle_direccion_y
    !
    ! Condiciones de frontera para la direcci\'on y
    !
    bucle_direccionx: do ii = 2, mi-1
       !***********************
       !Condiciones de frontera
       BC_o(1,ii)     = 1._DBL
       BN_o(1,ii)     = 0.0_DBL
       Ry_o(1,ii)     = 0.0_DBL
       au_o(ii,1)     = 1.e40_DBL !ACj(1)
       BC_o(nj+1,ii)  = 1._DBL
       BS_o(nj+1,ii)  = 0.0_DBL
       Ry_o(nj+1,ii)  = 0.0_DBL
       au_o(ii,nj+1)  = 1.e40_DBL !ACj(nj+1)       
    end do bucle_direccionx
  end subroutine ensambla_velu


  ! !*******************************************************************
  ! !
  ! ! ensambla_veluy
  ! !
  ! ! Subrutina que calcula los coeficientes de la matriz tridiagonal
  ! ! para la velocidad u barriendo la direcci\'on y
  ! !
  ! !*******************************************************************
  ! subroutine ensambla_veluy(yo,fey,d_xuo,d2_xuo,d_yvo,u_o,u_anto,v_o,&
  !      &temp_o,pres_o,gamma_uo,Ri_o,dt_o,du_o,au_o,rel_vo,&
  !      AIuy,ACuy,ADuy,Ruy)
  !   !**********************************
  !   ! Variables del flujo e inc'ognitas
  !   real(kind=DBL), dimension(mi,nj+1),   intent(in)    :: u_o
  !   real(kind=DBL), dimension(mi,nj+1),   intent(out)   :: au_o
  !   real(kind=DBL), dimension(mi,nj+1),   intent(in)    :: u_anto,gamma_uo,Ri_o
  !   real(kind=DBL), dimension(mi+1,nj),   intent(in)    :: v_o
  !   real(kind=DBL), dimension(mi+1,nj+1), intent(in)    :: temp_o,pres_o
  !   real(kind=DBL), dimension(mi,nj+1)                  :: fu
  !   !**********************************
  !   !Variables de la tridiagonal
  !   real(kind=DBL), dimension(nj+1),intent(out)         :: ACj,Ryj
  !   real(kind=DBL), dimension(nj),  intent(out)         :: AS,AN
  !   !*********************
  !   !Variables de interpolaci'on
  !   real(kind=DBL) :: ui,ud,vn,vs,di,dd,ds,dn,gamma_i,gamma_d,gamma_s,gamma_n
  !   real(kind=DBL) :: delta_x,delta_y,temp_int
  !   !*********************
  !   !Variables de la malla,volumen de control,incremento de tiempo y num Richardson
  !   real(kind=DBL), dimension(nj+1), intent(in) :: yo,fey
  !   real(kind=DBL), dimension(mi-1), intent(in) :: d_xuo,d2_xuo
  !   real(kind=DBL), dimension(nj-1), intent(in) :: d_yvo
  !   real(kind=DBL), intent(in) :: dt_o,rel_vo
  !   !********************
  !   !Variables auxiliares
  !   real(kind=DBL) :: alpha,beta

  !   !$OMP  PARALLEL DO DEFAULT(NONE)&
  !   !$OMP& PRIVATE(AS,ACj,AN,Ryj,ui,ud,vs,vn,di,dd,ds,dn,gamma_i,gamma_d,gamma_s,&
  !   !$OMP& gamma_n,delta_x,delta_y,temp_int,alpha,beta,info)&
  !   !$OMP& SHARED(gamma_uo,au_o,u_o,u_anto,pres_o,v_o,d_xuo,d2_xuo,&
  !   !$OMP& yo,fey,d_yvo,temp_o,rel_vo,Ri_o,dt_o)
  !   DO i = 2, mi-1
  !      !***********************
  !      !Condiciones de frontera
  !      ACj(1) = 1._DBL
  !      AN(1)  = 0.0_DBL
  !      Ryj(1) = 0.0_DBL
  !      au_o(i,1) =1.e40_DBL !ACj(1)
  !      !*************************
  !      !Llenado de la matriz en y
  !      DO j = 2, nj
  !         !**************************
  !         !Interpolaciones necesarias
  !         ui = (u_o(i-1,j)+u_o(i,j))/2._DBL
  !         ud = (u_o(i,j)+u_o(i+1,j))/2._DBL
  !         vs = v_o(i,j-1)+d_xuo(i-1)/d2_xuo(i)*(v_o(i+1,j-1)-v_o(i,j-1))
  !         vn = v_o(i,j)+d_xuo(i-1)/d2_xuo(i)*(v_o(i+1,j)-v_o(i,j))
  !         di = d_xuo(i-1)
  !         dd = d_xuo(i)
  !         ds = yo(j)-yo(j-1)
  !         dn = yo(j+1)-yo(j)
  !         gamma_i = 2._DBL * gamma_uo(i-1,j) * gamma_uo(i,j) / (gamma_uo(i-1,j) + gamma_uo(i,j))
  !         gamma_d = 2._DBL * gamma_uo(i+1,j) * gamma_uo(i,j) / (gamma_uo(i+1,j) + gamma_uo(i,j))
  !         gamma_s = 1._DBL/((1._DBL-fey(j)) / gamma_uo(i,j-1) + fey(j) / gamma_uo(i,j))
  !         gamma_n = 1._DBL/((1._DBL-fey(j+1)) / gamma_uo(i,j) + fey(j+1) / gamma_uo(i,j+1))
  !         delta_x = d2_xuo(i)/2._DBL
  !         delta_y = d_yvo(j-1)
  !         temp_int= (1._DBL-di/(2._DBL*delta_x))*temp_o(i,j)+di/(2._DBL*delta_x)*temp_o(i+1,j)
  !         !     (temp_o(i,j)+temp_o(i+1,j))/2._DBL
  !         !************************
  !         alpha  = gamma_i*delta_y/di*DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ui*di/gamma_i))**5._DBL)+&
  !              &DMAX1(0.0_DBL,ui*delta_y)
  !         beta   = gamma_d*delta_y/dd*DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ud*dd/gamma_d))**5._DBL)+&
  !              &DMAX1(0.0_DBL,-ud*delta_y)
  !         AS(j-1)=-(gamma_s*delta_x/ds*DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vs*ds/gamma_s))**5._DBL)+&
  !              &DMAX1(0.0_DBL,vs*delta_x))
  !         AN(j)  =-(gamma_n*delta_x/dn*DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vn*dn/gamma_n))**5._DBL)+&
  !              &DMAX1(0.0_DBL,-vn*delta_x))
  !         ACj(j) = (alpha+beta-AS(j-1)-AN(j)+delta_x*delta_y/dt_o) / rel_vo
  !         Ryj(j) = alpha*u_o(i-1,j)+beta*u_o(i+1,j)-delta_x*delta_y*Ri_o(i,j)*temp_int+&
  !              &delta_x*delta_y*u_anto(i,j)/dt_o+(pres_o(i,j)-pres_o(i+1,j))*delta_y+&
  !              &ACj(j)*(1._DBL-rel_vo)*u_o(i,j)
  !         au_o(i,j) = ACj(j) * rel_vo
  !      END DO
  !      !***********************
  !      !Condiciones de frontera
  !      ACj(nj+1) = 1._DBL
  !      AS(nj)    = 0.0_DBL
  !      Ryj(nj+1) = 0.0_DBL
  !      au_o(i,nj+1) = 1.e40_DBL !ACj(nj+1)
  !      !*************************
  !   end DO
  ! end subroutine ensambla_veluy


  ! !*******************************************************************
  ! !
  ! ! ensambla_velvy
  ! !
  ! ! Subrutina que calcula los coeficientes de la matriz tridiagonal
  ! ! para la velocidad v barriendo la direcci\'on y
  ! !
  ! !*******************************************************************
  ! subroutine ensambla_velvy(xo,fex,d_yvo,d2_yvo,d_xuo,u_o,v_o,v_anto,&
  !      &pres_o,gamma_vo,dt_o,dv_o,av_o,rel_vo,&
  !      )
  !   IMPLICIT NONE
  !   INTEGER :: i,j,k,info
  !   !**********************************
  !   ! Variables del flujo e inc'ognitas
  !   REAL(kind=DBL), DIMENSION(mi+1,nj), INTENT(inout) :: v_o
  !   REAL(kind=DBL), DIMENSION(mi+1,nj), INTENT(out)   :: av_o,dv_o
  !   REAL(kind=DBL), DIMENSION(mi+1,nj), INTENT(in)    :: v_anto,gamma_vo
  !   REAL(kind=DBL), DIMENSION(mi,nj+1), INTENT(in)    :: u_o
  !   REAL(kind=DBL), DIMENSION(mi+1,nj+1), INTENT(in)  :: pres_o
  !   REAL(kind=DBL), DIMENSION(mi+1,nj) :: fv
  !   !**********************************
  !   !Variables de la tridiagonal
  !   REAL(kind=DBL), DIMENSION(mi+1) :: ACi,Rxi
  !   REAL(kind=DBL), DIMENSION(mi)   :: AI,AD
  !   REAL(kind=DBL), DIMENSION(nj)   :: ACj,Ryj
  !   REAL(kind=DBL), DIMENSION(nj-1) :: AS,AN
  !   !*********************
  !   !Variables de interpolaci'on
  !   REAL(kind=DBL) :: ui,ud,vn,vs,di,dd,ds,dn,gamma_i,gamma_d,gamma_s,gamma_n,delta_x,delta_y
  !   !*********************
  !   !Variables de la malla,volumen de control,incremento de tiempo y num Richardson
  !   REAL(kind=DBL), DIMENSION(mi+1), INTENT(in) :: xo,fex
  !   REAL(kind=DBL), DIMENSION(mi-1), INTENT(in) :: d_xuo
  !   REAL(kind=DBL), DIMENSION(nj-1), INTENT(in) :: d_yvo,d2_yvo
  !   REAL(kind=DBL), INTENT(in) :: dt_o,rel_vo
  !   !********************
  !   !Variables auxiliares
  !   REAL(kind=DBL) :: alpha,beta
  !   !*************************
  !   !auxiliar para calcular dv
  !   fv = v_o
  !   !********************
  !   !TDMA en direccion x
  !   !$OMP  PARALLEL DO DEFAULT(NONE)&
  !   !$OMP& PRIVATE(AI,ACi,AD,Rxi,ui,ud,vs,vn,di,dd,ds,dn,gamma_i,gamma_d,gamma_s,gamma_n,&
  !   !$OMP& delta_x,delta_y,alpha,beta,info)&
  !   !$OMP& SHARED(v_o,fv,v_anto,pres_o,gamma_vo,av_o,u_o,d_xuo,d2_yvo,xo,fex,d_yvo,rel_vo,dt_o)
  !   DO j = 2, nj-1
  !      !***********************
  !      !Condiciones de frontera
  !      ACi(1) = 1._DBL
  !      AD(1)  = 0.0_DBL
  !      Rxi(1) = 0.0_DBL
  !      av_o(1,j) = 1.e40_DBL !ACi(1)
  !      !*************************
  !      !Llenado de la matriz en x
  !      DO i = 2, mi
  !         !**************************
  !         !Interpolaciones necesarias
  !         ud = u_o(i,j)+d_yvo(j-1)/d2_yvo(j)*(u_o(i,j+1)-u_o(i,j))
  !         ui = u_o(i-1,j)+d_yvo(j-1)/d2_yvo(j)*(u_o(i-1,j+1)-u_o(i-1,j))
  !         vn = (v_o(i,j+1)+v_o(i,j))/2._DBL
  !         vs = (v_o(i,j)+v_o(i,j-1))/2._DBL
  !         di = xo(i)-xo(i-1)
  !         dd = xo(i+1)-xo(i)
  !         ds = d_yvo(j-1)
  !         dn = d_yvo(j)
  !         gamma_i = 1._DBL/((1._DBL-fex(i)) / gamma_vo(i-1,j) + fex(i) / gamma_vo(i,j))
  !         gamma_d = 1._DBL/((1._DBL-fex(i+1)) / gamma_vo(i,j) + fex(i+1) / gamma_vo(i+1,j))
  !         gamma_s = 2._DBL * gamma_vo(i,j-1) * gamma_vo(i,j) / (gamma_vo(i,j-1) + gamma_vo(i,j))
  !         gamma_n = 2._DBL * gamma_vo(i,j+1) * gamma_vo(i,j) / (gamma_vo(i,j+1) + gamma_vo(i,j))
  !         delta_x = d_xuo(i-1)
  !         delta_y = d2_yvo(j)/2._DBL
  !         !************************
  !         AI(i-1) =-(gamma_i*delta_y/di*DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ui*di/gamma_i))**5._DBL)+&
  !              &DMAX1(0.0_DBL,ui*delta_y))
  !         AD(i)   =-(gamma_d*delta_y/dd*DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ud*dd/gamma_d))**5._DBL)+&
  !              &DMAX1(0.0_DBL,-ud*delta_y))
  !         alpha   = gamma_s*delta_x/ds*DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vs*ds/gamma_s))**5._DBL)+&
  !              &DMAX1(0.0_DBL,vs*delta_x)
  !         beta    = gamma_n*delta_x/dn*DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vn*dn/gamma_n))**5._DBL)+&
  !              &DMAX1(0.0_DBL,-vn*delta_x)
  !         ACi(i)  =(-AI(i-1)-AD(i)+alpha+beta+delta_x*delta_y/dt_o) / rel_vo
  !         Rxi(i)  = alpha*v_o(i,j-1)+beta*v_o(i,j+1)+delta_x*delta_y*v_anto(i,j)/dt_o+&
  !              &(pres_o(i,j)-pres_o(i,j+1))*delta_x+ACi(i)*(1._DBL-rel_vo)*v_o(i,j)
  !         av_o(i,j) = ACi(i) * rel_vo
  !      END DO
  !      !***********************
  !      !Condiciones de frontera
  !      ACi(mi+1) = 1._DBL
  !      AI(mi)    = 0.0_DBL !-1._DBL
  !      Rxi(mi+1) = 0.0_DBL
  !      av_o(mi+1,j) = 1.e40_DBL !ACi(mi+1)
  !      !***************************
  !      CALL tri(AI,ACi,AD,Rxi,mi+1)
  !      DO i = 1, mi+1
  !         v_o(i,j) = Rxi(i)
  !      END DO
  !   END DO
  !   !$OMP END PARALLEL DO
  !   !$OMP  PARALLEL DO DEFAULT(NONE)&
  !   !$OMP& PRIVATE(AS,ACj,AN,Ryj,ui,ud,vs,vn,di,dd,ds,dn,gamma_i,gamma_d,gamma_s,gamma_n,delta_x,delta_y,alpha,beta,info)&
  !   !$OMP& SHARED(v_o,v_anto,pres_o,gamma_vo,av_o,u_o,d_xuo,d2_yvo,xo,fex,d_yvo,rel_vo,dt_o)
  !   DO i = 2, mi
  !      !***********************
  !      !Condiciones de frontera
  !      ACj(1) = 1._DBL
  !      AN(1)  = 0.0_DBL
  !      Ryj(1) = 0.0_DBL
  !      av_o(i,1) = 1.e40_DBL !ACj(1)
  !      !*************************
  !      !Llenado de la matriz en x
  !      DO j = 2, nj-1
  !         !**************************
  !         !Interpolaciones necesarias
  !         ud = u_o(i,j)+d_yvo(j-1)/d2_yvo(j)*(u_o(i,j+1)-u_o(i,j))
  !         ui = u_o(i-1,j)+d_yvo(j-1)/d2_yvo(j)*(u_o(i-1,j+1)-u_o(i-1,j))
  !         vn = (v_o(i,j+1)+v_o(i,j))/2._DBL
  !         vs = (v_o(i,j)+v_o(i,j-1))/2._DBL
  !         di = xo(i)-xo(i-1)
  !         dd = xo(i+1)-xo(i)
  !         ds = d_yvo(j-1)
  !         dn = d_yvo(j)
  !         gamma_i = 1._DBL/((1._DBL-fex(i)) / gamma_vo(i-1,j) + fex(i) / gamma_vo(i,j))
  !         gamma_d = 1._DBL/((1._DBL-fex(i+1)) / gamma_vo(i,j) + fex(i+1) / gamma_vo(i+1,j))
  !         gamma_s = 2._DBL * gamma_vo(i,j-1) * gamma_vo(i,j) / (gamma_vo(i,j-1) + gamma_vo(i,j))
  !         gamma_n = 2._DBL * gamma_vo(i,j+1) * gamma_vo(i,j) / (gamma_vo(i,j+1) + gamma_vo(i,j))
  !         delta_x = d_xuo(i-1)
  !         delta_y = d2_yvo(j)/2._DBL
  !         !************************
  !         alpha   = gamma_i*delta_y/di*DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ui*di/gamma_i))**5._DBL)+&
  !              &DMAX1(0.0_DBL,ui*delta_y)
  !         beta    = gamma_d*delta_y/dd*DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(ud*dd/gamma_d))**5._DBL)+&
  !              &DMAX1(0.0_DBL,-ud*delta_y)
  !         AS(j-1) =-(gamma_s*delta_x/ds*DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vs*ds/gamma_s))**5._DBL)+&
  !              &DMAX1(0.0_DBL,vs*delta_x))
  !         AN(j)   =-(gamma_n*delta_x/dn*DMAX1(0.0_DBL,(1._DBL-0.1_DBL*dabs(vn*dn/gamma_n))**5._DBL)+&
  !              &DMAX1(0.0_DBL,-vn*delta_x))
  !         ACj(j)  =(-AS(j-1)-AN(j)+alpha+beta+delta_x*delta_y/dt_o) / rel_vo
  !         Ryj(j)  = alpha*v_o(i-1,j)+beta*v_o(i+1,j)+delta_x*delta_y*v_anto(i,j)/dt_o+&
  !              &(pres_o(i,j)-pres_o(i,j+1))*delta_x+ACj(j)*(1._DBL-rel_vo)*v_o(i,j)
  !         av_o(i,j) = ACj(j) * rel_vo
  !      END DO
  !      !***********************
  !      !Condiciones de frontera
  !      ACj(nj)  = 1._DBL
  !      AS(nj-1) = 0.0_DBL
  !      Ryj(nj)  = 0.0_DBL
  !      av_o(i,nj) = 1.e40_DBL !ACj(nj) 
  !      !*************************
  !   END DO
  !   !$OMP END PARALLEL DO

  ! end subroutine ensambla_vely
end MODULE ensamblaje
