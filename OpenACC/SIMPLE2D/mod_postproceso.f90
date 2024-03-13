!
! M\'odulo de postproceso
!
! Este m\'odulo contiene las subrutinas que calculan cantidades 
! de salida y archivos de postproceso
!
!
module postproceso
  use malla, only :  mi, nj, DBL
  implicit none
  
contains
  !*******************************************************************
  !
  ! nusselt_promedio
  !
  ! Subrutina que calcula el n\'umero de nusselt promedio
  ! en paredes verticales
  !
  !*******************************************************************
  subroutine nusselt_promedio_y(&
       &xpo,ypo,deltaxpo,deltaypo,&
       &temp_o,nusselt0_o,nusselt1_o,&
       &i_oo,i_1o&
       &)
    implicit none
    ! $acc routine
    !
    !-------------------------------------
    !
    ! Variables de malla, nusselt y temperatura
    ! 
    real(kind=DBL), dimension(mi+1,nj+1), intent(in) :: temp_o
    real(kind=DBL), dimension(mi+1), intent(in)      :: xpo
    real(kind=DBL), dimension(nj+1), intent(in)      :: ypo
    real(kind=DBL), dimension(mi),   intent(in)      :: deltaxpo
    real(kind=DBL), dimension(nj),   intent(in)      :: deltaypo
    real(kind=DBL), intent(out)                      :: nusselt0_o, nusselt1_o
    integer, intent(in)                              :: i_oo, i_1o
    !
    ! Variables de interpolaci\'on para derivadas
    !
    real(kind=DBL)                  :: a,b,c,dx1,dx2,dy1,dy2
    real(kind=DBL), dimension(mi+1) :: derivada
    !
    ! Variables de interpolaci\'on para integrales
    !
    real(kind=DBL) :: alpha,beta,gamma
    real(kind=DBL) :: x_o,x_1,x_2,y_o,y_1,y_2
    !
    ! Variables auxiliares
    !
    integer :: ii,jj
    !----------------------------------------
    !
    ! C\'alculo para jj = 1 usando interpolaci\'on cuadr\'atica
    !
    derivada = 0.0_DBL
    do ii = i_oo, i_1o
       jj = 1
       dy1 = temp_o(ii,jj+1)-temp_o(ii,jj)
       dy2 = temp_o(ii,jj+2)-temp_o(ii,jj+1)
       dx1 = ypo(jj+1)-ypo(jj)
       dx2 = ypo(jj+2)-ypo(jj+1)
       a = (dy2/dx2-dy1/dx1)/(dx1+dx2)
       b = dy2/dx2-a*(ypo(jj+1)+ypo(jj+2))
       c = temp_o(ii,jj)-a*ypo(jj)*ypo(jj)-b*ypo(jj)
       derivada(ii) = 2._DBL*a*ypo(jj)+b
    end do
    !
    ! Integra con interpolaci\'on cuadr\'atica
    !
    nusselt0_o = 0._DBL
    do ii = i_oo, i_1o, 2
       x_o   = xpo(ii)
       x_1   = xpo(ii+1)
       x_2   = xpo(ii+2)
       y_o   = derivada(ii)
       y_1   = derivada(ii+1)
       y_2   = derivada(ii+2)
       alpha = ((y_2-y_1)/(x_2-x_1)-(y_1-y_o)/(x_1-x_o))/(x_2-x_o)
       beta  = (y_2-y_1)/(x_2-x_1)-alpha*(x_1+x_2)
       gamma = y_o-alpha*x_o*x_o-beta*x_o
       nusselt0_o = nusselt0_o+&
            alpha/3._DBL*(x_2*x_2*x_2-x_o*x_o*x_o)+&
            &beta/2._DBL*(x_2*x_2-x_o*x_o)+&
            &gamma*(x_2-x_o)
    end do
    !
    ! Signo por la ley de Fourier
    !
    nusselt0_o=-nusselt0_o
    !
    !---------------------------------------------------------
    !
    ! C\'alculo para jj = nj usando interpolaci\'on cuadr\'atica
    !
    derivada = 0.0_DBL
    do ii = i_oo, i_1o
       jj = nj-1
       dy1 = temp_o(ii,jj+1)-temp_o(ii,jj)
       dy2 = temp_o(ii,jj+2)-temp_o(ii,jj+1)
       dx1 = ypo(jj+1)-ypo(jj)
       dx2 = ypo(jj+2)-ypo(jj+1)
       a = (dy2/dx2-dy1/dx1)/(dx1+dx2)
       b = dy2/dx2-a*(ypo(jj+1)+ypo(jj+2))
       c = temp_o(ii,jj)-a*ypo(jj)*ypo(jj)-b*ypo(jj)
       derivada(ii) = 2._DBL*a*ypo(jj)+b
    end do
    !
    ! Integra con interpolaci\'on cuadr\'atica
    !
    nusselt0_o = 0._DBL
    do ii = i_oo, i_1o, 2
       x_o   = xpo(ii)
       x_1   = xpo(ii+1)
       x_2   = xpo(ii+2)
       y_o   = derivada(ii)
       y_1   = derivada(ii+1)
       y_2   = derivada(ii+2)
       alpha = ((y_2-y_1)/(x_2-x_1)-(y_1-y_o)/(x_1-x_o))/(x_2-x_o)
       beta  = (y_2-y_1)/(x_2-x_1)-alpha*(x_1+x_2)
       gamma = y_o-alpha*x_o*x_o-beta*x_o
       nusselt1_o = nusselt1_o+&
            alpha/3._DBL*(x_2*x_2*x_2-x_o*x_o*x_o)+&
            &beta/2._DBL*(x_2*x_2-x_o*x_o)+&
            &gamma*(x_2-x_o)
    end do
    !
    ! Signo por la ley de Fourier
    !
    nusselt1_o=-nusselt1_o   
    !
  end subroutine nusselt_promedio_y
  !
end module postproceso
