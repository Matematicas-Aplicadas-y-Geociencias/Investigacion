  !
  ! M\'odulo de inversi\'on de matrices
  !
  ! Este m\'odulo contiene subrutinas para resolver problemas
  ! algebraicos de la forma
  !
  !    Ax=b
  !
module solucionador
  !
  use malla, only : DBL
  !
  implicit none
  !
contains
  !
  !*******************************************************************
  !
  ! tridiagonal
  !
  ! Subrutina que resuelve un sistema tridiagonal de ecuaciones
  ! algebraicas
  !
  !*******************************************************************
  subroutine tridiagonal(a,b,c,r,n)
    !
    !$acc routine seq
    !
    implicit none
    integer, intent(in) :: n
    real(kind=DBL), intent(in)    :: a(n),c(n)
    real(kind=DBL), intent(inout) :: b(n),r(n)
    ! double precision, intent(inout) :: u(n)
    integer :: ii
    
    ! eliminacion elementos bajo la matriz
    gausselim: do ii=2,n
       r(ii)=r(ii)-(a(ii)/b(ii-1))*r(ii-1)
       b(ii)=b(ii)-(a(ii)/b(ii-1))*c(ii-1)
    end do gausselim
    ! solucion para r
    r(n)=r(n)/b(n)
    sustatras: do ii=n-1,1,-1
       r(ii)=(r(ii)-c(ii)*r(ii+1))/b(ii)
    end do sustatras
    
  end subroutine tridiagonal
  
end module solucionador
