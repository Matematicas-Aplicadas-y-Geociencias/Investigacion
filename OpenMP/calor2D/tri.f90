!
!	flame_tri.f90
!	flame2stepsinterna
!
!	Created by J.E. Barrios on 9/17/06.
!
subroutine tri(a,b,c,r,u,n)
implicit none
integer, intent(in) :: n
double precision, intent(in) :: a(n),c(n)
double precision, intent(inout) :: b(n),r(n)
double precision, intent(inout) :: u(n)
integer :: i

u=0.d0
! eliminacion elementos bajo la matriz
gausselim: do i=2,n
  r(i)=r(i)-(a(i)/b(i-1))*r(i-1)
  b(i)=b(i)-(a(i)/b(i-1))*c(i-1)
end do gausselim
! solucion para u
u(n)=r(n)/b(n)
sustatras: do i=n-1,1,-1
  u(i)=(r(i)-c(i)*u(i+1))/b(i)
end do sustatras
end subroutine tri
