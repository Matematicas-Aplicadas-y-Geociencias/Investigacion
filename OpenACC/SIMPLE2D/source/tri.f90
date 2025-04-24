!
!	flame_tri.f90
!	flame2stepsinterna
!
!	Created by J.E. Barrios on 9/17/06.
!	Copyright 2006 __MyCompanyName__. All rights reserved.
!
subroutine tri(a,b,c,r,n)

implicit none
integer, intent(in) :: n
double precision, intent(in) :: a(n-1),c(n-1)
double precision, intent(inout) :: b(n),r(n)
integer :: i
! eliminacion elementos bajo la matriz
gausselim: do i=2,n
  r(i)=r(i)-(a(i-1)/b(i-1))*r(i-1)
  b(i)=b(i)-(a(i-1)/b(i-1))*c(i-1)
end do gausselim
! solucion para u
r(n)=r(n)/b(n)
sustatras: do i=n-1,1,-1
  r(i)=(r(i)-c(i)*r(i+1))/b(i)
end do sustatras
end subroutine tri
