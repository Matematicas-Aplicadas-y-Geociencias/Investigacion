module functions
    use mod_constantes
    implicit none

contains

    real(kind=DBL) function funA(x,y)
        implicit none
        real(kind=DBL), intent(in) :: x, y

        funA = 3.0_DBL * x * x + 4.5_DBL * y * y
    end function funA

    real(kind=DBL) function funB(x,y)
        implicit none
        real(kind=DBL), intent(in) :: x, y

        funB = 1 - (x**2+y**2)
    end function funB

    real(kind=DBL) function funC(x,y)
        implicit none
        real(kind=DBL), intent(in) :: x, y

        funC = y * cos(x) - 2
    end function funC

    real(kind=DBL) function funD(x,y)
        implicit none
        real(kind=DBL), intent(in) :: x, y

        funD = cos(x*y)
    end function funD

    real(kind=DBL) function funE(x,y)
        implicit none
        real(kind=DBL), intent(in) :: x, y

        funE = x*y*exp(x+y) ! No pasó el test
    end function funE

    real(kind=DBL) function funF(x,y)
        implicit none
        real(kind=DBL), intent(in) :: x, y

        funF = exp(x+y)
    end function funF

    real(kind=DBL) function funG(x,y)
        implicit none
        real(kind=DBL), intent(in) :: x, y

        funG = -y*log(x) ! No pasó el test
    end function funG

    real(kind=DBL) function funH(x,y)
        implicit none
        real(kind=DBL), intent(in) :: x, y

        funH = 0.d0*x*x +2.d0*x+ 0.d0*y*y+3.0_DBL * y
    end function funH

    real(kind=DBL) function funI(x,y)
        implicit none
        real(kind=DBL), intent(in) :: x, y

        funI = y * cos(x) + 2
    end function funI
end module functions

! real(kind=DBL) function fun(x,y)
!     implicit none

!     real(kind=DBL), intent(in) :: x, y

!     ! fun = 3.0_DBL * x * x + 4.5_DBL * y * y
!     ! fun = 1 - (x**2+y**2)
!     ! fun = y * cos(x) - 2
!     ! fun = cos(x*y)
!     ! fun = x*y*exp(x+y) ! No pasó el test
!     ! fun = exp(x+y)
!     ! fun = -y*log(x) ! No pasó el test
!     ! fun = 0.d0*x*x +2.d0*x+ 0.d0*y*y+3.0_DBL * y
!     ! fun = y * cos(x) + 2
! end function fun
