module functions
    use mod_constantes
    implicit none

contains

    real(kind=DBL) function funA(x,y)
        implicit none
        real(kind=DBL), intent(in) :: x, y

        fun = 3.0_DBL * x * x + 4.5_DBL * y * y
    end function fun

    real(kind=DBL) function funB(x,y)
        implicit none
        real(kind=DBL), intent(in) :: x, y

        fun = 1 - (x**2+y**2)
    end function fun

    real(kind=DBL) function funC(x,y)
        implicit none
        real(kind=DBL), intent(in) :: x, y

        fun = y * cos(x) - 2
    end function fun

    real(kind=DBL) function funD(x,y)
        implicit none
        real(kind=DBL), intent(in) :: x, y

        fun = cos(x*y)
    end function fun

    real(kind=DBL) function funE(x,y)
        implicit none
        real(kind=DBL), intent(in) :: x, y

        fun = x*y*exp(x+y) ! No pasó el test
    end function fun

    real(kind=DBL) function funF(x,y)
        implicit none
        real(kind=DBL), intent(in) :: x, y

        fun = exp(x+y)
    end function fun

    real(kind=DBL) function funG(x,y)
        implicit none
        real(kind=DBL), intent(in) :: x, y

        fun = -y*log(x) ! No pasó el test
    end function fun

    real(kind=DBL) function funH(x,y)
        implicit none
        real(kind=DBL), intent(in) :: x, y

        fun = 0.d0*x*x +2.d0*x+ 0.d0*y*y+3.0_DBL * y
    end function fun

    real(kind=DBL) function funI(x,y)
        implicit none
        real(kind=DBL), intent(in) :: x, y

        fun = y * cos(x) + 2
    end function fun
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
