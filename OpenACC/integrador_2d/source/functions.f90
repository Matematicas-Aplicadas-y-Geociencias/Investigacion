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

            funE = (0.5+x)*(0.5+y)*exp(x+y)
    end function funE

    real(kind=DBL) function funF(x,y)
        implicit none
        real(kind=DBL), intent(in) :: x, y

        funF = exp(x+y)
    end function funF

    real(kind=DBL) function funG(x,y)
        implicit none
        real(kind=DBL), intent(in) :: x, y

        funG = -y*log(x+0.5)
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

    subroutine get_fun_po(fun,xpo,ypo,fpo,ii)
        implicit none
        integer, intent(in) :: ii
        real(kind=DBL) :: fun
        real(kind=DBL), dimension(3), intent(in) :: xpo, ypo
        real(kind=DBL), dimension(5), intent(out) :: fpo
        real(kind=DBL) :: x1, x2, x3, y1, y2, y3
        real(kind=DBL) :: dx1, dx2, dy1, dy2

        x1 = xpo(ii)
        x2 = xpo(ii + 1)
        x3 = xpo(ii + 2)
        y1 = ypo(ii)
        y2 = ypo(ii + 1)
        y3 = ypo(ii + 2)
        ! -------------------------------------------------------------------------------
        dx2 = x2 - x1
        dx1 = x3 - x2
        dy2 = y2 - y1
        dy1 = y3 - y2
        ! -------------------------------------------------------------------------------
        fpo(1) = fun(0.0_DBL, 0.0_DBL)
        fpo(2) = fun(dx1, 0.0_DBL)
        fpo(3) = fun(-dx2, 0.0_DBL)
        fpo(4) = fun(0.0_DBL, dy1)
        fpo(5) = fun(0.0_DBL, -dy2)
    end subroutine get_fun_po
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
