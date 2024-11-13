real(kind=DBL) function fun(x,y)
    use mod_constantes

    implicit none

    real(kind=DBL), intent(in) :: x, y

    ! fun = x**2+y**2
    ! fun = 1 - (x**2+y**2)
    ! fun = y * cos(x) + 2
    ! fun = x*y*exp(x+y)
    ! fun = -x*log(y)
    fun = x + y
end function fun
