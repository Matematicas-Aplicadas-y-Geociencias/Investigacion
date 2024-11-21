real(kind=DBL) function fun(x,y)
    use mod_constantes

    implicit none

    real(kind=DBL), intent(in) :: x, y

    ! fun = 3.0_DBL * x*x - 4.5_DBL * y*y
    ! fun = 1 - (x**2+y**2)
    ! fun = y * cos(x) + 2
    ! fun = cos(x)
    ! fun = x*y*exp(x+y)
    ! fun = -x*log(y)
    fun = 0.d0*x*x +2.d0*x+ 0.d0*y*y+3.0_DBL * y
end function fun
