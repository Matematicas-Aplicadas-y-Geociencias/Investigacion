subroutine sp_integrador_2d(x0, x2, cx, cy, y0, y2, integral)

    use mod_constantes

    implicit none

    real(kind=DBL) :: fun
    real(kind=DBL), intent(out) :: integral
    real(kind=DBL), intent(in) :: x0, x2, cx, cy, y0, y2

    real(kind=DBL) :: a, b, c, d, e
    real(kind=DBL) :: dx1, dx2, dy1, dy2
    real(kind=DBL) :: f0,f1,f2,f3,f4
    real(kind=DBL) :: df0,df1,df2,df3,df4

    dx1 = cx - x0
    dx2 = x2 - cx
    dy1 = cy - y0
    dy2 = y2 - cy

    f1 = fun(dx1, 0.0_DBL)
    f2 = fun(-dx2, 0.0_DBL)
    f0 = fun(0.0_DBL, 0.0_DBL)
    f3 = fun(0.0_DBL, dy1)
    f4 = fun(0.0_DBL, -dy2)

    df1 = f1 - f0
    df2 = f2 - f0
    df0 = f0
    df3 = f3 - f0
    df4 = f4 - f0

    call get_coef(dx1,dx2,dy1,dy2,df0,df1,df2,df3,df4,a,b,c,d,e)

    call integrador_3d(a,b,c,d,e,x0,x2,y0,y2,integral)

end subroutine sp_integrador_2d
