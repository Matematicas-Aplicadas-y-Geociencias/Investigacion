subroutine sp_integrador_2d(integral)

    use mod_constantes

    implicit none

    real(kind=DBL), intent(out) :: integral

    real(kind=DBL) :: x0, x2, y0, y2
    real(kind=DBL) :: a, b, c, d, e
    real(kind=DBL) :: dx1, dx2, dy1, dy2
    real(kind=DBL) :: df1,df2,df3,df4,df5

    x0 = -1!1
    x2 = 1!2
    y0 = -1!3
    y2 = 1!4

    dx1 = 1
    dx2 = 1
    dy1 = 1
    dy2 = 1

    df1 = -1
    df2 = -1
    df3 = 1
    df4 = -1
    df5 = -1

    call get_coef(dx1,dx2,dy1,dy2,df1,df2,df3,df4,df5,a,b,c,d,e)

    call integrador_3d(a,b,c,d,e,x0,x2,y0,y2,integral)
    print *,integral
end subroutine sp_integrador_2d
