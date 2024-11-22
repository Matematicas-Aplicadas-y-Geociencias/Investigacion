subroutine integrador_2D(xpo, ypo, ii, valor_integral)
    ! TODO
    ! - Buscar mejores nombres para las variables
    use mod_constantes

    implicit none

    ! -------------------------------------------------------------------------------
    ! Declaracion de Variables ---------------------------------------------------------------------
    ! -------------------------------------------------------------------------------
    integer, intent(in) :: ii
    real(kind=DBL) :: fun
    real(kind=DBL), dimension(3), intent(in) :: xpo, ypo
    real(kind=DBL), intent(out) :: valor_integral

    real(kind=DBL) :: x1, x2, x3, y1, y2, y3
    real(kind=DBL) :: a, b, c, d, e
    real(kind=DBL) :: dx1, dx2, dy1, dy2
    real(kind=DBL) :: f0,f1,f2,f3,f4
    real(kind=DBL) :: df0,df1,df2,df3,df4
    real(kind=DBL) :: numerador, denominador
    real(kind=DBL) :: term_1, term_2, term_3, term_4, term_5
    ! -------------------------------------------------------------------------------
    ! Inicializacion de Variables ---------------------------------------------------
    ! -------------------------------------------------------------------------------
    x1 = xpo(ii)
    x2 = xpo(ii + 1)
    x3 = xpo(ii + 2)
    y1 = ypo(ii)
    y2 = ypo(ii + 1)
    y3 = ypo(ii + 2)

    dx2 = x2 - x1
    dx1 = x3 - x2
    dy2 = y2 - y1
    dy1 = y3 - y2

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
    ! print*, "DEBUG: funs = ", f1, f2, f0
    ! -------------------------------------------------------------------------------
    ! Calcular coeficientes ---------------------------------------------------------
    ! -------------------------------------------------------------------------------
    ! Coeficiente a
    numerador = dx2 * df1 + dx1 * df2
    denominador = dx1 * dx2 * (dx1 + dx2)
    a = numerador / denominador

    ! print*, "DEBUG: num den a = ", numerador, denominador, a

    ! Coeficiente c
    numerador = a * dx2 * dx2 - df2
    denominador = dx2
    c = numerador / denominador

    ! Coeficiente b
    numerador = dy2 * df3 + dy1 * df4
    denominador = dy1 * dy2 * (dy1 + dy2)
    b = numerador / denominador

    ! Coeficiente d
    numerador = b * dy2 * dy2 - df4
    denominador = dy2
    d = numerador / denominador

    e = df0

    print *, 'Coef a:', a
    print *, 'Coef b:', b
    print *, 'Coef c:', c
    print *, 'Coef d:', d
    print *, 'Coef e:', e
    ! -------------------------------------------------------------------------------
    ! Calcular el valor de la integral -------------------------------------------------------------
    ! -------------------------------------------------------------------------------
    term_1 = (a/3) * (x3 * x3 * x3 - x1 * x1 * x1) * (y3 - y1)
    term_2 = (b/3) * (x3 - x1) * (y3 * y3 *y3 - y1 * y1 * y1)
    term_3 = (c/2) * (x3 * x3 - x1 * x1) * (y3 - y1)
    term_4 = (d/2) * (x3 - x1) * (y3 * y3 - y1 * y1)
    term_5 = e * (x3 - x1) * (y3 - y1)

    valor_integral = term_1 + term_2 + term_3 + term_4 + term_5
    ! -------------------------------------------------------------------------------
end subroutine integrador_2D
