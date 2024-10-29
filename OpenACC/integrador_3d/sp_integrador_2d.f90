subroutine sp_integrador_2d(integral)

    use mod_constantes

    implicit none

    ! real(kind=DBL), dimension(NUMERO_PUNTOS), intent(in) :: x
    ! real(kind=DBL), dimension(NUMERO_PUNTOS), intent(in) :: y
    ! real(kind=DBL), dimension(NUMERO_PUNTOS), intent(in) :: f
    real(kind=DBL), intent(out) :: integral

    ! real(kind=DBL), dimension(NUMERO_PUNTOS) :: coef
    real(kind=DBL) :: a, b, c, d, e
    real(kind=DBL) :: dx1, dx2, dy1, dy2
    real(kind=DBL) :: df1,df2,df3,df4,df5
    ! real(kind=DBL) :: termino_A, termino_B, termino_C, termino_D, termino_E

    dx1 = 2.0
    dx2 = 3.0
    dy1 = 4.5
    dy2 = 2.3

    df1 = 8.0
    df2 = 7.0
    df3 = 0.0
    df4 = 5.6
    df5 = 1.1

    call get_coef(dx1,dx2,dy1,dy2,df1,df2,df3,df4,df5,a,b,c,d,e)

    integral = 3.8
end subroutine sp_integrador_2d
