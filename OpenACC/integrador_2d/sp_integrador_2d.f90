subroutine sp_integrador_2d(x, y, f, integral)

    use mod_constantes

    implicit none

    real(kind=DBL), dimension(NUMERO_PUNTOS), intent(in) :: x
    real(kind=DBL), dimension(NUMERO_PUNTOS), intent(in) :: y
    real(kind=DBL), dimension(NUMERO_PUNTOS), intent(in) :: f
    real(kind=DBL), intent(out) :: integral

    real(kind=DBL), dimension(NUMERO_PUNTOS) :: coef
    real(kind=DBL) :: a, b, c, d, e
    real(kind=DBL) :: termino_A, termino_B, termino_C, termino_D, termino_E

    coef = f

    call resolver_sel(x, y, coef)

    print *, "Coef: ", coef

    ! Calculo de la integral doble
    a = coef(1) / 3
    b = coef(2) / 3
    c = coef(3) / 2
    d = coef(4) / 2
    e = coef(5)

    ! termino_A = a *
end subroutine sp_integrador_2d
