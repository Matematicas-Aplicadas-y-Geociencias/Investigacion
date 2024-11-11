subroutine integrador_3d(a,b,c,d,e,x0,x2,y0,y2,integral)

    use mod_constantes

    implicit none

    real(kind=DBL), intent(in) :: a, b, c, d, e
    real(kind=DBL), intent(in) :: x0, x2, y0, y2
    real(kind=DBL), intent(out) :: integral

    real(kind=DBL) :: term_a, term_b, term_c, term_d, term_e

    term_a = (a/3) * (x2**3 - x0**3) * (y2 - y0)
    term_b = (b/3) * (x2 - x0) * (y2**3 - y0**3)
    term_c = (c/2) * (x2**2 - x0**2) * (y2 - y0)
    term_d = (d/2) * (x2 - x0) * (y2**2 - y0**2)
    term_e = e * (x2 - x0) * (y2 - y0)

    integral = term_a + term_b + term_c + term_d + term_e
end subroutine integrador_3d
