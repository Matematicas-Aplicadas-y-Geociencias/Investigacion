subroutine sp_integrador_2d(x, y, f, integral)

  use mod_constantes

  implicit none

  real(kind=DBL), dimension(NUMERO_PUNTOS), intent(in) :: x
  real(kind=DBL), dimension(NUMERO_PUNTOS), intent(in) :: y
  real(kind=DBL), dimension(NUMERO_PUNTOS), intent(in) :: f
  real(kind=DBL), intent(out) :: integral

  real(kind=DBL), dimension(NUMERO_PUNTOS) :: coef

  coef = f

  call resolver_sel(x, y, coef)

  print *, "Coef: ", coef
  
  integral = 3.3
  
end subroutine sp_integrador_2d
