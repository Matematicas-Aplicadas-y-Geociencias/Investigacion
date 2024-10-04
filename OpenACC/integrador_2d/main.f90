program main

  use mod_constantes

  implicit none

  real(kind=DBL) :: integral

  real(kind=DBL), dimension(5) :: x
  real(kind=DBL), dimension(5) :: y
  real(kind=DBL), dimension(5) :: f

  x = [5.0_DBL, 5.0_DBL, 5.0_DBL, 3.0_DBL, 8.0_DBL]
  y = [11.0_DBL, 7.0_DBL, 4.0_DBL, 7.0_DBL, 7.0_DBL]
  f = [5.4_DBL, 7.0_DBL, 1.2_DBL, 1.1_DBL, 6.0_DBL]

  call sp_integrador_2d(x, y, f, integral)

  print '(A, F6.4)', "El valor de la integral es: ", integral

end program main
