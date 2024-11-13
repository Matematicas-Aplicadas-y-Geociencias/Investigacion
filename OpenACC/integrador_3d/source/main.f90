program main

  use mod_constantes

  implicit none

  real(kind=DBL) :: x0, x2, cx, cy, y0, y2
  real(kind=DBL) :: integral

  x0 = 0.0_DBL
  x2 = 1.0_DBL
  cx = 0.5_DBL
  cy = 0.5_DBL
  y0 = 0.0_DBL
  y2 = 1.0_DBL

  call sp_integrador_2d(x0, x2, cx, cy, y0, y2, integral)

  print '(A, F8.4)', "El valor de la integral es: ", integral

end program main
