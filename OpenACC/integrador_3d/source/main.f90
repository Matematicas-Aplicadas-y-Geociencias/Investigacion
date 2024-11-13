program main

  use mod_constantes

  implicit none

  real(kind=DBL) :: x0, x2, cx, cy, y0, y2
  real(kind=DBL) :: integral

  x0 =-0.1_DBL
  x2 = 0.1_DBL
  cx = 0.0_DBL
  cy = 0.0_DBL
  y0 =-0.1_DBL
  y2 = 0.1_DBL

  call sp_integrador_2d(x0, x2, cx, cy, y0, y2, integral)

  print '(A, F8.4)', "El valor de la integral es: ", integral

end program main
