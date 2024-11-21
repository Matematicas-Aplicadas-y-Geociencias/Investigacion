program main

  use mod_constantes

  implicit none

  integer :: ii
  real(kind=DBL) :: x0, x2, cx, cy, y0, y2
  real(kind=DBL) :: integral
  real(kind=DBL), dimension(3) :: xpo, ypo

  ii = 1

  xpo = [-0.1_DBL, 0.0_DBL, 0.1_DBL]
  ypo = [-0.1_DBL, 0.0_DBL, 0.1_DBL]

  x0 =-0.1_DBL
  x2 = 0.1_DBL
  cx = 0.0_DBL
  cy = 0.0_DBL
  y0 =-0.1_DBL
  y2 = 0.1_DBL

  call integrador_2D(xpo, ypo, ii, integral)

  print *, "El valor de la integral es: ", integral

end program main
