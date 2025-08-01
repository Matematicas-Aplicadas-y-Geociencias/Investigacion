program main

  use mod_constantes
  use functions

  implicit none

  integer :: ii
  real(kind=DBL) :: valor_integral
  real(kind=DBL), dimension(3) :: xpo, ypo
  real(kind=DBL), dimension(5) :: fpo

  ii = 1

  xpo = [-0.05_DBL, 0.00_DBL, 0.05_DBL]
  ypo = [-0.05_DBL, 0.00_DBL, 0.05_DBL]

  call integrador_2D(funI, xpo, ypo, ii, valor_integral)
  print *, "El valor de la integral es: ", valor_integral

  call get_fun_po(funI, xpo, ypo, fpo, ii)
  call integrate_2D(xpo, ypo, fpo, ii, valor_integral)
  print *, "El valor de la integral es: ", valor_integral

end program main
