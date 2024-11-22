program main

  use mod_constantes

  implicit none

  integer :: ii
  real(kind=DBL) :: valor_integral
  real(kind=DBL), dimension(3) :: xpo, ypo

  ii = 1

  xpo = [-0.1_DBL, 0.0_DBL, 0.1_DBL]
  ypo = [-0.1_DBL, 0.0_DBL, 0.1_DBL]

  call integrador_2D(xpo, ypo, ii, valor_integral)

  print *, "El valor de la integral es: ", valor_integral

end program main
