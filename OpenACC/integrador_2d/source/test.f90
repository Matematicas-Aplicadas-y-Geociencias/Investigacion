program test

  use mod_constantes
  use asserts

  implicit none

  integer :: ii
  real(kind=DBL) :: valor_integral
  real(kind=DBL), dimension(3) :: xpo, ypo

  ii = 1

  xpo = [0.0_DBL, 0.5_DBL, 1.0_DBL]
  ypo = [0.0_DBL, 0.5_DBL, 1.0_DBL]

  call integrador_2D(xpo, ypo, ii, valor_integral)

  call assertFloat(2.13935_DBL, valor_integral)
end program test
