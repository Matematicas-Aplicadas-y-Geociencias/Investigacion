program test

  use mod_constantes
  use functions
  use asserts

  implicit none

  integer :: ii
  real(kind=DBL) :: valor_integral
  real(kind=DBL), dimension(3) :: xpo, ypo

  ii = 1

  xpo = [0.0_DBL, 0.05_DBL, 0.1_DBL]
  ypo = [0.0_DBL, 0.05_DBL, 0.1_DBL]

  call integrador_2D(funA, xpo, ypo, ii, valor_integral)
  call assertFloat(0.00025_DBL, valor_integral)

  call integrador_2D(funB, xpo, ypo, ii, valor_integral)
  call assertFloat(0.00993333_DBL, valor_integral)

  call integrador_2D(funC, xpo, ypo, ii, valor_integral)
  print *, valor_integral
  call assertFloat(-0.0195008_DBL, valor_integral)

  call integrador_2D(funD, xpo, ypo, ii, valor_integral)
  call assertFloat(0.00025_DBL, valor_integral)
end program test
