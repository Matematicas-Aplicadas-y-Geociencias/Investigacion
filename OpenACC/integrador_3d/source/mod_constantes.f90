module mod_constantes

  implicit none

  integer, parameter :: DBL=SELECTED_REAL_KIND(P=15,R=307)
  real(kind=DBL), parameter :: PI = 4.0_DBL * atan (1.0_DBL)

end module mod_constantes
