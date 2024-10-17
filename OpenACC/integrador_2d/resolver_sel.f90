subroutine resolver_sel(puntos, coef)

  use mod_constantes

  implicit none

  real(kind=DBL), dimension(NUMERO_PUNTOS), intent(in) :: puntos
  real(kind=DBL), dimension(NUMERO_PUNTOS), intent(inout) :: coef

  integer :: i, N, NRHS, LDA, LDB, INFO
  integer, dimension(NUMERO_PUNTOS) :: IPIV
  real(kind=DBL), dimension(NUMERO_PUNTOS, NUMERO_PUNTOS) :: matriz_sistema

  N = NUMERO_PUNTOS
  NRHS = 1
  LDA = N
  LDB = N
  matriz_sistema = 0.0_DBL

  do i = 1, NUMERO_PUNTOS - 2
    matriz_sistema(i, 1) = puntos(i) ** 2
    matriz_sistema(i, 2) = puntos(3) ** 2
    matriz_sistema(i, 3) = puntos(i)
  end do
  do i = 1, NUMERO_PUNTOS - 2
    matriz_sistema(i, 4) = puntos(3)
    matriz_sistema(i, 5) = 1
  end do

  ! Llamada a la rutina DGESV para resolver el sistema Ax = b
  call dgesv(N, NRHS, matriz_sistema, LDA, IPIV, coef, LDB, INFO)

  ! Verificar si la solución fue exitosa
  if (INFO /= 0) then
      print *, "Error: LAPACK dgesv failed with info =", info
  end if

end subroutine resolver_sel
