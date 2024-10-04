subroutine resolver_sel(x, y, coef)
   
  use mod_constantes

  implicit none
  
  integer :: i, N, NRHS, LDA, LDB, INFO
  integer, dimension(NUMERO_PUNTOS) :: IPIV

  real(kind=DBL), dimension(NUMERO_PUNTOS), intent(in) :: x
  real(kind=DBL), dimension(NUMERO_PUNTOS), intent(in) :: y
  real(kind=DBL), dimension(NUMERO_PUNTOS), intent(inout) :: coef


  real(kind=DBL), dimension(NUMERO_PUNTOS, NUMERO_PUNTOS) :: matriz_sistema

  N = NUMERO_PUNTOS
  NRHS = 1
  LDA = N
  LDB = N
  matriz_sistema = 0.0_DBL

  do i = 1, NUMERO_PUNTOS
    matriz_sistema(i, 1) = x(i) * x(i)
    matriz_sistema(i, 2) = y(i) * y(i)
    matriz_sistema(i, 3) = x(i)
    matriz_sistema(i, 4) = y(i)
    matriz_sistema(i, 5) = 1
  end do
    
  ! Llamada a la rutina DGESV para resolver el sistema Ax = b
  call dgesv(N, NRHS, matriz_sistema, LDA, IPIV, coef, LDB, INFO)

  ! Verificar si la solución fue exitosa
  if (INFO /= 0) then
      print *, "Error: LAPACK dgesv failed with info =", info
  end if

end subroutine resolver_sel
