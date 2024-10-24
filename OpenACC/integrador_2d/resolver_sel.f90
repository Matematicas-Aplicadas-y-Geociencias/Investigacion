subroutine resolver_sel(puntos_funcion, coeficientes_matriz)

    use mod_constantes

    implicit none

    real(kind=DBL), dimension(NUMERO_PUNTOS), intent(in) :: puntos_funcion
    real(kind=DBL), dimension(NUMERO_PUNTOS), intent(inout) :: coeficientes_matriz

    integer :: i, N, NRHS, LDA, LDB, INFO
    integer, dimension(NUMERO_PUNTOS) :: IPIV
    real(kind=DBL), dimension(NUMERO_PUNTOS, NUMERO_PUNTOS) :: matriz_sistema

    N = NUMERO_PUNTOS
    NRHS = 1
    LDA = N
    LDB = N
    matriz_sistema = 0.0_DBL

    do i = 1, NUMERO_PUNTOS - 3
        matriz_sistema(i, 1) = puntos_funcion(i) ** 2
        matriz_sistema(i, 2) = puntos_funcion(3) ** 2
    end do
    matriz_sistema(i, 3) = puntos_funcion(i)
    do i = 1, NUMERO_PUNTOS - 2
        matriz_sistema(i, 4) = puntos_funcion(3)
        matriz_sistema(i, 5) = 1
    end do

    ! Llamada a la rutina DGESV para resolver el sistema Ax = b
    call dgesv(N, NRHS, matriz_sistema, LDA, IPIV, coeficientes_matriz, LDB, INFO)

    ! Verificar si la solución fue exitosa
    if (INFO /= 0) then
        print *, "Error: LAPACK dgesv failed with info =", info
    end if

end subroutine resolver_sel
