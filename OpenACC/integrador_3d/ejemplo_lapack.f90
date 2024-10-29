program solve_system
    implicit none

    ! Declaración de variables
    integer :: n, lda, info, ipiv(2)
    double precision :: a(2,2), b(2)
    
    ! Inicializar matriz A y vector B
    a = reshape([ 1.0d0, 2.0d0, 3.0d0, 4.0d0 ], shape(a))
    b = [5.0d0, 6.0d0]

    ! Dimensiones de la matriz
    n = 2
    lda = 2

    ! Llamada a la rutina DGESV para resolver el sistema Ax = b
    call dgesv(n, 1, a, lda, ipiv, b, lda, info)

    ! Verificar si la solución fue exitosa
    if (info == 0) then
        print *, "Solution: ", b
    else
        print *, "Error: LAPACK dgesv failed with info =", info
    end if
end program solve_system
