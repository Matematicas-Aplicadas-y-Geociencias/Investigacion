program magma_lu_example
  use iso_c_binding
  implicit none

  interface
    subroutine magma_init() bind(C, name="magma_init")
    end subroutine magma_init

    subroutine magma_finalize() bind(C, name="magma_finalize")
    end subroutine magma_finalize

    subroutine magma_dgetrf(m, n, a, lda, ipiv, info) bind(C, name="magma_dgetrf")
      use iso_c_binding
      integer(c_int), value :: m, n, lda
      real(c_double) :: a(lda, *)
      integer(c_int) :: ipiv(*)
      integer(c_int) :: info
    end subroutine magma_dgetrf
  end interface


  integer(c_int) :: N, lda, info
  integer(c_int), allocatable :: ipiv(:)
  real(c_double), allocatable :: A(:,:)

  N = 3
  lda = N
  allocate(A(lda,N), ipiv(N))

  ! Matriz a factorizar
  A = reshape([2.0d0, 4.0d0, 6.0d0, &
               3.0d0, 7.0d0, 18.0d0, &
               1.0d0, 2.0d0, -1.0d0], shape(A))

  print *, A
  call magma_init()
  call magma_dgetrf(N, N, A, lda, ipiv, info)

  print *, 'Factorización LU completada.'
  print *, 'Info: ', info
  print *, 'Matriz A después de factorizar:'
  print *, A

  call magma_finalize()
end program magma_lu_example

