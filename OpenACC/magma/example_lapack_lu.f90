program magma_lu_example
  implicit none

  integer :: N, lda, info
  integer, allocatable :: ipiv(:)
  real(kind=8), allocatable :: A(:,:)

  N = 3
  lda = N
  allocate(A(lda,N), ipiv(N))

  ! Matriz a factorizar
  A = reshape([2.0d0, 4.0d0, 6.0d0, &
               3.0d0, 7.0d0, 18.0d0, &
               1.0d0, 2.0d0, -1.0d0], shape(A))

  print *, A
  call dgetrf(N, N, A, lda, ipiv, info)

  print *, 'Factorización LU completada.'
  print *, 'Info: ', info
  print *, 'Matriz A después de factorizar:'
  print *, A

end program magma_lu_example

