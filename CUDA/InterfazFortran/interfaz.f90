!==================================================================
!Interface to cusolverDn and CUDA C functions
!==================================================================
    
  module cuda_cusolver

  interface

     ! cudaMalloc
     integer (c_int) function cudaMalloc ( buffer, size ) bind (C, name="cudaMalloc" ) 
       use iso_c_binding
       implicit none
       type (c_ptr)  :: buffer
       integer (c_size_t), value :: size
     end function cudaMalloc

     ! cudaMemcpy 
     integer (c_int) function cudaMemcpy ( dst, src, count, kind ) bind (C, name="cudaMemcpy" )
       ! note: cudaMemcpyHostToDevice = 1
       ! note: cudaMemcpyDeviceToHost = 2
       use iso_c_binding
       type (C_PTR), value :: dst, src
       integer (c_size_t), value :: count, kind
     end function cudaMemcpy

     ! cudaFree
     integer (c_int) function cudaFree(buffer)  bind(C, name="cudaFree")
       use iso_c_binding
       implicit none
       type (C_PTR), value :: buffer
     end function cudaFree

     integer (c_int) function cudaMemGetInfo(fre, tot) bind(C, name="cudaMemGetInfo")
       use iso_c_binding
       implicit none
       type(c_ptr),value :: fre
       type(c_ptr),value :: tot
     end function cudaMemGetInfo


     integer(c_int) function cusolverDnCreate(cusolver_Hndl) bind(C,name="cusolverDnCreate")
     
     use iso_c_binding
     implicit none
     
     type(c_ptr)::cusolver_Hndl
     
     end function
     
     integer(c_int) function cusolverDnDestroy(cusolver_Hndl) bind(C,name="cusolverDnDestroy")
     
     use iso_c_binding
     implicit none
     
     type(c_ptr),value::cusolver_Hndl
     
     end function

    integer(c_int) function cusolverDnSgetrf_bufferSize(cusolver_Hndl,m,n,d_A,lda,Lwork) bind(C,name="cusolverDnSgetrf_bufferSize") 
      
      use iso_c_binding
      implicit none
    
      type(c_ptr),value::cusolver_Hndl
      integer(c_int),value::m
      integer(c_int),value::n
      type(c_ptr),value::d_A
      integer(c_int),value::lda
      type(c_ptr),value::Lwork
    end function

     integer(c_int) function cusolverDnSgetrf(cusolver_Hndl,m,n,d_A,lda,d_WS,d_Ipiv,d_devInfo) bind(C, name="cusolverDnSgetrf")

       use iso_c_binding
       implicit none
       type(c_ptr),value::cusolver_Hndl
       integer(c_int),value::m
       integer(c_int),value::n
       type(c_ptr),value::d_A
       integer(c_int),value::lda
       type(c_ptr),value::d_WS
       type(c_ptr),value::d_Ipiv
       type(c_ptr),value::d_devInfo
       
     end function 

     integer (c_int) function cusolverDnSgetrs(cusolver_Hndl,trans,n,nrhs,d_A,lda,d_Ipiv,d_B,ldb,d_devInfo) bind(C, name="cusolverDnSgetrs")

       use iso_c_binding
       implicit none
       type(c_ptr),value::cusolver_Hndl
       integer(c_int), value::trans
       integer(c_int), value::n
       integer(c_int), value::nrhs
       type(c_ptr),value::d_A
       integer(c_int), value::lda    
       type(c_ptr),value::d_Ipiv
       type(c_ptr),value::d_B
       integer(c_int),value::ldb
       type(c_ptr),value::d_devInfo
       
     end function
   
  end interface  
  
end module 

!=========================================
program main
!=========================================

! The purpose of this routine is to provide for GPU based solution of 
! large dense systems of equations in legacy FORTRAN Applications.

! This is the development version of the routine which does not yet work ...

!cudaGetDeviceCount(nDevices)
  
  use iso_c_binding
  use cuda_cusolver
  
  ! ------------- matrix definition & host CPU storage variables 
  
  integer(c_int) m,n,lda,ldb
  
  real, allocatable :: A(:,:)!,ATest(:,:)
  real, allocatable :: B(:,:)
  real, allocatable :: X(:,:)
  integer, allocatable::Ipiv(:)
  integer, target::Lwork
  real, allocatable :: IRN(:), JCN(:), val(:)


!-------------------- CPU equivalents of device variables
  integer devInfo
  integer(c_int) nrhs
  
  integer X_size
  integer*8 Ipiv_size,devInfo_size,Lwork_size
  integer*8 A_size,B_size,Workspace
  
! -------------------- handle to device
  
  type(c_ptr) :: cusolver_Hndl

  
! -------------------- pointers to device memory
  
  type(c_ptr)::d_A
  type(c_ptr)::d_B
  type(c_ptr)::d_Lwork
  type(c_ptr)::d_WS
  type(c_ptr)::d_Ipiv
  type(c_ptr)::d_devInfo

! --------------------- function result variables
  
  integer cusolver_stat
  integer A_mem_stat
  integer B_mem_stat
  integer X_mem_stat
  integer Ipiv_mem_stat
  integer WS_mem_stat
  integer devInfo_mem_stat
  integer Lwork_mem_stat
  
! -------------------- pointers to host CPU memory  
  
  type(c_ptr)::CPU_A_ptr
!  type(c_ptr)::CPU_ATest_ptr
  type(c_ptr)::CPU_B_ptr
  type(c_ptr)::CPU_X_ptr
  type(c_ptr)::CPU_Lwork_ptr
    
  target A,B,X,Ipiv !,ATest
  
  type(c_ptr)::cpfre,cptot
  integer*8,target::free,total
  integer res
  integer*8 cudaMemcpyDeviceToHost, cudaMemcpyHostToDevice
  integer*4 CUBLAS_OP_N, CUBLAS_OP_T
  parameter (cudaMemcpyHostToDevice=1)
  parameter (cudaMemcpyDeviceToHost=2)
  parameter (CUBLAS_OP_N=0)
  parameter (CUBLAS_OP_T=1)
  
  ! ================================================
  free = 0
  total = 0
  res = 1
  cpfre = c_loc(free)
  cptot = c_loc(total)
  
  res = cudaMemGetInfo(cpfre,cptot)
  if (res .ne. 0 ) then
      write (*,*)
      write (*, '(A, I2)') " cudaMemGetInfo error: ", res
      write (*,*)
      stop
  end if
  write (*, '(A, I12)') "  free mem: ", free
  write (*, '(A, I12)') " total mem: ", total
  
!  open(121, file='Input.dat')
!    read(121, *) lda, ldb, nnz
       
  lda = 3
  ldb = 3
  m=lda
  n=ldb
  nrhs = 1
    
  allocate(A(lda,ldb))
!  allocate(ATest(lda,ldb))
  allocate(B(n,1))
  allocate(X(n,1))
  allocate(Ipiv(lda))
  allocate(IRN(nnz))
  allocate(JCN(nnz))
  allocate(val(nnz))
!  allocate(Lwork)
  A_size = SIZEOF(A)
  B_size = SIZEOF(B)
  X_size = SIZEOF(X)
  Ipiv_size = SIZEOF(Ipiv)
  devInfo_size = SIZEOF(devInfo)
  Lwork_size = SIZEOF(Lwork)
  
  print*, val(nnz)
  
  !  "QR FACTORIZATION DENSE LINEAR SOLVER" 
  
  ! -------------- DEFINE [A] AND [B] --------------------------
  
  !           [A]         [x]  =    [B]
  !    | 1.0  2.0  3.0 | |1.0|    | 6.0|
  !    | 4.0  5.0  6.0 | |1.0| =  |15.0|
  !    | 2.0  1.0  1.0 | |1.0|    | 4.0|
  
  A(1,1) = 1.0;   A(1,2) = 2.0;  A(1,3) = 3.0
  A(2,1) = 4.0;   A(2,2) = 5.0;  A(2,3) = 6.0
  A(3,1) = 2.0;   A(3,2) = 1.0;  A(3,3) = 1.0
  
  B(1,1) = 6.0
  B(2,1) = 15.0
  B(3,1) = 4.0
  
!  A(:,:) = 0
!    print *, A(1,7962) 
  
! Define [A]
!  do i=1, nnz
!    read(121, *) IRN(i), JCN(i), val(i)
!    col = IRN(i)
!    row = JCN(i)
!    A(row,col)  = val(i)
!    print *, row," ",col, " ", val(i)
!  enddo
!  close(121)
  
! Define [B]  
!  open(122, file='b.dat')
!  do j=1, lda
!    read(122, *) B(j,1)
!  enddo
!  close(122)
      
  CPU_A_ptr = C_LOC(A)
!  CPU_ATest_ptr = C_LOC(ATest)
  CPU_B_ptr = C_LOC(B)
  CPU_X_ptr = C_LOC(X)
  CPU_Lwork_ptr = C_Loc(Lwork)
  
  A_size = SIZEOF(A)
  B_size = SIZEOF(B)
  X_size = SIZEOF(X)
  Ipiv_size = SIZEOF(Ipiv)
  devInfo_size = SIZEOF(devInfo)
  Lwork_size = SIZEOF(Lwork)
  
  
! Step 1: Create cudense handle ---------------
  cusolver_stat = cusolverDnCreate(cusolver_Hndl)  
  if (cusolver_stat .ne. 0 ) then
      write (*,*)
      write (*, '(A, I2)') " cusolverDnCreate error: ", cusolver_stat
      write (*,*)
      stop
  end if
  
! Step 2: copy A and B to Device
  
  A_mem_stat    = cudaMalloc(d_A,A_size)
  if (A_mem_stat .ne. 0 ) then
      write (*,*)
      write (*, '(A, I2)') " cudaMalloc 1 error: ", A_mem_stat
      write (*,*)
      stop
  end if  
  
  B_mem_stat    = cudaMalloc(d_B,B_size)  
  if (B_mem_stat .ne. 0 ) then
      write (*,*)
      write (*, '(A, I2)') " cudaMalloc 2 error: ", B_mem_stat
      write (*,*)
      stop
  end if
  
! ---------- also allocate space for other device based variables  
  
  Ipiv_mem_stat = cudaMalloc(d_Ipiv,Ipiv_size)
  if (Ipiv_mem_stat .ne. 0 ) then
      write (*,*)
      write (*, '(A, I2)') " cudaMalloc 3 error: ", Ipiv_mem_stat
      write (*,*)
      stop
  end if
  devInfo_mem_stat = cudaMalloc(d_devInfo,devInfo_size)
  if (devInfo_mem_stat .ne. 0 ) then
      write (*,*)
      write (*, '(A, I2)') " cudaMalloc 4 error: ", devInfo_mem_stat
      write (*,*)
      stop
  end if
  Lwork_mem_stat   = cudaMalloc(d_Lwork,Lwork_size)
  if (Lwork_mem_stat .ne. 0 ) then
      write (*,*)
      write (*, '(A, I2)') " cudaMalloc 5 error: ", Lwork_mem_stat
      write (*,*)
      stop
  end if
    
  ! ---------- copy A and B to Device
  
  A_mem_stat = cudaMemcpy(d_A,CPU_A_ptr,A_size,cudaMemcpyHostToDevice)
  if (A_mem_stat .ne. 0 ) then
      write (*,*)
      write (*, '(A, I2)') " cudaMemcpy 1 error: ", A_mem_stat
      write (*,*)
!      stop
  end if
  B_mem_stat = cudaMemcpy(d_B,CPU_B_ptr,B_size,cudaMemcpyHostToDevice)
  if (B_mem_stat .ne. 0 ) then
      write (*,*)
      write (*, '(A, I2)') " cudaMemcpy 2 error: ", B_mem_stat
      write (*,*)
!      stop
  end if
  
!! ---------- check transfer of A matrix by copy back to [ATest]
!  A_mem_stat = cudaMemcpy(CPU_ATest_ptr,d_A,A_size,cudaMemcpyDeviceToHost)  ! Test A matrix memory copy
!  if (A_mem_stat .ne. 0 ) then
!      write (*,*)
!      write (*, '(A, I2)') " cudaMemcpy 3 error: ", A_mem_stat
!      write (*,*)
!      stop
!  end if
  
! Step 3: query working space of Sgetrf (and allocate memory on device)
  
  Lwork = 5
  cusolver_stat =  cusolverDnSgetrf_bufferSize(cusolver_Hndl,m,n,d_A,lda,CPU_Lwork_ptr) 
  if (cusolver_stat .ne. 0 ) then
      write (*,*)
      write (*, '(A, I2)') " DnSgetrf_bufferSize error: ", cusolver_stat
      write (*,*)
!      stop
  end if
 
      write (*,*)
      write (*, '(A, I12)') " Lwork: ", Lwork
      write (*,*)

 
  Workspace = 4*Lwork
  WS_mem_stat = cudaMalloc(d_WS,Workspace)
  if (WS_mem_stat .ne. 0 ) then
      write (*,*)
      write (*, '(A, I2)') " cudaMalloc 6 error: ", WS_mem_stat
      write (*,*)
!      stop
  end if
  
! Step 4: compute LU factorization of [A] 
    
  cusolver_stat = cusolverDnSgetrf(cusolver_Hndl,m,n,d_A,lda,d_WS,d_Ipiv,d_devInfo) 
  if (cusolver_stat .ne. 0 ) then
      write (*,*)
      write (*, '(A, I2)') " cusolverDnSgetrf error: ", WS_mem_stat
      write (*,*)
!      stop
  end if

! Step 5: compute solution vector [X] for Right hand side [B]
  
  cusolver_stat = cusolverDnSgetrs(cusolver_Hndl,CUBLAS_OP_N,n,nrhs,d_A,lda,d_Ipiv,d_B,ldb,d_devInfo)  
  if (cusolver_stat .ne. 0 ) then
      write (*,*)
      write (*, '(A, I2)') " cusolverDnSgetrs error: ", WS_mem_stat
      write (*,*)
!      stop
  end if
  
! Step 6: copy solution vector stored in [B] on device into [X] vector on host
  
  X_mem_stat = cudaMemcpy(CPU_X_ptr,d_B,B_size,cudaMemcpyDeviceToHost)  
  if (X_mem_stat .ne. 0 ) then
      write (*,*)
      write (*, '(A, I2)') " cudaMemcpy 4 error: ", WS_mem_stat
      write (*,*)
!      stop
  end if
  
  do i = 1, n
     print *, X(i,1)
  enddo
  
! step 7: free memory on device and release CPU-side resources
  
  A_mem_Stat    = cudafree(d_A)
  B_mem_Stat    = cudafree(d_B)
  Ipiv_mem_stat = cudafree(d_Ipiv)
  WS_mem_stat   = cudafree(d_WS)
  Lwork_mem_stat = cudafree(d_Lwork)
  
  cusolver_stat = cusolverDnDestroy(cusolver_Hndl)
  
  ! Step 8: deallocate memory on host before exit
  
  deallocate(A)
!  deallocate(ATest)
  deallocate(B)
  deallocate(X)
  deallocate(Ipiv)

  
End Program
