Subroutine postprocess_vtk(xo,yo,zo,uo,vo,wo,presso,tempo,bo,archivoo)
  !
  use malla, only : mi, nj, lk, DBL
  !
  ! use constantes
  !
  implicit none
  INTEGER :: i,j,k
  REAL(kind=DBL), DIMENSION(mi+1), INTENT(in):: xo
  REAL(kind=DBL), DIMENSION(nj+1), INTENT(in):: yo
  REAL(kind=DBL), DIMENSION(lk+1), INTENT(in):: zo
  REAL(kind=DBL), DIMENSION(mi+1,nj+1,lk+1),   INTENT(in)    :: uo, vo, wo
  REAL(kind=DBL), DIMENSION(mi+1,nj+1,lk+1),   INTENT(in)    :: tempo,presso
  REAL(kind=DBL), DIMENSION(lk+1,mi+1,nj+1),   INTENT(in)    :: bo
  CHARACTER(len=64) :: archivoo
  !
  character(64)   :: mico,njco,zkco
  character(128)  :: npuntosc
  !************************************
  !
  ! Creaci\'on de cadenas de caracteres para el contenido de los archivos
  !
  write(mico,*) mi+1
  write(njco,*) nj+1
  write(zkco,*) lk+1
  write(npuntosc,*) (mi+1)*(nj+1)*(lk+1)
  !
  ! VTK
  open(78, file = trim(archivoo), access='stream', convert="big_endian")

  write(78) '# vtk DataFile Version 2.3'//new_line(' ')
  write(78) '3D Mesh'//new_line(' ')
  write(78) 'BINARY'//new_line(' ')
  write(78) 'DATASET STRUCTURED_GRID'//new_line(' ')
  write(78) 'DIMENSIONS '//trim(mico)//trim(njco)//trim(zkco)//new_line('a')
  write(78) 'POINTS '//trim(npuntosc)//' float',new_line('a')
  do k = 1, lk+1
     do j = 1, nj+1
        do i = 1, mi+1
           write(78) real(xo(i)),real(yo(j)),real(zo(k))
        enddo
     enddo
  end do
  write(78) new_line('a')//'POINT_DATA '//trim(npuntosc)
  write(78) 'SCALARS presion float',new_line('a')
  write(78) 'LOOKUP_TABLE default',new_line('a')
  do k = 1, lk+1
     do j =1, nj+1
        do i =1, mi+1
           write(78) real(presso(i,j,k))
        end do
     end do
  end do
  write(78) new_line('a')//'SCALARS temperatura float',new_line('a')
  write(78) 'LOOKUP_TABLE default',new_line('a')
  do k = 1, lk+1
     do j =1, nj+1
        do i =1, mi+1
           write(78) real(tempo(i,j,k))
        end do
     end do
  end do
  write(78) new_line('a')//'VECTORS velocidad float',new_line('a')
  do k = 1, lk+1
     do j =1, nj+1
        do i =1, mi+1
           write(78) real(uo(i,j,k)),real(vo(i,j,k)),real(wo(i,j,k))
        end do
     end do
  end do
  close(78)
! 100 FORMAT(3(f12.6));
! 110 FORMAT(A);
! 111 FORMAT(A,/);
! 120 FORMAT(A,I4,I4,I4);
! 130 FORMAT(A,I10,A);
! 140 FORMAT(A,I10);
  
END Subroutine postprocess_vtk

