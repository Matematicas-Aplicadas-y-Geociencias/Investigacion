Subroutine postprocess_vtk(xo,yo,uo,vo,presso,tempo,bo,archivoo)
  use constantes
  implicit none
  INTEGER :: i,j,k
  REAL(kind=DBL), DIMENSION(mi+1), INTENT(in)           :: xo
  REAL(kind=DBL), DIMENSION(nj+1), INTENT(in)           :: yo
  REAL(kind=DBL), DIMENSION(mi+1,nj+1),   INTENT(in)    :: uo, vo
  REAL(kind=DBL), DIMENSION(mi+1,nj+1),   INTENT(in)    :: tempo,presso
  REAL(kind=DBL), DIMENSION(mi+1,nj+1),   INTENT(in)    :: bo
  CHARACTER(46), INTENT(in)                             :: archivoo
  !************************************
  ! VTK
  open(78, file = trim(archivoo), form = 'formatted')

  write(78,110) '# vtk DataFile Version 2.3'
  write(78,110) '3D Mesh'
  write(78,111) 'ASCII'
  
  write(78,110) 'DATASET STRUCTURED_GRID'
  write(78,120) 'DIMENSIONS', mi+1, nj+1, 1
  write(78,130) 'POINTS', (mi+1)*(nj+1)*(1), ' float'

  do k = 1, 1
     do j = 1, nj+1
        do i = 1, mi+1
           write(78,100) xo(i), yo(j), 0.0
           !write(78,100) dx*float(i-1),dy*float(j-1)
        enddo
     enddo
  end do

  write(78,140) 'POINT_DATA', (mi+1)*(nj+1)*(1)

  write(78,110) 'SCALARS PRESS float'
  write(78,110) 'LOOKUP_TABLE default'
  do k = 1, 1
     do j =1, nj+1
        do i =1, mi+1
           write(78,100) presso(i,j)
        end do
     end do
  end do
  
  write(78,110) 'SCALARS TEMPER float'
  write(78,110) 'LOOKUP_TABLE default'
   
  do k = 1, 1
     do j =1, nj+1
        do i =1, mi+1
           write(78,100) tempo(i,j)
        end do
     end do
  end do
  
  write(78,110) 'VECTORS VELOCITY float'
  do k = 1, 1
     do j =1, nj+1
        do i =1, mi+1
           write(78,100) uo(i,j), vo(i,j), 0.0 
        end do
     end do
  end do
close(78)

100  FORMAT(3(f12.6));
110  FORMAT(A);
111  FORMAT(A,/);
120  FORMAT(A,I4,I4,I4);
130  FORMAT(A,I10,A);
140  FORMAT(A,I10);

END Subroutine postprocess_vtk

