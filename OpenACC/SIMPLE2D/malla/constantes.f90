Module constantes
  !
  implicit none
  !
  INTEGER, PARAMETER :: mi=101, nj=501, itermax=800000
  INTEGER, PARAMETER :: SGL=SELECTED_REAL_KIND(P=6,R=37)
  INTEGER, PARAMETER :: DBL=SELECTED_REAL_KIND(P=15,R=300)
  !
contains
  !
  !************************************************************
  !
  ! entero_caracter
  !
  ! Subrutina que devuelve una cadena a partir de un entero
  !
  !************************************************************
  !
  function entero_caracter(entero)
    
    implicit none
    
    character(6)        :: entero_caracter 
    integer, intent(in) :: entero

    integer             :: uni, dec, cen, mil, dmi
    character(1)        :: un,  de,  ce,  mi,  dm

    dmi = entero/10000
    mil = ( entero-dmi*10000 ) / 1000
    cen = ( entero-dmi*10000-mil*1000 ) / 100
    dec = ( entero-dmi*10000-mil*1000-cen*100 ) / 10
    uni = ( entero-dmi*10000-mil*1000-cen*100-dec*10 )

    write(un,16) uni; 16 format(I1)
    write(de,16) dec
    write(ce,16) cen
    write(mi,16) mil
    write(dm,16) dmi

    entero_caracter = dm//mi//ce//de//un

  end function entero_caracter
  !
    !************************************************************
  !
  ! postprocess_vtk
  !
  ! Subrutina de postproceso en formato vtk (archivos binarios)
  !
  !************************************************************
  !
  subroutine postproceso_vtk(&
       &xo,yo,uo,vo,presso,tempo,bo,archivoo&
       &)
    !use malla, only : mi, nj, DBL, mic, njc, zkc
    implicit none
    INTEGER :: i,j,k
    REAL(kind=DBL), DIMENSION(mi+1), INTENT(in)           :: xo
    REAL(kind=DBL), DIMENSION(nj+1), INTENT(in)           :: yo
    REAL(kind=DBL), DIMENSION(mi+1,nj+1),   INTENT(in)    :: uo, vo
    REAL(kind=DBL), DIMENSION(mi+1,nj+1),   INTENT(in)    :: tempo,presso
    REAL(kind=DBL), DIMENSION(mi+1,nj+1),   INTENT(in)    :: bo
    CHARACTER(46), INTENT(in)                             :: archivoo
    character(64)                                         :: mico,njco,zkco
    character(128)                                        :: npuntosc
    !
    ! Creaci\'on de cadenas de caracteres para el contenido de los archivos
    !
    write(mico,*) mi+1
    write(njco,*) nj+1
    write(zkco,*) 1
    write(npuntosc,*) (mi+1)*(nj+1)*1
    !
    !************************************
    ! VTK
    open(78, file = trim(archivoo), access='stream', convert="big_endian")

    write(78) '# vtk DataFile Version 2.3'//new_line(' ')
    write(78) '3D Mesh'//new_line(' ')
    write(78) 'BINARY'//new_line(' ')
    write(78) 'DATASET STRUCTURED_GRID'//new_line(' ')
    write(78) 'DIMENSIONS '//trim(mico)//trim(njco)//trim(zkco)//new_line('a')
    write(78) 'POINTS '//trim(npuntosc)//' float',new_line('a')
    do k = 1, 1
       do j = 1, nj+1
          do i = 1, mi+1
             write(78) real(xo(i)),real(yo(j)),0.0
          enddo
       enddo
    end do
    write(78) new_line('a')//'POINT_DATA '//trim(npuntosc)
    write(78) 'SCALARS PRESS float',new_line('a')
    write(78) 'LOOKUP_TABLE default',new_line('a')
    do k = 1, 1
       do j =1, nj+1
          do i =1, mi+1
             write(78) real(presso(i,j))
          end do
       end do
    end do
    write(78) new_line('a')//'SCALARS TEMPER float',new_line('a')
    write(78) 'LOOKUP_TABLE default',new_line('a')
    do k = 1, 1
       do j =1, nj+1
          do i =1, mi+1
             write(78) real(tempo(i,j))
          end do
       end do
    end do
    write(78) new_line('a')//'VECTORS VELOCITY float',new_line('a')
    do k = 1, 1
       do j =1, nj+1
          do i =1, mi+1
             write(78) real(uo(i,j)),real(vo(i,j)),0.0 
          end do
       end do
    end do
    close(78)
    !
    ! 100 FORMAT(3(f12.6));
    ! 110 FORMAT(A);
    ! 111 FORMAT(A,/);
    ! 120 FORMAT(A,I4,I4,I4);
    ! 130 FORMAT(A,I10,A);
    ! 140 FORMAT(A,I10);
    !
  end subroutine postproceso_vtk
  !
END MODULE constantes
