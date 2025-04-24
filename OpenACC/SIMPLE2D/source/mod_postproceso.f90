!
! M\'odulo de postproceso
!
! Este m\'odulo contiene las subrutinas que calculan cantidades 
! de salida y archivos de postproceso
!
! Autor: J.C. Cajas
!
module postproceso
  !
  use malla, only :  mi, nj, DBL
  !
  implicit none
  !
contains
  !
  !*******************************************************************
  !
  ! nusselt_promedio
  !
  ! Subrutina que calcula el n\'umero de nusselt promedio
  ! en paredes verticales
  !
  !*******************************************************************
  subroutine nusselt_promedio_y(&
       &xpo,ypo,deltaxpo,deltaypo,&
       &temp_o,nusselt0_o,nusselt1_o,&
       &i_oo,i_1o&
       &)
    implicit none
    ! $acc routine
    !
    !-------------------------------------
    !
    ! Variables de malla, nusselt y temperatura
    ! 
    real(kind=DBL), dimension(mi+1,nj+1), intent(in) :: temp_o
    real(kind=DBL), dimension(mi+1), intent(in)      :: xpo
    real(kind=DBL), dimension(nj+1), intent(in)      :: ypo
    real(kind=DBL), dimension(mi),   intent(in)      :: deltaxpo
    real(kind=DBL), dimension(nj),   intent(in)      :: deltaypo
    real(kind=DBL), intent(out)                      :: nusselt0_o, nusselt1_o
    integer, intent(in)                              :: i_oo, i_1o
    !
    ! Variables de interpolaci\'on para derivadas
    !
    real(kind=DBL)                  :: a,b,c,dx1,dx2,dy1,dy2
    real(kind=DBL), dimension(mi+1) :: derivada
    !
    ! Variables de interpolaci\'on para integrales
    !
    real(kind=DBL) :: alpha,beta,gamma
    real(kind=DBL) :: x_o,x_1,x_2,y_o,y_1,y_2
    !
    ! Variables auxiliares
    !
    integer :: ii,jj
    !----------------------------------------
    !
    ! C\'alculo para jj = 1 usando interpolaci\'on cuadr\'atica
    !
    derivada = 0.0_DBL
    jj = 1
    do ii = i_oo, i_1o
       ! jj = 1
       dy1 = temp_o(ii,jj+1)-temp_o(ii,jj)
       dy2 = temp_o(ii,jj+2)-temp_o(ii,jj+1)
       dx1 = ypo(jj+1)-ypo(jj)
       dx2 = ypo(jj+2)-ypo(jj+1)
       a = (dy2/dx2-dy1/dx1)/(dx1+dx2)
       b = dy2/dx2-a*(ypo(jj+1)+ypo(jj+2))
       c = temp_o(ii,jj)-a*ypo(jj)*ypo(jj)-b*ypo(jj)
       derivada(ii) = 2._DBL*a*ypo(jj)+b
    end do
    !
    ! Integra con interpolaci\'on cuadr\'atica
    !
    nusselt0_o = 0.0_DBL
    do ii = i_oo, i_1o, 2
       x_o   = xpo(ii)
       x_1   = xpo(ii+1)
       x_2   = xpo(ii+2)
       y_o   = derivada(ii)
       y_1   = derivada(ii+1)
       y_2   = derivada(ii+2)
       alpha = ((y_2-y_1)/(x_2-x_1)-(y_1-y_o)/(x_1-x_o))/(x_2-x_o)
       beta  = (y_2-y_1)/(x_2-x_1)-alpha*(x_1+x_2)
       gamma = y_o-alpha*x_o*x_o-beta*x_o
       nusselt0_o = nusselt0_o+&
            alpha/3._DBL*(x_2*x_2*x_2-x_o*x_o*x_o)+&
            &beta/2._DBL*(x_2*x_2-x_o*x_o)+&
            &gamma*(x_2-x_o)
    end do
    !
    ! Signo por la ley de Fourier
    !
    nusselt0_o=-nusselt0_o
    !
    !---------------------------------------------------------
    !
    ! C\'alculo para jj = nj usando interpolaci\'on cuadr\'atica
    !
    derivada = 0.0_DBL
    jj = nj-1
    do ii = i_oo, i_1o
       ! jj = nj-1
       dy1 = temp_o(ii,jj+1)-temp_o(ii,jj)
       dy2 = temp_o(ii,jj+2)-temp_o(ii,jj+1)
       dx1 = ypo(jj+1)-ypo(jj)
       dx2 = ypo(jj+2)-ypo(jj+1)
       a = (dy2/dx2-dy1/dx1)/(dx1+dx2)
       b = dy2/dx2-a*(ypo(jj+1)+ypo(jj+2))
       c = temp_o(ii,jj)-a*ypo(jj)*ypo(jj)-b*ypo(jj)
       derivada(ii) = 2._DBL*a*ypo(jj+2)+b
    end do
    !
    ! Integra con interpolaci\'on cuadr\'atica
    !
    nusselt1_o = 0.0_DBL
    do ii = i_oo, i_1o, 2
       x_o   = xpo(ii)
       x_1   = xpo(ii+1)
       x_2   = xpo(ii+2)
       y_o   = derivada(ii)
       y_1   = derivada(ii+1)
       y_2   = derivada(ii+2)
       alpha = ((y_2-y_1)/(x_2-x_1)-(y_1-y_o)/(x_1-x_o))/(x_2-x_o)
       beta  = (y_2-y_1)/(x_2-x_1)-alpha*(x_1+x_2)
       gamma = y_o-alpha*x_o*x_o-beta*x_o
       nusselt1_o = nusselt1_o+&
            alpha/3._DBL*(x_2*x_2*x_2-x_o*x_o*x_o)+&
            &beta/2._DBL*(x_2*x_2-x_o*x_o)+&
            &gamma*(x_2-x_o)
    end do
    !
    ! Signo por la ley de Fourier
    !
    nusselt1_o=-nusselt1_o   
    !
  end subroutine nusselt_promedio_y
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
    use malla, only : mi, nj, DBL, mic, njc, zkc
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
  !************************************************************
  !
  ! postprocess_bin
  !
  ! Subrutina de postproceso en formato binario
  !
  !************************************************************
  !
  subroutine postproceso_bin(xuo,yvo,xpo,ypo,&
       &uo,vo,presso,tempo,bo, &
       &Rxc                    &
       )
    use malla, only : mi, nj, DBL, mic, njc, zkc 
    implicit none
    integer :: i,j,k
    real(kind=DBL), DIMENSION(mi), INTENT(in)           :: xuo
    real(kind=DBL), DIMENSION(nj), INTENT(in)           :: yvo
    real(kind=DBL), DIMENSION(mi+1), INTENT(in)         :: xpo
    real(kind=DBL), DIMENSION(nj+1), INTENT(in)         :: ypo   
    real(kind=DBL), DIMENSION(mi,nj+1),   INTENT(in)    :: uo
    real(kind=DBL), DIMENSION(mi+1,nj),   INTENT(in)    :: vo
    real(kind=DBL), DIMENSION(mi+1,nj+1), INTENT(in)    :: tempo,presso
    real(kind=DBL), DIMENSION(mi+1,nj+1), INTENT(in)    :: bo
    character(6),   intent(in)                          :: Rxc
    ! character(64)                                       :: Rxc=repeat(' ',64)
    ! character(46),  INTENT(in)                          :: archivoo
    ! character(64)                                       :: mico,njco,zkco
    !
    ! Creaci\'on de cadenas de caracteres para el contenido de los archivos
    !
    ! write(mico,*) mi+1
    ! write(njco,*) nj+1
    ! write(zkco,*) 1
    ! Rxc = entero_caracter(ceiling(Rx))
    !********************************
    !*** Formato de escritura dat ***
    open(unit=2,file='out_n'//trim(njc)//'m'//trim(mic)//'_R'//trim(Rxc)//'u.bin',access='stream')
    !write(2) placa_min,placa_max,itera_total,ao
    do j = 1, nj+1
       do i = 1, mi
          write(2) xuo(i),ypo(j),uo(i,j)
       end do
    end do
    close(unit=2)
    ! -----------
    open(unit=3,file='out_n'//trim(njc)//'m'//trim(mic)//'_R'//trim(Rxc)//'v.bin',access='stream')
    ! WRITE(3) placa_min,placa_max,itera_total,ao
    DO j = 1, nj
       DO i = 1, mi+1
          WRITE(3) xpo(i),yvo(j),vo(i,j)
       END DO
    END DO
    CLOSE(unit=3)
    ! -----------
    OPEN(unit=4,file='out_n'//trim(njc)//'m'//trim(mic)//'_R'//trim(Rxc)//'p.bin',access='stream')
    ! WRITE(4,*) placa_min,placa_max,itera_total,ao
    DO j = 1, nj+1
       DO i = 1, mi+1
          WRITE(4) tempo(i,j),presso(i,j)
       END DO
    END DO
    CLOSE(unit=4)

  end subroutine postproceso_bin
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

end module postproceso
