!
! M\'odulo de postproceso
!
! Este m\'odulo contiene las subrutinas que calculan cantidades
! de salida y archivos de postproceso
!
! Autor: J.C. Cajas, K.A. Figueroa
!
module postproceso
  !
  use malla, only :  mi, nj, DBL
  !
  implicit none
  !
contains
  !
  !************************************************************
  !
  ! postprocesa_parametros
  !
  ! Subrutina de postproceso que muestra los parámetros empleados
  ! y verifica que el directorio esté creado
  !
  !************************************************************
  !
  subroutine postprocesa_parametros(&
       &Ra,&
       &Pr,&
       &dt,&
       &itermax,&
       &paq_itera,&
       &Rin,&
       &rel_pres,&
       &rel_vel,&
       &rel_tem,&
       &conv_u,&
       &conv_t,&
       &conv_p,&
       &conv_resi,&
       &conv_paso,&
       &simpmax,&
       &ecuamax,&
       &entrada_u,&
       &entrada_v,&
       &entrada_w,&
       &entrada_tp,&
       &directorio&
       &)
       !
    use malla, only : mi, nj, DBL, mic, njc
    implicit none
    INTEGER          :: itermax, paq_itera, simpmax, ecuamax
    REAL(kind=DBL)   :: Ra,Pr,dt,Rin,rel_pres,rel_vel,rel_tem
    REAL(kind=DBL)   :: conv_u,conv_p,conv_t,conv_resi,conv_paso
    CHARACTER(len=64):: entrada_u,entrada_v,entrada_w,entrada_tp,directorio
    !
    ! Se escribe la informaci\'on con la que se realiza la ejecuci\'on que
    ! produce los archivos de salida en el directorio nxxxmxxxRxxx/
    ! Esto ayuda a detectar errores de ejecuci\'on por la ausencia de este
    ! directorio
    !
    OPEN(unit=10, file=directorio)
      write (10,*) 'numero de Rayleigh                    ', Ra
      write (10,*) 'numero de Prandtl                     ', Pr
      write (10,*) 'incremento de tiempo                  ', dt
      write (10,*) 'iteraciones maximas                   ', itermax
      write (10,*) 'paquete de iteraciones                ', paq_itera
      write (10,*) 'numero de Richardson                  ', Rin
      write (10,*) 'relajacion de la presion              ', rel_pres
      write (10,*) 'relajacion de la velocidad            ', rel_vel
      write (10,*) 'relajacion de la temperatura          ', rel_tem
      write (10,*) 'convergencia de la velocidad          ', conv_u
      write (10,*) 'convergencia de la temperatura        ', conv_t
      write (10,*) 'convergencia de la presion            ', conv_p
      write (10,*) 'convergencia del residuo              ', conv_resi
      write (10,*) 'convergencia del paso de tiempo       ', conv_paso
      write (10,*) 'iteraciones maximas de SIMPLE         ', simpmax
      write (10,*) 'iteraciones maximas de las ecuaciones ', ecuamax
      write (10,*) 'archivo de entrada para u             ', entrada_u
      write (10,*) 'archivo de entrada para v             ', entrada_v
      write (10,*) 'archivo de entrada para w             ', entrada_w
      write (10,*) 'archivo de entrada para t y p         ', entrada_tp
    CLOSE(unit=10)
    !
    !
  end subroutine postprocesa_parametros
  !
  !************************************************************
  !
  ! postprocess_vtk
  !
  ! Subrutina de postproceso en formato vtk
  !
  !************************************************************
  !
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
    !
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
  !
  end MODULE postproceso
