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
    INTEGER :: itermax, paq_itera, simpmax, ecuamax
    REAL(kind=DBL)   :: Ra,Pr,dt,Rin,rel_pres,rel_vel,rel_tem
    REAL(kind=DBL)   :: conv_u,conv_p,conv_t,conv_resi,conv_paso
    CHARACTER(len=64):: entrada_u,entrada_v,entrada_w,entrada_tp,directorio
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
end MODULE postproceso
