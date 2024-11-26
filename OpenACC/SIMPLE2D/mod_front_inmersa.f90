!
! M\'odulo de frontera  inmersa
!
! Este m\'odulo contiene las subrutinas que permiten incluir s\'olidos en el dominio del fluido 
! Se utiliza la media geom\'etrica de los coeficientes de difusi\'on en las caras de los
! vol\'umenes de control como se describe en Patankar 1987.
!
! Autor: J.C. Cajas
!
module frontera_inmersa
  !
  use malla, only : mi, nj, DBL
  !
  implicit none
  !
contains
  !
  !*******************************************************************
  !
  ! definir_cuerpo
  !
  ! Subrutina que define si un punto en la malla est\'a dentro de un
  ! cuerpo s\'olido o no. Asigna los valores correspondientes a las
  ! constantes de difusi\'on gamma.
  !
  !*******************************************************************
  !
  subroutine definir_cuerpo(gamma_momeno, gamma_enero, dominiox, funcion)
    
  end subroutine definir_cuerpo
end module frontera_inmersa
