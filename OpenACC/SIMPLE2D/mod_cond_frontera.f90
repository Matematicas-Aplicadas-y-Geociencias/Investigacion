!
!
! M\'odulo de herramientas para la lectura de condiciones de frontera del programa SIMPLE 2D
!
! Este m\'odulo contiene herramientas que se utilizan para leer las condiciones de frontera
! en un dominio rectangular
!
MODULE cond_frontera
  !
  implicit none
  !
  use malla, only : mi, nj, DBL
  !
contains
  !
  !-----------------------------------------------------------------------------
  !
  ! Esta subrutina abre los archivos para leer las condiciones de frontera
  ! y crea los arreglos necesarios para imponerlas durante la construcci\'on
  ! de las matrices
  !
  subroutine lectura_cond_frontera(entrada_frontera,&
       &u_anto,v_anto,presso,temp_anto,&
       &xpo,ypo,xuo,yvo,&
       &deltaxpo,deltaypo,&
       &deltaxuo,deltayuo,&
       &deltaxvo,deltayvo,&
       &fexpo,feypo,fexuo,feyvo,&
       &ao,i_0,i_1,iter_inio)
end MODULE cond_frontera
