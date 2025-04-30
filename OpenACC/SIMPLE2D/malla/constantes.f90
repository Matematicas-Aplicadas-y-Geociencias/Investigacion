Module constantes
  !
  implicit none
  !
  INTEGER, PARAMETER :: mi=10001, nj=5001, itermax=800000
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
END MODULE constantes
