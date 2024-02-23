PROGRAM Calor2D
!
use constantes
implicit none !(type, external)
! external :: sgesv
integer :: ii, jj, kk
!
! Variables de tamaño para la placa
double precision :: a, b, deltax, deltay
!
! Parámetros físicos del problema
double precision :: cond_ter, temp_ini, temp_fin
double precision :: flux_aba, flux_arr, alpha
!
! Arreglos para incógnitas, término fuente y posición
double precision, dimension(mi,nj) :: temper,temp_ant
double precision, dimension(mi)    :: sourcex, xx
double precision, dimension(nj)    :: sourcey, yy
double precision, dimension(mi,nj) :: resultx
double precision, dimension(nj,mi) :: resulty,tempy
!
! Matriz a invertir
double precision, dimension(mi,nj) :: AI,AD
double precision, dimension(mi,nj) :: AC
double precision, dimension(nj,mi) :: BI,BD
double precision, dimension(nj,mi) :: BC
!
! Se crea la malla 2D
!
a = 1.d0
b = 5.d-1
do ii = 1, mi
    xx(ii) = a*dfloat(ii-1)/dfloat(mi-1)
!     print*, ii, xx(ii)
end do
do jj = 1, nj
    yy(jj) = b*dfloat(jj-1)/dfloat(nj-1)
end do
deltax = 1.d0/dfloat(mi-1)
deltay = 1.d0/dfloat(nj-1)
!
! Se definen los parámetros físicos del problema
!
cond_ter = 1.d0
temp_ini = 308.d0
temp_fin = 298.d0
flux_aba = 5.0d0
flux_arr = 5.0d0
alpha    = 0.5d0 ! par\'ametro de relajaci\'on
!
! Inicialización de los arreglos a utilizar
!
AI = 0.d0; AC = 0.d0; AD = 0.d0
BI = 0.d0; BC = 0.d0; BD = 0.d0
temper = (temp_fin+temp_ini)/2.d0
temp_ant = temper
resultx = 0.0d0
resulty = 0.0d0
!
! Abrimos la región de datos paralela
!
!$acc data copy( temper(1:mi,1:nj) ) &
!$acc & copyin( deltax,deltay,temp_ant(1:mi,1:nj),cond_ter,temp_ini,&
!$acc & temp_fin,flux_aba,flux_arr,alpha,resultx(1:mi,1:nj),resulty(1:nj,1:mi),&
!$acc & AI(1:mi,1:nj),AC(1:mi,1:nj),AD(1:mi,1:nj),&
!$acc & BI(1:nj,1:mi),BC(1:nj,1:mi),BD(1:nj,1:mi),tempy(1:nj,1:mi) )
!
! Bucle de pseudotiempo
! 
do kk = 1, 500
   !
   ! Inicia el ciclo que recorre la coordenada y resolviendo
   ! problemas 1D en la dirección de x
   
   !$acc parallel loop
   do jj = 2, nj-1
      !
      ! Ensamblamos matrices en direcci\'on x
      !
      call ensambla_tdmax(AI,AC,AD,resultx,deltax,deltay,temp_ant,cond_ter,temp_ini,temp_fin,jj)
      !
      ! Llamamos al TDMA
      call tri(AI(1:mi,jj),AC(1:mi,jj),AD(1:mi,jj),resultx(1:mi,jj),temper(1:mi,jj),mi)
      
    end do
    !
    ! Inicia el ciclo que recorre la coordenada x resolviendo
    ! problemas 1D en la dirección de y
    !
    !$acc parallel loop
    do ii = 2, mi-1
       !
       ! Ensamblamos matrices tridiagonales en direcci\'on y
       !
       call ensambla_tdmay(BI,BC,BD,resulty,deltax,deltay,temper,cond_ter,flux_aba,flux_arr,ii)
       !
       ! Llamamos al TDMA
       !
       call tri(BI(1:nj,ii),BC(1:nj,ii),BD(1:nj,ii),resulty(1:nj,ii),tempy(1:nj,ii),nj)
       do jj =1, nj
          temper(ii,jj) = tempy(jj,ii)
       end do
       
    end do
    !
    ! se aplica un esquema de relajaci\'on fija
    !
    temper = alpha*temper+(1.d0-alpha)*temp_ant
    !
    ! se actualiza la temperatura de la iteraci\'on anterior
    !
    temp_ant = temper
    !
    !print*, "Residuo: ", kk, maxval(temper-temp_ant)
    !     temp_ant = temper
    !
 end do
 !
 ! Cerramos la región de datos paralela
 !
 !$acc end data
 !
 !! Escritura de resultados
 !
 do ii = 1, mi
    do jj = 1, nj
       write(101,*) xx(ii), yy(jj), temper(ii,jj)
    end do
    write(101,*)
 end do
 !
 ! Obtenemos la diferencia entre la solución analítica y la solución numérica, así como el error porcentual
 !do ii =2, mi-1
 !    do jj = 2, nj-1
 !        !
 !        write(102,*) xx(ii), yy(jj), temper(ii,jj)-(-xx(ii)+308), (temper(ii,jj)-(-xx(ii)+308))/(-xx(ii)+308)
 !    end do
 !        write(102,*)
 !end do
 !
 !
 ! Obtenemos la diferencia entre la solución analítica y la solución numérica, así como el error porcentual
 ! do ii =2, mi-1
 !     do jj = 2, nj-1
 !         !
 !         write(102,*) xx(ii), yy(jj), temper(ii,jj)-(-xx(ii)+308), (temper(ii,jj)-(-xx(ii)+308))/(temp_fin-temp_ini)
 !     end do
 !         write(102,*)
 ! end do
 !
END PROGRAM Calor2D
