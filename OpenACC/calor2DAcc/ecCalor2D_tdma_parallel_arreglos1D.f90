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
! double precision, dimension(mi,nj) :: temper,temp_ant
double precision, dimension(mi*nj) :: tempx, tempy, t_ant
double precision, dimension(mi)    :: sourcex, xx
double precision, dimension(nj)    :: sourcey, yy
! double precision, dimension(mi,nj) :: resultx
! double precision, dimension(nj,mi) :: resulty
!
! Matriz a invertir
!
double precision, dimension(mi*nj) :: ax,bx,cx,rx
double precision, dimension(nj*mi) :: ay,by,cy,ry
!
! Se crea la malla 2D
!
a = 10.d0
b = 5.d0
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
cond_ter = 100.d0
temp_ini = 308.d0
temp_fin = 298.d0
flux_aba = 0.0d0
flux_arr = 0.0d0
alpha    = 0.5d0 ! par\'ametro de relajaci\'on
 !
 ! Inicialización de los arreglos a utilizar
 !
 ax = 0.d0
 bx = 0.d0
 cx = 0.d0
 rx = 0.d0
 !
 ay = 0.d0
 by = 0.d0
 cy = 0.d0
 ry = 0.d0
 !
 tempx = (temp_fin+temp_ini)/2.d0
 tempy = (temp_fin+temp_ini)/2.d0 
 t_ant = (temp_fin+temp_ini)/2.d0 
 !
 ! Abrimos la región de datos paralela
 !
 !$acc data copy( tempx(1:mi*nj),tempy(1:mi*nj) ) &
 !$acc & copyin( deltax,deltay,cond_ter,temp_ini,&
 !$acc & temp_fin,flux_aba,flux_arr,alpha, &
 !$acc & t_ant(1:mi*nj),rx(1:mi*nj),ry(1:mi*nj) ) &
 !$acc & create( ax(1:mi*nj),bx(1:mi*nj),cx(1:mi*nj),ay(1:mi*nj),by(1:mi*nj),cy(1:mi*nj) )
 !
 ! Bucle de pseudotiempo
 ! 
 do kk = 1, 5000
    !
    ! Inicia el ciclo que recorre la coordenada y resolviendo
    ! problemas 1D en la dirección de x
    
    !$acc parallel loop gang
    ! gang vector_length(128)
    TDMAX: do jj = 2, nj-1
       !
       ! Ensamblamos matrices en direcci\'on x
       !
       call ensambla_tdmax_1D(ax,bx,cx,rx,deltax,deltay,t_ant,cond_ter,temp_ini,temp_fin,jj)
    end do TDMAX
    !
    ! Llamamos al resolvedor
    !
    !$acc parallel loop gang
    ! gang vector_length(128) 
    SOLX: do jj = 2, nj-1
       call tri(ax((jj-1)*mi+1:(jj-1)*mi+mi),bx((jj-1)*mi+1:(jj-1)*mi+mi),&
            &cx((jj-1)*mi+1:(jj-1)*mi+mi),rx((jj-1)*mi+1:(jj-1)*mi+mi),mi)
       do ii =1, mi
          tempx((jj-1)*mi+ii) = rx((jj-1)*mi+ii)
       end do
    end do SOLX
    !
    ! Se reordena el arreglo tempx para resolver en direcci\'on y
    !
    !$acc parallel loop
    REORDEN: do ii = 1, mi
       do jj = 2, nj-1
          tempy((ii-1)*nj+jj)=tempx((jj-1)*mi+ii)
       end do
    end do REORDEN
    ! Inicia el ciclo que recorre la coordenada x resolviendo
    ! problemas 1D en la dirección de y
    !
    !$acc parallel loop gang
    ! gang vector_length(128)
    TDMAY: do ii = 2, mi-1
       !
       ! Ensamblamos matrices tridiagonales en direcci\'on y
       !
       call ensambla_tdmay_1D(ay,by,cy,ry,deltax,deltay,tempy,cond_ter,flux_aba,flux_arr,ii)
    end do TDMAY
    !$acc parallel loop gang
    ! gang vector_length(128)
    SOLY: do ii=2, mi-1
       !
       ! Llamamos al resolvedor
       !
       call tri(ay((ii-1)*nj+1:(ii-1)*nj+nj),by((ii-1)*nj+1:(ii-1)*nj+nj),&
            &cy((ii-1)*nj+1:(ii-1)*nj+nj),ry((ii-1)*nj+1:(ii-1)*nj+nj),nj)
       do jj =1, nj
          tempy((ii-1)*nj+jj) = ry((ii-1)*nj+jj)
       end do
    end do SOLY
    !
    ! se aplica un esquema de relajaci\'on fija
    !
    ! temper = alpha*temper+(1.d0-alpha)*temp_ant
    !
    ! se actualiza la temperatura de la iteraci\'on anterior
    !
    !$acc parallel loop gang
    ACTUALIZACION: do ii = 2, mi-1
       do jj = 1, nj
          t_ant((jj-1)*mi+ii) = tempy((ii-1)*nj+jj)
       end do
    end do ACTUALIZACION
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
       write(101,*) xx(ii), yy(jj), tempy((ii-1)*nj+jj)
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
