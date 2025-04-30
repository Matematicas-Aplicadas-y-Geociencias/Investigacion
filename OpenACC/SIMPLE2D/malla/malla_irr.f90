Program malla_irr
use constantes
Implicit none
Integer i,j,k,i_o,i_1
!*****************************************
!Variables de la malla, volumen de control
REAL(kind=DBL) xu(mi),yv(nj),ao,u(mi,nj+1)
REAL(kind=DBL) xu_eta(mi),yv_eta(nj),x_eta(mi+1),y_eta(nj+1)
REAL(kind=DBL) my,ky,mx,kx,tau,b,xc,h,beta
REAL(kind=DBL) x_placa_min, x_placa_max
character(len=3) :: njc,mic
!********************************
if (nj<100) then
  write(njc,170) int(nj);170 format(I2)
  njc='0'//njc
else
  write(njc,160) int(nj);160 format(I3)
endif
if (mi<100) then
  write(mic,170) int(mi)
  mic='0'//mic
else
  write(mic,160) int(mi)
endif
!*********************************
!generaci'on de la malla irregular
!puntos en x
ao=10.0_DBL
xc=8.5_DBL
x_placa_min=8.0_DBL
x_placa_max=9.0_DBL
do i=1,mi
  xu(i)=dfloat(i-1)/dfloat(mi-1)
end do
tau=2._DBL
b=1._DBL/(2._DBL*tau)*dlog((1._DBL+(dexp(tau)-1._DBL)*(xc/ao))/(1._DBL+(dexp(-tau)-1._DBL)*(xc/ao)))
do i =1,mi
  xu_eta(i)=xc*(1._DBL+dsinh(tau*(xu(i)-b))/dsinh(tau*b))
end do
! kx=10._DBL
! mx=1.0_DBL!empaquetamiento central en dirección x,menor coeficiente, mayor empaquetamiento
! do i=1,mi
!   xu_eta(i)=ao*(0.5_DBL+mx*(xu(i)-0.5_DBL)+0.5_DBL*(1-mx)*dsinh(kx*(xu(i)-0.5_DBL))/dsinh(0.5_DBL*kx))
! end do
!************************
x_eta(1)=(3._DBL*xu_eta(1)-xu_eta(2))/2._DBL
! x_eta(1)=xu_eta(1)
do i=1,mi-1
  x_eta(i+1)=(xu_eta(i+1)+xu_eta(i))/2._DBL
  if(x_eta(i)<x_placa_min)then
    i_o=i
  end if
  if(x_eta(i)<=x_placa_max)then
    i_1=i
  end if 
end do
x_eta(mi+1)=(3._DBL*xu_eta(mi)-xu_eta(mi-1))/2._DBL
! x_eta(mi+1)=xu_eta(mi)
i_o=i_o+1
! i_1 = mi + 1
!***********************
!***********************
!puntos en y
do j=1,nj
  yv(j)=dfloat(j-1)/dfloat(nj-1)
end do
!
! -----------------------------------------------------------------------------
!
ky=10._DBL
my=0.4_DBL !empaquetamiento en las paredes,menor coeficiente,mayor empaquetamiento
do j=1,nj
   yv_eta(j)=1.0_DBL*(0.5_DBL+my*(yv(j)-0.5_DBL)+0.5_DBL*(1-my)*derf(ky*(yv(j)-0.5_DBL))/derf(0.5_DBL*ky))
end do
y_eta(1)=(3._DBL*yv_eta(1)-yv_eta(2))/2._DBL
! y_eta(1)=yv_eta(1)
do j=1,nj-1
  y_eta(j+1)=(yv_eta(j+1)+yv_eta(j))/2._DBL
end do
y_eta(nj+1)=(3._DBL*yv_eta(nj)-yv_eta(nj-1))/2._DBL
!
!----------------------------------------------------
!
! Puntos en y para capa l\'imite
!
! h = 3.0_DBL
! beta = 1.01_DBL
! do j=1,nj
!    yv_eta(j) = h * ( beta+1.0_DBL - ( beta-1.0_DBL ) * ( (beta+1.0_DBL)/(beta-1.0_DBL) )**( 1.0_DBL-yv(j) ) ) &
!         &/ ( ( (beta+1.0_DBL)/(beta-1.0_DBL) )**( 1.0_DBL-yv(j) ) + 1.0_DBL )
! end do
! y_eta(1)=(3._DBL*yv_eta(1)-yv_eta(2))/2._DBL
! ! y_eta(1)=yv_eta(1)
! do j=1,nj-1
!   y_eta(j+1)=(yv_eta(j+1)+yv_eta(j))/2._DBL
! end do
! y_eta(nj+1)=(3._DBL*yv_eta(nj)-yv_eta(nj-1))/2._DBL
! y_eta(nj+1)=yv_eta(nj)
!************************
!velocidad inicial en x
do j=1,nj+1
  u(1,j)=0._DBL
  do i=2,mi
    u(i,j)=0._DBL
  end do
end do
u(1,1)=0._DBL
u(1,nj+1)=0._DBL
!**********************************
!**********************************
!archivo de escritura de la malla u
Open(unit=1,file='out_n'//njc//'m'//mic//'_Rxxxu.dat')
write(1,*) i_o,i_1,0,ao
do j=1,nj+1
  do i=1,mi
    write(1,24) xu_eta(i),y_eta(j),u(i,j);24 format (3D23.15)
  end do
end do
close(unit=1)
!**********************************
!archivo de escritura de la malla v
Open(unit=2,file='out_n'//njc//'m'//mic//'_Rxxxv.dat')
write(2,*) i_o,i_1,0,ao
do j=1,nj
  do i=1,mi+1
    write(2,24) x_eta(i),yv_eta(j),0._DBL
  end do
end do
close(unit=2)
!**********************************
!archivo de escritura de la malla p
Open(unit=3,file='out_n'//njc//'m'//mic//'_Rxxxp.dat')
write(3,*) i_o,i_1,0,ao
do j=1,nj+1
  do i=1,mi+1
    write(3,25) x_eta(i),y_eta(j),0._DBL,0._DBL;25 format (4D23.15)
  end do
end do
close(unit=3)
!**********************************
!*** Formato de escritura Tecplot ***
Open(unit=1,file='malla_irregularp.plt')
write(1,*)'TITLE     = "malla" '
write(1,*)'VARIABLES = "x"'
write(1,*)'"y"'
write(1,*)'ZONE T= "Matrix"'
write(1,*)'I=',mi,' J=',nj+1,' K=1,F=POINT'
23 format (2D23.15)
do j=1,nj+1
  do i=1,mi
    write(1,23) y_eta(j),ao-xu_eta(i)
  end do
end do
close(unit=1)
End program malla_irr
