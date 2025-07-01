Program malla_irr_3D
use constantes
Implicit none
Integer i,j,k,i_o,i_1
!*****************************************
!Variables de la malla, volumen de control
REAL(kind=DBL) xu(mi),yv(nj),zw(lk),ao,aox,aoy,aoz,as,u(mi,nj+1)
REAL(kind=DBL) xu_eta(mi),yv_eta(nj),zw_eta(lk),x_eta(mi+1),y_eta(nj+1),z_eta(lk+1)
REAL(kind=DBL) my,ky,mx,kx,mz,kz,tau,b,zc
REAL(kind=DBL) z_placa_min, z_placa_max
character(len=3) :: njc,mic,lkc
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
if (lk<100) then
  write(lkc,170) int(lk)
  lkc='0'//lkc
else
  write(lkc,160) int(lk)
endif
!***********************
!par'ametros de la malla
aox=1.00_DBL
aoy=1.00_DBL
aoz=1.00e0_DBL
as=0.0_DBL
zc=0.01_DBL
z_placa_min=4.5_DBL
z_placa_max=5.5_DBL
!********************************
!puntos en x
!generaci'on de la malla regular
do i = 1, nsolid-1
   xu_eta(i) = as * dfloat(i-1)/dfloat(nsolid-1)
end do
!*********************************
!generaci'on de la malla irregular
do i=nsolid,mi-nsolid+1
  xu(i)=dfloat(i-(nsolid))/dfloat(mi-2*nsolid+1)
end do
kx=5._DBL
mx=0.95_DBL !empaquetamiento en las paredes,menor coeficiente,mayor empaquetamiento
do i=nsolid,mi-nsolid+1
  xu_eta(i) = as + aox * (0.5_DBL+mx*(xu(i)-0.5_DBL)+0.5_DBL*(1-mx)*derf(kx*(xu(i)-0.5_DBL))/derf(0.5_DBL*kx))
end do
!********************************
!gerneraci'on de la malla regular
do i = mi-nsolid+2,mi
   xu_eta(i) = aox + as + as * dfloat(i+1-(mi-nsolid+2))/dfloat(nsolid-1)
end do
!**********************
!malla para la presi'on
x_eta(1)=(3._DBL*xu_eta(1)-xu_eta(2))/2._DBL
do i=1,mi-1
  x_eta(i+1)=(xu_eta(i+1)+xu_eta(i))/2._DBL
end do
x_eta(mi+1)=(3._DBL*xu_eta(mi)-xu_eta(mi-1))/2._DBL !xu_eta(mi)
!*******************************
!*******************************
!puntos en y
!generaci'on de la malla regular
do j = 1, nsolid-1
   yv_eta(j) = as * dfloat(j-1)/dfloat(nsolid-1)
end do
!*********************************
!generaci'on de la malla irregular
do j=nsolid,nj-nsolid+1
  yv(j)=dfloat(j-(nsolid))/dfloat(nj-2*nsolid+1)
end do
ky=5._DBL
my=0.95_DBL !empaquetamiento en las paredes,menor coeficiente,mayor empaquetamiento
do j=nsolid,nj-nsolid+1
  yv_eta(j) = as + aoy * (0.5_DBL+my*(yv(j)-0.5_DBL)+0.5_DBL*(1-my)*derf(ky*(yv(j)-0.5_DBL))/derf(0.5_DBL*ky))
end do
!********************************
!gerneraci'on de la malla regular
do j = nj-nsolid+2,nj
   yv_eta(j) = aoy + as + as * dfloat(j+1-(nj-nsolid+2))/dfloat(nsolid-1)
end do
!**********************
!malla para la presi'on
y_eta(1)=(3._DBL*yv_eta(1)-yv_eta(2))/2._DBL
do j=1,nj-1
  y_eta(j+1)=(yv_eta(j+1)+yv_eta(j))/2._DBL
end do
y_eta(nj+1)=(3._DBL*yv_eta(nj)-yv_eta(nj-1))/2._DBL
!*******************************
!*******************************
!puntos en z
!generaci'on de la malla regular
do k = 1, nsolid-1
   zw_eta(k) = as * dfloat(k-1)/dfloat(nsolid-1)
end do
!*********************************
!generaci'on de la malla irregular
do k=nsolid,lk-nsolid+1
  zw(k)=dfloat(k-(nsolid))/dfloat(lk-2*nsolid+1)
end do
kz=5._DBL
mz=0.15_DBL !empaquetamiento en las paredes,menor coeficiente,mayor empaquetamiento
do k=nsolid,lk-nsolid+1
  zw_eta(k) = as + aoz * (0.5_DBL+mz*(zw(k)-0.5_DBL)+0.5_DBL*(1-mz)*derf(kz*(zw(k)-0.5_DBL))/derf(0.5_DBL*kz))
end do
!********************************
!gerneraci'on de la malla regular
do k = lk-nsolid+2,lk
   zw_eta(k) = aoz + as + as * dfloat(k+1-(lk-nsolid+2))/dfloat(nsolid-1)
end do
! !******************************
! !malla empaquetada en el centro
! !puntos en z
! do k = 1, nsolid-1
!    zw_eta(k) = as * dfloat(k-1)/dfloat(nsolid-1)
! end do
! !*********************************
! !generaci'on de la malla irregular
! do k=nsolid,lk-nsolid+1
!   zw(k)=dfloat(k-(nsolid))/dfloat(lk-2*nsolid+1)
! end do
! tau=5._DBL
! b=1._DBL/(2._DBL*tau)*dlog((1._DBL+(dexp(tau)-1._DBL)*(zc/aoz))/(1._DBL+(dexp(-tau)-1._DBL)*(zc/aoz)))
! do k =nsolid,lk-nsolid+1
!   zw_eta(k) = as + zc*(1._DBL+dsinh(tau*(zw(k)-b))/dsinh(tau*b))
! end do
! !********************************
! !gerneraci'on de la malla regular
! do k = lk-nsolid+2,lk
!    zw_eta(k) = aoz + as + as * dfloat(k+1-(lk-nsolid+2))/dfloat(nsolid-1)
! end do
!**********************
!malla para la presi'on
z_eta(1)=(3._DBL*zw_eta(1)-zw_eta(2))/2._DBL
do k=1,lk-1
  z_eta(k+1)=(zw_eta(k+1)+zw_eta(k))/2._DBL
  if(z_eta(k)<z_placa_min)then
    i_o=k
  end if
  if(z_eta(k)<=z_placa_max)then
    i_1=k
  end if 
end do
z_eta(lk+1)=(3._DBL*zw_eta(lk)-zw_eta(lk-1))/2._DBL !zw_eta(lk)
if(z_eta(lk+1)<=z_placa_max)then
  i_1=lk+1
end if
i_o=i_o+1
!************************
!velocidad inicial en x
! do j=1,nj+1
!   u(1,j)=0._DBL
!   do i=2,mi
!     u(i,j)=0._DBL
!   end do
! end do
! u(1,1)=0._DBL
! u(1,nj+1)=0._DBL
!**********************************
!**********************************
!archivo de escritura de la malla u
Open(unit=1,file='out_n'//njc//'m'//mic//'k'//lkc//'_Rxxxu.dat')
write(1,*) i_o,i_1,0,aoz
do k=1,lk+1
  do j=1,nj+1
    do i=1,mi
      write(1,24) xu_eta(i),y_eta(j),z_eta(k), 0._DBL;24 format (4D23.15)
    end do
  end do
end do
close(unit=1)
!**********************************
!archivo de escritura de la malla v
Open(unit=2,file='out_n'//njc//'m'//mic//'k'//lkc//'_Rxxxv.dat')
write(2,*) i_o,i_1,0,aoz
do k=1,lk+1
  do j=1,nj
    do i=1,mi+1
      write(2,24) x_eta(i),yv_eta(j),z_eta(k),0._DBL
    end do
  end do
end do
close(unit=2)
!**********************************
!archivo de escritura de la malla w
Open(unit=3,file='out_n'//njc//'m'//mic//'k'//lkc//'_Rxxxw.dat')
write(3,*) i_o,i_1,0,aoz
do k=1,lk
  do j=1,nj+1
    do i=1,mi+1
      write(3,24) x_eta(i),y_eta(j),zw_eta(k), 0._DBL
    end do
  end do
end do
close(unit=3)
!**********************************
!archivo de escritura de la malla p
Open(unit=4,file='out_n'//njc//'m'//mic//'k'//lkc//'_Rxxxp.dat')
write(4,*) i_o,i_1,0,aoz
do k=1,lk+1
  do j=1,nj+1
    do i=1,mi+1
      write(4,25) x_eta(i),y_eta(j),z_eta(k),0._DBL,0._DBL;25 format (5D23.15)
    end do
  end do
end do
close(unit=4)
!**********************************
!*** Formato de escritura Tecplot ***
Open(unit=1,file='malla_irregularp.plt')
write(1,*)'TITLE     = "malla" '
write(1,*)'VARIABLES = "x"'
write(1,*)'"y"'
write(1,*)'"z"'
write(1,*)'ZONE T= "Matrix"'
write(1,*)'I=',mi+1,' J=',nj+1,' K=',lk+1,' F=POINT'
23 format (3D23.15)
do k=1,lk+1
  do j=1,nj+1
    do i=1,mi+1
      write(1,23) x_eta(i),y_eta(j),z_eta(k)
    end do
  end do
end do
close(unit=1)
End program malla_irr_3D
