Program malla_irr
  use constantes
  Implicit none
  Integer i,j,k,i_o,i_1
  !*****************************************
  !Variables de la malla, volumen de control
  REAL(kind=DBL)    :: xu(mi),yv(nj)
  real(kind=DBL)    :: u(mi,nj+1)
  real(kind=DBL)    :: uf(mi+1,nj+1), vf(mi+1,nj+1), temp(mi+1,nj+1), pres(mi+1,nj+1)
  real(kind=DBL)    :: ao,xc,bo,yc
  REAL(kind=DBL)    :: xu_eta(mi),yv_eta(nj),x_eta(mi+1),y_eta(nj+1)
  REAL(kind=DBL)    :: my,ky,mx,kx,tau,b,h,beta
  REAL(kind=DBL)    :: x_placa_min, x_placa_max
  character(len=20) :: njc,mic
  character(len=4)  :: tipox, tipoy
  character(len=46) :: archivo=repeat(' ',46)
  !
  !**************************************
  !
  ! Lectura de par\'ametros para la malla
  !
  open(unit=101,file='entrada-malla.dat')
  read(101,*) tipox
  read(101,*) ao
  read(101,*) xc
  read(101,*) tipoy
  read(101,*) bo
  read(101,*) yc
  close(unit=101)
  !
  !--------------------------------------------------------------
  !
  ! Definici'on de caracteres para nombres de archivos y mensajes
  !
  mic = entero_caracter(mi)
  njc = entero_caracter(nj)
  !
  print*,"DEBUG: ", tipox, tipoy, ao, xc, bo, yc, mic, njc
  !
  uf   = 0._DBL
  vf   = 0._DBL
  temp = 0._DBL
  pres = 0._DBL
  !*********************************
  ! generaci'on de la malla irregular
  ! puntos en x
  !
  do i = 1, mi
     xu(i)=dfloat(i-1)/dfloat(mi-1)
  end do
  !
  if( tipox == 'cent' )then
     !
     ! Malla para x compactada alrededor de xc
     !
     x_placa_min=8.0_DBL
     x_placa_max=9.0_DBL
     tau=2._DBL
     b=1._DBL/(2._DBL*tau)*dlog((1._DBL+(dexp(tau)-1._DBL)*(xc/ao))/(1._DBL+(dexp(-tau)-1._DBL)*(xc/ao)))
     do i =1,mi
        xu_eta(i)=xc*(1._DBL+dsinh(tau*(xu(i)-b))/dsinh(tau*b))
     end do
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
  else if( tipox == 'tubo' )then
     !
     ! Malla para x compactada en las paredes
     !
     ky=10._DBL
     my=0.7_DBL !empaquetamiento en las paredes,menor coeficiente,mayor empaquetamiento
     do j=1,mi
        xu_eta(j)=1.0_DBL*(0.5_DBL+my*(xu(j)-0.5_DBL)+0.5_DBL*(1-my)*derf(ky*(xu(j)-0.5_DBL))/derf(0.5_DBL*ky))
     end do
     x_eta(1)=(3._DBL*xu_eta(1)-xu_eta(2))/2._DBL
     ! y_eta(1)=yv_eta(1)
     do j=1,mi-1
        x_eta(j+1)=(xu_eta(j+1)+xu_eta(j))/2._DBL
     end do
     x_eta(mi+1)=(3._DBL*xu_eta(mi)-xu_eta(mi-1))/2._DBL
     !  
  end if
  !***********************
  !***********************
  !puntos en y
  do j=1,nj
     yv(j)=dfloat(j-1)/dfloat(nj-1)
  end do
  !
  if( tipoy == 'cent' )then
     !
     ! Malla para y compactada alrededor de yc
     !
     tau=2._DBL
     b=1._DBL/(2._DBL*tau)*dlog((1._DBL+(dexp(tau)-1._DBL)*(yc/bo))/(1._DBL+(dexp(-tau)-1._DBL)*(yc/bo)))
     do i =1,nj
        yv_eta(i)=yc*(1._DBL+dsinh(tau*(yv(i)-b))/dsinh(tau*b))
     end do
     !************************
     y_eta(1)=(3._DBL*yv_eta(1)-yv_eta(2))/2._DBL
     ! x_eta(1)=xu_eta(1)
     do i=1,nj-1
        y_eta(i+1)=(yv_eta(i+1)+yv_eta(i))/2._DBL
     end do
     y_eta(nj+1)=(3._DBL*yv_eta(nj)-yv_eta(nj-1))/2._DBL
     ! x_eta(mi+1)=xu_eta(mi)
  else if( tipoy == 'tubo' )then
     !
     ! Malla para y compactada en las paredes
     !
     ky=10._DBL
     my=0.7_DBL !empaquetamiento en las paredes,menor coeficiente,mayor empaquetamiento
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
  end if
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
  Open(unit=1,file='out_n'//trim(njc)//'m'//trim(mic)//'_Rxxxu.dat')
  write(1,*) i_o,i_1,0,ao
  do j=1,nj+1
     do i=1,mi
        write(1,24) xu_eta(i),y_eta(j),u(i,j);24 format (3D23.15)
     end do
  end do
  close(unit=1)
  !**********************************
  !archivo de escritura de la malla v
  Open(unit=2,file='out_n'//trim(njc)//'m'//trim(mic)//'_Rxxxv.dat')
  write(2,*) i_o,i_1,0,ao
  do j=1,nj
     do i=1,mi+1
        write(2,24) x_eta(i),yv_eta(j),0._DBL
     end do
  end do
  close(unit=2)
  !**********************************
  !archivo de escritura de la malla p
  Open(unit=3,file='out_n'//trim(njc)//'m'//trim(mic)//'_Rxxxp.dat')
  write(3,*) i_o,i_1,0,ao
  do j=1,nj+1
     do i=1,mi+1
        write(3,25) x_eta(i),y_eta(j),0._DBL,0._DBL;25 format (4D23.15)
     end do
  end do
  close(unit=3)
  !
  archivo = 'visualizar.vtk'
  call postproceso_vtk(x_eta,y_eta,uf,vf,pres,temp,temp,archivo)
  !
end program malla_irr
