module constantes
  integer, parameter :: mi = 1024, nj=512, nn = 1200000
  integer, parameter :: NUM_GANGSA = 512, VEC_LENGHTA = 1
  integer, parameter :: NUM_GANGSB = 1024, VEC_LENGHTB = 1
  real, parameter    :: pi = 3.1415926535
  INTEGER, PARAMETER :: DBL=SELECTED_REAL_KIND(P=15,R=300)
contains

    subroutine ensambla_tdmax_1D(aio,aco,ado,rxo,deltaxo,deltayo,t_anto,cond_tero,temp_inio,temp_fino,jjo)
      !$acc routine vector
      implicit none
      double precision, intent(in)                    :: deltaxo,deltayo,cond_tero
      double precision, intent(in)                    :: temp_inio,temp_fino
      double precision, dimension(mi*nj), intent(in)  :: t_anto
      
      double precision, dimension(mi*nj), intent(out) :: aio,ado
      double precision, dimension(mi*nj), intent(out) :: aco
      double precision, dimension(mi*nj), intent(out) :: rxo

      integer, intent(in) :: jjo
      integer             :: iin

      !
      ! Se definen las condiciones de frontera
      !
      aco((jjo-1)*mi+1)  = 1.0d0
      ado((jjo-1)*mi+1)  = 0.0d0
      rxo((jjo-1)*mi+1)  = temp_inio
      !
      aio((jjo-1)*mi+mi) = 0.0d0
      aco((jjo-1)*mi+mi) = 1.0d0
      rxo((jjo-1)*mi+mi) = temp_fino
      !
      ! Ensamblado de la matriz tridiagonal
      ! y del vector de resultados
      !
      !$acc loop vector
      do iin=2, mi-1
         aio((jjo-1)*mi+iin) =-1.0d0*cond_tero/(deltaxo*deltaxo)
         aco((jjo-1)*mi+iin) = 2.0d0*cond_tero*(1.d0/(deltaxo*deltaxo)+1.d0/(deltayo*deltayo))
         ado((jjo-1)*mi+iin) =-1.0d0*cond_tero/(deltaxo*deltaxo)
         rxo((jjo-1)*mi+iin) = cond_tero/(deltayo*deltayo)*t_anto((jjo)*mi+iin)+&
              &cond_tero/(deltayo*deltayo)*t_anto((jjo-2)*mi+iin)
      end do
      
    end subroutine ensambla_tdmax_1D
        !
    ! ---------------------------
    !
    subroutine ensambla_tdmay_1D(bio,bco,bdo,ryo,deltaxo,deltayo,tempero,cond_tero,flux_abao,flux_arro,iio)
      !$acc routine vector
      implicit none
      double precision, intent(in)                      :: deltaxo,deltayo,cond_tero
      double precision, intent(in)                      :: flux_arro,flux_abao
      double precision, dimension(mi*nj), intent(in)    :: tempero
      
      double precision, dimension(nj*mi), intent(out)   :: bio,bdo
      double precision, dimension(nj*mi), intent(out)   :: bco
      double precision, dimension(nj*mi), intent(out)   :: ryo
      integer :: iio,jjn
      !
      ! Se definen las condiciones de frontera
      !
      bco((iio-1)*nj+1)  =-1.0d0/deltayo ! 1.d0 
      bdo((iio-1)*nj+1)  = 1.0d0/deltayo ! 0.d0/deltayo 
      ryo((iio-1)*nj+1)  = flux_abao ! 308.d0 
      bio((iio-1)*nj+nj) =-1.d0/deltayo ! 0.d0   
      bco((iio-1)*nj+nj) = 1.0d0/deltayo ! 1.d0   
      ryo((iio-1)*nj+nj) = flux_arro ! 308.d0  !(temp_fin+temp_ini)/2.d0
      !
      ! Ensamblado de la matriz tridiagonal
      ! y del vector de resultados
      !
      !$acc loop vector
      do jjn=2, nj-1
         bio((iio-1)*nj+jjn) =-1.0d0*cond_tero/(deltayo*deltayo)
         bco((iio-1)*nj+jjn) = 2.0d0*cond_tero*(1.d0/(deltayo*deltayo)+1.d0/(deltaxo*deltaxo))
         bdo((iio-1)*nj+jjn) =-1.0d0*cond_tero/(deltayo*deltayo)
         ryo((iio-1)*nj+jjn) = cond_tero/(deltaxo*deltaxo)*tempero((jjn-1)*mi+iio+1)+&
              &cond_tero/(deltaxo*deltaxo)*tempero((jjn-1)*mi+iio-1)
      end do
    end subroutine ensambla_tdmay_1D
    !
    ! -------------------------------
    !
    subroutine ensambla_tdmax(AIo,ACo,ADo,resultxo,deltaxo,deltayo,temp_anto,cond_tero,temp_inio,temp_fino,jjo)
      !$acc routine vector
      implicit none
      double precision, intent(in)                    :: deltaxo,deltayo,cond_tero
      double precision, intent(in)                    :: temp_inio,temp_fino
      double precision, dimension(mi,nj), intent(in)  :: temp_anto
      
      double precision, dimension(mi,nj), intent(out) :: AIo,ADo
      double precision, dimension(mi,nj), intent(out) :: ACo
      double precision, dimension(mi,nj), intent(out) :: resultxo

      integer :: iin,jjo

      !
      ! Se definen las condiciones de frontera
      !
      ACo(1,jjo)       = 1.0d0
      ADo(1,jjo)       = 0.0d0
      resultxo(1,jjo)  = temp_inio
      AIo(mi,jjo)      = 0.d0
      ACo(mi,jjo)      = 1.0d0
      resultxo(mi,jjo) = temp_fino
      !
      ! Ensamblado de la matriz tridiagonal
      ! y del vector de resultados
      !
      !$acc loop vector
      do iin=2, mi-1
         AIo(iin,jjo)      =-1.0d0*cond_tero/(deltaxo*deltaxo)
         ACo(iin,jjo)      = 2.0d0*cond_tero*(1.d0/(deltaxo*deltaxo)+1.d0/(deltayo*deltayo))
         ADo(iin,jjo)      =-1.0d0*cond_tero/(deltaxo*deltaxo)
         resultxo(iin,jjo) = cond_tero/(deltayo*deltayo)*temp_anto(iin,jjo+1)+&
              &cond_tero/(deltayo*deltayo)*temp_anto(iin,jjo-1)
      end do
      
    end subroutine ensambla_tdmax
    !
    ! ---------------------------
    !
    subroutine ensambla_tdmay(BIo,BCo,BDo,resultyo,deltaxo,deltayo,tempero,cond_tero,flux_abao,flux_arro,iio)
      !$acc routine vector
      implicit none
      double precision, intent(in)                      :: deltaxo,deltayo,cond_tero
      double precision, intent(in)                      :: flux_arro,flux_abao
      double precision, dimension(mi,nj), intent(in)    :: tempero
      
      double precision, dimension(nj,mi), intent(out)   :: BIo,BDo
      double precision, dimension(nj,mi), intent(out)   :: BCo
      double precision, dimension(nj,mi), intent(out)   :: resultyo
      integer :: iio,jjn
      !
      ! Se definen las condiciones de frontera
      !
      BCo(1,iio)       =-1.0d0/deltayo
      BDo(1,iio)       = 1.0d0/deltayo
      resultyo(1,iio)  = flux_abao
      BIo(nj,iio)      =-1.d0/deltayo
      BCo(nj,iio)      = 1.0d0/deltayo
      resultyo(nj,iio) = flux_arro !(temp_fin+temp_ini)/2.d0
      !
      ! Ensamblado de la matriz tridiagonal
      ! y del vector de resultados
      !
      !$acc loop vector
      do jjn=2, nj-1
         BIo(jjn,iio)      =-1.0d0*cond_tero/(deltayo*deltayo)
         BCo(jjn,iio)      = 2.0d0*cond_tero*(1.d0/(deltayo*deltayo)+1.d0/(deltaxo*deltaxo))
         BDo(jjn,iio)      =-1.0d0*cond_tero/(deltayo*deltayo)
         resultyo(jjn,iio) = cond_tero/(deltaxo*deltaxo)*tempero(iio+1,jjn)+&
              &cond_tero/(deltaxo*deltaxo)*tempero(iio-1,jjn)
      end do
    end subroutine ensambla_tdmay
    !
    ! ---------------------------
    !    
    subroutine tri(a,b,c,r,n)
      !$acc routine
      implicit none
      integer, intent(in) :: n
      double precision, intent(in)    :: a(n),c(n)
      double precision, intent(inout) :: b(n),r(n)
      ! double precision, intent(inout) :: u(n)
      integer :: i
      
      ! eliminacion elementos bajo la matriz
      gausselim: do i=2,n
         r(i)=r(i)-(a(i)/b(i-1))*r(i-1)
         b(i)=b(i)-(a(i)/b(i-1))*c(i-1)
      end do gausselim
      ! solucion para r
      r(n)=r(n)/b(n)
      sustatras: do i=n-1,1,-1
         r(i)=(r(i)-c(i)*r(i+1))/b(i)
      end do sustatras

    end subroutine tri
    
end module
