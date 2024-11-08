subroutine get_coef(dx1,dx2,dy1,dy2,df1,df2,df3,df4,df5,a,b,c,d,e)

    use mod_constantes

    implicit none

    real(kind=DBL), intent(in):: dx1, dx2, dy1, dy2
    real(kind=DBL), intent(in):: df1, df2, df3, df4, df5
    real(kind=DBL), intent(out):: a, b, c, d, e

    real(kind=DBL) numerador, denominador

   ! Coeficiente a
   numerador = dx2 * df1 + dx1 * df2
   denominador = dx1 * dx2 * (dx1 + dx2)
   a = numerador / denominador

   ! Coeficiente c
   numerador = df1 + df2 - a * (dx1**2 + dx2**2)
   denominador = dx1 - dx2
   c = numerador / denominador

   ! Coeficiente b
   numerador = dy2 * df4 + dy1 * df5
   denominador = dy1 * dy2 * (dy1 + dy2)
   b = numerador / denominador

   ! Coeficiente d
   numerador = df4 + df5 - b * (dy1**2 + dy2**2)
   denominador = dy1 - dy2
   d = numerador / denominador

   e = df3

   print *,a
   print *,c
   print *,'-------------------------'
   print *,b
   print *,d

end subroutine get_coef
