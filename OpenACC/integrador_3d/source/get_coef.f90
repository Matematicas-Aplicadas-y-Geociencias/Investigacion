subroutine get_coef(dx1,dx2,dy1,dy2,df0,df1,df2,df3,df4,a,b,c,d,e)

    use mod_constantes

    implicit none

    real(kind=DBL), intent(in):: dx1, dx2, dy1, dy2
    real(kind=DBL), intent(in):: df0, df1, df2, df3, df4
    real(kind=DBL), intent(out):: a, b, c, d, e

    real(kind=DBL) :: numerador, denominador

   ! Coeficiente a
   numerador = dx2 * df1 + dx1 * df2
   denominador = dx1 * dx2 * (dx1 + dx2)
   a = numerador / denominador

   ! Coeficiente c
   numerador = dx2**2 * df1 - dx1**2 * df2
   ! denominador = dx1 - dx2 + epsilon(1.0d0)
   denominador = dx1 * dx2 * (dx1 + dx2)
   c = numerador / denominador

   ! Coeficiente b
   numerador = dy2 * df3 + dy1 * df4
   denominador = dy1 * dy2 * (dy1 + dy2)
   b = numerador / denominador

   ! Coeficiente d
   numerador = dy2**2 * df3 - dy1**2 * df4
   ! denominador = dy1 - dy2 + epsilon(1.0d0)
   denominador = dy1 * dy2 * (dy1 + dy2)
   d = numerador / denominador

   e = df0

   print '(A, F12.6)', 'Coef a:', a
   print '(A, F12.6)', 'Coef b:', b
   print '(A, F12.6)', 'Coef c:', c
   print '(A, F12.6)', 'Coef d:', d
   print '(A, F12.6)', 'Coef e:', e

end subroutine get_coef
