module asserts
    use mod_constantes
    implicit none
contains

    subroutine assertFloat(valor_esperado, valor_obtenido)

        implicit none

        real(kind=DBL), intent(in) :: valor_esperado, valor_obtenido
        logical :: are_equals
        character(len=20) :: assert_result

        are_equals = nearly_equals(valor_esperado, valor_obtenido)
        assert_result = check(are_equals)
        print *, assert_result

    end subroutine assertFloat

    logical function nearly_equals(valor_esperado, valor_obtenido)

        implicit none

        real(kind=DBL), intent(in) :: valor_esperado, valor_obtenido
        real(kind=DBL) :: abs_value_a, abs_value_b, diff, eps

        eps = 1.0E-5 ! epsilon(valor_obtenido)

        abs_value_a = abs(valor_esperado)
        abs_value_b = abs(valor_obtenido)
        diff = abs(valor_esperado - valor_obtenido)

        if (valor_esperado == valor_obtenido) then
            nearly_equals = .true.
        else if (valor_esperado * valor_obtenido == 0) then
            nearly_equals = diff < (eps * eps)
        else
            nearly_equals = diff / (abs_value_a + abs_value_b) < eps
        end if
    end function nearly_equals

    character(len=20) function check(result)

        implicit none

        logical, intent(in) :: result

        if (result .eqv. .true.) then
            check = "Test passed"
        else
            check = "Test failure"
        end if

    end function check
end module asserts
