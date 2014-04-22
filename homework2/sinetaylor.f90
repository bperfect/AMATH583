!==========================
function sinetaylor(x,n)
!==========================
    implicit none

    ! function arguments:
    real (kind=8), intent(in) :: x
    integer, intent(in) :: n
    real (kind=8) :: sinetaylor

    ! local variables:
    real (kind=8) :: term, fixedTerm, partial_sum
    integer :: j

    partial_sum = 0.
    term = 1.

    do j=1,n
        ! j'th term is  x**j / j!  which is the previous term times x/j:
        term = term*x/j
        IF (MOD(j,2)==0) THEN
            fixedTerm=0
        ELSE IF (MOD(j,4)==3) THEN
            fixedTerm=term*(-1)
        ELSE
            fixedTerm=term
        END IF           

        ! add this term to the partial sum:
        partial_sum = partial_sum + fixedTerm   
        enddo
     sinetaylor = partial_sum  ! this is the value returned
end function sinetaylor

