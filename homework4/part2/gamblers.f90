
module gamblers

    implicit none
    integer :: kwalks, max_steps
    save

contains

    subroutine walk(n1in, n2in, p, verbose, nsteps, winner)

    implicit none
    integer, intent(in) :: n1in,n2in
    real(kind=8), intent(in) :: p
    logical, intent(in) :: verbose
    integer, intent(out) :: nsteps, winner

    ! local variables
    real(kind=8) :: r
    integer :: nstep, n1, n2
    real(kind=8), dimension(max_steps) :: rand

    ! initialize n1 and n2 with input values
    ! need to change n1 and n2 internally during this walk but do not
    ! want to affect n1 and n2 in main program. 
    n1 = n1in
    n2 = n2in
    winner=0
    nsteps=0
    call random_number(rand)

    do nsteps=0,max_steps-1
        !Each iteration is one round of the game
        if (rand(nsteps)>p) then
            n1=n1+1
            n2=n2-1
        else
            n1=n1-1
            n2=n2+1
        endif

        if (verbose) then
        print *, 'In step ', nsteps, 'r= ',rand(nsteps), 'and n1= ',n1,'n2= ',n2
        endif

        !Check for win conditions
        if (n1==0) then 
            winner=2
            exit
        endif
        if (n2==0) then 
            winner=1
            exit
        endif
    enddo


    end subroutine walk

end module gamblers
