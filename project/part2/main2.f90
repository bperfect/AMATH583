
program main_heat2

    use heat_solvers, only: solve_heat_explicit,solve_heat_implicit
    use problem, only: k, u_true
    implicit none
    integer :: n, i, nsteps, method, j
    real(kind=8) :: error_max, dx, dt, pi, tfinal, t0, zero
    real(kind=8), dimension(:), allocatable :: x, u, ustar, u0

    pi = acos(-1.d0)
    zero=0

    open(unit=21, file='input_data.txt', status='old')
    read(21,*) n
    read(21,*) k
    read(21,*) tfinal
    read(21,*) nsteps
    read(21,*) method
    print '("n = ",i6)', n
    print '("k = ",i6)', k
    print '("tfinal = ",f7.4)', tfinal
    print '("nsteps = ",i6)', nsteps
    print '("method = ",i2)', method

    allocate(x(0:n+1), u(0:n+1), ustar(0:n+1), u0(0:n+1))

    t0 = 0.d0
    dt = (tfinal-t0)/float(nsteps)
    dx = pi / (n+1)
    do i=0,n+1
        x(i) = i*dx
        enddo


    do i=0,n+1
        u0(i) = sin(k*x(i))
        ustar(i) = u_true(x(i), tfinal)
    enddo
        
    !write initial solution to file too 
    open(unit=22, file='solution.txt', status='unknown')
    do i=0,n+1
        write(22,220) u0(i)
220     format(3e22.14)
    enddo

    do i=1,nsteps
        if (method == 1) then
            call solve_heat_explicit(x, u0, zero, dt, 1, u)
        else
            call solve_heat_implicit(x, u0, zero, dt, 1, u)
        endif
        !write the data
        do j=0,n+1
            if (abs(u(j))<1e-15) then !Prevent numerical errors from machine epsilon going from fortran to python
                u(j) = 0
            endif
            write(22,220) u(j)
            !print *, u(j)
        enddo
        u0 = u !Set the old u for the next iteration
    enddo
    close(22)


    error_max = 0.d0
    do i=0,n+1
        ustar(i) = u_true(x(i), tfinal)
        error_max = max(error_max, abs(u(i) - ustar(i)))
        enddo

    print '("error_max = ", e13.6)', error_max


end program main_heat2
