
module heat_solvers

contains

subroutine solve_heat_explicit(x, u0, t0, tfinal, nsteps, u)
    implicit none
    real(kind=8), intent(in) :: x(0:), u0(0:), t0, tfinal
    integer, intent(in) :: nsteps
    real(kind=8), intent(out) :: u(0:)
    integer :: n, i, nstep
    real(kind=8) :: dx, dt, dtdx2, um1, up1
    real(kind=8), allocatable :: uxx(:)

    n = size(x) - 2
    dx = x(2) - x(1)
    dt = (tfinal - t0) / float(nsteps)

    print 21, n, dx, dt
21  format("Solving heat equation with n = ",i6, " dx = ",e10.3," dt = ",e10.3)

    dtdx2 = dt / dx**2
    print 22, dtdx2
22  format("  dt/dx**2 = ",f7.3)
	if (dtdx2 > 0.5d0) then
		print *, "*** Warning: the explicit method is unstable"
		endif

    allocate(uxx(1:n))

    u = u0

    do nstep=1,nsteps

        do i=1,n
            uxx(i) = (u(i-1) - 2.d0*u(i) + u(i+1)) / dx**2
            enddo

        do i=1,n
            u(i) = u(i) + dt * uxx(i)
            enddo

        enddo

end subroutine solve_heat_explicit

subroutine solve_heat_implicit(x, u0, t0, tfinal, nsteps, u)
    implicit none
    real(kind=8), intent(in) :: x(:), u0(:), t0, tfinal
    integer, intent(in) :: nsteps
    real(kind=8), intent(out) :: u(:)
    integer :: n, i, nstep, INFO
    real(kind=8) :: dx, dt, dtdx2, um1, up1
    real(kind=8), allocatable :: uxx(:),DL(:),DU(:),D(:)
    
    n = size(x) - 2
    dx = x(2) - x(1)
    dt = (tfinal - t0) / float(nsteps)

    print 21, n, dx, dt
21  format("Solving heat equation using Crank-Nicholson with n = ",i6, " dx = ",e10.3," dt = ",e10.3)

    dtdx2 = dt / dx**2
    print 22, dtdx2
22  format("  dt/dx**2 = ",f7.3)

    allocate(uxx(1:n+2))
    allocate(DL(1:n+1))
    allocate(DU(1:n+1))
    allocate(D(1:n+2))
    u = u0
    
    do i=1,n+1
        DL(i) = -dtdx2/2
        DU(i) = -dtdx2/2
        D(i) = 1+dtdx2
    enddo
    D(n+2) = 1
    D(1) = 1
    !DU(1)=0
    !DL(n+1)=0
    !print *, D
    !print *, DU
    !print *, DL
    uxx=u0 !added this
    do nstep=1,nsteps
        print *, 'iter'
        uxx(1) = 0
        uxx(n+2) = 0
        u(1) = 0
        u(n+2)=0
        do i=2,n+1
            uxx(i) = u(i)+dtdx2*(-u(i)+u(i-1)/2+u(i+1)/2)
        enddo
        !print *, uxx
        call DGTSV(N, 1, DL, D, DU, uxx, N, INFO)
        !print *, uxx
        u=uxx
        u(1) = 0
        u(n+2)=0

    enddo

end subroutine solve_heat_implicit


end module heat_solvers
