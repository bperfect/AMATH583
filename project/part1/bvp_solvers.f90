
module bvp_solvers

contains

subroutine solve_bvp_direct(x, u_left, u_right, u)
    use problem, only: f
    implicit none
    real(kind=8), intent(in) :: x(0:)
    real(kind=8), intent(in) :: u_left, u_right
    real(kind=8), intent(out) :: u(0:)
    integer:: INFO, N, i
    real(kind=8), dimension(:), allocatable :: D, DL, DU!, B
    real(kind=8) :: dx

    ! decleare local variables and continue writing this code...
    N = SIZE(u)
    dx = x(2)-x(1)
    allocate(D(N))
    allocate(DL(N-1))
    allocate(DU(N-1))

    print *, 'Solving tridiagonal system with n= ', N-2

    do i=0,N-2
        D(i+1) = -2
        DU(i+1) = 1
        DL(i+1) = 1
        U(i) = -dx**2*f(x(i))
    enddo

    !Set the boundary conditions
    D(1)=1
    D(N)=1
    U(0) = u_left
    U(N-1) = u_right
    DU(1) = 0
    DL(N-1) = 0

    call DGTSV(N, 1, DL, D, DU, U, N, INFO)
    
    deallocate(D)
    deallocate(DL)
    deallocate(DU)
    
end subroutine solve_bvp_direct

end module bvp_solvers
