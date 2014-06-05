
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

!Solve the problem by splitting the domain in half
subroutine solve_bvp_split(x, u_left, u_right, u)
    use problem, only: f
    implicit none
    real(kind=8), intent(in) :: x(0:)
    real(kind=8), intent(in) :: u_left, u_right
    real(kind=8), intent(out) :: u(0:)
    integer:: INFO, pivot
    real(kind=8) :: dx, G0, G1, u_mid, z
    real(kind=8), dimension(:), allocatable :: V0, V1, u1, u2, x1, x2

    ! declare local variables and continue writing this code...
    dx = x(2)-x(1)
    
    !want to include the center point in each subproblem
    pivot = FLOOR(SIZE(x)/2.0)
    allocate(V0(SIZE(x)))
    allocate(V1(SIZE(x)))
    allocate(x1(FLOOR(SIZE(x)/2.0+1)))
    allocate(x2(pivot))
    allocate(u1(FLOOR(SIZE(x)/2.0+1)))
    allocate(u2(pivot))
    
    x1 = x(0:pivot)
    x2 = x(pivot+1:size(x))
    
    u1 = u(0:pivot)
    u2 = u(pivot+1:size(u))
    u_mid = 0
    call solve_bvp_direct(x1,u_left,u_mid,u1)
    call solve_bvp_direct(x2,u_mid,u_right,u2)
    G0 = (u1(SIZE(u1)-1)-2*u_mid+u2(2))+dx**2*f(x(pivot)) !solve the PDE at the middle
    V0 = (/ u1, u2(2:size(u2)) /) !Concatenate matrices together
    
    u1 = u(0:pivot)
    u2 = u(pivot+1:size(u))
    u_mid = 1
    call solve_bvp_direct(x1,u_left,u_mid,u1)
    call solve_bvp_direct(x2,u_mid,u_right,u2)
    G1 = (u1(SIZE(u1)-1)-2*u_mid+u2(2))+dx**2*f(x(pivot)) !solve the PDE at the middle
    V1 = (/ u1, u2(2:size(u2)) /) !Concatenate matrices together
    
    z=G1/(G1-G0)
    
    print *, 'Computed G0= ', G0, ' Computed G1= ', G1, 'z= ', z
    u = z*V0+(1-z)*V1
    
    deallocate(V0)
    deallocate(V1)
    deallocate(x1)
    deallocate(x2)
    deallocate(u1)
    deallocate(u2)

end subroutine solve_bvp_split

!perform the split domain simulation using openMP
subroutine solve_bvp_split_omp(x, u_left, u_right, u)
    use problem, only: f
    use omp_lib
    implicit none
    real(kind=8), intent(in) :: x(0:)
    real(kind=8), intent(in) :: u_left, u_right
    real(kind=8), intent(out) :: u(0:)
    integer:: pivot,thread
    real(kind=8) :: dx, G0, G1, u_mid1, u_mid2, z
    real(kind=8), dimension(:), allocatable :: V0, V1, u1, u2, u3, u4, x1, x2

    ! declare local variables and continue writing this code...
    dx = x(2)-x(1)
    !want to include the center point in each subproblem
    pivot = FLOOR(SIZE(x)/2.0)
    allocate(V0(SIZE(x)))
    allocate(V1(SIZE(x)))
    allocate(x1(FLOOR(SIZE(x)/2.0)+1))
    allocate(x2(pivot))
    allocate(u1(FLOOR(SIZE(x)/2.0+1)))
    allocate(u3(pivot))
    allocate(u2(FLOOR(SIZE(x)/2.0+1)))
    allocate(u4(pivot))
    
    !$ call omp_set_num_threads(4)
    u_mid1=0
    u_mid2=1
    
    !$omp parallel shared(u_left,u_right,pivot,x,u1,u2,u3,u4,u_mid1,u_mid2), private(x1,x2,thread)
    !$omp sections
    !$omp section
        x1 = x(0:pivot)
        thread=omp_get_thread_num()
        print 155, thread, x1(1),x1(pivot+1),u_mid1
155     format ('Thread ', i6, ' taking from ', f5.3, ' to ', f5.3, ' with u_mid = ', f3.1)
        call solve_bvp_direct(x1,u_left,u_mid1,u1)
    !$omp section
        x1 = x(0:pivot)
        thread=omp_get_thread_num()
        print 156, thread, x1(1),x1(pivot+1),u_mid2
156     format ('Thread ', i6, ' taking from ', f5.3, ' to ', f5.3, ' with u_mid = ', f3.1)
        call solve_bvp_direct(x1,u_left,u_mid2,u2)
    !$omp section
        x2 = x(pivot+1:size(x))
        thread=omp_get_thread_num()
        print 157, thread, x2(1),x2(pivot+1),u_mid1
157     format ('Thread ', i6, ' taking from ', f5.3, ' to ', f5.3, ' with u_mid = ', f3.1)
        call solve_bvp_direct(x2,u_mid1,u_right,u3)
    !$omp section
        x2 = x(pivot+1:size(x))
        thread=omp_get_thread_num()
        print 158, thread, x2(1),x2(pivot+1),u_mid2
158     format ('Thread ', i6, ' taking from ', f5.3, ' to ', f5.3, ' with u_mid = ', f3.1)
        call solve_bvp_direct(x2,u_mid2,u_right,u4)
    !$omp end sections
    !$omp end parallel

    
    G0 = (u1(SIZE(u1)-1)-2*u_mid1+u3(2))+dx**2*f(x(pivot)) !solve the PDE at the middle
    V0 = (/ u1, u3(2:size(u3)) /) !Concatenate matrices together
    G1 = (u2(SIZE(u2)-1)-2*u_mid2+u4(2))+dx**2*f(x(pivot)) !solve the PDE at the middle
    V1 = (/ u2, u4(2:size(u4)) /) !Concatenate matrices together

    z=G1/(G1-G0)
    print *, 'Computed G0= ', G0, ' Computed G1= ', G1, 'z= ', z
    u = z*V0+(1-z)*V1
    
    deallocate(V0)
    deallocate(V1)
    deallocate(x1)
    deallocate(x2)
    deallocate(u1)
    deallocate(u3)
    deallocate(u2)
    deallocate(u4)

end subroutine solve_bvp_split_omp

end module bvp_solvers
