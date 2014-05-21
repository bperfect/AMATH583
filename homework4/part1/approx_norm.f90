
module approx_norm

    implicit none
    integer :: nsamples, seed
    save

contains

    subroutine approx_norm1(a, anorm)

    use omp_lib
    implicit none
    real(kind=8), dimension(:,:), intent(in) :: a
    real(kind=8), intent(out) :: anorm

    
    real(kind=8), allocatable, dimension(:) :: x,ax,allx
    integer :: i,j,k,n,k1,k2,nthreads,istart,iend,thread_num,points_per_thread
    real(kind=8) :: xnorm, axnorm, ratio, max_ratio

    n = size(a,1) 
    allocate(x(n),ax(n))

    ! generate all the random numbers we'll ever need:
    allocate(allx(n*nsamples))
    call random_number(allx)

    max_ratio = 0.d0

! Start the parallel section
    !$omp parallel private(k,i,j,x,ax,xnorm,axnorm,ratio,istart,iend,thread_num)
   
!Figure out what each thread's task is 
    !$ nthreads = omp_get_num_threads()
    points_per_thread = (nsamples + nthreads - 1) / nthreads
    thread_num = 0     ! needed in serial mode
    !$ thread_num = omp_get_thread_num()    ! unique for each thread


    ! Determine start and end index for the set of points to be 
    ! handled by this thread:
    istart = thread_num * points_per_thread + 1
    iend = min((thread_num+1) * points_per_thread, nsamples)
    print *, "Thread ",thread_num, " will take k= ", istart, " through k= ", iend
!Perform the main work for each thread
    do k=istart,iend
        ! choose a random vector x by taking n elements from allx :
        k1 = (k-1)*n+1
        k2 = k*n
        x = allx(k1:k2)

        ! compute matrix-vector product ax:
        do i=1,n    
            ax(i) = 0.d0
            do j=1,n
                ax(i) = ax(i) + a(i,j)*x(j)
            enddo
        enddo
            ! compute 1-norm of x and ax
        xnorm = 0.d0
        axnorm = 0.d0
        do i=1,n
            xnorm = xnorm + abs(x(i))
            axnorm = axnorm + abs(ax(i))
        enddo

        ratio = axnorm/xnorm
        max_ratio = max(max_ratio, ratio)
    enddo
!Combine results to get a final answer
    !$omp critical
    anorm = max(anorm,max_ratio)
    !$omp end critical
    !$omp end parallel
    end subroutine approx_norm1

end module approx_norm
