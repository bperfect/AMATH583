
module approx_norm

    implicit none
    integer :: nsamples, seed, method
    save

contains

    subroutine approx_norm1(a, anorm)

    use random_util, only: init_random_seed
    use omp_lib

    implicit none
    real(kind=8), dimension(:,:), intent(in) :: a
    real(kind=8), intent(out) :: anorm
    integer :: n, i, j
    real(kind=8), allocatable, dimension(:) :: x,c
    real(kind=8) :: csum, xsum, norm
    
    anorm=0 !Set the norm to the lowest possible value
    n=size(a,1)
    call init_random_seed(seed) !Set up our ability to create random vectors 
    allocate(x(n))
    allocate(c(n))

    !#omp parallel do if (method==1) private(x, csum, xsum, c, norm)
    do j=1,nsamples
        call random_number(x)
        !#omp parallel if (method==2)
        c=matmul(a,x)
        csum=0
        xsum=0
        do i=1,n
            csum=csum+abs(c(i))
            xsum=xsum+abs(x(i))
        enddo
        norm=csum/xsum
        anorm=max(norm,anorm)
    enddo

    end subroutine approx_norm1

end module approx_norm
