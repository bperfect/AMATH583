
! $sinetaylor.f90

program taylor_main2

    implicit none                  
    real (kind=8) :: x, sine_true, y
    real (kind=8), external :: sinetaylor
    integer :: n

    n = 20               ! number of terms to use
    x = 1.0
    sine_true = sin(x)
    y = sinetaylor(x,n)   ! uses function below
    print *, "x = ",x
    print *, "n = ",n
    print *, "exp_true  = ",sine_true
    print *, "sinetaylor = ",y
    print *, "error     = ",y - sine_true

end program taylor_main2


