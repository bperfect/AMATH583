
program main2

    use omp_lib
    use random_util, only: init_random_seed
    use gamblers, only: max_steps, walk

    implicit none
    real(kind=8) :: p, q, p1, p2, frac1wins, frac2wins,t1,t2
    integer :: thread_num,seed, i, nthreads, n1, n2, nsteps, winner,kwalks,draw
    integer(kind=8) :: tclock1,tclock2,clock_rate
    real(kind=8):: elapsed_time
    integer, allocatable, dimension(:) :: nsteps_thread

    open(unit=21, file='input_data.txt', status='old')
    read(21,*) n1
    read(21,*) n2
    read(21,*) p
    read(21,*) max_steps
    read(21,*) seed
    read(21,*) nthreads
    read(21,*) kwalks

    print "('n1 = ',i6)", n1
    print "('n2 = ',i6)", n2
    print "('p = ',f8.4)", p
    print "('max_steps = ',i9)", max_steps
    print "('nthreads = ',i6)", nthreads
    print "('kwalks = ',i6)", kwalks

    allocate(nsteps_thread(nthreads))
    do i=1,nthreads
        nsteps_thread(i)=0
    enddo
    call init_random_seed(seed)
    draw=0
    p1=0
    p2=0
    
    !$ call omp_set_num_threads(nthreads)
    !---- check the time for the parallel portion:

    call system_clock(tclock1)
    call cpu_time(t1)   ! start cpu timer

    !$omp parallel do private(thread_num,winner,nsteps),&
           !$omp  reduction(+:draw,p1,p2)

    do i=1,kwalks
        call walk(n1, n2, p, .false., nsteps, winner)
        !$ thread_num = omp_get_thread_num()
        nsteps_thread(thread_num+1)=nsteps_thread(thread_num+1)+nsteps
        if (winner==0) then
            draw=draw+1
        elseif (winner==1) then
            p1=p1+1
        else
            p2=p2+1
        endif
    enddo

! The work is all done, now we end the timer and print out the relevant info


!Timing information
    call cpu_time(t2)   ! end cpu timer
    print 10, t2-t1
 10 format("CPU time = ",f12.8, " seconds")
    call system_clock(tclock2, clock_rate)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print 11, elapsed_time
 11 format("Elapsed time = ",f12.8, " seconds")

!Draw Warning information
    if (draw > 0) then
        print *, 'Warning: ', draw, 'walks out of ',kwalks, 'did not result in a win by either player'
    endif

!Win probability information
    print *, 'Player 1 won ', p1/kwalks, 'fraction of the time'
    print *, 'Player 2 won ', p2/kwalks, 'fraction of the time'
    q=1-p
    if (p == 0.5) then
        frac1wins=0.5
        frac2wins=0.5
    else
        frac1wins=(1-(q/p)**n1)/(1-(q/p)**(n1+n2))
        frac2wins=1-frac1wins
    endif
    print *, 'True probabilities are P1 = ',frac1wins, 'P2 = ', frac2wins


!Path length information    
    print *, 'The average path length is ', sum(nsteps_thread)/kwalks
   if (p==0.5) then
        i=n1*n2
   else 
        i=n1/(q-p)-(n1+n2)/(q-p)*(1-(q/p)**n1)/(1-(q/p)**(n1+n2))
   endif
    print *, 'True mean path length is ',i 

!Thread use information 
    do i=1,nthreads
        print *, 'Thread ',i-1,' took ',nsteps_thread(i), ' steps'
    enddo
    deallocate(nsteps_thread)
end program main2
