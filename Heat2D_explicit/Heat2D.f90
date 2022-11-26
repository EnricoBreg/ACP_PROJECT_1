PROGRAM Heat2D
    USE OMP_LIB
    IMPLICIT NONE
    !-------------------------------------------------------!
    INTEGER, PARAMETER :: IMAX=200                          ! max number of cells
    INTEGER, PARAMETER :: NMAX=1e6                          ! max number of time steps
    INTEGER            :: i, j, k, n, n_out
    REAL               :: xL, xR, dx, dx2                   ! domain boundaries on x-axis
    REAL               :: yT, yB, dy, dy2                   ! domain boundaries on y-axis
    REAL               :: xD                                ! x-coord where the discontinuity is
    REAL               :: kidx2, kidy2
    REAL               :: time, tend, dt                    ! initial and final time and timestep   
    REAL, ALLOCATABLE  :: x(:), y(:)                        ! coordinates array
    REAL, ALLOCATABLE  :: T(:,:), Tnew(:,:)                 ! Temperature IMAX*IMAX matrices
    REAL               :: kappa, CFL, beta                  ! Problem dipendent coefficient
    REAL               :: TL=100.0, TR=50.0                 ! boundary conditions
    REAL               :: t0, t1, totTime                   ! CPU time
    !                                                       !
    INTEGER, PARAMETER :: NCPU=4                            ! Total number of CPU to use
    INTEGER            :: NPROCS                            ! Total number of CPU avaiable
    !-------------------------------------------------------!
    !
    ! Prompt
    PRINT *, ' +----------------------------------------------------------------+ '
    PRINT *, ' |                   HEAT 2D WITH EXPLICIT SOLVER                 | '
    !
    IF(NCPU.GT.1) THEN
        PRINT *, ' +----------------------------------------------------------------+ '
        PRINT *, ' | Program compiled using OpenMP directives [PARALLEL]            | '
        PRINT *, ' +----------------------------------------------------------------+ '
    ELSE
        PRINT *, ' +----------------------------------------------------------------+ '
        PRINT *, ' | Program compiled using OpenMP directives [SERIAL]              | '
        PRINT *, ' +----------------------------------------------------------------+ '
    ENDIF
    !
    NPROCS = OMP_GET_NUM_PROCS()    ! get the total number of CPU avaiable
    CALL OMP_SET_NUM_THREADS(NCPU)  ! set the number of CPU used in parallel simulation
    !
    ! Allocate memory
    ALLOCATE( x(IMAX), y(IMAX) )
    ALLOCATE( T(0:IMAX+1,0:IMAX+1) )
    ALLOCATE( Tnew(0:IMAX+1,0:IMAX+1) )
    !
    ! Computational domain
    xL  = 0.0    
    xR  = 2.0
    yB  = 0.0
    yT  = 2.0
    xD  = 1.0                   ! location for discontinuity
    dx  = (xR-xL)/REAL(IMAX-1)
    dy  = (yT-yB)/REAL(IMAX-1)
    dx2 = dx**2
    dy2 = dy**2
    !
    ! Initialization of coords arrays
    x(1)    = xL
    x(IMAX) = xR
    y(1)    = yB
    y(IMAX) = yT
    DO i = 1, IMAX-1
        x(i+1) = x(i) + dx
        y(i+1) = y(i) + dy
    ENDDO
    !
    ! Boundary conditions
    TL    = 100.0
    TR    =  50.0
    kappa =   1.0  ! Heat conduction coefficient
    CFL   =   1.0  ! Courant-Friedrichs-Lewy coefficient
    beta  =  0.45  ! stability condition for explicit method
    !
    ! Initial conditions
    DO j = 1, IMAX
        IF(x(j).LE.xD) THEN
            T(:,j) = TL
        ELSE
            T(:,j) = TR
        ENDIF
    ENDDO
    T(:,0)      = TL
    T(:,IMAX+1) = TR
    Tnew        = T
    !
    ! Plot the initial condition
    CALL DataOutput(0,IMAX,x,y,T(1:IMAX,1:IMAX),0.0)
    !
    ! MAIN COMPUTATION: time loop
    !
    ! Get initial wct
    t0 = OMP_GET_WTIME()
    !
!$OMP PARALLEL PRIVATE(time, dt, n, kidx2, kidy2)
    time  = 0.0     
    tend  = 0.05
    dt    = CFL*beta*( (dx2*dy2)/( kappa*(dx2+dy2) ) )
    kidx2 = dt*kappa/dx2
    kidy2 = dt*kappa/dy2
    !
    DO n = 1, NMAX
        !
        ! time step check
        IF((time+dt).GT.tend) THEN
            dt    = tend - time
            kidx2 = dt*kappa/dx2
            kidy2 = dt*kappa/dy2
        ENDIF
        IF(time.GE.tend) THEN
            n_out = n
            EXIT
        ENDIF
        !
        ! Finite difference scheme for the solution of 2D heat eqn
        ! Explicit method
!$OMP DO
        DO j = 1, IMAX
            DO i = 1, IMAX
                Tnew(i,j) = T(i,j) + kidx2 * ( T(i+1,j) - 2.0*T(i,j) + T(i-1,j) ) &
                                   + kidy2 * ( T(i,j+1) - 2.0*T(i,j) + T(i,j-1) )
            ENDDO !i
        ENDDO !j
!$OMP ENDDO
        !
        ! Overwrite the current solution:
        ! - update the boundary
!$OMP SINGLE
        Tnew(0,:)      = Tnew(1,:)
        Tnew(IMAX+1,:) = Tnew(IMAX,:)
        Tnew(:,0)      = Tnew(:,1)
        Tnew(:,IMAX+1) = Tnew(:,IMAX)
!$OMP END SINGLE NOWAIT
        !
        ! - update current solution
!$OMP WORKSHARE
        T = Tnew
!$OMP END WORKSHARE NOWAIT
        !        
        ! Update time
        time = time + dt
        !
    ENDDO !n
    !
!$OMP END PARALLEL
    !
    ! Get final wct
    t1 = OMP_GET_WTIME()
    !
    ! Plot the results
    CALL DataOutput(n_out,IMAX,x,y,T(1:IMAX,1:IMAX),time)
    !
    ! Print computational time on screen
    PRINT *, ' +----------------------------------------------------------------+ '
    PRINT *, ' | - IMAX value: ', IMAX
    PRINT *, ' | - Total computational time: ', t1 - t0, ' seconds.'
    PRINT *, ' | - Number of processors in use: ', NCPU
    PRINT *, ' | - Total number of processors avaiable: ', NPROCS
    PRINT *, ' |                                                                | '
    PRINT *, ' |             Program finished successfully. Bye :-)             | '
    PRINT *, ' +----------------------------------------------------------------+ '
    !
    PAUSE
    !
    ! Empty memory
    DEALLOCATE( x, y, T, Tnew )
    !    
END PROGRAM Heat2D
    