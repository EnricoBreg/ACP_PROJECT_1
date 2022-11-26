PROGRAM Heat2D
    USE OMP_LIB
    IMPLICIT NONE
    !-------------------------------------------------!
    INTEGER             :: IMAX                       ! max number of cells
    INTEGER, PARAMETER  :: NMAX=1e6                   ! max number of time steps
    INTEGER             :: i, j, n, n_out          
    REAL                :: xL, xR, dx, dx2            ! domain boundaries on x-axis
    REAL                :: yT, yB, dy, dy2            ! domain boundaries on y-axis
    REAL                :: xD                         ! x-coord where the discontinuity is
    REAL                :: kidx2, kidy2
    REAL                :: time, tend, dt
    REAL, ALLOCATABLE   :: x(:), y(:)                 ! coordinates array
    REAL, ALLOCATABLE   :: T(:,:), Tnew(:,:)
    REAL                :: kappa, CFL, beta
    REAL                :: TL=100.0, TR=50.0          ! boundary conditions
    !
    CHARACTER(LEN=200)  :: TestName                   ! name for the test
    REAL                :: t0, t1
    INTEGER, PARAMETER  :: NCPU=1
    INTEGER             :: NPROCS
    !-------------------------------------------------!
    !
    PRINT *, ' +----------------------------------------------------------------+ '
    PRINT *, ' |              HEAT 2D WITH CONJUGATE GRADIENT METHOD            | '
    !
#ifdef PARALLEL_RUN !-----------------------------------------------------------------
    ! Initialize OpenMP
    !$ PRINT *, ' +----------------------------------------------------------------+ '
    !$ PRINT *, ' | Program compiled using OpenMP directives.                      | ' 
    !$ PRINT *, ' +----------------------------------------------------------------+ '
    NPROCS = OMP_GET_NUM_PROCS()    ! get the number of total CPU avaiable
    CALL OMP_SET_NUM_THREADS(NCPU)  ! set the number of CPU for computation
    !
    TestName = 'Heat2D-CG-parallel' ! set the name of the test
    !
#else !-------------------------------------------------------------------------------
    PRINT *, ' +----------------------------------------------------------------+ '
    PRINT *, ' | Running SERIAL simulation.                                     | '
    PRINT *, ' +----------------------------------------------------------------+ '
    !
    TestName = 'Heat2D-CG-serial'   ! set the name of test
    !
#endif !------------------------------------------------------------------------------
    !
    ! Computational domain
    IMAX = 500
    xL   = 0.0    
    xR   = 2.0
    yB   = 0.0
    yT   = 2.0
    xD   = 1.0                   ! location of the discontinuity
    dx   = (xR-xL)/REAL(IMAX-1)
    dy   = (yT-yB)/REAL(IMAX-1)
    !
    ! Boundary conditions
    TL = 100.0
    TR =  50.0
    !
    kappa = 1.0  ! heat conduction coefficient
    CFL   = 1.0  ! Courant-Friedrichs-Lewy coefficient
    beta  = 2.0  ! any value for implicit method is ok
    !
    dx2 = dx**2
    dy2 = dy**2
    !
    ! Allocate memory
    ALLOCATE( x(IMAX), y(IMAX)        )
    ALLOCATE( T(0:IMAX+1,0:IMAX+1)    )
    ALLOCATE( Tnew(0:IMAX+1,0:IMAX+1) )
    !
    ! Initialization of coords arrays
    x(1)    = xL
    y(1)    = yB
    x(IMAX) = xR
    y(IMAX) = yT
    DO i = 1, IMAX-1
        x(i+1) = x(i) + dx
        y(i+1) = y(i) + dy
    ENDDO
    !
    ! Set initial condition for temperature
    DO j = 1, IMAX
        IF(x(j).LE.xD) THEN
            T(:,j) = TL
        ELSE
            T(:,j) = TR
        ENDIF
    ENDDO
    T(:,0)      = TL
    T(:,IMAX+1) = TR
    !
    Tnew = T ! Init the new temperature
    !
    ! Plot initial condition
    CALL DataOutput(TestName,0,IMAX,x,y,T(1:IMAX,1:IMAX),0.0)
    !
    ! Get CPU time
    t0 = OMP_GET_WTIME()
    !   
    ! MAIN COMPUTATION: time loop
    time = 0.0     
    tend = 0.05
    dt   = CFL*beta*( (dx2*dy2)/( kappa*(dx2+dy2) ) )
    !
    DO n = 1, NMAX
        !
        ! time step check
        IF((time+dt).GT.tend) THEN
            dt = tend - time
        ENDIF
        IF(time.GE.tend) THEN
            n_out = n
            EXIT
        ENDIF
        !
        ! Implicit finite difference scheme for the solution of 2D heat eqn
#ifdef PARALLEL_RUN !---------------------------------------------------------------------
        CALL CG2D_v2(IMAX,Tnew,T,kappa,dt,dx2,dy2)  ! call parallel version of CG method
#else   !---------------------------------------------------------------------------------
        CALL CG2D(IMAX,Tnew,T,kappa,dt,dx2,dy2)     ! call standard version of CG method
#endif  !---------------------------------------------------------------------------------
        !
        T    = Tnew      ! overwrite the current solution
        time = time + dt ! update time
        !
    ENDDO !n
    !
    ! Get CPU time
    t1 = OMP_GET_WTIME()
    !
    ! Plot the results
    CALL DataOutput(TestName,n_out,IMAX,x,y,T(1:IMAX,1:IMAX),time)
    !
    ! Print computational time on screen
    PRINT *, ' +----------------------------------------------------------------+ '
    PRINT *, ' | - IMAX value: ', IMAX
    PRINT *, ' | - Total computational time: ', t1 - t0, ' seconds.'
#ifdef PARALLEL_RUN !----------------------------------------------------------------
    PRINT *, ' | - Number of processors in use: ', NCPU
    PRINT *, ' | - Total number of processors avaiable: ', NPROCS
#endif !-----------------------------------------------------------------------------
    PRINT *, ' |                                                                  '
    PRINT *, ' |             Program finished successfully. Bye :-)               '
    PRINT *, ' +----------------------------------------------------------------+ '
    !           
    PAUSE
    !
    ! Empty memory
    DEALLOCATE( x, y, T, Tnew )
    !    
END PROGRAM Heat2D
