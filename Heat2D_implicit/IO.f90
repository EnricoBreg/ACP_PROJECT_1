SUBROUTINE DataOutput(TestName,timestep,N,x,y,T,time)
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list
    INTEGER, INTENT(IN)            :: timestep, N
    REAL,    INTENT(IN)            :: x(N), y(N), T(N,N), time
    CHARACTER(LEN=200), INTENT(IN) :: TestName
    ! Local variables
    INTEGER                        :: i, j, DataUnit
    CHARACTER(LEN=10)              :: citer
    CHARACTER(LEN=200)             :: IoFileName
    !-------------------------------------------------------------------------!
    
    WRITE(*,*) ' | Plotting data at iteration = ', timestep
    
    WRITE(citer, '(I7.7)') timestep                       ! convert timestep integer to string
    IoFileName = TRIM(TestName)//'-'//TRIM(citer)//'.dat'         ! create file name
    DataUnit   = 100                                      ! to write on file
    
    ! open/create the file
    OPEN(UNIT=DataUnit, FILE=TRIM(IoFileName), STATUS='UNKNOWN', ACTION='WRITE')
    
    ! write problem dimension
    WRITE(DataUnit, *) N
    
    ! write x-coordinates
    DO i = 1, N
        WRITE(DataUnit, *) x(i)
    ENDDO
    
    ! write y-coordinates
    DO j = 1, N
        WRITE(DataUnit, *) y(j)
    ENDDO
    
    ! write solution
    DO j = 1, N
        DO i = 1, N
            WRITE(DataUnit, *) T(i,j)
        ENDDO !i
    ENDDO !j
    
    ! close file
    CLOSE(DataUnit)
    
END SUBROUTINE DataOutput    
