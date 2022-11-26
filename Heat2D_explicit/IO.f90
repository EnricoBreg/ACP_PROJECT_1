SUBROUTINE DataOutput(timestep,N,x,y,T,time)
    IMPLICIT NONE
    INTEGER            :: timestep, N
    REAL               :: x(N), y(N), T(N,N), time
    !
    INTEGER            :: i, j, DataUnit
    CHARACTER(LEN=10)  :: citer
    CHARACTER(LEN=200) :: IoFileName
    !
    INTENT(IN)         :: timestep, x, y, T, time
    !
    
    WRITE(*,*) ' | Plotting data at iteration = ', timestep
    
    WRITE(citer, '(I7.7)') timestep                     ! convert timestep integer to string
    IoFileName = 'Heat2D-'//TRIM(citer)//'.dat'         ! create file name
    DataUnit   = 100                                    ! to write on file
    
    ! open/create the file
    OPEN(UNIT=DataUnit, FILE=TRIM(IoFileName), STATUS='UNKNOWN', ACTION='WRITE')
    
    ! write on file
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
