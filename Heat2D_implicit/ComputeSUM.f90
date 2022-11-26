!-------------------------------------------------------------!
!     Compute the result of intrinsic function: SUM(v*u)      !
!-------------------------------------------------------------!
SUBROUTINE ComputeSUM(N,sum_vu,v,u)
    IMPLICIT NONE
    !--------------------------------------------------------!
    ! Argument list
    INTEGER     :: N
    REAL        :: sum_vu
    REAL        :: v(0:N+1,0:N+1), u(0:N+1,0:N+1)
    ! Local variables
    INTEGER     :: i, j
    !--------------------------------------------------------!
    !
    sum_vu = 0.0 ! init the variable
    !
!$OMP PARALLEL DO REDUCTION(+:sum_vu)
    DO j = 0, N+1
        DO i = 0, N+1
            sum_vu = sum_vu + v(i,j)*u(i,j)
        ENDDO !i
    ENDDO !j
!$OMP END PARALLEL DO
    !
END SUBROUTINE ComputeSUM