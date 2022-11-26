SUBROUTINE matop2D(N,Ap,p,kappa,dt,dx2,dy2)
    IMPLICIT NONE
    !--------------------------------------------------------!
    ! Argument list
    INTEGER     :: N
    REAL        :: Ap(0:N+1,0:N+1), p(0:N+1,0:N+1)
    REAL        :: kappa, dt, dx2, dy2
    ! Local variables
    INTEGER     :: i, j
    REAL        :: a, b, c, d, e
    !--------------------------------------------------------!
    !
!$OMP PARALLEL
    !
    a = -kappa*dt/dy2
    b = -kappa*dt/dx2
    c = -kappa*dt/dy2
    d = -kappa*dt/dx2
    e = 1.0 + 2.0*kappa*dt/dx2 + 2.0*kappa*dt/dy2
    !
!$OMP DO
    DO j = 1, N
        DO i = 1, N
            Ap(i,j) = a*p(i,j-1) + b*p(i+1,j) + c*p(i,j+1) + d*p(i-1,j) + e*p(i,j)
        ENDDO !i
    ENDDO !j
!$OMP END DO
    !
!$OMP SINGLE
    ! Update the boundaries (update the ghost cells)
    Ap(:,0)   = Ap(:,1)  ! first column 
    Ap(:,N+1) = Ap(:,N)  ! last column
    Ap(0,:)   = Ap(1,:)  ! first row
    Ap(N+1,:) = Ap(N,:)  ! last row
!$OMP END SINGLE NOWAIT
    !
!$OMP END PARALLEL
    !
END SUBROUTINE matop2D
