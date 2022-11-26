!    
!---------------------------------------------------------------------------------------!    
!                                       CG2D vers 1                                     !    
!---------------------------------------------------------------------------------------!    
SUBROUTINE CG2D(N,x,b,kappa,dt,dx2,dy2)
    !
    IMPLICIT NONE
    !--------------------------------------------------!
    ! Argument list
    INTEGER         :: N                                ! dimension of the problem
    REAL            :: x(0:N+1,0:N+1)                           
    REAL            :: b(0:N+1,0:N+1)                   ! rhs
    REAL            :: kappa, dt, dx2, dy2              ! problem dipendent variables
    ! Local variables
    INTEGER         :: k, KMAX                             
    REAL            :: Ax(0:N+1,0:N+1), r(0:N+1,0:N+1)
    REAL            :: Ap(0:N+1,0:N+1), p(0:N+1,0:N+1)
    REAL            :: alpha, alphak, lambda
    REAL, PARAMETER :: tol=1e-12
    !--------------------------------------------------!
    !
    KMAX   = N*N       ! maximum number of interations
    x      = b
    CALL matop2D(N,Ax,x,kappa,dt,dx2,dy2)
    r      = b - Ax
    p      = r
    alphak = SUM(r*r)
    !
    DO k = 1, KMAX
        !        
        IF(SQRT(alphak).LT.tol) THEN
            PRINT *, ' | CG iter: ', k, ' CG res = ', SQRT(alphak)
            RETURN
        ENDIF
        !
        CALL matop2D(N,Ap,p,kappa,dt,dx2,dy2)
        lambda = alphak/SUM(p*Ap)
        x      = x + lambda*p
        r      = r - lambda*Ap
        alpha  = SUM(r*r)
        p      = r + alpha/alphak*p
        alphak = alpha
        !
    ENDDO !k
    !
    IF(k.GE.KMAX) THEN
        PRINT *, ' ERROR. CG method did not converge! '
        PRINT *, ' CG res = ', SQRT(alphak)
        PAUSE
        STOP
    ENDIF
    !
END SUBROUTINE CG2D
