!---------------------------------------------------------------------------------------!    
!                                       CG2D vers 2                                     !
!---------------------------------------------------------------------------------------!    
!      Implementation of conjugate gradient method using SUM instrinsic function        !    
!---------------------------------------------------------------------------------------!       
SUBROUTINE CG2D_v2(N,x,b,kappa,dt,dx2,dy2)
    IMPLICIT NONE
    !--------------------------------------------------!
    ! Argument list
    INTEGER         :: N                                ! dimension of the problem
    REAL            :: x(0:N+1,0:N+1)                           
    REAL            :: b(0:N+1,0:N+1)                   ! rhs
    REAL            :: kappa, dt, dx2, dy2              ! problem dipendent variables
    ! Local variables
    INTEGER         :: i, j, k, KMAX                             
    REAL            :: Ax(0:N+1,0:N+1), r(0:N+1,0:N+1)
    REAL            :: Ap(0:N+1,0:N+1), p(0:N+1,0:N+1)  
    REAL            :: alpha, alphak, lambda
    REAL, PARAMETER :: tol=1e-12
    REAL            :: sum_rr, sum_pAp
    !--------------------------------------------------!
    !
    KMAX = N*N  ! CG maximum number of interations
!$OMP PARALLEL WORKSHARE
    x = b
!$OMP END PARALLEL WORKSHARE
    CALL matop2D(N,Ax,x,kappa,dt,dx2,dy2)
!$OMP PARALLEL WORKSHARE
    r = b - Ax
    p = r
!$OMP END PARALLEL WORKSHARE
    CALL ComputeSUM(N,sum_rr,r,r)
    alphak = sum_rr
    !
    DO k = 1, KMAX
        !
        IF(SQRT(alphak).LT.tol) THEN
            PRINT *, ' | CG iter: ', k, ' CG res = ', SQRT(alphak)
            RETURN
        ENDIF
        CALL matop2D(N,Ap,p,kappa,dt,dx2,dy2)
        CALL ComputeSUM(N,sum_pAp,p,Ap)
        lambda = alphak/sum_pAp
!$OMP PARALLEL WORKSHARE
        x      = x + lambda*p
        r      = r - lambda*Ap
!$OMP END PARALLEL WORKSHARE
        CALL ComputeSUM(N,sum_rr,r,r)
        alpha = sum_rr
!$OMP PARALLEL WORKSHARE
        p      = r + alpha/alphak*p
!$OMP END PARALLEL WORKSHARE
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
END SUBROUTINE CG2D_v2
