!
!---------------------------------------------------------------------------------------!    
!                                       CG2D vers 3                                     !
!---------------------------------------------------------------------------------------!    
!   Implementazione del metodo del Gradiente Coniugato con l'uso di un ciclo DO WHILE   !    
!                               invece che di un ciclo DO.                              !    
!---------------------------------------------------------------------------------------!    
SUBROUTINE CG2D_v3(N,x,b,kappa,dt,dx2,dy2)
    IMPLICIT NONE
    !---------------------------------------------------!
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
    REAL            :: sum_pAp, res
    !---------------------------------------------------!
    !
    KMAX = N*N  ! set maximum number of interations
    x = b
    CALL matop2D(N,Ax,x,kappa,dt,dx2,dy2)
    r = b - Ax
    p = r
    CALL ComputeSUM(N,alphak,r,r)
    !
    k   = 1
    res = SQRT(alphak)
    !
    DO WHILE( (k.LE.KMAX) .AND. (res.GT.tol) )
        !
        CALL matop2D(N,Ap,p,kappa,dt,dx2,dy2)  
        CALL ComputeSUM(N,sum_pAp,p,Ap)
        lambda = alphak/sum_pAp
        x = x + lambda*p
        r = r - lambda*Ap
        CALL ComputeSUM(N,alpha,r,r)
        p = r + alpha/alphak*p
        alphak = alpha       
        !
        res = SQRT(alphak) ! update the residual
        k   = k + 1        ! update the number of iters  
        !
    ENDDO !k
    !
    IF(res.LT.tol) THEN
        PRINT *, ' | CG iter = ', k, ' CG res = ', res
    ENDIF
    !
    IF(k.GE.KMAX) THEN
        PRINT *, ' | ERROR. CG method did not converge! '
        PRINT *, ' | CG res = ', SQRT(alphak)
        PAUSE
        STOP
    ENDIF
    !
END SUBROUTINE CG2D_v3    
