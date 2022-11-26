        !COMPILER-GENERATED INTERFACE MODULE: Thu Nov 10 15:15:07 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DATAOUTPUT__genmod
          INTERFACE 
            SUBROUTINE DATAOUTPUT(TIMESTEP,N,X,Y,T,TIME)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4), INTENT(IN) :: TIMESTEP
              REAL(KIND=4), INTENT(IN) :: X(N)
              REAL(KIND=4), INTENT(IN) :: Y(N)
              REAL(KIND=4), INTENT(IN) :: T(N,N)
              REAL(KIND=4), INTENT(IN) :: TIME
            END SUBROUTINE DATAOUTPUT
          END INTERFACE 
        END MODULE DATAOUTPUT__genmod
