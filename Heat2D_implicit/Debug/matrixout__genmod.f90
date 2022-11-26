        !COMPILER-GENERATED INTERFACE MODULE: Thu Nov 10 15:15:06 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MATRIXOUT__genmod
          INTERFACE 
            SUBROUTINE MATRIXOUT(M,N)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=4), INTENT(IN) :: M(N,N)
            END SUBROUTINE MATRIXOUT
          END INTERFACE 
        END MODULE MATRIXOUT__genmod
