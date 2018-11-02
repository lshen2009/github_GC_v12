MODULE gckpp_LinearAlgebra

  USE gckpp_Parameters
  USE gckpp_JacobianSP

  IMPLICIT NONE

CONTAINS
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE KppSolveIndirect( JVS, X )
  USE gckpp_Parameters
  USE gckpp_JacobianSP
      INTEGER  :: i, j
      REAL     :: JVS(LU_NONZERO), X(NVAR), sum

      DO i=1,NVAR
         DO j = LU_CROW(i), LU_DIAG(i)-1 
             X(i) = X(i) - JVS(j)*X(LU_ICOL(j));
         END DO  
      END DO

      DO i=NVAR,1,-1
        sum = X(i);
        DO j = LU_DIAG(i)+1, LU_CROW(i+1)-1
          sum = sum - JVS(j)*X(LU_ICOL(j));
        END DO
        X(i) = sum/JVS(LU_DIAG(i));
      END DO
END SUBROUTINE KppSolveIndirect

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE KppSolveTRIndirect( JVS, X )
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        Complex sparse solve transpose subroutine using indirect addressing
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  USE gckpp_Parameters
  USE gckpp_JacobianSP

      INTEGER       :: i, j
      REAL :: JVS(LU_NONZERO), X(NVAR)

      DO i=1,NVAR
        X(i) = X(i)/JVS(LU_DIAG(i))
        ! subtract all nonzero elements in row i of JVS from X
        DO j=LU_DIAG(i)+1,LU_CROW(i+1)-1
          X(LU_ICOL(j)) = X(LU_ICOL(j))-JVS(j)*X(i)
        END DO
      END DO

      DO i=NVAR, 1, -1
        ! subtract all nonzero elements in row i of JVS from X
        DO j=LU_CROW(i),LU_DIAG(i)-1
          X(LU_ICOL(j)) = X(LU_ICOL(j))-JVS(j)*X(i)
        END DO
      END DO
      
END SUBROUTINE KppSolveTRIndirect

SUBROUTINE KppSolve ( JVS, X )

  REAL          :: JVS(LU_NONZERO)
  REAL          :: X(NVAR)

  X(5) = X(5)-JVS(11)*X(3)-JVS(12)*X(4)
  X(5) = X(5)/JVS(13)
  X(4) = (X(4)-JVS(10)*X(5))/(JVS(9))
  X(3) = (X(3)-JVS(8)*X(5))/(JVS(7))
  X(2) = (X(2)-JVS(5)*X(4)-JVS(6)*X(5))/(JVS(4))
  X(1) = (X(1)-JVS(2)*X(3)-JVS(3)*X(5))/(JVS(1))
      
END SUBROUTINE KppSolve


END MODULE gckpp_LinearAlgebra
