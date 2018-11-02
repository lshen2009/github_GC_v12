program main
  USE gckpp_Parameters
  USE gckpp_JacobianSP
  USE initialize
  USE gckpp_LinearAlgebra

  IMPLICIT NONE

  REAL :: JVS(LU_NONZERO)
  REAL :: X(NVAR),Y1(NVAR),Y2(NVAR),Y3(NVAR)

  CALL initialize_1D(JVS)
  CALL initialize_1D(X)
  print *, X
  Y1=X
  CALL KppSolve ( JVS, Y1)
  print *,Y1

  Y2=X
  CALL KppSolveIndirect(JVS,Y2)
  print *,Y2

  Y3=X
  CALL KppSolveTRIndirect(JVS,Y3)
  print *,Y3
end program
