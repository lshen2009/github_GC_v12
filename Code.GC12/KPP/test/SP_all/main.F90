program main
  USE gckpp_Parameters
  USE gckpp_JacobianSP
  USE initialize
  USE gckpp_LinearAlgebra

  IMPLICIT NONE

  REAL(kind=dp) :: JVS(LU_NONZERO)
  REAL(kind=dp) :: X(NVAR),Y1(NVAR),Y2(NVAR),Y3(NVAR)

  CALL initialize_1D(JVS)
  CALL initialize_1D(X)
  Y1=X
  CALL KppSolve ( JVS, Y1)
  print *,Y1(11:16)

  Y2=X
  CALL KppSolveIndirect(JVS,Y2)
  print *,Y2(11:16)

end program
