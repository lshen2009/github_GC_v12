program main
  USE gckpp_Parameters
  USE gckpp_JacobianSP
  USE initialize
  USE gckpp_LinearAlgebra

  IMPLICIT NONE

  REAL(kind=dp) :: JVS(LU_NONZERO),JVS2(LU_NONZERO)
  REAL(kind=dp) :: JVS_orig(LU_NONZERO)
  REAL(kind=dp) :: X(NVAR)
  LOGICAL :: ind(LU_NONZERO)
  INTEGER:: IER,i,j,k,ind2(LU_NONZERO)
  REAL, ALLOCATABLE :: new(:)
  !INTEGER, ALLOCATABLE :: LU_IROW2(:),LU_ICOL2(:),LU_DIAG2(:)
  INTEGER :: LU_DIAG2(NVAR+1),LU_CROW2(NVAR+1),LU_ICOL2(LU_NONZERO)

  CALL initialize_1D(JVS_orig)
  CALL initialize_1D(X)
  IER=0

  JVS=JVS_orig
  CALL KppDecomp2(JVS,IER)
  JVS2=JVS_orig
  CALL KppDecomp3(JVS2,IER)
  CALL KppSolveIndirect(JVS,X)
end program
