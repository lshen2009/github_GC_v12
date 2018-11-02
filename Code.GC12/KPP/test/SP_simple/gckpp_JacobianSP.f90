! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Sparse Jacobian Data Structures File
! 
! Generated by KPP-2.2.4_gc symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : gckpp_JacobianSP.f90
! Time                 : Mon Jul 30 10:56:32 2018
! Working directory    : /n/regal/mickley_lab/lshen/Code.v11-02c_HGI/KPP/test
! Equation file        : gckpp.kpp
! Output root filename : gckpp
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE gckpp_JacobianSP

  PUBLIC
  SAVE


! Sparse Jacobian Data


  INTEGER, PARAMETER, DIMENSION(13) :: LU_IROW = (/ &
       1,  1,  1,  2,  2,  2,  3,  3,  4,  4,  5,  5, & ! index 1 - 12
       5 /)

  INTEGER, PARAMETER, DIMENSION(13) :: LU_ICOL = (/ &
       1,  3,  5,  2,  4,  5,  3,  5,  4,  5,  3,  4, & ! index 1 - 12
       5 /)

  INTEGER, PARAMETER, DIMENSION(6) :: LU_CROW = (/ &
       1,  4,  7,  9, 11, 14 /)

  INTEGER, PARAMETER, DIMENSION(6) :: LU_DIAG = (/ &
       1,  4,  7,  9, 13, 14 /)


END MODULE gckpp_JacobianSP

