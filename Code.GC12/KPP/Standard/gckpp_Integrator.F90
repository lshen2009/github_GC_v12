MODULE gckpp_Integrator

  USE gckpp_Parameters, ONLY: NVAR, NFIX, NSPEC
  USE gckpp_JacobianSP
  USE gckpp_Global
  IMPLICIT NONE
  PUBLIC
  SAVE
  
!~~~>  Statistics on the work performed by the Rosenbrock method
  INTEGER, PARAMETER :: Nfun=1, Njac=2, Nstp=3, Nacc=4, &
                        Nrej=5, Ndec=6, Nsol=7, Nsng=8, &
                        Ntexit=1, Nhexit=2, Nhnew = 3

CONTAINS

SUBROUTINE INTEGRATE_1( TIN, TOUT, LS_type,LS_NSEL, LS_NDEL,&
  Prate, Lrate, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )

   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: TIN  ! Start Time
   REAL(kind=dp), INTENT(IN) :: TOUT ! End Time
   
   ! Optional input parameters and statistics
   INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   REAL(kind=dp), INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   REAL(kind=dp), INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U
   INTEGER,		  INTENT(IN) :: LS_type,LS_NSEL,LS_NDEL
   REAL(kind=dp), INTENT(IN) :: Prate(NVAR),Lrate(NVAR)
   
   REAL(kind=dp) :: RCNTRL(20), RSTATUS(20),deltaT
   INTEGER       :: ICNTRL(20), ISTATUS(20), IERR
  
   
   INTEGER, SAVE :: Ntotal = 0

   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.0_dp
   ISTATUS(:) = 0
   RSTATUS(:) = 0.0_dp
   deltaT =  TOUT-TIN
   !~~~> fine-tune the integrator:
   ICNTRL(1) = 0	! 0 - non-autonomous, 1 - autonomous
   ICNTRL(2) = 0	! 0 - vector tolerances, 1 - scalars

   ! If optional parameters are given, and if they are >0, 
   ! then they overwrite default settings. 
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF
   		
		WHERE(Lrate<=(0.01/deltaT))
			VAR=VAR+deltaT*(Prate-Lrate*VAR)
		ELSEWHERE
			VAR=Prate/Lrate+(VAR-Prate/Lrate)*EXP(-Lrate*deltaT)
		END WHERE

   
   !~~~> Debug option: show no of steps
   ! Ntotal = Ntotal + ISTATUS(Nstp)
   ! PRINT*,'NSTEPS=',ISTATUS(Nstp),' (',Ntotal,')','  O3=', VAR(ind_O3)

   STEPMIN = RSTATUS(Nhexit)
   ! if optional parameters are given for output they 
   ! are updated with the return information
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE INTEGRATE_1

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE INTEGRATE_2( TIN, TOUT, LS_type,LS_NSEL, LS_NDEL,&
  Prate, Lrate, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )

   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: TIN  ! Start Time
   REAL(kind=dp), INTENT(IN) :: TOUT ! End Time
   
   ! Optional input parameters and statistics
   INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   REAL(kind=dp), INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   REAL(kind=dp), INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U
   INTEGER,		  INTENT(IN) :: LS_type,LS_NSEL,LS_NDEL
   REAL(kind=dp), INTENT(IN) :: Prate(NVAR),Lrate(NVAR)
   
   REAL(kind=dp) :: RCNTRL(20), RSTATUS(20),deltaT
   INTEGER       :: ICNTRL(20), ISTATUS(20), IERR
   REAL(kind=dp) :: VAR_selected(LS_NSEL),VAR_deleted(LS_NDEL),LS_P(LS_NDEL),LS_L(LS_NDEL)
   
   INTEGER, SAVE :: Ntotal = 0

   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.0_dp
   ISTATUS(:) = 0
   RSTATUS(:) = 0.0_dp
   deltaT =  TOUT-TIN
   !~~~> fine-tune the integrator:
   ICNTRL(1) = 0	! 0 - non-autonomous, 1 - autonomous
   ICNTRL(2) = 0	! 0 - vector tolerances, 1 - scalars

   ! If optional parameters are given, and if they are >0, 
   ! then they overwrite default settings. 
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF
   
        VAR_selected=VAR(select_ind_2)
	    VAR_deleted=VAR(delete_ind_2)
		LS_P=Prate(delete_ind_2)
		LS_L=Lrate(delete_ind_2)
        CALL Rosenbrock_2(LU_NSEL_2,VAR_selected,TIN,TOUT,ATOL,RTOL,&
          RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		  LU_NONZERO_2,LU_NSEL_2,LU_CROW_2,LU_DIAG_2,LU_IROW_2,LU_ICOL_2, LS_type)		
		WHERE(LS_L<=(0.01/deltaT))
			VAR_deleted=VAR_deleted+deltaT*(LS_P-LS_L*VAR_deleted)
		ELSEWHERE
			VAR_deleted=LS_P/LS_L+(VAR_deleted-LS_P/LS_L)*EXP(-LS_L*deltaT)
		END WHERE

	    VAR(select_ind_2)=VAR_selected
	    VAR(delete_ind_2)=VAR_deleted
   
   !~~~> Debug option: show no of steps
   ! Ntotal = Ntotal + ISTATUS(Nstp)
   ! PRINT*,'NSTEPS=',ISTATUS(Nstp),' (',Ntotal,')','  O3=', VAR(ind_O3)

   STEPMIN = RSTATUS(Nhexit)
   ! if optional parameters are given for output they 
   ! are updated with the return information
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE INTEGRATE_2

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE Rosenbrock_2(N,Y,Tstart,Tend, &
           AbsTol,RelTol,              &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		   LS_LU_NONZERO,LS_NSEL,LS_LU_CROW,LS_LU_DIAG,LS_LU_IROW,LS_LU_ICOL, LS_type)
  USE gckpp_Parameters
  USE gckpp_LinearAlgebra
  IMPLICIT NONE

   INTEGER,       INTENT(IN)    :: N,LS_LU_NONZERO,LS_NSEL,LS_type   
   INTEGER,INTENT(IN)::LS_LU_CROW(LS_NSEL+1),LS_LU_DIAG(LS_NSEL+1)
   INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO),LS_LU_ICOL(LS_LU_NONZERO)   
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
   REAL(kind=dp), INTENT(IN)    :: Tstart,Tend
   REAL(kind=dp), INTENT(IN)    :: AbsTol(N),RelTol(N)
   INTEGER,       INTENT(IN)    :: ICNTRL(20)
   REAL(kind=dp), INTENT(IN)    :: RCNTRL(20)
   INTEGER,       INTENT(INOUT) :: ISTATUS(20)
   REAL(kind=dp), INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT)   :: IERR
!~~~>  Parameters of the Rosenbrock method, up to 6 stages
   INTEGER ::  ros_S, rosMethod
   INTEGER, PARAMETER :: RS2=1, RS3=2, RS4=3, RD3=4, RD4=5, RG3=6
   REAL(kind=dp) :: ros_A(15), ros_C(15), ros_M(6), ros_E(6), &
                    ros_Alpha(6), ros_Gamma(6), ros_ELO
   LOGICAL :: ros_NewF(6)
   CHARACTER(LEN=12) :: ros_Name
!~~~>  Local variables
   REAL(kind=dp) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   REAL(kind=dp) :: Hmin, Hmax, Hstart
   REAL(kind=dp) :: Texit
   INTEGER       :: i, UplimTol, Max_no_steps
   LOGICAL       :: Autonomous, VectorTol
!~~~>   Parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp

!~~~>  Initialize statistics
   ISTATUS(1:8) = 0
   RSTATUS(1:3) = ZERO

!~~~>  Autonomous or time dependent ODE. Default is time dependent.
   Autonomous = .NOT.(ICNTRL(1) == 0)

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:N) and RelTol(1:N)
   IF (ICNTRL(2) == 0) THEN
      VectorTol = .TRUE.
      UplimTol  = N
   ELSE
      VectorTol = .FALSE.
      UplimTol  = 1
   END IF

!~~~>   Initialize the particular Rosenbrock method selected
   SELECT CASE (ICNTRL(3))
     CASE (1)
       CALL Ros2
     CASE (2)
       CALL Ros3
     CASE (3)
       CALL Ros4
     CASE (0,4)
       CALL Rodas3
     CASE (5)
       CALL Rodas4
     CASE (6)
       CALL Rang3
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(3)=',ICNTRL(3) 
       CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT

!~~~>   The maximum number of steps admitted
   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 200000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Unit roundoff (1+Roundoff>1)
   Roundoff = WLAMCH('E')

!~~~>  Lower bound on the step size: (positive value)
   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Upper bound on the step size: (positive value)
   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Starting step size: (positive value)
   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hold < FacMax
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Check if tolerances are reasonable
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO


!~~~>  CALL Rosenbrock method
   CALL ros_Integrator(Y, Tstart, Tend, Texit,   &
        AbsTol, RelTol,                          &
!  Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!  Error indicator
        IERR)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS !  SUBROUTINES internal to Rosenbrock
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 SUBROUTINE ros_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
   REAL(kind=dp), INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR
   
   IERR = Code
   PRINT * , &
     'Forced exit from Rosenbrock due to the following error:' 
     
   SELECT CASE (Code)
    CASE (-1)
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)    
      PRINT * , '--> Selected Rosenbrock method not implemented'
    CASE (-3)    
      PRINT * , '--> Hmin/Hmax/Hstart must be positive'
    CASE (-4)    
      PRINT * , '--> FacMin/FacMax/FacRej must be positive'
    CASE (-5) 
      PRINT * , '--> Improper tolerance values'
    CASE (-6) 
      PRINT * , '--> No of steps exceeds maximum bound'
    CASE (-7) 
      PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
    CASE (-8)    
      PRINT * , '--> Matrix is repeatedly singular'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT
   
   PRINT *, "T=", T, "and H=", H
     
 END SUBROUTINE ros_ErrorMsg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_Integrator (Y, Tstart, Tend, T,  &
        AbsTol, RelTol,                          &
!~~~> Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic Rosenbrock method
!      defined by ros_S (no of stages)
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

!~~~> Input: the initial condition at Tstart; Output: the solution at T
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
!~~~> Input: integration interval
   REAL(kind=dp), INTENT(IN) :: Tstart,Tend
!~~~> Output: time at which the solution is returned (T=Tend if success)
   REAL(kind=dp), INTENT(OUT) ::  T
!~~~> Input: tolerances
   REAL(kind=dp), INTENT(IN) ::  AbsTol(N), RelTol(N)
!~~~> Input: integration parameters
   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   REAL(kind=dp), INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   REAL(kind=dp), INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables
   REAL(kind=dp) :: Ynew(N), Fcn0(N), Fcn(N)
   REAL(kind=dp) :: K(N*ros_S), dFdT(N)
#ifdef FULL_ALGEBRA    
   REAL(kind=dp) :: Jac0(N,N), Ghimj(N,N)
#else
   REAL(kind=dp) :: Jac0(LS_LU_NONZERO), Ghimj(LS_LU_NONZERO)
#endif
   REAL(kind=dp) :: H, Hnew, HC, HG, Fac, Tau
   REAL(kind=dp) :: Err, Yerr(N)
   INTEGER :: Pivot(N), Direction, ioffset, j, istage
   LOGICAL :: RejectLastH, RejectMoreH, Singular
!~~~>  Local parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp
!~~~>  Locally called functions
!    REAL(kind=dp) WLAMCH
!    EXTERNAL WLAMCH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~>  Initial preparations
   T = Tstart
   RSTATUS(Nhexit) = ZERO
   H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   IF (ABS(H) <= 10.0_dp*Roundoff) H = DeltaMin

   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF
   H = Direction*H

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.

!~~~> Time loop begins below

TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) )

   IF ( ISTATUS(Nstp) > Max_no_steps ) THEN  ! Too many steps
      CALL ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
      CALL ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF

!~~~>  Limit H if necessary to avoid going beyond Tend
   H = MIN(H,ABS(Tend-T))

!~~~>   Compute the function at current time
   CALL FunTemplate_2(T,Y,Fcn0, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1

!~~~>  Compute the function derivative with respect to T
   IF (.NOT.Autonomous) THEN
      CALL ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
   END IF

!~~~>   Compute the Jacobian at current time
   CALL JacTemplate_2(T,Y,Jac0,LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)   
   ISTATUS(Njac) = ISTATUS(Njac) + 1

!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO

   CALL ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular)
   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       CALL ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S

      ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
       ioffset = N*(istage-1)

      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         !slim: CALL WCOPY(N,Fcn0,1,Fcn,1)
	 Fcn(1:N) = Fcn0(1:N)
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
         !slim: CALL WCOPY(N,Y,1,Ynew,1)
	 Ynew(1:N) = Y(1:N)
         DO j = 1, istage-1
           CALL WAXPY(N,ros_A((istage-1)*(istage-2)/2+j), &
            K(N*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL FunTemplate_2(Tau,Ynew,Fcn, LS_NSEL, LS_type)
         ISTATUS(Nfun) = ISTATUS(Nfun) + 1
       END IF ! if istage == 1 elseif ros_NewF(istage)
       !slim: CALL WCOPY(N,Fcn,1,K(ioffset+1),1)
       K(ioffset+1:ioffset+N) = Fcn(1:N)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL WAXPY(N,HC,K(N*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL WAXPY(N,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL ros_Solve(Ghimj, Pivot, K(ioffset+1))

   END DO Stage


!~~~>  Compute the new solution
   !slim: CALL WCOPY(N,Y,1,Ynew,1)
   Ynew(1:N) = Y(1:N)
   DO j=1,ros_S
         CALL WAXPY(N,ros_M(j),K(N*(j-1)+1),1,Ynew,1)
   END DO

!~~~>  Compute the error estimation
   !slim: CALL WSCAL(N,ZERO,Yerr,1)
   Yerr(1:N) = ZERO
   DO j=1,ros_S
        CALL WAXPY(N,ros_E(j),K(N*(j-1)+1),1,Yerr,1)
   END DO
   Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac

!~~~>  Check the error magnitude and adjust step size
   ISTATUS(Nstp) = ISTATUS(Nstp) + 1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
      ISTATUS(Nacc) = ISTATUS(Nacc) + 1
      !slim: CALL WCOPY(N,Ynew,1,Y,1)
      Y(1:N) = Ynew(1:N)
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H)
      END IF
      RSTATUS(Nhexit) = H
      RSTATUS(Nhnew)  = Hnew
      RSTATUS(Ntexit) = T
      RejectLastH = .FALSE.
      RejectMoreH = .FALSE.
      H = Hnew
      EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
   ELSE           !~~~> Reject step
      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (ISTATUS(Nacc) >= 1)  ISTATUS(Nrej) = ISTATUS(Nrej) + 1
   END IF ! Err <= 1

   END DO UntilAccepted

   
   END DO TimeLoop

!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful

  END SUBROUTINE ros_Integrator


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  REAL(kind=dp) FUNCTION ros_ErrorNorm ( Y, Ynew, Yerr, &
                               AbsTol, RelTol, VectorTol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

! Input arguments
   REAL(kind=dp), INTENT(IN) :: Y(N), Ynew(N), &
          Yerr(N), AbsTol(N), RelTol(N)
   LOGICAL, INTENT(IN) ::  VectorTol
! Local variables
   REAL(kind=dp) :: Err, Scale, Ymax
   INTEGER  :: i
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp

   Err = ZERO
   DO i=1,N
     Ymax = MAX(ABS(Y(i)),ABS(Ynew(i)))
     IF (VectorTol) THEN
       Scale = AbsTol(i)+RelTol(i)*Ymax
     ELSE
       Scale = AbsTol(1)+RelTol(1)*Ymax
     END IF
     Err = Err+(Yerr(i)/Scale)**2
   END DO
   Err  = SQRT(Err/N)

   ros_ErrorNorm = MAX(Err,1.0d-10)

  END FUNCTION ros_ErrorNorm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> The time partial derivative of the function by finite differences
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

!~~~> Input arguments
   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(N), Fcn0(N)
!~~~> Output arguments
   REAL(kind=dp), INTENT(OUT) :: dFdT(N)
!~~~> Local variables
   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL FunTemplate_2(T+Delta,Y,dFdT, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1
   CALL WAXPY(N,(-ONE),Fcn0,1,dFdT,1)
   CALL WSCAL(N,(ONE/Delta),dFdT,1)

  END SUBROUTINE ros_FunTimeDerivative


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_PrepareMatrix ( H, Direction, gam, &
             Jac0, Ghimj, Pivot, Singular )
! --- --- --- --- --- --- --- --- --- --- --- --- ---
!  Prepares the LHS matrix for stage calculations
!  1.  Construct Ghimj = 1/(H*ham) - Jac0
!      "(Gamma H) Inverse Minus Jacobian"
!  2.  Repeat LU decomposition of Ghimj until successful.
!       -half the step size if LU decomposition fails and retry
!       -exit after 5 consecutive fails
! --- --- --- --- --- --- --- --- --- --- --- --- ---
   IMPLICIT NONE

!~~~> Input arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) ::  Jac0(N,N)
#else
   REAL(kind=dp), INTENT(IN) ::  Jac0(LS_LU_NONZERO)
#endif   
   REAL(kind=dp), INTENT(IN) ::  gam
   INTEGER, INTENT(IN) ::  Direction
!~~~> Output arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(OUT) :: Ghimj(N,N)
#else
   REAL(kind=dp), INTENT(OUT) :: Ghimj(LS_LU_NONZERO)
#endif   
   LOGICAL, INTENT(OUT) ::  Singular
   INTEGER, INTENT(OUT) ::  Pivot(N)
!~~~> Inout arguments
   REAL(kind=dp), INTENT(INOUT) :: H   ! step size is decreased when LU fails
!~~~> Local variables
   INTEGER  :: i, ISING, Nconsecutive
   REAL(kind=dp) :: ghinv
   REAL(kind=dp), PARAMETER :: ONE  = 1.0_dp, HALF = 0.5_dp

   Nconsecutive = 0
   Singular = .TRUE.

   DO WHILE (Singular)

!~~~>    Construct Ghimj = 1/(H*gam) - Jac0
#ifdef FULL_ALGEBRA    
     !slim: CALL WCOPY(N*N,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(N*N,(-ONE),Ghimj,1)
     Ghimj = -Jac0
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(i,i) = Ghimj(i,i)+ghinv
     END DO
#else
     !slim: CALL WCOPY(LS_LU_NONZERO,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(LS_LU_NONZERO,(-ONE),Ghimj,1)
     Ghimj(1:LS_LU_NONZERO) = -Jac0(1:LS_LU_NONZERO)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(LS_LU_DIAG(i)) = Ghimj(LS_LU_DIAG(i))+ghinv
     END DO
#endif   
!~~~>    Compute LU decomposition
     CALL ros_Decomp( Ghimj, Pivot, ISING )
     IF (ISING == 0) THEN
!~~~>    If successful done
        Singular = .FALSE.
     ELSE ! ISING .ne. 0
!~~~>    If unsuccessful half the step size; if 5 consecutive fails then return
        ISTATUS(Nsng) = ISTATUS(Nsng) + 1
        Nconsecutive = Nconsecutive+1
        Singular = .TRUE.
        PRINT*,'Warning: LU Decomposition returned ISING = ',ISING
        IF (Nconsecutive <= 5) THEN ! Less than 5 consecutive failed decompositions
           H = H*HALF
        ELSE  ! More than 5 consecutive failed decompositions
           RETURN
        END IF  ! Nconsecutive
      END IF    ! ISING

   END DO ! WHILE Singular

  END SUBROUTINE ros_PrepareMatrix


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Decomp( A, Pivot, ISING )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the LU decomposition
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Inout variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(INOUT) :: A(N,N)
#else   
   REAL(kind=dp), INTENT(INOUT) :: A(LS_LU_NONZERO)
#endif
!~~~> Output variables
   INTEGER, INTENT(OUT) :: Pivot(N), ISING

#ifdef FULL_ALGEBRA    
   CALL  DGETRF( N, N, A, N, Pivot, ISING )
#else   
   CALL KppDecomp ( A, ISING,LS_NSEL,LS_LU_NONZERO,LS_LU_CROW,LS_LU_DIAG,LS_LU_ICOL)   
   Pivot(1) = 1
#endif
   ISTATUS(Ndec) = ISTATUS(Ndec) + 1

  END SUBROUTINE ros_Decomp


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Solve( A, Pivot, b)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the forward/backward substitution (using pre-computed LU decomposition)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Input variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) :: A(N,N)
   INTEGER :: ISING
#else   
   REAL(kind=dp), INTENT(IN) :: A(LS_LU_NONZERO)
#endif
   INTEGER, INTENT(IN) :: Pivot(N)
!~~~> InOut variables
   REAL(kind=dp), INTENT(INOUT) :: b(N)

#ifdef FULL_ALGEBRA    
   CALL  DGETRS( 'N', N , 1, A, N, Pivot, b, N, ISING )
   IF ( Info < 0 ) THEN
      PRINT*,"Error in DGETRS. ISING=",ISING
   END IF  
#else  
   	CALL KppSolve_2(LS_LU_NONZERO,LS_NSEL, A, b,LS_type )
#endif

   ISTATUS(Nsol) = ISTATUS(Nsol) + 1

  END SUBROUTINE ros_Solve



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   DOUBLE PRECISION g

    g = 1.0_dp + 1.0_dp/SQRT(2.0_dp)
    rosMethod = RS2
!~~~> Name of the method
    ros_Name = 'ROS-2'
!~~~> Number of stages
    ros_S = 2

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = (1.0_dp)/g
    ros_C(1) = (-2.0_dp)/g
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
!~~~> M_i = Coefficients for new step solution
    ros_M(1)= (3.0_dp)/(2.0_dp*g)
    ros_M(2)= (1.0_dp)/(2.0_dp*g)
! E_i = Coefficients for error estimator
    ros_E(1) = 1.0_dp/(2.0_dp*g)
    ros_E(2) = 1.0_dp/(2.0_dp*g)
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus one
    ros_ELO = 2.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.0_dp
    ros_Alpha(2) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = g
    ros_Gamma(2) =-g

 END SUBROUTINE Ros2


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   rosMethod = RS3
!~~~> Name of the method
   ros_Name = 'ROS-3'
!~~~> Number of stages
   ros_S = 3

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1)= 1.0_dp
   ros_A(2)= 1.0_dp
   ros_A(3)= 0.0_dp

   ros_C(1) = -0.10156171083877702091975600115545E+01_dp
   ros_C(2) =  0.40759956452537699824805835358067E+01_dp
   ros_C(3) =  0.92076794298330791242156818474003E+01_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) =  0.1E+01_dp
   ros_M(2) =  0.61697947043828245592553615689730E+01_dp
   ros_M(3) = -0.42772256543218573326238373806514_dp
! E_i = Coefficients for error estimator
   ros_E(1) =  0.5_dp
   ros_E(2) = -0.29079558716805469821718236208017E+01_dp
   ros_E(3) =  0.22354069897811569627360909276199_dp
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1)= 0.0_dp
   ros_Alpha(2)= 0.43586652150845899941601945119356_dp
   ros_Alpha(3)= 0.43586652150845899941601945119356_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1)= 0.43586652150845899941601945119356_dp
   ros_Gamma(2)= 0.24291996454816804366592249683314_dp
   ros_Gamma(3)= 0.21851380027664058511513169485832E+01_dp

  END SUBROUTINE Ros3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
!     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1990)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RS4
!~~~> Name of the method
   ros_Name = 'ROS-4'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.2000000000000000E+01_dp
   ros_A(2) = 0.1867943637803922E+01_dp
   ros_A(3) = 0.2344449711399156_dp
   ros_A(4) = ros_A(2)
   ros_A(5) = ros_A(3)
   ros_A(6) = 0.0_dp

   ros_C(1) =-0.7137615036412310E+01_dp
   ros_C(2) = 0.2580708087951457E+01_dp
   ros_C(3) = 0.6515950076447975_dp
   ros_C(4) =-0.2137148994382534E+01_dp
   ros_C(5) =-0.3214669691237626_dp
   ros_C(6) =-0.6949742501781779_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .TRUE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 0.2255570073418735E+01_dp
   ros_M(2) = 0.2870493262186792_dp
   ros_M(3) = 0.4353179431840180_dp
   ros_M(4) = 0.1093502252409163E+01_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) =-0.2815431932141155_dp
   ros_E(2) =-0.7276199124938920E-01_dp
   ros_E(3) =-0.1082196201495311_dp
   ros_E(4) =-0.1093502252409163E+01_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 4.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.1145640000000000E+01_dp
   ros_Alpha(3) = 0.6552168638155900_dp
   ros_Alpha(4) = ros_Alpha(3)
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5728200000000000_dp
   ros_Gamma(2) =-0.1769193891319233E+01_dp
   ros_Gamma(3) = 0.7592633437920482_dp
   ros_Gamma(4) =-0.1049021087100450_dp

  END SUBROUTINE Ros4

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- A STIFFLY-STABLE METHOD, 4 stages, order 3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RD3
!~~~> Name of the method
   ros_Name = 'RODAS-3'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.0_dp
   ros_A(2) = 2.0_dp
   ros_A(3) = 0.0_dp
   ros_A(4) = 2.0_dp
   ros_A(5) = 0.0_dp
   ros_A(6) = 1.0_dp

   ros_C(1) = 4.0_dp
   ros_C(2) = 1.0_dp
   ros_C(3) =-1.0_dp
   ros_C(4) = 1.0_dp
   ros_C(5) =-1.0_dp
   ros_C(6) =-(8.0_dp/3.0_dp)

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .FALSE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .TRUE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 2.0_dp
   ros_M(2) = 0.0_dp
   ros_M(3) = 1.0_dp
   ros_M(4) = 1.0_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) = 0.0_dp
   ros_E(2) = 0.0_dp
   ros_E(3) = 0.0_dp
   ros_E(4) = 1.0_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.0_dp
   ros_Alpha(3) = 1.0_dp
   ros_Alpha(4) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5_dp
   ros_Gamma(2) = 1.5_dp
   ros_Gamma(3) = 0.0_dp
   ros_Gamma(4) = 0.0_dp

  END SUBROUTINE Rodas3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1996)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RD4
!~~~> Name of the method
    ros_Name = 'RODAS-4'
!~~~> Number of stages
    ros_S = 6

!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.000_dp
    ros_Alpha(2) = 0.386_dp
    ros_Alpha(3) = 0.210_dp
    ros_Alpha(4) = 0.630_dp
    ros_Alpha(5) = 1.000_dp
    ros_Alpha(6) = 1.000_dp

!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = 0.2500000000000000_dp
    ros_Gamma(2) =-0.1043000000000000_dp
    ros_Gamma(3) = 0.1035000000000000_dp
    ros_Gamma(4) =-0.3620000000000023E-01_dp
    ros_Gamma(5) = 0.0_dp
    ros_Gamma(6) = 0.0_dp

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = 0.1544000000000000E+01_dp
    ros_A(2) = 0.9466785280815826_dp
    ros_A(3) = 0.2557011698983284_dp
    ros_A(4) = 0.3314825187068521E+01_dp
    ros_A(5) = 0.2896124015972201E+01_dp
    ros_A(6) = 0.9986419139977817_dp
    ros_A(7) = 0.1221224509226641E+01_dp
    ros_A(8) = 0.6019134481288629E+01_dp
    ros_A(9) = 0.1253708332932087E+02_dp
    ros_A(10) =-0.6878860361058950_dp
    ros_A(11) = ros_A(7)
    ros_A(12) = ros_A(8)
    ros_A(13) = ros_A(9)
    ros_A(14) = ros_A(10)
    ros_A(15) = 1.0_dp

    ros_C(1) =-0.5668800000000000E+01_dp
    ros_C(2) =-0.2430093356833875E+01_dp
    ros_C(3) =-0.2063599157091915_dp
    ros_C(4) =-0.1073529058151375_dp
    ros_C(5) =-0.9594562251023355E+01_dp
    ros_C(6) =-0.2047028614809616E+02_dp
    ros_C(7) = 0.7496443313967647E+01_dp
    ros_C(8) =-0.1024680431464352E+02_dp
    ros_C(9) =-0.3399990352819905E+02_dp
    ros_C(10) = 0.1170890893206160E+02_dp
    ros_C(11) = 0.8083246795921522E+01_dp
    ros_C(12) =-0.7981132988064893E+01_dp
    ros_C(13) =-0.3152159432874371E+02_dp
    ros_C(14) = 0.1631930543123136E+02_dp
    ros_C(15) =-0.6058818238834054E+01_dp

!~~~> M_i = Coefficients for new step solution
    ros_M(1) = ros_A(7)
    ros_M(2) = ros_A(8)
    ros_M(3) = ros_A(9)
    ros_M(4) = ros_A(10)
    ros_M(5) = 1.0_dp
    ros_M(6) = 1.0_dp

!~~~> E_i  = Coefficients for error estimator
    ros_E(1) = 0.0_dp
    ros_E(2) = 0.0_dp
    ros_E(3) = 0.0_dp
    ros_E(4) = 0.0_dp
    ros_E(5) = 0.0_dp
    ros_E(6) = 1.0_dp

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.
    ros_NewF(5) = .TRUE.
    ros_NewF(6) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 4.0_dp

  END SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! STIFFLY-STABLE W METHOD OF ORDER 3, WITH 4 STAGES
!
! J. RANG and L. ANGERMANN
! NEW ROSENBROCK W-METHODS OF ORDER 3
! FOR PARTIAL DIFFERENTIAL ALGEBRAIC
!        EQUATIONS OF INDEX 1                                             
! BIT Numerical Mathematics (2005) 45: 761-787
!  DOI: 10.1007/s10543-005-0035-y
! Table 4.1-4.2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RG3
!~~~> Name of the method
    ros_Name = 'RANG-3'
!~~~> Number of stages
    ros_S = 4

    ros_A(1) = 5.09052051067020d+00;
    ros_A(2) = 5.09052051067020d+00;
    ros_A(3) = 0.0d0;
    ros_A(4) = 4.97628111010787d+00;
    ros_A(5) = 2.77268164715849d-02;
    ros_A(6) = 2.29428036027904d-01;

    ros_C(1) = -1.16790812312283d+01;
    ros_C(2) = -1.64057326467367d+01;
    ros_C(3) = -2.77268164715850d-01;
    ros_C(4) = -8.38103960500476d+00;
    ros_C(5) = -8.48328409199343d-01;
    ros_C(6) =  2.87009860433106d-01;

    ros_M(1) =  5.22582761233094d+00;
    ros_M(2) = -5.56971148154165d-01;
    ros_M(3) =  3.57979469353645d-01;
    ros_M(4) =  1.72337398521064d+00;

    ros_E(1) = -5.16845212784040d+00;
    ros_E(2) = -1.26351942603842d+00;
    ros_E(3) = -1.11022302462516d-16;
    ros_E(4) =  2.22044604925031d-16;

    ros_Alpha(1) = 0.0d00;
    ros_Alpha(2) = 2.21878746765329d+00;
    ros_Alpha(3) = 2.21878746765329d+00;
    ros_Alpha(4) = 1.55392337535788d+00;

    ros_Gamma(1) =  4.35866521508459d-01;
    ros_Gamma(2) = -1.78292094614483d+00;
    ros_Gamma(3) = -2.46541900496934d+00;
    ros_Gamma(4) = -8.05529997906370d-01;

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 3.0_dp

  END SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END SUBROUTINE Rosenbrock_2


SUBROUTINE FunTemplate_2( T, Y, Ydot, LS_NSEL, LS_type )
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Function
 INTEGER,INTENT(IN)::LS_NSEL, LS_type
   REAL(kind=dp) :: T, Y(LS_NSEL)
   REAL(kind=dp) :: Ydot(LS_NSEL)
   REAL(kind=dp) :: Told
   Told = TIME
   TIME = T
   CALL Fun_2( Y, FIX, RCONST, Ydot, LS_NSEL)
   TIME = Told
END SUBROUTINE FunTemplate_2

SUBROUTINE JacTemplate_2( T, Y, Jcb, LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Jacobian
 USE gckpp_LinearAlgebra
 
 INTEGER,INTENT(IN)::LS_NSEL, LS_LU_NONZERO, LS_type
 INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO), LS_LU_ICOL(LS_LU_NONZERO)
    REAL(kind=dp) :: T, Y(LS_NSEL)
#ifdef FULL_ALGEBRA    
    REAL(kind=dp) :: JV(LS_LU_NONZERO), Jcb(LS_NSEL,LS_NSEL)
#else
    REAL(kind=dp) :: Jcb(LS_LU_NONZERO)
#endif   
    REAL(kind=dp) :: Told
#ifdef FULL_ALGEBRA    
    INTEGER :: i, j
#endif   

    Told = TIME
    TIME = T
#ifdef FULL_ALGEBRA        	
	CALL Jac_SP_2(Y, FIX, RCONST, JV,LS_NSEL,LS_LU_NONZERO)	
    DO j=1,LS_NSEL
      DO i=1,LS_NSEL
         Jcb(i,j) = 0.0_dp
      END DO
    END DO
    DO i=1,LS_LU_NONZERO
       Jcb(LS_LU_IROW(i),LS_LU_ICOL(i)) = JV(i)
    END DO
#else   
       CALL Jac_SP_2( Y, FIX, RCONST, Jcb,LS_NSEL,LS_LU_NONZERO)	 
#endif   
    TIME = Told
END SUBROUTINE JacTemplate_2



SUBROUTINE INTEGRATE_3( TIN, TOUT, LS_type,LS_NSEL, LS_NDEL,&
  Prate, Lrate, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )

   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: TIN  ! Start Time
   REAL(kind=dp), INTENT(IN) :: TOUT ! End Time
   
   ! Optional input parameters and statistics
   INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   REAL(kind=dp), INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   REAL(kind=dp), INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U
   INTEGER,		  INTENT(IN) :: LS_type,LS_NSEL,LS_NDEL
   REAL(kind=dp), INTENT(IN) :: Prate(NVAR),Lrate(NVAR)
   
   REAL(kind=dp) :: RCNTRL(20), RSTATUS(20),deltaT
   INTEGER       :: ICNTRL(20), ISTATUS(20), IERR
   REAL(kind=dp) :: VAR_selected(LS_NSEL),VAR_deleted(LS_NDEL),LS_P(LS_NDEL),LS_L(LS_NDEL)
   
   INTEGER, SAVE :: Ntotal = 0

   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.0_dp
   ISTATUS(:) = 0
   RSTATUS(:) = 0.0_dp
   deltaT =  TOUT-TIN
   !~~~> fine-tune the integrator:
   ICNTRL(1) = 0	! 0 - non-autonomous, 1 - autonomous
   ICNTRL(2) = 0	! 0 - vector tolerances, 1 - scalars

   ! If optional parameters are given, and if they are >0, 
   ! then they overwrite default settings. 
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF
   
        VAR_selected=VAR(select_ind_3)
	    VAR_deleted=VAR(delete_ind_3)
		LS_P=Prate(delete_ind_3)
		LS_L=Lrate(delete_ind_3)
        CALL Rosenbrock_3(LU_NSEL_3,VAR_selected,TIN,TOUT,ATOL,RTOL,&
          RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		  LU_NONZERO_3,LU_NSEL_3,LU_CROW_3,LU_DIAG_3,LU_IROW_3,LU_ICOL_3, LS_type)		
		WHERE(LS_L<=(0.01/deltaT))
			VAR_deleted=VAR_deleted+deltaT*(LS_P-LS_L*VAR_deleted)
		ELSEWHERE
			VAR_deleted=LS_P/LS_L+(VAR_deleted-LS_P/LS_L)*EXP(-LS_L*deltaT)
		END WHERE

	    VAR(select_ind_3)=VAR_selected
	    VAR(delete_ind_3)=VAR_deleted
   
   !~~~> Debug option: show no of steps
   ! Ntotal = Ntotal + ISTATUS(Nstp)
   ! PRINT*,'NSTEPS=',ISTATUS(Nstp),' (',Ntotal,')','  O3=', VAR(ind_O3)

   STEPMIN = RSTATUS(Nhexit)
   ! if optional parameters are given for output they 
   ! are updated with the return information
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE INTEGRATE_3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE Rosenbrock_3(N,Y,Tstart,Tend, &
           AbsTol,RelTol,              &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		   LS_LU_NONZERO,LS_NSEL,LS_LU_CROW,LS_LU_DIAG,LS_LU_IROW,LS_LU_ICOL, LS_type)
  USE gckpp_Parameters
  USE gckpp_LinearAlgebra
  IMPLICIT NONE

   INTEGER,       INTENT(IN)    :: N,LS_LU_NONZERO,LS_NSEL,LS_type   
   INTEGER,INTENT(IN)::LS_LU_CROW(LS_NSEL+1),LS_LU_DIAG(LS_NSEL+1)
   INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO),LS_LU_ICOL(LS_LU_NONZERO)   
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
   REAL(kind=dp), INTENT(IN)    :: Tstart,Tend
   REAL(kind=dp), INTENT(IN)    :: AbsTol(N),RelTol(N)
   INTEGER,       INTENT(IN)    :: ICNTRL(20)
   REAL(kind=dp), INTENT(IN)    :: RCNTRL(20)
   INTEGER,       INTENT(INOUT) :: ISTATUS(20)
   REAL(kind=dp), INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT)   :: IERR
!~~~>  Parameters of the Rosenbrock method, up to 6 stages
   INTEGER ::  ros_S, rosMethod
   INTEGER, PARAMETER :: RS2=1, RS3=2, RS4=3, RD3=4, RD4=5, RG3=6
   REAL(kind=dp) :: ros_A(15), ros_C(15), ros_M(6), ros_E(6), &
                    ros_Alpha(6), ros_Gamma(6), ros_ELO
   LOGICAL :: ros_NewF(6)
   CHARACTER(LEN=12) :: ros_Name
!~~~>  Local variables
   REAL(kind=dp) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   REAL(kind=dp) :: Hmin, Hmax, Hstart
   REAL(kind=dp) :: Texit
   INTEGER       :: i, UplimTol, Max_no_steps
   LOGICAL       :: Autonomous, VectorTol
!~~~>   Parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp

!~~~>  Initialize statistics
   ISTATUS(1:8) = 0
   RSTATUS(1:3) = ZERO

!~~~>  Autonomous or time dependent ODE. Default is time dependent.
   Autonomous = .NOT.(ICNTRL(1) == 0)

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:N) and RelTol(1:N)
   IF (ICNTRL(2) == 0) THEN
      VectorTol = .TRUE.
      UplimTol  = N
   ELSE
      VectorTol = .FALSE.
      UplimTol  = 1
   END IF

!~~~>   Initialize the particular Rosenbrock method selected
   SELECT CASE (ICNTRL(3))
     CASE (1)
       CALL Ros2
     CASE (2)
       CALL Ros3
     CASE (3)
       CALL Ros4
     CASE (0,4)
       CALL Rodas3
     CASE (5)
       CALL Rodas4
     CASE (6)
       CALL Rang3
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(3)=',ICNTRL(3) 
       CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT

!~~~>   The maximum number of steps admitted
   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 200000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Unit roundoff (1+Roundoff>1)
   Roundoff = WLAMCH('E')

!~~~>  Lower bound on the step size: (positive value)
   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Upper bound on the step size: (positive value)
   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Starting step size: (positive value)
   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hold < FacMax
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Check if tolerances are reasonable
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO


!~~~>  CALL Rosenbrock method
   CALL ros_Integrator(Y, Tstart, Tend, Texit,   &
        AbsTol, RelTol,                          &
!  Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!  Error indicator
        IERR)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS !  SUBROUTINES internal to Rosenbrock
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 SUBROUTINE ros_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
   REAL(kind=dp), INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR
   
   IERR = Code
   PRINT * , &
     'Forced exit from Rosenbrock due to the following error:' 
     
   SELECT CASE (Code)
    CASE (-1)
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)    
      PRINT * , '--> Selected Rosenbrock method not implemented'
    CASE (-3)    
      PRINT * , '--> Hmin/Hmax/Hstart must be positive'
    CASE (-4)    
      PRINT * , '--> FacMin/FacMax/FacRej must be positive'
    CASE (-5) 
      PRINT * , '--> Improper tolerance values'
    CASE (-6) 
      PRINT * , '--> No of steps exceeds maximum bound'
    CASE (-7) 
      PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
    CASE (-8)    
      PRINT * , '--> Matrix is repeatedly singular'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT
   
   PRINT *, "T=", T, "and H=", H
     
 END SUBROUTINE ros_ErrorMsg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_Integrator (Y, Tstart, Tend, T,  &
        AbsTol, RelTol,                          &
!~~~> Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic Rosenbrock method
!      defined by ros_S (no of stages)
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

!~~~> Input: the initial condition at Tstart; Output: the solution at T
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
!~~~> Input: integration interval
   REAL(kind=dp), INTENT(IN) :: Tstart,Tend
!~~~> Output: time at which the solution is returned (T=Tend if success)
   REAL(kind=dp), INTENT(OUT) ::  T
!~~~> Input: tolerances
   REAL(kind=dp), INTENT(IN) ::  AbsTol(N), RelTol(N)
!~~~> Input: integration parameters
   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   REAL(kind=dp), INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   REAL(kind=dp), INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables
   REAL(kind=dp) :: Ynew(N), Fcn0(N), Fcn(N)
   REAL(kind=dp) :: K(N*ros_S), dFdT(N)
#ifdef FULL_ALGEBRA    
   REAL(kind=dp) :: Jac0(N,N), Ghimj(N,N)
#else
   REAL(kind=dp) :: Jac0(LS_LU_NONZERO), Ghimj(LS_LU_NONZERO)
#endif
   REAL(kind=dp) :: H, Hnew, HC, HG, Fac, Tau
   REAL(kind=dp) :: Err, Yerr(N)
   INTEGER :: Pivot(N), Direction, ioffset, j, istage
   LOGICAL :: RejectLastH, RejectMoreH, Singular
!~~~>  Local parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp
!~~~>  Locally called functions
!    REAL(kind=dp) WLAMCH
!    EXTERNAL WLAMCH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~>  Initial preparations
   T = Tstart
   RSTATUS(Nhexit) = ZERO
   H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   IF (ABS(H) <= 10.0_dp*Roundoff) H = DeltaMin

   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF
   H = Direction*H

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.

!~~~> Time loop begins below

TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) )

   IF ( ISTATUS(Nstp) > Max_no_steps ) THEN  ! Too many steps
      CALL ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
      CALL ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF

!~~~>  Limit H if necessary to avoid going beyond Tend
   H = MIN(H,ABS(Tend-T))

!~~~>   Compute the function at current time
   CALL FunTemplate_3(T,Y,Fcn0, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1

!~~~>  Compute the function derivative with respect to T
   IF (.NOT.Autonomous) THEN
      CALL ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
   END IF

!~~~>   Compute the Jacobian at current time
   CALL JacTemplate_3(T,Y,Jac0,LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)   
   ISTATUS(Njac) = ISTATUS(Njac) + 1

!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO

   CALL ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular)
   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       CALL ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S

      ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
       ioffset = N*(istage-1)

      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         !slim: CALL WCOPY(N,Fcn0,1,Fcn,1)
	 Fcn(1:N) = Fcn0(1:N)
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
         !slim: CALL WCOPY(N,Y,1,Ynew,1)
	 Ynew(1:N) = Y(1:N)
         DO j = 1, istage-1
           CALL WAXPY(N,ros_A((istage-1)*(istage-2)/2+j), &
            K(N*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL FunTemplate_3(Tau,Ynew,Fcn, LS_NSEL, LS_type)
         ISTATUS(Nfun) = ISTATUS(Nfun) + 1
       END IF ! if istage == 1 elseif ros_NewF(istage)
       !slim: CALL WCOPY(N,Fcn,1,K(ioffset+1),1)
       K(ioffset+1:ioffset+N) = Fcn(1:N)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL WAXPY(N,HC,K(N*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL WAXPY(N,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL ros_Solve(Ghimj, Pivot, K(ioffset+1))

   END DO Stage


!~~~>  Compute the new solution
   !slim: CALL WCOPY(N,Y,1,Ynew,1)
   Ynew(1:N) = Y(1:N)
   DO j=1,ros_S
         CALL WAXPY(N,ros_M(j),K(N*(j-1)+1),1,Ynew,1)
   END DO

!~~~>  Compute the error estimation
   !slim: CALL WSCAL(N,ZERO,Yerr,1)
   Yerr(1:N) = ZERO
   DO j=1,ros_S
        CALL WAXPY(N,ros_E(j),K(N*(j-1)+1),1,Yerr,1)
   END DO
   Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac

!~~~>  Check the error magnitude and adjust step size
   ISTATUS(Nstp) = ISTATUS(Nstp) + 1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
      ISTATUS(Nacc) = ISTATUS(Nacc) + 1
      !slim: CALL WCOPY(N,Ynew,1,Y,1)
      Y(1:N) = Ynew(1:N)
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H)
      END IF
      RSTATUS(Nhexit) = H
      RSTATUS(Nhnew)  = Hnew
      RSTATUS(Ntexit) = T
      RejectLastH = .FALSE.
      RejectMoreH = .FALSE.
      H = Hnew
      EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
   ELSE           !~~~> Reject step
      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (ISTATUS(Nacc) >= 1)  ISTATUS(Nrej) = ISTATUS(Nrej) + 1
   END IF ! Err <= 1

   END DO UntilAccepted

   
   END DO TimeLoop

!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful

  END SUBROUTINE ros_Integrator


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  REAL(kind=dp) FUNCTION ros_ErrorNorm ( Y, Ynew, Yerr, &
                               AbsTol, RelTol, VectorTol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

! Input arguments
   REAL(kind=dp), INTENT(IN) :: Y(N), Ynew(N), &
          Yerr(N), AbsTol(N), RelTol(N)
   LOGICAL, INTENT(IN) ::  VectorTol
! Local variables
   REAL(kind=dp) :: Err, Scale, Ymax
   INTEGER  :: i
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp

   Err = ZERO
   DO i=1,N
     Ymax = MAX(ABS(Y(i)),ABS(Ynew(i)))
     IF (VectorTol) THEN
       Scale = AbsTol(i)+RelTol(i)*Ymax
     ELSE
       Scale = AbsTol(1)+RelTol(1)*Ymax
     END IF
     Err = Err+(Yerr(i)/Scale)**2
   END DO
   Err  = SQRT(Err/N)

   ros_ErrorNorm = MAX(Err,1.0d-10)

  END FUNCTION ros_ErrorNorm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> The time partial derivative of the function by finite differences
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

!~~~> Input arguments
   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(N), Fcn0(N)
!~~~> Output arguments
   REAL(kind=dp), INTENT(OUT) :: dFdT(N)
!~~~> Local variables
   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL FunTemplate_3(T+Delta,Y,dFdT, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1
   CALL WAXPY(N,(-ONE),Fcn0,1,dFdT,1)
   CALL WSCAL(N,(ONE/Delta),dFdT,1)

  END SUBROUTINE ros_FunTimeDerivative


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_PrepareMatrix ( H, Direction, gam, &
             Jac0, Ghimj, Pivot, Singular )
! --- --- --- --- --- --- --- --- --- --- --- --- ---
!  Prepares the LHS matrix for stage calculations
!  1.  Construct Ghimj = 1/(H*ham) - Jac0
!      "(Gamma H) Inverse Minus Jacobian"
!  2.  Repeat LU decomposition of Ghimj until successful.
!       -half the step size if LU decomposition fails and retry
!       -exit after 5 consecutive fails
! --- --- --- --- --- --- --- --- --- --- --- --- ---
   IMPLICIT NONE

!~~~> Input arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) ::  Jac0(N,N)
#else
   REAL(kind=dp), INTENT(IN) ::  Jac0(LS_LU_NONZERO)
#endif   
   REAL(kind=dp), INTENT(IN) ::  gam
   INTEGER, INTENT(IN) ::  Direction
!~~~> Output arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(OUT) :: Ghimj(N,N)
#else
   REAL(kind=dp), INTENT(OUT) :: Ghimj(LS_LU_NONZERO)
#endif   
   LOGICAL, INTENT(OUT) ::  Singular
   INTEGER, INTENT(OUT) ::  Pivot(N)
!~~~> Inout arguments
   REAL(kind=dp), INTENT(INOUT) :: H   ! step size is decreased when LU fails
!~~~> Local variables
   INTEGER  :: i, ISING, Nconsecutive
   REAL(kind=dp) :: ghinv
   REAL(kind=dp), PARAMETER :: ONE  = 1.0_dp, HALF = 0.5_dp

   Nconsecutive = 0
   Singular = .TRUE.

   DO WHILE (Singular)

!~~~>    Construct Ghimj = 1/(H*gam) - Jac0
#ifdef FULL_ALGEBRA    
     !slim: CALL WCOPY(N*N,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(N*N,(-ONE),Ghimj,1)
     Ghimj = -Jac0
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(i,i) = Ghimj(i,i)+ghinv
     END DO
#else
     !slim: CALL WCOPY(LS_LU_NONZERO,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(LS_LU_NONZERO,(-ONE),Ghimj,1)
     Ghimj(1:LS_LU_NONZERO) = -Jac0(1:LS_LU_NONZERO)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(LS_LU_DIAG(i)) = Ghimj(LS_LU_DIAG(i))+ghinv
     END DO
#endif   
!~~~>    Compute LU decomposition
     CALL ros_Decomp( Ghimj, Pivot, ISING )
     IF (ISING == 0) THEN
!~~~>    If successful done
        Singular = .FALSE.
     ELSE ! ISING .ne. 0
!~~~>    If unsuccessful half the step size; if 5 consecutive fails then return
        ISTATUS(Nsng) = ISTATUS(Nsng) + 1
        Nconsecutive = Nconsecutive+1
        Singular = .TRUE.
        PRINT*,'Warning: LU Decomposition returned ISING = ',ISING
        IF (Nconsecutive <= 5) THEN ! Less than 5 consecutive failed decompositions
           H = H*HALF
        ELSE  ! More than 5 consecutive failed decompositions
           RETURN
        END IF  ! Nconsecutive
      END IF    ! ISING

   END DO ! WHILE Singular

  END SUBROUTINE ros_PrepareMatrix


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Decomp( A, Pivot, ISING )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the LU decomposition
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Inout variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(INOUT) :: A(N,N)
#else   
   REAL(kind=dp), INTENT(INOUT) :: A(LS_LU_NONZERO)
#endif
!~~~> Output variables
   INTEGER, INTENT(OUT) :: Pivot(N), ISING

#ifdef FULL_ALGEBRA    
   CALL  DGETRF( N, N, A, N, Pivot, ISING )
#else   
   CALL KppDecomp ( A, ISING,LS_NSEL,LS_LU_NONZERO,LS_LU_CROW,LS_LU_DIAG,LS_LU_ICOL)   
   Pivot(1) = 1
#endif
   ISTATUS(Ndec) = ISTATUS(Ndec) + 1

  END SUBROUTINE ros_Decomp


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Solve( A, Pivot, b)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the forward/backward substitution (using pre-computed LU decomposition)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Input variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) :: A(N,N)
   INTEGER :: ISING
#else   
   REAL(kind=dp), INTENT(IN) :: A(LS_LU_NONZERO)
#endif
   INTEGER, INTENT(IN) :: Pivot(N)
!~~~> InOut variables
   REAL(kind=dp), INTENT(INOUT) :: b(N)

#ifdef FULL_ALGEBRA    
   CALL  DGETRS( 'N', N , 1, A, N, Pivot, b, N, ISING )
   IF ( Info < 0 ) THEN
      PRINT*,"Error in DGETRS. ISING=",ISING
   END IF  
#else  
   	CALL KppSolve_3(LS_LU_NONZERO,LS_NSEL, A, b,LS_type )
#endif

   ISTATUS(Nsol) = ISTATUS(Nsol) + 1

  END SUBROUTINE ros_Solve



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   DOUBLE PRECISION g

    g = 1.0_dp + 1.0_dp/SQRT(2.0_dp)
    rosMethod = RS2
!~~~> Name of the method
    ros_Name = 'ROS-2'
!~~~> Number of stages
    ros_S = 2

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = (1.0_dp)/g
    ros_C(1) = (-2.0_dp)/g
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
!~~~> M_i = Coefficients for new step solution
    ros_M(1)= (3.0_dp)/(2.0_dp*g)
    ros_M(2)= (1.0_dp)/(2.0_dp*g)
! E_i = Coefficients for error estimator
    ros_E(1) = 1.0_dp/(2.0_dp*g)
    ros_E(2) = 1.0_dp/(2.0_dp*g)
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus one
    ros_ELO = 2.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.0_dp
    ros_Alpha(2) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = g
    ros_Gamma(2) =-g

 END SUBROUTINE Ros2


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   rosMethod = RS3
!~~~> Name of the method
   ros_Name = 'ROS-3'
!~~~> Number of stages
   ros_S = 3

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1)= 1.0_dp
   ros_A(2)= 1.0_dp
   ros_A(3)= 0.0_dp

   ros_C(1) = -0.10156171083877702091975600115545E+01_dp
   ros_C(2) =  0.40759956452537699824805835358067E+01_dp
   ros_C(3) =  0.92076794298330791242156818474003E+01_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) =  0.1E+01_dp
   ros_M(2) =  0.61697947043828245592553615689730E+01_dp
   ros_M(3) = -0.42772256543218573326238373806514_dp
! E_i = Coefficients for error estimator
   ros_E(1) =  0.5_dp
   ros_E(2) = -0.29079558716805469821718236208017E+01_dp
   ros_E(3) =  0.22354069897811569627360909276199_dp
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1)= 0.0_dp
   ros_Alpha(2)= 0.43586652150845899941601945119356_dp
   ros_Alpha(3)= 0.43586652150845899941601945119356_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1)= 0.43586652150845899941601945119356_dp
   ros_Gamma(2)= 0.24291996454816804366592249683314_dp
   ros_Gamma(3)= 0.21851380027664058511513169485832E+01_dp

  END SUBROUTINE Ros3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
!     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1990)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RS4
!~~~> Name of the method
   ros_Name = 'ROS-4'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.2000000000000000E+01_dp
   ros_A(2) = 0.1867943637803922E+01_dp
   ros_A(3) = 0.2344449711399156_dp
   ros_A(4) = ros_A(2)
   ros_A(5) = ros_A(3)
   ros_A(6) = 0.0_dp

   ros_C(1) =-0.7137615036412310E+01_dp
   ros_C(2) = 0.2580708087951457E+01_dp
   ros_C(3) = 0.6515950076447975_dp
   ros_C(4) =-0.2137148994382534E+01_dp
   ros_C(5) =-0.3214669691237626_dp
   ros_C(6) =-0.6949742501781779_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .TRUE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 0.2255570073418735E+01_dp
   ros_M(2) = 0.2870493262186792_dp
   ros_M(3) = 0.4353179431840180_dp
   ros_M(4) = 0.1093502252409163E+01_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) =-0.2815431932141155_dp
   ros_E(2) =-0.7276199124938920E-01_dp
   ros_E(3) =-0.1082196201495311_dp
   ros_E(4) =-0.1093502252409163E+01_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 4.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.1145640000000000E+01_dp
   ros_Alpha(3) = 0.6552168638155900_dp
   ros_Alpha(4) = ros_Alpha(3)
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5728200000000000_dp
   ros_Gamma(2) =-0.1769193891319233E+01_dp
   ros_Gamma(3) = 0.7592633437920482_dp
   ros_Gamma(4) =-0.1049021087100450_dp

  END SUBROUTINE Ros4

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- A STIFFLY-STABLE METHOD, 4 stages, order 3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RD3
!~~~> Name of the method
   ros_Name = 'RODAS-3'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.0_dp
   ros_A(2) = 2.0_dp
   ros_A(3) = 0.0_dp
   ros_A(4) = 2.0_dp
   ros_A(5) = 0.0_dp
   ros_A(6) = 1.0_dp

   ros_C(1) = 4.0_dp
   ros_C(2) = 1.0_dp
   ros_C(3) =-1.0_dp
   ros_C(4) = 1.0_dp
   ros_C(5) =-1.0_dp
   ros_C(6) =-(8.0_dp/3.0_dp)

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .FALSE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .TRUE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 2.0_dp
   ros_M(2) = 0.0_dp
   ros_M(3) = 1.0_dp
   ros_M(4) = 1.0_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) = 0.0_dp
   ros_E(2) = 0.0_dp
   ros_E(3) = 0.0_dp
   ros_E(4) = 1.0_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.0_dp
   ros_Alpha(3) = 1.0_dp
   ros_Alpha(4) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5_dp
   ros_Gamma(2) = 1.5_dp
   ros_Gamma(3) = 0.0_dp
   ros_Gamma(4) = 0.0_dp

  END SUBROUTINE Rodas3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1996)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RD4
!~~~> Name of the method
    ros_Name = 'RODAS-4'
!~~~> Number of stages
    ros_S = 6

!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.000_dp
    ros_Alpha(2) = 0.386_dp
    ros_Alpha(3) = 0.210_dp
    ros_Alpha(4) = 0.630_dp
    ros_Alpha(5) = 1.000_dp
    ros_Alpha(6) = 1.000_dp

!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = 0.2500000000000000_dp
    ros_Gamma(2) =-0.1043000000000000_dp
    ros_Gamma(3) = 0.1035000000000000_dp
    ros_Gamma(4) =-0.3620000000000023E-01_dp
    ros_Gamma(5) = 0.0_dp
    ros_Gamma(6) = 0.0_dp

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = 0.1544000000000000E+01_dp
    ros_A(2) = 0.9466785280815826_dp
    ros_A(3) = 0.2557011698983284_dp
    ros_A(4) = 0.3314825187068521E+01_dp
    ros_A(5) = 0.2896124015972201E+01_dp
    ros_A(6) = 0.9986419139977817_dp
    ros_A(7) = 0.1221224509226641E+01_dp
    ros_A(8) = 0.6019134481288629E+01_dp
    ros_A(9) = 0.1253708332932087E+02_dp
    ros_A(10) =-0.6878860361058950_dp
    ros_A(11) = ros_A(7)
    ros_A(12) = ros_A(8)
    ros_A(13) = ros_A(9)
    ros_A(14) = ros_A(10)
    ros_A(15) = 1.0_dp

    ros_C(1) =-0.5668800000000000E+01_dp
    ros_C(2) =-0.2430093356833875E+01_dp
    ros_C(3) =-0.2063599157091915_dp
    ros_C(4) =-0.1073529058151375_dp
    ros_C(5) =-0.9594562251023355E+01_dp
    ros_C(6) =-0.2047028614809616E+02_dp
    ros_C(7) = 0.7496443313967647E+01_dp
    ros_C(8) =-0.1024680431464352E+02_dp
    ros_C(9) =-0.3399990352819905E+02_dp
    ros_C(10) = 0.1170890893206160E+02_dp
    ros_C(11) = 0.8083246795921522E+01_dp
    ros_C(12) =-0.7981132988064893E+01_dp
    ros_C(13) =-0.3152159432874371E+02_dp
    ros_C(14) = 0.1631930543123136E+02_dp
    ros_C(15) =-0.6058818238834054E+01_dp

!~~~> M_i = Coefficients for new step solution
    ros_M(1) = ros_A(7)
    ros_M(2) = ros_A(8)
    ros_M(3) = ros_A(9)
    ros_M(4) = ros_A(10)
    ros_M(5) = 1.0_dp
    ros_M(6) = 1.0_dp

!~~~> E_i  = Coefficients for error estimator
    ros_E(1) = 0.0_dp
    ros_E(2) = 0.0_dp
    ros_E(3) = 0.0_dp
    ros_E(4) = 0.0_dp
    ros_E(5) = 0.0_dp
    ros_E(6) = 1.0_dp

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.
    ros_NewF(5) = .TRUE.
    ros_NewF(6) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 4.0_dp

  END SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! STIFFLY-STABLE W METHOD OF ORDER 3, WITH 4 STAGES
!
! J. RANG and L. ANGERMANN
! NEW ROSENBROCK W-METHODS OF ORDER 3
! FOR PARTIAL DIFFERENTIAL ALGEBRAIC
!        EQUATIONS OF INDEX 1                                             
! BIT Numerical Mathematics (2005) 45: 761-787
!  DOI: 10.1007/s10543-005-0035-y
! Table 4.1-4.2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RG3
!~~~> Name of the method
    ros_Name = 'RANG-3'
!~~~> Number of stages
    ros_S = 4

    ros_A(1) = 5.09052051067020d+00;
    ros_A(2) = 5.09052051067020d+00;
    ros_A(3) = 0.0d0;
    ros_A(4) = 4.97628111010787d+00;
    ros_A(5) = 2.77268164715849d-02;
    ros_A(6) = 2.29428036027904d-01;

    ros_C(1) = -1.16790812312283d+01;
    ros_C(2) = -1.64057326467367d+01;
    ros_C(3) = -2.77268164715850d-01;
    ros_C(4) = -8.38103960500476d+00;
    ros_C(5) = -8.48328409199343d-01;
    ros_C(6) =  2.87009860433106d-01;

    ros_M(1) =  5.22582761233094d+00;
    ros_M(2) = -5.56971148154165d-01;
    ros_M(3) =  3.57979469353645d-01;
    ros_M(4) =  1.72337398521064d+00;

    ros_E(1) = -5.16845212784040d+00;
    ros_E(2) = -1.26351942603842d+00;
    ros_E(3) = -1.11022302462516d-16;
    ros_E(4) =  2.22044604925031d-16;

    ros_Alpha(1) = 0.0d00;
    ros_Alpha(2) = 2.21878746765329d+00;
    ros_Alpha(3) = 2.21878746765329d+00;
    ros_Alpha(4) = 1.55392337535788d+00;

    ros_Gamma(1) =  4.35866521508459d-01;
    ros_Gamma(2) = -1.78292094614483d+00;
    ros_Gamma(3) = -2.46541900496934d+00;
    ros_Gamma(4) = -8.05529997906370d-01;

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 3.0_dp

  END SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END SUBROUTINE Rosenbrock_3


SUBROUTINE FunTemplate_3( T, Y, Ydot, LS_NSEL, LS_type )
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Function
 INTEGER,INTENT(IN)::LS_NSEL, LS_type
   REAL(kind=dp) :: T, Y(LS_NSEL)
   REAL(kind=dp) :: Ydot(LS_NSEL)
   REAL(kind=dp) :: Told
   Told = TIME
   TIME = T
   CALL Fun_3( Y, FIX, RCONST, Ydot, LS_NSEL)
   TIME = Told
END SUBROUTINE FunTemplate_3

SUBROUTINE JacTemplate_3( T, Y, Jcb, LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Jacobian
 USE gckpp_LinearAlgebra
 
 INTEGER,INTENT(IN)::LS_NSEL, LS_LU_NONZERO, LS_type
 INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO), LS_LU_ICOL(LS_LU_NONZERO)
    REAL(kind=dp) :: T, Y(LS_NSEL)
#ifdef FULL_ALGEBRA    
    REAL(kind=dp) :: JV(LS_LU_NONZERO), Jcb(LS_NSEL,LS_NSEL)
#else
    REAL(kind=dp) :: Jcb(LS_LU_NONZERO)
#endif   
    REAL(kind=dp) :: Told
#ifdef FULL_ALGEBRA    
    INTEGER :: i, j
#endif   

    Told = TIME
    TIME = T
#ifdef FULL_ALGEBRA        	
	CALL Jac_SP_3(Y, FIX, RCONST, JV,LS_NSEL,LS_LU_NONZERO)	
    DO j=1,LS_NSEL
      DO i=1,LS_NSEL
         Jcb(i,j) = 0.0_dp
      END DO
    END DO
    DO i=1,LS_LU_NONZERO
       Jcb(LS_LU_IROW(i),LS_LU_ICOL(i)) = JV(i)
    END DO
#else   
       CALL Jac_SP_3( Y, FIX, RCONST, Jcb,LS_NSEL,LS_LU_NONZERO)	 
#endif   
    TIME = Told
END SUBROUTINE JacTemplate_3



SUBROUTINE INTEGRATE_4( TIN, TOUT, LS_type,LS_NSEL, LS_NDEL,&
  Prate, Lrate, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )

   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: TIN  ! Start Time
   REAL(kind=dp), INTENT(IN) :: TOUT ! End Time
   
   ! Optional input parameters and statistics
   INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   REAL(kind=dp), INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   REAL(kind=dp), INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U
   INTEGER,		  INTENT(IN) :: LS_type,LS_NSEL,LS_NDEL
   REAL(kind=dp), INTENT(IN) :: Prate(NVAR),Lrate(NVAR)
   
   REAL(kind=dp) :: RCNTRL(20), RSTATUS(20),deltaT
   INTEGER       :: ICNTRL(20), ISTATUS(20), IERR
   REAL(kind=dp) :: VAR_selected(LS_NSEL),VAR_deleted(LS_NDEL),LS_P(LS_NDEL),LS_L(LS_NDEL)
   
   INTEGER, SAVE :: Ntotal = 0

   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.0_dp
   ISTATUS(:) = 0
   RSTATUS(:) = 0.0_dp
   deltaT =  TOUT-TIN
   !~~~> fine-tune the integrator:
   ICNTRL(1) = 0	! 0 - non-autonomous, 1 - autonomous
   ICNTRL(2) = 0	! 0 - vector tolerances, 1 - scalars

   ! If optional parameters are given, and if they are >0, 
   ! then they overwrite default settings. 
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF
   
        VAR_selected=VAR(select_ind_4)
	    VAR_deleted=VAR(delete_ind_4)
		LS_P=Prate(delete_ind_4)
		LS_L=Lrate(delete_ind_4)
        CALL Rosenbrock_4(LU_NSEL_4,VAR_selected,TIN,TOUT,ATOL,RTOL,&
          RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		  LU_NONZERO_4,LU_NSEL_4,LU_CROW_4,LU_DIAG_4,LU_IROW_4,LU_ICOL_4, LS_type)		
		WHERE(LS_L<=(0.01/deltaT))
			VAR_deleted=VAR_deleted+deltaT*(LS_P-LS_L*VAR_deleted)
		ELSEWHERE
			VAR_deleted=LS_P/LS_L+(VAR_deleted-LS_P/LS_L)*EXP(-LS_L*deltaT)
		END WHERE

	    VAR(select_ind_4)=VAR_selected
	    VAR(delete_ind_4)=VAR_deleted
   
   !~~~> Debug option: show no of steps
   ! Ntotal = Ntotal + ISTATUS(Nstp)
   ! PRINT*,'NSTEPS=',ISTATUS(Nstp),' (',Ntotal,')','  O3=', VAR(ind_O3)

   STEPMIN = RSTATUS(Nhexit)
   ! if optional parameters are given for output they 
   ! are updated with the return information
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE INTEGRATE_4

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE Rosenbrock_4(N,Y,Tstart,Tend, &
           AbsTol,RelTol,              &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		   LS_LU_NONZERO,LS_NSEL,LS_LU_CROW,LS_LU_DIAG,LS_LU_IROW,LS_LU_ICOL, LS_type)
  USE gckpp_Parameters
  USE gckpp_LinearAlgebra
  IMPLICIT NONE

   INTEGER,       INTENT(IN)    :: N,LS_LU_NONZERO,LS_NSEL,LS_type   
   INTEGER,INTENT(IN)::LS_LU_CROW(LS_NSEL+1),LS_LU_DIAG(LS_NSEL+1)
   INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO),LS_LU_ICOL(LS_LU_NONZERO)   
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
   REAL(kind=dp), INTENT(IN)    :: Tstart,Tend
   REAL(kind=dp), INTENT(IN)    :: AbsTol(N),RelTol(N)
   INTEGER,       INTENT(IN)    :: ICNTRL(20)
   REAL(kind=dp), INTENT(IN)    :: RCNTRL(20)
   INTEGER,       INTENT(INOUT) :: ISTATUS(20)
   REAL(kind=dp), INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT)   :: IERR
!~~~>  Parameters of the Rosenbrock method, up to 6 stages
   INTEGER ::  ros_S, rosMethod
   INTEGER, PARAMETER :: RS2=1, RS3=2, RS4=3, RD3=4, RD4=5, RG3=6
   REAL(kind=dp) :: ros_A(15), ros_C(15), ros_M(6), ros_E(6), &
                    ros_Alpha(6), ros_Gamma(6), ros_ELO
   LOGICAL :: ros_NewF(6)
   CHARACTER(LEN=12) :: ros_Name
!~~~>  Local variables
   REAL(kind=dp) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   REAL(kind=dp) :: Hmin, Hmax, Hstart
   REAL(kind=dp) :: Texit
   INTEGER       :: i, UplimTol, Max_no_steps
   LOGICAL       :: Autonomous, VectorTol
!~~~>   Parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp

!~~~>  Initialize statistics
   ISTATUS(1:8) = 0
   RSTATUS(1:3) = ZERO

!~~~>  Autonomous or time dependent ODE. Default is time dependent.
   Autonomous = .NOT.(ICNTRL(1) == 0)

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:N) and RelTol(1:N)
   IF (ICNTRL(2) == 0) THEN
      VectorTol = .TRUE.
      UplimTol  = N
   ELSE
      VectorTol = .FALSE.
      UplimTol  = 1
   END IF

!~~~>   Initialize the particular Rosenbrock method selected
   SELECT CASE (ICNTRL(3))
     CASE (1)
       CALL Ros2
     CASE (2)
       CALL Ros3
     CASE (3)
       CALL Ros4
     CASE (0,4)
       CALL Rodas3
     CASE (5)
       CALL Rodas4
     CASE (6)
       CALL Rang3
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(3)=',ICNTRL(3) 
       CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT

!~~~>   The maximum number of steps admitted
   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 200000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Unit roundoff (1+Roundoff>1)
   Roundoff = WLAMCH('E')

!~~~>  Lower bound on the step size: (positive value)
   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Upper bound on the step size: (positive value)
   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Starting step size: (positive value)
   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hold < FacMax
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Check if tolerances are reasonable
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO


!~~~>  CALL Rosenbrock method
   CALL ros_Integrator(Y, Tstart, Tend, Texit,   &
        AbsTol, RelTol,                          &
!  Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!  Error indicator
        IERR)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS !  SUBROUTINES internal to Rosenbrock
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 SUBROUTINE ros_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
   REAL(kind=dp), INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR
   
   IERR = Code
   PRINT * , &
     'Forced exit from Rosenbrock due to the following error:' 
     
   SELECT CASE (Code)
    CASE (-1)
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)    
      PRINT * , '--> Selected Rosenbrock method not implemented'
    CASE (-3)    
      PRINT * , '--> Hmin/Hmax/Hstart must be positive'
    CASE (-4)    
      PRINT * , '--> FacMin/FacMax/FacRej must be positive'
    CASE (-5) 
      PRINT * , '--> Improper tolerance values'
    CASE (-6) 
      PRINT * , '--> No of steps exceeds maximum bound'
    CASE (-7) 
      PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
    CASE (-8)    
      PRINT * , '--> Matrix is repeatedly singular'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT
   
   PRINT *, "T=", T, "and H=", H
     
 END SUBROUTINE ros_ErrorMsg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_Integrator (Y, Tstart, Tend, T,  &
        AbsTol, RelTol,                          &
!~~~> Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic Rosenbrock method
!      defined by ros_S (no of stages)
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

!~~~> Input: the initial condition at Tstart; Output: the solution at T
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
!~~~> Input: integration interval
   REAL(kind=dp), INTENT(IN) :: Tstart,Tend
!~~~> Output: time at which the solution is returned (T=Tend if success)
   REAL(kind=dp), INTENT(OUT) ::  T
!~~~> Input: tolerances
   REAL(kind=dp), INTENT(IN) ::  AbsTol(N), RelTol(N)
!~~~> Input: integration parameters
   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   REAL(kind=dp), INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   REAL(kind=dp), INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables
   REAL(kind=dp) :: Ynew(N), Fcn0(N), Fcn(N)
   REAL(kind=dp) :: K(N*ros_S), dFdT(N)
#ifdef FULL_ALGEBRA    
   REAL(kind=dp) :: Jac0(N,N), Ghimj(N,N)
#else
   REAL(kind=dp) :: Jac0(LS_LU_NONZERO), Ghimj(LS_LU_NONZERO)
#endif
   REAL(kind=dp) :: H, Hnew, HC, HG, Fac, Tau
   REAL(kind=dp) :: Err, Yerr(N)
   INTEGER :: Pivot(N), Direction, ioffset, j, istage
   LOGICAL :: RejectLastH, RejectMoreH, Singular
!~~~>  Local parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp
!~~~>  Locally called functions
!    REAL(kind=dp) WLAMCH
!    EXTERNAL WLAMCH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~>  Initial preparations
   T = Tstart
   RSTATUS(Nhexit) = ZERO
   H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   IF (ABS(H) <= 10.0_dp*Roundoff) H = DeltaMin

   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF
   H = Direction*H

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.

!~~~> Time loop begins below

TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) )

   IF ( ISTATUS(Nstp) > Max_no_steps ) THEN  ! Too many steps
      CALL ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
      CALL ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF

!~~~>  Limit H if necessary to avoid going beyond Tend
   H = MIN(H,ABS(Tend-T))

!~~~>   Compute the function at current time
   CALL FunTemplate_4(T,Y,Fcn0, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1

!~~~>  Compute the function derivative with respect to T
   IF (.NOT.Autonomous) THEN
      CALL ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
   END IF

!~~~>   Compute the Jacobian at current time
   CALL JacTemplate_4(T,Y,Jac0,LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)   
   ISTATUS(Njac) = ISTATUS(Njac) + 1

!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO

   CALL ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular)
   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       CALL ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S

      ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
       ioffset = N*(istage-1)

      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         !slim: CALL WCOPY(N,Fcn0,1,Fcn,1)
	 Fcn(1:N) = Fcn0(1:N)
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
         !slim: CALL WCOPY(N,Y,1,Ynew,1)
	 Ynew(1:N) = Y(1:N)
         DO j = 1, istage-1
           CALL WAXPY(N,ros_A((istage-1)*(istage-2)/2+j), &
            K(N*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL FunTemplate_4(Tau,Ynew,Fcn, LS_NSEL, LS_type)
         ISTATUS(Nfun) = ISTATUS(Nfun) + 1
       END IF ! if istage == 1 elseif ros_NewF(istage)
       !slim: CALL WCOPY(N,Fcn,1,K(ioffset+1),1)
       K(ioffset+1:ioffset+N) = Fcn(1:N)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL WAXPY(N,HC,K(N*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL WAXPY(N,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL ros_Solve(Ghimj, Pivot, K(ioffset+1))

   END DO Stage


!~~~>  Compute the new solution
   !slim: CALL WCOPY(N,Y,1,Ynew,1)
   Ynew(1:N) = Y(1:N)
   DO j=1,ros_S
         CALL WAXPY(N,ros_M(j),K(N*(j-1)+1),1,Ynew,1)
   END DO

!~~~>  Compute the error estimation
   !slim: CALL WSCAL(N,ZERO,Yerr,1)
   Yerr(1:N) = ZERO
   DO j=1,ros_S
        CALL WAXPY(N,ros_E(j),K(N*(j-1)+1),1,Yerr,1)
   END DO
   Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac

!~~~>  Check the error magnitude and adjust step size
   ISTATUS(Nstp) = ISTATUS(Nstp) + 1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
      ISTATUS(Nacc) = ISTATUS(Nacc) + 1
      !slim: CALL WCOPY(N,Ynew,1,Y,1)
      Y(1:N) = Ynew(1:N)
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H)
      END IF
      RSTATUS(Nhexit) = H
      RSTATUS(Nhnew)  = Hnew
      RSTATUS(Ntexit) = T
      RejectLastH = .FALSE.
      RejectMoreH = .FALSE.
      H = Hnew
      EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
   ELSE           !~~~> Reject step
      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (ISTATUS(Nacc) >= 1)  ISTATUS(Nrej) = ISTATUS(Nrej) + 1
   END IF ! Err <= 1

   END DO UntilAccepted

   
   END DO TimeLoop

!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful

  END SUBROUTINE ros_Integrator


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  REAL(kind=dp) FUNCTION ros_ErrorNorm ( Y, Ynew, Yerr, &
                               AbsTol, RelTol, VectorTol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

! Input arguments
   REAL(kind=dp), INTENT(IN) :: Y(N), Ynew(N), &
          Yerr(N), AbsTol(N), RelTol(N)
   LOGICAL, INTENT(IN) ::  VectorTol
! Local variables
   REAL(kind=dp) :: Err, Scale, Ymax
   INTEGER  :: i
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp

   Err = ZERO
   DO i=1,N
     Ymax = MAX(ABS(Y(i)),ABS(Ynew(i)))
     IF (VectorTol) THEN
       Scale = AbsTol(i)+RelTol(i)*Ymax
     ELSE
       Scale = AbsTol(1)+RelTol(1)*Ymax
     END IF
     Err = Err+(Yerr(i)/Scale)**2
   END DO
   Err  = SQRT(Err/N)

   ros_ErrorNorm = MAX(Err,1.0d-10)

  END FUNCTION ros_ErrorNorm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> The time partial derivative of the function by finite differences
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

!~~~> Input arguments
   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(N), Fcn0(N)
!~~~> Output arguments
   REAL(kind=dp), INTENT(OUT) :: dFdT(N)
!~~~> Local variables
   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL FunTemplate_4(T+Delta,Y,dFdT, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1
   CALL WAXPY(N,(-ONE),Fcn0,1,dFdT,1)
   CALL WSCAL(N,(ONE/Delta),dFdT,1)

  END SUBROUTINE ros_FunTimeDerivative


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_PrepareMatrix ( H, Direction, gam, &
             Jac0, Ghimj, Pivot, Singular )
! --- --- --- --- --- --- --- --- --- --- --- --- ---
!  Prepares the LHS matrix for stage calculations
!  1.  Construct Ghimj = 1/(H*ham) - Jac0
!      "(Gamma H) Inverse Minus Jacobian"
!  2.  Repeat LU decomposition of Ghimj until successful.
!       -half the step size if LU decomposition fails and retry
!       -exit after 5 consecutive fails
! --- --- --- --- --- --- --- --- --- --- --- --- ---
   IMPLICIT NONE

!~~~> Input arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) ::  Jac0(N,N)
#else
   REAL(kind=dp), INTENT(IN) ::  Jac0(LS_LU_NONZERO)
#endif   
   REAL(kind=dp), INTENT(IN) ::  gam
   INTEGER, INTENT(IN) ::  Direction
!~~~> Output arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(OUT) :: Ghimj(N,N)
#else
   REAL(kind=dp), INTENT(OUT) :: Ghimj(LS_LU_NONZERO)
#endif   
   LOGICAL, INTENT(OUT) ::  Singular
   INTEGER, INTENT(OUT) ::  Pivot(N)
!~~~> Inout arguments
   REAL(kind=dp), INTENT(INOUT) :: H   ! step size is decreased when LU fails
!~~~> Local variables
   INTEGER  :: i, ISING, Nconsecutive
   REAL(kind=dp) :: ghinv
   REAL(kind=dp), PARAMETER :: ONE  = 1.0_dp, HALF = 0.5_dp

   Nconsecutive = 0
   Singular = .TRUE.

   DO WHILE (Singular)

!~~~>    Construct Ghimj = 1/(H*gam) - Jac0
#ifdef FULL_ALGEBRA    
     !slim: CALL WCOPY(N*N,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(N*N,(-ONE),Ghimj,1)
     Ghimj = -Jac0
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(i,i) = Ghimj(i,i)+ghinv
     END DO
#else
     !slim: CALL WCOPY(LS_LU_NONZERO,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(LS_LU_NONZERO,(-ONE),Ghimj,1)
     Ghimj(1:LS_LU_NONZERO) = -Jac0(1:LS_LU_NONZERO)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(LS_LU_DIAG(i)) = Ghimj(LS_LU_DIAG(i))+ghinv
     END DO
#endif   
!~~~>    Compute LU decomposition
     CALL ros_Decomp( Ghimj, Pivot, ISING )
     IF (ISING == 0) THEN
!~~~>    If successful done
        Singular = .FALSE.
     ELSE ! ISING .ne. 0
!~~~>    If unsuccessful half the step size; if 5 consecutive fails then return
        ISTATUS(Nsng) = ISTATUS(Nsng) + 1
        Nconsecutive = Nconsecutive+1
        Singular = .TRUE.
        PRINT*,'Warning: LU Decomposition returned ISING = ',ISING
        IF (Nconsecutive <= 5) THEN ! Less than 5 consecutive failed decompositions
           H = H*HALF
        ELSE  ! More than 5 consecutive failed decompositions
           RETURN
        END IF  ! Nconsecutive
      END IF    ! ISING

   END DO ! WHILE Singular

  END SUBROUTINE ros_PrepareMatrix


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Decomp( A, Pivot, ISING )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the LU decomposition
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Inout variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(INOUT) :: A(N,N)
#else   
   REAL(kind=dp), INTENT(INOUT) :: A(LS_LU_NONZERO)
#endif
!~~~> Output variables
   INTEGER, INTENT(OUT) :: Pivot(N), ISING

#ifdef FULL_ALGEBRA    
   CALL  DGETRF( N, N, A, N, Pivot, ISING )
#else   
   CALL KppDecomp ( A, ISING,LS_NSEL,LS_LU_NONZERO,LS_LU_CROW,LS_LU_DIAG,LS_LU_ICOL)   
   Pivot(1) = 1
#endif
   ISTATUS(Ndec) = ISTATUS(Ndec) + 1

  END SUBROUTINE ros_Decomp


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Solve( A, Pivot, b)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the forward/backward substitution (using pre-computed LU decomposition)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Input variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) :: A(N,N)
   INTEGER :: ISING
#else   
   REAL(kind=dp), INTENT(IN) :: A(LS_LU_NONZERO)
#endif
   INTEGER, INTENT(IN) :: Pivot(N)
!~~~> InOut variables
   REAL(kind=dp), INTENT(INOUT) :: b(N)

#ifdef FULL_ALGEBRA    
   CALL  DGETRS( 'N', N , 1, A, N, Pivot, b, N, ISING )
   IF ( Info < 0 ) THEN
      PRINT*,"Error in DGETRS. ISING=",ISING
   END IF  
#else  
   	CALL KppSolve_4(LS_LU_NONZERO,LS_NSEL, A, b,LS_type )
#endif

   ISTATUS(Nsol) = ISTATUS(Nsol) + 1

  END SUBROUTINE ros_Solve



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   DOUBLE PRECISION g

    g = 1.0_dp + 1.0_dp/SQRT(2.0_dp)
    rosMethod = RS2
!~~~> Name of the method
    ros_Name = 'ROS-2'
!~~~> Number of stages
    ros_S = 2

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = (1.0_dp)/g
    ros_C(1) = (-2.0_dp)/g
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
!~~~> M_i = Coefficients for new step solution
    ros_M(1)= (3.0_dp)/(2.0_dp*g)
    ros_M(2)= (1.0_dp)/(2.0_dp*g)
! E_i = Coefficients for error estimator
    ros_E(1) = 1.0_dp/(2.0_dp*g)
    ros_E(2) = 1.0_dp/(2.0_dp*g)
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus one
    ros_ELO = 2.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.0_dp
    ros_Alpha(2) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = g
    ros_Gamma(2) =-g

 END SUBROUTINE Ros2


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   rosMethod = RS3
!~~~> Name of the method
   ros_Name = 'ROS-3'
!~~~> Number of stages
   ros_S = 3

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1)= 1.0_dp
   ros_A(2)= 1.0_dp
   ros_A(3)= 0.0_dp

   ros_C(1) = -0.10156171083877702091975600115545E+01_dp
   ros_C(2) =  0.40759956452537699824805835358067E+01_dp
   ros_C(3) =  0.92076794298330791242156818474003E+01_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) =  0.1E+01_dp
   ros_M(2) =  0.61697947043828245592553615689730E+01_dp
   ros_M(3) = -0.42772256543218573326238373806514_dp
! E_i = Coefficients for error estimator
   ros_E(1) =  0.5_dp
   ros_E(2) = -0.29079558716805469821718236208017E+01_dp
   ros_E(3) =  0.22354069897811569627360909276199_dp
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1)= 0.0_dp
   ros_Alpha(2)= 0.43586652150845899941601945119356_dp
   ros_Alpha(3)= 0.43586652150845899941601945119356_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1)= 0.43586652150845899941601945119356_dp
   ros_Gamma(2)= 0.24291996454816804366592249683314_dp
   ros_Gamma(3)= 0.21851380027664058511513169485832E+01_dp

  END SUBROUTINE Ros3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
!     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1990)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RS4
!~~~> Name of the method
   ros_Name = 'ROS-4'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.2000000000000000E+01_dp
   ros_A(2) = 0.1867943637803922E+01_dp
   ros_A(3) = 0.2344449711399156_dp
   ros_A(4) = ros_A(2)
   ros_A(5) = ros_A(3)
   ros_A(6) = 0.0_dp

   ros_C(1) =-0.7137615036412310E+01_dp
   ros_C(2) = 0.2580708087951457E+01_dp
   ros_C(3) = 0.6515950076447975_dp
   ros_C(4) =-0.2137148994382534E+01_dp
   ros_C(5) =-0.3214669691237626_dp
   ros_C(6) =-0.6949742501781779_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .TRUE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 0.2255570073418735E+01_dp
   ros_M(2) = 0.2870493262186792_dp
   ros_M(3) = 0.4353179431840180_dp
   ros_M(4) = 0.1093502252409163E+01_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) =-0.2815431932141155_dp
   ros_E(2) =-0.7276199124938920E-01_dp
   ros_E(3) =-0.1082196201495311_dp
   ros_E(4) =-0.1093502252409163E+01_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 4.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.1145640000000000E+01_dp
   ros_Alpha(3) = 0.6552168638155900_dp
   ros_Alpha(4) = ros_Alpha(3)
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5728200000000000_dp
   ros_Gamma(2) =-0.1769193891319233E+01_dp
   ros_Gamma(3) = 0.7592633437920482_dp
   ros_Gamma(4) =-0.1049021087100450_dp

  END SUBROUTINE Ros4

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- A STIFFLY-STABLE METHOD, 4 stages, order 3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RD3
!~~~> Name of the method
   ros_Name = 'RODAS-3'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.0_dp
   ros_A(2) = 2.0_dp
   ros_A(3) = 0.0_dp
   ros_A(4) = 2.0_dp
   ros_A(5) = 0.0_dp
   ros_A(6) = 1.0_dp

   ros_C(1) = 4.0_dp
   ros_C(2) = 1.0_dp
   ros_C(3) =-1.0_dp
   ros_C(4) = 1.0_dp
   ros_C(5) =-1.0_dp
   ros_C(6) =-(8.0_dp/3.0_dp)

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .FALSE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .TRUE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 2.0_dp
   ros_M(2) = 0.0_dp
   ros_M(3) = 1.0_dp
   ros_M(4) = 1.0_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) = 0.0_dp
   ros_E(2) = 0.0_dp
   ros_E(3) = 0.0_dp
   ros_E(4) = 1.0_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.0_dp
   ros_Alpha(3) = 1.0_dp
   ros_Alpha(4) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5_dp
   ros_Gamma(2) = 1.5_dp
   ros_Gamma(3) = 0.0_dp
   ros_Gamma(4) = 0.0_dp

  END SUBROUTINE Rodas3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1996)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RD4
!~~~> Name of the method
    ros_Name = 'RODAS-4'
!~~~> Number of stages
    ros_S = 6

!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.000_dp
    ros_Alpha(2) = 0.386_dp
    ros_Alpha(3) = 0.210_dp
    ros_Alpha(4) = 0.630_dp
    ros_Alpha(5) = 1.000_dp
    ros_Alpha(6) = 1.000_dp

!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = 0.2500000000000000_dp
    ros_Gamma(2) =-0.1043000000000000_dp
    ros_Gamma(3) = 0.1035000000000000_dp
    ros_Gamma(4) =-0.3620000000000023E-01_dp
    ros_Gamma(5) = 0.0_dp
    ros_Gamma(6) = 0.0_dp

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = 0.1544000000000000E+01_dp
    ros_A(2) = 0.9466785280815826_dp
    ros_A(3) = 0.2557011698983284_dp
    ros_A(4) = 0.3314825187068521E+01_dp
    ros_A(5) = 0.2896124015972201E+01_dp
    ros_A(6) = 0.9986419139977817_dp
    ros_A(7) = 0.1221224509226641E+01_dp
    ros_A(8) = 0.6019134481288629E+01_dp
    ros_A(9) = 0.1253708332932087E+02_dp
    ros_A(10) =-0.6878860361058950_dp
    ros_A(11) = ros_A(7)
    ros_A(12) = ros_A(8)
    ros_A(13) = ros_A(9)
    ros_A(14) = ros_A(10)
    ros_A(15) = 1.0_dp

    ros_C(1) =-0.5668800000000000E+01_dp
    ros_C(2) =-0.2430093356833875E+01_dp
    ros_C(3) =-0.2063599157091915_dp
    ros_C(4) =-0.1073529058151375_dp
    ros_C(5) =-0.9594562251023355E+01_dp
    ros_C(6) =-0.2047028614809616E+02_dp
    ros_C(7) = 0.7496443313967647E+01_dp
    ros_C(8) =-0.1024680431464352E+02_dp
    ros_C(9) =-0.3399990352819905E+02_dp
    ros_C(10) = 0.1170890893206160E+02_dp
    ros_C(11) = 0.8083246795921522E+01_dp
    ros_C(12) =-0.7981132988064893E+01_dp
    ros_C(13) =-0.3152159432874371E+02_dp
    ros_C(14) = 0.1631930543123136E+02_dp
    ros_C(15) =-0.6058818238834054E+01_dp

!~~~> M_i = Coefficients for new step solution
    ros_M(1) = ros_A(7)
    ros_M(2) = ros_A(8)
    ros_M(3) = ros_A(9)
    ros_M(4) = ros_A(10)
    ros_M(5) = 1.0_dp
    ros_M(6) = 1.0_dp

!~~~> E_i  = Coefficients for error estimator
    ros_E(1) = 0.0_dp
    ros_E(2) = 0.0_dp
    ros_E(3) = 0.0_dp
    ros_E(4) = 0.0_dp
    ros_E(5) = 0.0_dp
    ros_E(6) = 1.0_dp

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.
    ros_NewF(5) = .TRUE.
    ros_NewF(6) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 4.0_dp

  END SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! STIFFLY-STABLE W METHOD OF ORDER 3, WITH 4 STAGES
!
! J. RANG and L. ANGERMANN
! NEW ROSENBROCK W-METHODS OF ORDER 3
! FOR PARTIAL DIFFERENTIAL ALGEBRAIC
!        EQUATIONS OF INDEX 1                                             
! BIT Numerical Mathematics (2005) 45: 761-787
!  DOI: 10.1007/s10543-005-0035-y
! Table 4.1-4.2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RG3
!~~~> Name of the method
    ros_Name = 'RANG-3'
!~~~> Number of stages
    ros_S = 4

    ros_A(1) = 5.09052051067020d+00;
    ros_A(2) = 5.09052051067020d+00;
    ros_A(3) = 0.0d0;
    ros_A(4) = 4.97628111010787d+00;
    ros_A(5) = 2.77268164715849d-02;
    ros_A(6) = 2.29428036027904d-01;

    ros_C(1) = -1.16790812312283d+01;
    ros_C(2) = -1.64057326467367d+01;
    ros_C(3) = -2.77268164715850d-01;
    ros_C(4) = -8.38103960500476d+00;
    ros_C(5) = -8.48328409199343d-01;
    ros_C(6) =  2.87009860433106d-01;

    ros_M(1) =  5.22582761233094d+00;
    ros_M(2) = -5.56971148154165d-01;
    ros_M(3) =  3.57979469353645d-01;
    ros_M(4) =  1.72337398521064d+00;

    ros_E(1) = -5.16845212784040d+00;
    ros_E(2) = -1.26351942603842d+00;
    ros_E(3) = -1.11022302462516d-16;
    ros_E(4) =  2.22044604925031d-16;

    ros_Alpha(1) = 0.0d00;
    ros_Alpha(2) = 2.21878746765329d+00;
    ros_Alpha(3) = 2.21878746765329d+00;
    ros_Alpha(4) = 1.55392337535788d+00;

    ros_Gamma(1) =  4.35866521508459d-01;
    ros_Gamma(2) = -1.78292094614483d+00;
    ros_Gamma(3) = -2.46541900496934d+00;
    ros_Gamma(4) = -8.05529997906370d-01;

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 3.0_dp

  END SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END SUBROUTINE Rosenbrock_4


SUBROUTINE FunTemplate_4( T, Y, Ydot, LS_NSEL, LS_type )
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Function
 INTEGER,INTENT(IN)::LS_NSEL, LS_type
   REAL(kind=dp) :: T, Y(LS_NSEL)
   REAL(kind=dp) :: Ydot(LS_NSEL)
   REAL(kind=dp) :: Told
   Told = TIME
   TIME = T
   CALL Fun_4( Y, FIX, RCONST, Ydot, LS_NSEL)
   TIME = Told
END SUBROUTINE FunTemplate_4

SUBROUTINE JacTemplate_4( T, Y, Jcb, LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Jacobian
 USE gckpp_LinearAlgebra
 
 INTEGER,INTENT(IN)::LS_NSEL, LS_LU_NONZERO, LS_type
 INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO), LS_LU_ICOL(LS_LU_NONZERO)
    REAL(kind=dp) :: T, Y(LS_NSEL)
#ifdef FULL_ALGEBRA    
    REAL(kind=dp) :: JV(LS_LU_NONZERO), Jcb(LS_NSEL,LS_NSEL)
#else
    REAL(kind=dp) :: Jcb(LS_LU_NONZERO)
#endif   
    REAL(kind=dp) :: Told
#ifdef FULL_ALGEBRA    
    INTEGER :: i, j
#endif   

    Told = TIME
    TIME = T
#ifdef FULL_ALGEBRA        	
	CALL Jac_SP_4(Y, FIX, RCONST, JV,LS_NSEL,LS_LU_NONZERO)	
    DO j=1,LS_NSEL
      DO i=1,LS_NSEL
         Jcb(i,j) = 0.0_dp
      END DO
    END DO
    DO i=1,LS_LU_NONZERO
       Jcb(LS_LU_IROW(i),LS_LU_ICOL(i)) = JV(i)
    END DO
#else   
       CALL Jac_SP_4( Y, FIX, RCONST, Jcb,LS_NSEL,LS_LU_NONZERO)	 
#endif   
    TIME = Told
END SUBROUTINE JacTemplate_4



SUBROUTINE INTEGRATE_5( TIN, TOUT, LS_type,LS_NSEL, LS_NDEL,&
  Prate, Lrate, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )

   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: TIN  ! Start Time
   REAL(kind=dp), INTENT(IN) :: TOUT ! End Time
   
   ! Optional input parameters and statistics
   INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   REAL(kind=dp), INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   REAL(kind=dp), INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U
   INTEGER,		  INTENT(IN) :: LS_type,LS_NSEL,LS_NDEL
   REAL(kind=dp), INTENT(IN) :: Prate(NVAR),Lrate(NVAR)
   
   REAL(kind=dp) :: RCNTRL(20), RSTATUS(20),deltaT
   INTEGER       :: ICNTRL(20), ISTATUS(20), IERR
   REAL(kind=dp) :: VAR_selected(LS_NSEL),VAR_deleted(LS_NDEL),LS_P(LS_NDEL),LS_L(LS_NDEL)
   
   INTEGER, SAVE :: Ntotal = 0

   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.0_dp
   ISTATUS(:) = 0
   RSTATUS(:) = 0.0_dp
   deltaT =  TOUT-TIN
   !~~~> fine-tune the integrator:
   ICNTRL(1) = 0	! 0 - non-autonomous, 1 - autonomous
   ICNTRL(2) = 0	! 0 - vector tolerances, 1 - scalars

   ! If optional parameters are given, and if they are >0, 
   ! then they overwrite default settings. 
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF
   
        VAR_selected=VAR(select_ind_5)
	    VAR_deleted=VAR(delete_ind_5)
		LS_P=Prate(delete_ind_5)
		LS_L=Lrate(delete_ind_5)
        CALL Rosenbrock_5(LU_NSEL_5,VAR_selected,TIN,TOUT,ATOL,RTOL,&
          RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		  LU_NONZERO_5,LU_NSEL_5,LU_CROW_5,LU_DIAG_5,LU_IROW_5,LU_ICOL_5, LS_type)		
		WHERE(LS_L<=(0.01/deltaT))
			VAR_deleted=VAR_deleted+deltaT*(LS_P-LS_L*VAR_deleted)
		ELSEWHERE
			VAR_deleted=LS_P/LS_L+(VAR_deleted-LS_P/LS_L)*EXP(-LS_L*deltaT)
		END WHERE

	    VAR(select_ind_5)=VAR_selected
	    VAR(delete_ind_5)=VAR_deleted
   
   !~~~> Debug option: show no of steps
   ! Ntotal = Ntotal + ISTATUS(Nstp)
   ! PRINT*,'NSTEPS=',ISTATUS(Nstp),' (',Ntotal,')','  O3=', VAR(ind_O3)

   STEPMIN = RSTATUS(Nhexit)
   ! if optional parameters are given for output they 
   ! are updated with the return information
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE INTEGRATE_5

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE Rosenbrock_5(N,Y,Tstart,Tend, &
           AbsTol,RelTol,              &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		   LS_LU_NONZERO,LS_NSEL,LS_LU_CROW,LS_LU_DIAG,LS_LU_IROW,LS_LU_ICOL, LS_type)
  USE gckpp_Parameters
  USE gckpp_LinearAlgebra
  IMPLICIT NONE

   INTEGER,       INTENT(IN)    :: N,LS_LU_NONZERO,LS_NSEL,LS_type   
   INTEGER,INTENT(IN)::LS_LU_CROW(LS_NSEL+1),LS_LU_DIAG(LS_NSEL+1)
   INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO),LS_LU_ICOL(LS_LU_NONZERO)   
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
   REAL(kind=dp), INTENT(IN)    :: Tstart,Tend
   REAL(kind=dp), INTENT(IN)    :: AbsTol(N),RelTol(N)
   INTEGER,       INTENT(IN)    :: ICNTRL(20)
   REAL(kind=dp), INTENT(IN)    :: RCNTRL(20)
   INTEGER,       INTENT(INOUT) :: ISTATUS(20)
   REAL(kind=dp), INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT)   :: IERR
!~~~>  Parameters of the Rosenbrock method, up to 6 stages
   INTEGER ::  ros_S, rosMethod
   INTEGER, PARAMETER :: RS2=1, RS3=2, RS4=3, RD3=4, RD4=5, RG3=6
   REAL(kind=dp) :: ros_A(15), ros_C(15), ros_M(6), ros_E(6), &
                    ros_Alpha(6), ros_Gamma(6), ros_ELO
   LOGICAL :: ros_NewF(6)
   CHARACTER(LEN=12) :: ros_Name
!~~~>  Local variables
   REAL(kind=dp) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   REAL(kind=dp) :: Hmin, Hmax, Hstart
   REAL(kind=dp) :: Texit
   INTEGER       :: i, UplimTol, Max_no_steps
   LOGICAL       :: Autonomous, VectorTol
!~~~>   Parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp

!~~~>  Initialize statistics
   ISTATUS(1:8) = 0
   RSTATUS(1:3) = ZERO

!~~~>  Autonomous or time dependent ODE. Default is time dependent.
   Autonomous = .NOT.(ICNTRL(1) == 0)

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:N) and RelTol(1:N)
   IF (ICNTRL(2) == 0) THEN
      VectorTol = .TRUE.
      UplimTol  = N
   ELSE
      VectorTol = .FALSE.
      UplimTol  = 1
   END IF

!~~~>   Initialize the particular Rosenbrock method selected
   SELECT CASE (ICNTRL(3))
     CASE (1)
       CALL Ros2
     CASE (2)
       CALL Ros3
     CASE (3)
       CALL Ros4
     CASE (0,4)
       CALL Rodas3
     CASE (5)
       CALL Rodas4
     CASE (6)
       CALL Rang3
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(3)=',ICNTRL(3) 
       CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT

!~~~>   The maximum number of steps admitted
   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 200000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Unit roundoff (1+Roundoff>1)
   Roundoff = WLAMCH('E')

!~~~>  Lower bound on the step size: (positive value)
   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Upper bound on the step size: (positive value)
   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Starting step size: (positive value)
   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hold < FacMax
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Check if tolerances are reasonable
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO


!~~~>  CALL Rosenbrock method
   CALL ros_Integrator(Y, Tstart, Tend, Texit,   &
        AbsTol, RelTol,                          &
!  Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!  Error indicator
        IERR)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS !  SUBROUTINES internal to Rosenbrock
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 SUBROUTINE ros_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
   REAL(kind=dp), INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR
   
   IERR = Code
   PRINT * , &
     'Forced exit from Rosenbrock due to the following error:' 
     
   SELECT CASE (Code)
    CASE (-1)
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)    
      PRINT * , '--> Selected Rosenbrock method not implemented'
    CASE (-3)    
      PRINT * , '--> Hmin/Hmax/Hstart must be positive'
    CASE (-4)    
      PRINT * , '--> FacMin/FacMax/FacRej must be positive'
    CASE (-5) 
      PRINT * , '--> Improper tolerance values'
    CASE (-6) 
      PRINT * , '--> No of steps exceeds maximum bound'
    CASE (-7) 
      PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
    CASE (-8)    
      PRINT * , '--> Matrix is repeatedly singular'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT
   
   PRINT *, "T=", T, "and H=", H
     
 END SUBROUTINE ros_ErrorMsg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_Integrator (Y, Tstart, Tend, T,  &
        AbsTol, RelTol,                          &
!~~~> Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic Rosenbrock method
!      defined by ros_S (no of stages)
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

!~~~> Input: the initial condition at Tstart; Output: the solution at T
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
!~~~> Input: integration interval
   REAL(kind=dp), INTENT(IN) :: Tstart,Tend
!~~~> Output: time at which the solution is returned (T=Tend if success)
   REAL(kind=dp), INTENT(OUT) ::  T
!~~~> Input: tolerances
   REAL(kind=dp), INTENT(IN) ::  AbsTol(N), RelTol(N)
!~~~> Input: integration parameters
   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   REAL(kind=dp), INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   REAL(kind=dp), INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables
   REAL(kind=dp) :: Ynew(N), Fcn0(N), Fcn(N)
   REAL(kind=dp) :: K(N*ros_S), dFdT(N)
#ifdef FULL_ALGEBRA    
   REAL(kind=dp) :: Jac0(N,N), Ghimj(N,N)
#else
   REAL(kind=dp) :: Jac0(LS_LU_NONZERO), Ghimj(LS_LU_NONZERO)
#endif
   REAL(kind=dp) :: H, Hnew, HC, HG, Fac, Tau
   REAL(kind=dp) :: Err, Yerr(N)
   INTEGER :: Pivot(N), Direction, ioffset, j, istage
   LOGICAL :: RejectLastH, RejectMoreH, Singular
!~~~>  Local parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp
!~~~>  Locally called functions
!    REAL(kind=dp) WLAMCH
!    EXTERNAL WLAMCH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~>  Initial preparations
   T = Tstart
   RSTATUS(Nhexit) = ZERO
   H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   IF (ABS(H) <= 10.0_dp*Roundoff) H = DeltaMin

   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF
   H = Direction*H

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.

!~~~> Time loop begins below

TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) )

   IF ( ISTATUS(Nstp) > Max_no_steps ) THEN  ! Too many steps
      CALL ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
      CALL ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF

!~~~>  Limit H if necessary to avoid going beyond Tend
   H = MIN(H,ABS(Tend-T))

!~~~>   Compute the function at current time
   CALL FunTemplate_5(T,Y,Fcn0, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1

!~~~>  Compute the function derivative with respect to T
   IF (.NOT.Autonomous) THEN
      CALL ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
   END IF

!~~~>   Compute the Jacobian at current time
   CALL JacTemplate_5(T,Y,Jac0,LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)   
   ISTATUS(Njac) = ISTATUS(Njac) + 1

!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO

   CALL ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular)
   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       CALL ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S

      ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
       ioffset = N*(istage-1)

      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         !slim: CALL WCOPY(N,Fcn0,1,Fcn,1)
	 Fcn(1:N) = Fcn0(1:N)
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
         !slim: CALL WCOPY(N,Y,1,Ynew,1)
	 Ynew(1:N) = Y(1:N)
         DO j = 1, istage-1
           CALL WAXPY(N,ros_A((istage-1)*(istage-2)/2+j), &
            K(N*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL FunTemplate_5(Tau,Ynew,Fcn, LS_NSEL, LS_type)
         ISTATUS(Nfun) = ISTATUS(Nfun) + 1
       END IF ! if istage == 1 elseif ros_NewF(istage)
       !slim: CALL WCOPY(N,Fcn,1,K(ioffset+1),1)
       K(ioffset+1:ioffset+N) = Fcn(1:N)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL WAXPY(N,HC,K(N*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL WAXPY(N,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL ros_Solve(Ghimj, Pivot, K(ioffset+1))

   END DO Stage


!~~~>  Compute the new solution
   !slim: CALL WCOPY(N,Y,1,Ynew,1)
   Ynew(1:N) = Y(1:N)
   DO j=1,ros_S
         CALL WAXPY(N,ros_M(j),K(N*(j-1)+1),1,Ynew,1)
   END DO

!~~~>  Compute the error estimation
   !slim: CALL WSCAL(N,ZERO,Yerr,1)
   Yerr(1:N) = ZERO
   DO j=1,ros_S
        CALL WAXPY(N,ros_E(j),K(N*(j-1)+1),1,Yerr,1)
   END DO
   Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac

!~~~>  Check the error magnitude and adjust step size
   ISTATUS(Nstp) = ISTATUS(Nstp) + 1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
      ISTATUS(Nacc) = ISTATUS(Nacc) + 1
      !slim: CALL WCOPY(N,Ynew,1,Y,1)
      Y(1:N) = Ynew(1:N)
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H)
      END IF
      RSTATUS(Nhexit) = H
      RSTATUS(Nhnew)  = Hnew
      RSTATUS(Ntexit) = T
      RejectLastH = .FALSE.
      RejectMoreH = .FALSE.
      H = Hnew
      EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
   ELSE           !~~~> Reject step
      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (ISTATUS(Nacc) >= 1)  ISTATUS(Nrej) = ISTATUS(Nrej) + 1
   END IF ! Err <= 1

   END DO UntilAccepted

   
   END DO TimeLoop

!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful

  END SUBROUTINE ros_Integrator


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  REAL(kind=dp) FUNCTION ros_ErrorNorm ( Y, Ynew, Yerr, &
                               AbsTol, RelTol, VectorTol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

! Input arguments
   REAL(kind=dp), INTENT(IN) :: Y(N), Ynew(N), &
          Yerr(N), AbsTol(N), RelTol(N)
   LOGICAL, INTENT(IN) ::  VectorTol
! Local variables
   REAL(kind=dp) :: Err, Scale, Ymax
   INTEGER  :: i
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp

   Err = ZERO
   DO i=1,N
     Ymax = MAX(ABS(Y(i)),ABS(Ynew(i)))
     IF (VectorTol) THEN
       Scale = AbsTol(i)+RelTol(i)*Ymax
     ELSE
       Scale = AbsTol(1)+RelTol(1)*Ymax
     END IF
     Err = Err+(Yerr(i)/Scale)**2
   END DO
   Err  = SQRT(Err/N)

   ros_ErrorNorm = MAX(Err,1.0d-10)

  END FUNCTION ros_ErrorNorm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> The time partial derivative of the function by finite differences
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

!~~~> Input arguments
   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(N), Fcn0(N)
!~~~> Output arguments
   REAL(kind=dp), INTENT(OUT) :: dFdT(N)
!~~~> Local variables
   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL FunTemplate_5(T+Delta,Y,dFdT, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1
   CALL WAXPY(N,(-ONE),Fcn0,1,dFdT,1)
   CALL WSCAL(N,(ONE/Delta),dFdT,1)

  END SUBROUTINE ros_FunTimeDerivative


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_PrepareMatrix ( H, Direction, gam, &
             Jac0, Ghimj, Pivot, Singular )
! --- --- --- --- --- --- --- --- --- --- --- --- ---
!  Prepares the LHS matrix for stage calculations
!  1.  Construct Ghimj = 1/(H*ham) - Jac0
!      "(Gamma H) Inverse Minus Jacobian"
!  2.  Repeat LU decomposition of Ghimj until successful.
!       -half the step size if LU decomposition fails and retry
!       -exit after 5 consecutive fails
! --- --- --- --- --- --- --- --- --- --- --- --- ---
   IMPLICIT NONE

!~~~> Input arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) ::  Jac0(N,N)
#else
   REAL(kind=dp), INTENT(IN) ::  Jac0(LS_LU_NONZERO)
#endif   
   REAL(kind=dp), INTENT(IN) ::  gam
   INTEGER, INTENT(IN) ::  Direction
!~~~> Output arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(OUT) :: Ghimj(N,N)
#else
   REAL(kind=dp), INTENT(OUT) :: Ghimj(LS_LU_NONZERO)
#endif   
   LOGICAL, INTENT(OUT) ::  Singular
   INTEGER, INTENT(OUT) ::  Pivot(N)
!~~~> Inout arguments
   REAL(kind=dp), INTENT(INOUT) :: H   ! step size is decreased when LU fails
!~~~> Local variables
   INTEGER  :: i, ISING, Nconsecutive
   REAL(kind=dp) :: ghinv
   REAL(kind=dp), PARAMETER :: ONE  = 1.0_dp, HALF = 0.5_dp

   Nconsecutive = 0
   Singular = .TRUE.

   DO WHILE (Singular)

!~~~>    Construct Ghimj = 1/(H*gam) - Jac0
#ifdef FULL_ALGEBRA    
     !slim: CALL WCOPY(N*N,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(N*N,(-ONE),Ghimj,1)
     Ghimj = -Jac0
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(i,i) = Ghimj(i,i)+ghinv
     END DO
#else
     !slim: CALL WCOPY(LS_LU_NONZERO,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(LS_LU_NONZERO,(-ONE),Ghimj,1)
     Ghimj(1:LS_LU_NONZERO) = -Jac0(1:LS_LU_NONZERO)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(LS_LU_DIAG(i)) = Ghimj(LS_LU_DIAG(i))+ghinv
     END DO
#endif   
!~~~>    Compute LU decomposition
     CALL ros_Decomp( Ghimj, Pivot, ISING )
     IF (ISING == 0) THEN
!~~~>    If successful done
        Singular = .FALSE.
     ELSE ! ISING .ne. 0
!~~~>    If unsuccessful half the step size; if 5 consecutive fails then return
        ISTATUS(Nsng) = ISTATUS(Nsng) + 1
        Nconsecutive = Nconsecutive+1
        Singular = .TRUE.
        PRINT*,'Warning: LU Decomposition returned ISING = ',ISING
        IF (Nconsecutive <= 5) THEN ! Less than 5 consecutive failed decompositions
           H = H*HALF
        ELSE  ! More than 5 consecutive failed decompositions
           RETURN
        END IF  ! Nconsecutive
      END IF    ! ISING

   END DO ! WHILE Singular

  END SUBROUTINE ros_PrepareMatrix


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Decomp( A, Pivot, ISING )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the LU decomposition
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Inout variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(INOUT) :: A(N,N)
#else   
   REAL(kind=dp), INTENT(INOUT) :: A(LS_LU_NONZERO)
#endif
!~~~> Output variables
   INTEGER, INTENT(OUT) :: Pivot(N), ISING

#ifdef FULL_ALGEBRA    
   CALL  DGETRF( N, N, A, N, Pivot, ISING )
#else   
   CALL KppDecomp ( A, ISING,LS_NSEL,LS_LU_NONZERO,LS_LU_CROW,LS_LU_DIAG,LS_LU_ICOL)   
   Pivot(1) = 1
#endif
   ISTATUS(Ndec) = ISTATUS(Ndec) + 1

  END SUBROUTINE ros_Decomp


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Solve( A, Pivot, b)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the forward/backward substitution (using pre-computed LU decomposition)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Input variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) :: A(N,N)
   INTEGER :: ISING
#else   
   REAL(kind=dp), INTENT(IN) :: A(LS_LU_NONZERO)
#endif
   INTEGER, INTENT(IN) :: Pivot(N)
!~~~> InOut variables
   REAL(kind=dp), INTENT(INOUT) :: b(N)

#ifdef FULL_ALGEBRA    
   CALL  DGETRS( 'N', N , 1, A, N, Pivot, b, N, ISING )
   IF ( Info < 0 ) THEN
      PRINT*,"Error in DGETRS. ISING=",ISING
   END IF  
#else  
   	CALL KppSolve_5(LS_LU_NONZERO,LS_NSEL, A, b,LS_type )
#endif

   ISTATUS(Nsol) = ISTATUS(Nsol) + 1

  END SUBROUTINE ros_Solve



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   DOUBLE PRECISION g

    g = 1.0_dp + 1.0_dp/SQRT(2.0_dp)
    rosMethod = RS2
!~~~> Name of the method
    ros_Name = 'ROS-2'
!~~~> Number of stages
    ros_S = 2

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = (1.0_dp)/g
    ros_C(1) = (-2.0_dp)/g
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
!~~~> M_i = Coefficients for new step solution
    ros_M(1)= (3.0_dp)/(2.0_dp*g)
    ros_M(2)= (1.0_dp)/(2.0_dp*g)
! E_i = Coefficients for error estimator
    ros_E(1) = 1.0_dp/(2.0_dp*g)
    ros_E(2) = 1.0_dp/(2.0_dp*g)
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus one
    ros_ELO = 2.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.0_dp
    ros_Alpha(2) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = g
    ros_Gamma(2) =-g

 END SUBROUTINE Ros2


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   rosMethod = RS3
!~~~> Name of the method
   ros_Name = 'ROS-3'
!~~~> Number of stages
   ros_S = 3

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1)= 1.0_dp
   ros_A(2)= 1.0_dp
   ros_A(3)= 0.0_dp

   ros_C(1) = -0.10156171083877702091975600115545E+01_dp
   ros_C(2) =  0.40759956452537699824805835358067E+01_dp
   ros_C(3) =  0.92076794298330791242156818474003E+01_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) =  0.1E+01_dp
   ros_M(2) =  0.61697947043828245592553615689730E+01_dp
   ros_M(3) = -0.42772256543218573326238373806514_dp
! E_i = Coefficients for error estimator
   ros_E(1) =  0.5_dp
   ros_E(2) = -0.29079558716805469821718236208017E+01_dp
   ros_E(3) =  0.22354069897811569627360909276199_dp
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1)= 0.0_dp
   ros_Alpha(2)= 0.43586652150845899941601945119356_dp
   ros_Alpha(3)= 0.43586652150845899941601945119356_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1)= 0.43586652150845899941601945119356_dp
   ros_Gamma(2)= 0.24291996454816804366592249683314_dp
   ros_Gamma(3)= 0.21851380027664058511513169485832E+01_dp

  END SUBROUTINE Ros3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
!     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1990)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RS4
!~~~> Name of the method
   ros_Name = 'ROS-4'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.2000000000000000E+01_dp
   ros_A(2) = 0.1867943637803922E+01_dp
   ros_A(3) = 0.2344449711399156_dp
   ros_A(4) = ros_A(2)
   ros_A(5) = ros_A(3)
   ros_A(6) = 0.0_dp

   ros_C(1) =-0.7137615036412310E+01_dp
   ros_C(2) = 0.2580708087951457E+01_dp
   ros_C(3) = 0.6515950076447975_dp
   ros_C(4) =-0.2137148994382534E+01_dp
   ros_C(5) =-0.3214669691237626_dp
   ros_C(6) =-0.6949742501781779_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .TRUE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 0.2255570073418735E+01_dp
   ros_M(2) = 0.2870493262186792_dp
   ros_M(3) = 0.4353179431840180_dp
   ros_M(4) = 0.1093502252409163E+01_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) =-0.2815431932141155_dp
   ros_E(2) =-0.7276199124938920E-01_dp
   ros_E(3) =-0.1082196201495311_dp
   ros_E(4) =-0.1093502252409163E+01_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 4.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.1145640000000000E+01_dp
   ros_Alpha(3) = 0.6552168638155900_dp
   ros_Alpha(4) = ros_Alpha(3)
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5728200000000000_dp
   ros_Gamma(2) =-0.1769193891319233E+01_dp
   ros_Gamma(3) = 0.7592633437920482_dp
   ros_Gamma(4) =-0.1049021087100450_dp

  END SUBROUTINE Ros4

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- A STIFFLY-STABLE METHOD, 4 stages, order 3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RD3
!~~~> Name of the method
   ros_Name = 'RODAS-3'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.0_dp
   ros_A(2) = 2.0_dp
   ros_A(3) = 0.0_dp
   ros_A(4) = 2.0_dp
   ros_A(5) = 0.0_dp
   ros_A(6) = 1.0_dp

   ros_C(1) = 4.0_dp
   ros_C(2) = 1.0_dp
   ros_C(3) =-1.0_dp
   ros_C(4) = 1.0_dp
   ros_C(5) =-1.0_dp
   ros_C(6) =-(8.0_dp/3.0_dp)

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .FALSE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .TRUE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 2.0_dp
   ros_M(2) = 0.0_dp
   ros_M(3) = 1.0_dp
   ros_M(4) = 1.0_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) = 0.0_dp
   ros_E(2) = 0.0_dp
   ros_E(3) = 0.0_dp
   ros_E(4) = 1.0_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.0_dp
   ros_Alpha(3) = 1.0_dp
   ros_Alpha(4) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5_dp
   ros_Gamma(2) = 1.5_dp
   ros_Gamma(3) = 0.0_dp
   ros_Gamma(4) = 0.0_dp

  END SUBROUTINE Rodas3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1996)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RD4
!~~~> Name of the method
    ros_Name = 'RODAS-4'
!~~~> Number of stages
    ros_S = 6

!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.000_dp
    ros_Alpha(2) = 0.386_dp
    ros_Alpha(3) = 0.210_dp
    ros_Alpha(4) = 0.630_dp
    ros_Alpha(5) = 1.000_dp
    ros_Alpha(6) = 1.000_dp

!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = 0.2500000000000000_dp
    ros_Gamma(2) =-0.1043000000000000_dp
    ros_Gamma(3) = 0.1035000000000000_dp
    ros_Gamma(4) =-0.3620000000000023E-01_dp
    ros_Gamma(5) = 0.0_dp
    ros_Gamma(6) = 0.0_dp

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = 0.1544000000000000E+01_dp
    ros_A(2) = 0.9466785280815826_dp
    ros_A(3) = 0.2557011698983284_dp
    ros_A(4) = 0.3314825187068521E+01_dp
    ros_A(5) = 0.2896124015972201E+01_dp
    ros_A(6) = 0.9986419139977817_dp
    ros_A(7) = 0.1221224509226641E+01_dp
    ros_A(8) = 0.6019134481288629E+01_dp
    ros_A(9) = 0.1253708332932087E+02_dp
    ros_A(10) =-0.6878860361058950_dp
    ros_A(11) = ros_A(7)
    ros_A(12) = ros_A(8)
    ros_A(13) = ros_A(9)
    ros_A(14) = ros_A(10)
    ros_A(15) = 1.0_dp

    ros_C(1) =-0.5668800000000000E+01_dp
    ros_C(2) =-0.2430093356833875E+01_dp
    ros_C(3) =-0.2063599157091915_dp
    ros_C(4) =-0.1073529058151375_dp
    ros_C(5) =-0.9594562251023355E+01_dp
    ros_C(6) =-0.2047028614809616E+02_dp
    ros_C(7) = 0.7496443313967647E+01_dp
    ros_C(8) =-0.1024680431464352E+02_dp
    ros_C(9) =-0.3399990352819905E+02_dp
    ros_C(10) = 0.1170890893206160E+02_dp
    ros_C(11) = 0.8083246795921522E+01_dp
    ros_C(12) =-0.7981132988064893E+01_dp
    ros_C(13) =-0.3152159432874371E+02_dp
    ros_C(14) = 0.1631930543123136E+02_dp
    ros_C(15) =-0.6058818238834054E+01_dp

!~~~> M_i = Coefficients for new step solution
    ros_M(1) = ros_A(7)
    ros_M(2) = ros_A(8)
    ros_M(3) = ros_A(9)
    ros_M(4) = ros_A(10)
    ros_M(5) = 1.0_dp
    ros_M(6) = 1.0_dp

!~~~> E_i  = Coefficients for error estimator
    ros_E(1) = 0.0_dp
    ros_E(2) = 0.0_dp
    ros_E(3) = 0.0_dp
    ros_E(4) = 0.0_dp
    ros_E(5) = 0.0_dp
    ros_E(6) = 1.0_dp

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.
    ros_NewF(5) = .TRUE.
    ros_NewF(6) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 4.0_dp

  END SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! STIFFLY-STABLE W METHOD OF ORDER 3, WITH 4 STAGES
!
! J. RANG and L. ANGERMANN
! NEW ROSENBROCK W-METHODS OF ORDER 3
! FOR PARTIAL DIFFERENTIAL ALGEBRAIC
!        EQUATIONS OF INDEX 1                                             
! BIT Numerical Mathematics (2005) 45: 761-787
!  DOI: 10.1007/s10543-005-0035-y
! Table 4.1-4.2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RG3
!~~~> Name of the method
    ros_Name = 'RANG-3'
!~~~> Number of stages
    ros_S = 4

    ros_A(1) = 5.09052051067020d+00;
    ros_A(2) = 5.09052051067020d+00;
    ros_A(3) = 0.0d0;
    ros_A(4) = 4.97628111010787d+00;
    ros_A(5) = 2.77268164715849d-02;
    ros_A(6) = 2.29428036027904d-01;

    ros_C(1) = -1.16790812312283d+01;
    ros_C(2) = -1.64057326467367d+01;
    ros_C(3) = -2.77268164715850d-01;
    ros_C(4) = -8.38103960500476d+00;
    ros_C(5) = -8.48328409199343d-01;
    ros_C(6) =  2.87009860433106d-01;

    ros_M(1) =  5.22582761233094d+00;
    ros_M(2) = -5.56971148154165d-01;
    ros_M(3) =  3.57979469353645d-01;
    ros_M(4) =  1.72337398521064d+00;

    ros_E(1) = -5.16845212784040d+00;
    ros_E(2) = -1.26351942603842d+00;
    ros_E(3) = -1.11022302462516d-16;
    ros_E(4) =  2.22044604925031d-16;

    ros_Alpha(1) = 0.0d00;
    ros_Alpha(2) = 2.21878746765329d+00;
    ros_Alpha(3) = 2.21878746765329d+00;
    ros_Alpha(4) = 1.55392337535788d+00;

    ros_Gamma(1) =  4.35866521508459d-01;
    ros_Gamma(2) = -1.78292094614483d+00;
    ros_Gamma(3) = -2.46541900496934d+00;
    ros_Gamma(4) = -8.05529997906370d-01;

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 3.0_dp

  END SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END SUBROUTINE Rosenbrock_5


SUBROUTINE FunTemplate_5( T, Y, Ydot, LS_NSEL, LS_type )
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Function
 INTEGER,INTENT(IN)::LS_NSEL, LS_type
   REAL(kind=dp) :: T, Y(LS_NSEL)
   REAL(kind=dp) :: Ydot(LS_NSEL)
   REAL(kind=dp) :: Told
   Told = TIME
   TIME = T
   CALL Fun_5( Y, FIX, RCONST, Ydot, LS_NSEL)
   TIME = Told
END SUBROUTINE FunTemplate_5

SUBROUTINE JacTemplate_5( T, Y, Jcb, LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Jacobian
 USE gckpp_LinearAlgebra
 
 INTEGER,INTENT(IN)::LS_NSEL, LS_LU_NONZERO, LS_type
 INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO), LS_LU_ICOL(LS_LU_NONZERO)
    REAL(kind=dp) :: T, Y(LS_NSEL)
#ifdef FULL_ALGEBRA    
    REAL(kind=dp) :: JV(LS_LU_NONZERO), Jcb(LS_NSEL,LS_NSEL)
#else
    REAL(kind=dp) :: Jcb(LS_LU_NONZERO)
#endif   
    REAL(kind=dp) :: Told
#ifdef FULL_ALGEBRA    
    INTEGER :: i, j
#endif   

    Told = TIME
    TIME = T
#ifdef FULL_ALGEBRA        	
	CALL Jac_SP_5(Y, FIX, RCONST, JV,LS_NSEL,LS_LU_NONZERO)	
    DO j=1,LS_NSEL
      DO i=1,LS_NSEL
         Jcb(i,j) = 0.0_dp
      END DO
    END DO
    DO i=1,LS_LU_NONZERO
       Jcb(LS_LU_IROW(i),LS_LU_ICOL(i)) = JV(i)
    END DO
#else   
       CALL Jac_SP_5( Y, FIX, RCONST, Jcb,LS_NSEL,LS_LU_NONZERO)	 
#endif   
    TIME = Told
END SUBROUTINE JacTemplate_5



SUBROUTINE INTEGRATE_6( TIN, TOUT, LS_type,LS_NSEL, LS_NDEL,&
  Prate, Lrate, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )

   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: TIN  ! Start Time
   REAL(kind=dp), INTENT(IN) :: TOUT ! End Time
   
   ! Optional input parameters and statistics
   INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   REAL(kind=dp), INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   REAL(kind=dp), INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U
   INTEGER,		  INTENT(IN) :: LS_type,LS_NSEL,LS_NDEL
   REAL(kind=dp), INTENT(IN) :: Prate(NVAR),Lrate(NVAR)
   
   REAL(kind=dp) :: RCNTRL(20), RSTATUS(20),deltaT
   INTEGER       :: ICNTRL(20), ISTATUS(20), IERR
   REAL(kind=dp) :: VAR_selected(LS_NSEL),VAR_deleted(LS_NDEL),LS_P(LS_NDEL),LS_L(LS_NDEL)
   
   INTEGER, SAVE :: Ntotal = 0

   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.0_dp
   ISTATUS(:) = 0
   RSTATUS(:) = 0.0_dp
   deltaT =  TOUT-TIN
   !~~~> fine-tune the integrator:
   ICNTRL(1) = 0	! 0 - non-autonomous, 1 - autonomous
   ICNTRL(2) = 0	! 0 - vector tolerances, 1 - scalars

   ! If optional parameters are given, and if they are >0, 
   ! then they overwrite default settings. 
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF
   
        VAR_selected=VAR(select_ind_6)
	    VAR_deleted=VAR(delete_ind_6)
		LS_P=Prate(delete_ind_6)
		LS_L=Lrate(delete_ind_6)
        CALL Rosenbrock_6(LU_NSEL_6,VAR_selected,TIN,TOUT,ATOL,RTOL,&
          RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		  LU_NONZERO_6,LU_NSEL_6,LU_CROW_6,LU_DIAG_6,LU_IROW_6,LU_ICOL_6, LS_type)		
		WHERE(LS_L<=(0.01/deltaT))
			VAR_deleted=VAR_deleted+deltaT*(LS_P-LS_L*VAR_deleted)
		ELSEWHERE
			VAR_deleted=LS_P/LS_L+(VAR_deleted-LS_P/LS_L)*EXP(-LS_L*deltaT)
		END WHERE

	    VAR(select_ind_6)=VAR_selected
	    VAR(delete_ind_6)=VAR_deleted
   
   !~~~> Debug option: show no of steps
   ! Ntotal = Ntotal + ISTATUS(Nstp)
   ! PRINT*,'NSTEPS=',ISTATUS(Nstp),' (',Ntotal,')','  O3=', VAR(ind_O3)

   STEPMIN = RSTATUS(Nhexit)
   ! if optional parameters are given for output they 
   ! are updated with the return information
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE INTEGRATE_6

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE Rosenbrock_6(N,Y,Tstart,Tend, &
           AbsTol,RelTol,              &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		   LS_LU_NONZERO,LS_NSEL,LS_LU_CROW,LS_LU_DIAG,LS_LU_IROW,LS_LU_ICOL, LS_type)
  USE gckpp_Parameters
  USE gckpp_LinearAlgebra
  IMPLICIT NONE

   INTEGER,       INTENT(IN)    :: N,LS_LU_NONZERO,LS_NSEL,LS_type   
   INTEGER,INTENT(IN)::LS_LU_CROW(LS_NSEL+1),LS_LU_DIAG(LS_NSEL+1)
   INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO),LS_LU_ICOL(LS_LU_NONZERO)   
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
   REAL(kind=dp), INTENT(IN)    :: Tstart,Tend
   REAL(kind=dp), INTENT(IN)    :: AbsTol(N),RelTol(N)
   INTEGER,       INTENT(IN)    :: ICNTRL(20)
   REAL(kind=dp), INTENT(IN)    :: RCNTRL(20)
   INTEGER,       INTENT(INOUT) :: ISTATUS(20)
   REAL(kind=dp), INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT)   :: IERR
!~~~>  Parameters of the Rosenbrock method, up to 6 stages
   INTEGER ::  ros_S, rosMethod
   INTEGER, PARAMETER :: RS2=1, RS3=2, RS4=3, RD3=4, RD4=5, RG3=6
   REAL(kind=dp) :: ros_A(15), ros_C(15), ros_M(6), ros_E(6), &
                    ros_Alpha(6), ros_Gamma(6), ros_ELO
   LOGICAL :: ros_NewF(6)
   CHARACTER(LEN=12) :: ros_Name
!~~~>  Local variables
   REAL(kind=dp) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   REAL(kind=dp) :: Hmin, Hmax, Hstart
   REAL(kind=dp) :: Texit
   INTEGER       :: i, UplimTol, Max_no_steps
   LOGICAL       :: Autonomous, VectorTol
!~~~>   Parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp

!~~~>  Initialize statistics
   ISTATUS(1:8) = 0
   RSTATUS(1:3) = ZERO

!~~~>  Autonomous or time dependent ODE. Default is time dependent.
   Autonomous = .NOT.(ICNTRL(1) == 0)

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:N) and RelTol(1:N)
   IF (ICNTRL(2) == 0) THEN
      VectorTol = .TRUE.
      UplimTol  = N
   ELSE
      VectorTol = .FALSE.
      UplimTol  = 1
   END IF

!~~~>   Initialize the particular Rosenbrock method selected
   SELECT CASE (ICNTRL(3))
     CASE (1)
       CALL Ros2
     CASE (2)
       CALL Ros3
     CASE (3)
       CALL Ros4
     CASE (0,4)
       CALL Rodas3
     CASE (5)
       CALL Rodas4
     CASE (6)
       CALL Rang3
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(3)=',ICNTRL(3) 
       CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT

!~~~>   The maximum number of steps admitted
   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 200000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Unit roundoff (1+Roundoff>1)
   Roundoff = WLAMCH('E')

!~~~>  Lower bound on the step size: (positive value)
   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Upper bound on the step size: (positive value)
   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Starting step size: (positive value)
   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hold < FacMax
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Check if tolerances are reasonable
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO


!~~~>  CALL Rosenbrock method
   CALL ros_Integrator(Y, Tstart, Tend, Texit,   &
        AbsTol, RelTol,                          &
!  Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!  Error indicator
        IERR)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS !  SUBROUTINES internal to Rosenbrock
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 SUBROUTINE ros_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
   REAL(kind=dp), INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR
   
   IERR = Code
   PRINT * , &
     'Forced exit from Rosenbrock due to the following error:' 
     
   SELECT CASE (Code)
    CASE (-1)
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)    
      PRINT * , '--> Selected Rosenbrock method not implemented'
    CASE (-3)    
      PRINT * , '--> Hmin/Hmax/Hstart must be positive'
    CASE (-4)    
      PRINT * , '--> FacMin/FacMax/FacRej must be positive'
    CASE (-5) 
      PRINT * , '--> Improper tolerance values'
    CASE (-6) 
      PRINT * , '--> No of steps exceeds maximum bound'
    CASE (-7) 
      PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
    CASE (-8)    
      PRINT * , '--> Matrix is repeatedly singular'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT
   
   PRINT *, "T=", T, "and H=", H
     
 END SUBROUTINE ros_ErrorMsg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_Integrator (Y, Tstart, Tend, T,  &
        AbsTol, RelTol,                          &
!~~~> Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic Rosenbrock method
!      defined by ros_S (no of stages)
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

!~~~> Input: the initial condition at Tstart; Output: the solution at T
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
!~~~> Input: integration interval
   REAL(kind=dp), INTENT(IN) :: Tstart,Tend
!~~~> Output: time at which the solution is returned (T=Tend if success)
   REAL(kind=dp), INTENT(OUT) ::  T
!~~~> Input: tolerances
   REAL(kind=dp), INTENT(IN) ::  AbsTol(N), RelTol(N)
!~~~> Input: integration parameters
   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   REAL(kind=dp), INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   REAL(kind=dp), INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables
   REAL(kind=dp) :: Ynew(N), Fcn0(N), Fcn(N)
   REAL(kind=dp) :: K(N*ros_S), dFdT(N)
#ifdef FULL_ALGEBRA    
   REAL(kind=dp) :: Jac0(N,N), Ghimj(N,N)
#else
   REAL(kind=dp) :: Jac0(LS_LU_NONZERO), Ghimj(LS_LU_NONZERO)
#endif
   REAL(kind=dp) :: H, Hnew, HC, HG, Fac, Tau
   REAL(kind=dp) :: Err, Yerr(N)
   INTEGER :: Pivot(N), Direction, ioffset, j, istage
   LOGICAL :: RejectLastH, RejectMoreH, Singular
!~~~>  Local parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp
!~~~>  Locally called functions
!    REAL(kind=dp) WLAMCH
!    EXTERNAL WLAMCH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~>  Initial preparations
   T = Tstart
   RSTATUS(Nhexit) = ZERO
   H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   IF (ABS(H) <= 10.0_dp*Roundoff) H = DeltaMin

   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF
   H = Direction*H

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.

!~~~> Time loop begins below

TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) )

   IF ( ISTATUS(Nstp) > Max_no_steps ) THEN  ! Too many steps
      CALL ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
      CALL ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF

!~~~>  Limit H if necessary to avoid going beyond Tend
   H = MIN(H,ABS(Tend-T))

!~~~>   Compute the function at current time
   CALL FunTemplate_6(T,Y,Fcn0, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1

!~~~>  Compute the function derivative with respect to T
   IF (.NOT.Autonomous) THEN
      CALL ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
   END IF

!~~~>   Compute the Jacobian at current time
   CALL JacTemplate_6(T,Y,Jac0,LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)   
   ISTATUS(Njac) = ISTATUS(Njac) + 1

!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO

   CALL ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular)
   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       CALL ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S

      ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
       ioffset = N*(istage-1)

      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         !slim: CALL WCOPY(N,Fcn0,1,Fcn,1)
	 Fcn(1:N) = Fcn0(1:N)
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
         !slim: CALL WCOPY(N,Y,1,Ynew,1)
	 Ynew(1:N) = Y(1:N)
         DO j = 1, istage-1
           CALL WAXPY(N,ros_A((istage-1)*(istage-2)/2+j), &
            K(N*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL FunTemplate_6(Tau,Ynew,Fcn, LS_NSEL, LS_type)
         ISTATUS(Nfun) = ISTATUS(Nfun) + 1
       END IF ! if istage == 1 elseif ros_NewF(istage)
       !slim: CALL WCOPY(N,Fcn,1,K(ioffset+1),1)
       K(ioffset+1:ioffset+N) = Fcn(1:N)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL WAXPY(N,HC,K(N*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL WAXPY(N,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL ros_Solve(Ghimj, Pivot, K(ioffset+1))

   END DO Stage


!~~~>  Compute the new solution
   !slim: CALL WCOPY(N,Y,1,Ynew,1)
   Ynew(1:N) = Y(1:N)
   DO j=1,ros_S
         CALL WAXPY(N,ros_M(j),K(N*(j-1)+1),1,Ynew,1)
   END DO

!~~~>  Compute the error estimation
   !slim: CALL WSCAL(N,ZERO,Yerr,1)
   Yerr(1:N) = ZERO
   DO j=1,ros_S
        CALL WAXPY(N,ros_E(j),K(N*(j-1)+1),1,Yerr,1)
   END DO
   Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac

!~~~>  Check the error magnitude and adjust step size
   ISTATUS(Nstp) = ISTATUS(Nstp) + 1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
      ISTATUS(Nacc) = ISTATUS(Nacc) + 1
      !slim: CALL WCOPY(N,Ynew,1,Y,1)
      Y(1:N) = Ynew(1:N)
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H)
      END IF
      RSTATUS(Nhexit) = H
      RSTATUS(Nhnew)  = Hnew
      RSTATUS(Ntexit) = T
      RejectLastH = .FALSE.
      RejectMoreH = .FALSE.
      H = Hnew
      EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
   ELSE           !~~~> Reject step
      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (ISTATUS(Nacc) >= 1)  ISTATUS(Nrej) = ISTATUS(Nrej) + 1
   END IF ! Err <= 1

   END DO UntilAccepted

   
   END DO TimeLoop

!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful

  END SUBROUTINE ros_Integrator


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  REAL(kind=dp) FUNCTION ros_ErrorNorm ( Y, Ynew, Yerr, &
                               AbsTol, RelTol, VectorTol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

! Input arguments
   REAL(kind=dp), INTENT(IN) :: Y(N), Ynew(N), &
          Yerr(N), AbsTol(N), RelTol(N)
   LOGICAL, INTENT(IN) ::  VectorTol
! Local variables
   REAL(kind=dp) :: Err, Scale, Ymax
   INTEGER  :: i
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp

   Err = ZERO
   DO i=1,N
     Ymax = MAX(ABS(Y(i)),ABS(Ynew(i)))
     IF (VectorTol) THEN
       Scale = AbsTol(i)+RelTol(i)*Ymax
     ELSE
       Scale = AbsTol(1)+RelTol(1)*Ymax
     END IF
     Err = Err+(Yerr(i)/Scale)**2
   END DO
   Err  = SQRT(Err/N)

   ros_ErrorNorm = MAX(Err,1.0d-10)

  END FUNCTION ros_ErrorNorm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> The time partial derivative of the function by finite differences
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

!~~~> Input arguments
   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(N), Fcn0(N)
!~~~> Output arguments
   REAL(kind=dp), INTENT(OUT) :: dFdT(N)
!~~~> Local variables
   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL FunTemplate_6(T+Delta,Y,dFdT, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1
   CALL WAXPY(N,(-ONE),Fcn0,1,dFdT,1)
   CALL WSCAL(N,(ONE/Delta),dFdT,1)

  END SUBROUTINE ros_FunTimeDerivative


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_PrepareMatrix ( H, Direction, gam, &
             Jac0, Ghimj, Pivot, Singular )
! --- --- --- --- --- --- --- --- --- --- --- --- ---
!  Prepares the LHS matrix for stage calculations
!  1.  Construct Ghimj = 1/(H*ham) - Jac0
!      "(Gamma H) Inverse Minus Jacobian"
!  2.  Repeat LU decomposition of Ghimj until successful.
!       -half the step size if LU decomposition fails and retry
!       -exit after 5 consecutive fails
! --- --- --- --- --- --- --- --- --- --- --- --- ---
   IMPLICIT NONE

!~~~> Input arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) ::  Jac0(N,N)
#else
   REAL(kind=dp), INTENT(IN) ::  Jac0(LS_LU_NONZERO)
#endif   
   REAL(kind=dp), INTENT(IN) ::  gam
   INTEGER, INTENT(IN) ::  Direction
!~~~> Output arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(OUT) :: Ghimj(N,N)
#else
   REAL(kind=dp), INTENT(OUT) :: Ghimj(LS_LU_NONZERO)
#endif   
   LOGICAL, INTENT(OUT) ::  Singular
   INTEGER, INTENT(OUT) ::  Pivot(N)
!~~~> Inout arguments
   REAL(kind=dp), INTENT(INOUT) :: H   ! step size is decreased when LU fails
!~~~> Local variables
   INTEGER  :: i, ISING, Nconsecutive
   REAL(kind=dp) :: ghinv
   REAL(kind=dp), PARAMETER :: ONE  = 1.0_dp, HALF = 0.5_dp

   Nconsecutive = 0
   Singular = .TRUE.

   DO WHILE (Singular)

!~~~>    Construct Ghimj = 1/(H*gam) - Jac0
#ifdef FULL_ALGEBRA    
     !slim: CALL WCOPY(N*N,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(N*N,(-ONE),Ghimj,1)
     Ghimj = -Jac0
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(i,i) = Ghimj(i,i)+ghinv
     END DO
#else
     !slim: CALL WCOPY(LS_LU_NONZERO,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(LS_LU_NONZERO,(-ONE),Ghimj,1)
     Ghimj(1:LS_LU_NONZERO) = -Jac0(1:LS_LU_NONZERO)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(LS_LU_DIAG(i)) = Ghimj(LS_LU_DIAG(i))+ghinv
     END DO
#endif   
!~~~>    Compute LU decomposition
     CALL ros_Decomp( Ghimj, Pivot, ISING )
     IF (ISING == 0) THEN
!~~~>    If successful done
        Singular = .FALSE.
     ELSE ! ISING .ne. 0
!~~~>    If unsuccessful half the step size; if 5 consecutive fails then return
        ISTATUS(Nsng) = ISTATUS(Nsng) + 1
        Nconsecutive = Nconsecutive+1
        Singular = .TRUE.
        PRINT*,'Warning: LU Decomposition returned ISING = ',ISING
        IF (Nconsecutive <= 5) THEN ! Less than 5 consecutive failed decompositions
           H = H*HALF
        ELSE  ! More than 5 consecutive failed decompositions
           RETURN
        END IF  ! Nconsecutive
      END IF    ! ISING

   END DO ! WHILE Singular

  END SUBROUTINE ros_PrepareMatrix


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Decomp( A, Pivot, ISING )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the LU decomposition
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Inout variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(INOUT) :: A(N,N)
#else   
   REAL(kind=dp), INTENT(INOUT) :: A(LS_LU_NONZERO)
#endif
!~~~> Output variables
   INTEGER, INTENT(OUT) :: Pivot(N), ISING

#ifdef FULL_ALGEBRA    
   CALL  DGETRF( N, N, A, N, Pivot, ISING )
#else   
   CALL KppDecomp ( A, ISING,LS_NSEL,LS_LU_NONZERO,LS_LU_CROW,LS_LU_DIAG,LS_LU_ICOL)   
   Pivot(1) = 1
#endif
   ISTATUS(Ndec) = ISTATUS(Ndec) + 1

  END SUBROUTINE ros_Decomp


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Solve( A, Pivot, b)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the forward/backward substitution (using pre-computed LU decomposition)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Input variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) :: A(N,N)
   INTEGER :: ISING
#else   
   REAL(kind=dp), INTENT(IN) :: A(LS_LU_NONZERO)
#endif
   INTEGER, INTENT(IN) :: Pivot(N)
!~~~> InOut variables
   REAL(kind=dp), INTENT(INOUT) :: b(N)

#ifdef FULL_ALGEBRA    
   CALL  DGETRS( 'N', N , 1, A, N, Pivot, b, N, ISING )
   IF ( Info < 0 ) THEN
      PRINT*,"Error in DGETRS. ISING=",ISING
   END IF  
#else  
   	CALL KppSolve_6(LS_LU_NONZERO,LS_NSEL, A, b,LS_type )
#endif

   ISTATUS(Nsol) = ISTATUS(Nsol) + 1

  END SUBROUTINE ros_Solve



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   DOUBLE PRECISION g

    g = 1.0_dp + 1.0_dp/SQRT(2.0_dp)
    rosMethod = RS2
!~~~> Name of the method
    ros_Name = 'ROS-2'
!~~~> Number of stages
    ros_S = 2

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = (1.0_dp)/g
    ros_C(1) = (-2.0_dp)/g
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
!~~~> M_i = Coefficients for new step solution
    ros_M(1)= (3.0_dp)/(2.0_dp*g)
    ros_M(2)= (1.0_dp)/(2.0_dp*g)
! E_i = Coefficients for error estimator
    ros_E(1) = 1.0_dp/(2.0_dp*g)
    ros_E(2) = 1.0_dp/(2.0_dp*g)
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus one
    ros_ELO = 2.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.0_dp
    ros_Alpha(2) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = g
    ros_Gamma(2) =-g

 END SUBROUTINE Ros2


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   rosMethod = RS3
!~~~> Name of the method
   ros_Name = 'ROS-3'
!~~~> Number of stages
   ros_S = 3

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1)= 1.0_dp
   ros_A(2)= 1.0_dp
   ros_A(3)= 0.0_dp

   ros_C(1) = -0.10156171083877702091975600115545E+01_dp
   ros_C(2) =  0.40759956452537699824805835358067E+01_dp
   ros_C(3) =  0.92076794298330791242156818474003E+01_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) =  0.1E+01_dp
   ros_M(2) =  0.61697947043828245592553615689730E+01_dp
   ros_M(3) = -0.42772256543218573326238373806514_dp
! E_i = Coefficients for error estimator
   ros_E(1) =  0.5_dp
   ros_E(2) = -0.29079558716805469821718236208017E+01_dp
   ros_E(3) =  0.22354069897811569627360909276199_dp
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1)= 0.0_dp
   ros_Alpha(2)= 0.43586652150845899941601945119356_dp
   ros_Alpha(3)= 0.43586652150845899941601945119356_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1)= 0.43586652150845899941601945119356_dp
   ros_Gamma(2)= 0.24291996454816804366592249683314_dp
   ros_Gamma(3)= 0.21851380027664058511513169485832E+01_dp

  END SUBROUTINE Ros3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
!     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1990)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RS4
!~~~> Name of the method
   ros_Name = 'ROS-4'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.2000000000000000E+01_dp
   ros_A(2) = 0.1867943637803922E+01_dp
   ros_A(3) = 0.2344449711399156_dp
   ros_A(4) = ros_A(2)
   ros_A(5) = ros_A(3)
   ros_A(6) = 0.0_dp

   ros_C(1) =-0.7137615036412310E+01_dp
   ros_C(2) = 0.2580708087951457E+01_dp
   ros_C(3) = 0.6515950076447975_dp
   ros_C(4) =-0.2137148994382534E+01_dp
   ros_C(5) =-0.3214669691237626_dp
   ros_C(6) =-0.6949742501781779_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .TRUE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 0.2255570073418735E+01_dp
   ros_M(2) = 0.2870493262186792_dp
   ros_M(3) = 0.4353179431840180_dp
   ros_M(4) = 0.1093502252409163E+01_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) =-0.2815431932141155_dp
   ros_E(2) =-0.7276199124938920E-01_dp
   ros_E(3) =-0.1082196201495311_dp
   ros_E(4) =-0.1093502252409163E+01_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 4.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.1145640000000000E+01_dp
   ros_Alpha(3) = 0.6552168638155900_dp
   ros_Alpha(4) = ros_Alpha(3)
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5728200000000000_dp
   ros_Gamma(2) =-0.1769193891319233E+01_dp
   ros_Gamma(3) = 0.7592633437920482_dp
   ros_Gamma(4) =-0.1049021087100450_dp

  END SUBROUTINE Ros4

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- A STIFFLY-STABLE METHOD, 4 stages, order 3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RD3
!~~~> Name of the method
   ros_Name = 'RODAS-3'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.0_dp
   ros_A(2) = 2.0_dp
   ros_A(3) = 0.0_dp
   ros_A(4) = 2.0_dp
   ros_A(5) = 0.0_dp
   ros_A(6) = 1.0_dp

   ros_C(1) = 4.0_dp
   ros_C(2) = 1.0_dp
   ros_C(3) =-1.0_dp
   ros_C(4) = 1.0_dp
   ros_C(5) =-1.0_dp
   ros_C(6) =-(8.0_dp/3.0_dp)

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .FALSE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .TRUE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 2.0_dp
   ros_M(2) = 0.0_dp
   ros_M(3) = 1.0_dp
   ros_M(4) = 1.0_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) = 0.0_dp
   ros_E(2) = 0.0_dp
   ros_E(3) = 0.0_dp
   ros_E(4) = 1.0_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.0_dp
   ros_Alpha(3) = 1.0_dp
   ros_Alpha(4) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5_dp
   ros_Gamma(2) = 1.5_dp
   ros_Gamma(3) = 0.0_dp
   ros_Gamma(4) = 0.0_dp

  END SUBROUTINE Rodas3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1996)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RD4
!~~~> Name of the method
    ros_Name = 'RODAS-4'
!~~~> Number of stages
    ros_S = 6

!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.000_dp
    ros_Alpha(2) = 0.386_dp
    ros_Alpha(3) = 0.210_dp
    ros_Alpha(4) = 0.630_dp
    ros_Alpha(5) = 1.000_dp
    ros_Alpha(6) = 1.000_dp

!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = 0.2500000000000000_dp
    ros_Gamma(2) =-0.1043000000000000_dp
    ros_Gamma(3) = 0.1035000000000000_dp
    ros_Gamma(4) =-0.3620000000000023E-01_dp
    ros_Gamma(5) = 0.0_dp
    ros_Gamma(6) = 0.0_dp

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = 0.1544000000000000E+01_dp
    ros_A(2) = 0.9466785280815826_dp
    ros_A(3) = 0.2557011698983284_dp
    ros_A(4) = 0.3314825187068521E+01_dp
    ros_A(5) = 0.2896124015972201E+01_dp
    ros_A(6) = 0.9986419139977817_dp
    ros_A(7) = 0.1221224509226641E+01_dp
    ros_A(8) = 0.6019134481288629E+01_dp
    ros_A(9) = 0.1253708332932087E+02_dp
    ros_A(10) =-0.6878860361058950_dp
    ros_A(11) = ros_A(7)
    ros_A(12) = ros_A(8)
    ros_A(13) = ros_A(9)
    ros_A(14) = ros_A(10)
    ros_A(15) = 1.0_dp

    ros_C(1) =-0.5668800000000000E+01_dp
    ros_C(2) =-0.2430093356833875E+01_dp
    ros_C(3) =-0.2063599157091915_dp
    ros_C(4) =-0.1073529058151375_dp
    ros_C(5) =-0.9594562251023355E+01_dp
    ros_C(6) =-0.2047028614809616E+02_dp
    ros_C(7) = 0.7496443313967647E+01_dp
    ros_C(8) =-0.1024680431464352E+02_dp
    ros_C(9) =-0.3399990352819905E+02_dp
    ros_C(10) = 0.1170890893206160E+02_dp
    ros_C(11) = 0.8083246795921522E+01_dp
    ros_C(12) =-0.7981132988064893E+01_dp
    ros_C(13) =-0.3152159432874371E+02_dp
    ros_C(14) = 0.1631930543123136E+02_dp
    ros_C(15) =-0.6058818238834054E+01_dp

!~~~> M_i = Coefficients for new step solution
    ros_M(1) = ros_A(7)
    ros_M(2) = ros_A(8)
    ros_M(3) = ros_A(9)
    ros_M(4) = ros_A(10)
    ros_M(5) = 1.0_dp
    ros_M(6) = 1.0_dp

!~~~> E_i  = Coefficients for error estimator
    ros_E(1) = 0.0_dp
    ros_E(2) = 0.0_dp
    ros_E(3) = 0.0_dp
    ros_E(4) = 0.0_dp
    ros_E(5) = 0.0_dp
    ros_E(6) = 1.0_dp

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.
    ros_NewF(5) = .TRUE.
    ros_NewF(6) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 4.0_dp

  END SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! STIFFLY-STABLE W METHOD OF ORDER 3, WITH 4 STAGES
!
! J. RANG and L. ANGERMANN
! NEW ROSENBROCK W-METHODS OF ORDER 3
! FOR PARTIAL DIFFERENTIAL ALGEBRAIC
!        EQUATIONS OF INDEX 1                                             
! BIT Numerical Mathematics (2005) 45: 761-787
!  DOI: 10.1007/s10543-005-0035-y
! Table 4.1-4.2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RG3
!~~~> Name of the method
    ros_Name = 'RANG-3'
!~~~> Number of stages
    ros_S = 4

    ros_A(1) = 5.09052051067020d+00;
    ros_A(2) = 5.09052051067020d+00;
    ros_A(3) = 0.0d0;
    ros_A(4) = 4.97628111010787d+00;
    ros_A(5) = 2.77268164715849d-02;
    ros_A(6) = 2.29428036027904d-01;

    ros_C(1) = -1.16790812312283d+01;
    ros_C(2) = -1.64057326467367d+01;
    ros_C(3) = -2.77268164715850d-01;
    ros_C(4) = -8.38103960500476d+00;
    ros_C(5) = -8.48328409199343d-01;
    ros_C(6) =  2.87009860433106d-01;

    ros_M(1) =  5.22582761233094d+00;
    ros_M(2) = -5.56971148154165d-01;
    ros_M(3) =  3.57979469353645d-01;
    ros_M(4) =  1.72337398521064d+00;

    ros_E(1) = -5.16845212784040d+00;
    ros_E(2) = -1.26351942603842d+00;
    ros_E(3) = -1.11022302462516d-16;
    ros_E(4) =  2.22044604925031d-16;

    ros_Alpha(1) = 0.0d00;
    ros_Alpha(2) = 2.21878746765329d+00;
    ros_Alpha(3) = 2.21878746765329d+00;
    ros_Alpha(4) = 1.55392337535788d+00;

    ros_Gamma(1) =  4.35866521508459d-01;
    ros_Gamma(2) = -1.78292094614483d+00;
    ros_Gamma(3) = -2.46541900496934d+00;
    ros_Gamma(4) = -8.05529997906370d-01;

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 3.0_dp

  END SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END SUBROUTINE Rosenbrock_6


SUBROUTINE FunTemplate_6( T, Y, Ydot, LS_NSEL, LS_type )
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Function
 INTEGER,INTENT(IN)::LS_NSEL, LS_type
   REAL(kind=dp) :: T, Y(LS_NSEL)
   REAL(kind=dp) :: Ydot(LS_NSEL)
   REAL(kind=dp) :: Told
   Told = TIME
   TIME = T
   CALL Fun_6( Y, FIX, RCONST, Ydot, LS_NSEL)
   TIME = Told
END SUBROUTINE FunTemplate_6

SUBROUTINE JacTemplate_6( T, Y, Jcb, LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Jacobian
 USE gckpp_LinearAlgebra
 
 INTEGER,INTENT(IN)::LS_NSEL, LS_LU_NONZERO, LS_type
 INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO), LS_LU_ICOL(LS_LU_NONZERO)
    REAL(kind=dp) :: T, Y(LS_NSEL)
#ifdef FULL_ALGEBRA    
    REAL(kind=dp) :: JV(LS_LU_NONZERO), Jcb(LS_NSEL,LS_NSEL)
#else
    REAL(kind=dp) :: Jcb(LS_LU_NONZERO)
#endif   
    REAL(kind=dp) :: Told
#ifdef FULL_ALGEBRA    
    INTEGER :: i, j
#endif   

    Told = TIME
    TIME = T
#ifdef FULL_ALGEBRA        	
	CALL Jac_SP_6(Y, FIX, RCONST, JV,LS_NSEL,LS_LU_NONZERO)	
    DO j=1,LS_NSEL
      DO i=1,LS_NSEL
         Jcb(i,j) = 0.0_dp
      END DO
    END DO
    DO i=1,LS_LU_NONZERO
       Jcb(LS_LU_IROW(i),LS_LU_ICOL(i)) = JV(i)
    END DO
#else   
       CALL Jac_SP_6( Y, FIX, RCONST, Jcb,LS_NSEL,LS_LU_NONZERO)	 
#endif   
    TIME = Told
END SUBROUTINE JacTemplate_6



SUBROUTINE INTEGRATE_7( TIN, TOUT, LS_type,LS_NSEL, LS_NDEL,&
  Prate, Lrate, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )

   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: TIN  ! Start Time
   REAL(kind=dp), INTENT(IN) :: TOUT ! End Time
   
   ! Optional input parameters and statistics
   INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   REAL(kind=dp), INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   REAL(kind=dp), INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U
   INTEGER,		  INTENT(IN) :: LS_type,LS_NSEL,LS_NDEL
   REAL(kind=dp), INTENT(IN) :: Prate(NVAR),Lrate(NVAR)
   
   REAL(kind=dp) :: RCNTRL(20), RSTATUS(20),deltaT
   INTEGER       :: ICNTRL(20), ISTATUS(20), IERR
   REAL(kind=dp) :: VAR_selected(LS_NSEL),VAR_deleted(LS_NDEL),LS_P(LS_NDEL),LS_L(LS_NDEL)
   
   INTEGER, SAVE :: Ntotal = 0

   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.0_dp
   ISTATUS(:) = 0
   RSTATUS(:) = 0.0_dp
   deltaT =  TOUT-TIN
   !~~~> fine-tune the integrator:
   ICNTRL(1) = 0	! 0 - non-autonomous, 1 - autonomous
   ICNTRL(2) = 0	! 0 - vector tolerances, 1 - scalars

   ! If optional parameters are given, and if they are >0, 
   ! then they overwrite default settings. 
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF
   
        VAR_selected=VAR(select_ind_7)
	    VAR_deleted=VAR(delete_ind_7)
		LS_P=Prate(delete_ind_7)
		LS_L=Lrate(delete_ind_7)
        CALL Rosenbrock_7(LU_NSEL_7,VAR_selected,TIN,TOUT,ATOL,RTOL,&
          RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		  LU_NONZERO_7,LU_NSEL_7,LU_CROW_7,LU_DIAG_7,LU_IROW_7,LU_ICOL_7, LS_type)		
		WHERE(LS_L<=(0.01/deltaT))
			VAR_deleted=VAR_deleted+deltaT*(LS_P-LS_L*VAR_deleted)
		ELSEWHERE
			VAR_deleted=LS_P/LS_L+(VAR_deleted-LS_P/LS_L)*EXP(-LS_L*deltaT)
		END WHERE

	    VAR(select_ind_7)=VAR_selected
	    VAR(delete_ind_7)=VAR_deleted
   
   !~~~> Debug option: show no of steps
   ! Ntotal = Ntotal + ISTATUS(Nstp)
   ! PRINT*,'NSTEPS=',ISTATUS(Nstp),' (',Ntotal,')','  O3=', VAR(ind_O3)

   STEPMIN = RSTATUS(Nhexit)
   ! if optional parameters are given for output they 
   ! are updated with the return information
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE INTEGRATE_7

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE Rosenbrock_7(N,Y,Tstart,Tend, &
           AbsTol,RelTol,              &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		   LS_LU_NONZERO,LS_NSEL,LS_LU_CROW,LS_LU_DIAG,LS_LU_IROW,LS_LU_ICOL, LS_type)
  USE gckpp_Parameters
  USE gckpp_LinearAlgebra
  IMPLICIT NONE

   INTEGER,       INTENT(IN)    :: N,LS_LU_NONZERO,LS_NSEL,LS_type   
   INTEGER,INTENT(IN)::LS_LU_CROW(LS_NSEL+1),LS_LU_DIAG(LS_NSEL+1)
   INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO),LS_LU_ICOL(LS_LU_NONZERO)   
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
   REAL(kind=dp), INTENT(IN)    :: Tstart,Tend
   REAL(kind=dp), INTENT(IN)    :: AbsTol(N),RelTol(N)
   INTEGER,       INTENT(IN)    :: ICNTRL(20)
   REAL(kind=dp), INTENT(IN)    :: RCNTRL(20)
   INTEGER,       INTENT(INOUT) :: ISTATUS(20)
   REAL(kind=dp), INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT)   :: IERR
!~~~>  Parameters of the Rosenbrock method, up to 6 stages
   INTEGER ::  ros_S, rosMethod
   INTEGER, PARAMETER :: RS2=1, RS3=2, RS4=3, RD3=4, RD4=5, RG3=6
   REAL(kind=dp) :: ros_A(15), ros_C(15), ros_M(6), ros_E(6), &
                    ros_Alpha(6), ros_Gamma(6), ros_ELO
   LOGICAL :: ros_NewF(6)
   CHARACTER(LEN=12) :: ros_Name
!~~~>  Local variables
   REAL(kind=dp) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   REAL(kind=dp) :: Hmin, Hmax, Hstart
   REAL(kind=dp) :: Texit
   INTEGER       :: i, UplimTol, Max_no_steps
   LOGICAL       :: Autonomous, VectorTol
!~~~>   Parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp

!~~~>  Initialize statistics
   ISTATUS(1:8) = 0
   RSTATUS(1:3) = ZERO

!~~~>  Autonomous or time dependent ODE. Default is time dependent.
   Autonomous = .NOT.(ICNTRL(1) == 0)

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:N) and RelTol(1:N)
   IF (ICNTRL(2) == 0) THEN
      VectorTol = .TRUE.
      UplimTol  = N
   ELSE
      VectorTol = .FALSE.
      UplimTol  = 1
   END IF

!~~~>   Initialize the particular Rosenbrock method selected
   SELECT CASE (ICNTRL(3))
     CASE (1)
       CALL Ros2
     CASE (2)
       CALL Ros3
     CASE (3)
       CALL Ros4
     CASE (0,4)
       CALL Rodas3
     CASE (5)
       CALL Rodas4
     CASE (6)
       CALL Rang3
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(3)=',ICNTRL(3) 
       CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT

!~~~>   The maximum number of steps admitted
   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 200000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Unit roundoff (1+Roundoff>1)
   Roundoff = WLAMCH('E')

!~~~>  Lower bound on the step size: (positive value)
   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Upper bound on the step size: (positive value)
   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Starting step size: (positive value)
   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hold < FacMax
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Check if tolerances are reasonable
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO


!~~~>  CALL Rosenbrock method
   CALL ros_Integrator(Y, Tstart, Tend, Texit,   &
        AbsTol, RelTol,                          &
!  Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!  Error indicator
        IERR)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS !  SUBROUTINES internal to Rosenbrock
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 SUBROUTINE ros_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
   REAL(kind=dp), INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR
   
   IERR = Code
   PRINT * , &
     'Forced exit from Rosenbrock due to the following error:' 
     
   SELECT CASE (Code)
    CASE (-1)
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)    
      PRINT * , '--> Selected Rosenbrock method not implemented'
    CASE (-3)    
      PRINT * , '--> Hmin/Hmax/Hstart must be positive'
    CASE (-4)    
      PRINT * , '--> FacMin/FacMax/FacRej must be positive'
    CASE (-5) 
      PRINT * , '--> Improper tolerance values'
    CASE (-6) 
      PRINT * , '--> No of steps exceeds maximum bound'
    CASE (-7) 
      PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
    CASE (-8)    
      PRINT * , '--> Matrix is repeatedly singular'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT
   
   PRINT *, "T=", T, "and H=", H
     
 END SUBROUTINE ros_ErrorMsg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_Integrator (Y, Tstart, Tend, T,  &
        AbsTol, RelTol,                          &
!~~~> Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic Rosenbrock method
!      defined by ros_S (no of stages)
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

!~~~> Input: the initial condition at Tstart; Output: the solution at T
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
!~~~> Input: integration interval
   REAL(kind=dp), INTENT(IN) :: Tstart,Tend
!~~~> Output: time at which the solution is returned (T=Tend if success)
   REAL(kind=dp), INTENT(OUT) ::  T
!~~~> Input: tolerances
   REAL(kind=dp), INTENT(IN) ::  AbsTol(N), RelTol(N)
!~~~> Input: integration parameters
   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   REAL(kind=dp), INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   REAL(kind=dp), INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables
   REAL(kind=dp) :: Ynew(N), Fcn0(N), Fcn(N)
   REAL(kind=dp) :: K(N*ros_S), dFdT(N)
#ifdef FULL_ALGEBRA    
   REAL(kind=dp) :: Jac0(N,N), Ghimj(N,N)
#else
   REAL(kind=dp) :: Jac0(LS_LU_NONZERO), Ghimj(LS_LU_NONZERO)
#endif
   REAL(kind=dp) :: H, Hnew, HC, HG, Fac, Tau
   REAL(kind=dp) :: Err, Yerr(N)
   INTEGER :: Pivot(N), Direction, ioffset, j, istage
   LOGICAL :: RejectLastH, RejectMoreH, Singular
!~~~>  Local parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp
!~~~>  Locally called functions
!    REAL(kind=dp) WLAMCH
!    EXTERNAL WLAMCH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~>  Initial preparations
   T = Tstart
   RSTATUS(Nhexit) = ZERO
   H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   IF (ABS(H) <= 10.0_dp*Roundoff) H = DeltaMin

   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF
   H = Direction*H

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.

!~~~> Time loop begins below

TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) )

   IF ( ISTATUS(Nstp) > Max_no_steps ) THEN  ! Too many steps
      CALL ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
      CALL ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF

!~~~>  Limit H if necessary to avoid going beyond Tend
   H = MIN(H,ABS(Tend-T))

!~~~>   Compute the function at current time
   CALL FunTemplate_7(T,Y,Fcn0, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1

!~~~>  Compute the function derivative with respect to T
   IF (.NOT.Autonomous) THEN
      CALL ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
   END IF

!~~~>   Compute the Jacobian at current time
   CALL JacTemplate_7(T,Y,Jac0,LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)   
   ISTATUS(Njac) = ISTATUS(Njac) + 1

!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO

   CALL ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular)
   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       CALL ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S

      ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
       ioffset = N*(istage-1)

      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         !slim: CALL WCOPY(N,Fcn0,1,Fcn,1)
	 Fcn(1:N) = Fcn0(1:N)
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
         !slim: CALL WCOPY(N,Y,1,Ynew,1)
	 Ynew(1:N) = Y(1:N)
         DO j = 1, istage-1
           CALL WAXPY(N,ros_A((istage-1)*(istage-2)/2+j), &
            K(N*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL FunTemplate_7(Tau,Ynew,Fcn, LS_NSEL, LS_type)
         ISTATUS(Nfun) = ISTATUS(Nfun) + 1
       END IF ! if istage == 1 elseif ros_NewF(istage)
       !slim: CALL WCOPY(N,Fcn,1,K(ioffset+1),1)
       K(ioffset+1:ioffset+N) = Fcn(1:N)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL WAXPY(N,HC,K(N*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL WAXPY(N,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL ros_Solve(Ghimj, Pivot, K(ioffset+1))

   END DO Stage


!~~~>  Compute the new solution
   !slim: CALL WCOPY(N,Y,1,Ynew,1)
   Ynew(1:N) = Y(1:N)
   DO j=1,ros_S
         CALL WAXPY(N,ros_M(j),K(N*(j-1)+1),1,Ynew,1)
   END DO

!~~~>  Compute the error estimation
   !slim: CALL WSCAL(N,ZERO,Yerr,1)
   Yerr(1:N) = ZERO
   DO j=1,ros_S
        CALL WAXPY(N,ros_E(j),K(N*(j-1)+1),1,Yerr,1)
   END DO
   Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac

!~~~>  Check the error magnitude and adjust step size
   ISTATUS(Nstp) = ISTATUS(Nstp) + 1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
      ISTATUS(Nacc) = ISTATUS(Nacc) + 1
      !slim: CALL WCOPY(N,Ynew,1,Y,1)
      Y(1:N) = Ynew(1:N)
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H)
      END IF
      RSTATUS(Nhexit) = H
      RSTATUS(Nhnew)  = Hnew
      RSTATUS(Ntexit) = T
      RejectLastH = .FALSE.
      RejectMoreH = .FALSE.
      H = Hnew
      EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
   ELSE           !~~~> Reject step
      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (ISTATUS(Nacc) >= 1)  ISTATUS(Nrej) = ISTATUS(Nrej) + 1
   END IF ! Err <= 1

   END DO UntilAccepted

   
   END DO TimeLoop

!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful

  END SUBROUTINE ros_Integrator


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  REAL(kind=dp) FUNCTION ros_ErrorNorm ( Y, Ynew, Yerr, &
                               AbsTol, RelTol, VectorTol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

! Input arguments
   REAL(kind=dp), INTENT(IN) :: Y(N), Ynew(N), &
          Yerr(N), AbsTol(N), RelTol(N)
   LOGICAL, INTENT(IN) ::  VectorTol
! Local variables
   REAL(kind=dp) :: Err, Scale, Ymax
   INTEGER  :: i
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp

   Err = ZERO
   DO i=1,N
     Ymax = MAX(ABS(Y(i)),ABS(Ynew(i)))
     IF (VectorTol) THEN
       Scale = AbsTol(i)+RelTol(i)*Ymax
     ELSE
       Scale = AbsTol(1)+RelTol(1)*Ymax
     END IF
     Err = Err+(Yerr(i)/Scale)**2
   END DO
   Err  = SQRT(Err/N)

   ros_ErrorNorm = MAX(Err,1.0d-10)

  END FUNCTION ros_ErrorNorm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> The time partial derivative of the function by finite differences
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

!~~~> Input arguments
   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(N), Fcn0(N)
!~~~> Output arguments
   REAL(kind=dp), INTENT(OUT) :: dFdT(N)
!~~~> Local variables
   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL FunTemplate_7(T+Delta,Y,dFdT, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1
   CALL WAXPY(N,(-ONE),Fcn0,1,dFdT,1)
   CALL WSCAL(N,(ONE/Delta),dFdT,1)

  END SUBROUTINE ros_FunTimeDerivative


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_PrepareMatrix ( H, Direction, gam, &
             Jac0, Ghimj, Pivot, Singular )
! --- --- --- --- --- --- --- --- --- --- --- --- ---
!  Prepares the LHS matrix for stage calculations
!  1.  Construct Ghimj = 1/(H*ham) - Jac0
!      "(Gamma H) Inverse Minus Jacobian"
!  2.  Repeat LU decomposition of Ghimj until successful.
!       -half the step size if LU decomposition fails and retry
!       -exit after 5 consecutive fails
! --- --- --- --- --- --- --- --- --- --- --- --- ---
   IMPLICIT NONE

!~~~> Input arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) ::  Jac0(N,N)
#else
   REAL(kind=dp), INTENT(IN) ::  Jac0(LS_LU_NONZERO)
#endif   
   REAL(kind=dp), INTENT(IN) ::  gam
   INTEGER, INTENT(IN) ::  Direction
!~~~> Output arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(OUT) :: Ghimj(N,N)
#else
   REAL(kind=dp), INTENT(OUT) :: Ghimj(LS_LU_NONZERO)
#endif   
   LOGICAL, INTENT(OUT) ::  Singular
   INTEGER, INTENT(OUT) ::  Pivot(N)
!~~~> Inout arguments
   REAL(kind=dp), INTENT(INOUT) :: H   ! step size is decreased when LU fails
!~~~> Local variables
   INTEGER  :: i, ISING, Nconsecutive
   REAL(kind=dp) :: ghinv
   REAL(kind=dp), PARAMETER :: ONE  = 1.0_dp, HALF = 0.5_dp

   Nconsecutive = 0
   Singular = .TRUE.

   DO WHILE (Singular)

!~~~>    Construct Ghimj = 1/(H*gam) - Jac0
#ifdef FULL_ALGEBRA    
     !slim: CALL WCOPY(N*N,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(N*N,(-ONE),Ghimj,1)
     Ghimj = -Jac0
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(i,i) = Ghimj(i,i)+ghinv
     END DO
#else
     !slim: CALL WCOPY(LS_LU_NONZERO,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(LS_LU_NONZERO,(-ONE),Ghimj,1)
     Ghimj(1:LS_LU_NONZERO) = -Jac0(1:LS_LU_NONZERO)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(LS_LU_DIAG(i)) = Ghimj(LS_LU_DIAG(i))+ghinv
     END DO
#endif   
!~~~>    Compute LU decomposition
     CALL ros_Decomp( Ghimj, Pivot, ISING )
     IF (ISING == 0) THEN
!~~~>    If successful done
        Singular = .FALSE.
     ELSE ! ISING .ne. 0
!~~~>    If unsuccessful half the step size; if 5 consecutive fails then return
        ISTATUS(Nsng) = ISTATUS(Nsng) + 1
        Nconsecutive = Nconsecutive+1
        Singular = .TRUE.
        PRINT*,'Warning: LU Decomposition returned ISING = ',ISING
        IF (Nconsecutive <= 5) THEN ! Less than 5 consecutive failed decompositions
           H = H*HALF
        ELSE  ! More than 5 consecutive failed decompositions
           RETURN
        END IF  ! Nconsecutive
      END IF    ! ISING

   END DO ! WHILE Singular

  END SUBROUTINE ros_PrepareMatrix


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Decomp( A, Pivot, ISING )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the LU decomposition
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Inout variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(INOUT) :: A(N,N)
#else   
   REAL(kind=dp), INTENT(INOUT) :: A(LS_LU_NONZERO)
#endif
!~~~> Output variables
   INTEGER, INTENT(OUT) :: Pivot(N), ISING

#ifdef FULL_ALGEBRA    
   CALL  DGETRF( N, N, A, N, Pivot, ISING )
#else   
   CALL KppDecomp ( A, ISING,LS_NSEL,LS_LU_NONZERO,LS_LU_CROW,LS_LU_DIAG,LS_LU_ICOL)   
   Pivot(1) = 1
#endif
   ISTATUS(Ndec) = ISTATUS(Ndec) + 1

  END SUBROUTINE ros_Decomp


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Solve( A, Pivot, b)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the forward/backward substitution (using pre-computed LU decomposition)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Input variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) :: A(N,N)
   INTEGER :: ISING
#else   
   REAL(kind=dp), INTENT(IN) :: A(LS_LU_NONZERO)
#endif
   INTEGER, INTENT(IN) :: Pivot(N)
!~~~> InOut variables
   REAL(kind=dp), INTENT(INOUT) :: b(N)

#ifdef FULL_ALGEBRA    
   CALL  DGETRS( 'N', N , 1, A, N, Pivot, b, N, ISING )
   IF ( Info < 0 ) THEN
      PRINT*,"Error in DGETRS. ISING=",ISING
   END IF  
#else  
   	CALL KppSolve_7(LS_LU_NONZERO,LS_NSEL, A, b,LS_type )
#endif

   ISTATUS(Nsol) = ISTATUS(Nsol) + 1

  END SUBROUTINE ros_Solve



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   DOUBLE PRECISION g

    g = 1.0_dp + 1.0_dp/SQRT(2.0_dp)
    rosMethod = RS2
!~~~> Name of the method
    ros_Name = 'ROS-2'
!~~~> Number of stages
    ros_S = 2

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = (1.0_dp)/g
    ros_C(1) = (-2.0_dp)/g
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
!~~~> M_i = Coefficients for new step solution
    ros_M(1)= (3.0_dp)/(2.0_dp*g)
    ros_M(2)= (1.0_dp)/(2.0_dp*g)
! E_i = Coefficients for error estimator
    ros_E(1) = 1.0_dp/(2.0_dp*g)
    ros_E(2) = 1.0_dp/(2.0_dp*g)
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus one
    ros_ELO = 2.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.0_dp
    ros_Alpha(2) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = g
    ros_Gamma(2) =-g

 END SUBROUTINE Ros2


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   rosMethod = RS3
!~~~> Name of the method
   ros_Name = 'ROS-3'
!~~~> Number of stages
   ros_S = 3

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1)= 1.0_dp
   ros_A(2)= 1.0_dp
   ros_A(3)= 0.0_dp

   ros_C(1) = -0.10156171083877702091975600115545E+01_dp
   ros_C(2) =  0.40759956452537699824805835358067E+01_dp
   ros_C(3) =  0.92076794298330791242156818474003E+01_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) =  0.1E+01_dp
   ros_M(2) =  0.61697947043828245592553615689730E+01_dp
   ros_M(3) = -0.42772256543218573326238373806514_dp
! E_i = Coefficients for error estimator
   ros_E(1) =  0.5_dp
   ros_E(2) = -0.29079558716805469821718236208017E+01_dp
   ros_E(3) =  0.22354069897811569627360909276199_dp
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1)= 0.0_dp
   ros_Alpha(2)= 0.43586652150845899941601945119356_dp
   ros_Alpha(3)= 0.43586652150845899941601945119356_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1)= 0.43586652150845899941601945119356_dp
   ros_Gamma(2)= 0.24291996454816804366592249683314_dp
   ros_Gamma(3)= 0.21851380027664058511513169485832E+01_dp

  END SUBROUTINE Ros3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
!     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1990)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RS4
!~~~> Name of the method
   ros_Name = 'ROS-4'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.2000000000000000E+01_dp
   ros_A(2) = 0.1867943637803922E+01_dp
   ros_A(3) = 0.2344449711399156_dp
   ros_A(4) = ros_A(2)
   ros_A(5) = ros_A(3)
   ros_A(6) = 0.0_dp

   ros_C(1) =-0.7137615036412310E+01_dp
   ros_C(2) = 0.2580708087951457E+01_dp
   ros_C(3) = 0.6515950076447975_dp
   ros_C(4) =-0.2137148994382534E+01_dp
   ros_C(5) =-0.3214669691237626_dp
   ros_C(6) =-0.6949742501781779_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .TRUE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 0.2255570073418735E+01_dp
   ros_M(2) = 0.2870493262186792_dp
   ros_M(3) = 0.4353179431840180_dp
   ros_M(4) = 0.1093502252409163E+01_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) =-0.2815431932141155_dp
   ros_E(2) =-0.7276199124938920E-01_dp
   ros_E(3) =-0.1082196201495311_dp
   ros_E(4) =-0.1093502252409163E+01_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 4.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.1145640000000000E+01_dp
   ros_Alpha(3) = 0.6552168638155900_dp
   ros_Alpha(4) = ros_Alpha(3)
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5728200000000000_dp
   ros_Gamma(2) =-0.1769193891319233E+01_dp
   ros_Gamma(3) = 0.7592633437920482_dp
   ros_Gamma(4) =-0.1049021087100450_dp

  END SUBROUTINE Ros4

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- A STIFFLY-STABLE METHOD, 4 stages, order 3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RD3
!~~~> Name of the method
   ros_Name = 'RODAS-3'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.0_dp
   ros_A(2) = 2.0_dp
   ros_A(3) = 0.0_dp
   ros_A(4) = 2.0_dp
   ros_A(5) = 0.0_dp
   ros_A(6) = 1.0_dp

   ros_C(1) = 4.0_dp
   ros_C(2) = 1.0_dp
   ros_C(3) =-1.0_dp
   ros_C(4) = 1.0_dp
   ros_C(5) =-1.0_dp
   ros_C(6) =-(8.0_dp/3.0_dp)

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .FALSE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .TRUE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 2.0_dp
   ros_M(2) = 0.0_dp
   ros_M(3) = 1.0_dp
   ros_M(4) = 1.0_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) = 0.0_dp
   ros_E(2) = 0.0_dp
   ros_E(3) = 0.0_dp
   ros_E(4) = 1.0_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.0_dp
   ros_Alpha(3) = 1.0_dp
   ros_Alpha(4) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5_dp
   ros_Gamma(2) = 1.5_dp
   ros_Gamma(3) = 0.0_dp
   ros_Gamma(4) = 0.0_dp

  END SUBROUTINE Rodas3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1996)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RD4
!~~~> Name of the method
    ros_Name = 'RODAS-4'
!~~~> Number of stages
    ros_S = 6

!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.000_dp
    ros_Alpha(2) = 0.386_dp
    ros_Alpha(3) = 0.210_dp
    ros_Alpha(4) = 0.630_dp
    ros_Alpha(5) = 1.000_dp
    ros_Alpha(6) = 1.000_dp

!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = 0.2500000000000000_dp
    ros_Gamma(2) =-0.1043000000000000_dp
    ros_Gamma(3) = 0.1035000000000000_dp
    ros_Gamma(4) =-0.3620000000000023E-01_dp
    ros_Gamma(5) = 0.0_dp
    ros_Gamma(6) = 0.0_dp

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = 0.1544000000000000E+01_dp
    ros_A(2) = 0.9466785280815826_dp
    ros_A(3) = 0.2557011698983284_dp
    ros_A(4) = 0.3314825187068521E+01_dp
    ros_A(5) = 0.2896124015972201E+01_dp
    ros_A(6) = 0.9986419139977817_dp
    ros_A(7) = 0.1221224509226641E+01_dp
    ros_A(8) = 0.6019134481288629E+01_dp
    ros_A(9) = 0.1253708332932087E+02_dp
    ros_A(10) =-0.6878860361058950_dp
    ros_A(11) = ros_A(7)
    ros_A(12) = ros_A(8)
    ros_A(13) = ros_A(9)
    ros_A(14) = ros_A(10)
    ros_A(15) = 1.0_dp

    ros_C(1) =-0.5668800000000000E+01_dp
    ros_C(2) =-0.2430093356833875E+01_dp
    ros_C(3) =-0.2063599157091915_dp
    ros_C(4) =-0.1073529058151375_dp
    ros_C(5) =-0.9594562251023355E+01_dp
    ros_C(6) =-0.2047028614809616E+02_dp
    ros_C(7) = 0.7496443313967647E+01_dp
    ros_C(8) =-0.1024680431464352E+02_dp
    ros_C(9) =-0.3399990352819905E+02_dp
    ros_C(10) = 0.1170890893206160E+02_dp
    ros_C(11) = 0.8083246795921522E+01_dp
    ros_C(12) =-0.7981132988064893E+01_dp
    ros_C(13) =-0.3152159432874371E+02_dp
    ros_C(14) = 0.1631930543123136E+02_dp
    ros_C(15) =-0.6058818238834054E+01_dp

!~~~> M_i = Coefficients for new step solution
    ros_M(1) = ros_A(7)
    ros_M(2) = ros_A(8)
    ros_M(3) = ros_A(9)
    ros_M(4) = ros_A(10)
    ros_M(5) = 1.0_dp
    ros_M(6) = 1.0_dp

!~~~> E_i  = Coefficients for error estimator
    ros_E(1) = 0.0_dp
    ros_E(2) = 0.0_dp
    ros_E(3) = 0.0_dp
    ros_E(4) = 0.0_dp
    ros_E(5) = 0.0_dp
    ros_E(6) = 1.0_dp

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.
    ros_NewF(5) = .TRUE.
    ros_NewF(6) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 4.0_dp

  END SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! STIFFLY-STABLE W METHOD OF ORDER 3, WITH 4 STAGES
!
! J. RANG and L. ANGERMANN
! NEW ROSENBROCK W-METHODS OF ORDER 3
! FOR PARTIAL DIFFERENTIAL ALGEBRAIC
!        EQUATIONS OF INDEX 1                                             
! BIT Numerical Mathematics (2005) 45: 761-787
!  DOI: 10.1007/s10543-005-0035-y
! Table 4.1-4.2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RG3
!~~~> Name of the method
    ros_Name = 'RANG-3'
!~~~> Number of stages
    ros_S = 4

    ros_A(1) = 5.09052051067020d+00;
    ros_A(2) = 5.09052051067020d+00;
    ros_A(3) = 0.0d0;
    ros_A(4) = 4.97628111010787d+00;
    ros_A(5) = 2.77268164715849d-02;
    ros_A(6) = 2.29428036027904d-01;

    ros_C(1) = -1.16790812312283d+01;
    ros_C(2) = -1.64057326467367d+01;
    ros_C(3) = -2.77268164715850d-01;
    ros_C(4) = -8.38103960500476d+00;
    ros_C(5) = -8.48328409199343d-01;
    ros_C(6) =  2.87009860433106d-01;

    ros_M(1) =  5.22582761233094d+00;
    ros_M(2) = -5.56971148154165d-01;
    ros_M(3) =  3.57979469353645d-01;
    ros_M(4) =  1.72337398521064d+00;

    ros_E(1) = -5.16845212784040d+00;
    ros_E(2) = -1.26351942603842d+00;
    ros_E(3) = -1.11022302462516d-16;
    ros_E(4) =  2.22044604925031d-16;

    ros_Alpha(1) = 0.0d00;
    ros_Alpha(2) = 2.21878746765329d+00;
    ros_Alpha(3) = 2.21878746765329d+00;
    ros_Alpha(4) = 1.55392337535788d+00;

    ros_Gamma(1) =  4.35866521508459d-01;
    ros_Gamma(2) = -1.78292094614483d+00;
    ros_Gamma(3) = -2.46541900496934d+00;
    ros_Gamma(4) = -8.05529997906370d-01;

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 3.0_dp

  END SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END SUBROUTINE Rosenbrock_7


SUBROUTINE FunTemplate_7( T, Y, Ydot, LS_NSEL, LS_type )
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Function
 INTEGER,INTENT(IN)::LS_NSEL, LS_type
   REAL(kind=dp) :: T, Y(LS_NSEL)
   REAL(kind=dp) :: Ydot(LS_NSEL)
   REAL(kind=dp) :: Told
   Told = TIME
   TIME = T
   CALL Fun_7( Y, FIX, RCONST, Ydot, LS_NSEL)
   TIME = Told
END SUBROUTINE FunTemplate_7

SUBROUTINE JacTemplate_7( T, Y, Jcb, LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Jacobian
 USE gckpp_LinearAlgebra
 
 INTEGER,INTENT(IN)::LS_NSEL, LS_LU_NONZERO, LS_type
 INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO), LS_LU_ICOL(LS_LU_NONZERO)
    REAL(kind=dp) :: T, Y(LS_NSEL)
#ifdef FULL_ALGEBRA    
    REAL(kind=dp) :: JV(LS_LU_NONZERO), Jcb(LS_NSEL,LS_NSEL)
#else
    REAL(kind=dp) :: Jcb(LS_LU_NONZERO)
#endif   
    REAL(kind=dp) :: Told
#ifdef FULL_ALGEBRA    
    INTEGER :: i, j
#endif   

    Told = TIME
    TIME = T
#ifdef FULL_ALGEBRA        	
	CALL Jac_SP_7(Y, FIX, RCONST, JV,LS_NSEL,LS_LU_NONZERO)	
    DO j=1,LS_NSEL
      DO i=1,LS_NSEL
         Jcb(i,j) = 0.0_dp
      END DO
    END DO
    DO i=1,LS_LU_NONZERO
       Jcb(LS_LU_IROW(i),LS_LU_ICOL(i)) = JV(i)
    END DO
#else   
       CALL Jac_SP_7( Y, FIX, RCONST, Jcb,LS_NSEL,LS_LU_NONZERO)	 
#endif   
    TIME = Told
END SUBROUTINE JacTemplate_7



SUBROUTINE INTEGRATE_8( TIN, TOUT, LS_type,LS_NSEL, LS_NDEL,&
  Prate, Lrate, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )

   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: TIN  ! Start Time
   REAL(kind=dp), INTENT(IN) :: TOUT ! End Time
   
   ! Optional input parameters and statistics
   INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   REAL(kind=dp), INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   REAL(kind=dp), INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U
   INTEGER,		  INTENT(IN) :: LS_type,LS_NSEL,LS_NDEL
   REAL(kind=dp), INTENT(IN) :: Prate(NVAR),Lrate(NVAR)
   
   REAL(kind=dp) :: RCNTRL(20), RSTATUS(20),deltaT
   INTEGER       :: ICNTRL(20), ISTATUS(20), IERR
   REAL(kind=dp) :: VAR_selected(LS_NSEL),VAR_deleted(LS_NDEL),LS_P(LS_NDEL),LS_L(LS_NDEL)
   
   INTEGER, SAVE :: Ntotal = 0

   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.0_dp
   ISTATUS(:) = 0
   RSTATUS(:) = 0.0_dp
   deltaT =  TOUT-TIN
   !~~~> fine-tune the integrator:
   ICNTRL(1) = 0	! 0 - non-autonomous, 1 - autonomous
   ICNTRL(2) = 0	! 0 - vector tolerances, 1 - scalars

   ! If optional parameters are given, and if they are >0, 
   ! then they overwrite default settings. 
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF
   
        VAR_selected=VAR(select_ind_8)
	    VAR_deleted=VAR(delete_ind_8)
		LS_P=Prate(delete_ind_8)
		LS_L=Lrate(delete_ind_8)
        CALL Rosenbrock_8(LU_NSEL_8,VAR_selected,TIN,TOUT,ATOL,RTOL,&
          RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		  LU_NONZERO_8,LU_NSEL_8,LU_CROW_8,LU_DIAG_8,LU_IROW_8,LU_ICOL_8, LS_type)		
		WHERE(LS_L<=(0.01/deltaT))
			VAR_deleted=VAR_deleted+deltaT*(LS_P-LS_L*VAR_deleted)
		ELSEWHERE
			VAR_deleted=LS_P/LS_L+(VAR_deleted-LS_P/LS_L)*EXP(-LS_L*deltaT)
		END WHERE

	    VAR(select_ind_8)=VAR_selected
	    VAR(delete_ind_8)=VAR_deleted
   
   !~~~> Debug option: show no of steps
   ! Ntotal = Ntotal + ISTATUS(Nstp)
   ! PRINT*,'NSTEPS=',ISTATUS(Nstp),' (',Ntotal,')','  O3=', VAR(ind_O3)

   STEPMIN = RSTATUS(Nhexit)
   ! if optional parameters are given for output they 
   ! are updated with the return information
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE INTEGRATE_8

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE Rosenbrock_8(N,Y,Tstart,Tend, &
           AbsTol,RelTol,              &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		   LS_LU_NONZERO,LS_NSEL,LS_LU_CROW,LS_LU_DIAG,LS_LU_IROW,LS_LU_ICOL, LS_type)
  USE gckpp_Parameters
  USE gckpp_LinearAlgebra
  IMPLICIT NONE

   INTEGER,       INTENT(IN)    :: N,LS_LU_NONZERO,LS_NSEL,LS_type   
   INTEGER,INTENT(IN)::LS_LU_CROW(LS_NSEL+1),LS_LU_DIAG(LS_NSEL+1)
   INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO),LS_LU_ICOL(LS_LU_NONZERO)   
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
   REAL(kind=dp), INTENT(IN)    :: Tstart,Tend
   REAL(kind=dp), INTENT(IN)    :: AbsTol(N),RelTol(N)
   INTEGER,       INTENT(IN)    :: ICNTRL(20)
   REAL(kind=dp), INTENT(IN)    :: RCNTRL(20)
   INTEGER,       INTENT(INOUT) :: ISTATUS(20)
   REAL(kind=dp), INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT)   :: IERR
!~~~>  Parameters of the Rosenbrock method, up to 6 stages
   INTEGER ::  ros_S, rosMethod
   INTEGER, PARAMETER :: RS2=1, RS3=2, RS4=3, RD3=4, RD4=5, RG3=6
   REAL(kind=dp) :: ros_A(15), ros_C(15), ros_M(6), ros_E(6), &
                    ros_Alpha(6), ros_Gamma(6), ros_ELO
   LOGICAL :: ros_NewF(6)
   CHARACTER(LEN=12) :: ros_Name
!~~~>  Local variables
   REAL(kind=dp) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   REAL(kind=dp) :: Hmin, Hmax, Hstart
   REAL(kind=dp) :: Texit
   INTEGER       :: i, UplimTol, Max_no_steps
   LOGICAL       :: Autonomous, VectorTol
!~~~>   Parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp

!~~~>  Initialize statistics
   ISTATUS(1:8) = 0
   RSTATUS(1:3) = ZERO

!~~~>  Autonomous or time dependent ODE. Default is time dependent.
   Autonomous = .NOT.(ICNTRL(1) == 0)

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:N) and RelTol(1:N)
   IF (ICNTRL(2) == 0) THEN
      VectorTol = .TRUE.
      UplimTol  = N
   ELSE
      VectorTol = .FALSE.
      UplimTol  = 1
   END IF

!~~~>   Initialize the particular Rosenbrock method selected
   SELECT CASE (ICNTRL(3))
     CASE (1)
       CALL Ros2
     CASE (2)
       CALL Ros3
     CASE (3)
       CALL Ros4
     CASE (0,4)
       CALL Rodas3
     CASE (5)
       CALL Rodas4
     CASE (6)
       CALL Rang3
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(3)=',ICNTRL(3) 
       CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT

!~~~>   The maximum number of steps admitted
   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 200000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Unit roundoff (1+Roundoff>1)
   Roundoff = WLAMCH('E')

!~~~>  Lower bound on the step size: (positive value)
   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Upper bound on the step size: (positive value)
   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Starting step size: (positive value)
   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hold < FacMax
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Check if tolerances are reasonable
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO


!~~~>  CALL Rosenbrock method
   CALL ros_Integrator(Y, Tstart, Tend, Texit,   &
        AbsTol, RelTol,                          &
!  Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!  Error indicator
        IERR)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS !  SUBROUTINES internal to Rosenbrock
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 SUBROUTINE ros_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
   REAL(kind=dp), INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR
   
   IERR = Code
   PRINT * , &
     'Forced exit from Rosenbrock due to the following error:' 
     
   SELECT CASE (Code)
    CASE (-1)
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)    
      PRINT * , '--> Selected Rosenbrock method not implemented'
    CASE (-3)    
      PRINT * , '--> Hmin/Hmax/Hstart must be positive'
    CASE (-4)    
      PRINT * , '--> FacMin/FacMax/FacRej must be positive'
    CASE (-5) 
      PRINT * , '--> Improper tolerance values'
    CASE (-6) 
      PRINT * , '--> No of steps exceeds maximum bound'
    CASE (-7) 
      PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
    CASE (-8)    
      PRINT * , '--> Matrix is repeatedly singular'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT
   
   PRINT *, "T=", T, "and H=", H
     
 END SUBROUTINE ros_ErrorMsg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_Integrator (Y, Tstart, Tend, T,  &
        AbsTol, RelTol,                          &
!~~~> Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic Rosenbrock method
!      defined by ros_S (no of stages)
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

!~~~> Input: the initial condition at Tstart; Output: the solution at T
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
!~~~> Input: integration interval
   REAL(kind=dp), INTENT(IN) :: Tstart,Tend
!~~~> Output: time at which the solution is returned (T=Tend if success)
   REAL(kind=dp), INTENT(OUT) ::  T
!~~~> Input: tolerances
   REAL(kind=dp), INTENT(IN) ::  AbsTol(N), RelTol(N)
!~~~> Input: integration parameters
   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   REAL(kind=dp), INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   REAL(kind=dp), INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables
   REAL(kind=dp) :: Ynew(N), Fcn0(N), Fcn(N)
   REAL(kind=dp) :: K(N*ros_S), dFdT(N)
#ifdef FULL_ALGEBRA    
   REAL(kind=dp) :: Jac0(N,N), Ghimj(N,N)
#else
   REAL(kind=dp) :: Jac0(LS_LU_NONZERO), Ghimj(LS_LU_NONZERO)
#endif
   REAL(kind=dp) :: H, Hnew, HC, HG, Fac, Tau
   REAL(kind=dp) :: Err, Yerr(N)
   INTEGER :: Pivot(N), Direction, ioffset, j, istage
   LOGICAL :: RejectLastH, RejectMoreH, Singular
!~~~>  Local parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp
!~~~>  Locally called functions
!    REAL(kind=dp) WLAMCH
!    EXTERNAL WLAMCH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~>  Initial preparations
   T = Tstart
   RSTATUS(Nhexit) = ZERO
   H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   IF (ABS(H) <= 10.0_dp*Roundoff) H = DeltaMin

   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF
   H = Direction*H

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.

!~~~> Time loop begins below

TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) )

   IF ( ISTATUS(Nstp) > Max_no_steps ) THEN  ! Too many steps
      CALL ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
      CALL ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF

!~~~>  Limit H if necessary to avoid going beyond Tend
   H = MIN(H,ABS(Tend-T))

!~~~>   Compute the function at current time
   CALL FunTemplate_8(T,Y,Fcn0, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1

!~~~>  Compute the function derivative with respect to T
   IF (.NOT.Autonomous) THEN
      CALL ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
   END IF

!~~~>   Compute the Jacobian at current time
   CALL JacTemplate_8(T,Y,Jac0,LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)   
   ISTATUS(Njac) = ISTATUS(Njac) + 1

!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO

   CALL ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular)
   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       CALL ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S

      ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
       ioffset = N*(istage-1)

      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         !slim: CALL WCOPY(N,Fcn0,1,Fcn,1)
	 Fcn(1:N) = Fcn0(1:N)
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
         !slim: CALL WCOPY(N,Y,1,Ynew,1)
	 Ynew(1:N) = Y(1:N)
         DO j = 1, istage-1
           CALL WAXPY(N,ros_A((istage-1)*(istage-2)/2+j), &
            K(N*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL FunTemplate_8(Tau,Ynew,Fcn, LS_NSEL, LS_type)
         ISTATUS(Nfun) = ISTATUS(Nfun) + 1
       END IF ! if istage == 1 elseif ros_NewF(istage)
       !slim: CALL WCOPY(N,Fcn,1,K(ioffset+1),1)
       K(ioffset+1:ioffset+N) = Fcn(1:N)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL WAXPY(N,HC,K(N*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL WAXPY(N,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL ros_Solve(Ghimj, Pivot, K(ioffset+1))

   END DO Stage


!~~~>  Compute the new solution
   !slim: CALL WCOPY(N,Y,1,Ynew,1)
   Ynew(1:N) = Y(1:N)
   DO j=1,ros_S
         CALL WAXPY(N,ros_M(j),K(N*(j-1)+1),1,Ynew,1)
   END DO

!~~~>  Compute the error estimation
   !slim: CALL WSCAL(N,ZERO,Yerr,1)
   Yerr(1:N) = ZERO
   DO j=1,ros_S
        CALL WAXPY(N,ros_E(j),K(N*(j-1)+1),1,Yerr,1)
   END DO
   Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac

!~~~>  Check the error magnitude and adjust step size
   ISTATUS(Nstp) = ISTATUS(Nstp) + 1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
      ISTATUS(Nacc) = ISTATUS(Nacc) + 1
      !slim: CALL WCOPY(N,Ynew,1,Y,1)
      Y(1:N) = Ynew(1:N)
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H)
      END IF
      RSTATUS(Nhexit) = H
      RSTATUS(Nhnew)  = Hnew
      RSTATUS(Ntexit) = T
      RejectLastH = .FALSE.
      RejectMoreH = .FALSE.
      H = Hnew
      EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
   ELSE           !~~~> Reject step
      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (ISTATUS(Nacc) >= 1)  ISTATUS(Nrej) = ISTATUS(Nrej) + 1
   END IF ! Err <= 1

   END DO UntilAccepted

   
   END DO TimeLoop

!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful

  END SUBROUTINE ros_Integrator


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  REAL(kind=dp) FUNCTION ros_ErrorNorm ( Y, Ynew, Yerr, &
                               AbsTol, RelTol, VectorTol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

! Input arguments
   REAL(kind=dp), INTENT(IN) :: Y(N), Ynew(N), &
          Yerr(N), AbsTol(N), RelTol(N)
   LOGICAL, INTENT(IN) ::  VectorTol
! Local variables
   REAL(kind=dp) :: Err, Scale, Ymax
   INTEGER  :: i
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp

   Err = ZERO
   DO i=1,N
     Ymax = MAX(ABS(Y(i)),ABS(Ynew(i)))
     IF (VectorTol) THEN
       Scale = AbsTol(i)+RelTol(i)*Ymax
     ELSE
       Scale = AbsTol(1)+RelTol(1)*Ymax
     END IF
     Err = Err+(Yerr(i)/Scale)**2
   END DO
   Err  = SQRT(Err/N)

   ros_ErrorNorm = MAX(Err,1.0d-10)

  END FUNCTION ros_ErrorNorm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> The time partial derivative of the function by finite differences
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

!~~~> Input arguments
   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(N), Fcn0(N)
!~~~> Output arguments
   REAL(kind=dp), INTENT(OUT) :: dFdT(N)
!~~~> Local variables
   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL FunTemplate_8(T+Delta,Y,dFdT, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1
   CALL WAXPY(N,(-ONE),Fcn0,1,dFdT,1)
   CALL WSCAL(N,(ONE/Delta),dFdT,1)

  END SUBROUTINE ros_FunTimeDerivative


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_PrepareMatrix ( H, Direction, gam, &
             Jac0, Ghimj, Pivot, Singular )
! --- --- --- --- --- --- --- --- --- --- --- --- ---
!  Prepares the LHS matrix for stage calculations
!  1.  Construct Ghimj = 1/(H*ham) - Jac0
!      "(Gamma H) Inverse Minus Jacobian"
!  2.  Repeat LU decomposition of Ghimj until successful.
!       -half the step size if LU decomposition fails and retry
!       -exit after 5 consecutive fails
! --- --- --- --- --- --- --- --- --- --- --- --- ---
   IMPLICIT NONE

!~~~> Input arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) ::  Jac0(N,N)
#else
   REAL(kind=dp), INTENT(IN) ::  Jac0(LS_LU_NONZERO)
#endif   
   REAL(kind=dp), INTENT(IN) ::  gam
   INTEGER, INTENT(IN) ::  Direction
!~~~> Output arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(OUT) :: Ghimj(N,N)
#else
   REAL(kind=dp), INTENT(OUT) :: Ghimj(LS_LU_NONZERO)
#endif   
   LOGICAL, INTENT(OUT) ::  Singular
   INTEGER, INTENT(OUT) ::  Pivot(N)
!~~~> Inout arguments
   REAL(kind=dp), INTENT(INOUT) :: H   ! step size is decreased when LU fails
!~~~> Local variables
   INTEGER  :: i, ISING, Nconsecutive
   REAL(kind=dp) :: ghinv
   REAL(kind=dp), PARAMETER :: ONE  = 1.0_dp, HALF = 0.5_dp

   Nconsecutive = 0
   Singular = .TRUE.

   DO WHILE (Singular)

!~~~>    Construct Ghimj = 1/(H*gam) - Jac0
#ifdef FULL_ALGEBRA    
     !slim: CALL WCOPY(N*N,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(N*N,(-ONE),Ghimj,1)
     Ghimj = -Jac0
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(i,i) = Ghimj(i,i)+ghinv
     END DO
#else
     !slim: CALL WCOPY(LS_LU_NONZERO,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(LS_LU_NONZERO,(-ONE),Ghimj,1)
     Ghimj(1:LS_LU_NONZERO) = -Jac0(1:LS_LU_NONZERO)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(LS_LU_DIAG(i)) = Ghimj(LS_LU_DIAG(i))+ghinv
     END DO
#endif   
!~~~>    Compute LU decomposition
     CALL ros_Decomp( Ghimj, Pivot, ISING )
     IF (ISING == 0) THEN
!~~~>    If successful done
        Singular = .FALSE.
     ELSE ! ISING .ne. 0
!~~~>    If unsuccessful half the step size; if 5 consecutive fails then return
        ISTATUS(Nsng) = ISTATUS(Nsng) + 1
        Nconsecutive = Nconsecutive+1
        Singular = .TRUE.
        PRINT*,'Warning: LU Decomposition returned ISING = ',ISING
        IF (Nconsecutive <= 5) THEN ! Less than 5 consecutive failed decompositions
           H = H*HALF
        ELSE  ! More than 5 consecutive failed decompositions
           RETURN
        END IF  ! Nconsecutive
      END IF    ! ISING

   END DO ! WHILE Singular

  END SUBROUTINE ros_PrepareMatrix


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Decomp( A, Pivot, ISING )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the LU decomposition
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Inout variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(INOUT) :: A(N,N)
#else   
   REAL(kind=dp), INTENT(INOUT) :: A(LS_LU_NONZERO)
#endif
!~~~> Output variables
   INTEGER, INTENT(OUT) :: Pivot(N), ISING

#ifdef FULL_ALGEBRA    
   CALL  DGETRF( N, N, A, N, Pivot, ISING )
#else   
   CALL KppDecomp ( A, ISING,LS_NSEL,LS_LU_NONZERO,LS_LU_CROW,LS_LU_DIAG,LS_LU_ICOL)   
   Pivot(1) = 1
#endif
   ISTATUS(Ndec) = ISTATUS(Ndec) + 1

  END SUBROUTINE ros_Decomp


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Solve( A, Pivot, b)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the forward/backward substitution (using pre-computed LU decomposition)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Input variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) :: A(N,N)
   INTEGER :: ISING
#else   
   REAL(kind=dp), INTENT(IN) :: A(LS_LU_NONZERO)
#endif
   INTEGER, INTENT(IN) :: Pivot(N)
!~~~> InOut variables
   REAL(kind=dp), INTENT(INOUT) :: b(N)

#ifdef FULL_ALGEBRA    
   CALL  DGETRS( 'N', N , 1, A, N, Pivot, b, N, ISING )
   IF ( Info < 0 ) THEN
      PRINT*,"Error in DGETRS. ISING=",ISING
   END IF  
#else  
   	CALL KppSolve_8(LS_LU_NONZERO,LS_NSEL, A, b,LS_type )
#endif

   ISTATUS(Nsol) = ISTATUS(Nsol) + 1

  END SUBROUTINE ros_Solve



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   DOUBLE PRECISION g

    g = 1.0_dp + 1.0_dp/SQRT(2.0_dp)
    rosMethod = RS2
!~~~> Name of the method
    ros_Name = 'ROS-2'
!~~~> Number of stages
    ros_S = 2

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = (1.0_dp)/g
    ros_C(1) = (-2.0_dp)/g
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
!~~~> M_i = Coefficients for new step solution
    ros_M(1)= (3.0_dp)/(2.0_dp*g)
    ros_M(2)= (1.0_dp)/(2.0_dp*g)
! E_i = Coefficients for error estimator
    ros_E(1) = 1.0_dp/(2.0_dp*g)
    ros_E(2) = 1.0_dp/(2.0_dp*g)
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus one
    ros_ELO = 2.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.0_dp
    ros_Alpha(2) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = g
    ros_Gamma(2) =-g

 END SUBROUTINE Ros2


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   rosMethod = RS3
!~~~> Name of the method
   ros_Name = 'ROS-3'
!~~~> Number of stages
   ros_S = 3

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1)= 1.0_dp
   ros_A(2)= 1.0_dp
   ros_A(3)= 0.0_dp

   ros_C(1) = -0.10156171083877702091975600115545E+01_dp
   ros_C(2) =  0.40759956452537699824805835358067E+01_dp
   ros_C(3) =  0.92076794298330791242156818474003E+01_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) =  0.1E+01_dp
   ros_M(2) =  0.61697947043828245592553615689730E+01_dp
   ros_M(3) = -0.42772256543218573326238373806514_dp
! E_i = Coefficients for error estimator
   ros_E(1) =  0.5_dp
   ros_E(2) = -0.29079558716805469821718236208017E+01_dp
   ros_E(3) =  0.22354069897811569627360909276199_dp
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1)= 0.0_dp
   ros_Alpha(2)= 0.43586652150845899941601945119356_dp
   ros_Alpha(3)= 0.43586652150845899941601945119356_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1)= 0.43586652150845899941601945119356_dp
   ros_Gamma(2)= 0.24291996454816804366592249683314_dp
   ros_Gamma(3)= 0.21851380027664058511513169485832E+01_dp

  END SUBROUTINE Ros3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
!     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1990)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RS4
!~~~> Name of the method
   ros_Name = 'ROS-4'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.2000000000000000E+01_dp
   ros_A(2) = 0.1867943637803922E+01_dp
   ros_A(3) = 0.2344449711399156_dp
   ros_A(4) = ros_A(2)
   ros_A(5) = ros_A(3)
   ros_A(6) = 0.0_dp

   ros_C(1) =-0.7137615036412310E+01_dp
   ros_C(2) = 0.2580708087951457E+01_dp
   ros_C(3) = 0.6515950076447975_dp
   ros_C(4) =-0.2137148994382534E+01_dp
   ros_C(5) =-0.3214669691237626_dp
   ros_C(6) =-0.6949742501781779_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .TRUE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 0.2255570073418735E+01_dp
   ros_M(2) = 0.2870493262186792_dp
   ros_M(3) = 0.4353179431840180_dp
   ros_M(4) = 0.1093502252409163E+01_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) =-0.2815431932141155_dp
   ros_E(2) =-0.7276199124938920E-01_dp
   ros_E(3) =-0.1082196201495311_dp
   ros_E(4) =-0.1093502252409163E+01_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 4.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.1145640000000000E+01_dp
   ros_Alpha(3) = 0.6552168638155900_dp
   ros_Alpha(4) = ros_Alpha(3)
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5728200000000000_dp
   ros_Gamma(2) =-0.1769193891319233E+01_dp
   ros_Gamma(3) = 0.7592633437920482_dp
   ros_Gamma(4) =-0.1049021087100450_dp

  END SUBROUTINE Ros4

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- A STIFFLY-STABLE METHOD, 4 stages, order 3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RD3
!~~~> Name of the method
   ros_Name = 'RODAS-3'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.0_dp
   ros_A(2) = 2.0_dp
   ros_A(3) = 0.0_dp
   ros_A(4) = 2.0_dp
   ros_A(5) = 0.0_dp
   ros_A(6) = 1.0_dp

   ros_C(1) = 4.0_dp
   ros_C(2) = 1.0_dp
   ros_C(3) =-1.0_dp
   ros_C(4) = 1.0_dp
   ros_C(5) =-1.0_dp
   ros_C(6) =-(8.0_dp/3.0_dp)

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .FALSE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .TRUE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 2.0_dp
   ros_M(2) = 0.0_dp
   ros_M(3) = 1.0_dp
   ros_M(4) = 1.0_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) = 0.0_dp
   ros_E(2) = 0.0_dp
   ros_E(3) = 0.0_dp
   ros_E(4) = 1.0_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.0_dp
   ros_Alpha(3) = 1.0_dp
   ros_Alpha(4) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5_dp
   ros_Gamma(2) = 1.5_dp
   ros_Gamma(3) = 0.0_dp
   ros_Gamma(4) = 0.0_dp

  END SUBROUTINE Rodas3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1996)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RD4
!~~~> Name of the method
    ros_Name = 'RODAS-4'
!~~~> Number of stages
    ros_S = 6

!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.000_dp
    ros_Alpha(2) = 0.386_dp
    ros_Alpha(3) = 0.210_dp
    ros_Alpha(4) = 0.630_dp
    ros_Alpha(5) = 1.000_dp
    ros_Alpha(6) = 1.000_dp

!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = 0.2500000000000000_dp
    ros_Gamma(2) =-0.1043000000000000_dp
    ros_Gamma(3) = 0.1035000000000000_dp
    ros_Gamma(4) =-0.3620000000000023E-01_dp
    ros_Gamma(5) = 0.0_dp
    ros_Gamma(6) = 0.0_dp

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = 0.1544000000000000E+01_dp
    ros_A(2) = 0.9466785280815826_dp
    ros_A(3) = 0.2557011698983284_dp
    ros_A(4) = 0.3314825187068521E+01_dp
    ros_A(5) = 0.2896124015972201E+01_dp
    ros_A(6) = 0.9986419139977817_dp
    ros_A(7) = 0.1221224509226641E+01_dp
    ros_A(8) = 0.6019134481288629E+01_dp
    ros_A(9) = 0.1253708332932087E+02_dp
    ros_A(10) =-0.6878860361058950_dp
    ros_A(11) = ros_A(7)
    ros_A(12) = ros_A(8)
    ros_A(13) = ros_A(9)
    ros_A(14) = ros_A(10)
    ros_A(15) = 1.0_dp

    ros_C(1) =-0.5668800000000000E+01_dp
    ros_C(2) =-0.2430093356833875E+01_dp
    ros_C(3) =-0.2063599157091915_dp
    ros_C(4) =-0.1073529058151375_dp
    ros_C(5) =-0.9594562251023355E+01_dp
    ros_C(6) =-0.2047028614809616E+02_dp
    ros_C(7) = 0.7496443313967647E+01_dp
    ros_C(8) =-0.1024680431464352E+02_dp
    ros_C(9) =-0.3399990352819905E+02_dp
    ros_C(10) = 0.1170890893206160E+02_dp
    ros_C(11) = 0.8083246795921522E+01_dp
    ros_C(12) =-0.7981132988064893E+01_dp
    ros_C(13) =-0.3152159432874371E+02_dp
    ros_C(14) = 0.1631930543123136E+02_dp
    ros_C(15) =-0.6058818238834054E+01_dp

!~~~> M_i = Coefficients for new step solution
    ros_M(1) = ros_A(7)
    ros_M(2) = ros_A(8)
    ros_M(3) = ros_A(9)
    ros_M(4) = ros_A(10)
    ros_M(5) = 1.0_dp
    ros_M(6) = 1.0_dp

!~~~> E_i  = Coefficients for error estimator
    ros_E(1) = 0.0_dp
    ros_E(2) = 0.0_dp
    ros_E(3) = 0.0_dp
    ros_E(4) = 0.0_dp
    ros_E(5) = 0.0_dp
    ros_E(6) = 1.0_dp

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.
    ros_NewF(5) = .TRUE.
    ros_NewF(6) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 4.0_dp

  END SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! STIFFLY-STABLE W METHOD OF ORDER 3, WITH 4 STAGES
!
! J. RANG and L. ANGERMANN
! NEW ROSENBROCK W-METHODS OF ORDER 3
! FOR PARTIAL DIFFERENTIAL ALGEBRAIC
!        EQUATIONS OF INDEX 1                                             
! BIT Numerical Mathematics (2005) 45: 761-787
!  DOI: 10.1007/s10543-005-0035-y
! Table 4.1-4.2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RG3
!~~~> Name of the method
    ros_Name = 'RANG-3'
!~~~> Number of stages
    ros_S = 4

    ros_A(1) = 5.09052051067020d+00;
    ros_A(2) = 5.09052051067020d+00;
    ros_A(3) = 0.0d0;
    ros_A(4) = 4.97628111010787d+00;
    ros_A(5) = 2.77268164715849d-02;
    ros_A(6) = 2.29428036027904d-01;

    ros_C(1) = -1.16790812312283d+01;
    ros_C(2) = -1.64057326467367d+01;
    ros_C(3) = -2.77268164715850d-01;
    ros_C(4) = -8.38103960500476d+00;
    ros_C(5) = -8.48328409199343d-01;
    ros_C(6) =  2.87009860433106d-01;

    ros_M(1) =  5.22582761233094d+00;
    ros_M(2) = -5.56971148154165d-01;
    ros_M(3) =  3.57979469353645d-01;
    ros_M(4) =  1.72337398521064d+00;

    ros_E(1) = -5.16845212784040d+00;
    ros_E(2) = -1.26351942603842d+00;
    ros_E(3) = -1.11022302462516d-16;
    ros_E(4) =  2.22044604925031d-16;

    ros_Alpha(1) = 0.0d00;
    ros_Alpha(2) = 2.21878746765329d+00;
    ros_Alpha(3) = 2.21878746765329d+00;
    ros_Alpha(4) = 1.55392337535788d+00;

    ros_Gamma(1) =  4.35866521508459d-01;
    ros_Gamma(2) = -1.78292094614483d+00;
    ros_Gamma(3) = -2.46541900496934d+00;
    ros_Gamma(4) = -8.05529997906370d-01;

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 3.0_dp

  END SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END SUBROUTINE Rosenbrock_8


SUBROUTINE FunTemplate_8( T, Y, Ydot, LS_NSEL, LS_type )
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Function
 INTEGER,INTENT(IN)::LS_NSEL, LS_type
   REAL(kind=dp) :: T, Y(LS_NSEL)
   REAL(kind=dp) :: Ydot(LS_NSEL)
   REAL(kind=dp) :: Told
   Told = TIME
   TIME = T
   CALL Fun_8( Y, FIX, RCONST, Ydot, LS_NSEL)
   TIME = Told
END SUBROUTINE FunTemplate_8

SUBROUTINE JacTemplate_8( T, Y, Jcb, LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Jacobian
 USE gckpp_LinearAlgebra
 
 INTEGER,INTENT(IN)::LS_NSEL, LS_LU_NONZERO, LS_type
 INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO), LS_LU_ICOL(LS_LU_NONZERO)
    REAL(kind=dp) :: T, Y(LS_NSEL)
#ifdef FULL_ALGEBRA    
    REAL(kind=dp) :: JV(LS_LU_NONZERO), Jcb(LS_NSEL,LS_NSEL)
#else
    REAL(kind=dp) :: Jcb(LS_LU_NONZERO)
#endif   
    REAL(kind=dp) :: Told
#ifdef FULL_ALGEBRA    
    INTEGER :: i, j
#endif   

    Told = TIME
    TIME = T
#ifdef FULL_ALGEBRA        	
	CALL Jac_SP_8(Y, FIX, RCONST, JV,LS_NSEL,LS_LU_NONZERO)	
    DO j=1,LS_NSEL
      DO i=1,LS_NSEL
         Jcb(i,j) = 0.0_dp
      END DO
    END DO
    DO i=1,LS_LU_NONZERO
       Jcb(LS_LU_IROW(i),LS_LU_ICOL(i)) = JV(i)
    END DO
#else   
       CALL Jac_SP_8( Y, FIX, RCONST, Jcb,LS_NSEL,LS_LU_NONZERO)	 
#endif   
    TIME = Told
END SUBROUTINE JacTemplate_8



SUBROUTINE INTEGRATE_9( TIN, TOUT, LS_type,LS_NSEL, LS_NDEL,&
  Prate, Lrate, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )

   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: TIN  ! Start Time
   REAL(kind=dp), INTENT(IN) :: TOUT ! End Time
   
   ! Optional input parameters and statistics
   INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   REAL(kind=dp), INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   REAL(kind=dp), INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U
   INTEGER,		  INTENT(IN) :: LS_type,LS_NSEL,LS_NDEL
   REAL(kind=dp), INTENT(IN) :: Prate(NVAR),Lrate(NVAR)
   
   REAL(kind=dp) :: RCNTRL(20), RSTATUS(20),deltaT
   INTEGER       :: ICNTRL(20), ISTATUS(20), IERR
   REAL(kind=dp) :: VAR_selected(LS_NSEL),VAR_deleted(LS_NDEL),LS_P(LS_NDEL),LS_L(LS_NDEL)
   
   INTEGER, SAVE :: Ntotal = 0

   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.0_dp
   ISTATUS(:) = 0
   RSTATUS(:) = 0.0_dp
   deltaT =  TOUT-TIN
   !~~~> fine-tune the integrator:
   ICNTRL(1) = 0	! 0 - non-autonomous, 1 - autonomous
   ICNTRL(2) = 0	! 0 - vector tolerances, 1 - scalars

   ! If optional parameters are given, and if they are >0, 
   ! then they overwrite default settings. 
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF
   
        VAR_selected=VAR(select_ind_9)
	    VAR_deleted=VAR(delete_ind_9)
		LS_P=Prate(delete_ind_9)
		LS_L=Lrate(delete_ind_9)
        CALL Rosenbrock_9(LU_NSEL_9,VAR_selected,TIN,TOUT,ATOL,RTOL,&
          RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		  LU_NONZERO_9,LU_NSEL_9,LU_CROW_9,LU_DIAG_9,LU_IROW_9,LU_ICOL_9, LS_type)		
		WHERE(LS_L<=(0.01/deltaT))
			VAR_deleted=VAR_deleted+deltaT*(LS_P-LS_L*VAR_deleted)
		ELSEWHERE
			VAR_deleted=LS_P/LS_L+(VAR_deleted-LS_P/LS_L)*EXP(-LS_L*deltaT)
		END WHERE

	    VAR(select_ind_9)=VAR_selected
	    VAR(delete_ind_9)=VAR_deleted
   
   !~~~> Debug option: show no of steps
   ! Ntotal = Ntotal + ISTATUS(Nstp)
   ! PRINT*,'NSTEPS=',ISTATUS(Nstp),' (',Ntotal,')','  O3=', VAR(ind_O3)

   STEPMIN = RSTATUS(Nhexit)
   ! if optional parameters are given for output they 
   ! are updated with the return information
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE INTEGRATE_9

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE Rosenbrock_9(N,Y,Tstart,Tend, &
           AbsTol,RelTol,              &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		   LS_LU_NONZERO,LS_NSEL,LS_LU_CROW,LS_LU_DIAG,LS_LU_IROW,LS_LU_ICOL, LS_type)
  USE gckpp_Parameters
  USE gckpp_LinearAlgebra
  IMPLICIT NONE

   INTEGER,       INTENT(IN)    :: N,LS_LU_NONZERO,LS_NSEL,LS_type   
   INTEGER,INTENT(IN)::LS_LU_CROW(LS_NSEL+1),LS_LU_DIAG(LS_NSEL+1)
   INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO),LS_LU_ICOL(LS_LU_NONZERO)   
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
   REAL(kind=dp), INTENT(IN)    :: Tstart,Tend
   REAL(kind=dp), INTENT(IN)    :: AbsTol(N),RelTol(N)
   INTEGER,       INTENT(IN)    :: ICNTRL(20)
   REAL(kind=dp), INTENT(IN)    :: RCNTRL(20)
   INTEGER,       INTENT(INOUT) :: ISTATUS(20)
   REAL(kind=dp), INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT)   :: IERR
!~~~>  Parameters of the Rosenbrock method, up to 6 stages
   INTEGER ::  ros_S, rosMethod
   INTEGER, PARAMETER :: RS2=1, RS3=2, RS4=3, RD3=4, RD4=5, RG3=6
   REAL(kind=dp) :: ros_A(15), ros_C(15), ros_M(6), ros_E(6), &
                    ros_Alpha(6), ros_Gamma(6), ros_ELO
   LOGICAL :: ros_NewF(6)
   CHARACTER(LEN=12) :: ros_Name
!~~~>  Local variables
   REAL(kind=dp) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   REAL(kind=dp) :: Hmin, Hmax, Hstart
   REAL(kind=dp) :: Texit
   INTEGER       :: i, UplimTol, Max_no_steps
   LOGICAL       :: Autonomous, VectorTol
!~~~>   Parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp

!~~~>  Initialize statistics
   ISTATUS(1:8) = 0
   RSTATUS(1:3) = ZERO

!~~~>  Autonomous or time dependent ODE. Default is time dependent.
   Autonomous = .NOT.(ICNTRL(1) == 0)

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:N) and RelTol(1:N)
   IF (ICNTRL(2) == 0) THEN
      VectorTol = .TRUE.
      UplimTol  = N
   ELSE
      VectorTol = .FALSE.
      UplimTol  = 1
   END IF

!~~~>   Initialize the particular Rosenbrock method selected
   SELECT CASE (ICNTRL(3))
     CASE (1)
       CALL Ros2
     CASE (2)
       CALL Ros3
     CASE (3)
       CALL Ros4
     CASE (0,4)
       CALL Rodas3
     CASE (5)
       CALL Rodas4
     CASE (6)
       CALL Rang3
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(3)=',ICNTRL(3) 
       CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT

!~~~>   The maximum number of steps admitted
   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 200000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Unit roundoff (1+Roundoff>1)
   Roundoff = WLAMCH('E')

!~~~>  Lower bound on the step size: (positive value)
   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Upper bound on the step size: (positive value)
   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Starting step size: (positive value)
   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hold < FacMax
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Check if tolerances are reasonable
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO


!~~~>  CALL Rosenbrock method
   CALL ros_Integrator(Y, Tstart, Tend, Texit,   &
        AbsTol, RelTol,                          &
!  Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!  Error indicator
        IERR)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS !  SUBROUTINES internal to Rosenbrock
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 SUBROUTINE ros_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
   REAL(kind=dp), INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR
   
   IERR = Code
   PRINT * , &
     'Forced exit from Rosenbrock due to the following error:' 
     
   SELECT CASE (Code)
    CASE (-1)
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)    
      PRINT * , '--> Selected Rosenbrock method not implemented'
    CASE (-3)    
      PRINT * , '--> Hmin/Hmax/Hstart must be positive'
    CASE (-4)    
      PRINT * , '--> FacMin/FacMax/FacRej must be positive'
    CASE (-5) 
      PRINT * , '--> Improper tolerance values'
    CASE (-6) 
      PRINT * , '--> No of steps exceeds maximum bound'
    CASE (-7) 
      PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
    CASE (-8)    
      PRINT * , '--> Matrix is repeatedly singular'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT
   
   PRINT *, "T=", T, "and H=", H
     
 END SUBROUTINE ros_ErrorMsg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_Integrator (Y, Tstart, Tend, T,  &
        AbsTol, RelTol,                          &
!~~~> Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic Rosenbrock method
!      defined by ros_S (no of stages)
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

!~~~> Input: the initial condition at Tstart; Output: the solution at T
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
!~~~> Input: integration interval
   REAL(kind=dp), INTENT(IN) :: Tstart,Tend
!~~~> Output: time at which the solution is returned (T=Tend if success)
   REAL(kind=dp), INTENT(OUT) ::  T
!~~~> Input: tolerances
   REAL(kind=dp), INTENT(IN) ::  AbsTol(N), RelTol(N)
!~~~> Input: integration parameters
   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   REAL(kind=dp), INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   REAL(kind=dp), INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables
   REAL(kind=dp) :: Ynew(N), Fcn0(N), Fcn(N)
   REAL(kind=dp) :: K(N*ros_S), dFdT(N)
#ifdef FULL_ALGEBRA    
   REAL(kind=dp) :: Jac0(N,N), Ghimj(N,N)
#else
   REAL(kind=dp) :: Jac0(LS_LU_NONZERO), Ghimj(LS_LU_NONZERO)
#endif
   REAL(kind=dp) :: H, Hnew, HC, HG, Fac, Tau
   REAL(kind=dp) :: Err, Yerr(N)
   INTEGER :: Pivot(N), Direction, ioffset, j, istage
   LOGICAL :: RejectLastH, RejectMoreH, Singular
!~~~>  Local parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp
!~~~>  Locally called functions
!    REAL(kind=dp) WLAMCH
!    EXTERNAL WLAMCH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~>  Initial preparations
   T = Tstart
   RSTATUS(Nhexit) = ZERO
   H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   IF (ABS(H) <= 10.0_dp*Roundoff) H = DeltaMin

   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF
   H = Direction*H

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.

!~~~> Time loop begins below

TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) )

   IF ( ISTATUS(Nstp) > Max_no_steps ) THEN  ! Too many steps
      CALL ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
      CALL ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF

!~~~>  Limit H if necessary to avoid going beyond Tend
   H = MIN(H,ABS(Tend-T))

!~~~>   Compute the function at current time
   CALL FunTemplate_9(T,Y,Fcn0, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1

!~~~>  Compute the function derivative with respect to T
   IF (.NOT.Autonomous) THEN
      CALL ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
   END IF

!~~~>   Compute the Jacobian at current time
   CALL JacTemplate_9(T,Y,Jac0,LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)   
   ISTATUS(Njac) = ISTATUS(Njac) + 1

!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO

   CALL ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular)
   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       CALL ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S

      ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
       ioffset = N*(istage-1)

      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         !slim: CALL WCOPY(N,Fcn0,1,Fcn,1)
	 Fcn(1:N) = Fcn0(1:N)
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
         !slim: CALL WCOPY(N,Y,1,Ynew,1)
	 Ynew(1:N) = Y(1:N)
         DO j = 1, istage-1
           CALL WAXPY(N,ros_A((istage-1)*(istage-2)/2+j), &
            K(N*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL FunTemplate_9(Tau,Ynew,Fcn, LS_NSEL, LS_type)
         ISTATUS(Nfun) = ISTATUS(Nfun) + 1
       END IF ! if istage == 1 elseif ros_NewF(istage)
       !slim: CALL WCOPY(N,Fcn,1,K(ioffset+1),1)
       K(ioffset+1:ioffset+N) = Fcn(1:N)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL WAXPY(N,HC,K(N*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL WAXPY(N,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL ros_Solve(Ghimj, Pivot, K(ioffset+1))

   END DO Stage


!~~~>  Compute the new solution
   !slim: CALL WCOPY(N,Y,1,Ynew,1)
   Ynew(1:N) = Y(1:N)
   DO j=1,ros_S
         CALL WAXPY(N,ros_M(j),K(N*(j-1)+1),1,Ynew,1)
   END DO

!~~~>  Compute the error estimation
   !slim: CALL WSCAL(N,ZERO,Yerr,1)
   Yerr(1:N) = ZERO
   DO j=1,ros_S
        CALL WAXPY(N,ros_E(j),K(N*(j-1)+1),1,Yerr,1)
   END DO
   Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac

!~~~>  Check the error magnitude and adjust step size
   ISTATUS(Nstp) = ISTATUS(Nstp) + 1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
      ISTATUS(Nacc) = ISTATUS(Nacc) + 1
      !slim: CALL WCOPY(N,Ynew,1,Y,1)
      Y(1:N) = Ynew(1:N)
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H)
      END IF
      RSTATUS(Nhexit) = H
      RSTATUS(Nhnew)  = Hnew
      RSTATUS(Ntexit) = T
      RejectLastH = .FALSE.
      RejectMoreH = .FALSE.
      H = Hnew
      EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
   ELSE           !~~~> Reject step
      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (ISTATUS(Nacc) >= 1)  ISTATUS(Nrej) = ISTATUS(Nrej) + 1
   END IF ! Err <= 1

   END DO UntilAccepted

   
   END DO TimeLoop

!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful

  END SUBROUTINE ros_Integrator


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  REAL(kind=dp) FUNCTION ros_ErrorNorm ( Y, Ynew, Yerr, &
                               AbsTol, RelTol, VectorTol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

! Input arguments
   REAL(kind=dp), INTENT(IN) :: Y(N), Ynew(N), &
          Yerr(N), AbsTol(N), RelTol(N)
   LOGICAL, INTENT(IN) ::  VectorTol
! Local variables
   REAL(kind=dp) :: Err, Scale, Ymax
   INTEGER  :: i
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp

   Err = ZERO
   DO i=1,N
     Ymax = MAX(ABS(Y(i)),ABS(Ynew(i)))
     IF (VectorTol) THEN
       Scale = AbsTol(i)+RelTol(i)*Ymax
     ELSE
       Scale = AbsTol(1)+RelTol(1)*Ymax
     END IF
     Err = Err+(Yerr(i)/Scale)**2
   END DO
   Err  = SQRT(Err/N)

   ros_ErrorNorm = MAX(Err,1.0d-10)

  END FUNCTION ros_ErrorNorm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> The time partial derivative of the function by finite differences
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

!~~~> Input arguments
   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(N), Fcn0(N)
!~~~> Output arguments
   REAL(kind=dp), INTENT(OUT) :: dFdT(N)
!~~~> Local variables
   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL FunTemplate_9(T+Delta,Y,dFdT, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1
   CALL WAXPY(N,(-ONE),Fcn0,1,dFdT,1)
   CALL WSCAL(N,(ONE/Delta),dFdT,1)

  END SUBROUTINE ros_FunTimeDerivative


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_PrepareMatrix ( H, Direction, gam, &
             Jac0, Ghimj, Pivot, Singular )
! --- --- --- --- --- --- --- --- --- --- --- --- ---
!  Prepares the LHS matrix for stage calculations
!  1.  Construct Ghimj = 1/(H*ham) - Jac0
!      "(Gamma H) Inverse Minus Jacobian"
!  2.  Repeat LU decomposition of Ghimj until successful.
!       -half the step size if LU decomposition fails and retry
!       -exit after 5 consecutive fails
! --- --- --- --- --- --- --- --- --- --- --- --- ---
   IMPLICIT NONE

!~~~> Input arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) ::  Jac0(N,N)
#else
   REAL(kind=dp), INTENT(IN) ::  Jac0(LS_LU_NONZERO)
#endif   
   REAL(kind=dp), INTENT(IN) ::  gam
   INTEGER, INTENT(IN) ::  Direction
!~~~> Output arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(OUT) :: Ghimj(N,N)
#else
   REAL(kind=dp), INTENT(OUT) :: Ghimj(LS_LU_NONZERO)
#endif   
   LOGICAL, INTENT(OUT) ::  Singular
   INTEGER, INTENT(OUT) ::  Pivot(N)
!~~~> Inout arguments
   REAL(kind=dp), INTENT(INOUT) :: H   ! step size is decreased when LU fails
!~~~> Local variables
   INTEGER  :: i, ISING, Nconsecutive
   REAL(kind=dp) :: ghinv
   REAL(kind=dp), PARAMETER :: ONE  = 1.0_dp, HALF = 0.5_dp

   Nconsecutive = 0
   Singular = .TRUE.

   DO WHILE (Singular)

!~~~>    Construct Ghimj = 1/(H*gam) - Jac0
#ifdef FULL_ALGEBRA    
     !slim: CALL WCOPY(N*N,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(N*N,(-ONE),Ghimj,1)
     Ghimj = -Jac0
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(i,i) = Ghimj(i,i)+ghinv
     END DO
#else
     !slim: CALL WCOPY(LS_LU_NONZERO,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(LS_LU_NONZERO,(-ONE),Ghimj,1)
     Ghimj(1:LS_LU_NONZERO) = -Jac0(1:LS_LU_NONZERO)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(LS_LU_DIAG(i)) = Ghimj(LS_LU_DIAG(i))+ghinv
     END DO
#endif   
!~~~>    Compute LU decomposition
     CALL ros_Decomp( Ghimj, Pivot, ISING )
     IF (ISING == 0) THEN
!~~~>    If successful done
        Singular = .FALSE.
     ELSE ! ISING .ne. 0
!~~~>    If unsuccessful half the step size; if 5 consecutive fails then return
        ISTATUS(Nsng) = ISTATUS(Nsng) + 1
        Nconsecutive = Nconsecutive+1
        Singular = .TRUE.
        PRINT*,'Warning: LU Decomposition returned ISING = ',ISING
        IF (Nconsecutive <= 5) THEN ! Less than 5 consecutive failed decompositions
           H = H*HALF
        ELSE  ! More than 5 consecutive failed decompositions
           RETURN
        END IF  ! Nconsecutive
      END IF    ! ISING

   END DO ! WHILE Singular

  END SUBROUTINE ros_PrepareMatrix


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Decomp( A, Pivot, ISING )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the LU decomposition
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Inout variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(INOUT) :: A(N,N)
#else   
   REAL(kind=dp), INTENT(INOUT) :: A(LS_LU_NONZERO)
#endif
!~~~> Output variables
   INTEGER, INTENT(OUT) :: Pivot(N), ISING

#ifdef FULL_ALGEBRA    
   CALL  DGETRF( N, N, A, N, Pivot, ISING )
#else   
   CALL KppDecomp ( A, ISING,LS_NSEL,LS_LU_NONZERO,LS_LU_CROW,LS_LU_DIAG,LS_LU_ICOL)   
   Pivot(1) = 1
#endif
   ISTATUS(Ndec) = ISTATUS(Ndec) + 1

  END SUBROUTINE ros_Decomp


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Solve( A, Pivot, b)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the forward/backward substitution (using pre-computed LU decomposition)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Input variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) :: A(N,N)
   INTEGER :: ISING
#else   
   REAL(kind=dp), INTENT(IN) :: A(LS_LU_NONZERO)
#endif
   INTEGER, INTENT(IN) :: Pivot(N)
!~~~> InOut variables
   REAL(kind=dp), INTENT(INOUT) :: b(N)

#ifdef FULL_ALGEBRA    
   CALL  DGETRS( 'N', N , 1, A, N, Pivot, b, N, ISING )
   IF ( Info < 0 ) THEN
      PRINT*,"Error in DGETRS. ISING=",ISING
   END IF  
#else  
   	CALL KppSolve_9(LS_LU_NONZERO,LS_NSEL, A, b,LS_type )
#endif

   ISTATUS(Nsol) = ISTATUS(Nsol) + 1

  END SUBROUTINE ros_Solve



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   DOUBLE PRECISION g

    g = 1.0_dp + 1.0_dp/SQRT(2.0_dp)
    rosMethod = RS2
!~~~> Name of the method
    ros_Name = 'ROS-2'
!~~~> Number of stages
    ros_S = 2

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = (1.0_dp)/g
    ros_C(1) = (-2.0_dp)/g
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
!~~~> M_i = Coefficients for new step solution
    ros_M(1)= (3.0_dp)/(2.0_dp*g)
    ros_M(2)= (1.0_dp)/(2.0_dp*g)
! E_i = Coefficients for error estimator
    ros_E(1) = 1.0_dp/(2.0_dp*g)
    ros_E(2) = 1.0_dp/(2.0_dp*g)
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus one
    ros_ELO = 2.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.0_dp
    ros_Alpha(2) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = g
    ros_Gamma(2) =-g

 END SUBROUTINE Ros2


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   rosMethod = RS3
!~~~> Name of the method
   ros_Name = 'ROS-3'
!~~~> Number of stages
   ros_S = 3

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1)= 1.0_dp
   ros_A(2)= 1.0_dp
   ros_A(3)= 0.0_dp

   ros_C(1) = -0.10156171083877702091975600115545E+01_dp
   ros_C(2) =  0.40759956452537699824805835358067E+01_dp
   ros_C(3) =  0.92076794298330791242156818474003E+01_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) =  0.1E+01_dp
   ros_M(2) =  0.61697947043828245592553615689730E+01_dp
   ros_M(3) = -0.42772256543218573326238373806514_dp
! E_i = Coefficients for error estimator
   ros_E(1) =  0.5_dp
   ros_E(2) = -0.29079558716805469821718236208017E+01_dp
   ros_E(3) =  0.22354069897811569627360909276199_dp
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1)= 0.0_dp
   ros_Alpha(2)= 0.43586652150845899941601945119356_dp
   ros_Alpha(3)= 0.43586652150845899941601945119356_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1)= 0.43586652150845899941601945119356_dp
   ros_Gamma(2)= 0.24291996454816804366592249683314_dp
   ros_Gamma(3)= 0.21851380027664058511513169485832E+01_dp

  END SUBROUTINE Ros3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
!     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1990)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RS4
!~~~> Name of the method
   ros_Name = 'ROS-4'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.2000000000000000E+01_dp
   ros_A(2) = 0.1867943637803922E+01_dp
   ros_A(3) = 0.2344449711399156_dp
   ros_A(4) = ros_A(2)
   ros_A(5) = ros_A(3)
   ros_A(6) = 0.0_dp

   ros_C(1) =-0.7137615036412310E+01_dp
   ros_C(2) = 0.2580708087951457E+01_dp
   ros_C(3) = 0.6515950076447975_dp
   ros_C(4) =-0.2137148994382534E+01_dp
   ros_C(5) =-0.3214669691237626_dp
   ros_C(6) =-0.6949742501781779_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .TRUE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 0.2255570073418735E+01_dp
   ros_M(2) = 0.2870493262186792_dp
   ros_M(3) = 0.4353179431840180_dp
   ros_M(4) = 0.1093502252409163E+01_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) =-0.2815431932141155_dp
   ros_E(2) =-0.7276199124938920E-01_dp
   ros_E(3) =-0.1082196201495311_dp
   ros_E(4) =-0.1093502252409163E+01_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 4.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.1145640000000000E+01_dp
   ros_Alpha(3) = 0.6552168638155900_dp
   ros_Alpha(4) = ros_Alpha(3)
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5728200000000000_dp
   ros_Gamma(2) =-0.1769193891319233E+01_dp
   ros_Gamma(3) = 0.7592633437920482_dp
   ros_Gamma(4) =-0.1049021087100450_dp

  END SUBROUTINE Ros4

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- A STIFFLY-STABLE METHOD, 4 stages, order 3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RD3
!~~~> Name of the method
   ros_Name = 'RODAS-3'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.0_dp
   ros_A(2) = 2.0_dp
   ros_A(3) = 0.0_dp
   ros_A(4) = 2.0_dp
   ros_A(5) = 0.0_dp
   ros_A(6) = 1.0_dp

   ros_C(1) = 4.0_dp
   ros_C(2) = 1.0_dp
   ros_C(3) =-1.0_dp
   ros_C(4) = 1.0_dp
   ros_C(5) =-1.0_dp
   ros_C(6) =-(8.0_dp/3.0_dp)

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .FALSE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .TRUE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 2.0_dp
   ros_M(2) = 0.0_dp
   ros_M(3) = 1.0_dp
   ros_M(4) = 1.0_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) = 0.0_dp
   ros_E(2) = 0.0_dp
   ros_E(3) = 0.0_dp
   ros_E(4) = 1.0_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.0_dp
   ros_Alpha(3) = 1.0_dp
   ros_Alpha(4) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5_dp
   ros_Gamma(2) = 1.5_dp
   ros_Gamma(3) = 0.0_dp
   ros_Gamma(4) = 0.0_dp

  END SUBROUTINE Rodas3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1996)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RD4
!~~~> Name of the method
    ros_Name = 'RODAS-4'
!~~~> Number of stages
    ros_S = 6

!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.000_dp
    ros_Alpha(2) = 0.386_dp
    ros_Alpha(3) = 0.210_dp
    ros_Alpha(4) = 0.630_dp
    ros_Alpha(5) = 1.000_dp
    ros_Alpha(6) = 1.000_dp

!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = 0.2500000000000000_dp
    ros_Gamma(2) =-0.1043000000000000_dp
    ros_Gamma(3) = 0.1035000000000000_dp
    ros_Gamma(4) =-0.3620000000000023E-01_dp
    ros_Gamma(5) = 0.0_dp
    ros_Gamma(6) = 0.0_dp

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = 0.1544000000000000E+01_dp
    ros_A(2) = 0.9466785280815826_dp
    ros_A(3) = 0.2557011698983284_dp
    ros_A(4) = 0.3314825187068521E+01_dp
    ros_A(5) = 0.2896124015972201E+01_dp
    ros_A(6) = 0.9986419139977817_dp
    ros_A(7) = 0.1221224509226641E+01_dp
    ros_A(8) = 0.6019134481288629E+01_dp
    ros_A(9) = 0.1253708332932087E+02_dp
    ros_A(10) =-0.6878860361058950_dp
    ros_A(11) = ros_A(7)
    ros_A(12) = ros_A(8)
    ros_A(13) = ros_A(9)
    ros_A(14) = ros_A(10)
    ros_A(15) = 1.0_dp

    ros_C(1) =-0.5668800000000000E+01_dp
    ros_C(2) =-0.2430093356833875E+01_dp
    ros_C(3) =-0.2063599157091915_dp
    ros_C(4) =-0.1073529058151375_dp
    ros_C(5) =-0.9594562251023355E+01_dp
    ros_C(6) =-0.2047028614809616E+02_dp
    ros_C(7) = 0.7496443313967647E+01_dp
    ros_C(8) =-0.1024680431464352E+02_dp
    ros_C(9) =-0.3399990352819905E+02_dp
    ros_C(10) = 0.1170890893206160E+02_dp
    ros_C(11) = 0.8083246795921522E+01_dp
    ros_C(12) =-0.7981132988064893E+01_dp
    ros_C(13) =-0.3152159432874371E+02_dp
    ros_C(14) = 0.1631930543123136E+02_dp
    ros_C(15) =-0.6058818238834054E+01_dp

!~~~> M_i = Coefficients for new step solution
    ros_M(1) = ros_A(7)
    ros_M(2) = ros_A(8)
    ros_M(3) = ros_A(9)
    ros_M(4) = ros_A(10)
    ros_M(5) = 1.0_dp
    ros_M(6) = 1.0_dp

!~~~> E_i  = Coefficients for error estimator
    ros_E(1) = 0.0_dp
    ros_E(2) = 0.0_dp
    ros_E(3) = 0.0_dp
    ros_E(4) = 0.0_dp
    ros_E(5) = 0.0_dp
    ros_E(6) = 1.0_dp

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.
    ros_NewF(5) = .TRUE.
    ros_NewF(6) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 4.0_dp

  END SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! STIFFLY-STABLE W METHOD OF ORDER 3, WITH 4 STAGES
!
! J. RANG and L. ANGERMANN
! NEW ROSENBROCK W-METHODS OF ORDER 3
! FOR PARTIAL DIFFERENTIAL ALGEBRAIC
!        EQUATIONS OF INDEX 1                                             
! BIT Numerical Mathematics (2005) 45: 761-787
!  DOI: 10.1007/s10543-005-0035-y
! Table 4.1-4.2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RG3
!~~~> Name of the method
    ros_Name = 'RANG-3'
!~~~> Number of stages
    ros_S = 4

    ros_A(1) = 5.09052051067020d+00;
    ros_A(2) = 5.09052051067020d+00;
    ros_A(3) = 0.0d0;
    ros_A(4) = 4.97628111010787d+00;
    ros_A(5) = 2.77268164715849d-02;
    ros_A(6) = 2.29428036027904d-01;

    ros_C(1) = -1.16790812312283d+01;
    ros_C(2) = -1.64057326467367d+01;
    ros_C(3) = -2.77268164715850d-01;
    ros_C(4) = -8.38103960500476d+00;
    ros_C(5) = -8.48328409199343d-01;
    ros_C(6) =  2.87009860433106d-01;

    ros_M(1) =  5.22582761233094d+00;
    ros_M(2) = -5.56971148154165d-01;
    ros_M(3) =  3.57979469353645d-01;
    ros_M(4) =  1.72337398521064d+00;

    ros_E(1) = -5.16845212784040d+00;
    ros_E(2) = -1.26351942603842d+00;
    ros_E(3) = -1.11022302462516d-16;
    ros_E(4) =  2.22044604925031d-16;

    ros_Alpha(1) = 0.0d00;
    ros_Alpha(2) = 2.21878746765329d+00;
    ros_Alpha(3) = 2.21878746765329d+00;
    ros_Alpha(4) = 1.55392337535788d+00;

    ros_Gamma(1) =  4.35866521508459d-01;
    ros_Gamma(2) = -1.78292094614483d+00;
    ros_Gamma(3) = -2.46541900496934d+00;
    ros_Gamma(4) = -8.05529997906370d-01;

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 3.0_dp

  END SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END SUBROUTINE Rosenbrock_9


SUBROUTINE FunTemplate_9( T, Y, Ydot, LS_NSEL, LS_type )
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Function
 INTEGER,INTENT(IN)::LS_NSEL, LS_type
   REAL(kind=dp) :: T, Y(LS_NSEL)
   REAL(kind=dp) :: Ydot(LS_NSEL)
   REAL(kind=dp) :: Told
   Told = TIME
   TIME = T
   CALL Fun_9( Y, FIX, RCONST, Ydot, LS_NSEL)
   TIME = Told
END SUBROUTINE FunTemplate_9

SUBROUTINE JacTemplate_9( T, Y, Jcb, LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Jacobian
 USE gckpp_LinearAlgebra
 
 INTEGER,INTENT(IN)::LS_NSEL, LS_LU_NONZERO, LS_type
 INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO), LS_LU_ICOL(LS_LU_NONZERO)
    REAL(kind=dp) :: T, Y(LS_NSEL)
#ifdef FULL_ALGEBRA    
    REAL(kind=dp) :: JV(LS_LU_NONZERO), Jcb(LS_NSEL,LS_NSEL)
#else
    REAL(kind=dp) :: Jcb(LS_LU_NONZERO)
#endif   
    REAL(kind=dp) :: Told
#ifdef FULL_ALGEBRA    
    INTEGER :: i, j
#endif   

    Told = TIME
    TIME = T
#ifdef FULL_ALGEBRA        	
	CALL Jac_SP_9(Y, FIX, RCONST, JV,LS_NSEL,LS_LU_NONZERO)	
    DO j=1,LS_NSEL
      DO i=1,LS_NSEL
         Jcb(i,j) = 0.0_dp
      END DO
    END DO
    DO i=1,LS_LU_NONZERO
       Jcb(LS_LU_IROW(i),LS_LU_ICOL(i)) = JV(i)
    END DO
#else   
       CALL Jac_SP_9( Y, FIX, RCONST, Jcb,LS_NSEL,LS_LU_NONZERO)	 
#endif   
    TIME = Told
END SUBROUTINE JacTemplate_9



SUBROUTINE INTEGRATE_10( TIN, TOUT, LS_type,LS_NSEL, LS_NDEL,&
  Prate, Lrate, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )

   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: TIN  ! Start Time
   REAL(kind=dp), INTENT(IN) :: TOUT ! End Time
   
   ! Optional input parameters and statistics
   INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   REAL(kind=dp), INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   REAL(kind=dp), INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U
   INTEGER,		  INTENT(IN) :: LS_type,LS_NSEL,LS_NDEL
   REAL(kind=dp), INTENT(IN) :: Prate(NVAR),Lrate(NVAR)
   
   REAL(kind=dp) :: RCNTRL(20), RSTATUS(20),deltaT
   INTEGER       :: ICNTRL(20), ISTATUS(20), IERR
   REAL(kind=dp) :: VAR_selected(LS_NSEL),VAR_deleted(LS_NDEL),LS_P(LS_NDEL),LS_L(LS_NDEL)
   
   INTEGER, SAVE :: Ntotal = 0

   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.0_dp
   ISTATUS(:) = 0
   RSTATUS(:) = 0.0_dp
   deltaT =  TOUT-TIN
   !~~~> fine-tune the integrator:
   ICNTRL(1) = 0	! 0 - non-autonomous, 1 - autonomous
   ICNTRL(2) = 0	! 0 - vector tolerances, 1 - scalars

   ! If optional parameters are given, and if they are >0, 
   ! then they overwrite default settings. 
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF
   
        VAR_selected=VAR(select_ind_10)
	    VAR_deleted=VAR(delete_ind_10)
		LS_P=Prate(delete_ind_10)
		LS_L=Lrate(delete_ind_10)
        CALL Rosenbrock_10(LU_NSEL_10,VAR_selected,TIN,TOUT,ATOL,RTOL,&
          RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		  LU_NONZERO_10,LU_NSEL_10,LU_CROW_10,LU_DIAG_10,LU_IROW_10,LU_ICOL_10, LS_type)		
		WHERE(LS_L<=(0.01/deltaT))
			VAR_deleted=VAR_deleted+deltaT*(LS_P-LS_L*VAR_deleted)
		ELSEWHERE
			VAR_deleted=LS_P/LS_L+(VAR_deleted-LS_P/LS_L)*EXP(-LS_L*deltaT)
		END WHERE

	    VAR(select_ind_10)=VAR_selected
	    VAR(delete_ind_10)=VAR_deleted
   
   !~~~> Debug option: show no of steps
   ! Ntotal = Ntotal + ISTATUS(Nstp)
   ! PRINT*,'NSTEPS=',ISTATUS(Nstp),' (',Ntotal,')','  O3=', VAR(ind_O3)

   STEPMIN = RSTATUS(Nhexit)
   ! if optional parameters are given for output they 
   ! are updated with the return information
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE INTEGRATE_10

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE Rosenbrock_10(N,Y,Tstart,Tend, &
           AbsTol,RelTol,              &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		   LS_LU_NONZERO,LS_NSEL,LS_LU_CROW,LS_LU_DIAG,LS_LU_IROW,LS_LU_ICOL, LS_type)
  USE gckpp_Parameters
  USE gckpp_LinearAlgebra
  IMPLICIT NONE

   INTEGER,       INTENT(IN)    :: N,LS_LU_NONZERO,LS_NSEL,LS_type   
   INTEGER,INTENT(IN)::LS_LU_CROW(LS_NSEL+1),LS_LU_DIAG(LS_NSEL+1)
   INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO),LS_LU_ICOL(LS_LU_NONZERO)   
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
   REAL(kind=dp), INTENT(IN)    :: Tstart,Tend
   REAL(kind=dp), INTENT(IN)    :: AbsTol(N),RelTol(N)
   INTEGER,       INTENT(IN)    :: ICNTRL(20)
   REAL(kind=dp), INTENT(IN)    :: RCNTRL(20)
   INTEGER,       INTENT(INOUT) :: ISTATUS(20)
   REAL(kind=dp), INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT)   :: IERR
!~~~>  Parameters of the Rosenbrock method, up to 6 stages
   INTEGER ::  ros_S, rosMethod
   INTEGER, PARAMETER :: RS2=1, RS3=2, RS4=3, RD3=4, RD4=5, RG3=6
   REAL(kind=dp) :: ros_A(15), ros_C(15), ros_M(6), ros_E(6), &
                    ros_Alpha(6), ros_Gamma(6), ros_ELO
   LOGICAL :: ros_NewF(6)
   CHARACTER(LEN=12) :: ros_Name
!~~~>  Local variables
   REAL(kind=dp) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   REAL(kind=dp) :: Hmin, Hmax, Hstart
   REAL(kind=dp) :: Texit
   INTEGER       :: i, UplimTol, Max_no_steps
   LOGICAL       :: Autonomous, VectorTol
!~~~>   Parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp

!~~~>  Initialize statistics
   ISTATUS(1:8) = 0
   RSTATUS(1:3) = ZERO

!~~~>  Autonomous or time dependent ODE. Default is time dependent.
   Autonomous = .NOT.(ICNTRL(1) == 0)

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:N) and RelTol(1:N)
   IF (ICNTRL(2) == 0) THEN
      VectorTol = .TRUE.
      UplimTol  = N
   ELSE
      VectorTol = .FALSE.
      UplimTol  = 1
   END IF

!~~~>   Initialize the particular Rosenbrock method selected
   SELECT CASE (ICNTRL(3))
     CASE (1)
       CALL Ros2
     CASE (2)
       CALL Ros3
     CASE (3)
       CALL Ros4
     CASE (0,4)
       CALL Rodas3
     CASE (5)
       CALL Rodas4
     CASE (6)
       CALL Rang3
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(3)=',ICNTRL(3) 
       CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT

!~~~>   The maximum number of steps admitted
   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 200000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Unit roundoff (1+Roundoff>1)
   Roundoff = WLAMCH('E')

!~~~>  Lower bound on the step size: (positive value)
   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Upper bound on the step size: (positive value)
   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Starting step size: (positive value)
   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hold < FacMax
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Check if tolerances are reasonable
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO


!~~~>  CALL Rosenbrock method
   CALL ros_Integrator(Y, Tstart, Tend, Texit,   &
        AbsTol, RelTol,                          &
!  Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!  Error indicator
        IERR)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS !  SUBROUTINES internal to Rosenbrock
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 SUBROUTINE ros_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
   REAL(kind=dp), INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR
   
   IERR = Code
   PRINT * , &
     'Forced exit from Rosenbrock due to the following error:' 
     
   SELECT CASE (Code)
    CASE (-1)
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)    
      PRINT * , '--> Selected Rosenbrock method not implemented'
    CASE (-3)    
      PRINT * , '--> Hmin/Hmax/Hstart must be positive'
    CASE (-4)    
      PRINT * , '--> FacMin/FacMax/FacRej must be positive'
    CASE (-5) 
      PRINT * , '--> Improper tolerance values'
    CASE (-6) 
      PRINT * , '--> No of steps exceeds maximum bound'
    CASE (-7) 
      PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
    CASE (-8)    
      PRINT * , '--> Matrix is repeatedly singular'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT
   
   PRINT *, "T=", T, "and H=", H
     
 END SUBROUTINE ros_ErrorMsg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_Integrator (Y, Tstart, Tend, T,  &
        AbsTol, RelTol,                          &
!~~~> Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic Rosenbrock method
!      defined by ros_S (no of stages)
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

!~~~> Input: the initial condition at Tstart; Output: the solution at T
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
!~~~> Input: integration interval
   REAL(kind=dp), INTENT(IN) :: Tstart,Tend
!~~~> Output: time at which the solution is returned (T=Tend if success)
   REAL(kind=dp), INTENT(OUT) ::  T
!~~~> Input: tolerances
   REAL(kind=dp), INTENT(IN) ::  AbsTol(N), RelTol(N)
!~~~> Input: integration parameters
   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   REAL(kind=dp), INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   REAL(kind=dp), INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables
   REAL(kind=dp) :: Ynew(N), Fcn0(N), Fcn(N)
   REAL(kind=dp) :: K(N*ros_S), dFdT(N)
#ifdef FULL_ALGEBRA    
   REAL(kind=dp) :: Jac0(N,N), Ghimj(N,N)
#else
   REAL(kind=dp) :: Jac0(LS_LU_NONZERO), Ghimj(LS_LU_NONZERO)
#endif
   REAL(kind=dp) :: H, Hnew, HC, HG, Fac, Tau
   REAL(kind=dp) :: Err, Yerr(N)
   INTEGER :: Pivot(N), Direction, ioffset, j, istage
   LOGICAL :: RejectLastH, RejectMoreH, Singular
!~~~>  Local parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp
!~~~>  Locally called functions
!    REAL(kind=dp) WLAMCH
!    EXTERNAL WLAMCH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~>  Initial preparations
   T = Tstart
   RSTATUS(Nhexit) = ZERO
   H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   IF (ABS(H) <= 10.0_dp*Roundoff) H = DeltaMin

   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF
   H = Direction*H

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.

!~~~> Time loop begins below

TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) )

   IF ( ISTATUS(Nstp) > Max_no_steps ) THEN  ! Too many steps
      CALL ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
      CALL ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF

!~~~>  Limit H if necessary to avoid going beyond Tend
   H = MIN(H,ABS(Tend-T))

!~~~>   Compute the function at current time
   CALL FunTemplate_10(T,Y,Fcn0, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1

!~~~>  Compute the function derivative with respect to T
   IF (.NOT.Autonomous) THEN
      CALL ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
   END IF

!~~~>   Compute the Jacobian at current time
   CALL JacTemplate_10(T,Y,Jac0,LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)   
   ISTATUS(Njac) = ISTATUS(Njac) + 1

!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO

   CALL ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular)
   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       CALL ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S

      ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
       ioffset = N*(istage-1)

      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         !slim: CALL WCOPY(N,Fcn0,1,Fcn,1)
	 Fcn(1:N) = Fcn0(1:N)
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
         !slim: CALL WCOPY(N,Y,1,Ynew,1)
	 Ynew(1:N) = Y(1:N)
         DO j = 1, istage-1
           CALL WAXPY(N,ros_A((istage-1)*(istage-2)/2+j), &
            K(N*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL FunTemplate_10(Tau,Ynew,Fcn, LS_NSEL, LS_type)
         ISTATUS(Nfun) = ISTATUS(Nfun) + 1
       END IF ! if istage == 1 elseif ros_NewF(istage)
       !slim: CALL WCOPY(N,Fcn,1,K(ioffset+1),1)
       K(ioffset+1:ioffset+N) = Fcn(1:N)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL WAXPY(N,HC,K(N*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL WAXPY(N,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL ros_Solve(Ghimj, Pivot, K(ioffset+1))

   END DO Stage


!~~~>  Compute the new solution
   !slim: CALL WCOPY(N,Y,1,Ynew,1)
   Ynew(1:N) = Y(1:N)
   DO j=1,ros_S
         CALL WAXPY(N,ros_M(j),K(N*(j-1)+1),1,Ynew,1)
   END DO

!~~~>  Compute the error estimation
   !slim: CALL WSCAL(N,ZERO,Yerr,1)
   Yerr(1:N) = ZERO
   DO j=1,ros_S
        CALL WAXPY(N,ros_E(j),K(N*(j-1)+1),1,Yerr,1)
   END DO
   Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac

!~~~>  Check the error magnitude and adjust step size
   ISTATUS(Nstp) = ISTATUS(Nstp) + 1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
      ISTATUS(Nacc) = ISTATUS(Nacc) + 1
      !slim: CALL WCOPY(N,Ynew,1,Y,1)
      Y(1:N) = Ynew(1:N)
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H)
      END IF
      RSTATUS(Nhexit) = H
      RSTATUS(Nhnew)  = Hnew
      RSTATUS(Ntexit) = T
      RejectLastH = .FALSE.
      RejectMoreH = .FALSE.
      H = Hnew
      EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
   ELSE           !~~~> Reject step
      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (ISTATUS(Nacc) >= 1)  ISTATUS(Nrej) = ISTATUS(Nrej) + 1
   END IF ! Err <= 1

   END DO UntilAccepted

   
   END DO TimeLoop

!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful

  END SUBROUTINE ros_Integrator


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  REAL(kind=dp) FUNCTION ros_ErrorNorm ( Y, Ynew, Yerr, &
                               AbsTol, RelTol, VectorTol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

! Input arguments
   REAL(kind=dp), INTENT(IN) :: Y(N), Ynew(N), &
          Yerr(N), AbsTol(N), RelTol(N)
   LOGICAL, INTENT(IN) ::  VectorTol
! Local variables
   REAL(kind=dp) :: Err, Scale, Ymax
   INTEGER  :: i
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp

   Err = ZERO
   DO i=1,N
     Ymax = MAX(ABS(Y(i)),ABS(Ynew(i)))
     IF (VectorTol) THEN
       Scale = AbsTol(i)+RelTol(i)*Ymax
     ELSE
       Scale = AbsTol(1)+RelTol(1)*Ymax
     END IF
     Err = Err+(Yerr(i)/Scale)**2
   END DO
   Err  = SQRT(Err/N)

   ros_ErrorNorm = MAX(Err,1.0d-10)

  END FUNCTION ros_ErrorNorm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> The time partial derivative of the function by finite differences
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

!~~~> Input arguments
   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(N), Fcn0(N)
!~~~> Output arguments
   REAL(kind=dp), INTENT(OUT) :: dFdT(N)
!~~~> Local variables
   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL FunTemplate_10(T+Delta,Y,dFdT, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1
   CALL WAXPY(N,(-ONE),Fcn0,1,dFdT,1)
   CALL WSCAL(N,(ONE/Delta),dFdT,1)

  END SUBROUTINE ros_FunTimeDerivative


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_PrepareMatrix ( H, Direction, gam, &
             Jac0, Ghimj, Pivot, Singular )
! --- --- --- --- --- --- --- --- --- --- --- --- ---
!  Prepares the LHS matrix for stage calculations
!  1.  Construct Ghimj = 1/(H*ham) - Jac0
!      "(Gamma H) Inverse Minus Jacobian"
!  2.  Repeat LU decomposition of Ghimj until successful.
!       -half the step size if LU decomposition fails and retry
!       -exit after 5 consecutive fails
! --- --- --- --- --- --- --- --- --- --- --- --- ---
   IMPLICIT NONE

!~~~> Input arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) ::  Jac0(N,N)
#else
   REAL(kind=dp), INTENT(IN) ::  Jac0(LS_LU_NONZERO)
#endif   
   REAL(kind=dp), INTENT(IN) ::  gam
   INTEGER, INTENT(IN) ::  Direction
!~~~> Output arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(OUT) :: Ghimj(N,N)
#else
   REAL(kind=dp), INTENT(OUT) :: Ghimj(LS_LU_NONZERO)
#endif   
   LOGICAL, INTENT(OUT) ::  Singular
   INTEGER, INTENT(OUT) ::  Pivot(N)
!~~~> Inout arguments
   REAL(kind=dp), INTENT(INOUT) :: H   ! step size is decreased when LU fails
!~~~> Local variables
   INTEGER  :: i, ISING, Nconsecutive
   REAL(kind=dp) :: ghinv
   REAL(kind=dp), PARAMETER :: ONE  = 1.0_dp, HALF = 0.5_dp

   Nconsecutive = 0
   Singular = .TRUE.

   DO WHILE (Singular)

!~~~>    Construct Ghimj = 1/(H*gam) - Jac0
#ifdef FULL_ALGEBRA    
     !slim: CALL WCOPY(N*N,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(N*N,(-ONE),Ghimj,1)
     Ghimj = -Jac0
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(i,i) = Ghimj(i,i)+ghinv
     END DO
#else
     !slim: CALL WCOPY(LS_LU_NONZERO,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(LS_LU_NONZERO,(-ONE),Ghimj,1)
     Ghimj(1:LS_LU_NONZERO) = -Jac0(1:LS_LU_NONZERO)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(LS_LU_DIAG(i)) = Ghimj(LS_LU_DIAG(i))+ghinv
     END DO
#endif   
!~~~>    Compute LU decomposition
     CALL ros_Decomp( Ghimj, Pivot, ISING )
     IF (ISING == 0) THEN
!~~~>    If successful done
        Singular = .FALSE.
     ELSE ! ISING .ne. 0
!~~~>    If unsuccessful half the step size; if 5 consecutive fails then return
        ISTATUS(Nsng) = ISTATUS(Nsng) + 1
        Nconsecutive = Nconsecutive+1
        Singular = .TRUE.
        PRINT*,'Warning: LU Decomposition returned ISING = ',ISING
        IF (Nconsecutive <= 5) THEN ! Less than 5 consecutive failed decompositions
           H = H*HALF
        ELSE  ! More than 5 consecutive failed decompositions
           RETURN
        END IF  ! Nconsecutive
      END IF    ! ISING

   END DO ! WHILE Singular

  END SUBROUTINE ros_PrepareMatrix


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Decomp( A, Pivot, ISING )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the LU decomposition
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Inout variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(INOUT) :: A(N,N)
#else   
   REAL(kind=dp), INTENT(INOUT) :: A(LS_LU_NONZERO)
#endif
!~~~> Output variables
   INTEGER, INTENT(OUT) :: Pivot(N), ISING

#ifdef FULL_ALGEBRA    
   CALL  DGETRF( N, N, A, N, Pivot, ISING )
#else   
   CALL KppDecomp ( A, ISING,LS_NSEL,LS_LU_NONZERO,LS_LU_CROW,LS_LU_DIAG,LS_LU_ICOL)   
   Pivot(1) = 1
#endif
   ISTATUS(Ndec) = ISTATUS(Ndec) + 1

  END SUBROUTINE ros_Decomp


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Solve( A, Pivot, b)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the forward/backward substitution (using pre-computed LU decomposition)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Input variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) :: A(N,N)
   INTEGER :: ISING
#else   
   REAL(kind=dp), INTENT(IN) :: A(LS_LU_NONZERO)
#endif
   INTEGER, INTENT(IN) :: Pivot(N)
!~~~> InOut variables
   REAL(kind=dp), INTENT(INOUT) :: b(N)

#ifdef FULL_ALGEBRA    
   CALL  DGETRS( 'N', N , 1, A, N, Pivot, b, N, ISING )
   IF ( Info < 0 ) THEN
      PRINT*,"Error in DGETRS. ISING=",ISING
   END IF  
#else  
   	CALL KppSolve_10(LS_LU_NONZERO,LS_NSEL, A, b,LS_type )
#endif

   ISTATUS(Nsol) = ISTATUS(Nsol) + 1

  END SUBROUTINE ros_Solve



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   DOUBLE PRECISION g

    g = 1.0_dp + 1.0_dp/SQRT(2.0_dp)
    rosMethod = RS2
!~~~> Name of the method
    ros_Name = 'ROS-2'
!~~~> Number of stages
    ros_S = 2

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = (1.0_dp)/g
    ros_C(1) = (-2.0_dp)/g
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
!~~~> M_i = Coefficients for new step solution
    ros_M(1)= (3.0_dp)/(2.0_dp*g)
    ros_M(2)= (1.0_dp)/(2.0_dp*g)
! E_i = Coefficients for error estimator
    ros_E(1) = 1.0_dp/(2.0_dp*g)
    ros_E(2) = 1.0_dp/(2.0_dp*g)
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus one
    ros_ELO = 2.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.0_dp
    ros_Alpha(2) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = g
    ros_Gamma(2) =-g

 END SUBROUTINE Ros2


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   rosMethod = RS3
!~~~> Name of the method
   ros_Name = 'ROS-3'
!~~~> Number of stages
   ros_S = 3

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1)= 1.0_dp
   ros_A(2)= 1.0_dp
   ros_A(3)= 0.0_dp

   ros_C(1) = -0.10156171083877702091975600115545E+01_dp
   ros_C(2) =  0.40759956452537699824805835358067E+01_dp
   ros_C(3) =  0.92076794298330791242156818474003E+01_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) =  0.1E+01_dp
   ros_M(2) =  0.61697947043828245592553615689730E+01_dp
   ros_M(3) = -0.42772256543218573326238373806514_dp
! E_i = Coefficients for error estimator
   ros_E(1) =  0.5_dp
   ros_E(2) = -0.29079558716805469821718236208017E+01_dp
   ros_E(3) =  0.22354069897811569627360909276199_dp
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1)= 0.0_dp
   ros_Alpha(2)= 0.43586652150845899941601945119356_dp
   ros_Alpha(3)= 0.43586652150845899941601945119356_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1)= 0.43586652150845899941601945119356_dp
   ros_Gamma(2)= 0.24291996454816804366592249683314_dp
   ros_Gamma(3)= 0.21851380027664058511513169485832E+01_dp

  END SUBROUTINE Ros3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
!     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1990)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RS4
!~~~> Name of the method
   ros_Name = 'ROS-4'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.2000000000000000E+01_dp
   ros_A(2) = 0.1867943637803922E+01_dp
   ros_A(3) = 0.2344449711399156_dp
   ros_A(4) = ros_A(2)
   ros_A(5) = ros_A(3)
   ros_A(6) = 0.0_dp

   ros_C(1) =-0.7137615036412310E+01_dp
   ros_C(2) = 0.2580708087951457E+01_dp
   ros_C(3) = 0.6515950076447975_dp
   ros_C(4) =-0.2137148994382534E+01_dp
   ros_C(5) =-0.3214669691237626_dp
   ros_C(6) =-0.6949742501781779_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .TRUE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 0.2255570073418735E+01_dp
   ros_M(2) = 0.2870493262186792_dp
   ros_M(3) = 0.4353179431840180_dp
   ros_M(4) = 0.1093502252409163E+01_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) =-0.2815431932141155_dp
   ros_E(2) =-0.7276199124938920E-01_dp
   ros_E(3) =-0.1082196201495311_dp
   ros_E(4) =-0.1093502252409163E+01_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 4.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.1145640000000000E+01_dp
   ros_Alpha(3) = 0.6552168638155900_dp
   ros_Alpha(4) = ros_Alpha(3)
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5728200000000000_dp
   ros_Gamma(2) =-0.1769193891319233E+01_dp
   ros_Gamma(3) = 0.7592633437920482_dp
   ros_Gamma(4) =-0.1049021087100450_dp

  END SUBROUTINE Ros4

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- A STIFFLY-STABLE METHOD, 4 stages, order 3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RD3
!~~~> Name of the method
   ros_Name = 'RODAS-3'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.0_dp
   ros_A(2) = 2.0_dp
   ros_A(3) = 0.0_dp
   ros_A(4) = 2.0_dp
   ros_A(5) = 0.0_dp
   ros_A(6) = 1.0_dp

   ros_C(1) = 4.0_dp
   ros_C(2) = 1.0_dp
   ros_C(3) =-1.0_dp
   ros_C(4) = 1.0_dp
   ros_C(5) =-1.0_dp
   ros_C(6) =-(8.0_dp/3.0_dp)

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .FALSE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .TRUE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 2.0_dp
   ros_M(2) = 0.0_dp
   ros_M(3) = 1.0_dp
   ros_M(4) = 1.0_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) = 0.0_dp
   ros_E(2) = 0.0_dp
   ros_E(3) = 0.0_dp
   ros_E(4) = 1.0_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.0_dp
   ros_Alpha(3) = 1.0_dp
   ros_Alpha(4) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5_dp
   ros_Gamma(2) = 1.5_dp
   ros_Gamma(3) = 0.0_dp
   ros_Gamma(4) = 0.0_dp

  END SUBROUTINE Rodas3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1996)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RD4
!~~~> Name of the method
    ros_Name = 'RODAS-4'
!~~~> Number of stages
    ros_S = 6

!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.000_dp
    ros_Alpha(2) = 0.386_dp
    ros_Alpha(3) = 0.210_dp
    ros_Alpha(4) = 0.630_dp
    ros_Alpha(5) = 1.000_dp
    ros_Alpha(6) = 1.000_dp

!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = 0.2500000000000000_dp
    ros_Gamma(2) =-0.1043000000000000_dp
    ros_Gamma(3) = 0.1035000000000000_dp
    ros_Gamma(4) =-0.3620000000000023E-01_dp
    ros_Gamma(5) = 0.0_dp
    ros_Gamma(6) = 0.0_dp

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = 0.1544000000000000E+01_dp
    ros_A(2) = 0.9466785280815826_dp
    ros_A(3) = 0.2557011698983284_dp
    ros_A(4) = 0.3314825187068521E+01_dp
    ros_A(5) = 0.2896124015972201E+01_dp
    ros_A(6) = 0.9986419139977817_dp
    ros_A(7) = 0.1221224509226641E+01_dp
    ros_A(8) = 0.6019134481288629E+01_dp
    ros_A(9) = 0.1253708332932087E+02_dp
    ros_A(10) =-0.6878860361058950_dp
    ros_A(11) = ros_A(7)
    ros_A(12) = ros_A(8)
    ros_A(13) = ros_A(9)
    ros_A(14) = ros_A(10)
    ros_A(15) = 1.0_dp

    ros_C(1) =-0.5668800000000000E+01_dp
    ros_C(2) =-0.2430093356833875E+01_dp
    ros_C(3) =-0.2063599157091915_dp
    ros_C(4) =-0.1073529058151375_dp
    ros_C(5) =-0.9594562251023355E+01_dp
    ros_C(6) =-0.2047028614809616E+02_dp
    ros_C(7) = 0.7496443313967647E+01_dp
    ros_C(8) =-0.1024680431464352E+02_dp
    ros_C(9) =-0.3399990352819905E+02_dp
    ros_C(10) = 0.1170890893206160E+02_dp
    ros_C(11) = 0.8083246795921522E+01_dp
    ros_C(12) =-0.7981132988064893E+01_dp
    ros_C(13) =-0.3152159432874371E+02_dp
    ros_C(14) = 0.1631930543123136E+02_dp
    ros_C(15) =-0.6058818238834054E+01_dp

!~~~> M_i = Coefficients for new step solution
    ros_M(1) = ros_A(7)
    ros_M(2) = ros_A(8)
    ros_M(3) = ros_A(9)
    ros_M(4) = ros_A(10)
    ros_M(5) = 1.0_dp
    ros_M(6) = 1.0_dp

!~~~> E_i  = Coefficients for error estimator
    ros_E(1) = 0.0_dp
    ros_E(2) = 0.0_dp
    ros_E(3) = 0.0_dp
    ros_E(4) = 0.0_dp
    ros_E(5) = 0.0_dp
    ros_E(6) = 1.0_dp

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.
    ros_NewF(5) = .TRUE.
    ros_NewF(6) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 4.0_dp

  END SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! STIFFLY-STABLE W METHOD OF ORDER 3, WITH 4 STAGES
!
! J. RANG and L. ANGERMANN
! NEW ROSENBROCK W-METHODS OF ORDER 3
! FOR PARTIAL DIFFERENTIAL ALGEBRAIC
!        EQUATIONS OF INDEX 1                                             
! BIT Numerical Mathematics (2005) 45: 761-787
!  DOI: 10.1007/s10543-005-0035-y
! Table 4.1-4.2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RG3
!~~~> Name of the method
    ros_Name = 'RANG-3'
!~~~> Number of stages
    ros_S = 4

    ros_A(1) = 5.09052051067020d+00;
    ros_A(2) = 5.09052051067020d+00;
    ros_A(3) = 0.0d0;
    ros_A(4) = 4.97628111010787d+00;
    ros_A(5) = 2.77268164715849d-02;
    ros_A(6) = 2.29428036027904d-01;

    ros_C(1) = -1.16790812312283d+01;
    ros_C(2) = -1.64057326467367d+01;
    ros_C(3) = -2.77268164715850d-01;
    ros_C(4) = -8.38103960500476d+00;
    ros_C(5) = -8.48328409199343d-01;
    ros_C(6) =  2.87009860433106d-01;

    ros_M(1) =  5.22582761233094d+00;
    ros_M(2) = -5.56971148154165d-01;
    ros_M(3) =  3.57979469353645d-01;
    ros_M(4) =  1.72337398521064d+00;

    ros_E(1) = -5.16845212784040d+00;
    ros_E(2) = -1.26351942603842d+00;
    ros_E(3) = -1.11022302462516d-16;
    ros_E(4) =  2.22044604925031d-16;

    ros_Alpha(1) = 0.0d00;
    ros_Alpha(2) = 2.21878746765329d+00;
    ros_Alpha(3) = 2.21878746765329d+00;
    ros_Alpha(4) = 1.55392337535788d+00;

    ros_Gamma(1) =  4.35866521508459d-01;
    ros_Gamma(2) = -1.78292094614483d+00;
    ros_Gamma(3) = -2.46541900496934d+00;
    ros_Gamma(4) = -8.05529997906370d-01;

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 3.0_dp

  END SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END SUBROUTINE Rosenbrock_10


SUBROUTINE FunTemplate_10( T, Y, Ydot, LS_NSEL, LS_type )
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Function
 INTEGER,INTENT(IN)::LS_NSEL, LS_type
   REAL(kind=dp) :: T, Y(LS_NSEL)
   REAL(kind=dp) :: Ydot(LS_NSEL)
   REAL(kind=dp) :: Told
   Told = TIME
   TIME = T
   CALL Fun_10( Y, FIX, RCONST, Ydot, LS_NSEL)
   TIME = Told
END SUBROUTINE FunTemplate_10

SUBROUTINE JacTemplate_10( T, Y, Jcb, LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Jacobian
 USE gckpp_LinearAlgebra
 
 INTEGER,INTENT(IN)::LS_NSEL, LS_LU_NONZERO, LS_type
 INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO), LS_LU_ICOL(LS_LU_NONZERO)
    REAL(kind=dp) :: T, Y(LS_NSEL)
#ifdef FULL_ALGEBRA    
    REAL(kind=dp) :: JV(LS_LU_NONZERO), Jcb(LS_NSEL,LS_NSEL)
#else
    REAL(kind=dp) :: Jcb(LS_LU_NONZERO)
#endif   
    REAL(kind=dp) :: Told
#ifdef FULL_ALGEBRA    
    INTEGER :: i, j
#endif   

    Told = TIME
    TIME = T
#ifdef FULL_ALGEBRA        	
	CALL Jac_SP_10(Y, FIX, RCONST, JV,LS_NSEL,LS_LU_NONZERO)	
    DO j=1,LS_NSEL
      DO i=1,LS_NSEL
         Jcb(i,j) = 0.0_dp
      END DO
    END DO
    DO i=1,LS_LU_NONZERO
       Jcb(LS_LU_IROW(i),LS_LU_ICOL(i)) = JV(i)
    END DO
#else   
       CALL Jac_SP_10( Y, FIX, RCONST, Jcb,LS_NSEL,LS_LU_NONZERO)	 
#endif   
    TIME = Told
END SUBROUTINE JacTemplate_10



SUBROUTINE INTEGRATE_11( TIN, TOUT, LS_type,LS_NSEL, LS_NDEL,&
  Prate, Lrate, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )

   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: TIN  ! Start Time
   REAL(kind=dp), INTENT(IN) :: TOUT ! End Time
   
   ! Optional input parameters and statistics
   INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   REAL(kind=dp), INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   REAL(kind=dp), INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U
   INTEGER,		  INTENT(IN) :: LS_type,LS_NSEL,LS_NDEL
   REAL(kind=dp), INTENT(IN) :: Prate(NVAR),Lrate(NVAR)
   
   REAL(kind=dp) :: RCNTRL(20), RSTATUS(20),deltaT
   INTEGER       :: ICNTRL(20), ISTATUS(20), IERR
   REAL(kind=dp) :: VAR_selected(LS_NSEL),VAR_deleted(LS_NDEL),LS_P(LS_NDEL),LS_L(LS_NDEL)
   
   INTEGER, SAVE :: Ntotal = 0

   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.0_dp
   ISTATUS(:) = 0
   RSTATUS(:) = 0.0_dp
   deltaT =  TOUT-TIN
   !~~~> fine-tune the integrator:
   ICNTRL(1) = 0	! 0 - non-autonomous, 1 - autonomous
   ICNTRL(2) = 0	! 0 - vector tolerances, 1 - scalars

   ! If optional parameters are given, and if they are >0, 
   ! then they overwrite default settings. 
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF
   
        VAR_selected=VAR(select_ind_11)
	    VAR_deleted=VAR(delete_ind_11)
		LS_P=Prate(delete_ind_11)
		LS_L=Lrate(delete_ind_11)
        CALL Rosenbrock_11(LU_NSEL_11,VAR_selected,TIN,TOUT,ATOL,RTOL,&
          RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		  LU_NONZERO_11,LU_NSEL_11,LU_CROW_11,LU_DIAG_11,LU_IROW_11,LU_ICOL_11, LS_type)		
		WHERE(LS_L<=(0.01/deltaT))
			VAR_deleted=VAR_deleted+deltaT*(LS_P-LS_L*VAR_deleted)
		ELSEWHERE
			VAR_deleted=LS_P/LS_L+(VAR_deleted-LS_P/LS_L)*EXP(-LS_L*deltaT)
		END WHERE

	    VAR(select_ind_11)=VAR_selected
	    VAR(delete_ind_11)=VAR_deleted
   
   !~~~> Debug option: show no of steps
   ! Ntotal = Ntotal + ISTATUS(Nstp)
   ! PRINT*,'NSTEPS=',ISTATUS(Nstp),' (',Ntotal,')','  O3=', VAR(ind_O3)

   STEPMIN = RSTATUS(Nhexit)
   ! if optional parameters are given for output they 
   ! are updated with the return information
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE INTEGRATE_11

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE Rosenbrock_11(N,Y,Tstart,Tend, &
           AbsTol,RelTol,              &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		   LS_LU_NONZERO,LS_NSEL,LS_LU_CROW,LS_LU_DIAG,LS_LU_IROW,LS_LU_ICOL, LS_type)
  USE gckpp_Parameters
  USE gckpp_LinearAlgebra
  IMPLICIT NONE

   INTEGER,       INTENT(IN)    :: N,LS_LU_NONZERO,LS_NSEL,LS_type   
   INTEGER,INTENT(IN)::LS_LU_CROW(LS_NSEL+1),LS_LU_DIAG(LS_NSEL+1)
   INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO),LS_LU_ICOL(LS_LU_NONZERO)   
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
   REAL(kind=dp), INTENT(IN)    :: Tstart,Tend
   REAL(kind=dp), INTENT(IN)    :: AbsTol(N),RelTol(N)
   INTEGER,       INTENT(IN)    :: ICNTRL(20)
   REAL(kind=dp), INTENT(IN)    :: RCNTRL(20)
   INTEGER,       INTENT(INOUT) :: ISTATUS(20)
   REAL(kind=dp), INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT)   :: IERR
!~~~>  Parameters of the Rosenbrock method, up to 6 stages
   INTEGER ::  ros_S, rosMethod
   INTEGER, PARAMETER :: RS2=1, RS3=2, RS4=3, RD3=4, RD4=5, RG3=6
   REAL(kind=dp) :: ros_A(15), ros_C(15), ros_M(6), ros_E(6), &
                    ros_Alpha(6), ros_Gamma(6), ros_ELO
   LOGICAL :: ros_NewF(6)
   CHARACTER(LEN=12) :: ros_Name
!~~~>  Local variables
   REAL(kind=dp) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   REAL(kind=dp) :: Hmin, Hmax, Hstart
   REAL(kind=dp) :: Texit
   INTEGER       :: i, UplimTol, Max_no_steps
   LOGICAL       :: Autonomous, VectorTol
!~~~>   Parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp

!~~~>  Initialize statistics
   ISTATUS(1:8) = 0
   RSTATUS(1:3) = ZERO

!~~~>  Autonomous or time dependent ODE. Default is time dependent.
   Autonomous = .NOT.(ICNTRL(1) == 0)

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:N) and RelTol(1:N)
   IF (ICNTRL(2) == 0) THEN
      VectorTol = .TRUE.
      UplimTol  = N
   ELSE
      VectorTol = .FALSE.
      UplimTol  = 1
   END IF

!~~~>   Initialize the particular Rosenbrock method selected
   SELECT CASE (ICNTRL(3))
     CASE (1)
       CALL Ros2
     CASE (2)
       CALL Ros3
     CASE (3)
       CALL Ros4
     CASE (0,4)
       CALL Rodas3
     CASE (5)
       CALL Rodas4
     CASE (6)
       CALL Rang3
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(3)=',ICNTRL(3) 
       CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT

!~~~>   The maximum number of steps admitted
   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 200000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Unit roundoff (1+Roundoff>1)
   Roundoff = WLAMCH('E')

!~~~>  Lower bound on the step size: (positive value)
   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Upper bound on the step size: (positive value)
   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Starting step size: (positive value)
   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hold < FacMax
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Check if tolerances are reasonable
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO


!~~~>  CALL Rosenbrock method
   CALL ros_Integrator(Y, Tstart, Tend, Texit,   &
        AbsTol, RelTol,                          &
!  Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!  Error indicator
        IERR)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS !  SUBROUTINES internal to Rosenbrock
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 SUBROUTINE ros_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
   REAL(kind=dp), INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR
   
   IERR = Code
   PRINT * , &
     'Forced exit from Rosenbrock due to the following error:' 
     
   SELECT CASE (Code)
    CASE (-1)
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)    
      PRINT * , '--> Selected Rosenbrock method not implemented'
    CASE (-3)    
      PRINT * , '--> Hmin/Hmax/Hstart must be positive'
    CASE (-4)    
      PRINT * , '--> FacMin/FacMax/FacRej must be positive'
    CASE (-5) 
      PRINT * , '--> Improper tolerance values'
    CASE (-6) 
      PRINT * , '--> No of steps exceeds maximum bound'
    CASE (-7) 
      PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
    CASE (-8)    
      PRINT * , '--> Matrix is repeatedly singular'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT
   
   PRINT *, "T=", T, "and H=", H
     
 END SUBROUTINE ros_ErrorMsg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_Integrator (Y, Tstart, Tend, T,  &
        AbsTol, RelTol,                          &
!~~~> Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic Rosenbrock method
!      defined by ros_S (no of stages)
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

!~~~> Input: the initial condition at Tstart; Output: the solution at T
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
!~~~> Input: integration interval
   REAL(kind=dp), INTENT(IN) :: Tstart,Tend
!~~~> Output: time at which the solution is returned (T=Tend if success)
   REAL(kind=dp), INTENT(OUT) ::  T
!~~~> Input: tolerances
   REAL(kind=dp), INTENT(IN) ::  AbsTol(N), RelTol(N)
!~~~> Input: integration parameters
   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   REAL(kind=dp), INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   REAL(kind=dp), INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables
   REAL(kind=dp) :: Ynew(N), Fcn0(N), Fcn(N)
   REAL(kind=dp) :: K(N*ros_S), dFdT(N)
#ifdef FULL_ALGEBRA    
   REAL(kind=dp) :: Jac0(N,N), Ghimj(N,N)
#else
   REAL(kind=dp) :: Jac0(LS_LU_NONZERO), Ghimj(LS_LU_NONZERO)
#endif
   REAL(kind=dp) :: H, Hnew, HC, HG, Fac, Tau
   REAL(kind=dp) :: Err, Yerr(N)
   INTEGER :: Pivot(N), Direction, ioffset, j, istage
   LOGICAL :: RejectLastH, RejectMoreH, Singular
!~~~>  Local parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp
!~~~>  Locally called functions
!    REAL(kind=dp) WLAMCH
!    EXTERNAL WLAMCH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~>  Initial preparations
   T = Tstart
   RSTATUS(Nhexit) = ZERO
   H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   IF (ABS(H) <= 10.0_dp*Roundoff) H = DeltaMin

   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF
   H = Direction*H

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.

!~~~> Time loop begins below

TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) )

   IF ( ISTATUS(Nstp) > Max_no_steps ) THEN  ! Too many steps
      CALL ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
      CALL ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF

!~~~>  Limit H if necessary to avoid going beyond Tend
   H = MIN(H,ABS(Tend-T))

!~~~>   Compute the function at current time
   CALL FunTemplate_11(T,Y,Fcn0, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1

!~~~>  Compute the function derivative with respect to T
   IF (.NOT.Autonomous) THEN
      CALL ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
   END IF

!~~~>   Compute the Jacobian at current time
   CALL JacTemplate_11(T,Y,Jac0,LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)   
   ISTATUS(Njac) = ISTATUS(Njac) + 1

!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO

   CALL ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular)
   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       CALL ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S

      ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
       ioffset = N*(istage-1)

      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         !slim: CALL WCOPY(N,Fcn0,1,Fcn,1)
	 Fcn(1:N) = Fcn0(1:N)
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
         !slim: CALL WCOPY(N,Y,1,Ynew,1)
	 Ynew(1:N) = Y(1:N)
         DO j = 1, istage-1
           CALL WAXPY(N,ros_A((istage-1)*(istage-2)/2+j), &
            K(N*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL FunTemplate_11(Tau,Ynew,Fcn, LS_NSEL, LS_type)
         ISTATUS(Nfun) = ISTATUS(Nfun) + 1
       END IF ! if istage == 1 elseif ros_NewF(istage)
       !slim: CALL WCOPY(N,Fcn,1,K(ioffset+1),1)
       K(ioffset+1:ioffset+N) = Fcn(1:N)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL WAXPY(N,HC,K(N*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL WAXPY(N,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL ros_Solve(Ghimj, Pivot, K(ioffset+1))

   END DO Stage


!~~~>  Compute the new solution
   !slim: CALL WCOPY(N,Y,1,Ynew,1)
   Ynew(1:N) = Y(1:N)
   DO j=1,ros_S
         CALL WAXPY(N,ros_M(j),K(N*(j-1)+1),1,Ynew,1)
   END DO

!~~~>  Compute the error estimation
   !slim: CALL WSCAL(N,ZERO,Yerr,1)
   Yerr(1:N) = ZERO
   DO j=1,ros_S
        CALL WAXPY(N,ros_E(j),K(N*(j-1)+1),1,Yerr,1)
   END DO
   Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac

!~~~>  Check the error magnitude and adjust step size
   ISTATUS(Nstp) = ISTATUS(Nstp) + 1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
      ISTATUS(Nacc) = ISTATUS(Nacc) + 1
      !slim: CALL WCOPY(N,Ynew,1,Y,1)
      Y(1:N) = Ynew(1:N)
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H)
      END IF
      RSTATUS(Nhexit) = H
      RSTATUS(Nhnew)  = Hnew
      RSTATUS(Ntexit) = T
      RejectLastH = .FALSE.
      RejectMoreH = .FALSE.
      H = Hnew
      EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
   ELSE           !~~~> Reject step
      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (ISTATUS(Nacc) >= 1)  ISTATUS(Nrej) = ISTATUS(Nrej) + 1
   END IF ! Err <= 1

   END DO UntilAccepted

   
   END DO TimeLoop

!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful

  END SUBROUTINE ros_Integrator


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  REAL(kind=dp) FUNCTION ros_ErrorNorm ( Y, Ynew, Yerr, &
                               AbsTol, RelTol, VectorTol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

! Input arguments
   REAL(kind=dp), INTENT(IN) :: Y(N), Ynew(N), &
          Yerr(N), AbsTol(N), RelTol(N)
   LOGICAL, INTENT(IN) ::  VectorTol
! Local variables
   REAL(kind=dp) :: Err, Scale, Ymax
   INTEGER  :: i
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp

   Err = ZERO
   DO i=1,N
     Ymax = MAX(ABS(Y(i)),ABS(Ynew(i)))
     IF (VectorTol) THEN
       Scale = AbsTol(i)+RelTol(i)*Ymax
     ELSE
       Scale = AbsTol(1)+RelTol(1)*Ymax
     END IF
     Err = Err+(Yerr(i)/Scale)**2
   END DO
   Err  = SQRT(Err/N)

   ros_ErrorNorm = MAX(Err,1.0d-10)

  END FUNCTION ros_ErrorNorm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> The time partial derivative of the function by finite differences
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

!~~~> Input arguments
   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(N), Fcn0(N)
!~~~> Output arguments
   REAL(kind=dp), INTENT(OUT) :: dFdT(N)
!~~~> Local variables
   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL FunTemplate_11(T+Delta,Y,dFdT, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1
   CALL WAXPY(N,(-ONE),Fcn0,1,dFdT,1)
   CALL WSCAL(N,(ONE/Delta),dFdT,1)

  END SUBROUTINE ros_FunTimeDerivative


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_PrepareMatrix ( H, Direction, gam, &
             Jac0, Ghimj, Pivot, Singular )
! --- --- --- --- --- --- --- --- --- --- --- --- ---
!  Prepares the LHS matrix for stage calculations
!  1.  Construct Ghimj = 1/(H*ham) - Jac0
!      "(Gamma H) Inverse Minus Jacobian"
!  2.  Repeat LU decomposition of Ghimj until successful.
!       -half the step size if LU decomposition fails and retry
!       -exit after 5 consecutive fails
! --- --- --- --- --- --- --- --- --- --- --- --- ---
   IMPLICIT NONE

!~~~> Input arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) ::  Jac0(N,N)
#else
   REAL(kind=dp), INTENT(IN) ::  Jac0(LS_LU_NONZERO)
#endif   
   REAL(kind=dp), INTENT(IN) ::  gam
   INTEGER, INTENT(IN) ::  Direction
!~~~> Output arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(OUT) :: Ghimj(N,N)
#else
   REAL(kind=dp), INTENT(OUT) :: Ghimj(LS_LU_NONZERO)
#endif   
   LOGICAL, INTENT(OUT) ::  Singular
   INTEGER, INTENT(OUT) ::  Pivot(N)
!~~~> Inout arguments
   REAL(kind=dp), INTENT(INOUT) :: H   ! step size is decreased when LU fails
!~~~> Local variables
   INTEGER  :: i, ISING, Nconsecutive
   REAL(kind=dp) :: ghinv
   REAL(kind=dp), PARAMETER :: ONE  = 1.0_dp, HALF = 0.5_dp

   Nconsecutive = 0
   Singular = .TRUE.

   DO WHILE (Singular)

!~~~>    Construct Ghimj = 1/(H*gam) - Jac0
#ifdef FULL_ALGEBRA    
     !slim: CALL WCOPY(N*N,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(N*N,(-ONE),Ghimj,1)
     Ghimj = -Jac0
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(i,i) = Ghimj(i,i)+ghinv
     END DO
#else
     !slim: CALL WCOPY(LS_LU_NONZERO,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(LS_LU_NONZERO,(-ONE),Ghimj,1)
     Ghimj(1:LS_LU_NONZERO) = -Jac0(1:LS_LU_NONZERO)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(LS_LU_DIAG(i)) = Ghimj(LS_LU_DIAG(i))+ghinv
     END DO
#endif   
!~~~>    Compute LU decomposition
     CALL ros_Decomp( Ghimj, Pivot, ISING )
     IF (ISING == 0) THEN
!~~~>    If successful done
        Singular = .FALSE.
     ELSE ! ISING .ne. 0
!~~~>    If unsuccessful half the step size; if 5 consecutive fails then return
        ISTATUS(Nsng) = ISTATUS(Nsng) + 1
        Nconsecutive = Nconsecutive+1
        Singular = .TRUE.
        PRINT*,'Warning: LU Decomposition returned ISING = ',ISING
        IF (Nconsecutive <= 5) THEN ! Less than 5 consecutive failed decompositions
           H = H*HALF
        ELSE  ! More than 5 consecutive failed decompositions
           RETURN
        END IF  ! Nconsecutive
      END IF    ! ISING

   END DO ! WHILE Singular

  END SUBROUTINE ros_PrepareMatrix


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Decomp( A, Pivot, ISING )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the LU decomposition
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Inout variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(INOUT) :: A(N,N)
#else   
   REAL(kind=dp), INTENT(INOUT) :: A(LS_LU_NONZERO)
#endif
!~~~> Output variables
   INTEGER, INTENT(OUT) :: Pivot(N), ISING

#ifdef FULL_ALGEBRA    
   CALL  DGETRF( N, N, A, N, Pivot, ISING )
#else   
   CALL KppDecomp ( A, ISING,LS_NSEL,LS_LU_NONZERO,LS_LU_CROW,LS_LU_DIAG,LS_LU_ICOL)   
   Pivot(1) = 1
#endif
   ISTATUS(Ndec) = ISTATUS(Ndec) + 1

  END SUBROUTINE ros_Decomp


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Solve( A, Pivot, b)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the forward/backward substitution (using pre-computed LU decomposition)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Input variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) :: A(N,N)
   INTEGER :: ISING
#else   
   REAL(kind=dp), INTENT(IN) :: A(LS_LU_NONZERO)
#endif
   INTEGER, INTENT(IN) :: Pivot(N)
!~~~> InOut variables
   REAL(kind=dp), INTENT(INOUT) :: b(N)

#ifdef FULL_ALGEBRA    
   CALL  DGETRS( 'N', N , 1, A, N, Pivot, b, N, ISING )
   IF ( Info < 0 ) THEN
      PRINT*,"Error in DGETRS. ISING=",ISING
   END IF  
#else  
   	CALL KppSolve_11(LS_LU_NONZERO,LS_NSEL, A, b,LS_type )
#endif

   ISTATUS(Nsol) = ISTATUS(Nsol) + 1

  END SUBROUTINE ros_Solve



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   DOUBLE PRECISION g

    g = 1.0_dp + 1.0_dp/SQRT(2.0_dp)
    rosMethod = RS2
!~~~> Name of the method
    ros_Name = 'ROS-2'
!~~~> Number of stages
    ros_S = 2

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = (1.0_dp)/g
    ros_C(1) = (-2.0_dp)/g
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
!~~~> M_i = Coefficients for new step solution
    ros_M(1)= (3.0_dp)/(2.0_dp*g)
    ros_M(2)= (1.0_dp)/(2.0_dp*g)
! E_i = Coefficients for error estimator
    ros_E(1) = 1.0_dp/(2.0_dp*g)
    ros_E(2) = 1.0_dp/(2.0_dp*g)
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus one
    ros_ELO = 2.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.0_dp
    ros_Alpha(2) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = g
    ros_Gamma(2) =-g

 END SUBROUTINE Ros2


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   rosMethod = RS3
!~~~> Name of the method
   ros_Name = 'ROS-3'
!~~~> Number of stages
   ros_S = 3

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1)= 1.0_dp
   ros_A(2)= 1.0_dp
   ros_A(3)= 0.0_dp

   ros_C(1) = -0.10156171083877702091975600115545E+01_dp
   ros_C(2) =  0.40759956452537699824805835358067E+01_dp
   ros_C(3) =  0.92076794298330791242156818474003E+01_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) =  0.1E+01_dp
   ros_M(2) =  0.61697947043828245592553615689730E+01_dp
   ros_M(3) = -0.42772256543218573326238373806514_dp
! E_i = Coefficients for error estimator
   ros_E(1) =  0.5_dp
   ros_E(2) = -0.29079558716805469821718236208017E+01_dp
   ros_E(3) =  0.22354069897811569627360909276199_dp
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1)= 0.0_dp
   ros_Alpha(2)= 0.43586652150845899941601945119356_dp
   ros_Alpha(3)= 0.43586652150845899941601945119356_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1)= 0.43586652150845899941601945119356_dp
   ros_Gamma(2)= 0.24291996454816804366592249683314_dp
   ros_Gamma(3)= 0.21851380027664058511513169485832E+01_dp

  END SUBROUTINE Ros3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
!     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1990)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RS4
!~~~> Name of the method
   ros_Name = 'ROS-4'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.2000000000000000E+01_dp
   ros_A(2) = 0.1867943637803922E+01_dp
   ros_A(3) = 0.2344449711399156_dp
   ros_A(4) = ros_A(2)
   ros_A(5) = ros_A(3)
   ros_A(6) = 0.0_dp

   ros_C(1) =-0.7137615036412310E+01_dp
   ros_C(2) = 0.2580708087951457E+01_dp
   ros_C(3) = 0.6515950076447975_dp
   ros_C(4) =-0.2137148994382534E+01_dp
   ros_C(5) =-0.3214669691237626_dp
   ros_C(6) =-0.6949742501781779_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .TRUE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 0.2255570073418735E+01_dp
   ros_M(2) = 0.2870493262186792_dp
   ros_M(3) = 0.4353179431840180_dp
   ros_M(4) = 0.1093502252409163E+01_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) =-0.2815431932141155_dp
   ros_E(2) =-0.7276199124938920E-01_dp
   ros_E(3) =-0.1082196201495311_dp
   ros_E(4) =-0.1093502252409163E+01_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 4.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.1145640000000000E+01_dp
   ros_Alpha(3) = 0.6552168638155900_dp
   ros_Alpha(4) = ros_Alpha(3)
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5728200000000000_dp
   ros_Gamma(2) =-0.1769193891319233E+01_dp
   ros_Gamma(3) = 0.7592633437920482_dp
   ros_Gamma(4) =-0.1049021087100450_dp

  END SUBROUTINE Ros4

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- A STIFFLY-STABLE METHOD, 4 stages, order 3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RD3
!~~~> Name of the method
   ros_Name = 'RODAS-3'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.0_dp
   ros_A(2) = 2.0_dp
   ros_A(3) = 0.0_dp
   ros_A(4) = 2.0_dp
   ros_A(5) = 0.0_dp
   ros_A(6) = 1.0_dp

   ros_C(1) = 4.0_dp
   ros_C(2) = 1.0_dp
   ros_C(3) =-1.0_dp
   ros_C(4) = 1.0_dp
   ros_C(5) =-1.0_dp
   ros_C(6) =-(8.0_dp/3.0_dp)

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .FALSE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .TRUE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 2.0_dp
   ros_M(2) = 0.0_dp
   ros_M(3) = 1.0_dp
   ros_M(4) = 1.0_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) = 0.0_dp
   ros_E(2) = 0.0_dp
   ros_E(3) = 0.0_dp
   ros_E(4) = 1.0_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.0_dp
   ros_Alpha(3) = 1.0_dp
   ros_Alpha(4) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5_dp
   ros_Gamma(2) = 1.5_dp
   ros_Gamma(3) = 0.0_dp
   ros_Gamma(4) = 0.0_dp

  END SUBROUTINE Rodas3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1996)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RD4
!~~~> Name of the method
    ros_Name = 'RODAS-4'
!~~~> Number of stages
    ros_S = 6

!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.000_dp
    ros_Alpha(2) = 0.386_dp
    ros_Alpha(3) = 0.210_dp
    ros_Alpha(4) = 0.630_dp
    ros_Alpha(5) = 1.000_dp
    ros_Alpha(6) = 1.000_dp

!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = 0.2500000000000000_dp
    ros_Gamma(2) =-0.1043000000000000_dp
    ros_Gamma(3) = 0.1035000000000000_dp
    ros_Gamma(4) =-0.3620000000000023E-01_dp
    ros_Gamma(5) = 0.0_dp
    ros_Gamma(6) = 0.0_dp

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = 0.1544000000000000E+01_dp
    ros_A(2) = 0.9466785280815826_dp
    ros_A(3) = 0.2557011698983284_dp
    ros_A(4) = 0.3314825187068521E+01_dp
    ros_A(5) = 0.2896124015972201E+01_dp
    ros_A(6) = 0.9986419139977817_dp
    ros_A(7) = 0.1221224509226641E+01_dp
    ros_A(8) = 0.6019134481288629E+01_dp
    ros_A(9) = 0.1253708332932087E+02_dp
    ros_A(10) =-0.6878860361058950_dp
    ros_A(11) = ros_A(7)
    ros_A(12) = ros_A(8)
    ros_A(13) = ros_A(9)
    ros_A(14) = ros_A(10)
    ros_A(15) = 1.0_dp

    ros_C(1) =-0.5668800000000000E+01_dp
    ros_C(2) =-0.2430093356833875E+01_dp
    ros_C(3) =-0.2063599157091915_dp
    ros_C(4) =-0.1073529058151375_dp
    ros_C(5) =-0.9594562251023355E+01_dp
    ros_C(6) =-0.2047028614809616E+02_dp
    ros_C(7) = 0.7496443313967647E+01_dp
    ros_C(8) =-0.1024680431464352E+02_dp
    ros_C(9) =-0.3399990352819905E+02_dp
    ros_C(10) = 0.1170890893206160E+02_dp
    ros_C(11) = 0.8083246795921522E+01_dp
    ros_C(12) =-0.7981132988064893E+01_dp
    ros_C(13) =-0.3152159432874371E+02_dp
    ros_C(14) = 0.1631930543123136E+02_dp
    ros_C(15) =-0.6058818238834054E+01_dp

!~~~> M_i = Coefficients for new step solution
    ros_M(1) = ros_A(7)
    ros_M(2) = ros_A(8)
    ros_M(3) = ros_A(9)
    ros_M(4) = ros_A(10)
    ros_M(5) = 1.0_dp
    ros_M(6) = 1.0_dp

!~~~> E_i  = Coefficients for error estimator
    ros_E(1) = 0.0_dp
    ros_E(2) = 0.0_dp
    ros_E(3) = 0.0_dp
    ros_E(4) = 0.0_dp
    ros_E(5) = 0.0_dp
    ros_E(6) = 1.0_dp

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.
    ros_NewF(5) = .TRUE.
    ros_NewF(6) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 4.0_dp

  END SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! STIFFLY-STABLE W METHOD OF ORDER 3, WITH 4 STAGES
!
! J. RANG and L. ANGERMANN
! NEW ROSENBROCK W-METHODS OF ORDER 3
! FOR PARTIAL DIFFERENTIAL ALGEBRAIC
!        EQUATIONS OF INDEX 1                                             
! BIT Numerical Mathematics (2005) 45: 761-787
!  DOI: 10.1007/s10543-005-0035-y
! Table 4.1-4.2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RG3
!~~~> Name of the method
    ros_Name = 'RANG-3'
!~~~> Number of stages
    ros_S = 4

    ros_A(1) = 5.09052051067020d+00;
    ros_A(2) = 5.09052051067020d+00;
    ros_A(3) = 0.0d0;
    ros_A(4) = 4.97628111010787d+00;
    ros_A(5) = 2.77268164715849d-02;
    ros_A(6) = 2.29428036027904d-01;

    ros_C(1) = -1.16790812312283d+01;
    ros_C(2) = -1.64057326467367d+01;
    ros_C(3) = -2.77268164715850d-01;
    ros_C(4) = -8.38103960500476d+00;
    ros_C(5) = -8.48328409199343d-01;
    ros_C(6) =  2.87009860433106d-01;

    ros_M(1) =  5.22582761233094d+00;
    ros_M(2) = -5.56971148154165d-01;
    ros_M(3) =  3.57979469353645d-01;
    ros_M(4) =  1.72337398521064d+00;

    ros_E(1) = -5.16845212784040d+00;
    ros_E(2) = -1.26351942603842d+00;
    ros_E(3) = -1.11022302462516d-16;
    ros_E(4) =  2.22044604925031d-16;

    ros_Alpha(1) = 0.0d00;
    ros_Alpha(2) = 2.21878746765329d+00;
    ros_Alpha(3) = 2.21878746765329d+00;
    ros_Alpha(4) = 1.55392337535788d+00;

    ros_Gamma(1) =  4.35866521508459d-01;
    ros_Gamma(2) = -1.78292094614483d+00;
    ros_Gamma(3) = -2.46541900496934d+00;
    ros_Gamma(4) = -8.05529997906370d-01;

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 3.0_dp

  END SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END SUBROUTINE Rosenbrock_11


SUBROUTINE FunTemplate_11( T, Y, Ydot, LS_NSEL, LS_type )
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Function
 INTEGER,INTENT(IN)::LS_NSEL, LS_type
   REAL(kind=dp) :: T, Y(LS_NSEL)
   REAL(kind=dp) :: Ydot(LS_NSEL)
   REAL(kind=dp) :: Told
   Told = TIME
   TIME = T
   CALL Fun_11( Y, FIX, RCONST, Ydot, LS_NSEL)
   TIME = Told
END SUBROUTINE FunTemplate_11

SUBROUTINE JacTemplate_11( T, Y, Jcb, LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Jacobian
 USE gckpp_LinearAlgebra
 
 INTEGER,INTENT(IN)::LS_NSEL, LS_LU_NONZERO, LS_type
 INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO), LS_LU_ICOL(LS_LU_NONZERO)
    REAL(kind=dp) :: T, Y(LS_NSEL)
#ifdef FULL_ALGEBRA    
    REAL(kind=dp) :: JV(LS_LU_NONZERO), Jcb(LS_NSEL,LS_NSEL)
#else
    REAL(kind=dp) :: Jcb(LS_LU_NONZERO)
#endif   
    REAL(kind=dp) :: Told
#ifdef FULL_ALGEBRA    
    INTEGER :: i, j
#endif   

    Told = TIME
    TIME = T
#ifdef FULL_ALGEBRA        	
	CALL Jac_SP_11(Y, FIX, RCONST, JV,LS_NSEL,LS_LU_NONZERO)	
    DO j=1,LS_NSEL
      DO i=1,LS_NSEL
         Jcb(i,j) = 0.0_dp
      END DO
    END DO
    DO i=1,LS_LU_NONZERO
       Jcb(LS_LU_IROW(i),LS_LU_ICOL(i)) = JV(i)
    END DO
#else   
       CALL Jac_SP_11( Y, FIX, RCONST, Jcb,LS_NSEL,LS_LU_NONZERO)	 
#endif   
    TIME = Told
END SUBROUTINE JacTemplate_11



SUBROUTINE INTEGRATE_12( TIN, TOUT, LS_type,LS_NSEL, LS_NDEL,&
  Prate, Lrate, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )

   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: TIN  ! Start Time
   REAL(kind=dp), INTENT(IN) :: TOUT ! End Time
   
   ! Optional input parameters and statistics
   INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   REAL(kind=dp), INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   REAL(kind=dp), INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U
   INTEGER,		  INTENT(IN) :: LS_type,LS_NSEL,LS_NDEL
   REAL(kind=dp), INTENT(IN) :: Prate(NVAR),Lrate(NVAR)
   
   REAL(kind=dp) :: RCNTRL(20), RSTATUS(20),deltaT
   INTEGER       :: ICNTRL(20), ISTATUS(20), IERR
   REAL(kind=dp) :: VAR_selected(LS_NSEL),VAR_deleted(LS_NDEL),LS_P(LS_NDEL),LS_L(LS_NDEL)
   
   INTEGER, SAVE :: Ntotal = 0

   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.0_dp
   ISTATUS(:) = 0
   RSTATUS(:) = 0.0_dp
   deltaT =  TOUT-TIN
   !~~~> fine-tune the integrator:
   ICNTRL(1) = 0	! 0 - non-autonomous, 1 - autonomous
   ICNTRL(2) = 0	! 0 - vector tolerances, 1 - scalars

   ! If optional parameters are given, and if they are >0, 
   ! then they overwrite default settings. 
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF
   
        VAR_selected=VAR(select_ind_12)
	    VAR_deleted=VAR(delete_ind_12)
		LS_P=Prate(delete_ind_12)
		LS_L=Lrate(delete_ind_12)
        CALL Rosenbrock_12(LU_NSEL_12,VAR_selected,TIN,TOUT,ATOL,RTOL,&
          RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		  LU_NONZERO_12,LU_NSEL_12,LU_CROW_12,LU_DIAG_12,LU_IROW_12,LU_ICOL_12, LS_type)		
		WHERE(LS_L<=(0.01/deltaT))
			VAR_deleted=VAR_deleted+deltaT*(LS_P-LS_L*VAR_deleted)
		ELSEWHERE
			VAR_deleted=LS_P/LS_L+(VAR_deleted-LS_P/LS_L)*EXP(-LS_L*deltaT)
		END WHERE

	    VAR(select_ind_12)=VAR_selected
	    VAR(delete_ind_12)=VAR_deleted
   
   !~~~> Debug option: show no of steps
   ! Ntotal = Ntotal + ISTATUS(Nstp)
   ! PRINT*,'NSTEPS=',ISTATUS(Nstp),' (',Ntotal,')','  O3=', VAR(ind_O3)

   STEPMIN = RSTATUS(Nhexit)
   ! if optional parameters are given for output they 
   ! are updated with the return information
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE INTEGRATE_12

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE Rosenbrock_12(N,Y,Tstart,Tend, &
           AbsTol,RelTol,              &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		   LS_LU_NONZERO,LS_NSEL,LS_LU_CROW,LS_LU_DIAG,LS_LU_IROW,LS_LU_ICOL, LS_type)
  USE gckpp_Parameters
  USE gckpp_LinearAlgebra
  IMPLICIT NONE

   INTEGER,       INTENT(IN)    :: N,LS_LU_NONZERO,LS_NSEL,LS_type   
   INTEGER,INTENT(IN)::LS_LU_CROW(LS_NSEL+1),LS_LU_DIAG(LS_NSEL+1)
   INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO),LS_LU_ICOL(LS_LU_NONZERO)   
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
   REAL(kind=dp), INTENT(IN)    :: Tstart,Tend
   REAL(kind=dp), INTENT(IN)    :: AbsTol(N),RelTol(N)
   INTEGER,       INTENT(IN)    :: ICNTRL(20)
   REAL(kind=dp), INTENT(IN)    :: RCNTRL(20)
   INTEGER,       INTENT(INOUT) :: ISTATUS(20)
   REAL(kind=dp), INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT)   :: IERR
!~~~>  Parameters of the Rosenbrock method, up to 6 stages
   INTEGER ::  ros_S, rosMethod
   INTEGER, PARAMETER :: RS2=1, RS3=2, RS4=3, RD3=4, RD4=5, RG3=6
   REAL(kind=dp) :: ros_A(15), ros_C(15), ros_M(6), ros_E(6), &
                    ros_Alpha(6), ros_Gamma(6), ros_ELO
   LOGICAL :: ros_NewF(6)
   CHARACTER(LEN=12) :: ros_Name
!~~~>  Local variables
   REAL(kind=dp) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   REAL(kind=dp) :: Hmin, Hmax, Hstart
   REAL(kind=dp) :: Texit
   INTEGER       :: i, UplimTol, Max_no_steps
   LOGICAL       :: Autonomous, VectorTol
!~~~>   Parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp

!~~~>  Initialize statistics
   ISTATUS(1:8) = 0
   RSTATUS(1:3) = ZERO

!~~~>  Autonomous or time dependent ODE. Default is time dependent.
   Autonomous = .NOT.(ICNTRL(1) == 0)

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:N) and RelTol(1:N)
   IF (ICNTRL(2) == 0) THEN
      VectorTol = .TRUE.
      UplimTol  = N
   ELSE
      VectorTol = .FALSE.
      UplimTol  = 1
   END IF

!~~~>   Initialize the particular Rosenbrock method selected
   SELECT CASE (ICNTRL(3))
     CASE (1)
       CALL Ros2
     CASE (2)
       CALL Ros3
     CASE (3)
       CALL Ros4
     CASE (0,4)
       CALL Rodas3
     CASE (5)
       CALL Rodas4
     CASE (6)
       CALL Rang3
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(3)=',ICNTRL(3) 
       CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT

!~~~>   The maximum number of steps admitted
   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 200000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Unit roundoff (1+Roundoff>1)
   Roundoff = WLAMCH('E')

!~~~>  Lower bound on the step size: (positive value)
   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Upper bound on the step size: (positive value)
   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Starting step size: (positive value)
   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hold < FacMax
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Check if tolerances are reasonable
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO


!~~~>  CALL Rosenbrock method
   CALL ros_Integrator(Y, Tstart, Tend, Texit,   &
        AbsTol, RelTol,                          &
!  Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!  Error indicator
        IERR)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS !  SUBROUTINES internal to Rosenbrock
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 SUBROUTINE ros_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
   REAL(kind=dp), INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR
   
   IERR = Code
   PRINT * , &
     'Forced exit from Rosenbrock due to the following error:' 
     
   SELECT CASE (Code)
    CASE (-1)
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)    
      PRINT * , '--> Selected Rosenbrock method not implemented'
    CASE (-3)    
      PRINT * , '--> Hmin/Hmax/Hstart must be positive'
    CASE (-4)    
      PRINT * , '--> FacMin/FacMax/FacRej must be positive'
    CASE (-5) 
      PRINT * , '--> Improper tolerance values'
    CASE (-6) 
      PRINT * , '--> No of steps exceeds maximum bound'
    CASE (-7) 
      PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
    CASE (-8)    
      PRINT * , '--> Matrix is repeatedly singular'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT
   
   PRINT *, "T=", T, "and H=", H
     
 END SUBROUTINE ros_ErrorMsg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_Integrator (Y, Tstart, Tend, T,  &
        AbsTol, RelTol,                          &
!~~~> Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic Rosenbrock method
!      defined by ros_S (no of stages)
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

!~~~> Input: the initial condition at Tstart; Output: the solution at T
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
!~~~> Input: integration interval
   REAL(kind=dp), INTENT(IN) :: Tstart,Tend
!~~~> Output: time at which the solution is returned (T=Tend if success)
   REAL(kind=dp), INTENT(OUT) ::  T
!~~~> Input: tolerances
   REAL(kind=dp), INTENT(IN) ::  AbsTol(N), RelTol(N)
!~~~> Input: integration parameters
   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   REAL(kind=dp), INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   REAL(kind=dp), INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables
   REAL(kind=dp) :: Ynew(N), Fcn0(N), Fcn(N)
   REAL(kind=dp) :: K(N*ros_S), dFdT(N)
#ifdef FULL_ALGEBRA    
   REAL(kind=dp) :: Jac0(N,N), Ghimj(N,N)
#else
   REAL(kind=dp) :: Jac0(LS_LU_NONZERO), Ghimj(LS_LU_NONZERO)
#endif
   REAL(kind=dp) :: H, Hnew, HC, HG, Fac, Tau
   REAL(kind=dp) :: Err, Yerr(N)
   INTEGER :: Pivot(N), Direction, ioffset, j, istage
   LOGICAL :: RejectLastH, RejectMoreH, Singular
!~~~>  Local parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp
!~~~>  Locally called functions
!    REAL(kind=dp) WLAMCH
!    EXTERNAL WLAMCH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~>  Initial preparations
   T = Tstart
   RSTATUS(Nhexit) = ZERO
   H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   IF (ABS(H) <= 10.0_dp*Roundoff) H = DeltaMin

   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF
   H = Direction*H

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.

!~~~> Time loop begins below

TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) )

   IF ( ISTATUS(Nstp) > Max_no_steps ) THEN  ! Too many steps
      CALL ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
      CALL ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF

!~~~>  Limit H if necessary to avoid going beyond Tend
   H = MIN(H,ABS(Tend-T))

!~~~>   Compute the function at current time
   CALL FunTemplate_12(T,Y,Fcn0, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1

!~~~>  Compute the function derivative with respect to T
   IF (.NOT.Autonomous) THEN
      CALL ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
   END IF

!~~~>   Compute the Jacobian at current time
   CALL JacTemplate_12(T,Y,Jac0,LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)   
   ISTATUS(Njac) = ISTATUS(Njac) + 1

!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO

   CALL ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular)
   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       CALL ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S

      ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
       ioffset = N*(istage-1)

      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         !slim: CALL WCOPY(N,Fcn0,1,Fcn,1)
	 Fcn(1:N) = Fcn0(1:N)
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
         !slim: CALL WCOPY(N,Y,1,Ynew,1)
	 Ynew(1:N) = Y(1:N)
         DO j = 1, istage-1
           CALL WAXPY(N,ros_A((istage-1)*(istage-2)/2+j), &
            K(N*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL FunTemplate_12(Tau,Ynew,Fcn, LS_NSEL, LS_type)
         ISTATUS(Nfun) = ISTATUS(Nfun) + 1
       END IF ! if istage == 1 elseif ros_NewF(istage)
       !slim: CALL WCOPY(N,Fcn,1,K(ioffset+1),1)
       K(ioffset+1:ioffset+N) = Fcn(1:N)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL WAXPY(N,HC,K(N*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL WAXPY(N,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL ros_Solve(Ghimj, Pivot, K(ioffset+1))

   END DO Stage


!~~~>  Compute the new solution
   !slim: CALL WCOPY(N,Y,1,Ynew,1)
   Ynew(1:N) = Y(1:N)
   DO j=1,ros_S
         CALL WAXPY(N,ros_M(j),K(N*(j-1)+1),1,Ynew,1)
   END DO

!~~~>  Compute the error estimation
   !slim: CALL WSCAL(N,ZERO,Yerr,1)
   Yerr(1:N) = ZERO
   DO j=1,ros_S
        CALL WAXPY(N,ros_E(j),K(N*(j-1)+1),1,Yerr,1)
   END DO
   Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac

!~~~>  Check the error magnitude and adjust step size
   ISTATUS(Nstp) = ISTATUS(Nstp) + 1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
      ISTATUS(Nacc) = ISTATUS(Nacc) + 1
      !slim: CALL WCOPY(N,Ynew,1,Y,1)
      Y(1:N) = Ynew(1:N)
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H)
      END IF
      RSTATUS(Nhexit) = H
      RSTATUS(Nhnew)  = Hnew
      RSTATUS(Ntexit) = T
      RejectLastH = .FALSE.
      RejectMoreH = .FALSE.
      H = Hnew
      EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
   ELSE           !~~~> Reject step
      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (ISTATUS(Nacc) >= 1)  ISTATUS(Nrej) = ISTATUS(Nrej) + 1
   END IF ! Err <= 1

   END DO UntilAccepted

   
   END DO TimeLoop

!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful

  END SUBROUTINE ros_Integrator


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  REAL(kind=dp) FUNCTION ros_ErrorNorm ( Y, Ynew, Yerr, &
                               AbsTol, RelTol, VectorTol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

! Input arguments
   REAL(kind=dp), INTENT(IN) :: Y(N), Ynew(N), &
          Yerr(N), AbsTol(N), RelTol(N)
   LOGICAL, INTENT(IN) ::  VectorTol
! Local variables
   REAL(kind=dp) :: Err, Scale, Ymax
   INTEGER  :: i
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp

   Err = ZERO
   DO i=1,N
     Ymax = MAX(ABS(Y(i)),ABS(Ynew(i)))
     IF (VectorTol) THEN
       Scale = AbsTol(i)+RelTol(i)*Ymax
     ELSE
       Scale = AbsTol(1)+RelTol(1)*Ymax
     END IF
     Err = Err+(Yerr(i)/Scale)**2
   END DO
   Err  = SQRT(Err/N)

   ros_ErrorNorm = MAX(Err,1.0d-10)

  END FUNCTION ros_ErrorNorm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> The time partial derivative of the function by finite differences
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

!~~~> Input arguments
   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(N), Fcn0(N)
!~~~> Output arguments
   REAL(kind=dp), INTENT(OUT) :: dFdT(N)
!~~~> Local variables
   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL FunTemplate_12(T+Delta,Y,dFdT, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1
   CALL WAXPY(N,(-ONE),Fcn0,1,dFdT,1)
   CALL WSCAL(N,(ONE/Delta),dFdT,1)

  END SUBROUTINE ros_FunTimeDerivative


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_PrepareMatrix ( H, Direction, gam, &
             Jac0, Ghimj, Pivot, Singular )
! --- --- --- --- --- --- --- --- --- --- --- --- ---
!  Prepares the LHS matrix for stage calculations
!  1.  Construct Ghimj = 1/(H*ham) - Jac0
!      "(Gamma H) Inverse Minus Jacobian"
!  2.  Repeat LU decomposition of Ghimj until successful.
!       -half the step size if LU decomposition fails and retry
!       -exit after 5 consecutive fails
! --- --- --- --- --- --- --- --- --- --- --- --- ---
   IMPLICIT NONE

!~~~> Input arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) ::  Jac0(N,N)
#else
   REAL(kind=dp), INTENT(IN) ::  Jac0(LS_LU_NONZERO)
#endif   
   REAL(kind=dp), INTENT(IN) ::  gam
   INTEGER, INTENT(IN) ::  Direction
!~~~> Output arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(OUT) :: Ghimj(N,N)
#else
   REAL(kind=dp), INTENT(OUT) :: Ghimj(LS_LU_NONZERO)
#endif   
   LOGICAL, INTENT(OUT) ::  Singular
   INTEGER, INTENT(OUT) ::  Pivot(N)
!~~~> Inout arguments
   REAL(kind=dp), INTENT(INOUT) :: H   ! step size is decreased when LU fails
!~~~> Local variables
   INTEGER  :: i, ISING, Nconsecutive
   REAL(kind=dp) :: ghinv
   REAL(kind=dp), PARAMETER :: ONE  = 1.0_dp, HALF = 0.5_dp

   Nconsecutive = 0
   Singular = .TRUE.

   DO WHILE (Singular)

!~~~>    Construct Ghimj = 1/(H*gam) - Jac0
#ifdef FULL_ALGEBRA    
     !slim: CALL WCOPY(N*N,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(N*N,(-ONE),Ghimj,1)
     Ghimj = -Jac0
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(i,i) = Ghimj(i,i)+ghinv
     END DO
#else
     !slim: CALL WCOPY(LS_LU_NONZERO,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(LS_LU_NONZERO,(-ONE),Ghimj,1)
     Ghimj(1:LS_LU_NONZERO) = -Jac0(1:LS_LU_NONZERO)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(LS_LU_DIAG(i)) = Ghimj(LS_LU_DIAG(i))+ghinv
     END DO
#endif   
!~~~>    Compute LU decomposition
     CALL ros_Decomp( Ghimj, Pivot, ISING )
     IF (ISING == 0) THEN
!~~~>    If successful done
        Singular = .FALSE.
     ELSE ! ISING .ne. 0
!~~~>    If unsuccessful half the step size; if 5 consecutive fails then return
        ISTATUS(Nsng) = ISTATUS(Nsng) + 1
        Nconsecutive = Nconsecutive+1
        Singular = .TRUE.
        PRINT*,'Warning: LU Decomposition returned ISING = ',ISING
        IF (Nconsecutive <= 5) THEN ! Less than 5 consecutive failed decompositions
           H = H*HALF
        ELSE  ! More than 5 consecutive failed decompositions
           RETURN
        END IF  ! Nconsecutive
      END IF    ! ISING

   END DO ! WHILE Singular

  END SUBROUTINE ros_PrepareMatrix


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Decomp( A, Pivot, ISING )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the LU decomposition
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Inout variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(INOUT) :: A(N,N)
#else   
   REAL(kind=dp), INTENT(INOUT) :: A(LS_LU_NONZERO)
#endif
!~~~> Output variables
   INTEGER, INTENT(OUT) :: Pivot(N), ISING

#ifdef FULL_ALGEBRA    
   CALL  DGETRF( N, N, A, N, Pivot, ISING )
#else   
   CALL KppDecomp ( A, ISING,LS_NSEL,LS_LU_NONZERO,LS_LU_CROW,LS_LU_DIAG,LS_LU_ICOL)   
   Pivot(1) = 1
#endif
   ISTATUS(Ndec) = ISTATUS(Ndec) + 1

  END SUBROUTINE ros_Decomp


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Solve( A, Pivot, b)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the forward/backward substitution (using pre-computed LU decomposition)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Input variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) :: A(N,N)
   INTEGER :: ISING
#else   
   REAL(kind=dp), INTENT(IN) :: A(LS_LU_NONZERO)
#endif
   INTEGER, INTENT(IN) :: Pivot(N)
!~~~> InOut variables
   REAL(kind=dp), INTENT(INOUT) :: b(N)

#ifdef FULL_ALGEBRA    
   CALL  DGETRS( 'N', N , 1, A, N, Pivot, b, N, ISING )
   IF ( Info < 0 ) THEN
      PRINT*,"Error in DGETRS. ISING=",ISING
   END IF  
#else  
   	CALL KppSolve_12(LS_LU_NONZERO,LS_NSEL, A, b,LS_type )
#endif

   ISTATUS(Nsol) = ISTATUS(Nsol) + 1

  END SUBROUTINE ros_Solve



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   DOUBLE PRECISION g

    g = 1.0_dp + 1.0_dp/SQRT(2.0_dp)
    rosMethod = RS2
!~~~> Name of the method
    ros_Name = 'ROS-2'
!~~~> Number of stages
    ros_S = 2

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = (1.0_dp)/g
    ros_C(1) = (-2.0_dp)/g
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
!~~~> M_i = Coefficients for new step solution
    ros_M(1)= (3.0_dp)/(2.0_dp*g)
    ros_M(2)= (1.0_dp)/(2.0_dp*g)
! E_i = Coefficients for error estimator
    ros_E(1) = 1.0_dp/(2.0_dp*g)
    ros_E(2) = 1.0_dp/(2.0_dp*g)
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus one
    ros_ELO = 2.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.0_dp
    ros_Alpha(2) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = g
    ros_Gamma(2) =-g

 END SUBROUTINE Ros2


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   rosMethod = RS3
!~~~> Name of the method
   ros_Name = 'ROS-3'
!~~~> Number of stages
   ros_S = 3

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1)= 1.0_dp
   ros_A(2)= 1.0_dp
   ros_A(3)= 0.0_dp

   ros_C(1) = -0.10156171083877702091975600115545E+01_dp
   ros_C(2) =  0.40759956452537699824805835358067E+01_dp
   ros_C(3) =  0.92076794298330791242156818474003E+01_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) =  0.1E+01_dp
   ros_M(2) =  0.61697947043828245592553615689730E+01_dp
   ros_M(3) = -0.42772256543218573326238373806514_dp
! E_i = Coefficients for error estimator
   ros_E(1) =  0.5_dp
   ros_E(2) = -0.29079558716805469821718236208017E+01_dp
   ros_E(3) =  0.22354069897811569627360909276199_dp
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1)= 0.0_dp
   ros_Alpha(2)= 0.43586652150845899941601945119356_dp
   ros_Alpha(3)= 0.43586652150845899941601945119356_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1)= 0.43586652150845899941601945119356_dp
   ros_Gamma(2)= 0.24291996454816804366592249683314_dp
   ros_Gamma(3)= 0.21851380027664058511513169485832E+01_dp

  END SUBROUTINE Ros3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
!     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1990)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RS4
!~~~> Name of the method
   ros_Name = 'ROS-4'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.2000000000000000E+01_dp
   ros_A(2) = 0.1867943637803922E+01_dp
   ros_A(3) = 0.2344449711399156_dp
   ros_A(4) = ros_A(2)
   ros_A(5) = ros_A(3)
   ros_A(6) = 0.0_dp

   ros_C(1) =-0.7137615036412310E+01_dp
   ros_C(2) = 0.2580708087951457E+01_dp
   ros_C(3) = 0.6515950076447975_dp
   ros_C(4) =-0.2137148994382534E+01_dp
   ros_C(5) =-0.3214669691237626_dp
   ros_C(6) =-0.6949742501781779_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .TRUE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 0.2255570073418735E+01_dp
   ros_M(2) = 0.2870493262186792_dp
   ros_M(3) = 0.4353179431840180_dp
   ros_M(4) = 0.1093502252409163E+01_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) =-0.2815431932141155_dp
   ros_E(2) =-0.7276199124938920E-01_dp
   ros_E(3) =-0.1082196201495311_dp
   ros_E(4) =-0.1093502252409163E+01_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 4.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.1145640000000000E+01_dp
   ros_Alpha(3) = 0.6552168638155900_dp
   ros_Alpha(4) = ros_Alpha(3)
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5728200000000000_dp
   ros_Gamma(2) =-0.1769193891319233E+01_dp
   ros_Gamma(3) = 0.7592633437920482_dp
   ros_Gamma(4) =-0.1049021087100450_dp

  END SUBROUTINE Ros4

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- A STIFFLY-STABLE METHOD, 4 stages, order 3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RD3
!~~~> Name of the method
   ros_Name = 'RODAS-3'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.0_dp
   ros_A(2) = 2.0_dp
   ros_A(3) = 0.0_dp
   ros_A(4) = 2.0_dp
   ros_A(5) = 0.0_dp
   ros_A(6) = 1.0_dp

   ros_C(1) = 4.0_dp
   ros_C(2) = 1.0_dp
   ros_C(3) =-1.0_dp
   ros_C(4) = 1.0_dp
   ros_C(5) =-1.0_dp
   ros_C(6) =-(8.0_dp/3.0_dp)

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .FALSE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .TRUE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 2.0_dp
   ros_M(2) = 0.0_dp
   ros_M(3) = 1.0_dp
   ros_M(4) = 1.0_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) = 0.0_dp
   ros_E(2) = 0.0_dp
   ros_E(3) = 0.0_dp
   ros_E(4) = 1.0_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.0_dp
   ros_Alpha(3) = 1.0_dp
   ros_Alpha(4) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5_dp
   ros_Gamma(2) = 1.5_dp
   ros_Gamma(3) = 0.0_dp
   ros_Gamma(4) = 0.0_dp

  END SUBROUTINE Rodas3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1996)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RD4
!~~~> Name of the method
    ros_Name = 'RODAS-4'
!~~~> Number of stages
    ros_S = 6

!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.000_dp
    ros_Alpha(2) = 0.386_dp
    ros_Alpha(3) = 0.210_dp
    ros_Alpha(4) = 0.630_dp
    ros_Alpha(5) = 1.000_dp
    ros_Alpha(6) = 1.000_dp

!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = 0.2500000000000000_dp
    ros_Gamma(2) =-0.1043000000000000_dp
    ros_Gamma(3) = 0.1035000000000000_dp
    ros_Gamma(4) =-0.3620000000000023E-01_dp
    ros_Gamma(5) = 0.0_dp
    ros_Gamma(6) = 0.0_dp

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = 0.1544000000000000E+01_dp
    ros_A(2) = 0.9466785280815826_dp
    ros_A(3) = 0.2557011698983284_dp
    ros_A(4) = 0.3314825187068521E+01_dp
    ros_A(5) = 0.2896124015972201E+01_dp
    ros_A(6) = 0.9986419139977817_dp
    ros_A(7) = 0.1221224509226641E+01_dp
    ros_A(8) = 0.6019134481288629E+01_dp
    ros_A(9) = 0.1253708332932087E+02_dp
    ros_A(10) =-0.6878860361058950_dp
    ros_A(11) = ros_A(7)
    ros_A(12) = ros_A(8)
    ros_A(13) = ros_A(9)
    ros_A(14) = ros_A(10)
    ros_A(15) = 1.0_dp

    ros_C(1) =-0.5668800000000000E+01_dp
    ros_C(2) =-0.2430093356833875E+01_dp
    ros_C(3) =-0.2063599157091915_dp
    ros_C(4) =-0.1073529058151375_dp
    ros_C(5) =-0.9594562251023355E+01_dp
    ros_C(6) =-0.2047028614809616E+02_dp
    ros_C(7) = 0.7496443313967647E+01_dp
    ros_C(8) =-0.1024680431464352E+02_dp
    ros_C(9) =-0.3399990352819905E+02_dp
    ros_C(10) = 0.1170890893206160E+02_dp
    ros_C(11) = 0.8083246795921522E+01_dp
    ros_C(12) =-0.7981132988064893E+01_dp
    ros_C(13) =-0.3152159432874371E+02_dp
    ros_C(14) = 0.1631930543123136E+02_dp
    ros_C(15) =-0.6058818238834054E+01_dp

!~~~> M_i = Coefficients for new step solution
    ros_M(1) = ros_A(7)
    ros_M(2) = ros_A(8)
    ros_M(3) = ros_A(9)
    ros_M(4) = ros_A(10)
    ros_M(5) = 1.0_dp
    ros_M(6) = 1.0_dp

!~~~> E_i  = Coefficients for error estimator
    ros_E(1) = 0.0_dp
    ros_E(2) = 0.0_dp
    ros_E(3) = 0.0_dp
    ros_E(4) = 0.0_dp
    ros_E(5) = 0.0_dp
    ros_E(6) = 1.0_dp

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.
    ros_NewF(5) = .TRUE.
    ros_NewF(6) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 4.0_dp

  END SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! STIFFLY-STABLE W METHOD OF ORDER 3, WITH 4 STAGES
!
! J. RANG and L. ANGERMANN
! NEW ROSENBROCK W-METHODS OF ORDER 3
! FOR PARTIAL DIFFERENTIAL ALGEBRAIC
!        EQUATIONS OF INDEX 1                                             
! BIT Numerical Mathematics (2005) 45: 761-787
!  DOI: 10.1007/s10543-005-0035-y
! Table 4.1-4.2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RG3
!~~~> Name of the method
    ros_Name = 'RANG-3'
!~~~> Number of stages
    ros_S = 4

    ros_A(1) = 5.09052051067020d+00;
    ros_A(2) = 5.09052051067020d+00;
    ros_A(3) = 0.0d0;
    ros_A(4) = 4.97628111010787d+00;
    ros_A(5) = 2.77268164715849d-02;
    ros_A(6) = 2.29428036027904d-01;

    ros_C(1) = -1.16790812312283d+01;
    ros_C(2) = -1.64057326467367d+01;
    ros_C(3) = -2.77268164715850d-01;
    ros_C(4) = -8.38103960500476d+00;
    ros_C(5) = -8.48328409199343d-01;
    ros_C(6) =  2.87009860433106d-01;

    ros_M(1) =  5.22582761233094d+00;
    ros_M(2) = -5.56971148154165d-01;
    ros_M(3) =  3.57979469353645d-01;
    ros_M(4) =  1.72337398521064d+00;

    ros_E(1) = -5.16845212784040d+00;
    ros_E(2) = -1.26351942603842d+00;
    ros_E(3) = -1.11022302462516d-16;
    ros_E(4) =  2.22044604925031d-16;

    ros_Alpha(1) = 0.0d00;
    ros_Alpha(2) = 2.21878746765329d+00;
    ros_Alpha(3) = 2.21878746765329d+00;
    ros_Alpha(4) = 1.55392337535788d+00;

    ros_Gamma(1) =  4.35866521508459d-01;
    ros_Gamma(2) = -1.78292094614483d+00;
    ros_Gamma(3) = -2.46541900496934d+00;
    ros_Gamma(4) = -8.05529997906370d-01;

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 3.0_dp

  END SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END SUBROUTINE Rosenbrock_12


SUBROUTINE FunTemplate_12( T, Y, Ydot, LS_NSEL, LS_type )
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Function
 INTEGER,INTENT(IN)::LS_NSEL, LS_type
   REAL(kind=dp) :: T, Y(LS_NSEL)
   REAL(kind=dp) :: Ydot(LS_NSEL)
   REAL(kind=dp) :: Told
   Told = TIME
   TIME = T
   CALL Fun_12( Y, FIX, RCONST, Ydot, LS_NSEL)
   TIME = Told
END SUBROUTINE FunTemplate_12

SUBROUTINE JacTemplate_12( T, Y, Jcb, LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Jacobian
 USE gckpp_LinearAlgebra
 
 INTEGER,INTENT(IN)::LS_NSEL, LS_LU_NONZERO, LS_type
 INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO), LS_LU_ICOL(LS_LU_NONZERO)
    REAL(kind=dp) :: T, Y(LS_NSEL)
#ifdef FULL_ALGEBRA    
    REAL(kind=dp) :: JV(LS_LU_NONZERO), Jcb(LS_NSEL,LS_NSEL)
#else
    REAL(kind=dp) :: Jcb(LS_LU_NONZERO)
#endif   
    REAL(kind=dp) :: Told
#ifdef FULL_ALGEBRA    
    INTEGER :: i, j
#endif   

    Told = TIME
    TIME = T
#ifdef FULL_ALGEBRA        	
	CALL Jac_SP_12(Y, FIX, RCONST, JV,LS_NSEL,LS_LU_NONZERO)	
    DO j=1,LS_NSEL
      DO i=1,LS_NSEL
         Jcb(i,j) = 0.0_dp
      END DO
    END DO
    DO i=1,LS_LU_NONZERO
       Jcb(LS_LU_IROW(i),LS_LU_ICOL(i)) = JV(i)
    END DO
#else   
       CALL Jac_SP_12( Y, FIX, RCONST, Jcb,LS_NSEL,LS_LU_NONZERO)	 
#endif   
    TIME = Told
END SUBROUTINE JacTemplate_12



SUBROUTINE INTEGRATE_13( TIN, TOUT, LS_type,LS_NSEL, LS_NDEL,&
  Prate, Lrate, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )

   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: TIN  ! Start Time
   REAL(kind=dp), INTENT(IN) :: TOUT ! End Time
   
   ! Optional input parameters and statistics
   INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   REAL(kind=dp), INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   REAL(kind=dp), INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U
   INTEGER,		  INTENT(IN) :: LS_type,LS_NSEL,LS_NDEL
   REAL(kind=dp), INTENT(IN) :: Prate(NVAR),Lrate(NVAR)
   
   REAL(kind=dp) :: RCNTRL(20), RSTATUS(20),deltaT
   INTEGER       :: ICNTRL(20), ISTATUS(20), IERR
   
   
   INTEGER, SAVE :: Ntotal = 0

   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.0_dp
   ISTATUS(:) = 0
   RSTATUS(:) = 0.0_dp
   deltaT =  TOUT-TIN
   !~~~> fine-tune the integrator:
   ICNTRL(1) = 0	! 0 - non-autonomous, 1 - autonomous
   ICNTRL(2) = 0	! 0 - vector tolerances, 1 - scalars

   ! If optional parameters are given, and if they are >0, 
   ! then they overwrite default settings. 
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF
   

        CALL Rosenbrock_13(LU_NSEL_13,VAR,TIN,TOUT,ATOL,RTOL,&
          RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		  LU_NONZERO_13,LU_NSEL_13,LU_CROW_13,LU_DIAG_13,LU_IROW_13,LU_ICOL_13, LS_type	
   
   !~~~> Debug option: show no of steps
   ! Ntotal = Ntotal + ISTATUS(Nstp)
   ! PRINT*,'NSTEPS=',ISTATUS(Nstp),' (',Ntotal,')','  O3=', VAR(ind_O3)

   STEPMIN = RSTATUS(Nhexit)
   ! if optional parameters are given for output they 
   ! are updated with the return information
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE INTEGRATE_13

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE Rosenbrock_13(N,Y,Tstart,Tend, &
           AbsTol,RelTol,              &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		   LS_LU_NONZERO,LS_NSEL,LS_LU_CROW,LS_LU_DIAG,LS_LU_IROW,LS_LU_ICOL, LS_type)
  USE gckpp_Parameters
  USE gckpp_LinearAlgebra
  IMPLICIT NONE

   INTEGER,       INTENT(IN)    :: N,LS_LU_NONZERO,LS_NSEL,LS_type   
   INTEGER,INTENT(IN)::LS_LU_CROW(LS_NSEL+1),LS_LU_DIAG(LS_NSEL+1)
   INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO),LS_LU_ICOL(LS_LU_NONZERO)   
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
   REAL(kind=dp), INTENT(IN)    :: Tstart,Tend
   REAL(kind=dp), INTENT(IN)    :: AbsTol(N),RelTol(N)
   INTEGER,       INTENT(IN)    :: ICNTRL(20)
   REAL(kind=dp), INTENT(IN)    :: RCNTRL(20)
   INTEGER,       INTENT(INOUT) :: ISTATUS(20)
   REAL(kind=dp), INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT)   :: IERR
!~~~>  Parameters of the Rosenbrock method, up to 6 stages
   INTEGER ::  ros_S, rosMethod
   INTEGER, PARAMETER :: RS2=1, RS3=2, RS4=3, RD3=4, RD4=5, RG3=6
   REAL(kind=dp) :: ros_A(15), ros_C(15), ros_M(6), ros_E(6), &
                    ros_Alpha(6), ros_Gamma(6), ros_ELO
   LOGICAL :: ros_NewF(6)
   CHARACTER(LEN=12) :: ros_Name
!~~~>  Local variables
   REAL(kind=dp) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   REAL(kind=dp) :: Hmin, Hmax, Hstart
   REAL(kind=dp) :: Texit
   INTEGER       :: i, UplimTol, Max_no_steps
   LOGICAL       :: Autonomous, VectorTol
!~~~>   Parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp

!~~~>  Initialize statistics
   ISTATUS(1:8) = 0
   RSTATUS(1:3) = ZERO

!~~~>  Autonomous or time dependent ODE. Default is time dependent.
   Autonomous = .NOT.(ICNTRL(1) == 0)

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:N) and RelTol(1:N)
   IF (ICNTRL(2) == 0) THEN
      VectorTol = .TRUE.
      UplimTol  = N
   ELSE
      VectorTol = .FALSE.
      UplimTol  = 1
   END IF

!~~~>   Initialize the particular Rosenbrock method selected
   SELECT CASE (ICNTRL(3))
     CASE (1)
       CALL Ros2
     CASE (2)
       CALL Ros3
     CASE (3)
       CALL Ros4
     CASE (0,4)
       CALL Rodas3
     CASE (5)
       CALL Rodas4
     CASE (6)
       CALL Rang3
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(3)=',ICNTRL(3) 
       CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT

!~~~>   The maximum number of steps admitted
   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 200000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Unit roundoff (1+Roundoff>1)
   Roundoff = WLAMCH('E')

!~~~>  Lower bound on the step size: (positive value)
   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Upper bound on the step size: (positive value)
   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Starting step size: (positive value)
   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hold < FacMax
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Check if tolerances are reasonable
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO


!~~~>  CALL Rosenbrock method
   CALL ros_Integrator(Y, Tstart, Tend, Texit,   &
        AbsTol, RelTol,                          &
!  Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!  Error indicator
        IERR)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS !  SUBROUTINES internal to Rosenbrock
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 SUBROUTINE ros_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
   REAL(kind=dp), INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR
   
   IERR = Code
   PRINT * , &
     'Forced exit from Rosenbrock due to the following error:' 
     
   SELECT CASE (Code)
    CASE (-1)
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)    
      PRINT * , '--> Selected Rosenbrock method not implemented'
    CASE (-3)    
      PRINT * , '--> Hmin/Hmax/Hstart must be positive'
    CASE (-4)    
      PRINT * , '--> FacMin/FacMax/FacRej must be positive'
    CASE (-5) 
      PRINT * , '--> Improper tolerance values'
    CASE (-6) 
      PRINT * , '--> No of steps exceeds maximum bound'
    CASE (-7) 
      PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
    CASE (-8)    
      PRINT * , '--> Matrix is repeatedly singular'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT
   
   PRINT *, "T=", T, "and H=", H
     
 END SUBROUTINE ros_ErrorMsg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_Integrator (Y, Tstart, Tend, T,  &
        AbsTol, RelTol,                          &
!~~~> Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic Rosenbrock method
!      defined by ros_S (no of stages)
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

!~~~> Input: the initial condition at Tstart; Output: the solution at T
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
!~~~> Input: integration interval
   REAL(kind=dp), INTENT(IN) :: Tstart,Tend
!~~~> Output: time at which the solution is returned (T=Tend if success)
   REAL(kind=dp), INTENT(OUT) ::  T
!~~~> Input: tolerances
   REAL(kind=dp), INTENT(IN) ::  AbsTol(N), RelTol(N)
!~~~> Input: integration parameters
   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   REAL(kind=dp), INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   REAL(kind=dp), INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables
   REAL(kind=dp) :: Ynew(N), Fcn0(N), Fcn(N)
   REAL(kind=dp) :: K(N*ros_S), dFdT(N)
#ifdef FULL_ALGEBRA    
   REAL(kind=dp) :: Jac0(N,N), Ghimj(N,N)
#else
   REAL(kind=dp) :: Jac0(LS_LU_NONZERO), Ghimj(LS_LU_NONZERO)
#endif
   REAL(kind=dp) :: H, Hnew, HC, HG, Fac, Tau
   REAL(kind=dp) :: Err, Yerr(N)
   INTEGER :: Pivot(N), Direction, ioffset, j, istage
   LOGICAL :: RejectLastH, RejectMoreH, Singular
!~~~>  Local parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp
!~~~>  Locally called functions
!    REAL(kind=dp) WLAMCH
!    EXTERNAL WLAMCH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~>  Initial preparations
   T = Tstart
   RSTATUS(Nhexit) = ZERO
   H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   IF (ABS(H) <= 10.0_dp*Roundoff) H = DeltaMin

   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF
   H = Direction*H

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.

!~~~> Time loop begins below

TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) )

   IF ( ISTATUS(Nstp) > Max_no_steps ) THEN  ! Too many steps
      CALL ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
      CALL ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF

!~~~>  Limit H if necessary to avoid going beyond Tend
   H = MIN(H,ABS(Tend-T))

!~~~>   Compute the function at current time
   CALL FunTemplate_13(T,Y,Fcn0, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1

!~~~>  Compute the function derivative with respect to T
   IF (.NOT.Autonomous) THEN
      CALL ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
   END IF

!~~~>   Compute the Jacobian at current time
   CALL JacTemplate_13(T,Y,Jac0,LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)   
   ISTATUS(Njac) = ISTATUS(Njac) + 1

!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO

   CALL ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular)
   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       CALL ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S

      ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
       ioffset = N*(istage-1)

      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         !slim: CALL WCOPY(N,Fcn0,1,Fcn,1)
	 Fcn(1:N) = Fcn0(1:N)
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
         !slim: CALL WCOPY(N,Y,1,Ynew,1)
	 Ynew(1:N) = Y(1:N)
         DO j = 1, istage-1
           CALL WAXPY(N,ros_A((istage-1)*(istage-2)/2+j), &
            K(N*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL FunTemplate_13(Tau,Ynew,Fcn, LS_NSEL, LS_type)
         ISTATUS(Nfun) = ISTATUS(Nfun) + 1
       END IF ! if istage == 1 elseif ros_NewF(istage)
       !slim: CALL WCOPY(N,Fcn,1,K(ioffset+1),1)
       K(ioffset+1:ioffset+N) = Fcn(1:N)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL WAXPY(N,HC,K(N*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL WAXPY(N,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL ros_Solve(Ghimj, Pivot, K(ioffset+1))

   END DO Stage


!~~~>  Compute the new solution
   !slim: CALL WCOPY(N,Y,1,Ynew,1)
   Ynew(1:N) = Y(1:N)
   DO j=1,ros_S
         CALL WAXPY(N,ros_M(j),K(N*(j-1)+1),1,Ynew,1)
   END DO

!~~~>  Compute the error estimation
   !slim: CALL WSCAL(N,ZERO,Yerr,1)
   Yerr(1:N) = ZERO
   DO j=1,ros_S
        CALL WAXPY(N,ros_E(j),K(N*(j-1)+1),1,Yerr,1)
   END DO
   Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac

!~~~>  Check the error magnitude and adjust step size
   ISTATUS(Nstp) = ISTATUS(Nstp) + 1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
      ISTATUS(Nacc) = ISTATUS(Nacc) + 1
      !slim: CALL WCOPY(N,Ynew,1,Y,1)
      Y(1:N) = Ynew(1:N)
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H)
      END IF
      RSTATUS(Nhexit) = H
      RSTATUS(Nhnew)  = Hnew
      RSTATUS(Ntexit) = T
      RejectLastH = .FALSE.
      RejectMoreH = .FALSE.
      H = Hnew
      EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
   ELSE           !~~~> Reject step
      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (ISTATUS(Nacc) >= 1)  ISTATUS(Nrej) = ISTATUS(Nrej) + 1
   END IF ! Err <= 1

   END DO UntilAccepted

   
   END DO TimeLoop

!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful

  END SUBROUTINE ros_Integrator


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  REAL(kind=dp) FUNCTION ros_ErrorNorm ( Y, Ynew, Yerr, &
                               AbsTol, RelTol, VectorTol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

! Input arguments
   REAL(kind=dp), INTENT(IN) :: Y(N), Ynew(N), &
          Yerr(N), AbsTol(N), RelTol(N)
   LOGICAL, INTENT(IN) ::  VectorTol
! Local variables
   REAL(kind=dp) :: Err, Scale, Ymax
   INTEGER  :: i
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp

   Err = ZERO
   DO i=1,N
     Ymax = MAX(ABS(Y(i)),ABS(Ynew(i)))
     IF (VectorTol) THEN
       Scale = AbsTol(i)+RelTol(i)*Ymax
     ELSE
       Scale = AbsTol(1)+RelTol(1)*Ymax
     END IF
     Err = Err+(Yerr(i)/Scale)**2
   END DO
   Err  = SQRT(Err/N)

   ros_ErrorNorm = MAX(Err,1.0d-10)

  END FUNCTION ros_ErrorNorm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> The time partial derivative of the function by finite differences
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

!~~~> Input arguments
   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(N), Fcn0(N)
!~~~> Output arguments
   REAL(kind=dp), INTENT(OUT) :: dFdT(N)
!~~~> Local variables
   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL FunTemplate_13(T+Delta,Y,dFdT, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1
   CALL WAXPY(N,(-ONE),Fcn0,1,dFdT,1)
   CALL WSCAL(N,(ONE/Delta),dFdT,1)

  END SUBROUTINE ros_FunTimeDerivative


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_PrepareMatrix ( H, Direction, gam, &
             Jac0, Ghimj, Pivot, Singular )
! --- --- --- --- --- --- --- --- --- --- --- --- ---
!  Prepares the LHS matrix for stage calculations
!  1.  Construct Ghimj = 1/(H*ham) - Jac0
!      "(Gamma H) Inverse Minus Jacobian"
!  2.  Repeat LU decomposition of Ghimj until successful.
!       -half the step size if LU decomposition fails and retry
!       -exit after 5 consecutive fails
! --- --- --- --- --- --- --- --- --- --- --- --- ---
   IMPLICIT NONE

!~~~> Input arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) ::  Jac0(N,N)
#else
   REAL(kind=dp), INTENT(IN) ::  Jac0(LS_LU_NONZERO)
#endif   
   REAL(kind=dp), INTENT(IN) ::  gam
   INTEGER, INTENT(IN) ::  Direction
!~~~> Output arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(OUT) :: Ghimj(N,N)
#else
   REAL(kind=dp), INTENT(OUT) :: Ghimj(LS_LU_NONZERO)
#endif   
   LOGICAL, INTENT(OUT) ::  Singular
   INTEGER, INTENT(OUT) ::  Pivot(N)
!~~~> Inout arguments
   REAL(kind=dp), INTENT(INOUT) :: H   ! step size is decreased when LU fails
!~~~> Local variables
   INTEGER  :: i, ISING, Nconsecutive
   REAL(kind=dp) :: ghinv
   REAL(kind=dp), PARAMETER :: ONE  = 1.0_dp, HALF = 0.5_dp

   Nconsecutive = 0
   Singular = .TRUE.

   DO WHILE (Singular)

!~~~>    Construct Ghimj = 1/(H*gam) - Jac0
#ifdef FULL_ALGEBRA    
     !slim: CALL WCOPY(N*N,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(N*N,(-ONE),Ghimj,1)
     Ghimj = -Jac0
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(i,i) = Ghimj(i,i)+ghinv
     END DO
#else
     !slim: CALL WCOPY(LS_LU_NONZERO,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(LS_LU_NONZERO,(-ONE),Ghimj,1)
     Ghimj(1:LS_LU_NONZERO) = -Jac0(1:LS_LU_NONZERO)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(LS_LU_DIAG(i)) = Ghimj(LS_LU_DIAG(i))+ghinv
     END DO
#endif   
!~~~>    Compute LU decomposition
     CALL ros_Decomp( Ghimj, Pivot, ISING )
     IF (ISING == 0) THEN
!~~~>    If successful done
        Singular = .FALSE.
     ELSE ! ISING .ne. 0
!~~~>    If unsuccessful half the step size; if 5 consecutive fails then return
        ISTATUS(Nsng) = ISTATUS(Nsng) + 1
        Nconsecutive = Nconsecutive+1
        Singular = .TRUE.
        PRINT*,'Warning: LU Decomposition returned ISING = ',ISING
        IF (Nconsecutive <= 5) THEN ! Less than 5 consecutive failed decompositions
           H = H*HALF
        ELSE  ! More than 5 consecutive failed decompositions
           RETURN
        END IF  ! Nconsecutive
      END IF    ! ISING

   END DO ! WHILE Singular

  END SUBROUTINE ros_PrepareMatrix


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Decomp( A, Pivot, ISING )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the LU decomposition
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Inout variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(INOUT) :: A(N,N)
#else   
   REAL(kind=dp), INTENT(INOUT) :: A(LS_LU_NONZERO)
#endif
!~~~> Output variables
   INTEGER, INTENT(OUT) :: Pivot(N), ISING

#ifdef FULL_ALGEBRA    
   CALL  DGETRF( N, N, A, N, Pivot, ISING )
#else   
   CALL KppDecomp ( A, ISING,LS_NSEL,LS_LU_NONZERO,LS_LU_CROW,LS_LU_DIAG,LS_LU_ICOL)   
   Pivot(1) = 1
#endif
   ISTATUS(Ndec) = ISTATUS(Ndec) + 1

  END SUBROUTINE ros_Decomp


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Solve( A, Pivot, b)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the forward/backward substitution (using pre-computed LU decomposition)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Input variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) :: A(N,N)
   INTEGER :: ISING
#else   
   REAL(kind=dp), INTENT(IN) :: A(LS_LU_NONZERO)
#endif
   INTEGER, INTENT(IN) :: Pivot(N)
!~~~> InOut variables
   REAL(kind=dp), INTENT(INOUT) :: b(N)

#ifdef FULL_ALGEBRA    
   CALL  DGETRS( 'N', N , 1, A, N, Pivot, b, N, ISING )
   IF ( Info < 0 ) THEN
      PRINT*,"Error in DGETRS. ISING=",ISING
   END IF  
#else  
   	CALL KppSolve_13(LS_LU_NONZERO,LS_NSEL, A, b,LS_type )
#endif

   ISTATUS(Nsol) = ISTATUS(Nsol) + 1

  END SUBROUTINE ros_Solve



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   DOUBLE PRECISION g

    g = 1.0_dp + 1.0_dp/SQRT(2.0_dp)
    rosMethod = RS2
!~~~> Name of the method
    ros_Name = 'ROS-2'
!~~~> Number of stages
    ros_S = 2

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = (1.0_dp)/g
    ros_C(1) = (-2.0_dp)/g
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
!~~~> M_i = Coefficients for new step solution
    ros_M(1)= (3.0_dp)/(2.0_dp*g)
    ros_M(2)= (1.0_dp)/(2.0_dp*g)
! E_i = Coefficients for error estimator
    ros_E(1) = 1.0_dp/(2.0_dp*g)
    ros_E(2) = 1.0_dp/(2.0_dp*g)
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus one
    ros_ELO = 2.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.0_dp
    ros_Alpha(2) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = g
    ros_Gamma(2) =-g

 END SUBROUTINE Ros2


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   rosMethod = RS3
!~~~> Name of the method
   ros_Name = 'ROS-3'
!~~~> Number of stages
   ros_S = 3

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1)= 1.0_dp
   ros_A(2)= 1.0_dp
   ros_A(3)= 0.0_dp

   ros_C(1) = -0.10156171083877702091975600115545E+01_dp
   ros_C(2) =  0.40759956452537699824805835358067E+01_dp
   ros_C(3) =  0.92076794298330791242156818474003E+01_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) =  0.1E+01_dp
   ros_M(2) =  0.61697947043828245592553615689730E+01_dp
   ros_M(3) = -0.42772256543218573326238373806514_dp
! E_i = Coefficients for error estimator
   ros_E(1) =  0.5_dp
   ros_E(2) = -0.29079558716805469821718236208017E+01_dp
   ros_E(3) =  0.22354069897811569627360909276199_dp
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1)= 0.0_dp
   ros_Alpha(2)= 0.43586652150845899941601945119356_dp
   ros_Alpha(3)= 0.43586652150845899941601945119356_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1)= 0.43586652150845899941601945119356_dp
   ros_Gamma(2)= 0.24291996454816804366592249683314_dp
   ros_Gamma(3)= 0.21851380027664058511513169485832E+01_dp

  END SUBROUTINE Ros3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
!     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1990)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RS4
!~~~> Name of the method
   ros_Name = 'ROS-4'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.2000000000000000E+01_dp
   ros_A(2) = 0.1867943637803922E+01_dp
   ros_A(3) = 0.2344449711399156_dp
   ros_A(4) = ros_A(2)
   ros_A(5) = ros_A(3)
   ros_A(6) = 0.0_dp

   ros_C(1) =-0.7137615036412310E+01_dp
   ros_C(2) = 0.2580708087951457E+01_dp
   ros_C(3) = 0.6515950076447975_dp
   ros_C(4) =-0.2137148994382534E+01_dp
   ros_C(5) =-0.3214669691237626_dp
   ros_C(6) =-0.6949742501781779_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .TRUE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 0.2255570073418735E+01_dp
   ros_M(2) = 0.2870493262186792_dp
   ros_M(3) = 0.4353179431840180_dp
   ros_M(4) = 0.1093502252409163E+01_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) =-0.2815431932141155_dp
   ros_E(2) =-0.7276199124938920E-01_dp
   ros_E(3) =-0.1082196201495311_dp
   ros_E(4) =-0.1093502252409163E+01_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 4.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.1145640000000000E+01_dp
   ros_Alpha(3) = 0.6552168638155900_dp
   ros_Alpha(4) = ros_Alpha(3)
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5728200000000000_dp
   ros_Gamma(2) =-0.1769193891319233E+01_dp
   ros_Gamma(3) = 0.7592633437920482_dp
   ros_Gamma(4) =-0.1049021087100450_dp

  END SUBROUTINE Ros4

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- A STIFFLY-STABLE METHOD, 4 stages, order 3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RD3
!~~~> Name of the method
   ros_Name = 'RODAS-3'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.0_dp
   ros_A(2) = 2.0_dp
   ros_A(3) = 0.0_dp
   ros_A(4) = 2.0_dp
   ros_A(5) = 0.0_dp
   ros_A(6) = 1.0_dp

   ros_C(1) = 4.0_dp
   ros_C(2) = 1.0_dp
   ros_C(3) =-1.0_dp
   ros_C(4) = 1.0_dp
   ros_C(5) =-1.0_dp
   ros_C(6) =-(8.0_dp/3.0_dp)

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .FALSE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .TRUE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 2.0_dp
   ros_M(2) = 0.0_dp
   ros_M(3) = 1.0_dp
   ros_M(4) = 1.0_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) = 0.0_dp
   ros_E(2) = 0.0_dp
   ros_E(3) = 0.0_dp
   ros_E(4) = 1.0_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.0_dp
   ros_Alpha(3) = 1.0_dp
   ros_Alpha(4) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5_dp
   ros_Gamma(2) = 1.5_dp
   ros_Gamma(3) = 0.0_dp
   ros_Gamma(4) = 0.0_dp

  END SUBROUTINE Rodas3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1996)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RD4
!~~~> Name of the method
    ros_Name = 'RODAS-4'
!~~~> Number of stages
    ros_S = 6

!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.000_dp
    ros_Alpha(2) = 0.386_dp
    ros_Alpha(3) = 0.210_dp
    ros_Alpha(4) = 0.630_dp
    ros_Alpha(5) = 1.000_dp
    ros_Alpha(6) = 1.000_dp

!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = 0.2500000000000000_dp
    ros_Gamma(2) =-0.1043000000000000_dp
    ros_Gamma(3) = 0.1035000000000000_dp
    ros_Gamma(4) =-0.3620000000000023E-01_dp
    ros_Gamma(5) = 0.0_dp
    ros_Gamma(6) = 0.0_dp

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = 0.1544000000000000E+01_dp
    ros_A(2) = 0.9466785280815826_dp
    ros_A(3) = 0.2557011698983284_dp
    ros_A(4) = 0.3314825187068521E+01_dp
    ros_A(5) = 0.2896124015972201E+01_dp
    ros_A(6) = 0.9986419139977817_dp
    ros_A(7) = 0.1221224509226641E+01_dp
    ros_A(8) = 0.6019134481288629E+01_dp
    ros_A(9) = 0.1253708332932087E+02_dp
    ros_A(10) =-0.6878860361058950_dp
    ros_A(11) = ros_A(7)
    ros_A(12) = ros_A(8)
    ros_A(13) = ros_A(9)
    ros_A(14) = ros_A(10)
    ros_A(15) = 1.0_dp

    ros_C(1) =-0.5668800000000000E+01_dp
    ros_C(2) =-0.2430093356833875E+01_dp
    ros_C(3) =-0.2063599157091915_dp
    ros_C(4) =-0.1073529058151375_dp
    ros_C(5) =-0.9594562251023355E+01_dp
    ros_C(6) =-0.2047028614809616E+02_dp
    ros_C(7) = 0.7496443313967647E+01_dp
    ros_C(8) =-0.1024680431464352E+02_dp
    ros_C(9) =-0.3399990352819905E+02_dp
    ros_C(10) = 0.1170890893206160E+02_dp
    ros_C(11) = 0.8083246795921522E+01_dp
    ros_C(12) =-0.7981132988064893E+01_dp
    ros_C(13) =-0.3152159432874371E+02_dp
    ros_C(14) = 0.1631930543123136E+02_dp
    ros_C(15) =-0.6058818238834054E+01_dp

!~~~> M_i = Coefficients for new step solution
    ros_M(1) = ros_A(7)
    ros_M(2) = ros_A(8)
    ros_M(3) = ros_A(9)
    ros_M(4) = ros_A(10)
    ros_M(5) = 1.0_dp
    ros_M(6) = 1.0_dp

!~~~> E_i  = Coefficients for error estimator
    ros_E(1) = 0.0_dp
    ros_E(2) = 0.0_dp
    ros_E(3) = 0.0_dp
    ros_E(4) = 0.0_dp
    ros_E(5) = 0.0_dp
    ros_E(6) = 1.0_dp

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.
    ros_NewF(5) = .TRUE.
    ros_NewF(6) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 4.0_dp

  END SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! STIFFLY-STABLE W METHOD OF ORDER 3, WITH 4 STAGES
!
! J. RANG and L. ANGERMANN
! NEW ROSENBROCK W-METHODS OF ORDER 3
! FOR PARTIAL DIFFERENTIAL ALGEBRAIC
!        EQUATIONS OF INDEX 1                                             
! BIT Numerical Mathematics (2005) 45: 761-787
!  DOI: 10.1007/s10543-005-0035-y
! Table 4.1-4.2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RG3
!~~~> Name of the method
    ros_Name = 'RANG-3'
!~~~> Number of stages
    ros_S = 4

    ros_A(1) = 5.09052051067020d+00;
    ros_A(2) = 5.09052051067020d+00;
    ros_A(3) = 0.0d0;
    ros_A(4) = 4.97628111010787d+00;
    ros_A(5) = 2.77268164715849d-02;
    ros_A(6) = 2.29428036027904d-01;

    ros_C(1) = -1.16790812312283d+01;
    ros_C(2) = -1.64057326467367d+01;
    ros_C(3) = -2.77268164715850d-01;
    ros_C(4) = -8.38103960500476d+00;
    ros_C(5) = -8.48328409199343d-01;
    ros_C(6) =  2.87009860433106d-01;

    ros_M(1) =  5.22582761233094d+00;
    ros_M(2) = -5.56971148154165d-01;
    ros_M(3) =  3.57979469353645d-01;
    ros_M(4) =  1.72337398521064d+00;

    ros_E(1) = -5.16845212784040d+00;
    ros_E(2) = -1.26351942603842d+00;
    ros_E(3) = -1.11022302462516d-16;
    ros_E(4) =  2.22044604925031d-16;

    ros_Alpha(1) = 0.0d00;
    ros_Alpha(2) = 2.21878746765329d+00;
    ros_Alpha(3) = 2.21878746765329d+00;
    ros_Alpha(4) = 1.55392337535788d+00;

    ros_Gamma(1) =  4.35866521508459d-01;
    ros_Gamma(2) = -1.78292094614483d+00;
    ros_Gamma(3) = -2.46541900496934d+00;
    ros_Gamma(4) = -8.05529997906370d-01;

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 3.0_dp

  END SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END SUBROUTINE Rosenbrock_13


SUBROUTINE FunTemplate_13( T, Y, Ydot, LS_NSEL, LS_type )
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Function
 INTEGER,INTENT(IN)::LS_NSEL, LS_type
   REAL(kind=dp) :: T, Y(LS_NSEL)
   REAL(kind=dp) :: Ydot(LS_NSEL)
   REAL(kind=dp) :: Told
   Told = TIME
   TIME = T
   CALL Fun_13( Y, FIX, RCONST, Ydot, LS_NSEL)
   TIME = Told
END SUBROUTINE FunTemplate_13

SUBROUTINE JacTemplate_13( T, Y, Jcb, LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Jacobian
 USE gckpp_LinearAlgebra
 
 INTEGER,INTENT(IN)::LS_NSEL, LS_LU_NONZERO, LS_type
 INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO), LS_LU_ICOL(LS_LU_NONZERO)
    REAL(kind=dp) :: T, Y(LS_NSEL)
#ifdef FULL_ALGEBRA    
    REAL(kind=dp) :: JV(LS_LU_NONZERO), Jcb(LS_NSEL,LS_NSEL)
#else
    REAL(kind=dp) :: Jcb(LS_LU_NONZERO)
#endif   
    REAL(kind=dp) :: Told
#ifdef FULL_ALGEBRA    
    INTEGER :: i, j
#endif   

    Told = TIME
    TIME = T
#ifdef FULL_ALGEBRA        	
	CALL Jac_SP_13(Y, FIX, RCONST, JV,LS_NSEL,LS_LU_NONZERO)	
    DO j=1,LS_NSEL
      DO i=1,LS_NSEL
         Jcb(i,j) = 0.0_dp
      END DO
    END DO
    DO i=1,LS_LU_NONZERO
       Jcb(LS_LU_IROW(i),LS_LU_ICOL(i)) = JV(i)
    END DO
#else   
       CALL Jac_SP_13( Y, FIX, RCONST, Jcb,LS_NSEL,LS_LU_NONZERO)	 
#endif   
    TIME = Told
END SUBROUTINE JacTemplate_13



SUBROUTINE INTEGRATE_14( TIN, TOUT, LS_type,LS_NSEL, LS_NDEL,&
  Prate, Lrate, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )

   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: TIN  ! Start Time
   REAL(kind=dp), INTENT(IN) :: TOUT ! End Time
   
   ! Optional input parameters and statistics
   INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   REAL(kind=dp), INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   REAL(kind=dp), INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U
   INTEGER,		  INTENT(IN) :: LS_type,LS_NSEL,LS_NDEL
   REAL(kind=dp), INTENT(IN) :: Prate(NVAR),Lrate(NVAR)
   
   REAL(kind=dp) :: RCNTRL(20), RSTATUS(20),deltaT
   INTEGER       :: ICNTRL(20), ISTATUS(20), IERR
   REAL(kind=dp) :: VAR_selected(LS_NSEL),VAR_deleted(LS_NDEL),LS_P(LS_NDEL),LS_L(LS_NDEL)
   
   INTEGER, SAVE :: Ntotal = 0

   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.0_dp
   ISTATUS(:) = 0
   RSTATUS(:) = 0.0_dp
   deltaT =  TOUT-TIN
   !~~~> fine-tune the integrator:
   ICNTRL(1) = 0	! 0 - non-autonomous, 1 - autonomous
   ICNTRL(2) = 0	! 0 - vector tolerances, 1 - scalars

   ! If optional parameters are given, and if they are >0, 
   ! then they overwrite default settings. 
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF
   
        VAR_selected=VAR(select_ind_14)
	    VAR_deleted=VAR(delete_ind_14)
		LS_P=Prate(delete_ind_14)
		LS_L=Lrate(delete_ind_14)
        CALL Rosenbrock_14(LU_NSEL_14,VAR_selected,TIN,TOUT,ATOL,RTOL,&
          RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		  LU_NONZERO_14,LU_NSEL_14,LU_CROW_14,LU_DIAG_14,LU_IROW_14,LU_ICOL_14, LS_type)		
		WHERE(LS_L<=(0.01/deltaT))
			VAR_deleted=VAR_deleted+deltaT*(LS_P-LS_L*VAR_deleted)
		ELSEWHERE
			VAR_deleted=LS_P/LS_L+(VAR_deleted-LS_P/LS_L)*EXP(-LS_L*deltaT)
		END WHERE

	    VAR(select_ind_14)=VAR_selected
	    VAR(delete_ind_14)=VAR_deleted
   
   !~~~> Debug option: show no of steps
   ! Ntotal = Ntotal + ISTATUS(Nstp)
   ! PRINT*,'NSTEPS=',ISTATUS(Nstp),' (',Ntotal,')','  O3=', VAR(ind_O3)

   STEPMIN = RSTATUS(Nhexit)
   ! if optional parameters are given for output they 
   ! are updated with the return information
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE INTEGRATE_14

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE Rosenbrock_14(N,Y,Tstart,Tend, &
           AbsTol,RelTol,              &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		   LS_LU_NONZERO,LS_NSEL,LS_LU_CROW,LS_LU_DIAG,LS_LU_IROW,LS_LU_ICOL, LS_type)
  USE gckpp_Parameters
  USE gckpp_LinearAlgebra
  IMPLICIT NONE

   INTEGER,       INTENT(IN)    :: N,LS_LU_NONZERO,LS_NSEL,LS_type   
   INTEGER,INTENT(IN)::LS_LU_CROW(LS_NSEL+1),LS_LU_DIAG(LS_NSEL+1)
   INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO),LS_LU_ICOL(LS_LU_NONZERO)   
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
   REAL(kind=dp), INTENT(IN)    :: Tstart,Tend
   REAL(kind=dp), INTENT(IN)    :: AbsTol(N),RelTol(N)
   INTEGER,       INTENT(IN)    :: ICNTRL(20)
   REAL(kind=dp), INTENT(IN)    :: RCNTRL(20)
   INTEGER,       INTENT(INOUT) :: ISTATUS(20)
   REAL(kind=dp), INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT)   :: IERR
!~~~>  Parameters of the Rosenbrock method, up to 6 stages
   INTEGER ::  ros_S, rosMethod
   INTEGER, PARAMETER :: RS2=1, RS3=2, RS4=3, RD3=4, RD4=5, RG3=6
   REAL(kind=dp) :: ros_A(15), ros_C(15), ros_M(6), ros_E(6), &
                    ros_Alpha(6), ros_Gamma(6), ros_ELO
   LOGICAL :: ros_NewF(6)
   CHARACTER(LEN=12) :: ros_Name
!~~~>  Local variables
   REAL(kind=dp) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   REAL(kind=dp) :: Hmin, Hmax, Hstart
   REAL(kind=dp) :: Texit
   INTEGER       :: i, UplimTol, Max_no_steps
   LOGICAL       :: Autonomous, VectorTol
!~~~>   Parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp

!~~~>  Initialize statistics
   ISTATUS(1:8) = 0
   RSTATUS(1:3) = ZERO

!~~~>  Autonomous or time dependent ODE. Default is time dependent.
   Autonomous = .NOT.(ICNTRL(1) == 0)

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:N) and RelTol(1:N)
   IF (ICNTRL(2) == 0) THEN
      VectorTol = .TRUE.
      UplimTol  = N
   ELSE
      VectorTol = .FALSE.
      UplimTol  = 1
   END IF

!~~~>   Initialize the particular Rosenbrock method selected
   SELECT CASE (ICNTRL(3))
     CASE (1)
       CALL Ros2
     CASE (2)
       CALL Ros3
     CASE (3)
       CALL Ros4
     CASE (0,4)
       CALL Rodas3
     CASE (5)
       CALL Rodas4
     CASE (6)
       CALL Rang3
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(3)=',ICNTRL(3) 
       CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT

!~~~>   The maximum number of steps admitted
   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 200000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Unit roundoff (1+Roundoff>1)
   Roundoff = WLAMCH('E')

!~~~>  Lower bound on the step size: (positive value)
   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Upper bound on the step size: (positive value)
   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Starting step size: (positive value)
   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hold < FacMax
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Check if tolerances are reasonable
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO


!~~~>  CALL Rosenbrock method
   CALL ros_Integrator(Y, Tstart, Tend, Texit,   &
        AbsTol, RelTol,                          &
!  Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!  Error indicator
        IERR)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS !  SUBROUTINES internal to Rosenbrock
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 SUBROUTINE ros_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
   REAL(kind=dp), INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR
   
   IERR = Code
   PRINT * , &
     'Forced exit from Rosenbrock due to the following error:' 
     
   SELECT CASE (Code)
    CASE (-1)
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)    
      PRINT * , '--> Selected Rosenbrock method not implemented'
    CASE (-3)    
      PRINT * , '--> Hmin/Hmax/Hstart must be positive'
    CASE (-4)    
      PRINT * , '--> FacMin/FacMax/FacRej must be positive'
    CASE (-5) 
      PRINT * , '--> Improper tolerance values'
    CASE (-6) 
      PRINT * , '--> No of steps exceeds maximum bound'
    CASE (-7) 
      PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
    CASE (-8)    
      PRINT * , '--> Matrix is repeatedly singular'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT
   
   PRINT *, "T=", T, "and H=", H
     
 END SUBROUTINE ros_ErrorMsg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_Integrator (Y, Tstart, Tend, T,  &
        AbsTol, RelTol,                          &
!~~~> Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic Rosenbrock method
!      defined by ros_S (no of stages)
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

!~~~> Input: the initial condition at Tstart; Output: the solution at T
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
!~~~> Input: integration interval
   REAL(kind=dp), INTENT(IN) :: Tstart,Tend
!~~~> Output: time at which the solution is returned (T=Tend if success)
   REAL(kind=dp), INTENT(OUT) ::  T
!~~~> Input: tolerances
   REAL(kind=dp), INTENT(IN) ::  AbsTol(N), RelTol(N)
!~~~> Input: integration parameters
   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   REAL(kind=dp), INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   REAL(kind=dp), INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables
   REAL(kind=dp) :: Ynew(N), Fcn0(N), Fcn(N)
   REAL(kind=dp) :: K(N*ros_S), dFdT(N)
#ifdef FULL_ALGEBRA    
   REAL(kind=dp) :: Jac0(N,N), Ghimj(N,N)
#else
   REAL(kind=dp) :: Jac0(LS_LU_NONZERO), Ghimj(LS_LU_NONZERO)
#endif
   REAL(kind=dp) :: H, Hnew, HC, HG, Fac, Tau
   REAL(kind=dp) :: Err, Yerr(N)
   INTEGER :: Pivot(N), Direction, ioffset, j, istage
   LOGICAL :: RejectLastH, RejectMoreH, Singular
!~~~>  Local parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp
!~~~>  Locally called functions
!    REAL(kind=dp) WLAMCH
!    EXTERNAL WLAMCH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~>  Initial preparations
   T = Tstart
   RSTATUS(Nhexit) = ZERO
   H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   IF (ABS(H) <= 10.0_dp*Roundoff) H = DeltaMin

   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF
   H = Direction*H

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.

!~~~> Time loop begins below

TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) )

   IF ( ISTATUS(Nstp) > Max_no_steps ) THEN  ! Too many steps
      CALL ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
      CALL ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF

!~~~>  Limit H if necessary to avoid going beyond Tend
   H = MIN(H,ABS(Tend-T))

!~~~>   Compute the function at current time
   CALL FunTemplate_14(T,Y,Fcn0, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1

!~~~>  Compute the function derivative with respect to T
   IF (.NOT.Autonomous) THEN
      CALL ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
   END IF

!~~~>   Compute the Jacobian at current time
   CALL JacTemplate_14(T,Y,Jac0,LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)   
   ISTATUS(Njac) = ISTATUS(Njac) + 1

!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO

   CALL ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular)
   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       CALL ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S

      ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
       ioffset = N*(istage-1)

      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         !slim: CALL WCOPY(N,Fcn0,1,Fcn,1)
	 Fcn(1:N) = Fcn0(1:N)
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
         !slim: CALL WCOPY(N,Y,1,Ynew,1)
	 Ynew(1:N) = Y(1:N)
         DO j = 1, istage-1
           CALL WAXPY(N,ros_A((istage-1)*(istage-2)/2+j), &
            K(N*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL FunTemplate_14(Tau,Ynew,Fcn, LS_NSEL, LS_type)
         ISTATUS(Nfun) = ISTATUS(Nfun) + 1
       END IF ! if istage == 1 elseif ros_NewF(istage)
       !slim: CALL WCOPY(N,Fcn,1,K(ioffset+1),1)
       K(ioffset+1:ioffset+N) = Fcn(1:N)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL WAXPY(N,HC,K(N*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL WAXPY(N,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL ros_Solve(Ghimj, Pivot, K(ioffset+1))

   END DO Stage


!~~~>  Compute the new solution
   !slim: CALL WCOPY(N,Y,1,Ynew,1)
   Ynew(1:N) = Y(1:N)
   DO j=1,ros_S
         CALL WAXPY(N,ros_M(j),K(N*(j-1)+1),1,Ynew,1)
   END DO

!~~~>  Compute the error estimation
   !slim: CALL WSCAL(N,ZERO,Yerr,1)
   Yerr(1:N) = ZERO
   DO j=1,ros_S
        CALL WAXPY(N,ros_E(j),K(N*(j-1)+1),1,Yerr,1)
   END DO
   Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac

!~~~>  Check the error magnitude and adjust step size
   ISTATUS(Nstp) = ISTATUS(Nstp) + 1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
      ISTATUS(Nacc) = ISTATUS(Nacc) + 1
      !slim: CALL WCOPY(N,Ynew,1,Y,1)
      Y(1:N) = Ynew(1:N)
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H)
      END IF
      RSTATUS(Nhexit) = H
      RSTATUS(Nhnew)  = Hnew
      RSTATUS(Ntexit) = T
      RejectLastH = .FALSE.
      RejectMoreH = .FALSE.
      H = Hnew
      EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
   ELSE           !~~~> Reject step
      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (ISTATUS(Nacc) >= 1)  ISTATUS(Nrej) = ISTATUS(Nrej) + 1
   END IF ! Err <= 1

   END DO UntilAccepted

   
   END DO TimeLoop

!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful

  END SUBROUTINE ros_Integrator


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  REAL(kind=dp) FUNCTION ros_ErrorNorm ( Y, Ynew, Yerr, &
                               AbsTol, RelTol, VectorTol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

! Input arguments
   REAL(kind=dp), INTENT(IN) :: Y(N), Ynew(N), &
          Yerr(N), AbsTol(N), RelTol(N)
   LOGICAL, INTENT(IN) ::  VectorTol
! Local variables
   REAL(kind=dp) :: Err, Scale, Ymax
   INTEGER  :: i
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp

   Err = ZERO
   DO i=1,N
     Ymax = MAX(ABS(Y(i)),ABS(Ynew(i)))
     IF (VectorTol) THEN
       Scale = AbsTol(i)+RelTol(i)*Ymax
     ELSE
       Scale = AbsTol(1)+RelTol(1)*Ymax
     END IF
     Err = Err+(Yerr(i)/Scale)**2
   END DO
   Err  = SQRT(Err/N)

   ros_ErrorNorm = MAX(Err,1.0d-10)

  END FUNCTION ros_ErrorNorm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> The time partial derivative of the function by finite differences
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

!~~~> Input arguments
   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(N), Fcn0(N)
!~~~> Output arguments
   REAL(kind=dp), INTENT(OUT) :: dFdT(N)
!~~~> Local variables
   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL FunTemplate_14(T+Delta,Y,dFdT, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1
   CALL WAXPY(N,(-ONE),Fcn0,1,dFdT,1)
   CALL WSCAL(N,(ONE/Delta),dFdT,1)

  END SUBROUTINE ros_FunTimeDerivative


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_PrepareMatrix ( H, Direction, gam, &
             Jac0, Ghimj, Pivot, Singular )
! --- --- --- --- --- --- --- --- --- --- --- --- ---
!  Prepares the LHS matrix for stage calculations
!  1.  Construct Ghimj = 1/(H*ham) - Jac0
!      "(Gamma H) Inverse Minus Jacobian"
!  2.  Repeat LU decomposition of Ghimj until successful.
!       -half the step size if LU decomposition fails and retry
!       -exit after 5 consecutive fails
! --- --- --- --- --- --- --- --- --- --- --- --- ---
   IMPLICIT NONE

!~~~> Input arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) ::  Jac0(N,N)
#else
   REAL(kind=dp), INTENT(IN) ::  Jac0(LS_LU_NONZERO)
#endif   
   REAL(kind=dp), INTENT(IN) ::  gam
   INTEGER, INTENT(IN) ::  Direction
!~~~> Output arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(OUT) :: Ghimj(N,N)
#else
   REAL(kind=dp), INTENT(OUT) :: Ghimj(LS_LU_NONZERO)
#endif   
   LOGICAL, INTENT(OUT) ::  Singular
   INTEGER, INTENT(OUT) ::  Pivot(N)
!~~~> Inout arguments
   REAL(kind=dp), INTENT(INOUT) :: H   ! step size is decreased when LU fails
!~~~> Local variables
   INTEGER  :: i, ISING, Nconsecutive
   REAL(kind=dp) :: ghinv
   REAL(kind=dp), PARAMETER :: ONE  = 1.0_dp, HALF = 0.5_dp

   Nconsecutive = 0
   Singular = .TRUE.

   DO WHILE (Singular)

!~~~>    Construct Ghimj = 1/(H*gam) - Jac0
#ifdef FULL_ALGEBRA    
     !slim: CALL WCOPY(N*N,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(N*N,(-ONE),Ghimj,1)
     Ghimj = -Jac0
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(i,i) = Ghimj(i,i)+ghinv
     END DO
#else
     !slim: CALL WCOPY(LS_LU_NONZERO,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(LS_LU_NONZERO,(-ONE),Ghimj,1)
     Ghimj(1:LS_LU_NONZERO) = -Jac0(1:LS_LU_NONZERO)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(LS_LU_DIAG(i)) = Ghimj(LS_LU_DIAG(i))+ghinv
     END DO
#endif   
!~~~>    Compute LU decomposition
     CALL ros_Decomp( Ghimj, Pivot, ISING )
     IF (ISING == 0) THEN
!~~~>    If successful done
        Singular = .FALSE.
     ELSE ! ISING .ne. 0
!~~~>    If unsuccessful half the step size; if 5 consecutive fails then return
        ISTATUS(Nsng) = ISTATUS(Nsng) + 1
        Nconsecutive = Nconsecutive+1
        Singular = .TRUE.
        PRINT*,'Warning: LU Decomposition returned ISING = ',ISING
        IF (Nconsecutive <= 5) THEN ! Less than 5 consecutive failed decompositions
           H = H*HALF
        ELSE  ! More than 5 consecutive failed decompositions
           RETURN
        END IF  ! Nconsecutive
      END IF    ! ISING

   END DO ! WHILE Singular

  END SUBROUTINE ros_PrepareMatrix


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Decomp( A, Pivot, ISING )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the LU decomposition
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Inout variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(INOUT) :: A(N,N)
#else   
   REAL(kind=dp), INTENT(INOUT) :: A(LS_LU_NONZERO)
#endif
!~~~> Output variables
   INTEGER, INTENT(OUT) :: Pivot(N), ISING

#ifdef FULL_ALGEBRA    
   CALL  DGETRF( N, N, A, N, Pivot, ISING )
#else   
   CALL KppDecomp ( A, ISING,LS_NSEL,LS_LU_NONZERO,LS_LU_CROW,LS_LU_DIAG,LS_LU_ICOL)   
   Pivot(1) = 1
#endif
   ISTATUS(Ndec) = ISTATUS(Ndec) + 1

  END SUBROUTINE ros_Decomp


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Solve( A, Pivot, b)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the forward/backward substitution (using pre-computed LU decomposition)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Input variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) :: A(N,N)
   INTEGER :: ISING
#else   
   REAL(kind=dp), INTENT(IN) :: A(LS_LU_NONZERO)
#endif
   INTEGER, INTENT(IN) :: Pivot(N)
!~~~> InOut variables
   REAL(kind=dp), INTENT(INOUT) :: b(N)

#ifdef FULL_ALGEBRA    
   CALL  DGETRS( 'N', N , 1, A, N, Pivot, b, N, ISING )
   IF ( Info < 0 ) THEN
      PRINT*,"Error in DGETRS. ISING=",ISING
   END IF  
#else  
   	CALL KppSolve_14(LS_LU_NONZERO,LS_NSEL, A, b,LS_type )
#endif

   ISTATUS(Nsol) = ISTATUS(Nsol) + 1

  END SUBROUTINE ros_Solve



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   DOUBLE PRECISION g

    g = 1.0_dp + 1.0_dp/SQRT(2.0_dp)
    rosMethod = RS2
!~~~> Name of the method
    ros_Name = 'ROS-2'
!~~~> Number of stages
    ros_S = 2

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = (1.0_dp)/g
    ros_C(1) = (-2.0_dp)/g
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
!~~~> M_i = Coefficients for new step solution
    ros_M(1)= (3.0_dp)/(2.0_dp*g)
    ros_M(2)= (1.0_dp)/(2.0_dp*g)
! E_i = Coefficients for error estimator
    ros_E(1) = 1.0_dp/(2.0_dp*g)
    ros_E(2) = 1.0_dp/(2.0_dp*g)
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus one
    ros_ELO = 2.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.0_dp
    ros_Alpha(2) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = g
    ros_Gamma(2) =-g

 END SUBROUTINE Ros2


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   rosMethod = RS3
!~~~> Name of the method
   ros_Name = 'ROS-3'
!~~~> Number of stages
   ros_S = 3

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1)= 1.0_dp
   ros_A(2)= 1.0_dp
   ros_A(3)= 0.0_dp

   ros_C(1) = -0.10156171083877702091975600115545E+01_dp
   ros_C(2) =  0.40759956452537699824805835358067E+01_dp
   ros_C(3) =  0.92076794298330791242156818474003E+01_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) =  0.1E+01_dp
   ros_M(2) =  0.61697947043828245592553615689730E+01_dp
   ros_M(3) = -0.42772256543218573326238373806514_dp
! E_i = Coefficients for error estimator
   ros_E(1) =  0.5_dp
   ros_E(2) = -0.29079558716805469821718236208017E+01_dp
   ros_E(3) =  0.22354069897811569627360909276199_dp
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1)= 0.0_dp
   ros_Alpha(2)= 0.43586652150845899941601945119356_dp
   ros_Alpha(3)= 0.43586652150845899941601945119356_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1)= 0.43586652150845899941601945119356_dp
   ros_Gamma(2)= 0.24291996454816804366592249683314_dp
   ros_Gamma(3)= 0.21851380027664058511513169485832E+01_dp

  END SUBROUTINE Ros3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
!     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1990)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RS4
!~~~> Name of the method
   ros_Name = 'ROS-4'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.2000000000000000E+01_dp
   ros_A(2) = 0.1867943637803922E+01_dp
   ros_A(3) = 0.2344449711399156_dp
   ros_A(4) = ros_A(2)
   ros_A(5) = ros_A(3)
   ros_A(6) = 0.0_dp

   ros_C(1) =-0.7137615036412310E+01_dp
   ros_C(2) = 0.2580708087951457E+01_dp
   ros_C(3) = 0.6515950076447975_dp
   ros_C(4) =-0.2137148994382534E+01_dp
   ros_C(5) =-0.3214669691237626_dp
   ros_C(6) =-0.6949742501781779_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .TRUE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 0.2255570073418735E+01_dp
   ros_M(2) = 0.2870493262186792_dp
   ros_M(3) = 0.4353179431840180_dp
   ros_M(4) = 0.1093502252409163E+01_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) =-0.2815431932141155_dp
   ros_E(2) =-0.7276199124938920E-01_dp
   ros_E(3) =-0.1082196201495311_dp
   ros_E(4) =-0.1093502252409163E+01_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 4.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.1145640000000000E+01_dp
   ros_Alpha(3) = 0.6552168638155900_dp
   ros_Alpha(4) = ros_Alpha(3)
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5728200000000000_dp
   ros_Gamma(2) =-0.1769193891319233E+01_dp
   ros_Gamma(3) = 0.7592633437920482_dp
   ros_Gamma(4) =-0.1049021087100450_dp

  END SUBROUTINE Ros4

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- A STIFFLY-STABLE METHOD, 4 stages, order 3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RD3
!~~~> Name of the method
   ros_Name = 'RODAS-3'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.0_dp
   ros_A(2) = 2.0_dp
   ros_A(3) = 0.0_dp
   ros_A(4) = 2.0_dp
   ros_A(5) = 0.0_dp
   ros_A(6) = 1.0_dp

   ros_C(1) = 4.0_dp
   ros_C(2) = 1.0_dp
   ros_C(3) =-1.0_dp
   ros_C(4) = 1.0_dp
   ros_C(5) =-1.0_dp
   ros_C(6) =-(8.0_dp/3.0_dp)

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .FALSE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .TRUE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 2.0_dp
   ros_M(2) = 0.0_dp
   ros_M(3) = 1.0_dp
   ros_M(4) = 1.0_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) = 0.0_dp
   ros_E(2) = 0.0_dp
   ros_E(3) = 0.0_dp
   ros_E(4) = 1.0_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.0_dp
   ros_Alpha(3) = 1.0_dp
   ros_Alpha(4) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5_dp
   ros_Gamma(2) = 1.5_dp
   ros_Gamma(3) = 0.0_dp
   ros_Gamma(4) = 0.0_dp

  END SUBROUTINE Rodas3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1996)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RD4
!~~~> Name of the method
    ros_Name = 'RODAS-4'
!~~~> Number of stages
    ros_S = 6

!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.000_dp
    ros_Alpha(2) = 0.386_dp
    ros_Alpha(3) = 0.210_dp
    ros_Alpha(4) = 0.630_dp
    ros_Alpha(5) = 1.000_dp
    ros_Alpha(6) = 1.000_dp

!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = 0.2500000000000000_dp
    ros_Gamma(2) =-0.1043000000000000_dp
    ros_Gamma(3) = 0.1035000000000000_dp
    ros_Gamma(4) =-0.3620000000000023E-01_dp
    ros_Gamma(5) = 0.0_dp
    ros_Gamma(6) = 0.0_dp

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = 0.1544000000000000E+01_dp
    ros_A(2) = 0.9466785280815826_dp
    ros_A(3) = 0.2557011698983284_dp
    ros_A(4) = 0.3314825187068521E+01_dp
    ros_A(5) = 0.2896124015972201E+01_dp
    ros_A(6) = 0.9986419139977817_dp
    ros_A(7) = 0.1221224509226641E+01_dp
    ros_A(8) = 0.6019134481288629E+01_dp
    ros_A(9) = 0.1253708332932087E+02_dp
    ros_A(10) =-0.6878860361058950_dp
    ros_A(11) = ros_A(7)
    ros_A(12) = ros_A(8)
    ros_A(13) = ros_A(9)
    ros_A(14) = ros_A(10)
    ros_A(15) = 1.0_dp

    ros_C(1) =-0.5668800000000000E+01_dp
    ros_C(2) =-0.2430093356833875E+01_dp
    ros_C(3) =-0.2063599157091915_dp
    ros_C(4) =-0.1073529058151375_dp
    ros_C(5) =-0.9594562251023355E+01_dp
    ros_C(6) =-0.2047028614809616E+02_dp
    ros_C(7) = 0.7496443313967647E+01_dp
    ros_C(8) =-0.1024680431464352E+02_dp
    ros_C(9) =-0.3399990352819905E+02_dp
    ros_C(10) = 0.1170890893206160E+02_dp
    ros_C(11) = 0.8083246795921522E+01_dp
    ros_C(12) =-0.7981132988064893E+01_dp
    ros_C(13) =-0.3152159432874371E+02_dp
    ros_C(14) = 0.1631930543123136E+02_dp
    ros_C(15) =-0.6058818238834054E+01_dp

!~~~> M_i = Coefficients for new step solution
    ros_M(1) = ros_A(7)
    ros_M(2) = ros_A(8)
    ros_M(3) = ros_A(9)
    ros_M(4) = ros_A(10)
    ros_M(5) = 1.0_dp
    ros_M(6) = 1.0_dp

!~~~> E_i  = Coefficients for error estimator
    ros_E(1) = 0.0_dp
    ros_E(2) = 0.0_dp
    ros_E(3) = 0.0_dp
    ros_E(4) = 0.0_dp
    ros_E(5) = 0.0_dp
    ros_E(6) = 1.0_dp

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.
    ros_NewF(5) = .TRUE.
    ros_NewF(6) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 4.0_dp

  END SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! STIFFLY-STABLE W METHOD OF ORDER 3, WITH 4 STAGES
!
! J. RANG and L. ANGERMANN
! NEW ROSENBROCK W-METHODS OF ORDER 3
! FOR PARTIAL DIFFERENTIAL ALGEBRAIC
!        EQUATIONS OF INDEX 1                                             
! BIT Numerical Mathematics (2005) 45: 761-787
!  DOI: 10.1007/s10543-005-0035-y
! Table 4.1-4.2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RG3
!~~~> Name of the method
    ros_Name = 'RANG-3'
!~~~> Number of stages
    ros_S = 4

    ros_A(1) = 5.09052051067020d+00;
    ros_A(2) = 5.09052051067020d+00;
    ros_A(3) = 0.0d0;
    ros_A(4) = 4.97628111010787d+00;
    ros_A(5) = 2.77268164715849d-02;
    ros_A(6) = 2.29428036027904d-01;

    ros_C(1) = -1.16790812312283d+01;
    ros_C(2) = -1.64057326467367d+01;
    ros_C(3) = -2.77268164715850d-01;
    ros_C(4) = -8.38103960500476d+00;
    ros_C(5) = -8.48328409199343d-01;
    ros_C(6) =  2.87009860433106d-01;

    ros_M(1) =  5.22582761233094d+00;
    ros_M(2) = -5.56971148154165d-01;
    ros_M(3) =  3.57979469353645d-01;
    ros_M(4) =  1.72337398521064d+00;

    ros_E(1) = -5.16845212784040d+00;
    ros_E(2) = -1.26351942603842d+00;
    ros_E(3) = -1.11022302462516d-16;
    ros_E(4) =  2.22044604925031d-16;

    ros_Alpha(1) = 0.0d00;
    ros_Alpha(2) = 2.21878746765329d+00;
    ros_Alpha(3) = 2.21878746765329d+00;
    ros_Alpha(4) = 1.55392337535788d+00;

    ros_Gamma(1) =  4.35866521508459d-01;
    ros_Gamma(2) = -1.78292094614483d+00;
    ros_Gamma(3) = -2.46541900496934d+00;
    ros_Gamma(4) = -8.05529997906370d-01;

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 3.0_dp

  END SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END SUBROUTINE Rosenbrock_14


SUBROUTINE FunTemplate_14( T, Y, Ydot, LS_NSEL, LS_type )
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Function
 INTEGER,INTENT(IN)::LS_NSEL, LS_type
   REAL(kind=dp) :: T, Y(LS_NSEL)
   REAL(kind=dp) :: Ydot(LS_NSEL)
   REAL(kind=dp) :: Told
   Told = TIME
   TIME = T
   CALL Fun_14( Y, FIX, RCONST, Ydot, LS_NSEL)
   TIME = Told
END SUBROUTINE FunTemplate_14

SUBROUTINE JacTemplate_14( T, Y, Jcb, LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Jacobian
 USE gckpp_LinearAlgebra
 
 INTEGER,INTENT(IN)::LS_NSEL, LS_LU_NONZERO, LS_type
 INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO), LS_LU_ICOL(LS_LU_NONZERO)
    REAL(kind=dp) :: T, Y(LS_NSEL)
#ifdef FULL_ALGEBRA    
    REAL(kind=dp) :: JV(LS_LU_NONZERO), Jcb(LS_NSEL,LS_NSEL)
#else
    REAL(kind=dp) :: Jcb(LS_LU_NONZERO)
#endif   
    REAL(kind=dp) :: Told
#ifdef FULL_ALGEBRA    
    INTEGER :: i, j
#endif   

    Told = TIME
    TIME = T
#ifdef FULL_ALGEBRA        	
	CALL Jac_SP_14(Y, FIX, RCONST, JV,LS_NSEL,LS_LU_NONZERO)	
    DO j=1,LS_NSEL
      DO i=1,LS_NSEL
         Jcb(i,j) = 0.0_dp
      END DO
    END DO
    DO i=1,LS_LU_NONZERO
       Jcb(LS_LU_IROW(i),LS_LU_ICOL(i)) = JV(i)
    END DO
#else   
       CALL Jac_SP_14( Y, FIX, RCONST, Jcb,LS_NSEL,LS_LU_NONZERO)	 
#endif   
    TIME = Told
END SUBROUTINE JacTemplate_14



SUBROUTINE INTEGRATE_15( TIN, TOUT, LS_type,LS_NSEL, LS_NDEL,&
  Prate, Lrate, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )

   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: TIN  ! Start Time
   REAL(kind=dp), INTENT(IN) :: TOUT ! End Time
   
   ! Optional input parameters and statistics
   INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   REAL(kind=dp), INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   REAL(kind=dp), INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U
   INTEGER,		  INTENT(IN) :: LS_type,LS_NSEL,LS_NDEL
   REAL(kind=dp), INTENT(IN) :: Prate(NVAR),Lrate(NVAR)
   
   REAL(kind=dp) :: RCNTRL(20), RSTATUS(20),deltaT
   INTEGER       :: ICNTRL(20), ISTATUS(20), IERR
   REAL(kind=dp) :: VAR_selected(LS_NSEL),VAR_deleted(LS_NDEL),LS_P(LS_NDEL),LS_L(LS_NDEL)
   
   INTEGER, SAVE :: Ntotal = 0

   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.0_dp
   ISTATUS(:) = 0
   RSTATUS(:) = 0.0_dp
   deltaT =  TOUT-TIN
   !~~~> fine-tune the integrator:
   ICNTRL(1) = 0	! 0 - non-autonomous, 1 - autonomous
   ICNTRL(2) = 0	! 0 - vector tolerances, 1 - scalars

   ! If optional parameters are given, and if they are >0, 
   ! then they overwrite default settings. 
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF
   
        VAR_selected=VAR(select_ind_15)
	    VAR_deleted=VAR(delete_ind_15)
		LS_P=Prate(delete_ind_15)
		LS_L=Lrate(delete_ind_15)
        CALL Rosenbrock_15(LU_NSEL_15,VAR_selected,TIN,TOUT,ATOL,RTOL,&
          RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		  LU_NONZERO_15,LU_NSEL_15,LU_CROW_15,LU_DIAG_15,LU_IROW_15,LU_ICOL_15, LS_type)		
		WHERE(LS_L<=(0.01/deltaT))
			VAR_deleted=VAR_deleted+deltaT*(LS_P-LS_L*VAR_deleted)
		ELSEWHERE
			VAR_deleted=LS_P/LS_L+(VAR_deleted-LS_P/LS_L)*EXP(-LS_L*deltaT)
		END WHERE

	    VAR(select_ind_15)=VAR_selected
	    VAR(delete_ind_15)=VAR_deleted
   
   !~~~> Debug option: show no of steps
   ! Ntotal = Ntotal + ISTATUS(Nstp)
   ! PRINT*,'NSTEPS=',ISTATUS(Nstp),' (',Ntotal,')','  O3=', VAR(ind_O3)

   STEPMIN = RSTATUS(Nhexit)
   ! if optional parameters are given for output they 
   ! are updated with the return information
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE INTEGRATE_15

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE Rosenbrock_15(N,Y,Tstart,Tend, &
           AbsTol,RelTol,              &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR, &
		   LS_LU_NONZERO,LS_NSEL,LS_LU_CROW,LS_LU_DIAG,LS_LU_IROW,LS_LU_ICOL, LS_type)
  USE gckpp_Parameters
  USE gckpp_LinearAlgebra
  IMPLICIT NONE

   INTEGER,       INTENT(IN)    :: N,LS_LU_NONZERO,LS_NSEL,LS_type   
   INTEGER,INTENT(IN)::LS_LU_CROW(LS_NSEL+1),LS_LU_DIAG(LS_NSEL+1)
   INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO),LS_LU_ICOL(LS_LU_NONZERO)   
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
   REAL(kind=dp), INTENT(IN)    :: Tstart,Tend
   REAL(kind=dp), INTENT(IN)    :: AbsTol(N),RelTol(N)
   INTEGER,       INTENT(IN)    :: ICNTRL(20)
   REAL(kind=dp), INTENT(IN)    :: RCNTRL(20)
   INTEGER,       INTENT(INOUT) :: ISTATUS(20)
   REAL(kind=dp), INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT)   :: IERR
!~~~>  Parameters of the Rosenbrock method, up to 6 stages
   INTEGER ::  ros_S, rosMethod
   INTEGER, PARAMETER :: RS2=1, RS3=2, RS4=3, RD3=4, RD4=5, RG3=6
   REAL(kind=dp) :: ros_A(15), ros_C(15), ros_M(6), ros_E(6), &
                    ros_Alpha(6), ros_Gamma(6), ros_ELO
   LOGICAL :: ros_NewF(6)
   CHARACTER(LEN=12) :: ros_Name
!~~~>  Local variables
   REAL(kind=dp) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   REAL(kind=dp) :: Hmin, Hmax, Hstart
   REAL(kind=dp) :: Texit
   INTEGER       :: i, UplimTol, Max_no_steps
   LOGICAL       :: Autonomous, VectorTol
!~~~>   Parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp

!~~~>  Initialize statistics
   ISTATUS(1:8) = 0
   RSTATUS(1:3) = ZERO

!~~~>  Autonomous or time dependent ODE. Default is time dependent.
   Autonomous = .NOT.(ICNTRL(1) == 0)

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:N) and RelTol(1:N)
   IF (ICNTRL(2) == 0) THEN
      VectorTol = .TRUE.
      UplimTol  = N
   ELSE
      VectorTol = .FALSE.
      UplimTol  = 1
   END IF

!~~~>   Initialize the particular Rosenbrock method selected
   SELECT CASE (ICNTRL(3))
     CASE (1)
       CALL Ros2
     CASE (2)
       CALL Ros3
     CASE (3)
       CALL Ros4
     CASE (0,4)
       CALL Rodas3
     CASE (5)
       CALL Rodas4
     CASE (6)
       CALL Rang3
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(3)=',ICNTRL(3) 
       CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT

!~~~>   The maximum number of steps admitted
   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 200000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Unit roundoff (1+Roundoff>1)
   Roundoff = WLAMCH('E')

!~~~>  Lower bound on the step size: (positive value)
   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Upper bound on the step size: (positive value)
   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Starting step size: (positive value)
   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hold < FacMax
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Check if tolerances are reasonable
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO


!~~~>  CALL Rosenbrock method
   CALL ros_Integrator(Y, Tstart, Tend, Texit,   &
        AbsTol, RelTol,                          &
!  Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!  Error indicator
        IERR)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS !  SUBROUTINES internal to Rosenbrock
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 SUBROUTINE ros_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
   REAL(kind=dp), INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR
   
   IERR = Code
   PRINT * , &
     'Forced exit from Rosenbrock due to the following error:' 
     
   SELECT CASE (Code)
    CASE (-1)
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)    
      PRINT * , '--> Selected Rosenbrock method not implemented'
    CASE (-3)    
      PRINT * , '--> Hmin/Hmax/Hstart must be positive'
    CASE (-4)    
      PRINT * , '--> FacMin/FacMax/FacRej must be positive'
    CASE (-5) 
      PRINT * , '--> Improper tolerance values'
    CASE (-6) 
      PRINT * , '--> No of steps exceeds maximum bound'
    CASE (-7) 
      PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
    CASE (-8)    
      PRINT * , '--> Matrix is repeatedly singular'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT
   
   PRINT *, "T=", T, "and H=", H
     
 END SUBROUTINE ros_ErrorMsg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_Integrator (Y, Tstart, Tend, T,  &
        AbsTol, RelTol,                          &
!~~~> Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic Rosenbrock method
!      defined by ros_S (no of stages)
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

!~~~> Input: the initial condition at Tstart; Output: the solution at T
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
!~~~> Input: integration interval
   REAL(kind=dp), INTENT(IN) :: Tstart,Tend
!~~~> Output: time at which the solution is returned (T=Tend if success)
   REAL(kind=dp), INTENT(OUT) ::  T
!~~~> Input: tolerances
   REAL(kind=dp), INTENT(IN) ::  AbsTol(N), RelTol(N)
!~~~> Input: integration parameters
   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   REAL(kind=dp), INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   REAL(kind=dp), INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables
   REAL(kind=dp) :: Ynew(N), Fcn0(N), Fcn(N)
   REAL(kind=dp) :: K(N*ros_S), dFdT(N)
#ifdef FULL_ALGEBRA    
   REAL(kind=dp) :: Jac0(N,N), Ghimj(N,N)
#else
   REAL(kind=dp) :: Jac0(LS_LU_NONZERO), Ghimj(LS_LU_NONZERO)
#endif
   REAL(kind=dp) :: H, Hnew, HC, HG, Fac, Tau
   REAL(kind=dp) :: Err, Yerr(N)
   INTEGER :: Pivot(N), Direction, ioffset, j, istage
   LOGICAL :: RejectLastH, RejectMoreH, Singular
!~~~>  Local parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp
!~~~>  Locally called functions
!    REAL(kind=dp) WLAMCH
!    EXTERNAL WLAMCH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~>  Initial preparations
   T = Tstart
   RSTATUS(Nhexit) = ZERO
   H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   IF (ABS(H) <= 10.0_dp*Roundoff) H = DeltaMin

   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF
   H = Direction*H

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.

!~~~> Time loop begins below

TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) )

   IF ( ISTATUS(Nstp) > Max_no_steps ) THEN  ! Too many steps
      CALL ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
      CALL ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF

!~~~>  Limit H if necessary to avoid going beyond Tend
   H = MIN(H,ABS(Tend-T))

!~~~>   Compute the function at current time
   CALL FunTemplate_15(T,Y,Fcn0, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1

!~~~>  Compute the function derivative with respect to T
   IF (.NOT.Autonomous) THEN
      CALL ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
   END IF

!~~~>   Compute the Jacobian at current time
   CALL JacTemplate_15(T,Y,Jac0,LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)   
   ISTATUS(Njac) = ISTATUS(Njac) + 1

!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO

   CALL ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular)
   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       CALL ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S

      ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
       ioffset = N*(istage-1)

      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         !slim: CALL WCOPY(N,Fcn0,1,Fcn,1)
	 Fcn(1:N) = Fcn0(1:N)
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
         !slim: CALL WCOPY(N,Y,1,Ynew,1)
	 Ynew(1:N) = Y(1:N)
         DO j = 1, istage-1
           CALL WAXPY(N,ros_A((istage-1)*(istage-2)/2+j), &
            K(N*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL FunTemplate_15(Tau,Ynew,Fcn, LS_NSEL, LS_type)
         ISTATUS(Nfun) = ISTATUS(Nfun) + 1
       END IF ! if istage == 1 elseif ros_NewF(istage)
       !slim: CALL WCOPY(N,Fcn,1,K(ioffset+1),1)
       K(ioffset+1:ioffset+N) = Fcn(1:N)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL WAXPY(N,HC,K(N*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL WAXPY(N,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL ros_Solve(Ghimj, Pivot, K(ioffset+1))

   END DO Stage


!~~~>  Compute the new solution
   !slim: CALL WCOPY(N,Y,1,Ynew,1)
   Ynew(1:N) = Y(1:N)
   DO j=1,ros_S
         CALL WAXPY(N,ros_M(j),K(N*(j-1)+1),1,Ynew,1)
   END DO

!~~~>  Compute the error estimation
   !slim: CALL WSCAL(N,ZERO,Yerr,1)
   Yerr(1:N) = ZERO
   DO j=1,ros_S
        CALL WAXPY(N,ros_E(j),K(N*(j-1)+1),1,Yerr,1)
   END DO
   Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac

!~~~>  Check the error magnitude and adjust step size
   ISTATUS(Nstp) = ISTATUS(Nstp) + 1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
      ISTATUS(Nacc) = ISTATUS(Nacc) + 1
      !slim: CALL WCOPY(N,Ynew,1,Y,1)
      Y(1:N) = Ynew(1:N)
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H)
      END IF
      RSTATUS(Nhexit) = H
      RSTATUS(Nhnew)  = Hnew
      RSTATUS(Ntexit) = T
      RejectLastH = .FALSE.
      RejectMoreH = .FALSE.
      H = Hnew
      EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
   ELSE           !~~~> Reject step
      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (ISTATUS(Nacc) >= 1)  ISTATUS(Nrej) = ISTATUS(Nrej) + 1
   END IF ! Err <= 1

   END DO UntilAccepted

   
   END DO TimeLoop

!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful

  END SUBROUTINE ros_Integrator


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  REAL(kind=dp) FUNCTION ros_ErrorNorm ( Y, Ynew, Yerr, &
                               AbsTol, RelTol, VectorTol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

! Input arguments
   REAL(kind=dp), INTENT(IN) :: Y(N), Ynew(N), &
          Yerr(N), AbsTol(N), RelTol(N)
   LOGICAL, INTENT(IN) ::  VectorTol
! Local variables
   REAL(kind=dp) :: Err, Scale, Ymax
   INTEGER  :: i
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp

   Err = ZERO
   DO i=1,N
     Ymax = MAX(ABS(Y(i)),ABS(Ynew(i)))
     IF (VectorTol) THEN
       Scale = AbsTol(i)+RelTol(i)*Ymax
     ELSE
       Scale = AbsTol(1)+RelTol(1)*Ymax
     END IF
     Err = Err+(Yerr(i)/Scale)**2
   END DO
   Err  = SQRT(Err/N)

   ros_ErrorNorm = MAX(Err,1.0d-10)

  END FUNCTION ros_ErrorNorm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> The time partial derivative of the function by finite differences
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

!~~~> Input arguments
   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(N), Fcn0(N)
!~~~> Output arguments
   REAL(kind=dp), INTENT(OUT) :: dFdT(N)
!~~~> Local variables
   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL FunTemplate_15(T+Delta,Y,dFdT, LS_NSEL, LS_type)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1
   CALL WAXPY(N,(-ONE),Fcn0,1,dFdT,1)
   CALL WSCAL(N,(ONE/Delta),dFdT,1)

  END SUBROUTINE ros_FunTimeDerivative


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_PrepareMatrix ( H, Direction, gam, &
             Jac0, Ghimj, Pivot, Singular )
! --- --- --- --- --- --- --- --- --- --- --- --- ---
!  Prepares the LHS matrix for stage calculations
!  1.  Construct Ghimj = 1/(H*ham) - Jac0
!      "(Gamma H) Inverse Minus Jacobian"
!  2.  Repeat LU decomposition of Ghimj until successful.
!       -half the step size if LU decomposition fails and retry
!       -exit after 5 consecutive fails
! --- --- --- --- --- --- --- --- --- --- --- --- ---
   IMPLICIT NONE

!~~~> Input arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) ::  Jac0(N,N)
#else
   REAL(kind=dp), INTENT(IN) ::  Jac0(LS_LU_NONZERO)
#endif   
   REAL(kind=dp), INTENT(IN) ::  gam
   INTEGER, INTENT(IN) ::  Direction
!~~~> Output arguments
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(OUT) :: Ghimj(N,N)
#else
   REAL(kind=dp), INTENT(OUT) :: Ghimj(LS_LU_NONZERO)
#endif   
   LOGICAL, INTENT(OUT) ::  Singular
   INTEGER, INTENT(OUT) ::  Pivot(N)
!~~~> Inout arguments
   REAL(kind=dp), INTENT(INOUT) :: H   ! step size is decreased when LU fails
!~~~> Local variables
   INTEGER  :: i, ISING, Nconsecutive
   REAL(kind=dp) :: ghinv
   REAL(kind=dp), PARAMETER :: ONE  = 1.0_dp, HALF = 0.5_dp

   Nconsecutive = 0
   Singular = .TRUE.

   DO WHILE (Singular)

!~~~>    Construct Ghimj = 1/(H*gam) - Jac0
#ifdef FULL_ALGEBRA    
     !slim: CALL WCOPY(N*N,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(N*N,(-ONE),Ghimj,1)
     Ghimj = -Jac0
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(i,i) = Ghimj(i,i)+ghinv
     END DO
#else
     !slim: CALL WCOPY(LS_LU_NONZERO,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(LS_LU_NONZERO,(-ONE),Ghimj,1)
     Ghimj(1:LS_LU_NONZERO) = -Jac0(1:LS_LU_NONZERO)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(LS_LU_DIAG(i)) = Ghimj(LS_LU_DIAG(i))+ghinv
     END DO
#endif   
!~~~>    Compute LU decomposition
     CALL ros_Decomp( Ghimj, Pivot, ISING )
     IF (ISING == 0) THEN
!~~~>    If successful done
        Singular = .FALSE.
     ELSE ! ISING .ne. 0
!~~~>    If unsuccessful half the step size; if 5 consecutive fails then return
        ISTATUS(Nsng) = ISTATUS(Nsng) + 1
        Nconsecutive = Nconsecutive+1
        Singular = .TRUE.
        PRINT*,'Warning: LU Decomposition returned ISING = ',ISING
        IF (Nconsecutive <= 5) THEN ! Less than 5 consecutive failed decompositions
           H = H*HALF
        ELSE  ! More than 5 consecutive failed decompositions
           RETURN
        END IF  ! Nconsecutive
      END IF    ! ISING

   END DO ! WHILE Singular

  END SUBROUTINE ros_PrepareMatrix


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Decomp( A, Pivot, ISING )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the LU decomposition
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Inout variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(INOUT) :: A(N,N)
#else   
   REAL(kind=dp), INTENT(INOUT) :: A(LS_LU_NONZERO)
#endif
!~~~> Output variables
   INTEGER, INTENT(OUT) :: Pivot(N), ISING

#ifdef FULL_ALGEBRA    
   CALL  DGETRF( N, N, A, N, Pivot, ISING )
#else   
   CALL KppDecomp ( A, ISING,LS_NSEL,LS_LU_NONZERO,LS_LU_CROW,LS_LU_DIAG,LS_LU_ICOL)   
   Pivot(1) = 1
#endif
   ISTATUS(Ndec) = ISTATUS(Ndec) + 1

  END SUBROUTINE ros_Decomp


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Solve( A, Pivot, b)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the forward/backward substitution (using pre-computed LU decomposition)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Input variables
#ifdef FULL_ALGEBRA    
   REAL(kind=dp), INTENT(IN) :: A(N,N)
   INTEGER :: ISING
#else   
   REAL(kind=dp), INTENT(IN) :: A(LS_LU_NONZERO)
#endif
   INTEGER, INTENT(IN) :: Pivot(N)
!~~~> InOut variables
   REAL(kind=dp), INTENT(INOUT) :: b(N)

#ifdef FULL_ALGEBRA    
   CALL  DGETRS( 'N', N , 1, A, N, Pivot, b, N, ISING )
   IF ( Info < 0 ) THEN
      PRINT*,"Error in DGETRS. ISING=",ISING
   END IF  
#else  
   	CALL KppSolve_15(LS_LU_NONZERO,LS_NSEL, A, b,LS_type )
#endif

   ISTATUS(Nsol) = ISTATUS(Nsol) + 1

  END SUBROUTINE ros_Solve



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   DOUBLE PRECISION g

    g = 1.0_dp + 1.0_dp/SQRT(2.0_dp)
    rosMethod = RS2
!~~~> Name of the method
    ros_Name = 'ROS-2'
!~~~> Number of stages
    ros_S = 2

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = (1.0_dp)/g
    ros_C(1) = (-2.0_dp)/g
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
!~~~> M_i = Coefficients for new step solution
    ros_M(1)= (3.0_dp)/(2.0_dp*g)
    ros_M(2)= (1.0_dp)/(2.0_dp*g)
! E_i = Coefficients for error estimator
    ros_E(1) = 1.0_dp/(2.0_dp*g)
    ros_E(2) = 1.0_dp/(2.0_dp*g)
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus one
    ros_ELO = 2.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.0_dp
    ros_Alpha(2) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = g
    ros_Gamma(2) =-g

 END SUBROUTINE Ros2


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   rosMethod = RS3
!~~~> Name of the method
   ros_Name = 'ROS-3'
!~~~> Number of stages
   ros_S = 3

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1)= 1.0_dp
   ros_A(2)= 1.0_dp
   ros_A(3)= 0.0_dp

   ros_C(1) = -0.10156171083877702091975600115545E+01_dp
   ros_C(2) =  0.40759956452537699824805835358067E+01_dp
   ros_C(3) =  0.92076794298330791242156818474003E+01_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) =  0.1E+01_dp
   ros_M(2) =  0.61697947043828245592553615689730E+01_dp
   ros_M(3) = -0.42772256543218573326238373806514_dp
! E_i = Coefficients for error estimator
   ros_E(1) =  0.5_dp
   ros_E(2) = -0.29079558716805469821718236208017E+01_dp
   ros_E(3) =  0.22354069897811569627360909276199_dp
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1)= 0.0_dp
   ros_Alpha(2)= 0.43586652150845899941601945119356_dp
   ros_Alpha(3)= 0.43586652150845899941601945119356_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1)= 0.43586652150845899941601945119356_dp
   ros_Gamma(2)= 0.24291996454816804366592249683314_dp
   ros_Gamma(3)= 0.21851380027664058511513169485832E+01_dp

  END SUBROUTINE Ros3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
!     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1990)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RS4
!~~~> Name of the method
   ros_Name = 'ROS-4'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.2000000000000000E+01_dp
   ros_A(2) = 0.1867943637803922E+01_dp
   ros_A(3) = 0.2344449711399156_dp
   ros_A(4) = ros_A(2)
   ros_A(5) = ros_A(3)
   ros_A(6) = 0.0_dp

   ros_C(1) =-0.7137615036412310E+01_dp
   ros_C(2) = 0.2580708087951457E+01_dp
   ros_C(3) = 0.6515950076447975_dp
   ros_C(4) =-0.2137148994382534E+01_dp
   ros_C(5) =-0.3214669691237626_dp
   ros_C(6) =-0.6949742501781779_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .TRUE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 0.2255570073418735E+01_dp
   ros_M(2) = 0.2870493262186792_dp
   ros_M(3) = 0.4353179431840180_dp
   ros_M(4) = 0.1093502252409163E+01_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) =-0.2815431932141155_dp
   ros_E(2) =-0.7276199124938920E-01_dp
   ros_E(3) =-0.1082196201495311_dp
   ros_E(4) =-0.1093502252409163E+01_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 4.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.1145640000000000E+01_dp
   ros_Alpha(3) = 0.6552168638155900_dp
   ros_Alpha(4) = ros_Alpha(3)
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5728200000000000_dp
   ros_Gamma(2) =-0.1769193891319233E+01_dp
   ros_Gamma(3) = 0.7592633437920482_dp
   ros_Gamma(4) =-0.1049021087100450_dp

  END SUBROUTINE Ros4

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- A STIFFLY-STABLE METHOD, 4 stages, order 3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RD3
!~~~> Name of the method
   ros_Name = 'RODAS-3'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.0_dp
   ros_A(2) = 2.0_dp
   ros_A(3) = 0.0_dp
   ros_A(4) = 2.0_dp
   ros_A(5) = 0.0_dp
   ros_A(6) = 1.0_dp

   ros_C(1) = 4.0_dp
   ros_C(2) = 1.0_dp
   ros_C(3) =-1.0_dp
   ros_C(4) = 1.0_dp
   ros_C(5) =-1.0_dp
   ros_C(6) =-(8.0_dp/3.0_dp)

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .FALSE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .TRUE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 2.0_dp
   ros_M(2) = 0.0_dp
   ros_M(3) = 1.0_dp
   ros_M(4) = 1.0_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) = 0.0_dp
   ros_E(2) = 0.0_dp
   ros_E(3) = 0.0_dp
   ros_E(4) = 1.0_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.0_dp
   ros_Alpha(3) = 1.0_dp
   ros_Alpha(4) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5_dp
   ros_Gamma(2) = 1.5_dp
   ros_Gamma(3) = 0.0_dp
   ros_Gamma(4) = 0.0_dp

  END SUBROUTINE Rodas3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1996)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RD4
!~~~> Name of the method
    ros_Name = 'RODAS-4'
!~~~> Number of stages
    ros_S = 6

!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.000_dp
    ros_Alpha(2) = 0.386_dp
    ros_Alpha(3) = 0.210_dp
    ros_Alpha(4) = 0.630_dp
    ros_Alpha(5) = 1.000_dp
    ros_Alpha(6) = 1.000_dp

!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = 0.2500000000000000_dp
    ros_Gamma(2) =-0.1043000000000000_dp
    ros_Gamma(3) = 0.1035000000000000_dp
    ros_Gamma(4) =-0.3620000000000023E-01_dp
    ros_Gamma(5) = 0.0_dp
    ros_Gamma(6) = 0.0_dp

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = 0.1544000000000000E+01_dp
    ros_A(2) = 0.9466785280815826_dp
    ros_A(3) = 0.2557011698983284_dp
    ros_A(4) = 0.3314825187068521E+01_dp
    ros_A(5) = 0.2896124015972201E+01_dp
    ros_A(6) = 0.9986419139977817_dp
    ros_A(7) = 0.1221224509226641E+01_dp
    ros_A(8) = 0.6019134481288629E+01_dp
    ros_A(9) = 0.1253708332932087E+02_dp
    ros_A(10) =-0.6878860361058950_dp
    ros_A(11) = ros_A(7)
    ros_A(12) = ros_A(8)
    ros_A(13) = ros_A(9)
    ros_A(14) = ros_A(10)
    ros_A(15) = 1.0_dp

    ros_C(1) =-0.5668800000000000E+01_dp
    ros_C(2) =-0.2430093356833875E+01_dp
    ros_C(3) =-0.2063599157091915_dp
    ros_C(4) =-0.1073529058151375_dp
    ros_C(5) =-0.9594562251023355E+01_dp
    ros_C(6) =-0.2047028614809616E+02_dp
    ros_C(7) = 0.7496443313967647E+01_dp
    ros_C(8) =-0.1024680431464352E+02_dp
    ros_C(9) =-0.3399990352819905E+02_dp
    ros_C(10) = 0.1170890893206160E+02_dp
    ros_C(11) = 0.8083246795921522E+01_dp
    ros_C(12) =-0.7981132988064893E+01_dp
    ros_C(13) =-0.3152159432874371E+02_dp
    ros_C(14) = 0.1631930543123136E+02_dp
    ros_C(15) =-0.6058818238834054E+01_dp

!~~~> M_i = Coefficients for new step solution
    ros_M(1) = ros_A(7)
    ros_M(2) = ros_A(8)
    ros_M(3) = ros_A(9)
    ros_M(4) = ros_A(10)
    ros_M(5) = 1.0_dp
    ros_M(6) = 1.0_dp

!~~~> E_i  = Coefficients for error estimator
    ros_E(1) = 0.0_dp
    ros_E(2) = 0.0_dp
    ros_E(3) = 0.0_dp
    ros_E(4) = 0.0_dp
    ros_E(5) = 0.0_dp
    ros_E(6) = 1.0_dp

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.
    ros_NewF(5) = .TRUE.
    ros_NewF(6) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 4.0_dp

  END SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! STIFFLY-STABLE W METHOD OF ORDER 3, WITH 4 STAGES
!
! J. RANG and L. ANGERMANN
! NEW ROSENBROCK W-METHODS OF ORDER 3
! FOR PARTIAL DIFFERENTIAL ALGEBRAIC
!        EQUATIONS OF INDEX 1                                             
! BIT Numerical Mathematics (2005) 45: 761-787
!  DOI: 10.1007/s10543-005-0035-y
! Table 4.1-4.2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RG3
!~~~> Name of the method
    ros_Name = 'RANG-3'
!~~~> Number of stages
    ros_S = 4

    ros_A(1) = 5.09052051067020d+00;
    ros_A(2) = 5.09052051067020d+00;
    ros_A(3) = 0.0d0;
    ros_A(4) = 4.97628111010787d+00;
    ros_A(5) = 2.77268164715849d-02;
    ros_A(6) = 2.29428036027904d-01;

    ros_C(1) = -1.16790812312283d+01;
    ros_C(2) = -1.64057326467367d+01;
    ros_C(3) = -2.77268164715850d-01;
    ros_C(4) = -8.38103960500476d+00;
    ros_C(5) = -8.48328409199343d-01;
    ros_C(6) =  2.87009860433106d-01;

    ros_M(1) =  5.22582761233094d+00;
    ros_M(2) = -5.56971148154165d-01;
    ros_M(3) =  3.57979469353645d-01;
    ros_M(4) =  1.72337398521064d+00;

    ros_E(1) = -5.16845212784040d+00;
    ros_E(2) = -1.26351942603842d+00;
    ros_E(3) = -1.11022302462516d-16;
    ros_E(4) =  2.22044604925031d-16;

    ros_Alpha(1) = 0.0d00;
    ros_Alpha(2) = 2.21878746765329d+00;
    ros_Alpha(3) = 2.21878746765329d+00;
    ros_Alpha(4) = 1.55392337535788d+00;

    ros_Gamma(1) =  4.35866521508459d-01;
    ros_Gamma(2) = -1.78292094614483d+00;
    ros_Gamma(3) = -2.46541900496934d+00;
    ros_Gamma(4) = -8.05529997906370d-01;

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 3.0_dp

  END SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END SUBROUTINE Rosenbrock_15


SUBROUTINE FunTemplate_15( T, Y, Ydot, LS_NSEL, LS_type )
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Function
 INTEGER,INTENT(IN)::LS_NSEL, LS_type
   REAL(kind=dp) :: T, Y(LS_NSEL)
   REAL(kind=dp) :: Ydot(LS_NSEL)
   REAL(kind=dp) :: Told
   Told = TIME
   TIME = T
   CALL Fun_15( Y, FIX, RCONST, Ydot, LS_NSEL)
   TIME = Told
END SUBROUTINE FunTemplate_15

SUBROUTINE JacTemplate_15( T, Y, Jcb, LS_NSEL, LS_LU_NONZERO, LS_LU_IROW, LS_LU_ICOL, LS_type)
 USE gckpp_Global, ONLY: FIX, RCONST, TIME
 USE gckpp_Jacobian
 USE gckpp_LinearAlgebra
 
 INTEGER,INTENT(IN)::LS_NSEL, LS_LU_NONZERO, LS_type
 INTEGER,INTENT(IN)::LS_LU_IROW(LS_LU_NONZERO), LS_LU_ICOL(LS_LU_NONZERO)
    REAL(kind=dp) :: T, Y(LS_NSEL)
#ifdef FULL_ALGEBRA    
    REAL(kind=dp) :: JV(LS_LU_NONZERO), Jcb(LS_NSEL,LS_NSEL)
#else
    REAL(kind=dp) :: Jcb(LS_LU_NONZERO)
#endif   
    REAL(kind=dp) :: Told
#ifdef FULL_ALGEBRA    
    INTEGER :: i, j
#endif   

    Told = TIME
    TIME = T
#ifdef FULL_ALGEBRA        	
	CALL Jac_SP_15(Y, FIX, RCONST, JV,LS_NSEL,LS_LU_NONZERO)	
    DO j=1,LS_NSEL
      DO i=1,LS_NSEL
         Jcb(i,j) = 0.0_dp
      END DO
    END DO
    DO i=1,LS_LU_NONZERO
       Jcb(LS_LU_IROW(i),LS_LU_ICOL(i)) = JV(i)
    END DO
#else   
       CALL Jac_SP_15( Y, FIX, RCONST, Jcb,LS_NSEL,LS_LU_NONZERO)	 
#endif   
    TIME = Told
END SUBROUTINE JacTemplate_15


END MODULE gckpp_Integrator