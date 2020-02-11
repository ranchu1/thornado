MODULE TwoMoment_PositivityLimiterModule_OrderV

  USE KindModule, ONLY: &
    DP, Zero, Half, One, &
    SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nNodesZ, &
    nDOFE, nDOFX
  USE TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_PL_Points, &
    Timer_PL_CellAverage
  USE ReferenceElementModule, ONLY: &
    nDOF_E, &
    nDOF_X1, &
    nDOF_X2, &
    nDOF_X3, &
    Weights_Q
  USE ReferenceElementModule_Lagrange, ONLY: &
    L_E_Dn,  L_E_Up, &
    L_X1_Dn, L_X1_Up, &
    L_X2_Dn, L_X2_Up, &
    L_X3_Dn, L_X3_Up
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply
  USE GeometryFieldsModuleE, ONLY: &
    nGE, iGE_Ep2
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    nCF
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3

  USE UtilitiesModule, ONLY: &
    WriteMatrix

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializePositivityLimiter_TwoMoment
  PUBLIC :: FinalizePositivityLimiter_TwoMoment
  PUBLIC :: ApplyPositivityLimiter_TwoMoment

  LOGICAL               :: UsePositivityLimiter
  LOGICAL               :: Verbose
  INTEGER               :: N_R
  INTEGER,    PARAMETER :: nPS = 9  ! Number of Positive Point Sets
  INTEGER               :: nPP(nPS) ! Number of Positive Points Per Set
  INTEGER               :: nPT      ! Total Number of Positive Points
  REAL(DP)              :: Min_1, Max_1, Min_2
  REAL(DP),   PARAMETER :: One_EPS = One - 1.0d1 * EPSILON( One )
  REAL(DP), ALLOCATABLE :: InterpMat(:,:)

CONTAINS


  SUBROUTINE InitializePositivityLimiter_TwoMoment &
    ( Min_1_Option, Max_1_Option, Min_2_Option, UsePositivityLimiter_Option, &
      Verbose_Option )

    REAL(DP), INTENT(in), OPTIONAL :: Min_1_Option
    REAL(DP), INTENT(in), OPTIONAL :: Max_1_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_2_Option
    LOGICAL,  INTENT(in), OPTIONAL :: UsePositivityLimiter_Option
    LOGICAL,  INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER :: i, iNodeZ, iOS

    PRINT*, "InitializePositivityLimiter_TwoMoment"

    IF( PRESENT( Min_1_Option ) )THEN
      Min_1 = Min_1_Option
    ELSE
      Min_1 = - HUGE( One )
    END IF

    IF( PRESENT( Max_1_Option ) )THEN
      Max_1 = Max_1_Option
    ELSE
      Max_1 = + HUGE( One )
    END IF

    IF( PRESENT( Min_2_Option ) )THEN
      Min_2 = Min_2_Option
    ELSE
      Min_2 = - HUGE( One )
    END IF

    IF( PRESENT( UsePositivityLimiter_Option ) )THEN
      UsePositivityLimiter = UsePositivityLimiter_Option
    ELSE
      UsePositivityLimiter = .TRUE.
    END IF

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .TRUE.
    END IF

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A2,A6,A)') '', 'INFO: ', 'InitializePositivityLimiter'
      WRITE(*,*)
      WRITE(*,'(A6,A,L1)') &
        '', 'Use Positivity Limiter: ', UsePositivityLimiter 
      WRITE(*,*)
      WRITE(*,'(A6,A12,ES23.15E3)') '', 'Min_1 = ', Min_1
      WRITE(*,'(A6,A12,ES23.15E3)') '', 'Max_1 = ', Max_1
      WRITE(*,'(A6,A12,ES23.15E3)') '', 'Min_2 = ', Min_2

    END IF

    nPP    = 0
    nPP(1) = PRODUCT( nNodesZ )

    DO i = 1, 4

      IF( nNodesZ(i) > 1 )THEN

        nPP(2*i:2*i+1) = PRODUCT( nNodesZ, MASK = [1,2,3,4] .NE. i )

      END IF

    END DO

    nPT = SUM( nPP )

    PRINT*, "nPP = ", nPP
    PRINT*, "nPT = ", nPT

    ALLOCATE( InterpMat(nPT,nDOFZ) )

    InterpMat = Zero
    DO iNodeZ = 1, nDOFZ

      InterpMat(iNodeZ,iNodeZ) = One

      IF( SUM( nPP(2:3) ) > 0 )THEN

        iOS = nPP(1)
        InterpMat(iOS+1:iOS+nDOF_E,iNodeZ) = L_E_Dn(1:nDOF_E,iNodeZ)

        iOS = iOS + nPP(2)
        InterpMat(iOS+1:iOS+nDOF_E,iNodeZ) = L_E_Up(1:nDOF_E,iNodeZ)

      END IF

      IF( SUM( nPP(4:5) ) > 0 )THEN

        iOS = SUM( nPP(1:3) )
        InterpMat(iOS+1:iOS+nDOF_X1,iNodeZ) = L_X1_Dn(1:nDOF_X1,iNodeZ)

        iOS = iOS + nPP(4)
        InterpMat(iOS+1:iOS+nDOF_X1,iNodeZ) = L_X1_Up(1:nDOF_X1,iNodeZ)

      END IF

      IF( SUM( nPP(6:7) ) > 0 )THEN

        iOS = SUM( nPP(1:5) )
        InterpMat(iOS+1:iOS+nDOF_X2,iNodeZ) = L_X2_Dn(1:nDOF_X2,iNodeZ)

        iOS = iOS + nPP(6)
        InterpMat(iOS+1:iOS+nDOF_X2,iNodeZ) = L_X2_Up(1:nDOF_X2,iNodeZ)

      END IF

      IF( SUM( nPP(8:9) ) > 0 )THEN

        iOS = SUM( nPP(1:7) )
        InterpMat(iOS+1:iOS+nDOF_X3,iNodeZ) = L_X3_Dn(1:nDOF_X3,iNodeZ)

        iOS = iOS + nPP(8)
        InterpMat(iOS+1:iOS+nDOF_X3,iNodeZ) = L_X3_Up(1:nDOF_X3,iNodeZ)

      END IF

    END DO

  END SUBROUTINE InitializePositivityLimiter_TwoMoment


  SUBROUTINE FinalizePositivityLimiter_TwoMoment

    PRINT*, "FinalizePositivityLimiter_TwoMoment"

    DEALLOCATE( InterpMat )

  END SUBROUTINE FinalizePositivityLimiter_TwoMoment


  SUBROUTINE ApplyPositivityLimiter_TwoMoment &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE (1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
                  iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      U_F(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
                  iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(inout) :: &
      U_R(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                  iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    LOGICAL  :: RecomputePointValues
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS, iP
    INTEGER  :: iNodeZ, iNodeE, iNodeX
    REAL(DP) :: Min_K, Max_K, Theta_1, Theta_2, Theta_P
    REAL(DP) :: Gamma, Gamma_Min
    REAL(DP) :: Tau_Q(nDOFZ, &
                      iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                      iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: N_Q(nDOFZ, &
                    iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                    iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                    nSpecies), &
                N_P(nPT  , &
                    iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                    iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                    nSpecies), &
                N_K(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                    iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                    nSpecies)
    REAL(DP) :: G1_Q(nDOFZ, &
                       iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                       iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                       nSpecies), &
                G1_P(nPT  , &
                       iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                       iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                       nSpecies), &
                G1_K(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                       iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                       nSpecies)
    REAL(DP) :: G2_Q(nDOFZ, &
                     iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                     nSpecies), &
                G2_P(nPT  , &
                     iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                     nSpecies), &
                G2_K(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                     nSpecies)
    REAL(DP) :: G3_Q(nDOFZ, &
                     iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                     nSpecies), &
                G3_P(nPT  , &
                     iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                     nSpecies), &
                G3_K(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                     nSpecies)

    IF( nDOFZ == 1 ) RETURN

    PRINT*, "      ApplyPositivityLimiter_TwoMoment"

    N_R = nSpecies * PRODUCT( iZ_E0 - iZ_B0 + 1 )

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        Tau_Q(iNodeZ,iZ1,iZ2,iZ3,iZ4) &
          = GE(iNodeE,iZ1,iGE_Ep2) * GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm)

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO

    DO iS = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      N_Q (:,iZ1,iZ2,iZ3,iZ4,iS) = U_R(:,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS)
      G1_Q(:,iZ1,iZ2,iZ3,iZ4,iS) = U_R(:,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS)
      G2_Q(:,iZ1,iZ2,iZ3,iZ4,iS) = U_R(:,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS)
      G3_Q(:,iZ1,iZ2,iZ3,iZ4,iS) = U_R(:,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStart( Timer_PL_Points )

    CALL ComputePointValues( iZ_B0, iZ_E0, N_Q , N_P  )
    CALL ComputePointValues( iZ_B0, iZ_E0, G1_Q, G1_P )
    CALL ComputePointValues( iZ_B0, iZ_E0, G2_Q, G2_P )
    CALL ComputePointValues( iZ_B0, iZ_E0, G3_Q, G3_P )

    CALL TimersStop( Timer_PL_Points )

    CALL TimersStart( Timer_PL_CellAverage )

    CALL ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, N_Q , N_K  )
    CALL ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, G1_Q, G1_K )
    CALL ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, G2_Q, G2_K )
    CALL ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, G3_Q, G3_K )

    CALL TimersStop( Timer_PL_CellAverage )

    ! --- Ensure Bounded Density ---

    RecomputePointValues = .FALSE.

    DO iS = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      Min_K = Min_1
      Max_K = Max_1

      DO iP = 1, nPT

        Min_K = MIN( Min_K, N_P(iP,iZ1,iZ2,iZ3,iZ4,iS) )
        Max_K = MAX( Max_K, N_P(iP,iZ1,iZ2,iZ3,iZ4,iS) )

      END DO

      IF( Min_K < Min_1 .OR. Max_K > Max_K )THEN

        Theta_1 &
          = MIN( One, &
                 ABS( ( Min_1 - N_K(iZ1,iZ2,iZ3,iZ4,iS) ) &
                      / ( Min_K - N_K(iZ1,iZ2,iZ3,iZ4,iS)+SqrtTiny ) ), &
                 ABS( ( Max_1 - N_K(iZ1,iZ2,iZ3,iZ4,iS) ) &
                      / ( Max_K - N_K(iZ1,iZ2,iZ3,iZ4,iS)+SqrtTiny ) ) )

        Theta_1 = One_EPS * Theta_1

        PRINT*, "Min_K, Max_K = ", Min_K, Max_K
        PRINT*, "Theta_1      = ", Theta_1

        N_Q(:,iZ1,iZ2,iZ3,iZ4,iS) &
          = Theta_1 * N_Q(:,iZ1,iZ2,iZ3,iZ4,iS) &
            + ( One - Theta_1 ) * N_K(iZ1,iZ2,iZ3,iZ4,iS)

        RecomputePointValues = .TRUE.

      END IF

    END DO
    END DO
    END DO
    END DO
    END DO

    IF( RecomputePointValues )THEN

      CALL ComputePointValues( iZ_B0, iZ_E0, N_Q , N_P  )

    END IF

    ! --- Ensure Positive "Gamma" ---

    DO iS = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      Gamma_Min = Min_2
      Theta_2   = One

      DO iP = 1, nPT

        Gamma &
          = GammaFun &
              ( N_P (iP,iZ1,iZ2,iZ3,iZ4,iS), &
                G1_P(iP,iZ1,iZ2,iZ3,iZ4,iS), &
                G2_P(iP,iZ1,iZ2,iZ3,iZ4,iS), &
                G3_P(iP,iZ1,iZ2,iZ3,iZ4,iS) )

        Gamma_Min = MIN( Gamma, Gamma_Min )

        IF( Gamma_Min < Min_2 )THEN

          CALL SolveTheta_Bisection &
                 ( N_P (iP,iZ1,iZ2,iZ3,iZ4,iS), &
                   G1_P(iP,iZ1,iZ2,iZ3,iZ4,iS), &
                   G2_P(iP,iZ1,iZ2,iZ3,iZ4,iS), &
                   G3_P(iP,iZ1,iZ2,iZ3,iZ4,iS), &
                   N_K (   iZ1,iZ2,iZ3,iZ4,iS), &
                   G1_K(   iZ1,iZ2,iZ3,iZ4,iS), &
                   G2_K(   iZ1,iZ2,iZ3,iZ4,iS), &
                   G3_K(   iZ1,iZ2,iZ3,iZ4,iS), &
                   Theta_P )

          Theta_2 = MIN( Theta_2, Theta_P )

        END IF

      END DO

      IF( Gamma_Min < Min_2 )THEN

        ! --- Limit Towards Cell Average ---

        PRINT*, "Gamma_Min = ", Gamma_Min
        PRINT*, "Theta_2   = ", Theta_2

        Theta_2 = One_EPS * Theta_2

        N_Q (:,iZ1,iZ2,iZ3,iZ4,iS) &
          = Theta_2 * N_Q (:,iZ1,iZ2,iZ3,iZ4,iS) &
            + ( One - Theta_2 ) * N_K (iZ1,iZ2,iZ3,iZ4,iS)

        G1_Q(:,iZ1,iZ2,iZ3,iZ4,iS) &
          = Theta_2 * G1_Q(:,iZ1,iZ2,iZ3,iZ4,iS) &
            + ( One - Theta_2 ) * G1_K(iZ1,iZ2,iZ3,iZ4,iS)

        G2_Q(:,iZ1,iZ2,iZ3,iZ4,iS) &
          = Theta_2 * G2_Q(:,iZ1,iZ2,iZ3,iZ4,iS) &
            + ( One - Theta_2 ) * G2_K(iZ1,iZ2,iZ3,iZ4,iS)

        G3_Q(:,iZ1,iZ2,iZ3,iZ4,iS) &
          = Theta_2 * G3_Q(:,iZ1,iZ2,iZ3,iZ4,iS) &
            + ( One - Theta_2 ) * G3_K(iZ1,iZ2,iZ3,iZ4,iS)

      END IF

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iS = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      U_R(:,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS) = N_Q (:,iZ1,iZ2,iZ3,iZ4,iS)
      U_R(:,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS) = G1_Q(:,iZ1,iZ2,iZ3,iZ4,iS)
      U_R(:,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS) = G2_Q(:,iZ1,iZ2,iZ3,iZ4,iS)
      U_R(:,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS) = G3_Q(:,iZ1,iZ2,iZ3,iZ4,iS)

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ApplyPositivityLimiter_TwoMoment


  SUBROUTINE ComputePointValues( iZ_B0, iZ_E0, U_Q, U_P )

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4)
    REAL(DP), INTENT(in)  :: &
      U_Q(nDOFZ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP), INTENT(out) :: &
      U_P(nPT  ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nPT, N_R, nDOFZ, One, InterpMat, nPT, &
             U_Q, nDOFZ, Zero, U_P, nPT )

  END SUBROUTINE ComputePointValues


  SUBROUTINE ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, U_Q, U_K )

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4)
    REAL(DP), INTENT(in)  :: &
      Tau_Q(nDOFZ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                  iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in)  :: &
      U_Q  (nDOFZ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                  iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP), INTENT(out) :: &
      U_K  (      iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                  iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iS

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      U_K(iZ1,iZ2,iZ3,iZ4,iS) &
        = SUM( Weights_Q(:) * Tau_Q(:,iZ1,iZ2,iZ3,iZ4) &
                 * U_Q(:,iZ1,iZ2,iZ3,iZ4,iS) ) &
          / SUM( Weights_Q(:) * Tau_Q(:,iZ1,iZ2,iZ3,iZ4) )

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ComputeCellAverage


  REAL(DP) FUNCTION GammaFun( N, G1, G2, G3 )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: N, G1, G2, G3

    GammaFun = N - SQRT( G1**2 + G2**2 + G3**2 )

    RETURN
  END FUNCTION GammaFun


  SUBROUTINE SolveTheta_Bisection &
    ( N_P, G1_P, G2_P, G3_P, N_K, G1_K, G2_K, G3_K, Theta )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: N_P, G1_P, G2_P, G3_P
    REAL(DP), INTENT(in)  :: N_K, G1_K, G2_K, G3_K
    REAL(DP), INTENT(out) :: Theta

    INTEGER,  PARAMETER :: ITERATION_MAX = 12
    REAL(DP), PARAMETER :: dx_min = 1.0d-3

    LOGICAL  :: CONVERGED
    INTEGER  :: ITERATION
    REAL(DP) :: x_a, x_b, x_c, dx
    REAL(DP) :: f_a, f_b, f_c

    x_a = Zero
    f_a = GammaFun( N_K, G1_K, G2_K, G3_K ) - Min_2

    x_b = One
    f_b = GammaFun( N_P, G1_P, G2_P, G3_P ) - Min_2

    dx = One

    ITERATION = 0
    CONVERGED = .FALSE.
    DO WHILE ( .NOT. CONVERGED .AND. ITERATION < ITERATION_MAX )

      ITERATION = ITERATION + 1

      dx = Half * dx
      x_c = x_a + dx

      f_c = GammaFun &
              ( x_c * N_P  + ( One - x_c ) * N_K,  &
                x_c * G1_P + ( One - x_c ) * G1_K, &
                x_c * G2_P + ( One - x_c ) * G2_K, &
                x_c * G3_P + ( One - x_c ) * G3_K ) - Min_2

      IF( f_a * f_c < Zero )THEN

        x_b = x_c
        f_b = f_c

      ELSE

        x_a = x_c
        f_a = f_c

      END IF

      IF( dx < dx_min ) CONVERGED = .TRUE.

    END DO

    IF( ITERATION >= ITERATION_MAX )THEN
      Theta = Zero
    ELSE
      Theta = x_a
    END IF

  END SUBROUTINE SolveTheta_Bisection


END MODULE TwoMoment_PositivityLimiterModule_OrderV