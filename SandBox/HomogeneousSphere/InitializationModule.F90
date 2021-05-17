MODULE InitializationModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two, TwoPi, SqrtTiny
  USE UnitsModule, ONLY: &
    Gram, Centimeter, &
    Kilometer, Kelvin, &
    MeV, BoltzmannConstant
  USE UtilitiesModule, ONLY: &
    Locate, &
    Interpolate1D_Linear, &
    NodeNumberX
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nX, nNodesX, &
    iX_B0, iX_E0, &
    iE_B0, iE_E0, &
    nDOF, nDOFX, nDOFE
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE ReferenceElementModule, ONLY: &
    NodeNumbersX, &
    NodeNumberTable
  USE MeshModule, ONLY: &
    MeshE, MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    CoordinateSystem, &
    uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uAF, iAF_P, iAF_T, iAF_Ye, iAF_S, iAF_E, &
    iAF_Me, iAF_Mp, iAF_Mn, iAF_Xp, iAF_Xn, &
    iAF_Xa, iAF_Xh, iAF_Gm
  USE RadiationFieldsModule, ONLY: &
    nSpecies, iNuE, iNuE_Bar, &
    uPR, nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3, &
    uCR, nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE EquationOfStateModule_TABLE, ONLY: &
    ApplyEquationOfState_TABLE, &
    ComputeThermodynamicStates_Primitive_TABLE
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeNeutrinoOpacities_EC_Points, &
    ComputeEquilibriumDistributions_Points
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputeConserved_TwoMoment

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  REAL(DP), PARAMETER :: OuterD = 1.0d8 * Gram / Centimeter**3
  REAL(DP), PARAMETER :: OuterT = 0.2_DP * MeV
  REAL(DP), PARAMETER :: OuterY = 0.4643_DP
  REAL(DP), PARAMETER :: R0 = 1.0d2 * Kilometer

  PUBLIC :: InitializeFields_HomogeneousSphere

CONTAINS

  SUBROUTINE InitializeFields_HomogeneousSphere( CentralConditions )

    CHARACTER(2), INTENT(in) :: CentralConditions

    CALL InitializeFluidFields_HomogeneousSphere( CentralConditions )

    CALL InitializeRadiationFields_HomogeneousSphere

  END SUBROUTINE InitializeFields_HomogeneousSphere


  SUBROUTINE InitializeFluidFields_HomogeneousSphere( CentralConditions_Option )

    CHARACTER(2), INTENT(in), OPTIONAL :: CentralConditions_Option

    CHARACTER(2) :: CentralConditions
    INTEGER      :: iX1, iX2, iX3, iE
    INTEGER      :: iNodeX1, iNodeX2, iNodeX3, iNodeE
    INTEGER      :: iNodeX, iNode
    REAL(DP)     :: X1, E
    REAL(DP)     :: Ratio, Width
    REAL(DP)     :: InnerD, InnerT, InnerY

    CentralConditions = '02'

    IF( PRESENT( CentralConditions_Option ) ) &
      CentralConditions = CentralConditions_Option

    Width  = 1.0d1 * Kilometer

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)

    ! --- Initialize Fluid Fields ---

    DO iX3 = 1, nX(3)
    DO iX2 = 1, nX(2)
    DO iX1 = 1, nX(1)

      DO iNodeX3 = 1, nNodesX(3)
        DO iNodeX2 = 1, nNodesX(2)
          DO iNodeX1 = 1, nNodesX(1)

            iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

            Ratio = Zero
            !Ratio = One / ( EXP( 1.0d1 / Width * (X1 - R0) ) + One )
            !Ratio = One / ( (X1/R0)**80 + One )
            IF( X1 <= R0 ) Ratio = One

            SELECT CASE ( CentralConditions )

              CASE ( '01' )

                InnerD = 1.0d14 * Gram / Centimeter**3
                InnerT = 21.0_DP * MeV
                InnerY = 0.25_DP

              CASE ( ' A' )

                InnerD = 1.0d13 * Gram / Centimeter**3
                InnerT = 16.0_DP * MeV
                InnerY = 0.14_DP

              CASE ( '03' )

                InnerD = 1.0d12 * Gram / Centimeter**3
                InnerT = 8.0_DP * MeV
                InnerY = 0.12_DP

              CASE ( ' C' )

                InnerD = 1.0d11 * Gram / Centimeter**3
                InnerT = 8.0_DP * MeV
                InnerY = 0.15_DP

              CASE ( '05' )

                InnerD = 1.0d10 * Gram / Centimeter**3
                InnerT = 3.0_DP * MeV
                InnerY = 0.26_DP

              CASE ( ' B' )

                InnerD = 3.0d11 * Gram / Centimeter**3
                InnerT = 4.0_DP * MeV
                InnerY = 0.20_DP

              CASE DEFAULT

                InnerD = 1.0d13 * Gram / Centimeter**3
                InnerT = 16.0_DP * MeV
                InnerY = 0.14_DP

            END SELECT

            uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
              = OuterD + ( InnerD - OuterD ) * Ratio

            uAF(iNodeX,iX1,iX2,iX3,iAF_T) &
              = OuterT + ( InnerT - OuterT ) * Ratio

            uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
              = OuterY + ( InnerY - OuterY ) * Ratio

          END DO
        END DO
      END DO

      CALL ComputeThermodynamicStates_Primitive_TABLE &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E ), &
               uAF(:,iX1,iX2,iX3,iAF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne) )

      uPF(:,iX1,iX2,iX3,iPF_V1) = Zero
      uPF(:,iX1,iX2,iX3,iPF_V2) = Zero
      uPF(:,iX1,iX2,iX3,iPF_V3) = Zero

      CALL ApplyEquationOfState_TABLE &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_P ), &
               uAF(:,iX1,iX2,iX3,iAF_S ), uAF(:,iX1,iX2,iX3,iAF_E ), &
               uAF(:,iX1,iX2,iX3,iAF_Me), uAF(:,iX1,iX2,iX3,iAF_Mp), &
               uAF(:,iX1,iX2,iX3,iAF_Mn), uAF(:,iX1,iX2,iX3,iAF_Xp), &
               uAF(:,iX1,iX2,iX3,iAF_Xn), uAF(:,iX1,iX2,iX3,iAF_Xa), &
               uAF(:,iX1,iX2,iX3,iAF_Xh), uAF(:,iX1,iX2,iX3,iAF_Gm) )

      ! --- Conserved ---

      uCF(:,iX1,iX2,iX3,iCF_D ) = uPF(:,iX1,iX2,iX3,iPF_D)
      uCF(:,iX1,iX2,iX3,iCF_S1) = Zero
      uCF(:,iX1,iX2,iX3,iCF_S2) = Zero
      uCF(:,iX1,iX2,iX3,iCF_S3) = Zero
      uCF(:,iX1,iX2,iX3,iCF_E ) = uPF(:,iX1,iX2,iX3,iPF_E)
      uCF(:,iX1,iX2,iX3,iCF_Ne) = uPF(:,iX1,iX2,iX3,iPF_Ne)

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFluidFields_HomogeneousSphere


  SUBROUTINE InitializeRadiationFields_HomogeneousSphere

    INTEGER  :: iE, iX1, iX2, iX3, iS
    INTEGER  :: iNode, iNodeE
    REAL(DP) :: Gm_dd_11(nDOF)
    REAL(DP) :: Gm_dd_22(nDOF)
    REAL(DP) :: Gm_dd_33(nDOF)

    INTEGER  :: i, nR, nE, iR
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    REAL(DP), ALLOCATABLE :: R_P(:), D_P(:), T_P(:), Y_P(:)
    REAL(DP), ALLOCATABLE :: Chi(:,:,:), fEQ(:,:,:), E_Nu(:)
    REAL(DP), ALLOCATABLE :: D_Nu_P(:,:,:), I1_Nu_P(:,:,:)

    REAL(DP) :: UnitR = Kilometer, &
                UnitD = Gram / Centimeter**3, &
                UnitT = Kelvin

    nR = ( iX_E0(1) - iX_B0(1) + 1 ) * nDOFX 
    nE = ( iE_E0 - iE_B0 + 1 ) * nDOFE

    ALLOCATE( E_Nu(nE) )
    ALLOCATE( R_P(nR) )
    ALLOCATE( D_P(nR) )
    ALLOCATE( T_P(nR) )
    ALLOCATE( Y_P(nR) )
    ALLOCATE( Chi(nE,nR,nSpecies) )
    ALLOCATE( fEQ(nE,nR,nSpecies) )
    ALLOCATE( D_Nu_P (nE,nR,nSpecies) )
    ALLOCATE( I1_Nu_P(nE,nR,nSpecies) )

    ! --- Compute Neutrino Energies ---
 
    i = 1
    DO iE = iE_B0, iE_E0
    DO iNodeE = 1, nDOFE

      E_Nu(i) = NodeCoordinate( MeshE, iE, iNodeE )

      i = i + 1

    END DO
    END DO

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
      DO iNodeX = 1, nDOFX

         iR = nDOFX * (iX1-1) + iNodeX
         iNodeX1 = NodeNumberTableX(1,iNodeX)

         R_P(iR) = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
         D_P(iR) = uPF(iNodeX,iX1,iX2,iX3,iPF_D )
         T_P(iR) = uAF(iNodeX,iX1,iX2,iX3,iAF_T )
         Y_P(iR) = uAF(iNodeX,iX1,iX2,iX3,iAF_Ye)

      END DO
    END DO
    END DO
    END DO

    ! --- Compute Neutrino Absorption Opacities ---

    DO iS = 1, nSpecies

      CALL ComputeNeutrinoOpacities_EC_Points &
             ( 1, nE, 1, nR, E_Nu, D_P, T_P, Y_P, iS, Chi(:,:,iS) )

      CALL ComputeEquilibriumDistributions_Points &
             ( 1, nE, 1, nR, E_Nu, D_P, T_P, Y_P, fEQ(:,:,iS), iS )

    END DO

    DO iS = 1, nSpecies
    DO iR = 1, nR
    DO iE = 1, nE

      ! --- Use Chi and fEQ to compute analytic solution ---
      IF( R_P(iR) <= R0 ) THEN

        CALL ComputeHomogeneousSphereSolution &
               ( Chi(iE,iR,iS), fEQ(iE,iR,iS), R_P(iR), &
                 D_Nu_P(iE,iR,iS), I1_Nu_P(iE,iR,iS) )

      ELSE

        CALL ComputeHomogeneousSphereSolution &
               ( Chi(iE,1,iS), fEQ(iE,1,iS), R_P(iR), &
                 D_Nu_P(iE,iR,iS), I1_Nu_P(iE,iR,iS) )
      END IF

      !WRITE(*,'(A,2I4,5ES12.3)') 'iE, iR, J, H',iE, iR, R_P(iR), &
      !     D_Nu_P(iE,iR,iS), I1_Nu_P(iE,iR,iS), &
      !     Chi(iE,iR,iS), fEQ(iE,iR,iS)

    END DO
    END DO
    END DO

    DO iS = 1, nSpecies
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      Gm_dd_11 &
        = uGF(NodeNumbersX,iX1,iX2,iX3,iGF_Gm_dd_11)

      Gm_dd_22 &
        = uGF(NodeNumbersX,iX1,iX2,iX3,iGF_Gm_dd_22)

      Gm_dd_33 &
        = uGF(NodeNumbersX,iX1,iX2,iX3,iGF_Gm_dd_33)

      DO iE = iE_B0, iE_E0
      DO iNode = 1, nDOF

        iNodeE = ( iE - 1 ) * nDOFE + NodeNumberTable(1,iNode)
        iR     = ( iX1- 1 ) * nDOFX + NodeNumberTable(2,iNode)

        uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS) &
          = D_Nu_P(iNodeE,iR,iS)

        uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,iS) &
          = I1_Nu_P(iNodeE,iR,iS)

        uPR(iNode,iE,iX1,iX2,iX3,iPR_I2,iS) &
          = Zero

        uPR(iNode,iE,iX1,iX2,iX3,iPR_I3,iS) &
          = Zero
        
      END DO

        CALL ComputeConserved_TwoMoment &
               ( uPR(:,iE,iX1,iX2,iX3,iPR_D, iS), &
                 uPR(:,iE,iX1,iX2,iX3,iPR_I1,iS), &
                 uPR(:,iE,iX1,iX2,iX3,iPR_I2,iS), &
                 uPR(:,iE,iX1,iX2,iX3,iPR_I3,iS), &
                 uCR(:,iE,iX1,iX2,iX3,iCR_N, iS), &
                 uCR(:,iE,iX1,iX2,iX3,iCR_G1,iS), &
                 uCR(:,iE,iX1,iX2,iX3,iCR_G2,iS), &
                 uCR(:,iE,iX1,iX2,iX3,iCR_G3,iS), &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33 )

      END DO

    END DO
    END DO
    END DO
    END DO

    DEALLOCATE( E_Nu, R_P, D_P, T_P, Y_P, Chi, fEQ, D_Nu_P, I1_Nu_P )

  END SUBROUTINE InitializeRadiationFields_HomogeneousSphere


  SUBROUTINE ComputeHomogeneousSphereSolution( Chi, f0, R, D, I )

    REAL(DP), INTENT(in)  :: Chi, f0, R
    REAL(DP), INTENT(out) :: D, I

    INTEGER, PARAMETER :: nMu = 8192
    INTEGER            :: iMu
    REAL(DP)           :: Mu(nMu), Distribution(nMu)

    DO iMu = 1, nMu

      Mu(iMu) = - One + Two * DBLE(iMu-1)/DBLE(nMu-1)

      Distribution(iMu) = f_A( R, Mu(iMu), f0, Chi )

    END DO

    D = Half * TRAPEZ( nMu, Mu, Distribution )
    I = Half * TRAPEZ( nMu, Mu, Distribution * Mu )

  END SUBROUTINE ComputeHomogeneousSphereSolution


  REAL(DP) FUNCTION f_A( R, Mu, f0, Chi )

    REAL(DP), INTENT(in) :: R, Mu, f0, Chi

    REAL(DP) :: s

    IF( R < R0 )THEN
      s = ( R * Mu + R0 * SQRT( One - ( R / R0 )**2 * ( One - Mu**2 ) ) )
    ELSE
      IF( ( Mu >= SQRT( One - ( R0 / R )**2 ) ) &
          .and. ( Mu <= One ) )THEN
        s = ( Two * R0 * SQRT( One - ( R / R0 )**2 * ( One - Mu**2 ) ) )
      ELSE
        s = Zero
      END IF
    END IF

    f_A = f0 * ( One - EXP( - Chi * s ) )

    RETURN
  END FUNCTION f_A


  REAL(DP) FUNCTION TRAPEZ( n, x, y )

    INTEGER,  INTENT(in) :: n
    REAL(dp), INTENT(in) :: x(n), y(n)

    INTEGER :: i

    TRAPEZ = 0.0_dp
    DO i = 1, n - 1
      TRAPEZ = TRAPEZ + 0.5_dp * ( x(i+1) - x(i) ) * ( y(i) + y(i+1) )
    END DO

    RETURN
  END FUNCTION TRAPEZ

END MODULE InitializationModule
