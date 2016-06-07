MODULE UnitsModule

  USE KindModule, ONLY: &
    DP
  USE PhysicalConstantsModule, ONLY: &
    SpeedOfLightMKS, &
    GravitationalConstantMKS, &
    BoltzmannConstantMKS, &
    ElectronVoltMKS

  IMPLICIT NONE
  PRIVATE

  REAL(DP), PUBLIC, PARAMETER :: &
    SpeedOfLight          = 1.0_DP, &
    GravitationalConstant = 1.0_DP, &
    BoltzmannConstant     = 1.0_DP, &
    PlanckConstant        = 1.0_DP

  ! --- Lenght ---

  REAL(DP), PUBLIC, PARAMETER :: &
    Meter      = 1.0_DP, &
    Centimeter = 1.0e-2_DP * Meter, &
    Kilometer  = 1.0e+3_DP * Meter

  ! --- Time ---

  REAL(DP), PUBLIC, PARAMETER :: &
    Second      = SpeedOfLightMKS / SpeedOfLight * Meter, &
    Millisecond = 1.0e-3_DP * Second

  ! --- Mass ---

  REAL(DP), PUBLIC, PARAMETER :: &
    Kilogram  = GravitationalConstantMKS / GravitationalConstant &
                  * Meter**3 / Second**2, &
    Gram      = 1.0e-3_DP * Kilogram, &
    SolarMass = 1.98892e30_DP * Kilogram

  ! --- Other Units of Measure ---

  REAL(DP), PUBLIC, PARAMETER :: &
    Joule        = Kilogram * ( Meter / Second )**2, &
    Erg          = Gram * ( Centimeter / Second )**2, &
    Bethe        = 1.0e51_DP * Erg, &
    ElectronVolt = ElectronVoltMKS * Joule, &
    MeV          = 1.0e6_DP * ElectronVolt, &
    Kelvin       = BoltzmannConstantMKS / BoltzmannConstant * Joule, &
    Newton       = Joule / Meter, &
    Dyne         = Erg / Centimeter


  ! --- Units Displayed During Execution and for IO ---

  CHARACTER(16), PRIVATE, PARAMETER :: &
    DisplayLabel_Null          = '', &
    DisplayLabel_Length        = 'km', &
    DisplayLabel_Time          = 'ms', &
    DisplayLabel_Mass          = 'M_sun', &
    DisplayLabel_MassDensity   = 'g/cm^3', &
    DisplayLabel_Energy        = 'MeV', &
    DisplayLabel_EnergyDensity = 'erg/cm^3', &
    DisplayLabel_Pressure      = 'erg/cm^3', &
    DisplayLabel_Temperature   = 'K'
  REAL(DP), PRIVATE, PARAMETER :: &
    DisplayUnit_Length        = Kilometer, &
    DisplayUnit_Time          = Millisecond, &
    DisplayUnit_Mass          = SolarMass, &
    DisplayUnit_MassDensity   = Erg / Centimeter**3, &
    DisplayUnit_Energy        = MeV, &
    DisplayUnit_EnergyDensity = Erg / Centimeter**3, &
    DisplayUnit_Pressure      = Erg / Centimeter**3, &
    DisplayUnit_Temperature   = Kelvin

  TYPE, PRIVATE :: UnitsDisplayType
    LOGICAL  :: &
      Active = .FALSE.
    CHARACTER(16) :: &
      LengthLabel        = DisplayLabel_Null, &
      TimeLabel          = DisplayLabel_Null, &
      MassLabel          = DisplayLabel_Null, &
      MassDensityLabel   = DisplayLabel_Null, &
      EnergyLabel        = DisplayLabel_Null, &
      EnergyDensityLabel = DisplayLabel_Null, &
      PressureLabel      = DisplayLabel_Null, &
      TemperatureLabel   = DisplayLabel_Null
    REAL(DP) :: &
      LengthUnit        = 1.0_DP, &
      TimeUnit          = 1.0_DP, &
      MassUnit          = 1.0_DP, &
      MassDensityUnit   = 1.0_DP, &
      EnergyUnit        = 1.0_DP, &
      EnergyDensityUnit = 1.0_DP, &
      PressureUnit      = 1.0_DP, &
      TemperatureUnit   = 1.0_DP
  END type UnitsDisplayType

  TYPE(UnitsDisplayType), PUBLIC :: UnitsDisplay

  PUBLIC :: ActivateUnitsDisplay

CONTAINS


  SUBROUTINE ActivateUnitsDisplay

    UnitsDisplay % Active = .TRUE.

    UnitsDisplay % LengthLabel        = DisplayLabel_Length
    UnitsDisplay % TimeLabel          = DisplayLabel_Time
    UnitsDisplay % MassLabel          = DisplayLabel_Mass
    UnitsDisplay % MassDensityLabel   = DisplayLabel_MassDensity
    UnitsDisplay % EnergyLabel        = DisplayLabel_Energy
    UnitsDisplay % EnergyDensityLabel = DisplayLabel_EnergyDensity
    UnitsDisplay % PressureLabel      = DisplayLabel_Pressure
    UnitsDisplay % TemperatureLabel   = DisplayLabel_Temperature

    UnitsDisplay % LengthUnit        = DisplayUnit_Length
    UnitsDisplay % TimeUnit          = DisplayUnit_Time
    UnitsDisplay % MassUnit          = DisplayUnit_Mass
    UnitsDisplay % MassDensityUnit   = DisplayUnit_MassDensity
    UnitsDisplay % EnergyUnit        = DisplayUnit_Energy
    UnitsDisplay % EnergyDensityUnit = DisplayUnit_EnergyDensity
    UnitsDisplay % PressureUnit      = DisplayUnit_Pressure
    UnitsDisplay % TemperatureUnit   = DisplayUnit_Temperature

  END SUBROUTINE ActivateUnitsDisplay


END MODULE UnitsModule