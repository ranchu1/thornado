MICROPHYSICS           = WEAKLIB
MOMENT_CLOSURE         = NN_BRONSON
#MOMENT_CLOSURE         = MAXIMUM_ENTROPY_CB
#MOMENT_CLOSURE         = LEVERMORE
#MOMENT_CLOSURE         = MINERBO
#MOMENT_CLOSURE         = KERSHAW_BL
NEUTRINO_MATTER_SOLVER  = EMAB
TWOMOMENT_ORDER         = ORDER_1


THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_Dependencies

VPATH += $(THORNADO_DIR)/SandBox/HomogeneousSphere
VPATH += $(THORNADO_DIR)/SandBox/Interface_FLASH
VPATH += $(THORNADO_DIR)/SandBox/TwoMoment_OrderV

#DEFINES += $(PP)-DTHORNADO_DEBUG

.DEFAULT_GOAL := all

all: HomogeneousSphere1D

HomogeneousSphere: \
	$(weaklib) \
	$(thornado) \
	InitializationModule.o \
	TimeSteppingModule_Flash.o \
	HomogeneousSphere.o
	$(FLINKER) $(FLAGS) -o HomogeneousSphere_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	InitializationModule.o \
	TimeSteppingModule_Flash.o \
	HomogeneousSphere.o \
	$(LIBRARIES)

HomogeneousSphere1D: \
	$(weaklib) \
	$(thornado) \
	InitializationModule.o \
	TwoMoment_TimersModule_OrderV.o \
	TwoMoment_TroubledCellIndicatorModule.o \
        TwoMoment_SlopeLimiterModule_OrderV.o \
	TimeSteppingModule_Flash.o \
	HomogeneousSphere1D.o
	$(FLINKER) $(FLAGS) -o HomogeneousSphere1D_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	InitializationModule.o \
	TwoMoment_TimersModule_OrderV.o \
	TwoMoment_TroubledCellIndicatorModule.o \
	TwoMoment_SlopeLimiterModule_OrderV.o \
	TimeSteppingModule_Flash.o \
	HomogeneousSphere1D.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld
	rm -f *.cub *.ptx

clobber: clean
	rm -f HomogeneousSphere_$(MACHINE)
	rm -f HomogeneousSphere1D_$(MACHINE)

HomogeneousSphere.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  UnitsModule.o \
  ProgramInitializationModule.o \
  TimersModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleX_Lagrange.o \
  ReferenceElementModuleE.o \
  ReferenceElementModuleE_Lagrange.o \
  ReferenceElementModule.o \
  ReferenceElementModule_Lagrange.o \
  GeometryFieldsModule.o \
  GeometryComputationModule.o \
  GeometryFieldsModuleE.o \
  GeometryComputationModuleE.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  EquationOfStateModule_TABLE.o \
  OpacityModule_TABLE.o \
  NeutrinoOpacitiesModule.o \
  NeutrinoOpacitiesComputationModule.o \
  TimeSteppingModule_Flash.o \
  InitializationModule.o \
  InputOutputModuleHDF.o \
  TwoMoment_ClosureModule.o \
  TwoMoment_PositivityLimiterModule.o \
  HomogeneousSphere.F90

HomogeneousSphere1D.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  UnitsModule.o \
  ProgramInitializationModule.o \
  TimersModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleX_Lagrange.o \
  ReferenceElementModuleE.o \
  ReferenceElementModuleE_Lagrange.o \
  ReferenceElementModule.o \
  ReferenceElementModule_Lagrange.o \
  GeometryFieldsModule.o \
  GeometryComputationModule.o \
  GeometryFieldsModuleE.o \
  GeometryComputationModuleE.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  EquationOfStateModule_TABLE.o \
  OpacityModule_TABLE.o \
  TimeSteppingModule_Flash.o \
  InitializationModule.o \
  InputOutputModuleHDF.o \
  Euler_UtilitiesModule.o \
  TwoMoment_ClosureModule.o \
  TwoMoment_SlopeLimiterModule_OrderV.o \
  TwoMoment_PositivityLimiterModule.o \
  HomogeneousSphere1D.F90

InitializationModule.o: \
  KindModule.o \
  UnitsModule.o \
  UtilitiesModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModule.o \
  MeshModule.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  EquationOfStateModule_TABLE.o \
  NeutrinoOpacitiesComputationModule.o \
  TwoMoment_UtilitiesModule.o \
  InitializationModule.F90

TimeSteppingModule_Flash.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  TimersModule.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  TwoMoment_DiscretizationModule_Streaming.o \
  TwoMoment_DiscretizationModule_Collisions_Neutrinos.o \
  TwoMoment_SlopeLimiterModule_OrderV.o \
  TwoMoment_PositivityLimiterModule.o \
  PhysicalConstantsModule.o \
  GeometryFieldsModuleE.o \
  GeometryFieldsModule.o \
  TimeSteppingModule_Flash.F90

