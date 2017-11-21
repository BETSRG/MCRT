PROGRAM Main_MonteCarlo
    !Program and Modules created by Bereket Nigusse, Fall 2004 for MAE 5823
    !Program and Modules updated and modified November 2012 by
    !John Holman, Rachel Spitler, and Sudha Sikha for MAE 5823

   USE Global
   USE EnclosureGeometry
   USE EnergyBundleLocation
   USE IntersectionEnergy_Surface
   USE EnergyAbsorbed_Reflected
   USE Distribution_Factors
   USE EnergyBalance
   USE OutPut

   IMPLICIT NONE
   INTEGER		 ::I,J,k,IOS,Index, logfileint
!
!  Initialize the CPU time
   CALL CPU_TIME(TIME1)
!
!  Assign number of Energy Bundles emitted per surface
!   JDS 11-9-2006 Replace this fixed number with an input from
!   a file, set below.
!   NBundles = 1000000

! JDS 11/08/2006 Use generic filenames:
!                                      For input:             input.vs3
!                                      For output:            MCoutput.txt
!                                      For temperatures in K: input.TK
!                                      For number of bundles: parameters.txt
    Open (Unit=2, file='input.vs3',status='unknown',Action='Read',IOSTAT=IOS)
	Open (Unit=3, file='MCoutput.txt',status='unknown',IOSTAT = IOS)    !Diffuse bundles and distribution factors
 	Open (Unit=4, file='logfile.out',status='unknown',IOSTAT = IOS)     !Ray emission, reflection, and absorption points
 	Open (Unit=7, File='input.TK',status='unknown',IOSTAT = IOS)        !Surface temperatures
    OPEN (Unit=6, File='SpecularDF.out',status='unknown', IOSTAT=IOS)   !Total Specular bundles and distribution factors
    OPEN (Unit=9, File='SpecReflecDF.out',status='unknown', IOSTAT=IOS) !Reflected and rereflected Specular bundles and distribution factors
    OPEN (Unit=10, File='SpecWRDF.out',status='unknown', IOSTAT=IOS)    !Non-Reflected AKA absorbed on first intersection Specular bundles and distribution factors
    OPEN (Unit=11, File='DebugFile.txt',status='unknown', IOSTAT=IOS)   !Lists rays that are not finding intersection points and whether they're reflected or not

!   JDS 11-9-2006 Read Nbundles from file so that it can be changed
!   without recompiling.
!   JDS 11-20-2006 Read flag to determine whether or not logfile should
!   be written for use with RTVT

    Open (Unit=8, File='parameters.txt',status='old',IOSTAT = IOS)
    Read (8,*)NBundles
    read (8,*) logfileint
    if (logfileint == 1) then
       WriteLogFile=.true.
       else
       WriteLogFile=.false.
       end if
    Close(Unit = 8)
!
!
    Call CalculateGeometry()
    Call InitializeSeed()
    Call AllocateArrays()
    Call InitializeArrays()
!
    Do SIndex = 1, NSurf
      CALL Calculate_SurfaceEquation()
      CALL Calculate_Area_Surfaces()
      CALL TangentVectors()
    END DO
!
!  Initialize the logical variable for the first emitted energy bundle
   Reflected = .False.
!
    DO SIndexR = 1, NSurf
      SIndex  = SIndexR
!
!  The counter only counts the emitted and absorbed energy
              Ntrials = 0
 If (SurfaceType(SIndex) .EQ. 'SDE') then

    SpIndex = 1 !JH: This begins the Specular ray tracing

 Specular:    DO
!  Calculating source locations for each energy bundle
             Call EnergySourceLocation()
!
!  Calculate the direction of the emitted energy bundle
             CALL DirectionEmittedEnergy()
!
!  Check the intersection points and determine the correct one
             CALL CheckingIntersection()

!  Determine whether the energy bundle is absorbed or reflected
             CALL AbsorptionReflection()

!  If the number of absorbed energy bundles is the same as the number
!  of emitted energy bundles, exit the specular emission loop
!
             If(NTrials==NBundles) Exit Specular
 END DO Specular

 SpIndex = 0 !JH: This begins the diffuse ray tracing for SDE surfaces

Ntrialsd=0
 Diffuse1:    DO
!  Calculating source locations for each energy bundle
             Call EnergySourceLocation()
!
!  Calculate the direction of emitted energy bundle
             CALL DirectionEmittedEnergy()
!
!  Check the intersection points and determine the correct one
             CALL CheckingIntersection()

!  Detrmine whether the enrgy bundle is absorbed or reflected
             CALL AbsorptionReflection()

!  If the number of absorbed energy bundles is the same as the
!  number of emitted energy bundles, exit the diffuse emission loop
!
             If(NTrialsd==NBundles) Exit Diffuse1
 END DO Diffuse1

 ElseIF (SurfaceType(SIndex) .EQ. "DIF" .OR. SurfaceType(SIndex) .EQ. "SDR") THEN
     NTrialsd=0
    !JH: This is the standard diffuse ray tracing
     Diffuse2:    DO

!  Calculating source locations for each energy bundle
             Call EnergySourceLocation()
!
!  Calculate the direction of emitted energy bundle
             CALL DirectionEmittedEnergy()
!
!  Check the intersection points and determine the correct one
             CALL CheckingIntersection()

!  Determine whether the enrgy bundle is absorbed or reflected
             CALL AbsorptionReflection()

!  If the number of absorbed energy bundles is the same as the
!  number of emitted energy bundles, exit the diffuse emission loop
!
             If(NTrialsd==NBundles) Exit Diffuse2
     END DO Diffuse2
 End if
    END DO
!
!  Calculate the radiation distribution factor
    Call Rad_Distribution_Factors()
!
!  Calculate the heat balance of the enclosure
    Call Radiation_Balance
!
!  Calculate the CPU Time
    CALL CPU_TIME(TIME2)

!  Write Results to a file
    Call Print_ViewFactor_HeatFlux

!
     Close(Unit = 3)
     Close(Unit = 4)
     Close(Unit = 7)
     CLOSE(Unit = 6)
     CLOSE(Unit =10)
     CLOSE(Unit =11)
     STOP
END PROGRAM Main_MonteCarlo
