PROGRAM MainMonteCarlo
    ! Program and Modules created by Bereket Nigusse, Fall 2004 for MAE 5823
    ! Program and Modules updated and modified November 2012 by
    ! John Holman, Rachel Spitler, and Sudha Sikha for MAE 5823
    ! Matt Mitchell for MAE 5823, December 2017

    USE Global
    USE EnclosureGeometry
    USE EnergyBundleLocation
    USE IntersectionEnergySurface
    USE EnergyAbsorbedReflected
    USE DistributionFactors
    USE EnergyBalance
    USE Output

    IMPLICIT NONE
    INTEGER :: I, J, K, IOS, Index, logfileint

    ! Initialize the CPU time
    CALL CPU_TIME(TIME1)

    OPEN(Unit = 2, file = 'input.vs3', status = 'unknown', Action = 'READ', IOSTAT = IOS)
    OPEN(Unit = 3, file = 'MCoutput.txt', status = 'unknown', IOSTAT = IOS)     ! Diffuse bundles and distribution factors
    OPEN(Unit = 4, file = 'logfile.dat', status = 'unknown', IOSTAT = IOS)      ! Ray emission, reflection, and absorption points
    OPEN(Unit = 7, File = 'input.TK', status = 'unknown', IOSTAT = IOS)         ! Surface temperatures
    OPEN(Unit = 8, File = 'parameters.txt', status = 'old', IOSTAT = IOS)       ! Geometry and ray data for RTVT
    OPEN(Unit = 12, File = 'MCOutput.csv', status = 'unknown', IOSTAT = IOS)     ! csv file with diffuse distribution factors

    ! Read simulation parameters
    READ(8, *) NBundles
    READ(8, *) logfileint
    CLOSE(Unit = 8)

    IF (logfileint == 1) THEN
        WriteLogFile = .true.
    ELSE
        WriteLogFile = .false.
    ENDIF

    WRITE(*, *) "Loading Geometry"
    CALL CalculateGeometry()

    WRITE(*, *) "Initializing Variables"
    CALL InitializeSeed()
    CALL AllocateAndInitArrays()

    WRITE(*, *) "Calculating Surface Areas"
    DO SIndex = 1, NSurf
        CALL CalculateSurfaceEquation()
        CALL CalculateAreaSurfaces()
        CALL TangentVectors()
    END DO

    WRITE(*, *) "Evaluating Surface Energy Bundles"

    DO SIndexRef = 1, NSurf

        ! This surface only reflects, so we don't need to compute emitted values
        IF (SType(SIndexRef) == "SDRO" ) THEN
            CYCLE
        ELSE

            ! Run for all bundles
            DO BIndex = 1, NBundles
                SIndex = SIndexRef
                Reflected = .False.
                RayAbsorbed = .false.
                ReflecCount = 0
                ! Run until absorbed
                DO
                    !  Calculating source locations for each energy bundle
                    CALL EnergySourceLocation()

                    !  Calculate the direction of the emitted energy bundle
                    CALL DirectionEmittedEnergy()

                    ! Check the intersection points and determine the correct one
                    CALL CheckingIntersection()

                    ! Determine whether the energy bundle is absorbed or reflected
                    CALL AbsorptionReflection()

                    IF (RayAbsorbed) THEN
                        EXIT
                    END IF
                END DO
            END DO
        END IF

        ! Update progress bar
        CALL Progress(SIndexRef, NSurf)

    END DO

    ! Calculate the radiation distribution factor
    WRITE(*, *) "Calculating Distribution Factors"
    CALL RadDistributionFactors()

    ! Calculate the heat balance of the enclosure
    WRITE(*, *) "Evaluating Radiation Balance"
    CALL RadiationBalance()

    WRITE(*, *) "Simulaton Complete"

    ! Calculate the CPU Time
    CALL CPU_TIME(TIME2)

    ! Write Results to a file
    CALL PrintViewFactorHeatFlux

    CLOSE(UNIT = 2)
    CLOSE(Unit = 3)
    CLOSE(Unit = 4)
    CLOSE(Unit = 7)
    CLOSE(Unit = 12)
    STOP
END PROGRAM MainMonteCarlo
