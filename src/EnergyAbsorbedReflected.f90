MODULE EnergyAbsorbedReflected


USE Global
USE EnclosureGeometry
USE EnergyBundleLocation

IMPLICIT NONE

CONTAINS

SUBROUTINE AbsorptionReflection()

!******************************************************************************
!
! PURPOSE:            Checking whether the energy bundle absorbed or reflected
!
!
!******************************************************************************
    IMPLICIT NONE
    INTEGER     :: I, J, K, IOS, count
    REAL(Prec2) :: R_incident
    REAL(Prec2) :: XPVal

    !   R_absorbed      Random number generated is used to verify whether the
    !                   intercepted energy is absorbed or reflected by comparing
    !                   it with surface absorptance

    !    JDS 11-10-2006 added all the " IF (WriteLogFile) THEN" blocks to control whether
    !    or not a log file is written.  Also changed FORMAT statements to remove commas
    !    so that RTVT could actually read the file.

    R_incident = Rand(6)
    ReflectedSpec = .false.

    ! Write point data to RTVT file
    IF (WriteLogFile) THEN
        IF (ReflecCount == 0) THEN
101         FORMAT(A1, 1(' ', I2, 3(' ', f6.3)))
            WRITE(4, 101, ADVANCE = 'NO')'P', SIndex, XLS(SIndex), YLS(SIndex), ZLS(SIndex)
        END IF
    END IF

    IF (R_incident .lt. DiffReflec(SInter)) THEN
        ! Reflect Diffusely
        Reflected = .true.
        ReflecCount = ReflecCount + 1
    ELSE IF ((DiffReflec(SInter) .lt. R_incident) .and. (R_incident .lt. (DiffReflec(SInter) + SpecReflec(SInter)))) THEN
        ! Reflect Specularly
        Reflected = .true.
        ReflectedSpec = .true.
        ReflecCount = ReflecCount + 1
    ELSE
        ! Absorb
        RayAbsorbed = .true.
    END IF

    ! Write point data to RTVT file
    IF (WriteLogFile) THEN
        IF (RayAbsorbed) THEN
111         FORMAT(1(' ', I2, 3(' ', f6.3)), ' ' '0')
            XPVal = XP(SIndex, SInter)
            WRITE(4, 111, ADVANCE = 'YES') SInter, XPVal, YP(SIndex, SInter), ZP(SIndex, SInter)
        ELSE
121         FORMAT(1(' ', I2, 3(' ', f6.3)))
            WRITE(4, 121, ADVANCE = 'NO') SInter, XP(SIndex, SInter), YP(SIndex, SInter), ZP(SIndex, SInter)
        END IF
    ENDIF

    PrevSurf = SIndex
    SIndex = SInter

END SUBROUTINE AbsorptionReflection
END MODULE EnergyAbsorbedReflected
