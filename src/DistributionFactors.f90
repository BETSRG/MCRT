MODULE DistributionFactors

USE Global
USE EnclosureGeometry
USE EnergyBundleLocation
USE IntersectionEnergySurface
USE EnergyAbsorbedReflected

IMPLICIT NONE
CONTAINS

SUBROUTINE RadDistributionFactors
!******************************************************************************
!
! PURPOSE:            Calculating the radiation distribution factor
!
!
!
!******************************************************************************
    IMPLICIT NONE
    INTEGER :: I, J, K, L, N_SCMB, IOS
    INTEGER, ALLOCATABLE, DIMENSION(:) :: NumEmitted
    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: NAEnergy_cmb
    REAL(Prec2), ALLOCATABLE, DIMENSION(:) :: Area_cmb_temp, Emit_cmb_temp

    ALLOCATE(NumEmitted(NSurf), STAT = IOS)
    ALLOCATE(RAD_D_F(NSurf, NSurf), STAT = IOS)
    ALLOCATE(NAEnergy_cmb(NSurf, NSurf), STAT = IOS)
    ALLOCATE(Area_cmb_temp(NSurf), Emit_cmb_temp(NSurf), STAT = IOS)

    ! Populate array for orignial surfaces
    DO I = 1, NSurf
        DO J = 1, NSurf
            ! Diffuse rays
            IF (TCOUNTA(I) .EQ. 0) THEN
                RAD_D_F(I, J) = 0
            ELSE
                RAD_D_F(I, J) = REAL(NAEnergy(I, J)) / REAL(TCOUNTA(J))
            ENDIF
        END DO
    END DO

    ! Now combine surfaces

    ! Identify number of surface combinations
    DO I = 1, NSurf
        DO J = 1, NSurf
            IF (I == CMB(J))THEN
                N_SCMB = N_SCMB + 1
            ENDIF
        END DO
    END DO

    ! Number of surfaces after combined
    NSurf_cmb = NSurf - N_SCMB

    !ALLOCATE(RAD_D_F_cmb(NSurf_cmb, NSurf_cmb), RAD_D_S_cmb(NSurf_cmb, NSurf_cmb), RAD_D_R_cmb(NSurf_cmb, NSurf_cmb), RAD_D_WR_cmb(NSurf_cmb, NSurf_cmb), STAT = IOS)
    ALLOCATE(RAD_D_F_cmb(NSurf_cmb, NSurf_cmb), STAT = IOS)
    ALLOCATE(Area_cmb(NSurf_cmb), Emit_cmb(NSurf_cmb), STAT = IOS)

    ! Copy over to arrays we can edit
    NAEnergy_cmb = NAEnergy

    ! Initialize arrays
    RAD_D_F_cmb = 0

    ! Combine count
    DO I = 1, NSurf

        IF (CMB(I) .gt. 0) THEN
            NumEmitted(I) = 0
        ELSE
            NumEmitted(I) = TCOUNTA(I)
        ENDIF

        DO J = 1, NSurf
            IF (I == CMB(J)) THEN
                NumEmitted(I) = NumEmitted(I) + TCOUNTA(J)
            ENDIF
        END DO
    END DO

    ! Combine columns
    DO I = 1, NSurf
        DO J = 1, NSurf
            IF (CMB(J) .gt. 0 ) THEN
                ! Diffuse rays
                NAEnergy_cmb(I, CMB(J)) = NAEnergy_cmb(I, CMB(J)) + NAEnergy_cmb(I, J)
                NAEnergy_cmb(I, J) = 0
            ENDIF
        END DO
    END DO

    ! Combine rows
    DO I = 1, NSurf
        DO J = 1, NSurf
            IF (CMB(I) .gt. 0 ) THEN
                ! Diffuse rays
                NAEnergy_cmb(CMB(I), J) = NAEnergy_cmb(CMB(I), J) + NAEnergy_cmb(I, J)
                NAEnergy_cmb(I, J) = 0
            ENDIF
        END DO
    END DO

    ! Copy to new reduced arrays
    K = 0
    DO I = 1, NSurf
        IF (NumEmitted(I) == 0) THEN
            CYCLE
        ELSE
            K = K + 1
            L = 0
            DO J = 1, NSurf
                IF (NumEmitted(J) == 0) THEN
                    CYCLE
                ELSE
                    L = L + 1

                    ! Diffuse rays
                    IF (NumEmitted(I) == 0) THEN
                        RAD_D_F_cmb(K, L) = 0
                    ELSE
                        RAD_D_F_cmb(K, L) = REAL(NAEnergy_cmb(I, J)) / REAL(NumEmitted(I))
                    END IF
                END IF
            END DO
        END IF
    END DO

    ! Combined surface areas
    ! Combined surface emittances
    ! Use area weighting
    Area_cmb = 0
    Area_cmb_temp = 0
    Emit_cmb = 0
    Emit_cmb_temp = 0

    DO I = 1, NSurf
        IF (CMB(I) == 0) THEN
            Area_cmb_temp(I) = Area(I)
            Emit_cmb_temp(I) = Emit(I) * Area(I)
        ELSE
            Area_cmb_temp(CMB(I)) = Area_cmb_temp(CMB(I)) + Area(I)
            Emit_cmb_temp(CMB(I)) = Emit_cmb_temp(CMB(I)) + Emit(I) * Area(I)
        END IF
    END DO

    J = 0
    DO I = 1, NSurf
        IF (Area_cmb_temp(I) == 0) THEN
            CYCLE
        ELSE
            J= J + 1
            Area_cmb(J) = Area_cmb_temp(I)
            Emit_cmb(J) = Emit_cmb_temp(I) / Area_cmb(J)
        END IF
    END DO

    RETURN
END SUBROUTINE RadDistributionFactors

END MODULE DistributionFactors
