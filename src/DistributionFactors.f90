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

    ! Delete me
    REAL(Prec2) :: Num, Den

    ALLOCATE(NumEmitted(NSurf), STAT = IOS)
    ALLOCATE(RAD_D_F(NSurf, NSurf), STAT = IOS)
    ALLOCATE(NAEnergy_cmb(NSurf, NSurf), STAT = IOS)
    ALLOCATE(Area_cmb_temp(NSurf), Emit_cmb_temp(NSurf), STAT = IOS)

    RAD_D_F = 0

    ! Populate array for orignial surfaces
    DO I = 1, NSurf
        DO J = 1, NSurf
            ! Diffuse rays
            IF (TCOUNTA(I) .EQ. 0) THEN
                RAD_D_F(I, J) = 0
            ELSE
                RAD_D_F(I, J) = REAL(NAEnergy(I, J)) / REAL(TCOUNTA(I))
            ENDIF
        END DO
    END DO

    ! Delete this section
    OPEN(Unit = 21, File = 'NAEnergy.csv', status = 'unknown', IOSTAT = IOS)
401 FORMAT(I3, 100(', ', I3))
    DO I = 1, NSurf
        WRITE(21, 401)NAEnergy(I, 1), (NAEnergy(I, J), J = 2, NSurf)
    END DO
    CLOSE(Unit = 21)

    ! Delete this section
    OPEN(Unit = 22, File = 'RAD_D_F.csv', status = 'unknown', IOSTAT = IOS)
402 FORMAT(f10.6, 100(', ', f10.6))
    DO I = 1, NSurf
        WRITE(22, 402)RAD_D_F(I, 1), (RAD_D_F(I, J), J = 2, NSurf)
    END DO
    CLOSE(Unit = 22)

    ! Now combine surfaces

    N_SCMB = 0
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
    NumEmitted = 0

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

    ! Delete this section
    OPEN(Unit = 23, File = 'NAEnergy_cmb_col.csv', status = 'unknown', IOSTAT = IOS)
    DO I = 1, NSurf
        WRITE(23, 401)NAEnergy_cmb(I, 1), (NAEnergy_cmb(I, J), J = 2, NSurf)
    END DO
    CLOSE(Unit = 23)

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

    ! Delete this section
    OPEN(Unit = 24, File = 'NAEnergy_cmb_row.csv', status = 'unknown', IOSTAT = IOS)
    DO I = 1, NSurf
        WRITE(24, 401)NAEnergy_cmb(I, 1), (NAEnergy_cmb(I, J), J = 2, NSurf)
    END DO
    CLOSE(Unit = 24)

    ! Copy to new reduced arrays
    K = 0
    DO I = 1, NSurf
        IF (CMB(I) .gt. 0) THEN
            CYCLE
        ELSE
            K = K + 1
            L = 0
            DO J = 1, NSurf
                IF (CMB(J) .gt. 0) THEN
                    CYCLE
                ELSE
                    L = L + 1

                    ! Diffuse rays
                    IF (NumEmitted(I) == 0) THEN
                        RAD_D_F_cmb(K, L) = 0
                    ELSE
                        Num = REAL(NAEnergy_cmb(I, J))
                        Den = REAL(NumEmitted(I))
                        RAD_D_F_cmb(K, L) = REAL(NAEnergy_cmb(I, J)) / REAL(NumEmitted(I))
                    END IF
                END IF
            END DO
        END IF
    END DO

    ! Delete this section
    OPEN(Unit = 25, File = 'RAD_D_F_cmb.csv', status = 'unknown', IOSTAT = IOS)
    DO I = 1, NSurf_cmb
        WRITE(25, 402)RAD_D_F_cmb(I, 1), (RAD_D_F_cmb(I, J), J = 2, NSurf_cmb)
    END DO
    CLOSE(Unit = 25)

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
