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
    INTEGER :: I, J, K, L, IOS
    INTEGER, ALLOCATABLE, DIMENSION(:) :: NumEmitted, NumReflected, NumReReflec, NumEmitted_cmb, NumReflected_cmb, NumReReflect_cmb
    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: NAEnergy_cmb, NAEnergyS_cmb, NAEnergyR_cmb, NAEnergyWR_cmb
    ALLOCATE(NumEmitted(NSurf), STAT = IOS)

    ALLOCATE(RAD_D_F(NSurf, NSurf), RAD_D_S(NSurf, NSurf), RAD_D_R(NSurf, NSurf), RAD_D_WR(NSurf, NSurf), STAT = IOS)
    ALLOCATE(NAEnergy_cmb(NSurf, NSurf), NAEnergyS_cmb(NSurf, NSurf), NAEnergyR_cmb(NSurf, NSurf), NAEnergyWR_cmb(NSurf, NSurf), STAT = IOS)

    ! Populate array for orignial surfaces
    DO I = 1, NSurf
        DO J = 1, NSurf
            ! Diffuse rays
            IF (TCOUNTA(I) .EQ. 0) THEN
                RAD_D_F(I, J) = 0
            ELSE
                RAD_D_F(I, J) = REAL(NAEnergy(I, J)) / REAL(TCOUNTA(J))
            ENDIF

            ! Specular rays
            IF (TSpecA(I) .EQ. 0) THEN
                RAD_D_S(I, J) = 0
            ELSE
                RAD_D_S(I, J) = REAL(NAEnergyS(I, J)) / REAL(TSpecA(I))
            ENDIF

            ! Reflected specular rays
            IF ((TSpecR(I) + TSpecRR(I)) .EQ. 0) THEN
                RAD_D_R(I, J) = 0
            ELSE
                RAD_D_R(I, J) = REAL(NAEnergyR(I, J)) / (REAL(TSpecR(I)) + REAL(TSpecRR(I)))
            ENDIF

            ! Non-Reflected (those absorbed at the first intersection point) specular rays
            IF ((TSpecA(I) - TSpecR(I)) .EQ. 0) THEN
                RAD_D_WR(I, J) = 0
            ELSE
                RAD_D_WR(I, J) = REAL(NAEnergyWR(I, J)) / (REAL(TSpecA(I)) - REAL(TSpecR(I)))
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

    ALLOCATE(RAD_D_F_cmb(NSurf_cmb, NSurf_cmb), RAD_D_S_cmb(NSurf_cmb, NSurf_cmb), RAD_D_R_cmb(NSurf_cmb, NSurf_cmb), RAD_D_WR_cmb(NSurf_cmb, NSurf_cmb), STAT = IOS)

    ! Copy over to arrays we can edit
    NAEnergy_cmb = NAEnergy
    NAEnergyS_cmb = NAEnergyS
    NAEnergyR_cmb = NAEnergyR
    NAEnergyWR_cmb = NAEnergyWR

    ! Initialize arrays
    DO I = 1, NSurf_cmb
      DO J = 1, NSurf_cmb
          RAD_D_F_cmb(I, J) = 0
          RAD_D_S_cmb(I, J) = 0
          RAD_D_R_cmb(I, J) = 0
          RAD_D_WR_cmb(I, J) = 0
      END DO
    END DO

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

                ! Specular rays
                NAEnergyS_cmb(I, CMB(J)) = NAEnergyS_cmb(I, CMB(J)) + NAEnergyS_cmb(I, J)
                NAEnergyS_cmb(I, J) = 0

                ! Reflected specular rays
                NAEnergyR_cmb(I, CMB(J)) = NAEnergyR_cmb(I, CMB(J)) + NAEnergyR_cmb(I, J)
                NAEnergyR_cmb(I, J) = 0

                ! Non-Reflected (those absorbed at the first intersection point) specular rays
                NAEnergyWR_cmb(I, CMB(J)) = NAEnergyWR_cmb(I, CMB(J)) + NAEnergyWR_cmb(I, J)
                NAEnergyWR_cmb(I, J) = 0
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

                ! Specular rays
                NAEnergyS_cmb(CMB(I), J) = NAEnergyS_cmb(CMB(I), J) + NAEnergyS_cmb(I, J)
                NAEnergyS_cmb(I, J) = 0

                ! Reflected specular rays
                NAEnergyR_cmb(CMB(I), J) = NAEnergyR_cmb(CMB(I), J) + NAEnergyR_cmb(I, J)
                NAEnergyR_cmb(I, J) = 0

                ! Non-Reflected (those absorbed at the first intersection point) specular rays
                NAEnergyWR_cmb(CMB(I), J) = NAEnergyWR_cmb(CMB(I), J) + NAEnergyWR_cmb(I, J)
                NAEnergyWR_cmb(I, J) = 0
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

                    ! Specular rays
                    IF (TSpecA(I) == 0) THEN
                        RAD_D_S_cmb(K, L) = 0
                    ELSE
                        RAD_D_S_cmb(K, L) = REAL(NAEnergyS_cmb(I, J)) / REAL(TSpecA(I))
                    END IF

                    ! Reflected specular rays
                    IF ((TSpecR(I) + TSpecRR(I)) == 0) THEN
                        RAD_D_R_cmb(K, L) = 0
                    ELSE
                        RAD_D_R_cmb(K, L) = REAL(NAEnergyR_cmb(I, J)) / (REAL(TSpecR(I)) + REAL(TSpecRR(I)))
                    END IF

                    ! Non-Reflected (those absorbed at the first intersection point) specular rays
                    IF ((TSpecA(I) - TSpecR(I)) == 0) THEN
                        RAD_D_WR_cmb(K, L) = 0
                    ELSE
                        RAD_D_WR_cmb(K, L) = REAL(NAEnergyWR_cmb(I, J)) / (REAL(TSpecA(I)) - REAL(TSpecR(I)))
                    END IF
                END IF
            END DO
        END IF
    END DO

    RETURN
END SUBROUTINE RadDistributionFactors

END MODULE DistributionFactors
