MODULE Distribution_Factors

USE Global
USE EnclosureGeometry
USE EnergyBundleLocation
USE IntersectionEnergy_Surface
USE EnergyAbsorbed_Reflected

IMPLICIT NONE
CONTAINS

SUBROUTINE Rad_Distribution_Factors
!******************************************************************************
!
! PURPOSE:            Calculating the radiation distribution factor
!
!
!
!******************************************************************************
    IMPLICIT NONE
    INTEGER :: I, J, K, L, M, Index, IOS
    INTEGER, ALLOCATABLE, DIMENSION(:) :: NumEmitted, NumReflected, NumReReflec, NumEmitted_cmb, NumReflected_cmb, NumReReflect_cmb

    !   NumEmitted          =    Number of total energy bundles emitted from a given surface
    !   NumReflected        =    Number of total energy bundles reflected from a given surface
    !   NumReReflec         =    Number of energy bundles rereflected rom a given surface
    !   NumEmitted_cmb      =    Number of total energy bundles emitted from a given surface after surface combination
    !   NumReflected_cmb    =    Number of energy bundles reflected ifrom a given surface after surface combination
    !   NumReReflect_cmb    = Number of energy bundles rereflected rom a given surface after surface combination

    ALLOCATE (NumEmitted(NSurf), NumReflected(NSurf), NumReReflec(NSurf), STAT = IOS)

    !   Identify number of surface combinations
    DO J = 1, NSurf
        DO M = 1, NSurf
            IF (J == CMB(M))THEN
                N_SCMB = N_SCMB + 1
            ENDIF
        END DO
    END DO

    NSurf_cmb = NSurf - N_SCMB        ! Number of Surfaces after combined

    ALLOCATE(NumEmitted_cmb(NSurf_cmb), NumReflected_cmb(NSurf_cmb), NumReReflect_cmb(NSurf_cmb), STAT = IOS)
    ALLOCATE(RAD_D_F(NSurf_cmb, NSurf_cmb), RAD_D_S(NSurf_cmb, NSurf_cmb), RAD_D_R(NSurf_cmb, NSurf_cmb), RAD_D_WR(NSurf_cmb, NSurf_cmb), STAT = IOS)

    DO I = 1, NSurf
        NumEmitted(I) = 0
        NumReflected(I) = 0
        NumReReflec(I) = 0
    END DO

    DO I = 1, NSurf_cmb
        NumEmitted(I) = 0
        NumReflected(I) = 0
        NumReReflec(I) = 0
    END DO

    ! Combine count
    DO J = 1, NSurf
        DO I = 1, NSurf
            IF (J == CMB(I)) THEN
                NumEmitted_cmb(J) = NumEmitted_cmb(J) + TCOUNTA(I)
            ENDIF
        END DO

        IF (CMB(J) .gt. 0) THEN
            NumEmitted(J) = 0
        ELSE
            NumEmitted(J) = TCOUNTA(J)
        ENDIF
    END DO

    ! Now combine rows


    !DO M = 1, NSurf_cmb
    !    DO J = 1, NSurf_cmb
    !        NAEnergyCMB(M, J) = 0
    !    END DO
    !
    !    NumEmitted_cmb(M) = 0
    !    NumReflected_cmb(M) = 0
    !    NumReReflect_cmb(M) = 0
    !END DO

    DO I = 1, NSurf  !Distribution Factors for Diffuse Rays
        DO Index = 1, NSurf
            IF (REAL(NumEmitted(I)) .EQ. 0) THEN
                RAD_D_F(I, Index) = 0.0000
            ELSE
                RAD_D_F(I, Index) = NAEnergy(I, Index) / REAL(NumEmitted(I))
            ENDIF
        END DO
    END DO

    DO I = 1, NSurf      !Distribution Factors for Specular Rays
        DO Index = 1, NSurf
            IF (REAL(TSpecA(I)) .EQ. 0) THEN
                RAD_D_S(I, Index) = 0.0000
            ELSE
                RAD_D_S(I, Index) = NAEnergyS(I, Index) / REAL(TSpecA(I))
            ENDIF
        END DO
    END DO

    DO I = 1, NSurf    !Distribution Factors for Reflected Specular Rays
        DO Index = 1, NSurf
            IF ((REAL(TSpecR(I)) + REAL(TSpecRR(I))) .EQ. 0) THEN
                RAD_D_R(I, Index) = 0.0000
            ELSE
                RAD_D_R(I, Index) = NAEnergyR(I, Index) / (REAL(TSpecR(I)) + REAL(TSpecRR(I))) !RS: NAEnergyR is reflected energy
            ENDIF
        END DO
    END DO

    DO I = 1, NSurf    !Distribution Factors for Non-Reflected (those absorbed at the first intersection point) Specular Rays
        DO Index = 1, NSurf
            IF ((REAL(TSpecA(I)) - REAL(TSpecR(I))) .EQ. 0) THEN
                RAD_D_WR(I, Index) = 0.0000
            ELSE
                RAD_D_WR(I, Index) = NAEnergyWR(I, Index) / (REAL(TSpecA(I)) - REAL(TSpecR(I))) !RS: NAEnergyWR is non - reflected energy
            ENDIF
        END DO
    END DO

    RETURN
END SUBROUTINE Rad_Distribution_Factors

END MODULE Distribution_Factors
