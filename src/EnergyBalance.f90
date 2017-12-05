 MODULE EnergyBalance

USE Global
USE EnclosureGeometry
USE EnergyBundleLocation
USE IntersectionEnergySurface
USE EnergyAbsorbedReflected
USE Distribution_Factors

IMPLICIT NONE
CONTAINS

SUBROUTINE RadiationBalance
!******************************************************************************
!
! PURPOSE:          Calculating the net radiation flux at each surface using
!                   the gray view factor or the radiation distribution factor
!
!
!******************************************************************************
    IMPLICIT NONE
    INTEGER  :: I, J, LWL, UPL, IOS
    INTEGER, ALLOCATABLE, DIMENSION(:) :: EB
    REAL(Prec2) :: SIGMA, EBSUM, T

    SIGMA = 5.67E-8    ! Stephan Boltzmann constant
    !   EBSUM   =   Is product sum of emissivities and balck body emissive power
    !               For each surface
    !   LWL     =   The lower surface index for which the temperatures to read is
    !               applicable
    !   UPL     =   The upper surface index for which the temperatures to read is
    !               applicable
    !   T       =   Temperature of the surfaces, K

    ALLOCATE(Ts(NSurf), EB(NSurf), QFLUX(NSurf), Q(NSurf), STAT = IOS)

    ! Read and assign surface temperatures
    DO I = 1, NSurf
        READ(7, *)LWL, UPL, T
        IF(LWL == 0)EXIT
        DO J = LWL, UPL
            Ts(J) = T
        END DO
    END DO

    DO J = 1, NSurf
        EB(J) = SIGMA * (Ts(J)**4)
    END DO

    DO I = 1, NSurf
        EBSUM = 0.0

        DO J = 1, NSurf
            EBSUM = EBSUM + RAD_D_F(I, J) * EB(J)
        END DO

        QFLUX(I) = EMIT(I) * EB(I) - EMIT(I) * EBSUM
        Q(I) = Area(I) * QFlux(I)
    END DO

END SUBROUTINE RadiationBalance

END MODULE EnergyBalance
