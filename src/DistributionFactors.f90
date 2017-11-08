  MODULE Distribution_Factors
!
!
     USE Global
     USE EnclosureGeometry
     USE EnergyBundleLocation
     USE IntersectionEnergy_Surface
     USE EnergyAbsorbed_Reflected
!
!
     IMPLICIT NONE
  CONTAINS
!
!
     SUBROUTINE Rad_Distribution_Factors
!******************************************************************************
!
! PURPOSE:                        Calculating the radiation distribution factor
!
!
!
!******************************************************************************
        IMPLICIT NONE
        INTEGER                :: I, J, k, l, m, Index, IOS, NEACMB, NAreaCMB, N_C_S_CMB
        INTEGER, ALLOCATABLE, DIMENSION(:) :: NTA, NTR, NTRR, NTRcmb, NTRRcmb, CMBCOUNT
        INTEGER, ALLOCATABLE, DIMENSION(:) :: CMBSURFS, ICOMBSURF, COMBSURF
        INTEGER, ALLOCATABLE, DIMENSION(:, :) :: NAEnergyDummy
!
!        NTA                        =        Number of total energy bundles absorbed in the enclosure
!                                        for an energy emitted from a given surface
!        NTR                        =   Number of total energy bundles reflected in the enclosure
!                                        for energy bundles emitted from a given surface
!
!        NTRR                =   Number of energy bundles re-reflected in the enclosure
!                                        for energy bundles emitted from a given surface
!        NTAcmb                =        Number of total energy bundles absorbed in the enclosure
!                                        for energy bundles emitted from a given surface after
!                                        surface combination
!        NTRcmb                =        Number of energy bundles reflected in the enclosure
!                                        for energy bundles emitted from a given surface after
!                                        surface combination

        ALLOCATE (NTA(NSurf), NTR(NSurf), NTRR(NSurf), COMBSURF(NSurf), &
                  NAEnergyDummy(NSurf, NSurf), STAT=IOS)
!
!   Identify number of surface combinations
        DO J = 1, NSurf
           DO m = 1, NSurf
           IF (J == CMB(m)) Then
              N_SCMB = N_SCMB + 1
           ELSE
           ENDIF
           END DO
        END DO
!
        NSurfcmb = NSurf - N_SCMB ! Number of Surfaces after combined
!
        ALLOCATE (NTAcmb(NSurfcmb), NTRcmb(NSurfcmb), NTRRcmb(NSurfcmb), &
                  NAEnergyCMB(NSurfcmb, NSurfcmb), CMBCOUNT(NSurfcmb), &
                  ICOMBSURF(N_SCMB), CMBSURFS(N_SCMB), AreaCMB(NSurfcmb), STAT=IOS)
!
        DO I = 1, NSurf
           NTA(I) = 0
           NTR(I) = 0
           NTRR(I) = 0
        END DO
!
        DO I = 1, NSurf
           NTA(I) = TCOUNTA(I)
        END DO

        DO m = 1, NSurfcmb
           DO J = 1, NSurfcmb
              NAEnergyCMB(m, J) = 0
           END DO
           NTAcmb(m) = 0
           NTRcmb(m) = 0
           NTRRcmb(m) = 0
        END DO

        DO I = 1, NSurf !Distribution Factors for Diffuse Rays
           DO Index = 1, NSurf
           IF (Real(NTA(I)) .EQ. 0) THEN
              RAD_D_F(I, Index) = 0.0000
           ELSE
              RAD_D_F(I, Index) = NAEnergy(I, Index)/Real(NTA(I))
           END IF
           END DO
        END DO

        DO I = 1, NSurf !Distribution Factors for Specular Rays
           DO Index = 1, NSurf
              IF (Real(TSpecA(I)) .EQ. 0) THEN
                 RAD_D_S(I, Index) = 0.0000
              ELSE
                 RAD_D_S(I, Index) = NAEnergyS(I, Index)/Real(TSpecA(I))
              END IF
           END DO
        END DO

        DO I = 1, NSurf !Distribution Factors for Reflected Specular Rays
           DO Index = 1, NSurf
              IF ((Real(TSpecR(I)) + Real(TSpecRR(I))) .EQ. 0) THEN
                 RAD_D_R(I, Index) = 0.0000
              ELSE
                 RAD_D_R(I, Index) = NAEnergyR(I, Index)/(Real(TSpecR(I)) + Real(TSpecRR(I))) !RS: NAEnergyR is reflected energy
              END IF
           END DO
        END DO

        DO I = 1, NSurf !Distribution Factors for Non-Reflected (those absorbed at the first intersection point) Specular Rays
           DO Index = 1, NSurf
              IF ((REAL(TSpecA(I)) - Real(TSpecR(I))) .EQ. 0) THEN
                 RAD_D_WR(I, Index) = 0.0000
              ELSE
                 RAD_D_WR(I, Index) = NAEnergyWR(I, Index)/(REAL(TSPecA(I)) - Real(TSpecR(I))) !RS: NAEnergyWR is non-reflected energy
              END IF
           END DO
        END DO

        Return
     END SUBROUTINE Rad_Distribution_Factors
!
  END MODULE Distribution_Factors

