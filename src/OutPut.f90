
MODULE OutPut

    USE Global
    USE EnclosureGeometry
    USE EnergyBundleLocation
    USE IntersectionEnergy_Surface
    USE EnergyAbsorbed_Reflected
    USE Distribution_Factors
    USE EnergyBalance

    IMPLICIT NONE
    CONTAINS

    SUBROUTINE Print_ViewFactor_HeatFlux()
    !******************************************************************************
    !
    ! PURPOSE:          Prints View Factors, Radiation Heat Flux and Heat Transfer
    !                   Rate at Each Surface
    !
    !
    !******************************************************************************
    IMPLICIT NONE
    INTEGER        :: I, J, k, Index

    !  WRITE the Title of the Program and Output data
    WRITE(3, 101)'Monte Carlo Method', 'PURPOSE:', 'Calculates The View &
                 Factors Using Monte Carlo Method', 'and', 'The Net Radiation &
                 Heat Flux at Each Surface'
101 FORMAT(//, 15x, A30, ///, 14x, A25, //, 14x, A52, /, 36x, A3, /, 11x, A50, //)

    DO k = 1, NSurf
        WRITE(3, 1001)NAEnergy(k, :), TCOUNTA(k)
        WRITE(6, 1001)NAEnergyS(k, :), TSpecA(k) !JH !Writing the number of total specular rays absorbed at each surface
        WRITE(9, 1001)NAEnergyR(k, :), TSpecR(k) !Writing the number of reflected specular rays absorbed at each surface
        WRITE(10, 1001)NAEnergyWR(k, :), (TSpecA(k) - TSpecR(k))   !Writing the number of specular rays absorbed on first contact at each surface
    END DO

1001 FORMAT(2x, 100(x, I8), I10)

     WRITE(3, 1002)
     WRITE(6, 1002)
     WRITE(9, 1002)
     WRITE(10, 1002)

1002 FORMAT(//)

     DO Index = 1, NSurf
        WRITE(3, 102)(RAD_D_F(Index, J), J = 1, NSurfcmb)   !Diffuse distribution factors
        WRITE(6, 102)(RAD_D_S(Index, J), J = 1, NSurfcmb)   !Total specular distribution factors
        WRITE(9, 102)(RAD_D_R(Index, J), J = 1, NSurfcmb)   !Reflected specular distribution factors
        WRITE(10, 102)(RAD_D_WR(Index, J), J = 1, NSurfcmb) !Absorbed at first intersection specular distribution factors
     END DO

     !Writing the rest of the outputs for MCOutput.txt
102  FORMAT(4x, 100(2x, f8.6))

     WRITE(3, 103)'Index', 'SURF_NAME', 'Temperature', 'Emissivity', 'Heat Flux', 'Heat Transfer Rate'

103  FORMAT(///, 8x, A5, 2x, A10, 6x, A12, 2x, A12, 4x, A12, 8x, A19)

     DO Index = 1, NSurf
        WRITE(3, 104)Index, SURF_NAME(Index), TS(Index), EMIT(Index), QFLUX(Index), Q(Index)
104     FORMAT(7x, I3, 8x, A12, 4x, F7.2, 8x, F5.2, 8x, ES12.3, 10x, ES12.3)
     END DO

     WRITE(*, 107)'Elapsed Time:', TIME2 - TIME1, 's'
     WRITE(3, 107)'Elapsed Time:', TIME2 - TIME1, 's'

107  FORMAT(//, 8x, A14, 1x, F14.2, x, A1)

    END SUBROUTINE Print_ViewFactor_HeatFlux
 END MODULE OutPut
