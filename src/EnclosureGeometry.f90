
MODULE EnclosureGeometry
!******************************************************************************
!
!  MODULE:        EnclosureGeometry
!
!  PURPOSE:       Reads the enclosure Geometry (vertex and vertices coordinates
!                 data) from a file for use in the program for surface equation
!                 determination
!
!******************************************************************************

USE Global
USE StringUtility
IMPLICIT NONE

CONTAINS

SUBROUTINE CalculateGeometry()

    IMPLICIT NONE
    INTEGER :: I, VertIndex, SurfIndex, TypeIndex, IOS
    CHARACTER (Len = 12) :: ReadStr
    REAL(prec2) TotalReflec

    NVertex = 0
    NSurf = 0

    ! Count verticies and surfaces so we can allocate arrays
    DO
        ReadStr = ''
        READ (2, *, IOSTAT = IOS) ReadStr
        IF (StrLowCase(TRIM(ReadStr)) == "v") THEN
            NVertex = NVertex + 1
        ELSE IF (StrLowCase(TRIM(ReadStr)) == "s") THEN
            NSurf = NSurf + 1
        END IF

        IF (IS_IOSTAT_END(IOS)) THEN
            EXIT
        ENDIF
    END DO

    REWIND(2)

    !  Allocate the size of the array
    ALLOCATE(V(NVertex), XS(NVertex), YS(NVertex), ZS(NVertex), STAT = IOS)
    ALLOCATE(SNumber(NSurf), SVertex(NSurf, NSurf), SType(NSurf), BaseP(NSurf), CMB(NSurf), Emit(NSurf), SurfName(NSurf), STAT = IOS)
    ALLOCATE(DirectionX(NSurf), DirectionY(NSurf), DirectionZ(NSurf), SpecReflec(NSurf), DiffReflec(NSurf), FracSpecEmit(NSurf), FracSpecReflec(NSurf))

    ! Initialize arrays
    DirectionX = 0
    DirectionY = 0
    DirectionZ = 0
    SpecReflec = 0
    DiffReflec = 0
    FracSpecEmit = 0
    FracSpecReflec = 0

201 FORMAT(A1, ' ', I2, 3(' ', f6.3))
203 FORMAT(A1, 5('  ', I2), 1('  ',f6.3), 1('  ',I2), 1('  ',f6.3), A15)

    ! Now, read all data into arrays
    DO
        ReadStr = ''
        READ (2, *, IOSTAT = IOS) ReadStr
        IF (StrLowCase(TRIM(ReadStr)) == "v") THEN

            ! Read in vertex information
            BACKSPACE(2)
            READ (2, *) ReadStr, VertIndex, XS(VertIndex), YS(VertIndex), ZS(VertIndex)
            V(VertIndex) = VertIndex

            ! Write vertex data for RTVT
            IF (WriteLogFile) THEN
                WRITE(4, 201, ADVANCE = 'YES') TRIM(ReadStr), VertIndex, XS(VertIndex), YS(VertIndex), ZS(VertIndex)
            END IF

        ELSE IF (StrLowCase(TRIM(ReadStr)) == "s") THEN

            ! Read in surface information
            BACKSPACE(2)
            READ (2, *) ReadStr, SurfIndex, (SVertex(SurfIndex, I), I = 1, 4), BaseP(SurfIndex), CMB(SurfIndex), Emit(SurfIndex), SurfName(SurfIndex)
            SNumber(SurfIndex) = SurfIndex

            ! Write surface data for RTVT
            IF (WriteLogFile) THEN
                WRITE(4, 203) TRIM(ReadStr), SurfIndex, (SVertex(SurfIndex, I), I = 1, 4), BaseP(SurfIndex), CMB(SurfIndex), Emit(SurfIndex), SurfName(SurfIndex)
            END IF

        ELSE IF (StrLowCase(TRIM(ReadStr)) == "t") THEN

            ! Read in surface type information
            BACKSPACE(2)
            READ(2, *) ReadStr, TypeIndex, SType(TypeIndex)
            BACKSPACE(2)
            IF (SType(TypeIndex) == "SDE") THEN
                READ(2, *) ReadStr, SNumber(TypeIndex), SType(TypeIndex), DirectionX(TypeIndex), DirectionY(TypeIndex), DirectionZ(TypeIndex), FracSpecEmit(TypeIndex), FracSpecReflec(TypeIndex)
            ELSE IF (SType(TypeIndex) == "SDRO") THEN
                READ(2, *) ReadStr, SNumber(TypeIndex), SType(TypeIndex), FracSpecReflec(TypeIndex) ! Reading in specular and diffuse reflection
            ELSE
                READ(2, *) ReadStr, SNumber(TypeIndex), SType(TypeIndex)
            END IF
        END IF

        IF (IS_IOSTAT_END(IOS)) THEN
            EXIT
        ENDIF
    END DO

    ! If no surface types were provided, set them to DIF here
    ! Update reflectance/emittance values and fractions here
    DO I = 1, NSurf
        IF (SType(I) == "SDE") THEN
            TotalReflec = 1 - Emit(I)
            DiffReflec(I) = TotalReflec * (1 - FracSpecReflec(I))
            SpecReflec(I) = TotalReflec * FracSpecReflec(I)
        ELSE IF (SType(I) == "SDRO") THEN
            TotalReflec = 1 - Emit(I)
            DiffReflec(I) = TotalReflec * (1 - FracSpecReflec(I))
            SpecReflec(I) = TotalReflec * FracSpecReflec(I)
        ELSE
            SType(I) = "DIF"
            DiffReflec(I) = 1 - Emit(I)
        END IF
    END DO

    CLOSE(Unit = 2)

END Subroutine  CalculateGeometry

SUBROUTINE CalculateSurfaceEquation()
!******************************************************************************
!
!  SUBROUTINE:    CalculateSurfaceEquation
!
!  PURPOSE:       Determines the coefficients of the surface equation using
!                 surface normal vector a point on the surface. The equation
!                 is of the form Ax + By + Cz + D  = 0
!
!******************************************************************************

!  Calculating the normal vector of the surfaces in the enclosure and the
!  coefficients of the surface equation. The equations is determined in
!  Cartesian coordinate system

    IMPLICIT NONE
    INTEGER :: I, J, K, M, IOS
    INTEGER, DIMENSION (:) :: VS(4)
    REAL(Prec2),  Dimension (4) :: X, Y, Z
    REAL(Prec2), Dimension (:, :) :: V_x(SIndex, 2), V_y(SIndex, 2), V_z(SIndex, 2)

    !  V_x(SIndex, 2)   Vectors on a surface used for normal vector determination
    !  V_y(SIndex, 2)   Vectors on a surface used for normal vector determination
    !  V_z(SIndex, 2)   Vectors on a surface used for normal vector determination
    !  X               x  - coordinate of a vertix
    !  Y               y  - coordinate of a vertix
    !  Z               z  - coordinate of a vertix

    ALLOCATE (SPlane(NSurf), NormalV(NSurf, 3), Width(NSurf), Length(NSurf), Height(NSurf), NormalUV(NSurf, 3), PolygonIndex(NSurf), STAT = IOS)

    !   Assign the vertices of a surfaces their corresponding vertices
    DO J = 1, 4
        VS(J) = SVertex(SIndex, J)
    END DO

    DO J = 1, 4
        IF(VS(4) .ne. 0 .or. J < 4)THEN
            X(J) = XS(VS(J))
            Y(J) = YS(VS(J))
            Z(J) = ZS(VS(J))
        ELSEIF(VS(4) .eq. 0)THEN
            X(4) = XS(VS(1))
            Y(4) = YS(VS(1))
            Z(4) = ZS(VS(1))
       ENDIF
    END DO

    IF(VS(4) == 0)THEN
       PolygonIndex(SIndex) = 3
    ELSE
       PolygonIndex(SIndex) = 4
    ENDIF

    DO I = 1, 2
        V_x(SIndex, I) = X(I + 1) - X(I)
        V_y(SIndex, I) = Y(I + 1) - Y(I)
        V_z(SIndex, I) = Z(I + 1) - Z(I)
    END DO

    CALL SurfaceNormal(V_x, V_y, V_z)

!    Allocate size of the array for coefficients of surface equation
    ALLOCATE (A(NSurf), B(NSurf), C(NSurf), D(NSurf), STAT = IOS)

    DO J = 1, 4
         VS(J) = SVertex(SIndex, J)
         IF(VS(4) .eq. 0)THEN
         ELSE
              X(J) = XS(VS(J))
             Y(J) = YS(VS(J))
             Z(J) = ZS(VS(J))
         ENDIF
    END DO

    ! Calculates the coefficients of the surface equation
    A(SIndex) = NormalUV(SIndex, 1)
    B(SIndex) = NormalUV(SIndex, 2)
    C(SIndex) = NormalUV(SIndex, 3)
    D(SIndex) = - (X(1) * A(SIndex) + Y(1) * B(SIndex) + Z(1) * C(SIndex))

END SUBROUTINE CalculateSurfaceEquation

SUBROUTINE SurfaceNormal(Vx, Vy, Vz)
!******************************************************************************
!
!  PURPOSE:        Determine normal unit vector of the surfaces in the enclosure
!
!
!******************************************************************************
    IMPLICIT NONE
    INTEGER :: I, J, K
    REAL(Prec2) :: NV(SIndex), Vector(3) !Norm_V,
    REAL(Prec2), Dimension (:, :) :: Vx(SIndex, 2), Vy(SIndex, 2), Vz(SIndex, 2)

    ! Norm_V         Magnitude of a vector
    ! NV(SIndex)     Magnitude of a normal vector of a surface SIndex
    ! Vector(3)      Coefficients of a normal vector

    ! Calculates the cross product of the vectors on a surface to determine the
    ! Surface Normal vector
    NormalV(SIndex, 1) = Vy(SIndex, 1) * Vz(SIndex, 2) - Vz(SIndex, 1) * Vy(SIndex, 2)
    NormalV(SIndex, 2) = Vz(SIndex, 1) * Vx(SIndex, 2) - Vx(SIndex, 1) * Vz(SIndex, 2)
    NormalV(SIndex, 3) = Vx(SIndex, 1) * Vy(SIndex, 2) - Vy(SIndex, 1) * Vx(SIndex, 2)

    DO K = 1, 3
        Vector(K) = NormalV(SIndex, K)
    END DO

    ! JDS 11-8-06 attempt to eliminate Norm_V linking problem
    ! NV(SIndex) = Norm_V(Vector)
    NV(Sindex) = sqrt(DOT_PRODUCT(Vector, Vector))

    ! Converts/Normalizes the normal vector to get the unit vector
    DO J = 1, 3
        NormalUV(SIndex, J) = Vector(J) / NV(SIndex)
    END DO

END SUBROUTINE SurfaceNormal

SUBROUTINE CalculateAreaSurfaces()
!******************************************************************************
!
!  PURPOSE:        Determine areas of the surfaces in the enclosure
!
!
!******************************************************************************

    IMPLICIT NONE
    INTEGER :: I, J, IOS
    INTEGER, DIMENSION (:) :: VS(4)
    REAL(Prec2), DIMENSION(:) :: X(4), Y(4), Z(4)
    REAL(prec2), ALLOCATABLE, DIMENSION(:, :) :: LR, LT
    REAL(prec2), ALLOCATABLE, DIMENSION(:) :: S

    !   LR            Length and width of a rectangular surface in the enclosure
    !   LT            The three sides of a triangular surface in the enclosure
    !   S             A parameter used to calculate area for triangular surfaces
    !                 using the Heron's formula s = (LT(1) + LT(2) + LT(3))/2
    !   VS            Vertices of a surface
    !   X, Y & Z      Are coordinates of a vertex


    !   Assign the surfaces their corresponding vertices and coordinates and
    !   and calculate areas of rectangular and triangular polygons

    ALLOCATE(LR(NSurf, 2), LT(NSurf, 3), S(NSurf), Area(NSurf), STAT = IOS)

    IF(PolygonIndex(SIndex) == 4)THEN
        DO J = 1, 4
            VS(J) = SVertex(SIndex, J)
            X(J) = XS(VS(J))
            Y(J) = YS(VS(J))
            Z(J) = ZS(VS(J))
        END DO

        DO I = 1, 2
            LR(SIndex, I) = sqrt((X(I + 1) - X(I))**2 + (Y(I + 1) - Y(I))**2 + (Z(I + 1) - Z(I))**2)
        END DO

        Area(SIndex) = LR(SIndex, 1) * LR(SIndex, 2)
!
    ELSEIF(PolygonIndex(SIndex) == 3)THEN
        DO J = 1, 4
            VS(J) = SVertex(SIndex, J)
            IF(J < 4)THEN
                X(J) = XS(VS(J))
                Y(J) = YS(VS(J))
                Z(J) = ZS(VS(J))
            ELSEIF(J == 4)THEN
                X(4) = XS(VS(1))
                Y(4) = YS(VS(1))
                Z(4) = ZS(VS(1))
            ENDIF
        END DO

        DO J = 1, 3
            LT(SIndex, J) = SQRT((X(J + 1) - X(J))**2 + (Y(J + 1) - Y(J))**2 + (Z(J + 1) - Z(J))**2)
        END DO

        S(SIndex) = (LT(SIndex, 1) + LT(SIndex, 2) + LT(SIndex, 3)) / 2
        Area(SIndex) = SQRT(S(SIndex) * (S(SIndex) - LT(SIndex, 1)) * (S(SIndex) - LT(SIndex, 2)) * (S(SIndex) - LT(SIndex, 3)))
    ENDIF
END SUBROUTINE CalculateAreaSurfaces

SUBROUTINE CrossProduct(Vec1, Vec2, Vec)
!******************************************************************************
!
!  PURPOSE:        Calculates the crossProduct of two vectors
!
!
!******************************************************************************

    REAL(Prec2) :: Vec1(3), Vec2(3), Vec(3)
    Vec(1) = Vec1(2) * Vec2(3) - Vec1(3) * Vec2(2)
    Vec(2) = Vec1(3) * Vec2(1) - Vec1(1) * Vec2(3)
    Vec(3) = Vec1(1) * Vec2(2) - Vec1(2) * Vec2(1)
END SUBROUTINE CrossProduct

Function Norm_V(V)
!******************************************************************************
!
!  PURPOSE:        Calculates the magnitude of a vector
!
!
!******************************************************************************
    IMPLICIT NONE
    REAL(Prec2) :: V(3), Norm_V
    !   V(3)        the vector whose magnitude is to be determined
    !   Norm_V      is the magnitude of the vector V

    Norm_V = 0.0d0
    Norm_V = SQRT(DOT_PRODUCT(V, V))
END Function Norm_V

SUBROUTINE AllocateAndInitArrays()
!******************************************************************************
!
!  PURPOSE:        Allocates the arrays
!
!
!******************************************************************************

    IMPLICIT NONE
    INTEGER :: I, J, IOS

    ALLOCATE(NAEnergy(NSurf, NSurf))
    !ALLOCATE(TCOUNTA(NSurf), TCOUNTR(NSurf), TCOUNTRR(NSurf), NTOTAL(NSurf), STAT = IOS)
    ALLOCATE(TCOUNTA(NSurf), STAT = IOS)
    ALLOCATE(XLS(NSurf), YLS(NSurf), ZLS(NSurf), STAT = IOS)
    ALLOCATE(XP(NSurf, NSurf), YP(NSurf, NSurf), ZP(NSurf, NSurf), Intersection(NSurf, NSurf), STAT = IOS)
    ALLOCATE(Xo(NSurf), Yo(NSurf), Zo(NSurf), Intersects(NSurf), STAT = IOS)
    ALLOCATE (EmittedUV(NSurf, 3), STAT = IOS )

    NAEnergy = 0
    TCOUNTA = 0
    XLS = 0
    YLS = 0
    ZLS = 0
    XP = 0
    YP = 0
    ZP = 0
    Intersection = 0
    Xo = 0
    Yo = 0
    Zo = 0
    EmittedUv = 0
    !ALLOCATE(TSpecA(NSurf), TSpecR(NSurf), TSpecRR(NSurf), NAEnergyS(NSurf, NSurf), NAEnergyR(NSurf, NSurf), NAEnergyWR(NSurf, NSurf))

   !Setting Specular Counter arrays to 0
    !DO I = 1, NSurf  !JH
    !    TSpecA(I) = 0
    !    TSpecR(I) = 0
    !    TSpecRR(I) = 0
    !    DO J = 1, NSurf
    !        NAEnergyS(I, J) = 0
    !        NAEnergyR(I, J) = 0
    !        NAEnergyWR(I, J) = 0
    !    END DO
    !END DO
END SUBROUTINE AllocateAndInitArrays

END MODULE EnclosureGeometry
