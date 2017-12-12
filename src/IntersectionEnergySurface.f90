  MODULE IntersectionEnergySurface
!****************************************************************************
!
!  MODULE:        IntersectionEnergySurface
!
!  PURPOSE:       Determines the point of intersection of the emitted energy &
!                 the surfaces in the enclosure
!
!****************************************************************************

  USE Global
  USE EnclosureGeometry
  USE EnergyBundleLocation
  USE EnergyAbsorbedReflected

  IMPLICIT NONE
  CONTAINS

! Checking intersection point of emitted ray and surfaces in the enclosure
! the emitted ray navigates through the equation of surfaces

SUBROUTINE CheckingIntersection()

!****************************************************************************
!
!  SUBROUTINE:  CheckingIntersection
!
!  PURPOSE:     Determines the point of intersection between the emitted
!               energy ray and the surfaces
!
!  CALLS:       Subroutines IntersectionPoints & SingleOutIntersection
!
!****************************************************************************

    IMPLICIT NONE
    INTEGER :: I, J, K, Index, IOS, InterCount

    CALL IntersectionPoints()
    CALL SingleOutIntersection()

END SUBROUTINE CheckingIntersection

SUBROUTINE IntersectionPoints()
!****************************************************************************
!
!  SUBROUTINE:  IntersectionPoints
!
!  PURPOSE:     Determines all possible points of intersection for the
!               surfaces in the enclosure
!
!****************************************************************************

    IMPLICIT NONE
    INTEGER :: I, J, K, Index, SCount, IOS, InterCount
    INTEGER, DIMENSION (:) :: VS(4)
    REAL(Prec2), DIMENSION(:) :: WV(3), UNV(3), EUV(3), W_V(3)
    REAL(Prec2) :: UNV_DOT_WV, UNV_DOT_EUV
    REAL(Prec2),  Dimension (:) :: X(4), Y(4), Z(4)

    !    SI             Scalar multiplier of emitted energy unit vector to locate
    !                   the intersection point
    !    UNV            Unit vector normal to the surfaces
    !    EUV            Unit vector in the direction of the emitted energy
    !    WV             A vector from a point on a surface intersection with the ray
    !                   to the source point the surface emitting the energy
    !    W_V            Unit vector in the direction of the emitted energy
    !    UNV_DOT_WV     Dot product of UNV and WV vectors
    !    UNV_DOT_EUV    Dot product of UNV and EUV vectors

    ALLOCATE(SI(NSurf), STAT = IOS)
    !  Assign surfaces their corresponding vertices and coordinates
    DO Index = 1, NSurf
        DO J = 1, 4
            VS(J) = SVertex(Index, J)
            IF(VS(4) .ne. 0 )THEN
                X(J) = XS(VS(J))
                Y(J) = YS(VS(J))
                Z(J) = ZS(VS(J))
            ELSEIF(J .lt. 4)THEN
                X(J) = XS(VS(J))
                Y(J) = YS(VS(J))
                Z(J) = ZS(VS(J))
            ELSE
            ENDIF
        END DO

    !  Determine a vector between a point on a surface considered for intersection
    !  and the emitted energy source point

    IF(Index .ne. SIndex) THEN
        WV(1) = - (XLS(SIndex) - X(1))
        WV(2) = - (YLS(SIndex) - Y(1))
        WV(3) = - (ZLS(SIndex) - Z(1))

        !  Determine the dot product of the surfaces unit vector and vector WV
        DO I = 1, 3
            UNV(I) = NormalUV(Index, I)
            W_V(I) = WV(I)
            EUV(I) = EmittedUV(SIndex, I)
        END DO

        UNV_DOT_WV  = DOT_PRODUCT(UNV, W_V)
        UNV_DOT_EUV = DOT_PRODUCT(UNV, EUV)
        SI(Index) = UNV_DOT_WV / UNV_DOT_EUV

        IF (UNV_DOT_EUV .EQ. 0) THEN
            SI(Index) = 0.0  !RS: In the case of division by 0
        ENDIF
    ELSE
        SI(Index) = 0.0
    ENDIF
        DO I = 1, 3
            UNV(I) = 0.0
            W_V(I) = 0.0
            EUV(I) = 0.0
        END DO
    END DO

END SUBROUTINE IntersectionPoints

SUBROUTINE SingleOutIntersection()
!******************************************************************************
!
!  SUBROUTINE:  SingleOutIntersection
!
!  PURPOSE:     Selects the exact intersection points from the possible
!               intersection points
!  USES:        Subroutine IntersectionTriangle(Scount) &
!               IntersectionRectangle(Scount)
!
!******************************************************************************
!
    IMPLICIT NONE
    INTEGER :: I, J, K, Index, Scount, IOS, InterCount
    INTEGER, DIMENSION (:) :: VS(4)
    REAL(Prec2), ALLOCATABLE, DIMENSION(:) :: SIINTER
    REAL(Prec2) SIMIN, SIMAX
    REAL(Prec2) XPVal

    !   SIMIN        the closest intersection distance
    !   SIMAX        Maximum real number
    !   Assign the maximum REAL number to SIMAX

    SIMAX  = 10000000000000000.0

    Allocate(SIINTER(NSurf), STAT = IOS)

    !  Calculates the vector position of the intersection point
    DO Index = 1, NSurf
        IF (Index .ne. SIndex) THEN
            XPVal = XLS(SIndex) + SI(Index) * EmittedUV(SIndex, 1)
            XP(SIndex, Index) = XPVal
            YP(SIndex, Index) = YLS(SIndex) + SI(Index) * EmittedUV(SIndex, 2)
            ZP(SIndex, Index) = ZLS(SIndex) + SI(Index) * EmittedUV(SIndex, 3)
!
            IF (SI(Index) > 0.0) THEN
                Intersection(SIndex, Index) = 1  !0 means no intersection, 1 means there is Inter.
            ELSE
                Intersection(SIndex, Index) = 0
            ENDIF
        ELSE
            Intersection(SIndex, Index) = 0
            Intersects(SIndex) = .FALSE.    !RS: Setting the intersection flag to false for cases when it's the emission surface
        ENDIF
    END DO

    DO Scount = 1, NSurf
        IF(PolygonIndex(Scount) .eq. 4 .and. Intersection(SIndex, Scount) == 1)THEN
            CALL IntersectionRectangle(Scount)
        ELSEIF(PolygonIndex(Scount) .eq. 3 .and. Intersection(SIndex, Scount) == 1)THEN
            CALL IntersectionTriangle(Scount)
        ENDIF

        !  Eliminate intersection point on the back side of emission
        IF(SI(Scount) > 0.0 .and. Intersection(SIndex, Scount) == 1)THen
            SIINTER(Scount) = SI(Scount)
        ELSE
            SIINTER(Scount) = SIMAX
        ENDIF
    END DO

    !  Assign the minimum distance from intersection point
    SIMIN = MINVAL(SIINTER)

    !  Determine intersection by selecting the closest point
    DO I = 1, Nsurf
        IF (Intersects(I))THEN
            IF(SIINTER(I) == SIMIN) THEN
                SInter = I
            ENDIF
        ENDIF
    END DO

END SUBROUTINE SingleOutIntersection

SUBROUTINE IntersectionRectangle(Index)
!******************************************************************************
!
!  SUBROUTINE:  IntersectionRectangle
!
!  PURPOSE:     Finds intersection point (IF any) for rectangular surface
!               JDS: Should also work for any trapezoidal or convex 4-sided
!               polygon
!
!
!******************************************************************************
!
!       Modifications:
!       24 November 2012 - JDS: clean up internal documentation whilst trying to
!                               figure out what is going on!
!
!       Input variables:
!       Index   =  index of surface that is being tested for possible intersection
!       Note: Current ray information is stored in Global variables:
!             Sindex: emitting (or reflecting) surface index
!             Intersection(i, j) = 1 IF the ray emitted from the ith surface intersects the plane of
!             the jth surface; else = 0
!             (JDS: IF this only applies to the current ray, why is it stored in an array?
!             We shouldn't even call this subroutine IF it doesn't intersect.)
!             XP, YP, ZP hold x, y, z coordinates of intersection on the plane, previously determined
!
!    UNV      = Unit normal vector of the rectangular surface
!    V_Int    = Vector from one vertex to the intersection (on plane of surface) point
!    V_edge   = Vector along the edges of the surfaces defined in consistent
!               direction
!    VcpS     = Cross product vector between the edges and intersection vector
!    VcpN     = Dot product of VcpS and the surface unit normal vector

    IMPLICIT NONE
    INTEGER :: I, J, K, Index, SCount, IOS, count
    INTEGER, DIMENSION (:) :: VS(4)
    REAL(Prec2), DIMENSION(:, :) :: VcpS(NSurf, 3), VcpN(NSurf, 4)
    REAL(Prec2), DIMENSION(:) ::  V(3), X(4), Y(4), Z(4), V_edge(3), V_Int(3), Vcp(3), UNV(3), Vedge1(3), Vedge2(3), Vedge3(3), Vedge4(3)
    REAl(Prec2) SIMIN

    !  checks whether the point of intersection of the surface's plane is within the
    !  enclosure
    !  Assign surface its corresponding vertices
    ! (JDS: Shouldn't this be done once globally?)

    DO J = 1, 4
        VS(J) = SVertex(Index, J)
        X(J) = XS(VS(J))
        Y(J) = YS(VS(J))
        Z(J) = ZS(VS(J))
    END DO

    !  Determine a vector for the surface edges using the vertices of the surfaces
    ! (JDS: Shouldn't this be done once globally?)

    IF(Index .ne. SIndex) THEN
        DO J = 1, 4
            IF (J < 4 )THEN
                V_edge(1) = (X(J + 1) - X(J))
                V_edge(2) = (Y(J + 1) - Y(J))
                V_edge(3) = (Z(J + 1) - Z(J))
            ELSEIF(J == 4)THEN
                V_edge(1) = (X(1) - X(4))
                V_edge(2) = (Y(1) - Y(4))
                V_edge(3) = (Z(1) - Z(4))
            ENDIF

            ! Determine a vector from a vertex on the surface to the intersection point on
            ! the plane of the same surface
            V_Int(1) = XP(SIndex, Index) - X(J)
            V_Int(2) = YP(SIndex, Index) - Y(J)
            V_Int(3) = ZP(SIndex, Index) - Z(J)

            CALL CrossProduct(V_edge, V_Int, Vcp)

            DO I = 1, 3
               UNV(I) = NormalUV(Index, I)
            END DO

            VcpN(Index, J) = DOT_PRODUCT(Vcp, UNV)

            DO I = 1, 3
                VcpS(Index, I) = Vcp(I)
            END DO
        END DO
    ENDIF

    !  Eliminate intersection point outside the surface domain
    IF(VcpN(Index, 1)> 0.0 .and. VcpN(Index, 4) > 0.0 .and. VcpN(Index, 2) > 0.0 .and. VcpN(Index, 3) > 0.0) THEN
        SInter = Index
        Intersects(Index) = .True.

    !  Save the intersection point coordinates

       Xo(SInter) = XP(SIndex, Index)
       Yo(SInter) = YP(SIndex, Index)
       Zo(SInter) = ZP(SIndex, Index)

    ! JDS: One possible problem - If intersection is on vertex or edge, it will be "false"

    ELSE
        Intersects(Index) = .false.
        Intersection(SIndex, Index) = 0
    ENDIF
END SUBROUTINE IntersectionRectangle

SUBROUTINE IntersectionTriangle(Index)
!******************************************************************************
!
!  SUBROUTINE:    IntersectionTriangle
!
!  PURPOSE:       Selects the exact intersection points for triangular surfaces
!
!
!******************************************************************************
!
!   UNV      = Unit normal vector of the surfaces
!   V_Int    = Vector from the vertices to the intersection point
!   V_edge   = Vector along the edges of the surfaces defined in consistent
!              direction
!   VcpS     = Cross product vector between the edges and intersection vector

    IMPLICIT NONE
    INTEGER :: I, J, K, Index, SCount, IOS, count
    INTEGER, DIMENSION(:) :: VS(4)
    REAL(Prec2), DIMENSION(:, :):: VcpS(NSurf, 3), VcpN(NSurf, 4)
    REAL(Prec2), DIMENSION(:)::V(3), X(4), Y(4), Z(4), V_edge(3), V_Int(3), Vcp(3), UNV(3)

    !  check whether the point of intersection of the surfaces is within the enclosure
    DO J = 1, 3
        VS(J) = SVertex(Index, J)
        X(J) = XS(VS(J))
        Y(J) = YS(VS(J))
        Z(J) = ZS(VS(J))
    END DO

    !  Determine a vector for the surface edges using the vertices of the surfaces
    IF(Index .ne. SIndex .and.  Intersection(SIndex, Index) == 1) THEN
        DO J = 1, 3
            IF (J < 3 )THEN
                V_edge(1) = (X(J + 1) - X(J))
                V_edge(2) = (Y(J + 1) - Y(J))
                V_edge(3) = (Z(J + 1) - Z(J))
            ELSEIF(J == 3)THEN
                V_edge(1) = (X(1) - X(3))
                V_edge(2) = (Y(1) - Y(3))
                V_edge(3) = (Z(1) - Z(3))
            ENDIF

            !  Determine a vector from a vertex on the surface to the intersection point on
            !  the plane of the same surface
            V_Int(1) = XP(SIndex, Index) - X(J)
            V_Int(2) = YP(SIndex, Index) - Y(J)
            V_Int(3) = ZP(SIndex, Index) - Z(J)

            CALL CrossProduct(V_edge, V_Int, Vcp)

            DO I = 1, 3
                UNV(I) = NormalUV(Index, I)
            END DO

            VcpN(Index, J) = DOT_PRODUCT(Vcp, UNV)

            DO I = 1, 3
                VcpS(Index, I) = Vcp(I)
            END DO
        END DO
    ELSE
    ENDIF

    !  Eliminate intersection point outside the surface domain
    IF(VcpN(Index, 1) > 0.0 .and. VcpN(Index, 2) > 0.0 .and. VcpN(Index, 3) > 0.0 .and. Intersection(SIndex, Index) == 1) THEN
        SInter = Index
        Intersects(Index) = .True.
        Intersection(SIndex, Index) = 1

        !  Save the intersection point coordinates
        Xo(SInter) = XP(SIndex, Index)
        Yo(SInter) = YP(SIndex, Index)
        Zo(SInter) = ZP(SIndex, Index)
    ELSE
        Intersects(Index) = .false.
        Intersection(SIndex, Index) = 0
    ENDIF
END SUBROUTINE IntersectionTriangle

END MODULE IntersectionEnergySurface
