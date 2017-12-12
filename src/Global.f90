MODULE Global
!	OKlahoma State University
!	School of Mechanical And Aerospace Engineering
!
!
!	PURPOSE:  Global Data for Program Monte Carlo Method
!

IMPLICIT NONE
SAVE

    INTEGER, PARAMETER :: Prec = SELECTED_REAL_KIND(P = 12)
    INTEGER, PARAMETER :: Prec2 = SELECTED_REAL_KIND(P = 12)

    INTEGER :: NSurf                       ! Number of Surfaces
    INTEGER :: NSurf_cmb                   ! Number of Surfaces after combination
    INTEGER :: SIndex                      ! Surface counting Index
    INTEGER :: SIndexRef                   ! Surface counting Index Reference
    INTEGER :: NVertex                     ! Number of Vertices
    INTEGER :: NBundles                    ! Number of Energy Bundles Emitted

    INTEGER , ALLOCATABLE, DIMENSION(:, :)  :: SVertex      ! Vertices of A surface
    INTEGER , ALLOCATABLE, DIMENSION(:)     :: SNumber      ! Index of a surface
    INTEGER , ALLOCATABLE, DIMENSION(:)     :: V            ! vertex Index
    INTEGER , ALLOCATABLE, DIMENSION(:)     :: SPlane       ! Plane of a Surface (x, y, z)
    INTEGER                                 :: SInter       ! Index of Intercepted Surface
    INTEGER , ALLOCATABLE, DIMENSION(:, :)  :: NAEnergy     ! Absorbed Energy Counter
    INTEGER , ALLOCATABLE, DIMENSION(:) :: TCOUNTA       ! Number of absorbed energy bundle
    INTEGER , ALLOCATABLE, DIMENSION(:) :: TCOUNTR       ! Number of reflected energy bundle
    INTEGER , ALLOCATABLE, DIMENSION(:) :: TCOUNTRR      ! Number of rereflected energy bundle
    INTEGER , ALLOCATABLE, DIMENSION(:) :: NTOTAL        ! Total Number of Energy bundles emitted
    INTEGER , ALLOCATABLE, DIMENSION(:) :: NTACMB        ! Total Number of Energy bundles emitted
                                                         ! after surface combinations

    INTEGER, ALLOCATABLE, DIMENSION(:) :: TSpecA         !Total Number of specular bundles absorbed on first bounce
    INTEGER, ALLOCATABLE, DIMENSION(:) :: TSpecR         !Total Number of specular bundles reflected
    INTEGER, ALLOCATABLE, DIMENSION(:) :: TSpecRR        !Total Number of specular bundles rereflected

    INTEGER , ALLOCATABLE, DIMENSION(:, :)  :: Intersection     ! Surface Intersection Index
    INTEGER , ALLOCATABLE, DIMENSION(:)     :: PolygonIndex     ! 3 is Triangle, 4 is Rectangle
    INTEGER , ALLOCATABLE, DIMENSION(:)     :: CMB              ! Index for surfaces to be combined

    REAL(Prec2), ALLOCATABLE, DIMENSION(:)  :: Emit         ! Emissivities of surfaces
    REAL(Prec2), ALLOCATABLE, DIMENSION(:)  :: Emit_cmb     ! Emissivities of combined surfaces
    REAL(Prec2), ALLOCATABLE, DIMENSION(:)  :: TS           ! surface Temperature, K
    REAL(Prec2), ALLOCATABLE, DIMENSION(:)  :: BaseP        ! Reference Point
    REAL(Prec2)                             :: Rand(7)      ! Random number (0 - 1)
    REAL(Prec2)                             :: TIME1        ! Starting Time in s
    REAL(Prec2)                             :: TIME2        ! Finishing Time in s

    CHARACTER (LEN = 12), ALLOCATABLE, DIMENSION(:) :: SurfName   ! Name of Surfaces
    LOGICAL :: Reflected
    LOGICAL :: ReflectedSpec
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: Intersects ! Surface Intersection Flag
    LOGICAL :: WriteLogFile                          ! Flag to indicate whether log file should be written

    REAL(Prec2), ALLOCATABLE, DIMENSION(:, :) :: XP     ! Intersection Point x - coordinates
    REAL(Prec2), ALLOCATABLE, DIMENSION(:, :) :: YP     ! Intersection Point y - coordinates
    REAL(Prec2), ALLOCATABLE, DIMENSION(:, :) :: ZP     ! Intersection Point z - coordinates
    REAL(Prec2), ALLOCATABLE, DIMENSION(:)   :: SI      ! Scalar Vector Multiplier
    REAL(Prec2), ALLOCATABLE, DIMENSION(:)   :: SIPOS   ! Scalar Vector Multiplier

    REAL(Prec2), ALLOCATABLE, DIMENSION(:) ::   XLS     ! X coordinate of Source Location
    REAL(Prec2), ALLOCATABLE, DIMENSION(:) ::   YLS     ! Y coordinate of Source Location
    REAL(Prec2), ALLOCATABLE, DIMENSION(:) ::   ZLS     ! Z coordinate of Source Location
    REAL(Prec2), ALLOCATABLE, DIMENSION(:) ::   QFLUX   ! Net radiation flux at each surface
    REAL(Prec2), ALLOCATABLE, DIMENSION(:) ::   Q       ! Net radiation heat transfer at each
                                                        ! surface
    REAL(Prec2), ALLOCATABLE, DIMENSION(:) ::   XS      ! x - coordinate of a vertex
    REAL(Prec2), ALLOCATABLE, DIMENSION(:) ::   YS      ! y - coordinate of a vertex
    REAL(Prec2), ALLOCATABLE, DIMENSION(:) ::   ZS      ! z - coordinate of a vertex

    REAL(Prec2), ALLOCATABLE, DIMENSION(:) ::   Xo      ! x - coordinate of intersection point
    REAL(Prec2), ALLOCATABLE, DIMENSION(:) ::   Yo      ! y - coordinate of intersection point
    REAL(Prec2), ALLOCATABLE, DIMENSION(:) ::   Zo      ! z - coordinate of intersection point

    REAL(Prec2), ALLOCATABLE, DIMENSION(:, :) ::  NormalV    ! Normal Vectors of surfaces
    REAL(Prec2), ALLOCATABLE, DIMENSION(:, :) ::  NormalUV   ! Normal Unit Vectors of surfaces
    REAL(Prec2), ALLOCATABLE, DIMENSION(:, :) ::  EmittedUV  ! Unit Vector of emitted energy
    REAL(Prec2), ALLOCATABLE, DIMENSION(:, :) ::  Tan_V1     ! Unit Vector tangent to the source S
    REAL(Prec2), ALLOCATABLE, DIMENSION(:, :) ::  Tan_V2     ! Unit Vector tangent to the source S

    REAL(Prec2), ALLOCATABLE, DIMENSION(:, :) ::  RAD_D_F     ! Diffuse Radiation Distribution Factor
    REAL(Prec2), ALLOCATABLE, DIMENSION(:, :) ::  RAD_D_F_cmb ! Diffuse Radiation Distribution Factor for combined surfaces

    REAL(Prec2), ALLOCATABLE, DIMENSION(:) ::   Area         ! Area of a surface
    REAL(Prec2), ALLOCATABLE, DIMENSION(:) ::   Area_cmb     ! Area of a combined surfaces

    REAL(Prec2), ALLOCATABLE, DIMENSION(:) ::   A            ! Coefficient of X in Surface equation
    REAL(Prec2), ALLOCATABLE, DIMENSION(:) ::   B            ! Coefficient of Y in Surface equation
    REAL(Prec2), ALLOCATABLE, DIMENSION(:) ::   C            ! Coefficient of Z in Surface equation
    REAL(Prec2), ALLOCATABLE, DIMENSION(:) ::   D            ! Constant in Surface equation

    CHARACTER(LEN = 4), ALLOCATABLE, DIMENSION(:) ::  SType ! Surface Type Array
    REAL(Prec2), ALLOCATABLE, DIMENSION(:) ::  DirectionX   ! X Vector Coordinates for SDE type
    REAL(Prec2), ALLOCATABLE, DIMENSION(:) ::  DirectionY   ! Y Vector Coordinates for SDE type
    REAL(Prec2), ALLOCATABLE, DIMENSION(:) ::  DirectionZ   ! Z Vector Coordinates for SDE type
    REAL(Prec2), ALLOCATABLE, DIMENSION(:) ::  SpecReflec   ! Specular Reflectance
    REAL(Prec2), ALLOCATABLE, DIMENSION(:) ::  DiffReflec   ! Diffuse Reflectance

    REAL(Prec2), ALLOCATABLE, DIMENSION(:) ::  FracSpecEmit  !- Fraction of total rays emitted which are specular
    REAL(Prec2), ALLOCATABLE, DIMENSION(:) ::  FracSpecReflec     !- Fraction of total rays reflected which are specular

    INTEGER :: PrevSurf
    INTEGER :: BIndex       !- Bundle index
    LOGICAL :: RayAbsorbed  !- Has bundle been absorbed flag
    INTEGER :: ReflecCount

 END MODULE Global
