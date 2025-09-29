!WRF:MODEL_LAYER:PHYSICS
!
! SNOWPACK Land Surface Model Driver
! WRF interface driver for SNOWPACK physics integration
!
MODULE module_sf_snowpackdrv

    ! Module-level configuration management
    LOGICAL, SAVE :: config_initialized = .false.
    CHARACTER(LEN=256), SAVE :: config_file_path = './io.ini'
    
    ! C interface for SNOWPACK configuration
    INTERFACE
        SUBROUTINE initialize_snowpack_config_c(ini_file_path) BIND(C, name='initialize_snowpack_config_c')
            USE ISO_C_BINDING
            CHARACTER(KIND=C_CHAR), INTENT(IN) :: ini_file_path(*)
        END SUBROUTINE initialize_snowpack_config_c
    END INTERFACE

CONTAINS

!-------------------------------------------------------------------
! Initialize SNOWPACK configuration (called once)
!-------------------------------------------------------------------
SUBROUTINE initialize_snowpack_module_config()
    IMPLICIT NONE
    
    IF (config_initialized) RETURN
    
    ! Use default io.ini path
    config_file_path = './io.ini'
    
    ! Initialize SNOWPACK configuration with the file
    CALL initialize_snowpack_config_c(TRIM(config_file_path)//C_NULL_CHAR)
    
    config_initialized = .true.
    
    ! Diagnostic message
    CALL wrf_message('*** SNOWPACK LSM: Configuration initialized from ' // TRIM(config_file_path) // ' ***')
    
END SUBROUTINE initialize_snowpack_module_config

!-------------------------------------------------------------------
SUBROUTINE snowpacklsm(ITIMESTEP,XLAT,XLONG,DZ8W,DT,DX, & ! IN : Time/Space-related
               IVGTYP,   ISLTYP,    XLAND,     XICE,XICE_THRESHOLD,       & ! IN : Vegetation/Soil characteristics  
               T_PHY,   QV_CURR,   U_PHY,    V_PHY,   SWDOWN,      GLW,   & ! IN : Atmospheric state
               P8W,RAINBL,                                                & ! IN : Forcing
               TSK,      HFX,      QFX,        LH,   GRDFLX,    SMSTAV, & ! IN/OUT LSM eqv
               SMSTOT,SFCRUNOFF, UDRUNOFF,    ALBEDO,    SNOWC,     SMOIS, & ! IN/OUT LSM eqv
               SH2O,     TSLB,     SNOW,     SNOWH,   CANWAT,    ACSNOM, & ! IN/OUT LSM eqv
               ACSNOW,    EMISS,     QSFC,                                 & ! IN/OUT LSM eqv
               Z0,      ZNT,                                           & ! IN/OUT LSM eqv
               ids,ide,  jds,jde,  kds,kde,                    &
               ims,ime,  jms,jme,  kms,kme,                    &
               its,ite,  jts,jte,  kts,kte     )
!-------------------------------------------------------------------
    USE module_sf_snowpack
    USE module_snowpack_types
    USE ISO_C_BINDING
    IMPLICIT NONE
!-------------------------------------------------------------------

! Parameters
    INTEGER, PARAMETER :: nsoil = 4  ! Number of soil layers (standard WRF Noah)

! IN only
    INTEGER,                                         INTENT(IN   ) ::  ITIMESTEP ! timestep number
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  XLAT      ! latitude [rad]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  XLONG     ! longitude [rad]
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  DZ8W      ! thickness of atmo layers [m]
    REAL,                                            INTENT(IN   ) ::  DT        ! timestep [s]
    REAL,                                            INTENT(IN   ) ::  DX        ! horizontal grid spacing [m]
    INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  IVGTYP    ! vegetation type
    INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  ISLTYP    ! soil type
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  XLAND     ! land mask (1 for land, 2 for water)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  XICE      ! sea ice flag
    REAL,                                            INTENT(IN   ) ::  XICE_THRESHOLD! sea ice threshold

! Atmospheric state
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  T_PHY     ! 3D temperature [K]
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  QV_CURR   ! 3D water vapor [kg/kg]
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  U_PHY     ! 3D u-wind [m/s]
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  V_PHY     ! 3D v-wind [m/s]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SWDOWN    ! downward SW flux at surface [W/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  GLW       ! downward LW flux at surface [W/m2]
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  P8W       ! 3D pressure [Pa]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  RAINBL    ! precipitation rate [mm/s]

! INOUT (with generic LSM equivalent)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  TSK       ! surface radiative temperature [K]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  HFX       ! sensible heat flux [W/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  QFX       ! latent heat flux [kg/s/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  LH        ! latent heat flux [W/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  GRDFLX    ! ground/snow heat flux [W/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SMSTAV    ! soil moisture availability
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SMSTOT    ! total soil water [mm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SFCRUNOFF ! accumulated surface runoff [m]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  UDRUNOFF  ! accumulated subsurface runoff [m]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ALBEDO    ! total grid albedo
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SNOWC     ! snow coverage fraction
    REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(INOUT) ::  SMOIS     ! volumetric soil moisture [m3/m3]
    REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(INOUT) ::  SH2O      ! volumetric liquid soil moisture [m3/m3]
    REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(INOUT) ::  TSLB      ! soil temperature [K]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SNOW      ! snow water equivalent [mm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SNOWH     ! physical snow depth [m]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  CANWAT    ! total canopy water [mm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACSNOM    ! accumulated snow melt [mm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACSNOW    ! accumulated snowfall [mm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  EMISS     ! surface bulk emissivity
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  QSFC      ! bulk surface specific humidity
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  Z0        ! combined z0 sent to coupled model
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ZNT       ! thermal roughness length


    INTEGER, INTENT(IN   )   ::     ids,ide, jds,jde, kds,kde
    INTEGER, INTENT(IN   )   ::     ims,ime, jms,jme, kms,kme
    INTEGER, INTENT(IN   )   ::     its,ite, jts,jte, kts,kte

! Local Variables
    INTEGER :: I, J, K
    TYPE(SnowpackInterface) :: interface_data

!-------------------------------------------------------------------
!  Call SNOWPACK physics (simplified interface without CRYOWRF arrays)
!-------------------------------------------------------------------

    ! Initialize SNOWPACK configuration on first call
    CALL initialize_snowpack_module_config()
    
    ! Diagnostic message
    IF (itimestep == 1) THEN
        CALL wrf_message('*** SNOWPACK LSM: Simplified interface active ***')
    END IF

    ! Grid loop with structured data approach
    DO J = jts,jte
    DO I = its,ite
        ! Skip water points
        IF((XLAND(I,J)-1.5) >= 0.) CYCLE
        
        ! Initialize interface structure
        CALL initialize_snowpack_interface(interface_data, itimestep, dt, dx)
        
        ! Pack WRF data into structured interface
        CALL pack_wrf_to_interface(interface_data, &
            ! Atmospheric inputs
            T_PHY(I,1,J), QV_CURR(I,1,J), U_PHY(I,1,J), V_PHY(I,1,J), P8W(I,1,J), &
            SWDOWN(I,J), GLW(I,J), RAINBL(I,J), &
            ! Grid information  
            XLAT(I,J), XLONG(I,J), DZ8W(I,1,J), &
            ! Surface state (input/output)
            TSK(I,J), HFX(I,J), QFX(I,J), LH(I,J), ALBEDO(I,J), SNOWC(I,J), &
            SNOW(I,J), SNOWH(I,J), QSFC(I,J), Z0(I,J), ZNT(I,J), &
            ! Grid indices
            I, J, 1)
        
        ! Call SNOWPACK physics with structured interface
        CALL snowpack_physics_structured(interface_data)
        
        ! Unpack results back to WRF variables
        CALL unpack_interface_to_wrf(interface_data, &
            TSK(I,J), HFX(I,J), QFX(I,J), LH(I,J), ALBEDO(I,J), SNOWC(I,J), &
            SNOW(I,J), SNOWH(I,J), QSFC(I,J), Z0(I,J), ZNT(I,J))
        
        ! Note: CRYOWRF layer arrays not populated in simplified interface
        ! Layer diagnostics available in interface_data%layers for future use
        
        ! Update accumulated variables from budget tracking
        ACSNOM(I,J) = interface_data%budgets%accum_snowmelt
        ACSNOW(I,J) = interface_data%budgets%accum_snowfall
        
        ! Set standard LSM outputs
        GRDFLX(I,J) = 0.06     ! Geothermal heat flux [W/m2]
        SMSTAV(I,J) = 1.0      ! Soil moisture availability
        SMSTOT(I,J) = 100.0    ! Total soil water [mm]
        SFCRUNOFF(I,J) = SFCRUNOFF(I,J) + 0.0  ! No surface runoff from SNOWPACK yet
        UDRUNOFF(I,J) = UDRUNOFF(I,J) + 0.0    ! No subsurface runoff from SNOWPACK yet
        CANWAT(I,J) = 0.0      ! No canopy water in SNOWPACK
        EMISS(I,J) = 0.98      ! Standard snow/soil emissivity
        
        ! Set soil properties (simplified - SNOWPACK handles snow layers)
        SMOIS(I,1,J) = 0.3  ! Volumetric soil moisture
        SH2O(I,1,J) = 0.3   ! Liquid soil moisture
        TSLB(I,1,J) = TSK(I,J)  ! Soil temperature = surface temperature
        
    END DO
    END DO

END SUBROUTINE snowpacklsm

END MODULE module_sf_snowpackdrv