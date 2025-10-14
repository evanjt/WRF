!WRF:MODEL_LAYER:PHYSICS
!
! SNOWPACK Land Surface Model Driver
! WRF interface driver for SNOWPACK physics integration
!
MODULE module_sf_snowpackdrv
    USE ISO_C_BINDING, ONLY: C_CHAR, C_NULL_CHAR, C_DOUBLE, C_INT, C_FLOAT
#ifdef _OPENMP
    USE omp_lib
#endif

    INTERFACE
        FUNCTION snowpack_get_measurement_height_floor() BIND(C, name='snowpack_get_measurement_height_floor')
            USE ISO_C_BINDING
            REAL(C_DOUBLE) :: snowpack_get_measurement_height_floor
        END FUNCTION
        FUNCTION snowpack_get_measurement_height_scale() BIND(C, name='snowpack_get_measurement_height_scale')
            USE ISO_C_BINDING
            REAL(C_DOUBLE) :: snowpack_get_measurement_height_scale
        END FUNCTION
        FUNCTION snowpack_get_snow_roughness_length() BIND(C, name='snowpack_get_snow_roughness_length')
            USE ISO_C_BINDING
            REAL(C_DOUBLE) :: snowpack_get_snow_roughness_length
        END FUNCTION
        FUNCTION snowpack_get_snow_emissivity() BIND(C, name='snowpack_get_snow_emissivity')
            USE ISO_C_BINDING
            REAL(C_DOUBLE) :: snowpack_get_snow_emissivity
        END FUNCTION
        FUNCTION snowpack_get_soil_emissivity() BIND(C, name='snowpack_get_soil_emissivity')
            USE ISO_C_BINDING
            REAL(C_DOUBLE) :: snowpack_get_soil_emissivity
        END FUNCTION
        SUBROUTINE snowpack_set_start_from_file_c(flag) BIND(C, name='snowpack_set_start_from_file')
            USE ISO_C_BINDING
            INTEGER(C_INT), VALUE :: flag
        END SUBROUTINE snowpack_set_start_from_file_c
        SUBROUTINE snowpack_save_state_c(i_grid, j_grid) BIND(C, name='snowpack_save_state')
            USE ISO_C_BINDING
            INTEGER(C_INT), VALUE :: i_grid, j_grid
        END SUBROUTINE snowpack_save_state_c
    END INTERFACE

    ! Module-level configuration management
    LOGICAL, SAVE :: config_initialized = .false.
    LOGICAL, SAVE :: config_logged = .false.
    LOGICAL, SAVE :: coupling_enabled = .TRUE.
    LOGICAL, SAVE :: coupling_disabled_logged = .FALSE.
    LOGICAL, SAVE :: multi_layer_soil_enabled = .FALSE.
    CHARACTER(LEN=256), SAVE :: config_file_path = './io.ini'

    ! Per-grid accumulation for Meteo averaging / throttling
    LOGICAL, SAVE :: accum_initialized = .FALSE.
    INTEGER, SAVE :: snowpack_compute_counter = 1
    REAL, ALLOCATABLE, SAVE :: acc_ta(:,:), acc_qv(:,:), acc_u(:,:), acc_v(:,:)
    REAL, ALLOCATABLE, SAVE :: acc_sw(:,:), acc_lw(:,:), acc_p(:,:), acc_rain(:,:)
    INTEGER, ALLOCATABLE, SAVE :: acc_samples(:,:)
    
    ! C interface for SNOWPACK configuration and layer extraction
    INTERFACE
        SUBROUTINE initialize_snowpack_config_c(ini_file_path) BIND(C, name='initialize_snowpack_config_c')
            USE ISO_C_BINDING
            CHARACTER(KIND=C_CHAR), INTENT(IN) :: ini_file_path(*)
        END SUBROUTINE initialize_snowpack_config_c
        
        SUBROUTINE initialize_wrf_simulation_time_c(start_year, start_month, start_day, &
            start_hour, start_minute) BIND(C, name='initialize_wrf_simulation_time_c')
            USE ISO_C_BINDING
            INTEGER(C_INT), VALUE :: start_year, start_month, start_day, start_hour, start_minute
        END SUBROUTINE initialize_wrf_simulation_time_c
        
        SUBROUTINE extract_snowpack_layers_c(i_grid, j_grid, &
            layer_temps, layer_thick, layer_voli, layer_volw, layer_volv, &
            layer_rg, layer_rb, layer_dd, layer_sp, layer_cdot, layer_meta, &
            layer_depd, layer_graintype, layer_mk, n_layers) &
            BIND(C, name='extract_snowpack_layers_c')
            USE ISO_C_BINDING
            INTEGER(C_INT), INTENT(IN) :: i_grid, j_grid
            REAL(C_FLOAT), INTENT(OUT) :: layer_temps(100), layer_thick(100)
            REAL(C_FLOAT), INTENT(OUT) :: layer_voli(100), layer_volw(100), layer_volv(100)
            REAL(C_FLOAT), INTENT(OUT) :: layer_rg(100), layer_rb(100), layer_dd(100), layer_sp(100)
            REAL(C_FLOAT), INTENT(OUT) :: layer_cdot(100), layer_meta(100), layer_depd(100)
            REAL(C_FLOAT), INTENT(OUT) :: layer_graintype(100), layer_mk(100)
            INTEGER(C_INT), INTENT(OUT) :: n_layers
        END SUBROUTINE extract_snowpack_layers_c

    END INTERFACE

CONTAINS

SUBROUTINE snowpack_tolower(str)
    CHARACTER(LEN=*), INTENT(INOUT) :: str
    INTEGER :: idx, ich

    DO idx = 1, LEN(str)
        ich = IACHAR(str(idx:idx))
        IF (ich >= IACHAR('A') .AND. ich <= IACHAR('Z')) THEN
            str(idx:idx) = CHAR(ich + 32)
        END IF
    END DO
END SUBROUTINE snowpack_tolower

!-------------------------------------------------------------------
! Initialize SNOWPACK configuration (called once)
!-------------------------------------------------------------------
SUBROUTINE initialize_snowpack_module_config()
    USE module_configure  ! For nl_get_* functions
    IMPLICIT NONE
    
    ! Local variables
    CHARACTER(LEN=256) :: snowpack_config_path
    INTEGER :: domain_id
    
    LOGICAL :: already_configured
    LOGICAL :: enable_coupling_local
    LOGICAL :: multi_layer_local

    already_configured = .FALSE.

!$OMP CRITICAL(snowpack_config_init)
    IF (config_initialized) THEN
        already_configured = .TRUE.
    ELSE
        ! Get SNOWPACK configuration file path from namelist
        domain_id = 1  ! Default to domain 1 for initialization
        CALL nl_get_snowpack_config_file(domain_id, snowpack_config_path)

        ! Use namelist-specified path or default
        IF (LEN_TRIM(snowpack_config_path) > 0) THEN
            config_file_path = TRIM(snowpack_config_path)
        ELSE
            config_file_path = './io.ini'
        ENDIF

        CALL nl_get_snpack_enable_coupling(domain_id, enable_coupling_local)
        coupling_enabled = enable_coupling_local

        CALL nl_get_snpack_multi_layer_soil(domain_id, multi_layer_local)
        multi_layer_soil_enabled = multi_layer_local

        ! Initialize SNOWPACK configuration with the file
        CALL initialize_snowpack_config_c(TRIM(config_file_path)//C_NULL_CHAR)

        ! Note: WRF time initialization will be done later in snowpacklsm()
        ! when namelist is fully available (CRYOWRF pattern)

        config_initialized = .TRUE.

        ! Diagnostic message
        CALL wrf_message('*** SNOWPACK LSM: Initialized with config file: ' // TRIM(config_file_path) // ' ***')
        IF (.NOT. coupling_enabled) THEN
            CALL wrf_message('*** SNOWPACK LSM: Coupling disabled via snpack_enable_coupling ***')
        END IF
        IF (multi_layer_soil_enabled) THEN
            CALL wrf_message('*** SNOWPACK LSM: Multi-layer soil export enabled ***')
        END IF
    END IF
!$OMP END CRITICAL(snowpack_config_init)

    IF (already_configured) RETURN
    
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
               SN_T,SN_THICK,SN_VOLI,SN_VOLW,SN_VOLV,SN_RG,SN_RB,SN_DD,SN_SP,SN_NLAYERS,SN_CDOT,SN_META,SN_DEPD,SN_GRAINTYPE,SN_MK, & ! OPTIONAL CRYOWRF layer arrays
               M_BUDG_PRECIP,M_BUDG_SUBLIM,M_BUDG_MELT,M_BUDG_SWE,M_BUDG_REFREEZE, & ! OPTIONAL CRYOWRF mass budget
               E_BUDG_LW_IN,E_BUDG_LW_OUT,E_BUDG_SW_IN,E_BUDG_SW_OUT, & ! OPTIONAL CRYOWRF energy budget
               E_BUDG_SENSIBLE,E_BUDG_LATENT,E_BUDG_GFLUX,E_BUDG_RAINE,E_BUDG_DU,  &
               ids,ide,  jds,jde,  kds,kde,                    &
               ims,ime,  jms,jme,  kms,kme,                    &
               its,ite,  jts,jte,  kts,kte     )
!-------------------------------------------------------------------
    USE module_sf_snowpack
    USE module_snowpack_types
    USE module_configure  ! For nl_get_* functions
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

! CRYOWRF layer arrays (OPTIONAL - matching CRYOWRF behavior)
    REAL,    DIMENSION( ims:ime, 1:100, jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_T        ! layer temperature [K]
    REAL,    DIMENSION( ims:ime, 1:100, jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_THICK    ! layer thickness [m]
    REAL,    DIMENSION( ims:ime, 1:100, jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_VOLI     ! ice volume fraction [%]
    REAL,    DIMENSION( ims:ime, 1:100, jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_VOLW     ! water volume fraction [%]
    REAL,    DIMENSION( ims:ime, 1:100, jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_VOLV     ! air volume fraction [%]
    REAL,    DIMENSION( ims:ime, 1:100, jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_RG       ! grain radius [mm]
    REAL,    DIMENSION( ims:ime, 1:100, jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_RB       ! bond radius [mm]
    REAL,    DIMENSION( ims:ime, 1:100, jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_DD       ! dendricity [-]
    REAL,    DIMENSION( ims:ime, 1:100, jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_SP       ! sphericity [-]
    REAL,    DIMENSION( ims:ime, 1:100, jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_CDOT     ! stress rate [-]
    REAL,    DIMENSION( ims:ime, 1:100, jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_META     ! metamorphism tracker [-]
    REAL,    DIMENSION( ims:ime, 1:100, jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_DEPD     ! deposition date (Julian)
    REAL,    DIMENSION( ims:ime, 1:100, jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_GRAINTYPE ! grain classification
    REAL,    DIMENSION( ims:ime, 1:100, jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_MK       ! grain marker
    INTEGER, DIMENSION( ims:ime,        jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_NLAYERS  ! number of layers

! CRYOWRF mass budget arrays (OPTIONAL - matching CRYOWRF behavior)
    REAL,    DIMENSION( ims:ime,        jms:jme ), INTENT(INOUT), OPTIONAL ::  M_BUDG_PRECIP    ! mass budget precip [kg m-2]
    REAL,    DIMENSION( ims:ime,        jms:jme ), INTENT(INOUT), OPTIONAL ::  M_BUDG_SUBLIM    ! mass budget sublim [kg m-2]
    REAL,    DIMENSION( ims:ime,        jms:jme ), INTENT(INOUT), OPTIONAL ::  M_BUDG_MELT      ! mass budget melt [kg m-2]
    REAL,    DIMENSION( ims:ime,        jms:jme ), INTENT(INOUT), OPTIONAL ::  M_BUDG_SWE       ! mass budget SWE change [kg m-2]
    REAL,    DIMENSION( ims:ime,        jms:jme ), INTENT(INOUT), OPTIONAL ::  M_BUDG_REFREEZE  ! mass budget refreeze [kg m-2]

! CRYOWRF energy budget arrays (OPTIONAL - matching CRYOWRF behavior)
    REAL,    DIMENSION( ims:ime,        jms:jme ), INTENT(INOUT), OPTIONAL ::  E_BUDG_LW_IN     ! energy budget LW in  [J m-2]
    REAL,    DIMENSION( ims:ime,        jms:jme ), INTENT(INOUT), OPTIONAL ::  E_BUDG_LW_OUT    ! energy budget LW out [J m-2]
    REAL,    DIMENSION( ims:ime,        jms:jme ), INTENT(INOUT), OPTIONAL ::  E_BUDG_SW_IN     ! energy budget SW in  [J m-2]
    REAL,    DIMENSION( ims:ime,        jms:jme ), INTENT(INOUT), OPTIONAL ::  E_BUDG_SW_OUT    ! energy budget SW out [J m-2]
    REAL,    DIMENSION( ims:ime,        jms:jme ), INTENT(INOUT), OPTIONAL ::  E_BUDG_SENSIBLE  ! energy budget sensible [J m-2]
    REAL,    DIMENSION( ims:ime,        jms:jme ), INTENT(INOUT), OPTIONAL ::  E_BUDG_LATENT    ! energy budget latent   [J m-2]
    REAL,    DIMENSION( ims:ime,        jms:jme ), INTENT(INOUT), OPTIONAL ::  E_BUDG_GFLUX     ! energy budget ground flux [J m-2]
    REAL,    DIMENSION( ims:ime,        jms:jme ), INTENT(INOUT), OPTIONAL ::  E_BUDG_RAINE     ! energy budget rain       [J m-2]
    REAL,    DIMENSION( ims:ime,        jms:jme ), INTENT(INOUT), OPTIONAL ::  E_BUDG_DU        ! internal energy change   [J m-2]

    INTEGER, INTENT(IN   )   ::     ids,ide, jds,jde, kds,kde
    INTEGER, INTENT(IN   )   ::     ims,ime, jms,jme, kms,kme
    INTEGER, INTENT(IN   )   ::     its,ite, jts,jte, kts,kte

! Local Variables
    INTEGER :: I, J, K
    INTEGER :: layer_count_local
    INTEGER :: write_interval_steps
    INTEGER :: layer_limit
    INTEGER :: no_snpack_lay_to_sav
    INTEGER :: n_layers_local
    INTEGER :: copy_count
    LOGICAL :: should_write_layers
    LOGICAL :: cell_active
    LOGICAL :: save_state_this_step
    REAL, PARAMETER :: SN_LAYER_SENTINEL = -999.0
    REAL(C_FLOAT) :: layer_temp_local(100)
    REAL(C_FLOAT) :: layer_thick_local(100)
    REAL(C_FLOAT) :: layer_voli_local(100)
    REAL(C_FLOAT) :: layer_volw_local(100)
    REAL(C_FLOAT) :: layer_volv_local(100)
    REAL(C_FLOAT) :: layer_rg_local(100)
    REAL(C_FLOAT) :: layer_rb_local(100)
    REAL(C_FLOAT) :: layer_dd_local(100)
    REAL(C_FLOAT) :: layer_sp_local(100)
    REAL(C_FLOAT) :: layer_cdot_local(100)
    REAL(C_FLOAT) :: layer_meta_local(100)
    REAL(C_FLOAT) :: layer_depd_local(100)
    REAL(C_FLOAT) :: layer_graintype_local(100)
    REAL(C_FLOAT) :: layer_mk_local(100)
    TYPE(SnowpackInterface) :: interface_data
    CHARACTER(LEN=256) :: snpack_mode_lower
    CHARACTER(LEN=256) :: snpack_mode_trim
    
    ! Time initialization variables (CRYOWRF pattern)
    INTEGER(C_INT), SAVE :: simulation_start_year, simulation_start_month, simulation_start_day
    INTEGER(C_INT), SAVE :: simulation_start_hour, simulation_start_minute
    LOGICAL, SAVE :: time_initialized = .FALSE.
    INTEGER :: domain_id
    
    ! SNOWPACK namelist parameters
    LOGICAL :: sn_start_from_file
    REAL :: snpack_dt, snpack_write_dt
    CHARACTER(LEN=256) :: snpack_mode
    INTEGER :: thread_id
    INTEGER :: num_threads
    REAL :: ta_avg_val, qv_avg_val, u_avg_val, v_avg_val
    REAL :: sw_avg_val, lw_avg_val, p_avg_val, rain_total
    REAL :: dt_effective, met_height

!-------------------------------------------------------------------
!  Call SNOWPACK physics (simplified interface without CRYOWRF arrays)
!-------------------------------------------------------------------

    ! Initialize SNOWPACK configuration on first call
    CALL initialize_snowpack_module_config()
    IF (.NOT. coupling_enabled) THEN
        IF (.NOT. coupling_disabled_logged .AND. itimestep == 1) THEN
            CALL wrf_message('SNOWPACK LSM: Coupling disabled via snpack_enable_coupling; skipping SNOWPACK physics')
            coupling_disabled_logged = .TRUE.
        END IF
        RETURN
    END IF
    
    ! Get domain ID from grid dimensions (WRF pattern)
    IF (ids == 1) THEN
        domain_id = 1  ! Parent domain
    ELSE
        ! For nested domains, use a simple mapping
        domain_id = 1  ! Default to domain 1 for now
    ENDIF
    
    ! Get SNOWPACK namelist parameters for this domain
    CALL nl_get_sn_start_from_file(domain_id, sn_start_from_file)
    CALL nl_get_snpack_dt(domain_id, snpack_dt)
    CALL nl_get_snpack_write_dt(domain_id, snpack_write_dt)
    CALL nl_get_snpack_mode(domain_id, snpack_mode)
    CALL nl_get_no_snpack_lay_to_sav(domain_id, no_snpack_lay_to_sav)

    snpack_mode_lower = snpack_mode
    CALL snowpack_tolower(snpack_mode_lower)
    IF (LEN_TRIM(snpack_mode_lower) == 0) THEN
        snpack_mode_trim = 'land'
    ELSE
        snpack_mode_trim = snpack_mode_lower
    END IF

    layer_limit = MAX(1, MIN(100, no_snpack_lay_to_sav))
    IF (snpack_write_dt > 0.0) THEN
        write_interval_steps = MAX(1, NINT(snpack_write_dt / dt))
    ELSE
        write_interval_steps = 1
    END IF
    should_write_layers = (MOD(itimestep, write_interval_steps) == 0)
    save_state_this_step = should_write_layers .AND. sn_start_from_file

    CALL snowpack_set_start_from_file_c(MERGE(1_C_INT, 0_C_INT, sn_start_from_file))

    ! Initialize accumulation buffers (per grid point) on first entry
#ifdef _OPENMP
    !$OMP CRITICAL(snowpack_accum_init)
#endif
    IF (.NOT. accum_initialized) THEN
        ! Mirror CRYOWRF's guard that only seeds the accumulators once (Coupler.cpp:820ff)
        ALLOCATE(acc_ta(ims:ime, jms:jme), acc_qv(ims:ime, jms:jme), acc_u(ims:ime, jms:jme), acc_v(ims:ime, jms:jme))
        ALLOCATE(acc_sw(ims:ime, jms:jme), acc_lw(ims:ime, jms:jme), acc_p(ims:ime, jms:jme), acc_rain(ims:ime, jms:jme))
        ALLOCATE(acc_samples(ims:ime, jms:jme))
        acc_ta = 0.0; acc_qv = 0.0; acc_u = 0.0; acc_v = 0.0
        acc_sw = 0.0; acc_lw = 0.0; acc_p = 0.0; acc_rain = 0.0
        acc_samples = 0
        accum_initialized = .TRUE.
    END IF
#ifdef _OPENMP
    !$OMP END CRITICAL(snowpack_accum_init)
#endif

    ! Determine how often Snowpack should run (relative to WRF timestep)
    snowpack_compute_counter = MAX(1, NINT((snpack_dt*60.0) / dt))
    
    ! Initialize WRF simulation time on first call (CRYOWRF pattern)
    IF (.NOT. time_initialized) THEN
        ! Initialize variables to detect if nl_get functions fail
        simulation_start_year   = -999
        simulation_start_month  = -999
        simulation_start_day    = -999
        simulation_start_hour   = -999
        simulation_start_minute = -999
        
        ! Namelist values will be stored in the variables above
        ! Get actual namelist values using WRF's nl_get_* functions (like CRYOWRF)
        CALL nl_get_start_year   ( domain_id, simulation_start_year   )
        
        CALL nl_get_start_month  ( domain_id, simulation_start_month  )
        
        CALL nl_get_start_day    ( domain_id, simulation_start_day    )
        
        CALL nl_get_start_hour   ( domain_id, simulation_start_hour   )
        
        CALL nl_get_start_minute ( domain_id, simulation_start_minute )
          
        ! Validate the retrieved values
        IF (simulation_start_year < 1900 .OR. simulation_start_year > 2100 .OR. &
            simulation_start_month < 1 .OR. simulation_start_month > 12 .OR. &
            simulation_start_day < 1 .OR. simulation_start_day > 31 .OR. &
            simulation_start_hour < 0 .OR. simulation_start_hour > 23 .OR. &
            simulation_start_minute < 0 .OR. simulation_start_minute > 59) THEN
            
            CALL wrf_error_fatal('SNOWPACK-FATAL: Invalid namelist time values - nl_get_* functions failed')
        ENDIF
        
        CALL initialize_wrf_simulation_time_c(simulation_start_year, simulation_start_month, &
                                              simulation_start_day, simulation_start_hour, &
                                              simulation_start_minute)
        time_initialized = .TRUE.
    ENDIF
    
    ! Diagnostic message showing namelist parameters
    IF (itimestep == 1 .AND. .NOT. config_logged) THEN
        CALL wrf_message('*** SNOWPACK LSM: Simplified interface active ***')
        WRITE(*,'(A,L1)') 'SNOWPACK-INFO: sn_start_from_file = ', sn_start_from_file
        WRITE(*,'(A,F8.1)') 'SNOWPACK-INFO: snpack_dt = ', snpack_dt
        WRITE(*,'(A,F8.1)') 'SNOWPACK-INFO: snpack_write_dt = ', snpack_write_dt
        WRITE(*,'(A,A)') 'SNOWPACK-INFO: snpack_mode = ', TRIM(snpack_mode_trim)
        config_logged = .TRUE.
    END IF

    ! Safety check: Verify loop bounds are reasonable
    IF (its < 0 .OR. ite < 0 .OR. jts < 0 .OR. jte < 0) THEN
        WRITE(*,'(A,4(I0,1X))') 'SNOWPACK-ERROR: Negative loop bounds detected! its,ite,jts,jte = ', its, ite, jts, jte
        RETURN
    ENDIF

    IF (its > ite .OR. jts > jte) THEN
        WRITE(*,'(A,4(I0,1X))') 'SNOWPACK-ERROR: Invalid loop bounds (start > end)! its,ite,jts,jte = ', its, ite, jts, jte
        RETURN
    ENDIF

    IF (ite > 10000 .OR. jte > 10000) THEN
        WRITE(*,'(A,4(I0,1X))') 'SNOWPACK-ERROR: Loop bounds too large! its,ite,jts,jte = ', its, ite, jts, jte
        RETURN
    ENDIF

#ifdef _OPENMP
    thread_id = omp_get_thread_num()
    num_threads = omp_get_num_threads()
#else
    thread_id = 0
    num_threads = 1
#endif

    IF (itimestep == 1 .AND. thread_id == 0) THEN
        WRITE(*,'(A,I4,A,4I10,A,I4)') 'SNOWPACK-TRACE tid=', thread_id, &
     &        ' its,ite,jts,jte=', its, ite, jts, jte, ' nthreads=', num_threads
    ENDIF

    ! Grid loop with structured data approach
    DO J = jts,jte
    DO I = its,ite
        ! Skip water points
        IF((XLAND(I,J)-1.5) >= 0.) CYCLE

        ! Debug: Print first few grid points to understand coordinate values
        IF (I <= 2 .AND. J <= 2) THEN
            WRITE(*,'(A,I10,I10)') 'SNOWPACK-DEBUG: Processing grid point I,J = ', I, J
        ENDIF

        cell_active = .TRUE.
        IF (snpack_mode_trim == 'antarctica') THEN
            IF (ISLTYP(I,J) /= 8) cell_active = .FALSE.
        ELSEIF (snpack_mode_trim == 'land') THEN
            IF (ISLTYP(I,J) == 14) cell_active = .FALSE.
        END IF

        IF (.NOT. cell_active) THEN
            acc_samples(I,J) = 0
            acc_ta(I,J) = 0.0
            acc_qv(I,J) = 0.0
            acc_u(I,J)  = 0.0
            acc_v(I,J)  = 0.0
            acc_sw(I,J) = 0.0
            acc_lw(I,J) = 0.0
            acc_p(I,J)  = 0.0
            acc_rain(I,J) = 0.0
            CYCLE
        END IF

        ! Accumulate forcing samples before handing them to Snowpack.
        ! Mirrors the buffering logic in the legacy coupler (CRYOWRF/src/coupler/main_coupler/Coupler.cpp:820-906).
        acc_samples(I,J) = acc_samples(I,J) + 1
        acc_ta(I,J)   = acc_ta(I,J)   + T_PHY(I,1,J)
        acc_qv(I,J)   = acc_qv(I,J)   + QV_CURR(I,1,J)
        acc_u(I,J)    = acc_u(I,J)    + U_PHY(I,1,J)
        acc_v(I,J)    = acc_v(I,J)    + V_PHY(I,1,J)
        acc_sw(I,J)   = acc_sw(I,J)   + SWDOWN(I,J)
        acc_lw(I,J)   = acc_lw(I,J)   + GLW(I,J)
        acc_p(I,J)    = acc_p(I,J)    + P8W(I,1,J)
        acc_rain(I,J) = acc_rain(I,J) + RAINBL(I,J)

        IF (acc_samples(I,J) < snowpack_compute_counter) CYCLE

        ! Compute averaged forcings for Snowpack call
        ta_avg_val = acc_ta(I,J)   / REAL(acc_samples(I,J))
        qv_avg_val = acc_qv(I,J)   / REAL(acc_samples(I,J))
        u_avg_val  = acc_u(I,J)    / REAL(acc_samples(I,J))
        v_avg_val  = acc_v(I,J)    / REAL(acc_samples(I,J))
        sw_avg_val = acc_sw(I,J)   / REAL(acc_samples(I,J))
        lw_avg_val = acc_lw(I,J)   / REAL(acc_samples(I,J))
        p_avg_val  = acc_p(I,J)    / REAL(acc_samples(I,J))
        rain_total = acc_rain(I,J)

        dt_effective = dt * REAL(acc_samples(I,J))   ! aggregate WRF steps into the Snowpack physics step
        ! CRYOWRF (src/wrf/phys/module_sf_snowpacklsm.F:527) samples meteo at half the first layer thickness.
        ! Reinstate the 0.5*DZ8W scaling so Snowpack sees the same effective measurement height as the
        ! reference implementation; keep the 1 m floor to avoid degeneracy over shallow domains.
        met_height = MAX(REAL(snowpack_get_measurement_height_floor(), KIND=KIND(met_height)), &
                          REAL(snowpack_get_measurement_height_scale(), KIND=KIND(met_height)) * DZ8W(I,1,J))

        ! Reset accumulators for next cycle
        acc_ta(I,J) = 0.0
        acc_qv(I,J) = 0.0
        acc_u(I,J)  = 0.0
        acc_v(I,J)  = 0.0
        acc_sw(I,J) = 0.0
        acc_lw(I,J) = 0.0
        acc_p(I,J)  = 0.0
        acc_rain(I,J) = 0.0
        acc_samples(I,J) = 0

        ! Initialize interface structure
        CALL initialize_snowpack_interface(interface_data, itimestep, dt, dx)

        ! Pack WRF data into structured interface
        CALL pack_wrf_to_interface(interface_data, &
            ! Atmospheric inputs
            ta_avg_val, qv_avg_val, u_avg_val, v_avg_val, p_avg_val, &
            sw_avg_val, lw_avg_val, rain_total, &
            ! Grid information  
            XLAT(I,J), XLONG(I,J), met_height, &
            ! Model configuration
            dt_effective, ITIMESTEP, &
            ! Surface state (input/output)
            TSK(I,J), HFX(I,J), QFX(I,J), LH(I,J), ALBEDO(I,J), SNOWC(I,J), &
            SNOW(I,J), SNOWH(I,J), QSFC(I,J), Z0(I,J), ZNT(I,J), &
            ! Grid indices
            I, J, 1)

        IF (itimestep == 1 .AND. I == its .AND. J == jts) THEN
            WRITE(*,'(A,I4,A,2I6,A,F8.2,A,F8.2,A,F7.3,A,F7.3)') &
     &             'SNOWPACK-PACK tid=', thread_id, ' ij=', I, J, &
     &             ' lat=', interface_data%atmos%latitude, &
     &             ' lon=', interface_data%atmos%longitude, &
     &             ' T=', interface_data%atmos%temperature, &
     &             ' RH=', interface_data%atmos%humidity
            WRITE(*,'(A,F8.2,A,F8.2,A,F8.3)') 'SNOWPACK-PACK wind_u=', &
     &             interface_data%atmos%wind_u, ' wind_v=', interface_data%atmos%wind_v, &
     &             ' precip=', interface_data%atmos%precipitation
        END IF

        ! Call SNOWPACK physics with structured interface (serialize core)
        !$OMP CRITICAL(snowpack_bridge)
        CALL snowpack_physics_structured(interface_data)

        IF (itimestep == 1 .AND. I == its .AND. J == jts) THEN
            WRITE(*,'(A,I4,A,2I6,A,F8.2,A,F8.2,A,F8.2)') &
     &           'SNOWPACK-RETURN tid=', thread_id, ' ij=', I, J, &
     &           ' Tsurf=', interface_data%surface%surface_temp, &
     &           ' SWE=', interface_data%surface%snow_swe, &
     &           ' depth=', interface_data%surface%snow_depth
            WRITE(*,'(A,F8.2,A,F8.2)') 'SNOWPACK-RETURN HFX=', &
     &           interface_data%surface%heat_sensible, ' QFX=', interface_data%surface%heat_latent_kg
        END IF

        ! Unpack results back to WRF variables
        CALL unpack_interface_to_wrf(interface_data, &
            TSK(I,J), HFX(I,J), QFX(I,J), LH(I,J), ALBEDO(I,J), SNOWC(I,J), &
            SNOW(I,J), SNOWH(I,J), QSFC(I,J), Z0(I,J), ZNT(I,J))

        ! Apply dynamic roughness calculation following CRYOWRF approach
        ! CRYOWRF: module_sf_snowpacklsm.F:742-750
        IF (SNOWH(I,J) > 0.0) THEN
            ! Snow-covered surface: use snow roughness length
            Z0(I,J) = REAL(snowpack_get_snow_roughness_length(), KIND=KIND(Z0))
            ZNT(I,J) = Z0(I,J)
        ELSE
            ! Snow-free surface: use underlying surface roughness (preserve WRF value)
            ! Z0 and ZNT already contain WRF-calculated values for bare ground
        END IF
        
        ! Extract CRYOWRF layer arrays if present (matching CRYOWRF behavior)
        IF (PRESENT(SN_T)) THEN
            IF (should_write_layers) THEN
                CALL extract_snowpack_layers_c(I, J, layer_temp_local, layer_thick_local, layer_voli_local, &
                    layer_volw_local, layer_volv_local, layer_rg_local, layer_rb_local, layer_dd_local, &
                    layer_sp_local, layer_cdot_local, layer_meta_local, layer_depd_local, &
                    layer_graintype_local, layer_mk_local, n_layers_local)

                IF (PRESENT(SN_NLAYERS)) SN_NLAYERS(I,J) = n_layers_local

                copy_count = MIN(layer_limit, n_layers_local)

                IF (copy_count > 0) THEN
                    SN_T(I,1:copy_count,J)     = layer_temp_local(1:copy_count)
                    IF (PRESENT(SN_THICK)) SN_THICK(I,1:copy_count,J) = layer_thick_local(1:copy_count)
                    IF (PRESENT(SN_VOLI))  SN_VOLI(I,1:copy_count,J)  = layer_voli_local(1:copy_count)
                    IF (PRESENT(SN_VOLW))  SN_VOLW(I,1:copy_count,J)  = layer_volw_local(1:copy_count)
                    IF (PRESENT(SN_VOLV))  SN_VOLV(I,1:copy_count,J)  = layer_volv_local(1:copy_count)
                    IF (PRESENT(SN_RG))    SN_RG(I,1:copy_count,J)    = layer_rg_local(1:copy_count)
                    IF (PRESENT(SN_RB))    SN_RB(I,1:copy_count,J)    = layer_rb_local(1:copy_count)
                    IF (PRESENT(SN_DD))    SN_DD(I,1:copy_count,J)    = layer_dd_local(1:copy_count)
                    IF (PRESENT(SN_SP))    SN_SP(I,1:copy_count,J)    = layer_sp_local(1:copy_count)
                    IF (PRESENT(SN_CDOT))  SN_CDOT(I,1:copy_count,J)  = layer_cdot_local(1:copy_count)
                    IF (PRESENT(SN_META))  SN_META(I,1:copy_count,J)  = layer_meta_local(1:copy_count)
                    IF (PRESENT(SN_DEPD))  SN_DEPD(I,1:copy_count,J)  = layer_depd_local(1:copy_count)
                    IF (PRESENT(SN_GRAINTYPE)) SN_GRAINTYPE(I,1:copy_count,J) = layer_graintype_local(1:copy_count)
                    IF (PRESENT(SN_MK))    SN_MK(I,1:copy_count,J)    = layer_mk_local(1:copy_count)
                ELSE
                    SN_T(I,1:layer_limit,J) = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_THICK)) SN_THICK(I,1:layer_limit,J) = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_VOLI))  SN_VOLI(I,1:layer_limit,J)  = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_VOLW))  SN_VOLW(I,1:layer_limit,J)  = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_VOLV))  SN_VOLV(I,1:layer_limit,J)  = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_RG))    SN_RG(I,1:layer_limit,J)    = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_RB))    SN_RB(I,1:layer_limit,J)    = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_DD))    SN_DD(I,1:layer_limit,J)    = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_SP))    SN_SP(I,1:layer_limit,J)    = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_CDOT))  SN_CDOT(I,1:layer_limit,J)  = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_META))  SN_META(I,1:layer_limit,J)  = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_DEPD))  SN_DEPD(I,1:layer_limit,J)  = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_GRAINTYPE)) SN_GRAINTYPE(I,1:layer_limit,J) = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_MK))    SN_MK(I,1:layer_limit,J)    = SN_LAYER_SENTINEL
                END IF

                IF (copy_count < layer_limit) THEN
                    SN_T(I,copy_count+1:layer_limit,J) = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_THICK)) SN_THICK(I,copy_count+1:layer_limit,J) = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_VOLI))  SN_VOLI(I,copy_count+1:layer_limit,J)  = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_VOLW))  SN_VOLW(I,copy_count+1:layer_limit,J)  = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_VOLV))  SN_VOLV(I,copy_count+1:layer_limit,J)  = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_RG))    SN_RG(I,copy_count+1:layer_limit,J)    = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_RB))    SN_RB(I,copy_count+1:layer_limit,J)    = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_DD))    SN_DD(I,copy_count+1:layer_limit,J)    = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_SP))    SN_SP(I,copy_count+1:layer_limit,J)    = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_CDOT))  SN_CDOT(I,copy_count+1:layer_limit,J)  = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_META))  SN_META(I,copy_count+1:layer_limit,J)  = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_DEPD))  SN_DEPD(I,copy_count+1:layer_limit,J)  = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_GRAINTYPE)) SN_GRAINTYPE(I,copy_count+1:layer_limit,J) = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_MK))    SN_MK(I,copy_count+1:layer_limit,J)    = SN_LAYER_SENTINEL
                END IF

                IF (layer_limit < SIZE(SN_T,2)) THEN
                    SN_T(I,layer_limit+1:SIZE(SN_T,2),J) = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_THICK)) SN_THICK(I,layer_limit+1:SIZE(SN_THICK,2),J) = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_VOLI))  SN_VOLI(I,layer_limit+1:SIZE(SN_VOLI,2),J)   = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_VOLW))  SN_VOLW(I,layer_limit+1:SIZE(SN_VOLW,2),J)   = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_VOLV))  SN_VOLV(I,layer_limit+1:SIZE(SN_VOLV,2),J)   = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_RG))    SN_RG(I,layer_limit+1:SIZE(SN_RG,2),J)       = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_RB))    SN_RB(I,layer_limit+1:SIZE(SN_RB,2),J)       = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_DD))    SN_DD(I,layer_limit+1:SIZE(SN_DD,2),J)       = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_SP))    SN_SP(I,layer_limit+1:SIZE(SN_SP,2),J)       = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_CDOT))  SN_CDOT(I,layer_limit+1:SIZE(SN_CDOT,2),J)   = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_META))  SN_META(I,layer_limit+1:SIZE(SN_META,2),J)   = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_DEPD))  SN_DEPD(I,layer_limit+1:SIZE(SN_DEPD,2),J)   = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_GRAINTYPE)) SN_GRAINTYPE(I,layer_limit+1:SIZE(SN_GRAINTYPE,2),J) = SN_LAYER_SENTINEL
                    IF (PRESENT(SN_MK))    SN_MK(I,layer_limit+1:SIZE(SN_MK,2),J)       = SN_LAYER_SENTINEL
                END IF
            END IF
        ELSEIF (PRESENT(SN_NLAYERS)) THEN
            IF (should_write_layers) SN_NLAYERS(I,J) = interface_data%layers%n_layers
        END IF
        
        ! Extract CRYOWRF budget arrays if present (matching CRYOWRF behavior)
        IF (PRESENT(M_BUDG_PRECIP))    M_BUDG_PRECIP(I,J)   = M_BUDG_PRECIP(I,J)   + interface_data%budgets%mass_precip
        IF (PRESENT(M_BUDG_SUBLIM))    M_BUDG_SUBLIM(I,J)   = M_BUDG_SUBLIM(I,J)   + interface_data%budgets%mass_sublim
        IF (PRESENT(M_BUDG_MELT))      M_BUDG_MELT(I,J)     = M_BUDG_MELT(I,J)     + interface_data%budgets%mass_melt
        IF (PRESENT(M_BUDG_SWE))       M_BUDG_SWE(I,J)      = M_BUDG_SWE(I,J)      + interface_data%budgets%mass_swe
        IF (PRESENT(M_BUDG_REFREEZE))  M_BUDG_REFREEZE(I,J) = M_BUDG_REFREEZE(I,J) + interface_data%budgets%mass_refreeze

        IF (PRESENT(E_BUDG_LW_IN))     E_BUDG_LW_IN(I,J)    = E_BUDG_LW_IN(I,J)    + interface_data%budgets%energy_lw_in
        IF (PRESENT(E_BUDG_LW_OUT))    E_BUDG_LW_OUT(I,J)   = E_BUDG_LW_OUT(I,J)   + interface_data%budgets%energy_lw_out
        IF (PRESENT(E_BUDG_SW_IN))     E_BUDG_SW_IN(I,J)    = E_BUDG_SW_IN(I,J)    + interface_data%budgets%energy_sw_in
        IF (PRESENT(E_BUDG_SW_OUT))    E_BUDG_SW_OUT(I,J)   = E_BUDG_SW_OUT(I,J)   + interface_data%budgets%energy_sw_out
        IF (PRESENT(E_BUDG_SENSIBLE))  E_BUDG_SENSIBLE(I,J) = E_BUDG_SENSIBLE(I,J) + interface_data%budgets%energy_sensible
        IF (PRESENT(E_BUDG_LATENT))    E_BUDG_LATENT(I,J)   = E_BUDG_LATENT(I,J)   + interface_data%budgets%energy_latent
        IF (PRESENT(E_BUDG_GFLUX))     E_BUDG_GFLUX(I,J)    = E_BUDG_GFLUX(I,J)    + interface_data%budgets%energy_ground
        IF (PRESENT(E_BUDG_RAINE))     E_BUDG_RAINE(I,J)    = E_BUDG_RAINE(I,J)    + interface_data%budgets%energy_rain
        IF (PRESENT(E_BUDG_DU))        E_BUDG_DU(I,J)       = E_BUDG_DU(I,J)       + interface_data%budgets%energy_total

        !$OMP END CRITICAL(snowpack_bridge)

        IF (save_state_this_step) CALL snowpack_save_state_c(I, J)

        ! Additional CRYOWRF budget components (extracted but not passed to WRF standard interface)
        ! These are available for diagnostics and future CRYOWRF compatibility:
        ! - interface_data%budgets%mass_erosion, interface_data%budgets%mass_deposition
        ! - interface_data%budgets%energy_sw_in, interface_data%budgets%energy_sw_out
        ! - interface_data%budgets%energy_lw_out, interface_data%budgets%energy_precip
        
        ! Update accumulated variables from budget tracking
        ACSNOM(I,J) = interface_data%budgets%accum_snowmelt
        ACSNOW(I,J) = interface_data%budgets%accum_snowfall
        
        ! Set standard LSM outputs with proper SNOWPACK extraction
        ! CRYOWRF SOURCE: Energy balance and soil moisture from CRYOWRF land surface scheme
        GRDFLX(I,J) = interface_data%budgets%energy_ground  ! Ground heat flux from SNOWPACK [W/m2]
        SMSTAV(I,J) = interface_data%soil%moisture_avail   ! Soil moisture availability from SNOWPACK
        SMSTOT(I,J) = interface_data%soil%moisture_total   ! Total soil water from SNOWPACK [mm]

        ! Runoff calculations from SNOWPACK mass balance
        SFCRUNOFF(I,J) = SFCRUNOFF(I,J) + interface_data%budgets%runoff_surface  ! Surface runoff from SNOWPACK
        UDRUNOFF(I,J) = UDRUNOFF(I,J) + interface_data%budgets%runoff_subsurface ! Subsurface runoff from SNOWPACK
        CANWAT(I,J) = 0.0      ! No canopy water in SNOWPACK

        ! Set emissivity based on surface conditions following CRYOWRF approach
        ! CRYOWRF: module_sf_snowpacklsm.F:772 uses 0.98 for snow/soil
        IF (SNOWH(I,J) > 0.0) THEN
            EMISS(I,J) = REAL(snowpack_get_snow_emissivity(), KIND=KIND(EMISS))
        ELSE
            EMISS(I,J) = REAL(snowpack_get_soil_emissivity(), KIND=KIND(EMISS))
        END IF

        ! Set soil properties from SNOWPACK ground interface
        ! CRYOWRF SOURCE: Soil layer initialization from CRYOWRF ground interface module
        IF (multi_layer_soil_enabled) THEN
            layer_count_local = MIN(nsoil, interface_data%soil%layer_count)

            IF (layer_count_local > 0) THEN
                DO K = 1, layer_count_local
                    SMOIS(I,K,J) = interface_data%soil%moisture_vol_layers(K)
                    SH2O(I,K,J)  = interface_data%soil%moisture_liq_layers(K)
                    TSLB(I,K,J)  = interface_data%soil%temperature_layers(K)
                END DO
                DO K = layer_count_local + 1, nsoil
                    SMOIS(I,K,J) = SMOIS(I,layer_count_local,J)
                    SH2O(I,K,J)  = SH2O(I,layer_count_local,J)
                    TSLB(I,K,J)  = TSLB(I,layer_count_local,J)
                END DO
            ELSE
                DO K = 1, nsoil
                    SMOIS(I,K,J) = interface_data%soil%moisture_volumetric
                    SH2O(I,K,J)  = interface_data%soil%moisture_liquid
                    TSLB(I,K,J)  = interface_data%soil%temperature
                END DO
            END IF
        ELSE
            DO K = 1, nsoil
                SMOIS(I,K,J) = interface_data%soil%moisture_volumetric
                SH2O(I,K,J)  = interface_data%soil%moisture_liquid
                TSLB(I,K,J)  = interface_data%soil%temperature
            END DO
        END IF
        
    END DO
    END DO

END SUBROUTINE snowpacklsm

!-------------------------------------------------------------------
! Basic 2m temperature and humidity calculations with stability
! Following CRYOWRF approach: module_sf_snowpacklsm.F:781-782
!-------------------------------------------------------------------
SUBROUTINE calculate_2m_variables(TSK, HFX, QFX, PSFC, Z0, U10, V10, &
                                   T2, Q2, I, J)
    IMPLICIT NONE

    ! Input variables
    REAL, INTENT(IN) :: TSK          ! Surface temperature [K]
    REAL, INTENT(IN) :: HFX          ! Sensible heat flux [W/m2]
    REAL, INTENT(IN) :: QFX          ! Moisture flux [kg/m2/s]
    REAL, INTENT(IN) :: PSFC         ! Surface pressure [Pa]
    REAL, INTENT(IN) :: Z0           ! Roughness length [m]
    REAL, INTENT(IN) :: U10, V10     ! 10m wind components [m/s]
    INTEGER, INTENT(IN) :: I, J       ! Grid indices (for debugging)

    ! Output variables
    REAL, INTENT(OUT) :: T2           ! 2m temperature [K]
    REAL, INTENT(OUT) :: Q2           ! 2m humidity [kg/kg]

    ! Local variables
    REAL :: wind_10m, rho_air, cp_air
    REAL :: profile_constant_T, profile_constant_Q
    REAL :: stability_correction, log_term

    ! Constants
    REAL, PARAMETER :: von_karman = 0.4
    REAL, PARAMETER :: cp = 1004.67   ! Specific heat of air [J/kg/K]
    REAL, PARAMETER :: R_d = 287.04   ! Gas constant for dry air [J/kg/K]

    ! Calculate 10m wind speed
    wind_10m = SQRT(U10*U10 + V10*V10)

    ! Air density at surface
    rho_air = PSFC / (R_d * TSK)

    ! Avoid division by zero
    IF (wind_10m < 0.1) wind_10m = 0.1

    ! Profile constants following CRYOWRF approach
    ! CRYOWRF: module_sf_snowpacklsm.F:775-776
    profile_constant_T = HFX / (rho_air * cp * wind_10m * von_karman)
    profile_constant_Q = QFX / (rho_air * wind_10m * von_karman)

    ! Logarithmic profile term with simple stability correction
    log_term = LOG(2.0 / Z0)
    IF (log_term < 0.7) log_term = 0.7  ! Minimum value as in CRYOWRF

    ! Simplified stability correction (neutral conditions assumed)
    stability_correction = log_term

    ! Calculate 2m variables
    ! CRYOWRF: module_sf_snowpacklsm.F:781-782
    T2 = TSK - profile_constant_T * stability_correction
    Q2 = QFX / (wind_10m * von_karman)  ! Simplified calculation

    ! Ensure physical limits
    IF (Q2 < 0.0) Q2 = 0.0
    IF (Q2 > 0.05) Q2 = 0.05  ! Maximum specific humidity

    ! Debug output for first few grid points
    IF (I <= 2 .AND. J <= 2) THEN
        WRITE(*,'(A,I3,A,I3,A,F8.3,A,F8.3)') &
            "SNOWPACK-DEBUG [2m_calc]: Grid(", I, ",", J, &
            ") T2=", T2, " Q2=", Q2*1000
    END IF

END SUBROUTINE calculate_2m_variables

END MODULE module_sf_snowpackdrv
