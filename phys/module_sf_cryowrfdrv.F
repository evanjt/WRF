!WRF:MODEL_LAYER:PHYSICS
!
! SNOWPACK Land Surface Model Driver
! WRF interface driver for SNOWPACK physics integration
!
MODULE module_sf_snowpackdrv
    USE ISO_C_BINDING, ONLY: C_CHAR, C_NULL_CHAR

    ! Module-level configuration management
    LOGICAL, SAVE :: config_initialized = .false.
    CHARACTER(LEN=256), SAVE :: config_file_path = './io.ini'
    
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
            layer_rg, layer_rb, layer_dd, layer_sp, n_layers) &
            BIND(C, name='extract_snowpack_layers_c')
            USE ISO_C_BINDING
            INTEGER(C_INT), INTENT(IN) :: i_grid, j_grid
            REAL(C_FLOAT), INTENT(OUT) :: layer_temps(100), layer_thick(100)
            REAL(C_FLOAT), INTENT(OUT) :: layer_voli(100), layer_volw(100), layer_volv(100)
            REAL(C_FLOAT), INTENT(OUT) :: layer_rg(100), layer_rb(100), layer_dd(100), layer_sp(100)
            INTEGER(C_INT), INTENT(OUT) :: n_layers
        END SUBROUTINE extract_snowpack_layers_c

    END INTERFACE

CONTAINS

!-------------------------------------------------------------------
! Initialize SNOWPACK configuration (called once)
!-------------------------------------------------------------------
SUBROUTINE initialize_snowpack_module_config()
    USE module_configure  ! For nl_get_* functions
    IMPLICIT NONE
    
    ! Local variables
    CHARACTER(LEN=256) :: snowpack_config_path
    INTEGER :: domain_id
    
    IF (config_initialized) RETURN
    
    ! Get SNOWPACK configuration file path from namelist
    domain_id = 1  ! Default to domain 1 for initialization
    CALL nl_get_snowpack_config_file(domain_id, snowpack_config_path)
    
    ! Use namelist-specified path or default
    IF (LEN_TRIM(snowpack_config_path) > 0) THEN
        config_file_path = TRIM(snowpack_config_path)
    ELSE
        config_file_path = './io.ini'
    ENDIF
    
    ! Initialize SNOWPACK configuration with the file
    CALL initialize_snowpack_config_c(TRIM(config_file_path)//C_NULL_CHAR)
    
    ! Note: WRF time initialization will be done later in snowpacklsm() 
    ! when namelist is fully available (CRYOWRF pattern)
    
    config_initialized = .true.
    
    ! Diagnostic message
    CALL wrf_message('*** SNOWPACK LSM: Initialized with config file: ' // TRIM(config_file_path) // ' ***')
    
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
               SN_T,SN_THICK,SN_VOLI,SN_VOLW,SN_VOLV,SN_RG,SN_RB,SN_DD,SN_SP,SN_NLAYERS, & ! OPTIONAL CRYOWRF layer arrays
               M_BUDG_PRECIP,M_BUDG_SUBLIM,M_BUDG_MELT,               & ! OPTIONAL CRYOWRF mass budget
               E_BUDG_LW_IN,E_BUDG_SENSIBLE,E_BUDG_LATENT,            & ! OPTIONAL CRYOWRF energy budget
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
    INTEGER, DIMENSION( ims:ime,        jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_NLAYERS  ! number of layers

! CRYOWRF mass budget arrays (OPTIONAL - matching CRYOWRF behavior)
    REAL,    DIMENSION( ims:ime,        jms:jme ), INTENT(INOUT), OPTIONAL ::  M_BUDG_PRECIP   ! mass budget precip [kg m-2]
    REAL,    DIMENSION( ims:ime,        jms:jme ), INTENT(INOUT), OPTIONAL ::  M_BUDG_SUBLIM   ! mass budget sublim [kg m-2]
    REAL,    DIMENSION( ims:ime,        jms:jme ), INTENT(INOUT), OPTIONAL ::  M_BUDG_MELT     ! mass budget melt [kg m-2]

! CRYOWRF energy budget arrays (OPTIONAL - matching CRYOWRF behavior)
    REAL,    DIMENSION( ims:ime,        jms:jme ), INTENT(INOUT), OPTIONAL ::  E_BUDG_LW_IN    ! energy budget LW in [W m-2]
    REAL,    DIMENSION( ims:ime,        jms:jme ), INTENT(INOUT), OPTIONAL ::  E_BUDG_SENSIBLE ! energy budget sensible [W m-2]
    REAL,    DIMENSION( ims:ime,        jms:jme ), INTENT(INOUT), OPTIONAL ::  E_BUDG_LATENT   ! energy budget latent [W m-2]

    INTEGER, INTENT(IN   )   ::     ids,ide, jds,jde, kds,kde
    INTEGER, INTENT(IN   )   ::     ims,ime, jms,jme, kms,kme
    INTEGER, INTENT(IN   )   ::     its,ite, jts,jte, kts,kte

! Local Variables
    INTEGER :: I, J, K
    TYPE(SnowpackInterface) :: interface_data
    
    ! Time initialization variables (CRYOWRF pattern)
    INTEGER(C_INT), SAVE :: simulation_start_year, simulation_start_month, simulation_start_day
    INTEGER(C_INT), SAVE :: simulation_start_hour, simulation_start_minute
    LOGICAL, SAVE :: time_initialized = .FALSE.
    INTEGER :: domain_id
    
    ! SNOWPACK namelist parameters
    LOGICAL :: sn_start_from_file
    REAL :: snpack_dt, snpack_write_dt
    CHARACTER(LEN=256) :: snpack_mode

!-------------------------------------------------------------------
!  Call SNOWPACK physics (simplified interface without CRYOWRF arrays)
!-------------------------------------------------------------------

    ! Initialize SNOWPACK configuration on first call
    CALL initialize_snowpack_module_config()
    
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
    IF (itimestep == 1) THEN
        CALL wrf_message('*** SNOWPACK LSM: Simplified interface active ***')
        WRITE(*,'(A,L1)') 'SNOWPACK-INFO: sn_start_from_file = ', sn_start_from_file
        WRITE(*,'(A,F8.1)') 'SNOWPACK-INFO: snpack_dt = ', snpack_dt
        WRITE(*,'(A,F8.1)') 'SNOWPACK-INFO: snpack_write_dt = ', snpack_write_dt
        WRITE(*,'(A,A)') 'SNOWPACK-INFO: snpack_mode = ', TRIM(snpack_mode)
    END IF

    ! Debug: Print loop bounds to identify coordinate corruption
    WRITE(*,'(A,I10,I10,I10,I10)') 'SNOWPACK-DEBUG: Loop bounds its,ite,jts,jte = ', its, ite, jts, jte

    ! Safety check: Verify loop bounds are reasonable
    IF (its < 0 .OR. ite < 0 .OR. jts < 0 .OR. jte < 0) THEN
        WRITE(*,'(A,I10,I10,I10,I10)') 'SNOWPACK-ERROR: Negative loop bounds detected! its,ite,jts,jte = ', its, ite, jts, jte
        RETURN
    ENDIF

    IF (its > ite .OR. jts > jte) THEN
        WRITE(*,'(A,I10,I10,I10,I10)') 'SNOWPACK-ERROR: Invalid loop bounds (start > end)! its,ite,jts,jte = ', its, ite, jts, jte
        RETURN
    ENDIF

    IF (ite > 10000 .OR. jte > 10000) THEN
        WRITE(*,'(A,I10,I10,I10,I10)') 'SNOWPACK-ERROR: Loop bounds too large! its,ite,jts,jte = ', its, ite, jts, jte
        RETURN
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
        
        ! Initialize interface structure
        CALL initialize_snowpack_interface(interface_data, itimestep, dt, dx)
        
        ! Pack WRF data into structured interface
        CALL pack_wrf_to_interface(interface_data, &
            ! Atmospheric inputs
            T_PHY(I,1,J), QV_CURR(I,1,J), U_PHY(I,1,J), V_PHY(I,1,J), P8W(I,1,J), &
            SWDOWN(I,J), GLW(I,J), RAINBL(I,J), &
            ! Grid information  
            XLAT(I,J), XLONG(I,J), DZ8W(I,1,J), &
            ! Model configuration
            DT, ITIMESTEP, &
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

        ! Apply dynamic roughness calculation following CRYOWRF approach
        ! CRYOWRF: module_sf_snowpacklsm.F:742-750
        IF (SNOWH(I,J) > 0.0) THEN
            ! Snow-covered surface: use snow roughness length
            ! TODO: Add mode selection (antarctica vs land) for different roughness values
            Z0(I,J) = 0.002   ! Snow roughness length [m]
            ZNT(I,J) = 0.002   ! Thermal roughness equals momentum roughness
        ELSE
            ! Snow-free surface: use underlying surface roughness (preserve WRF value)
            ! Z0 and ZNT already contain WRF-calculated values for bare ground
        END IF
        
        ! Extract CRYOWRF layer arrays if present (matching CRYOWRF behavior)
        IF (PRESENT(SN_T)) THEN
            CALL extract_snowpack_layers_c(I, J, &
                SN_T(I,:,J), SN_THICK(I,:,J), SN_VOLI(I,:,J), SN_VOLW(I,:,J), SN_VOLV(I,:,J), &
                SN_RG(I,:,J), SN_RB(I,:,J), SN_DD(I,:,J), SN_SP(I,:,J), SN_NLAYERS(I,J))
        END IF
        
        ! Extract CRYOWRF budget arrays if present (matching CRYOWRF behavior)
        IF (PRESENT(M_BUDG_PRECIP)) THEN
            M_BUDG_PRECIP(I,J) = interface_data%budgets%mass_precip     ! Precipitation mass [kg/m²]
            M_BUDG_SUBLIM(I,J) = interface_data%budgets%mass_sublim    ! Sublimation mass [kg/m²]
            M_BUDG_MELT(I,J) = interface_data%budgets%mass_melt        ! Melt mass [kg/m²]
        END IF

        IF (PRESENT(E_BUDG_LW_IN)) THEN
            E_BUDG_LW_IN(I,J) = interface_data%budgets%energy_lw_in     ! Longwave radiation in [J/m²]
            E_BUDG_SENSIBLE(I,J) = interface_data%budgets%energy_sensible ! Sensible heat [J/m²]
            E_BUDG_LATENT(I,J) = interface_data%budgets%energy_latent    ! Latent heat [J/m²]
        END IF

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
            EMISS(I,J) = 0.98  ! Snow emissivity
        ELSE
            EMISS(I,J) = 0.96  ! Bare soil emissivity (typical value)
        END IF

        ! Set soil properties from SNOWPACK ground interface
        ! CRYOWRF SOURCE: Soil layer initialization from CRYOWRF ground interface module
        SMOIS(I,1,J) = interface_data%soil%moisture_volumetric  ! Volumetric soil moisture from SNOWPACK
        SH2O(I,1,J) = interface_data%soil%moisture_liquid      ! Liquid soil moisture from SNOWPACK
        TSLB(I,1,J) = interface_data%soil%temperature         ! Soil temperature from SNOWPACK
        
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