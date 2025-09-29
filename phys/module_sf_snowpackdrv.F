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
            INTEGER(C_INT), INTENT(IN) :: start_year, start_month, start_day, start_hour, start_minute
        END SUBROUTINE initialize_wrf_simulation_time_c
        
        SUBROUTINE extract_snowpack_layers_c(i_grid, j_grid, &
            layer_temps, layer_thick, layer_voli, layer_volw, layer_volv, &
            layer_rg, layer_rb, layer_dd, layer_sp, n_layers) &
            BIND(C, name='extract_snowpack_layers_c')
            USE ISO_C_BINDING
            INTEGER(C_INT), INTENT(IN) :: i_grid, j_grid
            REAL(C_FLOAT), INTENT(OUT) :: layer_temps(50), layer_thick(50)
            REAL(C_FLOAT), INTENT(OUT) :: layer_voli(50), layer_volw(50), layer_volv(50)
            REAL(C_FLOAT), INTENT(OUT) :: layer_rg(50), layer_rb(50), layer_dd(50), layer_sp(50)
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
    
    ! Local variables for WRF namelist time values (CRYOWRF pattern)
    INTEGER :: simulation_start_year, simulation_start_month, simulation_start_day
    INTEGER :: simulation_start_hour, simulation_start_minute
    INTEGER :: domain_id
    
    IF (config_initialized) RETURN
    
    ! Use default io.ini path
    config_file_path = './io.ini'
    
    ! Initialize SNOWPACK configuration with the file
    CALL initialize_snowpack_config_c(TRIM(config_file_path)//C_NULL_CHAR)
    
    ! Note: WRF time initialization will be done later in snowpacklsm() 
    ! when namelist is fully available (CRYOWRF pattern)
    
    config_initialized = .true.
    
    ! Diagnostic message
    CALL wrf_message('*** SNOWPACK LSM: Initialized with WRF namelist time (CRYOWRF pattern) ***')
    
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
    REAL,    DIMENSION( ims:ime, 1:50, jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_T        ! layer temperature [K]
    REAL,    DIMENSION( ims:ime, 1:50, jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_THICK    ! layer thickness [m]
    REAL,    DIMENSION( ims:ime, 1:50, jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_VOLI     ! ice volume fraction [%]
    REAL,    DIMENSION( ims:ime, 1:50, jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_VOLW     ! water volume fraction [%]
    REAL,    DIMENSION( ims:ime, 1:50, jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_VOLV     ! air volume fraction [%]
    REAL,    DIMENSION( ims:ime, 1:50, jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_RG       ! grain radius [mm]
    REAL,    DIMENSION( ims:ime, 1:50, jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_RB       ! bond radius [mm]
    REAL,    DIMENSION( ims:ime, 1:50, jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_DD       ! dendricity [-]
    REAL,    DIMENSION( ims:ime, 1:50, jms:jme ), INTENT(INOUT), OPTIONAL ::  SN_SP       ! sphericity [-]
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
    INTEGER, SAVE :: simulation_start_year, simulation_start_month, simulation_start_day
    INTEGER, SAVE :: simulation_start_hour, simulation_start_minute
    LOGICAL, SAVE :: time_initialized = .FALSE.
    INTEGER :: domain_id

!-------------------------------------------------------------------
!  Call SNOWPACK physics (simplified interface without CRYOWRF arrays)
!-------------------------------------------------------------------

    ! Initialize SNOWPACK configuration on first call
    CALL initialize_snowpack_module_config()
    
    ! Initialize WRF simulation time on first call (CRYOWRF pattern)
    IF (.NOT. time_initialized) THEN
        ! Use proper domain ID from WRF (CRYOWRF pattern)
        domain_id = INEST  ! Use the actual domain ID passed by WRF
        
        WRITE(*,'(A,I3,A,I6)') 'SNOWPACK-DEBUG: Starting time initialization for domain ', &
                               domain_id, ' at timestep ', itimestep
        
        ! Initialize variables to detect if nl_get functions fail
        simulation_start_year   = -999
        simulation_start_month  = -999
        simulation_start_day    = -999
        simulation_start_hour   = -999
        simulation_start_minute = -999
        
        WRITE(*,'(A,5I12)') 'SNOWPACK-DEBUG: Initialized values (should be -999): ', &
                            simulation_start_year, simulation_start_month, simulation_start_day, &
                            simulation_start_hour, simulation_start_minute
        
        ! Get actual namelist values using WRF's nl_get_* functions (like CRYOWRF)
        WRITE(*,'(A)') 'SNOWPACK-DEBUG: About to call nl_get_start_year...'
        CALL nl_get_start_year   ( domain_id, simulation_start_year   )
        WRITE(*,'(A,I12)') 'SNOWPACK-DEBUG: After nl_get_start_year, value = ', simulation_start_year
        
        WRITE(*,'(A)') 'SNOWPACK-DEBUG: About to call nl_get_start_month...'
        CALL nl_get_start_month  ( domain_id, simulation_start_month  )
        WRITE(*,'(A,I12)') 'SNOWPACK-DEBUG: After nl_get_start_month, value = ', simulation_start_month
        
        WRITE(*,'(A)') 'SNOWPACK-DEBUG: About to call nl_get_start_day...'
        CALL nl_get_start_day    ( domain_id, simulation_start_day    )
        WRITE(*,'(A,I12)') 'SNOWPACK-DEBUG: After nl_get_start_day, value = ', simulation_start_day
        
        WRITE(*,'(A)') 'SNOWPACK-DEBUG: About to call nl_get_start_hour...'
        CALL nl_get_start_hour   ( domain_id, simulation_start_hour   )
        WRITE(*,'(A,I12)') 'SNOWPACK-DEBUG: After nl_get_start_hour, value = ', simulation_start_hour
        
        WRITE(*,'(A)') 'SNOWPACK-DEBUG: About to call nl_get_start_minute...'
        CALL nl_get_start_minute ( domain_id, simulation_start_minute )
        WRITE(*,'(A,I12)') 'SNOWPACK-DEBUG: After nl_get_start_minute, value = ', simulation_start_minute
        
        WRITE(*,'(A,I3,A,5I12)') 'SNOWPACK-DEBUG: Domain ', domain_id, &
                                ' final retrieved WRF namelist time: ', &
                                simulation_start_year, simulation_start_month, simulation_start_day, &
                                simulation_start_hour, simulation_start_minute
        
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
        
        ! Extract CRYOWRF layer arrays if present (matching CRYOWRF behavior)
        IF (PRESENT(SN_T)) THEN
            CALL extract_snowpack_layers_c(I, J, &
                SN_T(I,:,J), SN_THICK(I,:,J), SN_VOLI(I,:,J), SN_VOLW(I,:,J), SN_VOLV(I,:,J), &
                SN_RG(I,:,J), SN_RB(I,:,J), SN_DD(I,:,J), SN_SP(I,:,J), SN_NLAYERS(I,J))
        END IF
        
        ! Extract CRYOWRF budget arrays if present (matching CRYOWRF behavior)
        IF (PRESENT(M_BUDG_PRECIP)) THEN
            M_BUDG_PRECIP(I,J) = interface_data%budgets%mass_precip
            M_BUDG_SUBLIM(I,J) = interface_data%budgets%mass_sublim
            M_BUDG_MELT(I,J) = interface_data%budgets%mass_melt
        END IF
        
        IF (PRESENT(E_BUDG_LW_IN)) THEN
            E_BUDG_LW_IN(I,J) = interface_data%budgets%energy_lw_in
            E_BUDG_SENSIBLE(I,J) = interface_data%budgets%energy_sensible
            E_BUDG_LATENT(I,J) = interface_data%budgets%energy_latent
        END IF
        
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