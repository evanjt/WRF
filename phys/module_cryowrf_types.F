!WRF:MODEL_LAYER:PHYSICS
!
! SNOWPACK Data Structure Definitions
!
MODULE module_snowpack_types

USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_PTR, C_NULL_PTR, C_LOC
IMPLICIT NONE

INTERFACE
    FUNCTION snowpack_get_snow_roughness_length() BIND(C, name='snowpack_get_snow_roughness_length')
        USE ISO_C_BINDING
        REAL(C_DOUBLE) :: snowpack_get_snow_roughness_length
    END FUNCTION
    FUNCTION snowpack_get_default_soil_temperature() BIND(C, name='snowpack_get_default_soil_temperature')
        USE ISO_C_BINDING
        REAL(C_DOUBLE) :: snowpack_get_default_soil_temperature
    END FUNCTION
    FUNCTION snowpack_get_default_soil_moisture_avail() BIND(C, name='snowpack_get_default_soil_moisture_avail')
        USE ISO_C_BINDING
        REAL(C_DOUBLE) :: snowpack_get_default_soil_moisture_avail
    END FUNCTION
    FUNCTION snowpack_get_default_soil_moisture_total() BIND(C, name='snowpack_get_default_soil_moisture_total')
        USE ISO_C_BINDING
        REAL(C_DOUBLE) :: snowpack_get_default_soil_moisture_total
    END FUNCTION
    FUNCTION snowpack_get_default_soil_moisture_vol() BIND(C, name='snowpack_get_default_soil_moisture_vol')
        USE ISO_C_BINDING
        REAL(C_DOUBLE) :: snowpack_get_default_soil_moisture_vol
    END FUNCTION
    FUNCTION snowpack_get_default_soil_moisture_liq() BIND(C, name='snowpack_get_default_soil_moisture_liq')
        USE ISO_C_BINDING
        REAL(C_DOUBLE) :: snowpack_get_default_soil_moisture_liq
    END FUNCTION
    FUNCTION snowpack_get_default_soil_density() BIND(C, name='snowpack_get_default_soil_density')
        USE ISO_C_BINDING
        REAL(C_DOUBLE) :: snowpack_get_default_soil_density
    END FUNCTION
    FUNCTION snowpack_get_default_soil_conductivity() BIND(C, name='snowpack_get_default_soil_conductivity')
        USE ISO_C_BINDING
        REAL(C_DOUBLE) :: snowpack_get_default_soil_conductivity
    END FUNCTION
    FUNCTION snowpack_get_default_soil_heat_capacity() BIND(C, name='snowpack_get_default_soil_heat_capacity')
        USE ISO_C_BINDING
        REAL(C_DOUBLE) :: snowpack_get_default_soil_heat_capacity
    END FUNCTION
END INTERFACE

! Maximum number of snow layers for arrays
INTEGER, PARAMETER :: MAX_SNOW_LAYERS = 100
INTEGER, PARAMETER :: MAX_SOIL_LAYERS = 4

!-------------------------------------------------------------------
! C-interoperable structs matching C++ definitions in snowpack_bridge_structs.h
! These are used to call the C++ bridge directly via BIND(C)
!-------------------------------------------------------------------
TYPE, BIND(C) :: MeteoInput_C
    REAL(C_DOUBLE) :: temp_air
    REAL(C_DOUBLE) :: humidity
    REAL(C_DOUBLE) :: wind_speed
    REAL(C_DOUBLE) :: wind_dir
    REAL(C_DOUBLE) :: shortwave_in
    REAL(C_DOUBLE) :: longwave_in
    REAL(C_DOUBLE) :: precipitation
    REAL(C_DOUBLE) :: pressure
    REAL(C_DOUBLE) :: height
    REAL(C_DOUBLE) :: dt
    INTEGER(C_INT)  :: i_grid
    INTEGER(C_INT)  :: j_grid
    REAL(C_DOUBLE) :: wrf_lat
    REAL(C_DOUBLE) :: wrf_lon
END TYPE MeteoInput_C

TYPE, BIND(C) :: SnowpackOutput_C
    REAL(C_DOUBLE) :: snow_swe
    REAL(C_DOUBLE) :: snow_depth
    REAL(C_DOUBLE) :: surface_temp
    REAL(C_DOUBLE) :: heat_flux_sensible
    REAL(C_DOUBLE) :: heat_flux_latent
    REAL(C_DOUBLE) :: albedo
    REAL(C_DOUBLE) :: snow_coverage
    REAL(C_DOUBLE) :: friction_velocity
    REAL(C_DOUBLE) :: stability_param
    REAL(C_DOUBLE) :: roughness_mom
    REAL(C_DOUBLE) :: roughness_heat
    REAL(C_DOUBLE) :: latent_flux_kg
    REAL(C_DOUBLE) :: moisture_flux

    ! Soil properties from SNOWPACK ground interface
    ! CRYOWRF SOURCE: module_sf_snowpacklsm.F:364-366 updates soil moisture
    REAL(C_DOUBLE) :: soil_moisture_volumetric
    REAL(C_DOUBLE) :: soil_temperature
    REAL(C_DOUBLE) :: soil_density
    REAL(C_DOUBLE) :: soil_conductivity
    REAL(C_DOUBLE) :: soil_heat_capacity
    REAL(C_DOUBLE) :: soil_moisture_liquid
    REAL(C_DOUBLE) :: soil_moisture_avail
    REAL(C_DOUBLE) :: soil_moisture_total
    INTEGER(C_INT) :: soil_layer_count
    REAL(C_DOUBLE) :: soil_temp_layers(MAX_SOIL_LAYERS)
    REAL(C_DOUBLE) :: soil_moisture_vol_layers(MAX_SOIL_LAYERS)
    REAL(C_DOUBLE) :: soil_moisture_liq_layers(MAX_SOIL_LAYERS)
END TYPE SnowpackOutput_C

TYPE, BIND(C) :: SnowpackLayerData_C
    INTEGER(C_INT) :: n_layers
    REAL(C_DOUBLE) :: layer_temp(100)
    REAL(C_DOUBLE) :: layer_thick(100)
    REAL(C_DOUBLE) :: layer_vol_ice(100)
    REAL(C_DOUBLE) :: layer_vol_water(100)
    REAL(C_DOUBLE) :: layer_vol_air(100)
    REAL(C_DOUBLE) :: layer_grain_radius(100)
    REAL(C_DOUBLE) :: layer_bond_radius(100)
    REAL(C_DOUBLE) :: layer_dendricity(100)
    REAL(C_DOUBLE) :: layer_sphericity(100)
END TYPE SnowpackLayerData_C

TYPE, BIND(C) :: BudgetData_C
    REAL(C_DOUBLE) :: mass_precip
    REAL(C_DOUBLE) :: mass_sublim
    REAL(C_DOUBLE) :: mass_melt
    REAL(C_DOUBLE) :: mass_swe
    REAL(C_DOUBLE) :: mass_refreeze
    REAL(C_DOUBLE) :: energy_lw_in
    REAL(C_DOUBLE) :: energy_lw_out
    REAL(C_DOUBLE) :: energy_sw_in
    REAL(C_DOUBLE) :: energy_sw_out
    REAL(C_DOUBLE) :: energy_sensible
    REAL(C_DOUBLE) :: energy_latent
    REAL(C_DOUBLE) :: energy_ground_flux
    REAL(C_DOUBLE) :: energy_rain
    REAL(C_DOUBLE) :: energy_total
END TYPE BudgetData_C

! External C++ bridge function declaration
EXTERNAL :: snowpack_physics_c



!-------------------------------------------------------------------
! Atmospheric forcing data structure
!-------------------------------------------------------------------
TYPE :: AtmosForcing
    ! Basic meteorological variables
    REAL    :: temperature     ! Air temperature [K]
    REAL    :: humidity        ! Specific humidity [kg/kg] 
    REAL    :: wind_u          ! U-component wind [m/s]
    REAL    :: wind_v          ! V-component wind [m/s]
    REAL    :: pressure        ! Air pressure [Pa]
    
    ! Radiation and precipitation
    REAL    :: shortwave_down  ! Downward SW radiation [W/m²]
    REAL    :: longwave_down   ! Downward LW radiation [W/m²]
    REAL    :: precipitation   ! Precipitation rate [mm/s]
    
    ! Grid information
    REAL    :: height          ! Measurement height [m]
    REAL    :: latitude        ! Grid latitude [rad]
    REAL    :: longitude       ! Grid longitude [rad]
END TYPE AtmosForcing

!-------------------------------------------------------------------
! Snow layer data structure (for detailed diagnostics)
!-------------------------------------------------------------------
TYPE :: SnowLayerData
    INTEGER :: n_layers                                 ! Number of active layers
    REAL    :: layer_temp(MAX_SNOW_LAYERS)             ! Layer temperatures [K]
    REAL    :: layer_thick(MAX_SNOW_LAYERS)            ! Layer thickness [m]
    REAL    :: vol_ice(MAX_SNOW_LAYERS)                ! Ice volume fraction [%]
    REAL    :: vol_water(MAX_SNOW_LAYERS)              ! Water volume fraction [%]
    REAL    :: vol_air(MAX_SNOW_LAYERS)                ! Air volume fraction [%]
    REAL    :: grain_radius(MAX_SNOW_LAYERS)           ! Grain radius [mm]
    REAL    :: bond_radius(MAX_SNOW_LAYERS)            ! Bond radius [mm]
    REAL    :: dendricity(MAX_SNOW_LAYERS)             ! Dendricity [-]
    REAL    :: sphericity(MAX_SNOW_LAYERS)             ! Sphericity [-]
END TYPE SnowLayerData

!-------------------------------------------------------------------
! Surface state and fluxes
!-------------------------------------------------------------------
TYPE :: SurfaceState
    ! Basic snow state
    REAL    :: snow_swe        ! Snow water equivalent [mm]
    REAL    :: snow_depth      ! Physical snow depth [m]
    REAL    :: snow_coverage   ! Snow coverage fraction [0-1]
    REAL    :: surface_temp    ! Surface temperature [K]
    REAL    :: albedo          ! Surface albedo [0-1]
    REAL    :: emissivity      ! Surface emissivity [0-1]
    
    ! Surface fluxes
    REAL    :: heat_sensible   ! Sensible heat flux [W/m²]
    REAL    :: heat_latent     ! Latent heat flux [W/m²]
    REAL    :: heat_latent_kg  ! Latent heat flux [kg/m²/s]
    REAL    :: heat_ground     ! Ground heat flux [W/m²]
    REAL    :: moisture_flux   ! Surface moisture flux [kg/m²/s]
    
    ! Roughness
    REAL    :: roughness_mom   ! Momentum roughness length [m]
    REAL    :: roughness_heat  ! Heat roughness length [m]

    ! Atmospheric variables (CRYOWRF compatibility)
    REAL    :: friction_velocity ! Friction velocity [m/s]
    REAL    :: stability_param  ! Atmospheric stability parameter [dimensionless]
END TYPE SurfaceState

!-------------------------------------------------------------------
! Mass and energy budget tracking
!-------------------------------------------------------------------
TYPE :: PhysicsBudgets
    ! Mass budget components [kg/m²] - accumulated over timestep
    REAL    :: mass_precip     ! Mass from precipitation
    REAL    :: mass_sublim     ! Mass from sublimation (negative = loss)
    REAL    :: mass_evap       ! Mass from evaporation (negative = loss)
    REAL    :: mass_melt       ! Mass from melting (negative = loss)
    REAL    :: mass_swe        ! Snow water equivalent [kg/m²]
    REAL    :: mass_refreeze   ! Mass from refreezing (positive = gain)
    REAL    :: mass_erosion    ! Mass from wind erosion (negative = loss)
    REAL    :: mass_deposition ! Mass from wind deposition
    
    ! Energy budget components [J/m²] - accumulated over timestep
    REAL    :: energy_sw_in    ! Shortwave radiation input
    REAL    :: energy_sw_out   ! Shortwave radiation output (reflected)
    REAL    :: energy_lw_in    ! Longwave radiation input
    REAL    :: energy_lw_out   ! Longwave radiation output
    REAL    :: energy_sensible ! Sensible heat exchange
    REAL    :: energy_latent   ! Latent heat exchange
    REAL    :: energy_ground   ! Ground/soil heat exchange
    REAL    :: energy_precip   ! Energy from precipitation
    REAL    :: energy_rain     ! Energy from rain
    REAL    :: energy_total    ! Total energy balance
    
    ! Long-term accumulators [mm] - persistent across timesteps
    REAL    :: accum_snowfall  ! Accumulated snowfall since start
    REAL    :: accum_snowmelt  ! Accumulated snow melt since start
    REAL    :: accum_sublim    ! Accumulated sublimation since start

    ! Runoff components [mm] - accumulated over timestep (CRYOWRF compatibility)
    REAL    :: runoff_surface  ! Surface runoff from snow melt
    REAL    :: runoff_subsurface ! Subsurface runoff from soil
END TYPE PhysicsBudgets

!-------------------------------------------------------------------
! Soil state and properties
!-------------------------------------------------------------------
TYPE :: SoilState
    ! Soil moisture variables
    REAL    :: moisture_avail     ! Soil moisture availability [0-1]
    REAL    :: moisture_total     ! Total soil water [mm]
    REAL    :: moisture_volumetric ! Volumetric soil moisture [m³/m³]
    REAL    :: moisture_liquid    ! Liquid soil moisture [m³/m³]
    REAL    :: temperature        ! Soil temperature [K]
    REAL    :: density            ! Soil density [kg/m³]
    REAL    :: conductivity       ! Soil thermal conductivity [W/m/K]
    REAL    :: heat_capacity      ! Soil heat capacity [J/kg/K]
    INTEGER :: layer_count        ! Number of soil layers provided
    REAL    :: moisture_vol_layers(MAX_SOIL_LAYERS)
    REAL    :: moisture_liq_layers(MAX_SOIL_LAYERS)
    REAL    :: temperature_layers(MAX_SOIL_LAYERS)
END TYPE SoilState

!-------------------------------------------------------------------
! Model configuration and control
!-------------------------------------------------------------------
TYPE :: ModelConfig
    ! Timestep information
    INTEGER :: timestep_number ! Current timestep number
    REAL    :: dt              ! Model timestep [s]
    REAL    :: dx              ! Horizontal grid spacing [m]
    
    ! Physics options
    LOGICAL :: enable_detailed_layers  ! Extract detailed layer diagnostics
    LOGICAL :: enable_mass_budgets     ! Track mass budget components
    LOGICAL :: enable_energy_budgets   ! Track energy budget components
    LOGICAL :: enable_time_series      ! Output time series
    
    ! Grid indices for current call
    INTEGER :: i_start, i_end   ! I-dimension bounds
    INTEGER :: j_start, j_end   ! J-dimension bounds
    INTEGER :: i_grid, j_grid   ! Current grid point indices
END TYPE ModelConfig

!-------------------------------------------------------------------
! Combined interface structure for clean function calls
!-------------------------------------------------------------------
TYPE :: SnowpackInterface
    TYPE(AtmosForcing)   :: atmos      ! Atmospheric forcing
    TYPE(SurfaceState)   :: surface    ! Surface state and fluxes
    TYPE(SoilState)      :: soil       ! Soil state and properties
    TYPE(SnowLayerData)  :: layers     ! Detailed layer information
    TYPE(PhysicsBudgets) :: budgets    ! Mass/energy budget tracking
    TYPE(ModelConfig)    :: config     ! Model configuration
END TYPE SnowpackInterface

!-------------------------------------------------------------------
! Helper routines for data structure management
!-------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------
SUBROUTINE initialize_snowpack_interface(snowpack_data, itimestep, dt, dx)
!-------------------------------------------------------------------
    TYPE(SnowpackInterface), INTENT(INOUT) :: snowpack_data
    INTEGER, INTENT(IN) :: itimestep
    REAL, INTENT(IN) :: dt, dx
    
    ! Initialize configuration
    snowpack_data%config%timestep_number = itimestep
    snowpack_data%config%dt = dt
    snowpack_data%config%dx = dx
    
    ! Set default physics options (can be configured via namelist later)
    snowpack_data%config%enable_detailed_layers = .TRUE.
    snowpack_data%config%enable_mass_budgets = .TRUE.
    snowpack_data%config%enable_energy_budgets = .TRUE.
    snowpack_data%config%enable_time_series = .FALSE.
    
    ! Initialize layer data
    snowpack_data%layers%n_layers = 0
    snowpack_data%layers%layer_temp = 0.0
    snowpack_data%layers%layer_thick = 0.0

    ! Initialize budgets (important - don't reset accumulators!)
    snowpack_data%budgets%mass_precip = 0.0
    snowpack_data%budgets%mass_sublim = 0.0
    snowpack_data%budgets%energy_sw_in = 0.0
    snowpack_data%budgets%energy_lw_in = 0.0
    ! Note: accum_* variables should persist between calls

    ! Initialize atmospheric variables (CRYOWRF compatibility)
    snowpack_data%surface%friction_velocity = 0.0
    snowpack_data%surface%stability_param = 0.0
    ! Initial placeholder: overwritten after the first Snowpack call
    snowpack_data%surface%roughness_mom = REAL(snowpack_get_snow_roughness_length(), KIND=KIND(1.0))
    snowpack_data%surface%roughness_heat = snowpack_data%surface%roughness_mom
    snowpack_data%surface%heat_latent_kg = 0.0
    snowpack_data%surface%moisture_flux = 0.0

    ! Initialize soil variables with defaults for interface compatibility
    ! NOTE: These values are NOT used - real SNOWPACK soil data is handled elsewhere
    ! SNOWPACK soil data is extracted in bridge.cpp and passed through module_sf_cryowrf.F
    ! CRYOWRF SOURCE: module_sf_snowpacklsm.F:364-366 should extract SNOWPACK soil data
    ! Initial placeholders: each field is replaced with Snowpack output at runtime
    snowpack_data%soil%moisture_avail = REAL(snowpack_get_default_soil_moisture_avail(), KIND=KIND(1.0))
    snowpack_data%soil%moisture_total = REAL(snowpack_get_default_soil_moisture_total(), KIND=KIND(1.0))
    snowpack_data%soil%moisture_volumetric = REAL(snowpack_get_default_soil_moisture_vol(), KIND=KIND(1.0))
    snowpack_data%soil%moisture_liquid = REAL(snowpack_get_default_soil_moisture_liq(), KIND=KIND(1.0))
    snowpack_data%soil%temperature = REAL(snowpack_get_default_soil_temperature(), KIND=KIND(1.0))
    snowpack_data%soil%density = REAL(snowpack_get_default_soil_density(), KIND=KIND(1.0))
    snowpack_data%soil%conductivity = REAL(snowpack_get_default_soil_conductivity(), KIND=KIND(1.0))
    snowpack_data%soil%heat_capacity = REAL(snowpack_get_default_soil_heat_capacity(), KIND=KIND(1.0))
    snowpack_data%soil%layer_count = 0
    snowpack_data%soil%moisture_vol_layers = 0.0
    snowpack_data%soil%moisture_liq_layers = 0.0
    snowpack_data%soil%temperature_layers = REAL(snowpack_get_default_soil_temperature(), KIND=KIND(1.0))
    
END SUBROUTINE initialize_snowpack_interface

!-------------------------------------------------------------------
SUBROUTINE pack_wrf_to_interface(snowpack_data, &
    ! Atmospheric inputs
    t_phy, qv_curr, u_phy, v_phy, p8w, swdown, glw, rainbl, &
    ! Grid information  
    xlat, xlong, meas_height, &
    ! Model configuration
    dt, itimestep, &
    ! Surface state (input/output)
    tsk, hfx, qfx, lh, albedo, snowc, snow, snowh, qsfc, z0, znt, &
    ! Grid indices
    i, j, k)
!-------------------------------------------------------------------
    TYPE(SnowpackInterface), INTENT(INOUT) :: snowpack_data
    
    ! WRF variables
    REAL, INTENT(IN) :: t_phy, qv_curr, u_phy, v_phy, p8w
    REAL, INTENT(IN) :: swdown, glw, rainbl
    REAL, INTENT(IN) :: xlat, xlong, meas_height
    REAL, INTENT(IN) :: dt
    INTEGER, INTENT(IN) :: itimestep
    REAL, INTENT(INOUT) :: tsk, hfx, qfx, lh, albedo, snowc
    REAL, INTENT(INOUT) :: snow, snowh, qsfc, z0, znt
    INTEGER, INTENT(IN) :: i, j, k
    
    ! Pack atmospheric forcing
    ! CRITICAL: Convert QV (mixing ratio) to RH (relative humidity) following CRYOWRF approach
    ! CRYOWRF: module_sf_snowpacklsm.F lines 529, 575-578
    REAL :: q_ml, es_sat, qs_sat, rh_calc

    ! Convert mixing ratio to specific humidity: Q = QV/(1+QV)
    q_ml = qv_curr / (1.0 + qv_curr)

    ! Calculate saturation vapor pressure and specific humidity (simplified Magnus formula)
    ! CRYOWRF uses Qsat_noah function - this is equivalent for typical atmospheric conditions
    es_sat = 611.2 * EXP(17.67 * (t_phy - 273.15) / (t_phy - 29.65))  ! Saturation vapor pressure [Pa]
    qs_sat = 0.622 * es_sat / (p8w - 0.378 * es_sat)  ! Saturation specific humidity [kg/kg]

    ! Calculate relative humidity: RH = Q / Q_sat, clamped to [0,1]
    rh_calc = MIN(1.0, MAX(q_ml / qs_sat, 0.0))

    snowpack_data%atmos%temperature = t_phy
    snowpack_data%atmos%humidity = rh_calc  ! Properly calculated RH, not raw QV!
    snowpack_data%atmos%wind_u = u_phy
    snowpack_data%atmos%wind_v = v_phy
    snowpack_data%atmos%pressure = p8w
    snowpack_data%atmos%shortwave_down = swdown
    snowpack_data%atmos%longwave_down = glw
    snowpack_data%atmos%precipitation = rainbl
    snowpack_data%atmos%height = meas_height
    snowpack_data%atmos%latitude = xlat
    snowpack_data%atmos%longitude = xlong
    
    ! Debug first few grid points
    IF (i <= 2 .AND. j <= 2) THEN
        WRITE(*,'(A,I3,A,I3,A,F10.5,A,F10.5)') &
            "SNOWPACK-DEBUG [pack_wrf_to_interface]: Grid(", i, ",", j, &
            ") Packed lat=", xlat, " lon=", xlong
    ENDIF
    
    ! Pack current surface state
    snowpack_data%surface%snow_swe = snow
    snowpack_data%surface%snow_depth = snowh
    snowpack_data%surface%snow_coverage = snowc
    snowpack_data%surface%surface_temp = tsk
    snowpack_data%surface%albedo = albedo
    snowpack_data%surface%heat_sensible = hfx
    snowpack_data%surface%heat_latent = lh
    snowpack_data%surface%heat_latent_kg = qfx
    snowpack_data%surface%moisture_flux = qsfc
    snowpack_data%surface%roughness_mom = MAX(1.0E-4, z0)
    snowpack_data%surface%roughness_heat = MAX(1.0E-4, znt)
    
    ! Pack model configuration
    snowpack_data%config%dt = dt
    snowpack_data%config%timestep_number = itimestep
    snowpack_data%config%i_grid = i
    snowpack_data%config%j_grid = j
    
END SUBROUTINE pack_wrf_to_interface

!-------------------------------------------------------------------
SUBROUTINE unpack_interface_to_wrf(snowpack_data, &
    ! Surface state outputs
    tsk, hfx, qfx, lh, albedo, snowc, snow, snowh, qsfc, z0, znt)
!-------------------------------------------------------------------
    TYPE(SnowpackInterface), INTENT(IN) :: snowpack_data
    
    ! WRF output variables
    REAL, INTENT(OUT) :: tsk, hfx, qfx, lh, albedo, snowc
    REAL, INTENT(OUT) :: snow, snowh, qsfc, z0, znt
    
    ! Unpack surface state
    tsk = snowpack_data%surface%surface_temp
    hfx = snowpack_data%surface%heat_sensible
    qfx = snowpack_data%surface%heat_latent_kg
    lh = snowpack_data%surface%heat_latent
    albedo = snowpack_data%surface%albedo
    snowc = snowpack_data%surface%snow_coverage
    snow = snowpack_data%surface%snow_swe
    snowh = snowpack_data%surface%snow_depth
    qsfc = snowpack_data%surface%moisture_flux
    z0 = snowpack_data%surface%roughness_mom
    znt = snowpack_data%surface%roughness_heat

    ! Note: CRYOWRF atmospheric variables (friction_velocity, stability_param)
    ! are available in snowpack_data%surface% but not used in standard WRF LSM interface
    
END SUBROUTINE unpack_interface_to_wrf

!-------------------------------------------------------------------
SUBROUTINE snowpack_physics_structured(snowpack_data)
!-------------------------------------------------------------------
! Structured interface to SNOWPACK physics with layer extraction
!-------------------------------------------------------------------
    USE, INTRINSIC :: iso_c_binding
#ifdef _OPENMP
    USE omp_lib
#endif
    TYPE(SnowpackInterface), INTENT(INOUT) :: snowpack_data
    REAL(C_DOUBLE), PARAMETER :: RAD2DEG = 180.0_C_DOUBLE / 3.1415926535897932384626433832795_C_DOUBLE

    ! Interface for C++ bridge function with proper BIND(C) to prevent name mangling
    INTERFACE
        SUBROUTINE snowpack_physics_c(meteo_ptr, output_ptr, layer_ptr, budget_ptr) BIND(C, name='snowpack_physics_c')
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE :: meteo_ptr
            TYPE(C_PTR), VALUE :: output_ptr
            TYPE(C_PTR), VALUE :: layer_ptr
            TYPE(C_PTR), VALUE :: budget_ptr
        END SUBROUTINE snowpack_physics_c
    END INTERFACE

    
    ! Local variables for C interface (use C-interoperable derived types)
    TYPE(MeteoInput_C), TARGET :: meteo_c
    TYPE(SnowpackOutput_C), TARGET :: output_c
    TYPE(SnowpackLayerData_C), TARGET :: layers_c
    TYPE(BudgetData_C), TARGET :: budgets_c
    INTEGER :: k
    INTEGER :: thread_id
    INTEGER(C_INTPTR_T) :: meteo_addr

#ifdef _OPENMP
    thread_id = omp_get_thread_num()
#else
    thread_id = 0
#endif
    
    ! Pack data into C-interoperable MeteoInput_C
    meteo_c%temp_air = DBLE(snowpack_data%atmos%temperature)
    meteo_c%humidity = DBLE(snowpack_data%atmos%humidity)
    meteo_c%wind_speed = SQRT(DBLE(snowpack_data%atmos%wind_u)**2 + DBLE(snowpack_data%atmos%wind_v)**2)
    meteo_c%wind_dir = ATAN2(DBLE(snowpack_data%atmos%wind_v), DBLE(snowpack_data%atmos%wind_u)) * RAD2DEG
    meteo_c%shortwave_in = DBLE(snowpack_data%atmos%shortwave_down)
    meteo_c%longwave_in = DBLE(snowpack_data%atmos%longwave_down)
    meteo_c%precipitation = DBLE(snowpack_data%atmos%precipitation)
    meteo_c%pressure = DBLE(snowpack_data%atmos%pressure)
    meteo_c%height = DBLE(snowpack_data%atmos%height)
    meteo_c%dt = DBLE(snowpack_data%config%dt)
    meteo_c%i_grid = snowpack_data%config%i_grid
    meteo_c%j_grid = snowpack_data%config%j_grid
    meteo_c%wrf_lat = DBLE(snowpack_data%atmos%latitude)
    meteo_c%wrf_lon = DBLE(snowpack_data%atmos%longitude)

    ! Initialize optional outputs to safe defaults
    output_c%snow_swe = 0.0
    output_c%snow_depth = 0.0
    output_c%surface_temp = 0.0
    output_c%heat_flux_sensible = 0.0
    output_c%heat_flux_latent = 0.0
    output_c%albedo = 0.0
    output_c%snow_coverage = 0.0
    output_c%friction_velocity = 0.0
    output_c%stability_param = 0.0
    output_c%roughness_mom = 0.0
    output_c%roughness_heat = 0.0
    output_c%latent_flux_kg = 0.0
    output_c%moisture_flux = 0.0

    layers_c%n_layers = 0
    layers_c%layer_temp = 0.0_C_DOUBLE
    layers_c%layer_thick = 0.0_C_DOUBLE
    layers_c%layer_vol_ice = 0.0_C_DOUBLE
    layers_c%layer_vol_water = 0.0_C_DOUBLE
    layers_c%layer_vol_air = 0.0_C_DOUBLE
    layers_c%layer_grain_radius = 0.0_C_DOUBLE
    layers_c%layer_bond_radius = 0.0_C_DOUBLE
    layers_c%layer_dendricity = 0.0_C_DOUBLE
    layers_c%layer_sphericity = 0.0_C_DOUBLE

    budgets_c%mass_precip = 0.0_C_DOUBLE
    budgets_c%mass_sublim = 0.0_C_DOUBLE
    budgets_c%mass_melt = 0.0_C_DOUBLE
    budgets_c%mass_swe = 0.0_C_DOUBLE
    budgets_c%mass_refreeze = 0.0_C_DOUBLE
    budgets_c%energy_lw_in = 0.0_C_DOUBLE
    budgets_c%energy_lw_out = 0.0_C_DOUBLE
    budgets_c%energy_sw_in = 0.0_C_DOUBLE
    budgets_c%energy_sw_out = 0.0_C_DOUBLE
    budgets_c%energy_sensible = 0.0_C_DOUBLE
    budgets_c%energy_latent = 0.0_C_DOUBLE
    budgets_c%energy_ground_flux = 0.0_C_DOUBLE
    budgets_c%energy_rain = 0.0_C_DOUBLE
    budgets_c%energy_total = 0.0_C_DOUBLE

    IF (snowpack_data%config%timestep_number == 1 .AND. snowpack_data%config%i_grid <= 4 .AND. snowpack_data%config%j_grid <= 4) THEN
        meteo_addr = TRANSFER(C_LOC(meteo_c), 0_C_INTPTR_T)
        WRITE(*,'(A,I4,A,2I6,A,F8.2,A,F8.2,A,F7.3,A,F7.3)') &
             'SNOWPACK-STRUCT pack tid=', thread_id, ' ij=', &
             snowpack_data%config%i_grid, snowpack_data%config%j_grid, &
             ' T=', snowpack_data%atmos%temperature, &
             ' wind=', meteo_c%wind_speed, ' RH=', meteo_c%humidity, &
             ' height=', meteo_c%height
        WRITE(*,'(A,I0)') 'SNOWPACK-STRUCT meteo_ptr=', meteo_addr
    END IF

    ! Call the C++ bridge using the C-interoperable derived types
    CALL snowpack_physics_c(C_LOC(meteo_c), C_LOC(output_c), C_LOC(layers_c), C_LOC(budgets_c))

    IF (snowpack_data%config%timestep_number == 1 .AND. snowpack_data%config%i_grid <= 4 .AND. snowpack_data%config%j_grid <= 4) THEN
        WRITE(*,'(A,I4,A,2I6,A,F8.2,A,F8.2,A,F8.2)') 'SNOWPACK-STRUCT return tid=', thread_id, &
     &       ' ij=', snowpack_data%config%i_grid, snowpack_data%config%j_grid, &
     &       ' Tsurf=', output_c%surface_temp, &
     &       ' SWE=', output_c%snow_swe, &
     &       ' depth=', output_c%snow_depth
    END IF

    ! Unpack surface results back to interface structure
    snowpack_data%surface%snow_swe = REAL(output_c%snow_swe)
    snowpack_data%surface%snow_depth = REAL(output_c%snow_depth)
    snowpack_data%surface%surface_temp = REAL(output_c%surface_temp)
    snowpack_data%surface%heat_sensible = REAL(output_c%heat_flux_sensible)
    snowpack_data%surface%heat_latent = REAL(output_c%heat_flux_latent)
    snowpack_data%surface%albedo = REAL(output_c%albedo)
    snowpack_data%surface%snow_coverage = REAL(output_c%snow_coverage)
    snowpack_data%surface%friction_velocity = REAL(output_c%friction_velocity)
    snowpack_data%surface%stability_param = REAL(output_c%stability_param)
    snowpack_data%surface%roughness_mom = MAX(1.0E-4, REAL(output_c%roughness_mom))
    snowpack_data%surface%roughness_heat = MAX(1.0E-4, REAL(output_c%roughness_heat))
    snowpack_data%surface%heat_latent_kg = REAL(output_c%latent_flux_kg)
    snowpack_data%surface%moisture_flux = REAL(output_c%moisture_flux)
    IF (ABS(snowpack_data%surface%moisture_flux) < 1.0E-12) THEN
        snowpack_data%surface%moisture_flux = snowpack_data%surface%heat_latent_kg
    END IF

    ! Store extracted layer data
    snowpack_data%layers%n_layers = MIN(INT(layers_c%n_layers), MAX_SNOW_LAYERS)
    DO k = 1, snowpack_data%layers%n_layers
        snowpack_data%layers%layer_temp(k) = REAL(layers_c%layer_temp(k))
        snowpack_data%layers%layer_thick(k) = REAL(layers_c%layer_thick(k))
        snowpack_data%layers%vol_ice(k) = REAL(layers_c%layer_vol_ice(k))
        snowpack_data%layers%vol_water(k) = REAL(layers_c%layer_vol_water(k))
        snowpack_data%layers%vol_air(k) = REAL(layers_c%layer_vol_air(k))
        snowpack_data%layers%grain_radius(k) = REAL(layers_c%layer_grain_radius(k))
        snowpack_data%layers%bond_radius(k) = REAL(layers_c%layer_bond_radius(k))
        snowpack_data%layers%dendricity(k) = REAL(layers_c%layer_dendricity(k))
        snowpack_data%layers%sphericity(k) = REAL(layers_c%layer_sphericity(k))
    END DO

    ! Clear remaining layers
    DO k = snowpack_data%layers%n_layers + 1, MAX_SNOW_LAYERS
        snowpack_data%layers%layer_temp(k) = 0.0
        snowpack_data%layers%layer_thick(k) = 0.0
        snowpack_data%layers%vol_ice(k) = 0.0
        snowpack_data%layers%vol_water(k) = 0.0
        snowpack_data%layers%vol_air(k) = 0.0
        snowpack_data%layers%grain_radius(k) = 0.0
        snowpack_data%layers%bond_radius(k) = 0.0
        snowpack_data%layers%dendricity(k) = 0.0
        snowpack_data%layers%sphericity(k) = 0.0
    END DO

    ! Store budget data
    snowpack_data%budgets%mass_precip = REAL(budgets_c%mass_precip)
    snowpack_data%budgets%mass_sublim = REAL(budgets_c%mass_sublim)
    snowpack_data%budgets%mass_melt = REAL(budgets_c%mass_melt)
    snowpack_data%budgets%mass_swe = REAL(budgets_c%mass_swe)
    snowpack_data%budgets%mass_refreeze = REAL(budgets_c%mass_refreeze)
    snowpack_data%budgets%energy_lw_in = REAL(budgets_c%energy_lw_in)
    snowpack_data%budgets%energy_lw_out = REAL(budgets_c%energy_lw_out)
    snowpack_data%budgets%energy_sw_in = REAL(budgets_c%energy_sw_in)
    snowpack_data%budgets%energy_sw_out = REAL(budgets_c%energy_sw_out)
    snowpack_data%budgets%energy_sensible = REAL(budgets_c%energy_sensible)
    snowpack_data%budgets%energy_latent = REAL(budgets_c%energy_latent)
    snowpack_data%budgets%energy_ground = REAL(budgets_c%energy_ground_flux)
    snowpack_data%budgets%energy_rain = REAL(budgets_c%energy_rain)
    snowpack_data%budgets%energy_total = REAL(budgets_c%energy_total)

    ! Store soil state
    snowpack_data%soil%moisture_volumetric = REAL(output_c%soil_moisture_volumetric)
    snowpack_data%soil%temperature = REAL(output_c%soil_temperature)
    snowpack_data%soil%density = REAL(output_c%soil_density)
    snowpack_data%soil%conductivity = REAL(output_c%soil_conductivity)
    snowpack_data%soil%heat_capacity = REAL(output_c%soil_heat_capacity)
    snowpack_data%soil%moisture_liquid = REAL(output_c%soil_moisture_liquid)
    snowpack_data%soil%moisture_avail = REAL(output_c%soil_moisture_avail)
    snowpack_data%soil%moisture_total = REAL(output_c%soil_moisture_total)
    snowpack_data%soil%layer_count = MIN(output_c%soil_layer_count, MAX_SOIL_LAYERS)
    DO k = 1, snowpack_data%soil%layer_count
        snowpack_data%soil%temperature_layers(k) = REAL(output_c%soil_temp_layers(k))
        snowpack_data%soil%moisture_vol_layers(k) = REAL(output_c%soil_moisture_vol_layers(k))
        snowpack_data%soil%moisture_liq_layers(k) = REAL(output_c%soil_moisture_liq_layers(k))
    END DO
    DO k = snowpack_data%soil%layer_count + 1, MAX_SOIL_LAYERS
        snowpack_data%soil%temperature_layers(k) = snowpack_data%soil%temperature
        snowpack_data%soil%moisture_vol_layers(k) = snowpack_data%soil%moisture_volumetric
        snowpack_data%soil%moisture_liq_layers(k) = snowpack_data%soil%moisture_liquid
    END DO

    ! Calculate additional budget components if enabled
    IF (snowpack_data%config%enable_mass_budgets) THEN
        CALL update_mass_budgets_structured(snowpack_data)
    END IF

    IF (snowpack_data%config%enable_energy_budgets) THEN
        CALL update_energy_budgets_structured(snowpack_data)
    END IF

    ! Extract soil state from SNOWPACK ground interface
    ! SNOWPACK provides complete soil physics - this is why we use it instead of NoahMP
    ! CRYOWRF SOURCE: module_sf_snowpacklsm.F:364-366 should extract SNOWPACK soil data
    ! Soil properties extracted from SNOWPACK ElementData structure in bridge.cpp

END SUBROUTINE snowpack_physics_structured

!-------------------------------------------------------------------
SUBROUTINE update_mass_budgets_structured(snowpack_data)
!-------------------------------------------------------------------
    TYPE(SnowpackInterface), INTENT(INOUT) :: snowpack_data

    REAL :: dt_hours, melt_water

    dt_hours = snowpack_data%config%dt / 3600.0

    ! Precipitation contribution
    snowpack_data%budgets%mass_precip = snowpack_data%atmos%precipitation * dt_hours

    ! Accumulate snowfall if temperature is cold enough
    IF (snowpack_data%atmos%temperature < 273.15) THEN
        snowpack_data%budgets%accum_snowfall = snowpack_data%budgets%accum_snowfall + &
                                               snowpack_data%budgets%mass_precip
    END IF

    ! Simple sublimation estimate from latent heat flux
    IF (snowpack_data%surface%heat_latent_kg < 0.0) THEN
        snowpack_data%budgets%mass_sublim = snowpack_data%surface%heat_latent_kg * dt_hours
        snowpack_data%budgets%accum_sublim = snowpack_data%budgets%accum_sublim + &
                                             ABS(snowpack_data%budgets%mass_sublim)
    END IF

    ! Calculate runoff from snow melt following CRYOWRF approach
    ! CRYOWRF: module_sf_snowpacklsm.F:716-717 accumulates melt water
    IF (snowpack_data%budgets%mass_melt < 0.0) THEN  ! Negative mass_melt indicates melting
        melt_water = ABS(snowpack_data%budgets%mass_melt)
        ! Split melt into surface and subsurface runoff (simplified)
        snowpack_data%budgets%runoff_surface = melt_water * 0.7      ! 70% surface runoff
        snowpack_data%budgets%runoff_subsurface = melt_water * 0.3   ! 30% subsurface runoff
    ELSE
        snowpack_data%budgets%runoff_surface = 0.0
        snowpack_data%budgets%runoff_subsurface = 0.0
    END IF

END SUBROUTINE update_mass_budgets_structured

!-------------------------------------------------------------------
SUBROUTINE update_energy_budgets_structured(snowpack_data)
!-------------------------------------------------------------------
    TYPE(SnowpackInterface), INTENT(INOUT) :: snowpack_data
    
    REAL :: dt_sec
    REAL :: lw_out_flux, ground_flux, rain_flux, total_flux
    
    dt_sec = snowpack_data%config%dt
    lw_out_flux = snowpack_data%budgets%energy_lw_out
    ground_flux = snowpack_data%budgets%energy_ground
    rain_flux = snowpack_data%budgets%energy_rain
    total_flux = snowpack_data%budgets%energy_total
    
    ! Energy budget components [J/m²] over timestep
    snowpack_data%budgets%energy_sw_in = snowpack_data%atmos%shortwave_down * dt_sec
    snowpack_data%budgets%energy_sw_out = snowpack_data%atmos%shortwave_down * &
                                          snowpack_data%surface%albedo * dt_sec
    snowpack_data%budgets%energy_lw_in = snowpack_data%atmos%longwave_down * dt_sec
    snowpack_data%budgets%energy_sensible = snowpack_data%surface%heat_sensible * dt_sec
    snowpack_data%budgets%energy_latent = snowpack_data%surface%heat_latent * dt_sec
    snowpack_data%budgets%energy_lw_out = lw_out_flux * dt_sec
    snowpack_data%budgets%energy_ground = ground_flux * dt_sec
    snowpack_data%budgets%energy_rain = rain_flux * dt_sec
    snowpack_data%budgets%energy_total = total_flux * dt_sec
    
END SUBROUTINE update_energy_budgets_structured

!-------------------------------------------------------------------
SUBROUTINE extract_layer_diagnostics_structured(snowpack_data)
!-------------------------------------------------------------------
    ! This routine is now obsolete - layer extraction happens directly in snowpack_physics_structured
    ! Kept for backward compatibility but does nothing
    TYPE(SnowpackInterface), INTENT(INOUT) :: snowpack_data
    ! Layer data is now populated directly from SNOWPACK in snowpack_physics_structured
END SUBROUTINE extract_layer_diagnostics_structured

END MODULE module_snowpack_types
