!WRF:MODEL_LAYER:PHYSICS
!
! SNOWPACK Data Structure Definitions
!
MODULE module_snowpack_types

IMPLICIT NONE

! Maximum number of snow layers for arrays
INTEGER, PARAMETER :: MAX_SNOW_LAYERS = 50

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
    
    ! Long-term accumulators [mm] - persistent across timesteps
    REAL    :: accum_snowfall  ! Accumulated snowfall since start
    REAL    :: accum_snowmelt  ! Accumulated snow melt since start
    REAL    :: accum_sublim    ! Accumulated sublimation since start
END TYPE PhysicsBudgets

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
END TYPE ModelConfig

!-------------------------------------------------------------------
! Combined interface structure for clean function calls
!-------------------------------------------------------------------
TYPE :: SnowpackInterface
    TYPE(AtmosForcing)   :: atmos      ! Atmospheric forcing
    TYPE(SurfaceState)   :: surface    ! Surface state and fluxes
    TYPE(SnowLayerData)  :: layers     ! Detailed layer information
    TYPE(PhysicsBudgets) :: budgets    ! Mass/energy budget tracking
    TYPE(ModelConfig)    :: config     ! Model configuration
END TYPE SnowpackInterface

!-------------------------------------------------------------------
! Helper routines for data structure management
!-------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------
SUBROUTINE initialize_snowpack_interface(interface_data, itimestep, dt, dx)
!-------------------------------------------------------------------
    TYPE(SnowpackInterface), INTENT(INOUT) :: interface_data
    INTEGER, INTENT(IN) :: itimestep
    REAL, INTENT(IN) :: dt, dx
    
    ! Initialize configuration
    interface_data%config%timestep_number = itimestep
    interface_data%config%dt = dt
    interface_data%config%dx = dx
    
    ! Set default physics options (can be configured via namelist later)
    interface_data%config%enable_detailed_layers = .TRUE.
    interface_data%config%enable_mass_budgets = .TRUE.
    interface_data%config%enable_energy_budgets = .TRUE.
    interface_data%config%enable_time_series = .FALSE.
    
    ! Initialize layer data
    interface_data%layers%n_layers = 0
    interface_data%layers%layer_temp = 0.0
    interface_data%layers%layer_thick = 0.0
    
    ! Initialize budgets (important - don't reset accumulators!)
    interface_data%budgets%mass_precip = 0.0
    interface_data%budgets%mass_sublim = 0.0
    interface_data%budgets%energy_sw_in = 0.0
    interface_data%budgets%energy_lw_in = 0.0
    ! Note: accum_* variables should persist between calls
    
END SUBROUTINE initialize_snowpack_interface

!-------------------------------------------------------------------
SUBROUTINE pack_wrf_to_interface(interface_data, &
    ! Atmospheric inputs
    t_phy, qv_curr, u_phy, v_phy, p8w, swdown, glw, rainbl, &
    ! Grid information  
    xlat, xlong, dz8w, &
    ! Surface state (input/output)
    tsk, hfx, qfx, lh, albedo, snowc, snow, snowh, qsfc, z0, znt, &
    ! Grid indices
    i, j, k)
!-------------------------------------------------------------------
    TYPE(SnowpackInterface), INTENT(INOUT) :: interface_data
    
    ! WRF variables
    REAL, INTENT(IN) :: t_phy, qv_curr, u_phy, v_phy, p8w
    REAL, INTENT(IN) :: swdown, glw, rainbl
    REAL, INTENT(IN) :: xlat, xlong, dz8w
    REAL, INTENT(INOUT) :: tsk, hfx, qfx, lh, albedo, snowc
    REAL, INTENT(INOUT) :: snow, snowh, qsfc, z0, znt
    INTEGER, INTENT(IN) :: i, j, k
    
    ! Pack atmospheric forcing
    interface_data%atmos%temperature = t_phy
    interface_data%atmos%humidity = qv_curr  
    interface_data%atmos%wind_u = u_phy
    interface_data%atmos%wind_v = v_phy
    interface_data%atmos%pressure = p8w
    interface_data%atmos%shortwave_down = swdown
    interface_data%atmos%longwave_down = glw
    interface_data%atmos%precipitation = rainbl
    interface_data%atmos%height = dz8w
    interface_data%atmos%latitude = xlat
    interface_data%atmos%longitude = xlong
    
    ! Pack current surface state
    interface_data%surface%snow_swe = snow
    interface_data%surface%snow_depth = snowh
    interface_data%surface%snow_coverage = snowc
    interface_data%surface%surface_temp = tsk
    interface_data%surface%albedo = albedo
    interface_data%surface%heat_sensible = hfx
    interface_data%surface%heat_latent = lh
    interface_data%surface%heat_latent_kg = qfx
    interface_data%surface%moisture_flux = qsfc
    interface_data%surface%roughness_mom = z0
    interface_data%surface%roughness_heat = znt
    
END SUBROUTINE pack_wrf_to_interface

!-------------------------------------------------------------------
SUBROUTINE unpack_interface_to_wrf(interface_data, &
    ! Surface state outputs
    tsk, hfx, qfx, lh, albedo, snowc, snow, snowh, qsfc, z0, znt)
!-------------------------------------------------------------------
    TYPE(SnowpackInterface), INTENT(IN) :: interface_data
    
    ! WRF output variables
    REAL, INTENT(OUT) :: tsk, hfx, qfx, lh, albedo, snowc
    REAL, INTENT(OUT) :: snow, snowh, qsfc, z0, znt
    
    ! Unpack surface state
    tsk = interface_data%surface%surface_temp
    hfx = interface_data%surface%heat_sensible
    qfx = interface_data%surface%heat_latent_kg
    lh = interface_data%surface%heat_latent
    albedo = interface_data%surface%albedo
    snowc = interface_data%surface%snow_coverage
    snow = interface_data%surface%snow_swe
    snowh = interface_data%surface%snow_depth
    qsfc = interface_data%surface%moisture_flux
    z0 = interface_data%surface%roughness_mom
    znt = interface_data%surface%roughness_heat
    
END SUBROUTINE unpack_interface_to_wrf

!-------------------------------------------------------------------
SUBROUTINE snowpack_physics_structured(interface_data)
!-------------------------------------------------------------------
! Structured interface to SNOWPACK physics
!-------------------------------------------------------------------
    TYPE(SnowpackInterface), INTENT(INOUT) :: interface_data
    
    ! Local variables for C interface
    REAL(8) :: temp_air, humidity, wind_speed, wind_dir
    REAL(8) :: shortwave_in, longwave_in, precipitation, pressure, height, dt
    INTEGER :: i_grid, j_grid
    REAL(8) :: snow_swe, snow_depth, surface_temp
    REAL(8) :: heat_flux_sensible, heat_flux_latent, albedo, snow_coverage
    
    ! Pack data for C interface
    temp_air = DBLE(interface_data%atmos%temperature)
    humidity = DBLE(interface_data%atmos%humidity)
    wind_speed = SQRT(interface_data%atmos%wind_u**2 + interface_data%atmos%wind_v**2)
    wind_dir = ATAN2(interface_data%atmos%wind_v, interface_data%atmos%wind_u) * 180.0 / 3.14159
    shortwave_in = DBLE(interface_data%atmos%shortwave_down)
    longwave_in = DBLE(interface_data%atmos%longwave_down)  
    precipitation = DBLE(interface_data%atmos%precipitation)
    pressure = DBLE(interface_data%atmos%pressure)
    height = DBLE(interface_data%atmos%height)
    dt = DBLE(interface_data%config%dt)
    
    i_grid = 1
    j_grid = 1
    
    snow_swe = DBLE(interface_data%surface%snow_swe)
    snow_depth = DBLE(interface_data%surface%snow_depth)
    
    ! Call SNOWPACK physics
    CALL snowpack_physics_c(temp_air, humidity, wind_speed, wind_dir, &
                           shortwave_in, longwave_in, precipitation, pressure, height, dt, &
                           i_grid, j_grid, &
                           snow_swe, snow_depth, surface_temp, &
                           heat_flux_sensible, heat_flux_latent, albedo, snow_coverage)
    
    ! Unpack results back to interface structure
    interface_data%surface%snow_swe = REAL(snow_swe)
    interface_data%surface%snow_depth = REAL(snow_depth) 
    interface_data%surface%surface_temp = REAL(surface_temp)
    interface_data%surface%heat_sensible = REAL(heat_flux_sensible)
    interface_data%surface%heat_latent = REAL(heat_flux_latent)
    interface_data%surface%albedo = REAL(albedo)
    interface_data%surface%snow_coverage = REAL(snow_coverage)
    
    ! Update mass budgets if enabled
    IF (interface_data%config%enable_mass_budgets) THEN
        CALL update_mass_budgets_structured(interface_data)
    END IF
    
    ! Update energy budgets if enabled
    IF (interface_data%config%enable_energy_budgets) THEN
        CALL update_energy_budgets_structured(interface_data)
    END IF
    
    ! Extract detailed layer data if enabled
    IF (interface_data%config%enable_detailed_layers) THEN
        CALL extract_layer_diagnostics_structured(interface_data)
    END IF
    
END SUBROUTINE snowpack_physics_structured

!-------------------------------------------------------------------
SUBROUTINE update_mass_budgets_structured(interface_data)
!-------------------------------------------------------------------
    TYPE(SnowpackInterface), INTENT(INOUT) :: interface_data
    
    REAL :: dt_hours
    
    dt_hours = interface_data%config%dt / 3600.0
    
    ! Precipitation contribution
    interface_data%budgets%mass_precip = interface_data%atmos%precipitation * dt_hours
    
    ! Accumulate snowfall if temperature is cold enough
    IF (interface_data%atmos%temperature < 273.15) THEN
        interface_data%budgets%accum_snowfall = interface_data%budgets%accum_snowfall + &
                                               interface_data%budgets%mass_precip
    END IF
    
    ! Simple sublimation estimate from latent heat flux
    IF (interface_data%surface%heat_latent_kg < 0.0) THEN
        interface_data%budgets%mass_sublim = interface_data%surface%heat_latent_kg * dt_hours
        interface_data%budgets%accum_sublim = interface_data%budgets%accum_sublim + &
                                             ABS(interface_data%budgets%mass_sublim)
    END IF
    
END SUBROUTINE update_mass_budgets_structured

!-------------------------------------------------------------------
SUBROUTINE update_energy_budgets_structured(interface_data)
!-------------------------------------------------------------------
    TYPE(SnowpackInterface), INTENT(INOUT) :: interface_data
    
    REAL :: dt_sec
    
    dt_sec = interface_data%config%dt
    
    ! Energy budget components [J/m²] over timestep
    interface_data%budgets%energy_sw_in = interface_data%atmos%shortwave_down * dt_sec
    interface_data%budgets%energy_sw_out = interface_data%atmos%shortwave_down * &
                                          interface_data%surface%albedo * dt_sec
    interface_data%budgets%energy_lw_in = interface_data%atmos%longwave_down * dt_sec
    interface_data%budgets%energy_sensible = interface_data%surface%heat_sensible * dt_sec
    interface_data%budgets%energy_latent = interface_data%surface%heat_latent * dt_sec
    
END SUBROUTINE update_energy_budgets_structured

!-------------------------------------------------------------------
SUBROUTINE extract_layer_diagnostics_structured(interface_data)
!-------------------------------------------------------------------
    TYPE(SnowpackInterface), INTENT(INOUT) :: interface_data
    
    ! For now, just set basic layer info based on total depth
    IF (interface_data%surface%snow_depth > 0.01) THEN
        ! Simple single-layer approximation (will be replaced with actual SNOWPACK data)
        interface_data%layers%n_layers = 1
        interface_data%layers%layer_temp(1) = interface_data%surface%surface_temp
        interface_data%layers%layer_thick(1) = interface_data%surface%snow_depth
        interface_data%layers%vol_ice(1) = 30.0  ! Approximate ice fraction [%]
        interface_data%layers%vol_water(1) = 5.0  ! Approximate water fraction [%]
        interface_data%layers%vol_air(1) = 65.0   ! Approximate air fraction [%]
        interface_data%layers%grain_radius(1) = 0.5  ! Approximate grain size [mm]
        interface_data%layers%bond_radius(1) = 0.1   ! Approximate bond size [mm]
        interface_data%layers%dendricity(1) = 0.5    ! Approximate dendricity
        interface_data%layers%sphericity(1) = 0.7    ! Approximate sphericity
    ELSE
        interface_data%layers%n_layers = 0
    END IF
    
END SUBROUTINE extract_layer_diagnostics_structured

END MODULE module_snowpack_types