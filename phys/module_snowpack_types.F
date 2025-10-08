!WRF:MODEL_LAYER:PHYSICS
!
! SNOWPACK Data Structure Definitions
!
MODULE module_snowpack_types

IMPLICIT NONE

! Maximum number of snow layers for arrays
INTEGER, PARAMETER :: MAX_SNOW_LAYERS = 100

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

    ! Initialize atmospheric variables (CRYOWRF compatibility)
    interface_data%surface%friction_velocity = 0.0
    interface_data%surface%stability_param = 0.0

    ! Initialize soil variables
    interface_data%soil%moisture_avail = 1.0      ! Maximum availability
    interface_data%soil%moisture_total = 100.0    ! 100mm total water
    interface_data%soil%moisture_volumetric = 0.3 ! 30% volumetric moisture
    interface_data%soil%moisture_liquid = 0.3     ! All moisture is liquid
    interface_data%soil%temperature = 273.15      ! Freezing point
    
END SUBROUTINE initialize_snowpack_interface

!-------------------------------------------------------------------
SUBROUTINE pack_wrf_to_interface(interface_data, &
    ! Atmospheric inputs
    t_phy, qv_curr, u_phy, v_phy, p8w, swdown, glw, rainbl, &
    ! Grid information  
    xlat, xlong, dz8w, &
    ! Model configuration
    dt, itimestep, &
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

    interface_data%atmos%temperature = t_phy
    interface_data%atmos%humidity = rh_calc  ! Properly calculated RH, not raw QV!
    interface_data%atmos%wind_u = u_phy
    interface_data%atmos%wind_v = v_phy
    interface_data%atmos%pressure = p8w
    interface_data%atmos%shortwave_down = swdown
    interface_data%atmos%longwave_down = glw
    interface_data%atmos%precipitation = rainbl
    interface_data%atmos%height = dz8w
    interface_data%atmos%latitude = xlat
    interface_data%atmos%longitude = xlong
    
    ! Debug first few grid points
    IF (i <= 2 .AND. j <= 2) THEN
        WRITE(*,'(A,I3,A,I3,A,F10.5,A,F10.5)') &
            "SNOWPACK-DEBUG [pack_wrf_to_interface]: Grid(", i, ",", j, &
            ") Packed lat=", xlat, " lon=", xlong
    ENDIF
    
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
    
    ! Pack model configuration
    interface_data%config%dt = dt
    interface_data%config%timestep_number = itimestep
    interface_data%config%i_grid = i
    interface_data%config%j_grid = j
    
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

    ! Note: CRYOWRF atmospheric variables (friction_velocity, stability_param)
    ! are available in interface_data%surface% but not used in standard WRF LSM interface
    
END SUBROUTINE unpack_interface_to_wrf

!-------------------------------------------------------------------
SUBROUTINE snowpack_physics_structured(interface_data)
!-------------------------------------------------------------------
! Structured interface to SNOWPACK physics with layer extraction
!-------------------------------------------------------------------
    USE, INTRINSIC :: iso_c_binding
    TYPE(SnowpackInterface), INTENT(INOUT) :: interface_data
    
    ! C interface declaration (avoiding circular dependency)
    INTERFACE
        SUBROUTINE snowpack_physics_layers_c( &
            temp_air, humidity, wind_speed, wind_dir, &
            shortwave_in, longwave_in, precipitation, pressure, height, dt, &
            i_grid, j_grid, wrf_lat, wrf_lon, &
            snow_swe, snow_depth, surface_temp, &
            heat_flux_sensible, heat_flux_latent, albedo, snow_coverage, &
            friction_velocity, stability_param, &
            n_layers, &
            layer_temp, layer_thick, layer_vol_ice, layer_vol_water, layer_vol_air, &
            layer_grain_radius, layer_bond_radius, layer_dendricity, layer_sphericity, &
            mass_precip, mass_sublim, mass_melt, mass_swe, mass_refreeze, &
            energy_lw_in, energy_lw_out, energy_sw_in, energy_sw_out, &
            energy_sensible, energy_latent, energy_ground_flux, energy_rain, energy_total) BIND(C, name="snowpack_physics_layers")
          
          USE, INTRINSIC :: iso_c_binding
          REAL(C_DOUBLE), VALUE, INTENT(IN) :: temp_air, humidity, wind_speed, wind_dir
          REAL(C_DOUBLE), VALUE, INTENT(IN) :: shortwave_in, longwave_in, precipitation, pressure, height, dt
          INTEGER(C_INT), VALUE, INTENT(IN) :: i_grid, j_grid
          REAL(C_DOUBLE), VALUE, INTENT(IN) :: wrf_lat, wrf_lon
          
          REAL(C_DOUBLE), INTENT(INOUT) :: snow_swe, snow_depth, surface_temp
          REAL(C_DOUBLE), INTENT(OUT) :: heat_flux_sensible, heat_flux_latent, albedo, snow_coverage
          REAL(C_DOUBLE), INTENT(OUT) :: friction_velocity, stability_param
          
          INTEGER(C_INT), INTENT(OUT) :: n_layers
          REAL(C_DOUBLE), INTENT(OUT) :: layer_temp(100), layer_thick(100)
          REAL(C_DOUBLE), INTENT(OUT) :: layer_vol_ice(100), layer_vol_water(100), layer_vol_air(100)
          REAL(C_DOUBLE), INTENT(OUT) :: layer_grain_radius(100), layer_bond_radius(100)
          REAL(C_DOUBLE), INTENT(OUT) :: layer_dendricity(100), layer_sphericity(100)
          REAL(C_DOUBLE), INTENT(OUT) :: mass_precip, mass_sublim, mass_melt, mass_swe, mass_refreeze
          REAL(C_DOUBLE), INTENT(OUT) :: energy_lw_in, energy_lw_out, energy_sw_in, energy_sw_out
          REAL(C_DOUBLE), INTENT(OUT) :: energy_sensible, energy_latent, energy_ground_flux, energy_rain, energy_total
          
        END SUBROUTINE snowpack_physics_layers_c
    END INTERFACE
    
    ! Local variables for C interface
    REAL(8) :: temp_air, humidity, wind_speed, wind_dir
    REAL(8) :: shortwave_in, longwave_in, precipitation, pressure, height, dt
    INTEGER :: i_grid, j_grid
    REAL(8) :: snow_swe, snow_depth, surface_temp
    REAL(8) :: heat_flux_sensible, heat_flux_latent, albedo, snow_coverage
    REAL(8) :: friction_velocity, stability_param
    
    ! Layer arrays for detailed extraction
    INTEGER :: n_layers
    REAL(8), DIMENSION(100) :: layer_temp, layer_thick
    REAL(8), DIMENSION(100) :: layer_vol_ice, layer_vol_water, layer_vol_air
    REAL(8), DIMENSION(100) :: layer_grain_radius, layer_bond_radius
    REAL(8), DIMENSION(100) :: layer_dendricity, layer_sphericity
    
    ! Budget tracking variables
    REAL(8) :: mass_precip, mass_sublim, mass_melt, mass_swe, mass_refreeze
    REAL(8) :: energy_lw_in, energy_lw_out, energy_sw_in, energy_sw_out
    REAL(8) :: energy_sensible, energy_latent, energy_ground_flux, energy_rain, energy_total
    
    INTEGER :: k
    
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
    
    ! Use actual grid indices passed in from WRF
    i_grid = interface_data%config%i_grid
    j_grid = interface_data%config%j_grid
    
    snow_swe = DBLE(interface_data%surface%snow_swe)
    snow_depth = DBLE(interface_data%surface%snow_depth)
    
    ! Call SNOWPACK physics with layer extraction
    CALL snowpack_physics_layers_c(temp_air, humidity, wind_speed, wind_dir, &
                                   shortwave_in, longwave_in, precipitation, pressure, height, dt, &
                                   i_grid, j_grid, DBLE(interface_data%atmos%latitude), DBLE(interface_data%atmos%longitude), &
                                   snow_swe, snow_depth, surface_temp, &
                                   heat_flux_sensible, heat_flux_latent, albedo, snow_coverage, &
                                   friction_velocity, stability_param, &
                                   n_layers, &
                                   layer_temp, layer_thick, layer_vol_ice, layer_vol_water, layer_vol_air, &
                                   layer_grain_radius, layer_bond_radius, layer_dendricity, layer_sphericity, &
                                   mass_precip, mass_sublim, mass_melt, mass_swe, mass_refreeze, &
                                   energy_lw_in, energy_lw_out, energy_sw_in, energy_sw_out, &
                                   energy_sensible, energy_latent, energy_ground_flux, energy_rain, energy_total)
    
    ! Unpack surface results back to interface structure
    interface_data%surface%snow_swe = REAL(snow_swe)
    interface_data%surface%snow_depth = REAL(snow_depth) 
    interface_data%surface%surface_temp = REAL(surface_temp)
    interface_data%surface%heat_sensible = REAL(heat_flux_sensible)
    interface_data%surface%heat_latent = REAL(heat_flux_latent)
    interface_data%surface%albedo = REAL(albedo)
    interface_data%surface%snow_coverage = REAL(snow_coverage)
    interface_data%surface%friction_velocity = REAL(friction_velocity)
    interface_data%surface%stability_param = REAL(stability_param)
    
    ! Store extracted layer data
    interface_data%layers%n_layers = n_layers
    DO k = 1, MIN(n_layers, MAX_SNOW_LAYERS)
        interface_data%layers%layer_temp(k) = REAL(layer_temp(k))
        interface_data%layers%layer_thick(k) = REAL(layer_thick(k))
        interface_data%layers%vol_ice(k) = REAL(layer_vol_ice(k))
        interface_data%layers%vol_water(k) = REAL(layer_vol_water(k))
        interface_data%layers%vol_air(k) = REAL(layer_vol_air(k))
        interface_data%layers%grain_radius(k) = REAL(layer_grain_radius(k))
        interface_data%layers%bond_radius(k) = REAL(layer_bond_radius(k))
        interface_data%layers%dendricity(k) = REAL(layer_dendricity(k))
        interface_data%layers%sphericity(k) = REAL(layer_sphericity(k))
    END DO
    
    ! Clear unused layers
    DO k = n_layers+1, MAX_SNOW_LAYERS
        interface_data%layers%layer_temp(k) = 0.0
        interface_data%layers%layer_thick(k) = 0.0
        interface_data%layers%vol_ice(k) = 0.0
        interface_data%layers%vol_water(k) = 0.0
        interface_data%layers%vol_air(k) = 0.0
        interface_data%layers%grain_radius(k) = 0.0
        interface_data%layers%bond_radius(k) = 0.0
        interface_data%layers%dendricity(k) = 0.0
        interface_data%layers%sphericity(k) = 0.0
    END DO
    
    ! Store mass budget data
    interface_data%budgets%mass_precip = REAL(mass_precip)
    interface_data%budgets%mass_sublim = REAL(mass_sublim)
    interface_data%budgets%mass_melt = REAL(mass_melt)
    
    ! Store energy budget data
    interface_data%budgets%energy_lw_in = REAL(energy_lw_in)
    interface_data%budgets%energy_sensible = REAL(energy_sensible)
    interface_data%budgets%energy_latent = REAL(energy_latent)
    
    ! Calculate additional budget components if enabled
    IF (interface_data%config%enable_mass_budgets) THEN
        CALL update_mass_budgets_structured(interface_data)
    END IF
    
    IF (interface_data%config%enable_energy_budgets) THEN
        CALL update_energy_budgets_structured(interface_data)
    END IF

    ! Update soil state based on surface conditions
    interface_data%soil%temperature = interface_data%surface%surface_temp  ! Soil temp follows surface temp

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
    ! This routine is now obsolete - layer extraction happens directly in snowpack_physics_structured
    ! Kept for backward compatibility but does nothing
    TYPE(SnowpackInterface), INTENT(INOUT) :: interface_data
    ! Layer data is now populated directly from SNOWPACK in snowpack_physics_structured
END SUBROUTINE extract_layer_diagnostics_structured

END MODULE module_snowpack_types