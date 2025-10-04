!WRF:MODEL_LAYER:PHYSICS
!
! SNOWPACK Interface Module
! Provides ISO_C_BINDING interface between WRF and SNOWPACK physics
!
      MODULE module_sf_snowpack

      USE, INTRINSIC :: iso_c_binding
      IMPLICIT NONE

      PRIVATE
      PUBLIC :: snowpack_driver

      INTERFACE
        SUBROUTINE snowpack_physics_c( &
            temp_air, humidity, wind_speed, wind_dir, &
            shortwave_in, longwave_in, precipitation, pressure, height, dt, &
            i_grid, j_grid, &
            snow_swe, snow_depth, surface_temp, &
            heat_flux_sensible, heat_flux_latent, albedo, snow_coverage) BIND(C, name="snowpack_physics")
          
          USE, INTRINSIC :: iso_c_binding
          REAL(C_DOUBLE), VALUE, INTENT(IN) :: temp_air, humidity, wind_speed, wind_dir
          REAL(C_DOUBLE), VALUE, INTENT(IN) :: shortwave_in, longwave_in, precipitation, pressure, height, dt
          INTEGER(C_INT), VALUE, INTENT(IN) :: i_grid, j_grid
          
          REAL(C_DOUBLE), INTENT(INOUT) :: snow_swe, snow_depth, surface_temp
          REAL(C_DOUBLE), INTENT(OUT) :: heat_flux_sensible, heat_flux_latent, albedo, snow_coverage
          
        END SUBROUTINE snowpack_physics_c
        
        SUBROUTINE snowpack_physics_layers_c( &
            temp_air, humidity, wind_speed, wind_dir, &
            shortwave_in, longwave_in, precipitation, pressure, height, dt, &
            i_grid, j_grid, &
            snow_swe, snow_depth, surface_temp, &
            heat_flux_sensible, heat_flux_latent, albedo, snow_coverage, &
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

          REAL(C_DOUBLE), INTENT(INOUT) :: snow_swe, snow_depth, surface_temp
          REAL(C_DOUBLE), INTENT(OUT) :: heat_flux_sensible, heat_flux_latent, albedo, snow_coverage

          INTEGER(C_INT), INTENT(OUT) :: n_layers
          REAL(C_DOUBLE), DIMENSION(50), INTENT(OUT) :: layer_temp, layer_thick
          REAL(C_DOUBLE), DIMENSION(50), INTENT(OUT) :: layer_vol_ice, layer_vol_water, layer_vol_air
          REAL(C_DOUBLE), DIMENSION(50), INTENT(OUT) :: layer_grain_radius, layer_bond_radius
          REAL(C_DOUBLE), DIMENSION(50), INTENT(OUT) :: layer_dendricity, layer_sphericity

          REAL(C_DOUBLE), INTENT(OUT) :: mass_precip, mass_sublim, mass_melt, mass_swe, mass_refreeze
          REAL(C_DOUBLE), INTENT(OUT) :: energy_lw_in, energy_lw_out, energy_sw_in, energy_sw_out
          REAL(C_DOUBLE), INTENT(OUT) :: energy_sensible, energy_latent, energy_ground_flux, energy_rain, energy_total

        END SUBROUTINE snowpack_physics_layers_c
      END INTERFACE

      CONTAINS

      SUBROUTINE snowpack_driver( &
        itimestep, dt, &
        t_phy, qv_curr, u_phy, v_phy, p8w, dz8w, &
        swdown, glw, rainbl, &
        tsk, hfx, qfx, lh, albedo, snowc, &
        snow, snowh, qsfc, &
        sn_T, sn_thick, sn_volI, sn_volW, sn_volV, &
        sn_rg, sn_rb, sn_dd, sn_sp, sn_nlayers, &
        m_budg_precip, m_budg_sublim, m_budg_melt, m_budg_swe, m_budg_refreeze, &
        e_budg_lw_in, e_budg_lw_out, e_budg_sw_in, e_budg_sw_out, &
        e_budg_sensible, e_budg_latent, e_budg_ground_flux, e_budg_rain, e_budg_total, &
        snowpack_acsnom, snowpack_acsnow, &
        ids,ide, jds,jde, kds,kde, &
        ims,ime, jms,jme, kms,kme, &
        its,ite, jts,jte, kts,kte )
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: itimestep
        REAL,    INTENT(IN) :: dt

        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: t_phy
        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: qv_curr
        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: u_phy
        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: v_phy
        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: p8w
        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: dz8w

        REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: swdown
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: glw
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: rainbl

        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: tsk
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: hfx
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: qfx
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: lh
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: albedo
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: snowc

        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: snow
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: snowh
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: qsfc

        REAL, DIMENSION(ims:ime, 1:50, jms:jme), INTENT(INOUT) :: sn_T
        REAL, DIMENSION(ims:ime, 1:50, jms:jme), INTENT(INOUT) :: sn_thick
        REAL, DIMENSION(ims:ime, 1:50, jms:jme), INTENT(INOUT) :: sn_volI
        REAL, DIMENSION(ims:ime, 1:50, jms:jme), INTENT(INOUT) :: sn_volW
        REAL, DIMENSION(ims:ime, 1:50, jms:jme), INTENT(INOUT) :: sn_volV
        REAL, DIMENSION(ims:ime, 1:50, jms:jme), INTENT(INOUT) :: sn_rg
        REAL, DIMENSION(ims:ime, 1:50, jms:jme), INTENT(INOUT) :: sn_rb
        REAL, DIMENSION(ims:ime, 1:50, jms:jme), INTENT(INOUT) :: sn_dd
        REAL, DIMENSION(ims:ime, 1:50, jms:jme), INTENT(INOUT) :: sn_sp
        INTEGER, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: sn_nlayers

        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: m_budg_precip
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: m_budg_sublim
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: m_budg_melt
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: m_budg_swe
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: m_budg_refreeze
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: e_budg_lw_in
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: e_budg_lw_out
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: e_budg_sw_in
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: e_budg_sw_out
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: e_budg_sensible
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: e_budg_latent
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: e_budg_ground_flux
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: e_budg_rain
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: e_budg_total

        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: snowpack_acsnom
        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: snowpack_acsnow

        INTEGER, INTENT(IN) :: ids,ide, jds,jde, kds,kde
        INTEGER, INTENT(IN) :: ims,ime, jms,jme, kms,kme
        INTEGER, INTENT(IN) :: its,ite, jts,jte, kts,kte

        INTEGER :: i, j, k
        REAL :: z_ml, t_ml, q_ml, p_ml, u_ml, v_ml
        REAL :: wind_speed, wind_dir, rh, precip_rate
        REAL :: heat_sens, heat_lat, alb_out, snow_cov
        REAL :: es, qs
        REAL(C_DOUBLE) :: snow_swe_c, snow_depth_c, surface_temp_c
        REAL(C_DOUBLE) :: heat_sens_c, heat_lat_c, albedo_c, snow_cov_c

        INTEGER(C_INT) :: n_layers_c
        REAL(C_DOUBLE), DIMENSION(50) :: layer_temp_c, layer_thick_c
        REAL(C_DOUBLE), DIMENSION(50) :: layer_voli_c, layer_volw_c, layer_volv_c
        REAL(C_DOUBLE), DIMENSION(50) :: layer_rg_c, layer_rb_c, layer_dd_c, layer_sp_c
        REAL(C_DOUBLE) :: mass_precip_c, mass_sublim_c, mass_melt_c, mass_swe_c, mass_refreeze_c
        REAL(C_DOUBLE) :: energy_lw_in_c, energy_lw_out_c, energy_sw_in_c, energy_sw_out_c
        REAL(C_DOUBLE) :: energy_sens_c, energy_lat_c, energy_gflux_c, energy_rain_c, energy_total_c

        REAL, PARAMETER :: svp1 = 0.6112, svp2 = 17.67
        REAL, PARAMETER :: svp3 = 29.65, svpt0 = 273.15, ep_2 = 0.622

!-------------------------------------------------------------------
! Grid loop - process each grid point  
!-------------------------------------------------------------------
    
    WRITE(*,*) 'SNOWPACK-DEBUG [Fortran/module_sf_snowpack.F]: snowpack_driver called, processing grid range i=',its,'-',ite,' j=',jts,'-',jte
    
    ILOOP: DO i = its, ite
    JLOOP: DO j = jts, jte

      ! Extract atmospheric conditions at lowest model level
      z_ml  = 0.5 * dz8w(i,1,j)                      ! Height of first model level [m]
      t_ml  = t_phy(i,1,j)                           ! Temperature [K]
      q_ml  = qv_curr(i,1,j) / (1.0 + qv_curr(i,1,j)) ! Specific humidity [kg/kg]
      p_ml  = p8w(i,1,j)                             ! Pressure [Pa]  
      u_ml  = u_phy(i,1,j)                           ! U wind [m/s]
      v_ml  = v_phy(i,1,j)                           ! V wind [m/s]

      ! Calculate wind speed and direction
      wind_speed = SQRT(u_ml*u_ml + v_ml*v_ml)     ! Wind speed [m/s]
      wind_dir = ATAN2(u_ml, v_ml) * 180.0/3.14159 ! Wind direction [deg]
      
      ! Calculate relative humidity from specific humidity
      es = svp1 * EXP(svp2 * (t_ml - svpt0) / (t_ml - svp3))  ! Saturation vapor pressure
      qs = ep_2 * es / (p_ml/100.0 - es)                       ! Saturation mixing ratio
      rh = MIN(1.0, MAX(0.0, q_ml / qs))                      ! Relative humidity [0-1]

      ! Convert precipitation from [mm/dt] to [mm/s]  
      precip_rate = rainbl(i,j) / dt

      ! Convert to C_DOUBLE for interface
      snow_swe_c = REAL(snow(i,j), C_DOUBLE)
      snow_depth_c = REAL(snowh(i,j), C_DOUBLE)
      surface_temp_c = REAL(tsk(i,j), C_DOUBLE)

      ! Call SNOWPACK physics through C interface with layer diagnostics
      WRITE(*,*) 'SNOWPACK-DEBUG [Fortran/module_sf_snowpack.F]: Calling C interface for grid (',i,',',j,'), temp=',t_ml
      CALL snowpack_physics_layers_c( &
        REAL(t_ml, C_DOUBLE), REAL(rh, C_DOUBLE), &
        REAL(wind_speed, C_DOUBLE), REAL(wind_dir, C_DOUBLE), &
        REAL(swdown(i,j), C_DOUBLE), REAL(glw(i,j), C_DOUBLE), &
        REAL(precip_rate, C_DOUBLE), REAL(p_ml, C_DOUBLE), REAL(z_ml, C_DOUBLE), &
        REAL(dt, C_DOUBLE), INT(i, C_INT), INT(j, C_INT), &
        snow_swe_c, snow_depth_c, surface_temp_c, &
        heat_sens_c, heat_lat_c, albedo_c, snow_cov_c, &
        n_layers_c, &
        layer_temp_c, layer_thick_c, &
        layer_voli_c, layer_volw_c, layer_volv_c, &
        layer_rg_c, layer_rb_c, layer_dd_c, layer_sp_c, &
        mass_precip_c, mass_sublim_c, mass_melt_c, mass_swe_c, mass_refreeze_c, &
        energy_lw_in_c, energy_lw_out_c, energy_sw_in_c, energy_sw_out_c, &
        energy_sens_c, energy_lat_c, energy_gflux_c, energy_rain_c, energy_total_c)

      ! Update WRF variables with SNOWPACK results
      snow(i,j)   = REAL(snow_swe_c)                ! Update SWE
      snowh(i,j)  = REAL(snow_depth_c)              ! Update snow depth
      tsk(i,j)    = REAL(surface_temp_c)            ! Update surface temperature
      hfx(i,j)    = REAL(heat_sens_c)               ! Sensible heat flux [W/m2]
      lh(i,j)     = REAL(heat_lat_c)                ! Latent heat flux [W/m2]
      qfx(i,j)    = REAL(heat_lat_c) / 2.5e6        ! Moisture flux [kg/m2/s] (approx)
      albedo(i,j) = REAL(albedo_c)                  ! Surface albedo [0-1]
      snowc(i,j)  = REAL(snow_cov_c)                ! Snow coverage fraction [0-1]

      ! Store layer diagnostics in WRF CRYOWRF arrays
      sn_nlayers(i,j) = INT(n_layers_c)

      ! CRYOWRF only writes layer data at periodic intervals (snpack_write_dt)
      IF (MOD(itimestep, INT(21600.0/dt)) == 0) THEN
          DO k = 1, MIN(INT(n_layers_c), 50)
              sn_T(i,k,j)     = REAL(layer_temp_c(k))
              sn_thick(i,k,j) = REAL(layer_thick_c(k))
              sn_volI(i,k,j)  = REAL(layer_voli_c(k))
              sn_volW(i,k,j)  = REAL(layer_volw_c(k))
              sn_volV(i,k,j)  = REAL(layer_volv_c(k))
              sn_rg(i,k,j)    = REAL(layer_rg_c(k))
              sn_rb(i,k,j)    = REAL(layer_rb_c(k))
              sn_dd(i,k,j)    = REAL(layer_dd_c(k))
              sn_sp(i,k,j)    = REAL(layer_sp_c(k))
          END DO

          ! Initialize remaining layers to missing value
          DO k = INT(n_layers_c) + 1, 50
              sn_T(i,k,j)     = -999.0
              sn_thick(i,k,j) = -999.0
              sn_volI(i,k,j)  = -999.0
              sn_volW(i,k,j)  = -999.0
              sn_volV(i,k,j)  = -999.0
              sn_rg(i,k,j)    = -999.0
              sn_rb(i,k,j)    = -999.0
              sn_dd(i,k,j)    = -999.0
              sn_sp(i,k,j)    = -999.0
          END DO
      END IF

      ! Store budget data (accumulated each timestep for mass, instantaneous for energy)
      m_budg_precip(i,j)    = m_budg_precip(i,j) + REAL(mass_precip_c * dt)
      m_budg_sublim(i,j)    = m_budg_sublim(i,j) + REAL(mass_sublim_c * dt)
      m_budg_melt(i,j)      = m_budg_melt(i,j) + REAL(mass_melt_c * dt)
      m_budg_swe(i,j)       = REAL(mass_swe_c)
      m_budg_refreeze(i,j)  = m_budg_refreeze(i,j) + REAL(mass_refreeze_c * dt)
      e_budg_lw_in(i,j)     = REAL(energy_lw_in_c)
      e_budg_lw_out(i,j)    = REAL(energy_lw_out_c)
      e_budg_sw_in(i,j)     = REAL(energy_sw_in_c)
      e_budg_sw_out(i,j)    = REAL(energy_sw_out_c)
      e_budg_sensible(i,j)  = REAL(energy_sens_c)
      e_budg_latent(i,j)    = REAL(energy_lat_c)
      e_budg_ground_flux(i,j) = REAL(energy_gflux_c)
      e_budg_rain(i,j)      = REAL(energy_rain_c)
      e_budg_total(i,j)     = REAL(energy_total_c)

      ! Update time-integrated accumulation variables
      snowpack_acsnom(i,j) = m_budg_melt(i,j)
      snowpack_acsnow(i,j) = m_budg_precip(i,j)

    END DO JLOOP
    END DO ILOOP

    WRITE(*,*) 'SNOWPACK-DEBUG [Fortran/module_sf_snowpack.F]: snowpack_driver completed successfully'

  END SUBROUTINE snowpack_driver

END MODULE module_sf_snowpack