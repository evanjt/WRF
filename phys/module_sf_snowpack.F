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
            shortwave_in, longwave_in, precipitation, pressure, height, &
            snow_swe, snow_depth, surface_temp, &
            heat_flux_sensible, heat_flux_latent, albedo, snow_coverage, &
            dt, i_grid, j_grid) BIND(C, name="snowpack_physics")
          
          USE, INTRINSIC :: iso_c_binding
          REAL(C_DOUBLE), VALUE, INTENT(IN) :: temp_air
          REAL(C_DOUBLE), VALUE, INTENT(IN) :: humidity
          REAL(C_DOUBLE), VALUE, INTENT(IN) :: wind_speed
          REAL(C_DOUBLE), VALUE, INTENT(IN) :: wind_dir
          REAL(C_DOUBLE), VALUE, INTENT(IN) :: shortwave_in
          REAL(C_DOUBLE), VALUE, INTENT(IN) :: longwave_in
          REAL(C_DOUBLE), VALUE, INTENT(IN) :: precipitation
          REAL(C_DOUBLE), VALUE, INTENT(IN) :: pressure
          REAL(C_DOUBLE), VALUE, INTENT(IN) :: height
          REAL(C_DOUBLE), VALUE, INTENT(IN) :: dt
          INTEGER(C_INT), VALUE, INTENT(IN) :: i_grid
          INTEGER(C_INT), VALUE, INTENT(IN) :: j_grid
          
          REAL(C_DOUBLE), INTENT(INOUT) :: snow_swe
          REAL(C_DOUBLE), INTENT(INOUT) :: snow_depth
          REAL(C_DOUBLE), INTENT(INOUT) :: surface_temp
          
          REAL(C_DOUBLE), INTENT(OUT) :: heat_flux_sensible
          REAL(C_DOUBLE), INTENT(OUT) :: heat_flux_latent
          REAL(C_DOUBLE), INTENT(OUT) :: albedo
          REAL(C_DOUBLE), INTENT(OUT) :: snow_coverage
          
        END SUBROUTINE snowpack_physics_c
      END INTERFACE

      CONTAINS

      SUBROUTINE snowpack_driver( &
        itimestep, dt, &
        t_phy, qv_curr, u_phy, v_phy, p8w, dz8w, &
        swdown, glw, rainbl, &
        tsk, hfx, qfx, lh, albedo, snowc, &
        snow, snowh, qsfc, &
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

        INTEGER, INTENT(IN) :: ids,ide, jds,jde, kds,kde
        INTEGER, INTENT(IN) :: ims,ime, jms,jme, kms,kme
        INTEGER, INTENT(IN) :: its,ite, jts,jte, kts,kte

        INTEGER :: i, j
        REAL :: z_ml, t_ml, q_ml, p_ml, u_ml, v_ml
        REAL :: wind_speed, wind_dir, rh, precip_rate
        REAL :: heat_sens, heat_lat, alb_out, snow_cov
        REAL :: es, qs
        REAL(C_DOUBLE) :: snow_swe_c, snow_depth_c, surface_temp_c
        REAL(C_DOUBLE) :: heat_sens_c, heat_lat_c, albedo_c, snow_cov_c
        
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

      ! Call SNOWPACK physics through C interface
      WRITE(*,*) 'SNOWPACK-DEBUG [Fortran/module_sf_snowpack.F]: Calling C interface for grid (',i,',',j,'), temp=',t_ml
      CALL snowpack_physics_c( &
        REAL(t_ml, C_DOUBLE), REAL(rh, C_DOUBLE), &
        REAL(wind_speed, C_DOUBLE), REAL(wind_dir, C_DOUBLE), &
        REAL(swdown(i,j), C_DOUBLE), REAL(glw(i,j), C_DOUBLE), &
        REAL(precip_rate, C_DOUBLE), REAL(p_ml, C_DOUBLE), REAL(z_ml, C_DOUBLE), &
        snow_swe_c, snow_depth_c, surface_temp_c, &
        heat_sens_c, heat_lat_c, albedo_c, snow_cov_c, &
        REAL(dt, C_DOUBLE), INT(i, C_INT), INT(j, C_INT))

      ! Update WRF variables with SNOWPACK results
      snow(i,j)   = REAL(snow_swe_c)                ! Update SWE
      snowh(i,j)  = REAL(snow_depth_c)              ! Update snow depth
      tsk(i,j)    = REAL(surface_temp_c)            ! Update surface temperature
      hfx(i,j)    = REAL(heat_sens_c)               ! Sensible heat flux [W/m2]  
      lh(i,j)     = REAL(heat_lat_c)                ! Latent heat flux [W/m2]
      qfx(i,j)    = REAL(heat_lat_c) / 2.5e6        ! Moisture flux [kg/m2/s] (approx)
      albedo(i,j) = REAL(albedo_c)                  ! Surface albedo [0-1]
      snowc(i,j)  = REAL(snow_cov_c)                ! Snow coverage fraction [0-1]

    END DO JLOOP
    END DO ILOOP

    WRITE(*,*) 'SNOWPACK-DEBUG [Fortran/module_sf_snowpack.F]: snowpack_driver completed successfully'

  END SUBROUTINE snowpack_driver

END MODULE module_sf_snowpack