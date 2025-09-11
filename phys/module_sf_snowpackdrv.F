!WRF:MODEL_LAYER:PHYSICS
!
! SNOWPACK Land Surface Model Driver
! WRF interface driver for SNOWPACK physics integration
!
MODULE module_sf_snowpackdrv

CONTAINS

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
    INTEGER :: I, J

!-------------------------------------------------------------------
!  Call SNOWPACK physics - this replaces the entire CRYOWRF coupling system
!-------------------------------------------------------------------

    ! Diagnostic message to confirm SNOWPACK LSM is active
    IF (itimestep == 1) THEN
        CALL wrf_message('*** SNOWPACK LSM ACTIVATED: Running modern SNOWPACK physics ***')
    END IF

    ! Call basic SNOWPACK driver - CRYOWRF arrays will come from Registry later
    CALL snowpack_driver( &
        itimestep, dt, &
        t_phy, qv_curr, u_phy, v_phy, p8w, dz8w, &   ! Atmospheric state
        swdown, glw, rainbl, &                        ! Forcing  
        tsk, hfx, qfx, lh, albedo, snowc, &           ! Surface fluxes
        snow, snowh, qsfc, &                          ! Snow states
        ids,ide, jds,jde, kds,kde, &                  ! Dimensions  
        ims,ime, jms,jme, kms,kme, &
        its,ite, jts,jte, kts,kte )

    ! Set some standard LSM outputs that SNOWPACK doesn't directly provide
    DO J = jts,jte
    DO I = its,ite
        ! Skip water points
        IF((XLAND(I,J)-1.5) >= 0.) CYCLE
        
        ! Rough approximations for variables SNOWPACK doesn't compute
        GRDFLX(I,J) = 0.06  ! Geothermal heat flux [W/m2] - typical value
        SMSTAV(I,J) = 1.0   ! Soil moisture availability - assume available
        SMSTOT(I,J) = 100.0 ! Total soil water [mm] - placeholder
        SFCRUNOFF(I,J) = SFCRUNOFF(I,J) + 0.0  ! No surface runoff from SNOWPACK yet
        UDRUNOFF(I,J) = UDRUNOFF(I,J) + 0.0    ! No subsurface runoff from SNOWPACK yet
        CANWAT(I,J) = 0.0   ! No canopy water in basic SNOWPACK
        EMISS(I,J) = 0.98   ! Standard snow/soil emissivity
        Z0(I,J) = 0.001     ! Roughness length [m] - typical for snow
        ZNT(I,J) = Z0(I,J)  ! Thermal roughness = momentum roughness
        
        ! Accumulate snowfall and melt (simplified)
        IF(RAINBL(I,J) > 0.0 .AND. T_PHY(I,1,J) < 273.15) THEN
            ACSNOW(I,J) = ACSNOW(I,J) + RAINBL(I,J) * DT
        END IF
        
        ! Simple snow melt accumulation
        IF(SNOWH(I,J) < 0.001 .AND. TSK(I,J) > 273.15) THEN
            ACSNOM(I,J) = ACSNOM(I,J) + MAX(0.0, (TSK(I,J) - 273.15) * 0.1 * DT)
        END IF
        
        ! Set soil properties (simplified - SNOWPACK handles snow layers)
        SMOIS(I,1,J) = 0.3  ! Volumetric soil moisture - placeholder
        SH2O(I,1,J) = 0.3   ! Liquid soil moisture - placeholder  
        TSLB(I,1,J) = TSK(I,J)  ! Soil temperature = surface temperature (simplified)
        
    END DO
    END DO

END SUBROUTINE snowpacklsm

END MODULE module_sf_snowpackdrv