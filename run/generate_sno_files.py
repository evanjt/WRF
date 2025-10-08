#!/usr/bin/env python3
"""
Generate SNOWPACK .sno files from WRF initial conditions
Based on CRYOWRF's make_sno_for_land_ERA5.py

Creates initial snowpack state files for each grid point following 
CRYOWRF naming convention: snpack_domain_j_i.sno
"""

import os
import sys
import datetime as dt
import numpy as np

def create_sno_file(i, j, domain, latitude, longitude, elevation, snow_depth, 
                   soil_temp, output_dir="./snowpack_states"):
    """
    Create a .sno file for a single grid point
    
    Args:
        i, j: Grid indices (WRF convention)
        domain: Domain number (1, 2, 3, etc.)
        latitude, longitude: Grid point coordinates [degrees]
        elevation: Grid point elevation [m]
        snow_depth: Initial snow depth [m]
        soil_temp: Soil temperature [K]
        output_dir: Directory for output files
    """
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate filename following CRYOWRF convention
    filename = os.path.join(output_dir, f"snpack_{domain}_{j}_{i}.sno")
    
    # Station ID following CRYOWRF convention
    station_id = f"snpack_{domain}_{j}_{i}"
    
    # Current date for ProfileDate
    current_date = dt.datetime.now()
    
    # Determine number of snow layers
    if snow_depth > 0.001:  # Minimum snow depth threshold
        n_snow_layers = 10
        # Divide snow into 10 equal layers
        layer_thickness = snow_depth / 10.0
    else:
        n_snow_layers = 0
        layer_thickness = 0.0
    
    # Soil layer configuration (4 layers)
    soil_thicknesses = [1.89, 0.72, 0.21, 0.07]  # [m]
    n_soil_layers = 4
    
    # Default snow properties
    snow_density = 300.0  # [kg/m³] - reasonable fresh snow density
    snow_vol_ice = snow_density / 918.0  # Ice density
    snow_vol_air = 1.0 - snow_vol_ice
    snow_temperature = min(soil_temp, 273.15)  # Snow temp <= 0°C
    
    # Default soil properties
    soil_vol_solid = 0.7    # Soil volume fraction
    soil_vol_air = 0.3      # Air volume fraction
    soil_vol_water = 0.0    # Assume dry soil initially
    soil_vol_ice = 0.0      # No frozen soil initially
    soil_density = 2400.0   # [kg/m³]
    soil_conductivity = 2.0 # [W/m/K]
    soil_heat_capacity = 2000.0 # [J/m³/K]
    
    # Open file for writing
    with open(filename, 'w') as f:
        # SMET header
        f.write("SMET 1.1 ASCII\\n")
        f.write("[HEADER]\\n")
        f.write(f"station_id       = {station_id}\\n")
        f.write("station_name     = WRF:SNOWPACK\\n")
        f.write(f"latitude         = {latitude:.4f}\\n")
        f.write(f"longitude        = {longitude:.4f}\\n")
        f.write(f"altitude         = {elevation:.2f}\\n")
        f.write("nodata           = -999\\n")
        f.write(f"source           = WRF_to_SNOWPACK; Generated {current_date.strftime('%Y-%m-%d')}\\n")
        f.write(f"ProfileDate      = {current_date.isoformat()}\\n")
        f.write(f"HS_Last          = {snow_depth:.6f}\\n")
        f.write("SlopeAngle       = 0.00\\n")
        f.write("SlopeAzi         = 0.00\\n")
        f.write(f"nSoilLayerData   = {n_soil_layers}\\n")
        f.write(f"nSnowLayerData   = {n_snow_layers}\\n")
        f.write("SoilAlbedo       = 0.2\\n")
        f.write("BareSoil_z0      = 0.02\\n")
        f.write("CanopyHeight     = 0.00\\n")
        f.write("CanopyLeafAreaIndex = 0.00\\n")
        f.write("CanopyDirectThroughfall = 1.00\\n")
        f.write("WindScalingFactor = 1.00\\n")
        f.write(f"ErosionLevel     = {n_snow_layers}\\n")
        f.write("TimeCountDeltaHS = 0.000000\\n")
        f.write("fields           = timestamp Layer_Thick  T  Vol_Frac_I  Vol_Frac_W  Vol_Frac_V  Vol_Frac_S Rho_S Conduc_S HeatCapac_S  rg  rb  dd  sp  mk mass_hoar ne CDot metamo\\n")
        f.write("[DATA]\\n")
        
        # Soil layers (4 layers from bottom to top)
        for layer in range(n_soil_layers):
            f.write(f"{current_date.isoformat()}   ")
            f.write(f"{soil_thicknesses[layer]:6.3f}   ")  # Layer thickness [m]
            f.write(f"{soil_temp:7.3f}   ")               # Temperature [K]
            f.write(f"{soil_vol_ice:5.3f}   ")            # Ice volume fraction
            f.write(f"{soil_vol_water:5.3f}   ")          # Water volume fraction
            f.write(f"{soil_vol_air:5.3f}   ")            # Air volume fraction
            f.write(f"{soil_vol_solid:5.3f}   ")          # Soil volume fraction
            f.write(f"{soil_density:8.3f}   ")            # Soil density [kg/m³]
            f.write(f"{soil_conductivity:6.3f}   ")       # Heat conductivity [W/m/K]
            f.write(f"{soil_heat_capacity:8.3f}   ")      # Heat capacity [J/m³/K]
            f.write("1000.000   ")                        # Grain radius [mm]
            f.write("0.000   ")                           # Bond radius [mm]
            f.write("0.000   ")                           # Dendricity
            f.write("1.000   ")                           # Sphericity
            f.write("1   ")                               # Grain marker
            f.write("0.0   ")                             # Mass hoar
            f.write("1   ")                               # Number of elements
            f.write("0.0   ")                             # Stress rate
            f.write("0.0\\n")                             # Metamorphism
        
        # Snow layers (if present)
        if n_snow_layers > 0:
            for layer in range(n_snow_layers):
                f.write(f"{current_date.isoformat()}   ")
                f.write(f"{layer_thickness:6.3f}   ")     # Layer thickness [m]
                f.write(f"{snow_temperature:7.3f}   ")    # Temperature [K]
                f.write(f"{snow_vol_ice:5.3f}   ")        # Ice volume fraction
                f.write("0.000   ")                       # Water volume fraction
                f.write(f"{snow_vol_air:5.3f}   ")        # Air volume fraction
                f.write("0.000   ")                       # Soil volume fraction
                f.write("0.000   ")                       # Density (not used for snow)
                f.write("0.000   ")                       # Conductivity (calculated by SNOWPACK)
                f.write("0.000   ")                       # Heat capacity (calculated)
                f.write("0.400   ")                       # Grain radius [mm]
                f.write("0.100   ")                       # Bond radius [mm]
                f.write("0.000   ")                       # Dendricity
                f.write("1.000   ")                       # Sphericity
                f.write("1   ")                           # Grain marker
                f.write("0.0   ")                         # Mass hoar
                f.write("1   ")                           # Number of elements
                f.write("0.0   ")                         # Stress rate  
                f.write("0.0\\n")                         # Metamorphism
    
    print(f"Created .sno file: {filename}")

def generate_simple_test_files(domain=1, nx=5, ny=5):
    """
    Generate simple test .sno files for a small domain
    
    Args:
        domain: Domain number
        nx, ny: Grid dimensions
    """
    print(f"Generating test .sno files for {nx}x{ny} grid, domain {domain}")
    
    # Simple test grid
    for j in range(1, ny + 1):  # WRF uses 1-based indexing in file names
        for i in range(1, nx + 1):
            # Simple test values
            latitude = 45.0 + j * 0.1
            longitude = -110.0 + i * 0.1  
            elevation = 1500.0 + i * 10.0 + j * 5.0
            snow_depth = 0.1 if (i + j) % 3 == 0 else 0.0  # Some points have snow
            soil_temp = 275.0  # 2°C soil temperature
            
            create_sno_file(i, j, domain, latitude, longitude, elevation, 
                           snow_depth, soil_temp)
    
    print(f"Generated {nx * ny} .sno files in ./snowpack_states/")

if __name__ == "__main__":
    print("SNOWPACK .sno file generator for WRF integration")
    print("Based on CRYOWRF methodology")
    print()
    
    if len(sys.argv) > 1 and sys.argv[1] == "--test":
        # Generate test files
        generate_simple_test_files()
    else:
        print("Usage:")
        print("  python generate_sno_files.py --test    # Generate test files")
        print()
        print("To generate from actual WRF files, modify the script to read wrfinput_d0X files")
        print("See CRYOWRF's make_sno_for_land_ERA5.py for full implementation")