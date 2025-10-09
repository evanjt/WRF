# CRYOWRF SNOWPACK & MeteoIO Build Makefile
# Builds SNOWPACK and MeteoIO libraries for WRF integration

# Build configuration
CC ?= gcc
CXX ?= g++
F77 ?= gfortran
FC ?= gfortran

# Paths
BASE_DIR ?= /home/wrfuser/wrf-build
CRYOWRF_DIR := $(shell pwd)
BUILD_DIR := $(CRYOWRF_DIR)/build
METEOIO_DIR := $(CRYOWRF_DIR)/meteoio
SNOWPACK_DIR := $(CRYOWRF_DIR)/snowpack

# Installation directories
METEOIO_INSTALL_DIR := $(CRYOWRF_DIR)/meteoio/lib
SNOWPACK_INSTALL_DIR := $(CRYOWRF_DIR)/snowpack/lib

# Compiler flags
CFLAGS := -O3 -fPIC $(CFLAGS)
CXXFLAGS := -O3 -fPIC -std=c++11 $(CXXFLAGS)
FFLAGS := -O3 -fPIC -fallow-argument-mismatch -fallow-invalid-boz $(FFLAGS)

# NetCDF paths (from WRF build environment)
NETCDF_DIR ?= $(BASE_DIR)/libs/netcdf
NETCDF_LIB := -L$(NETCDF_DIR)/lib -lnetcdf -lnetcdff
NETCDF_INCLUDE := -I$(NETCDF_DIR)/include

# Target libraries
METEOIO_LIB := $(METEOIO_INSTALL_DIR)/libmeteoio.a
SNOWPACK_LIB := $(SNOWPACK_INSTALL_DIR)/libsnowpack.a

# Default target
all: $(SNOWPACK_LIB)

# Help target
help:
	@echo "CRYOWRF SNOWPACK & MeteoIO Build System"
	@echo "========================================"
	@echo ""
	@echo "Targets:"
	@echo "  all        - Build both MeteoIO and SNOWPACK (default)"
	@echo "  meteoio    - Build MeteoIO library only"
	@echo "  snowpack   - Build SNOWPACK library only"
	@echo "  clean      - Clean all build artifacts"
	@echo "  install    - Install libraries to standard locations"
	@echo "  help       - Show this help message"
	@echo ""
	@echo "Environment Variables:"
	@echo "  BASE_DIR   - WRF build directory (default: /home/wrfuser/wrf-build)"
	@echo "  NETCDF_DIR - NetCDF installation directory"
	@echo "  CC, CXX, FC - Compilers"
	@echo ""

# MeteoIO build target
meteoio: $(METEOIO_LIB)

$(METEOIO_LIB):
	@echo "Building MeteoIO..."
	@mkdir -p $(BUILD_DIR)/meteoio
	@mkdir -p $(METEOIO_INSTALL_DIR)
	@if [ ! -d "$(METEOIO_DIR)" ]; then \
		echo "ERROR: MeteoIO source directory not found: $(METEOIO_DIR)"; \
		echo "Please clone MeteoIO to $(METEOIO_DIR)"; \
		exit 1; \
	fi
	cd $(BUILD_DIR)/meteoio && \
	cmake $(METEOIO_DIR) \
		-DCMAKE_INSTALL_PREFIX=$(CRYOWRF_DIR)/meteoio \
		-DCMAKE_C_COMPILER=$(CC) \
		-DCMAKE_CXX_COMPILER=$(CXX) \
		-DCMAKE_Fortran_COMPILER=$(FC) \
		-DCMAKE_C_FLAGS="$(CFLAGS)" \
		-DCMAKE_CXX_FLAGS="$(CXXFLAGS)" \
		-DCMAKE_Fortran_FLAGS="$(FFLAGS)" \
		-DNETCDF_INCLUDE_DIR=$(NETCDF_DIR)/include \
		-DNETCDF_LIBRARY=$(NETCDF_DIR)/lib \
		-DBUILD_PLUGIN=OFF \
		-DBUILD_EXAMPLES=OFF \
		-DBUILD_TESTS=OFF \
		-DENABLE_DBO=OFF \
		-DENABLE_IMIS=OFF \
		-DENABLE_MYSQL=OFF \
		-DENABLE_ODBC=OFF \
		-DENABLE_ORACLE=OFF \
		-DENABLE_POSTGRESQL=OFF \
		-DENABLE_SQLITE=OFF \
		-DENABLE_SMET=ON \
		-DENABLE_GSN=ON \
		-DENABLE_BORMA=ON \
		-DENABLE_GRASS=ON \
		-DENABLE_GEOTOP=ON \
		-DENABLE_PROJ=ON \
		-DENABLE_CURL=ON
	$(MAKE) -C $(BUILD_DIR)/meteoio -j$(shell nproc)
	$(MAKE) -C $(BUILD_DIR)/meteoio install
	@echo "MeteoIO built successfully: $(METEOIO_LIB)"

# SNOWPACK build target
snowpack: $(SNOWPACK_LIB)

$(SNOWPACK_LIB): $(METEOIO_LIB)
	@echo "Building SNOWPACK..."
	@mkdir -p $(BUILD_DIR)/snowpack
	@mkdir -p $(SNOWPACK_INSTALL_DIR)
	@if [ ! -d "$(SNOWPACK_DIR)" ]; then \
		echo "ERROR: SNOWPACK source directory not found: $(SNOWPACK_DIR)"; \
		echo "Please clone SNOWPACK to $(SNOWPACK_DIR)"; \
		exit 1; \
	fi
	cd $(BUILD_DIR)/snowpack && \
	cmake $(SNOWPACK_DIR) \
		-DCMAKE_INSTALL_PREFIX=$(CRYOWRF_DIR)/snowpack \
		-DCMAKE_C_COMPILER=$(CC) \
		-DCMAKE_CXX_COMPILER=$(CXX) \
		-DCMAKE_Fortran_COMPILER=$(FC) \
		-DCMAKE_C_FLAGS="$(CFLAGS)" \
		-DCMAKE_CXX_FLAGS="$(CXXFLAGS)" \
		-DCMAKE_Fortran_FLAGS="$(FFLAGS)" \
		-DNETCDF_INCLUDE_DIR=$(NETCDF_DIR)/include \
		-DNETCDF_LIBRARY=$(NETCDF_DIR)/lib \
		-DENABLE_CRYOWRF=ON \
		-DBUILD_PLUGIN=OFF \
		-DBUILD_EXAMPLES=OFF \
		-DBUILD_TOOLS=OFF \
		-DBUILD_TESTS=OFF \
		-DOPERATIONAL=ON \
		-DMAKE_METEOIO_LIB=ON \
		-DMETEOIO_BASE=$(CRYOWRF_DIR)/meteoio \
		-DMAKE_FORTRAN=ON
	$(MAKE) -C $(BUILD_DIR)/snowpack -j$(shell nproc)
	$(MAKE) -C $(BUILD_DIR)/snowpack install
	@echo "SNOWPACK built successfully: $(SNOWPACK_LIB)"

# Clean target
clean:
	@echo "Cleaning build artifacts..."
	rm -rf $(BUILD_DIR)
	rm -f $(METEOIO_LIB)
	rm -f $(SNOWPACK_LIB)
	@echo "Clean completed"

# Install target (for potential system-wide installation)
install: $(SNOWPACK_LIB)
	@echo "Installing CRYOWRF libraries..."
	@mkdir -p $(BASE_DIR)/lib
	@mkdir -p $(BASE_DIR)/include
	cp $(METEOIO_LIB) $(BASE_DIR)/lib/ 2>/dev/null || echo "MeteoIO lib not found"
	cp $(SNOWPACK_LIB) $(BASE_DIR)/lib/ 2>/dev/null || echo "SNOWPACK lib not found"
	cp -r $(CRYOWRF_DIR)/meteoio/include/* $(BASE_DIR)/include/ 2>/dev/null || echo "MeteoIO headers not found"
	cp -r $(CRYOWRF_DIR)/snowpack/include/* $(BASE_DIR)/include/ 2>/dev/null || echo "SNOWPACK headers not found"
	@echo "Installation completed"

# Status target
status:
	@echo "CRYOWRF Build Status"
	@echo "===================="
	@echo "Base Directory: $(BASE_DIR)"
	@echo "NetCDF Directory: $(NETCDF_DIR)"
	@echo "CRYOWRF Directory: $(CRYOWRF_DIR)"
	@echo ""
	@echo "MeteoIO Source: $(METEOIO_DIR)"
	@echo "MeteoIO Library: $(METEOIO_LIB)"
	@echo "SNOWPACK Source: $(SNOWPACK_DIR)"
	@echo "SNOWPACK Library: $(SNOWPACK_LIB)"
	@echo ""
	@if [ -f "$(METEOIO_LIB)" ]; then \
		echo "MeteoIO Status: ✓ Built"; \
		ls -la $(METEOIO_LIB); \
	else \
		echo "MeteoIO Status: ✗ Not built"; \
	fi
	@if [ -f "$(SNOWPACK_LIB)" ]; then \
		echo "SNOWPACK Status: ✓ Built"; \
		ls -la $(SNOWPACK_LIB); \
	else \
		echo "SNOWPACK Status: ✗ Not built"; \
	fi

# Phony targets
.PHONY: all meteoio snowpack clean install status help