# WRF Build Environment
# Ubuntu 24.04 with custom library builds
FROM ubuntu:24.04 AS base

LABEL description="WRF atmospheric model build environment"

# Build configuration
ENV DEBIAN_FRONTEND=noninteractive \
    WRF_VERSION=v4.7+ \
    BASE_DIR=/home/wrfuser/wrf-build \
    WRF_CONFIGURE_OPTION=34 \
    WPS_CONFIGURE_OPTION=3

# Library paths configuration
ENV DIR=${BASE_DIR}/libs \
    NETCDF=${BASE_DIR}/libs/netcdf \
    JASPERLIB=${BASE_DIR}/libs/grib2/lib \
    JASPERINC=${BASE_DIR}/libs/grib2/include \
    LD_LIBRARY_PATH="${BASE_DIR}/libs/lib:${BASE_DIR}/libs/netcdf/lib:${BASE_DIR}/libs/grib2/lib" \
    PATH=".:/home/wrfuser/.local/bin:${BASE_DIR}/libs/netcdf/bin:${BASE_DIR}/libs/bin:${PATH}"

# Compiler environment
ENV CC=gcc \
    CXX=g++ \
    F77=gfortran \
    FC=gfortran \
    CFLAGS='-O3 -fPIC' \
    CXXFLAGS='-O3 -fPIC' \
    FFLAGS='-O3 -fPIC -fallow-argument-mismatch -fallow-invalid-boz' \
    FCFLAGS='-O3 -fPIC -fallow-argument-mismatch -fallow-invalid-boz' \
    CPP='gcc -E' \
    CXXCPP='g++ -E'

# Install all WRF dependencies for Ubuntu 24.04
RUN apt-get update && apt-get install -y \
    # Build essentials
    build-essential \
    gfortran \
    gcc \
    g++ \
    make \
    cmake \
    m4 \
    perl \
    csh \
    tcsh \
    pkg-config \
    # Network tools
    curl \
    wget \
    git \
    python3 \
    python3-pip \
    # NetCDF, HDF5, and GRIB2 dependencies
    libhdf5-dev \
    libhdf5-serial-dev \
    zlib1g-dev \
    libjpeg-dev \
    # Build tools for custom libraries
    autoconf \
    libtool \
    bzip2 \
    # Optional dependencies
    libproj-dev \
    libopenmpi-dev \
    openmpi-bin \
    # Utilities
    vim \
    time \
    libtirpc-dev \
    && rm -rf /var/lib/apt/lists/*

# Create normal user and directories
RUN useradd -m -s /bin/bash wrfuser && \
    mkdir -p ${BASE_DIR} ${LIB_DIR} ${INCLUDE_DIR} ${BIN_DIR} ${RUN_DIR} && \
    chown -R wrfuser:wrfuser /home/wrfuser

USER wrfuser
WORKDIR ${BASE_DIR}

# Stage 1: Build WRF libraries from source
FROM base AS libraries-build
USER wrfuser
WORKDIR ${BASE_DIR}

# Create libs directory
RUN cd ${BASE_DIR} && \
    mkdir -p libs

# Build MPICH
RUN cd ${BASE_DIR} && \
    unset F90 && \
    unset F90FLAGS && \
    export CFLAGS="$CFLAGS -fcommon" && \
    wget https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compile_tutorial/tar_files/mpich-3.0.4.tar.gz && \
    tar -xf mpich-3.0.4.tar.gz && \
    cd mpich-3.0.4 && \
    ./configure --prefix=$DIR --disable-wrapper-rpath && \
    make -j 32 2>&1 && \
    make install && \
    cd .. && \
    rm -rf mpich*

# Build zlib
RUN cd ${BASE_DIR} && \
    wget https://www2.mmm.ucar.edu/people/duda/files/mpas/sources/zlib-1.2.11.tar.gz && \
    tar xzvf zlib-1.2.11.tar.gz && \
    cd zlib-1.2.11 && \
    ./configure --prefix=$DIR/grib2 && \
    make -j 32 && \
    make install && \
    cd .. && \
    rm -rf zlib*

# Build HDF5
RUN cd ${BASE_DIR} && \
    export LDFLAGS="-L$DIR/lib -L$DIR/grib2/lib" && \
    export CPPFLAGS="-I$DIR/include -I$DIR/grib2/include" && \
    wget https://www2.mmm.ucar.edu/people/duda/files/mpas/sources/hdf5-1.10.5.tar.bz2 && \
    tar -xf hdf5-1.10.5.tar.bz2 && \
    cd hdf5-1.10.5 && \
    ./configure --prefix=$DIR --with-zlib=$DIR/grib2 --enable-fortran --enable-shared && \
    make -j 32 && \
    make install && \
    cd .. && \
    rm -rf hdf5*

# Build NetCDF-C
RUN cd ${BASE_DIR} && \
    export LDFLAGS="-L$DIR/lib -L$DIR/grib2/lib" && \
    export CPPFLAGS="-I$DIR/include -I$DIR/grib2/include" && \
    wget https://github.com/Unidata/netcdf-c/archive/v4.7.2.tar.gz && \
    tar -xf v4.7.2.tar.gz && \
    cd netcdf-c-4.7.2 && \
    ./configure --enable-shared --enable-netcdf4 --disable-filter-testing --disable-dap --prefix=$DIR/netcdf && \
    make -j 32 && \
    make install && \
    cd .. && \
    rm -rf v4.7.2.tar.gz netcdf-c*

# Build NetCDF-Fortran
RUN cd ${BASE_DIR} && \
    export LIBS="-lnetcdf -lhdf5_hl -lhdf5 -lz" && \
    export LDFLAGS="-L$DIR/lib -L$DIR/netcdf/lib -L$DIR/grib2/lib" && \
    export CPPFLAGS="-I$DIR/include -I$DIR/netcdf/include -I$DIR/grib2/include" && \
    wget https://github.com/Unidata/netcdf-fortran/archive/v4.5.2.tar.gz && \
    tar -xf v4.5.2.tar.gz && \
    cd netcdf-fortran-4.5.2 && \
    ./configure --enable-shared --prefix=$DIR/netcdf && \
    make -j 32 && \
    make install && \
    cd .. && \
    rm -rf netcdf-fortran* v4.5.2.tar.gz

# Build libpng
RUN cd ${BASE_DIR} && \
    wget https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compile_tutorial/tar_files/libpng-1.2.50.tar.gz && \
    tar xzvf libpng-1.2.50.tar.gz && \
    cd libpng-1.2.50 && \
    ./configure --prefix=$DIR/grib2 && \
    make -j 32 && \
    make install && \
    cd .. && \
    rm -rf libpng*

# Build jasper
RUN cd ${BASE_DIR} && \
    wget https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compile_tutorial/tar_files/jasper-1.900.1.tar.gz && \
    tar xzvf jasper-1.900.1.tar.gz && \
    cd jasper-1.900.1 && \
    ./configure --prefix=$DIR/grib2 && \
    make -j 32 && \
    make install && \
    cd .. && \
    rm -rf jasper*

# Stage 2: Copy WRF source code (optimized for caching)
FROM libraries-build AS wrf-source

# Copy entire WRF source tree (dockerignore filters out unnecessary files)
COPY --chown=wrfuser:wrfuser . ${BASE_DIR}/WRF/

# Verify WRF source structure and dependencies
RUN cd ${BASE_DIR}/WRF && \
    test -f configure && \
    test -d Registry && \
    test -d phys

# Stage 3: Configure and build WRF 
FROM wrf-source AS wrf-build-attempt
RUN cd ${BASE_DIR}/WRF && \
    export JASPER=$DIR/grib2 && \
    export JASPER_ROOT=$DIR/grib2 && \
    printf "${WRF_CONFIGURE_OPTION}\n1\n" | ./configure && \
    ./compile em_real -j 20

# Stage 3.5: Analysis stage where we can examine the log even if build failed  
FROM wrf-build-attempt AS wrf-build

# Stage 4: Build WPS
FROM wrf-build AS wps-build
RUN cd ${BASE_DIR} && \
    git clone https://github.com/wrf-model/WPS.git && \
    cd WPS && \
    export WRF_DIR=${BASE_DIR}/WRF && \
    printf "${WPS_CONFIGURE_OPTION}\n" | ./configure && \
    ./compile > log.compile 2>&1

# Final stage: Runtime environment
FROM wps-build AS final

# Create convenient links in run directory
RUN ln -sf ${BASE_DIR}/WRF/main/wrf.exe ${RUN_DIR}/ 2>/dev/null || true && \
    ln -sf ${BASE_DIR}/WRF/main/real.exe ${RUN_DIR}/ 2>/dev/null || true && \
    ln -sf ${BASE_DIR}/WRF/main/ndown.exe ${RUN_DIR}/ 2>/dev/null || true && \
    ln -sf ${BASE_DIR}/WRF/main/tc.exe ${RUN_DIR}/ 2>/dev/null || true && \
    ln -sf ${BASE_DIR}/WPS/geogrid/src/geogrid.exe ${RUN_DIR}/ 2>/dev/null || true && \
    ln -sf ${BASE_DIR}/WPS/ungrib/src/ungrib.exe ${RUN_DIR}/ 2>/dev/null || true && \
    ln -sf ${BASE_DIR}/WPS/metgrid/src/metgrid.exe ${RUN_DIR}/ 2>/dev/null || true

# Set up working directory and entrypoint
WORKDIR /home/wrfuser
ENTRYPOINT ["/bin/bash"]

# Container metadata
LABEL description="WRF atmospheric model build environment" \
      wrf_version="4.7+" \
      ubuntu_version="24.04"
