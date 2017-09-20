#!/usr/bin/env bash
BUILD_CONFIG=Release

# use globs to take into account various possible suffixes: m, u, d
PYTHON_LIBRARY=`ls -d ${PREFIX}/lib/libpython* | head -n 1`
PYTHON_INCLUDE=`ls -d ${PREFIX}/include/python* | head -n 1`

if [ `uname` = 'Darwin' ]; then
    export CXXFLAGS="${CXXFLAGS} -stdlib=libstdc++"
fi

# we live inside the vtktudoss source
SRC_DIR=$RECIPE_DIR/..
cd $SRC_DIR

mkdir -p build
cd build

# Compile the code from the Samples/SampleCode directory
cmake ../ -G "Ninja" \
    -DCMAKE_BUILD_TYPE=${BUILD_CONFIG} \
    -DCMAKE_INSTALL_PREFIX:PATH="${PREFIX}" \
    -DCMAKE_INSTALL_RPATH:PATH="${PREFIX}/lib" \
    -DCMAKE_CXX_FLAGS="${CXXFLAGS}" \
    ${MACOSX_DEPLOYMENT_TARGET:+-DCMAKE_OSX_DEPLOYMENT_TARGET='10.9'} \
    -DWRAP_PYTHON:BOOL=ON \
    -DINSTALL_PYTHON_MODULE_DIR:PATH="${SP_DIR}" \
    -DPYTHON_INCLUDE_DIR:PATH=${PYTHON_INCLUDE} \
    -DPYTHON_LIBRARY:FILEPATH=${PYTHON_LIBRARY}

ninja install
