set BUILD_CONFIG=Release
set PYTHON_LIBRARY=%PREFIX%\libs\python%PY_VER:~0,1%%PY_VER:~2,1%.lib

cmake . -G "NMake Makefiles" ^
    -Wno-dev ^
    -DPYTHON_INCLUDE_DIR:PATH="%PREFIX%/include" ^
    -DPYTHON_LIBRARY:FILEPATH="%PYTHON_LIBRARY%" ^
    -DINSTALL_PYTHON_MODULE_DIR:PATH="%SP_DIR%" ^
    -DWRAP_PYTHON:BOOL=ON ^
    -DCMAKE_BUILD_TYPE=%BUILD_CONFIG% ^
    -DCMAKE_INSTALL_PREFIX="%LIBRARY_PREFIX%"
if errorlevel 1 exit 1

nmake install
if errorlevel 1 exit 1
