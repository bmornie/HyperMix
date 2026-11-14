@echo off
setlocal

echo Configuring project...

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

if errorlevel 1 exit /b 1

echo.
echo Building Python extension...

cmake --build build --config Release

if errorlevel 1 exit /b 1

echo.
echo Build completed successfully!