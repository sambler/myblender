echo No explicit msvc version requested, autodetecting version.

call "%~dp0\detect_msvc2013.cmd"
if %ERRORLEVEL% EQU 0 goto DetectionComplete

call "%~dp0\detect_msvc2015.cmd"
if %ERRORLEVEL% EQU 0 goto DetectionComplete

call "%~dp0\detect_msvc2017.cmd"
if %ERRORLEVEL% EQU 0 goto DetectionComplete

echo Compiler Detection failed. Use verbose switch for more information. 
exit /b 1

:DetectionComplete
echo Compiler Detection successfull, detected VS%BUILD_VS_YEAR%
exit /b 0