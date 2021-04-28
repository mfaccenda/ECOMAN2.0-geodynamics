@echo off
set APPDIR=%~dp0
set CMD_LINE_ARGS=%1
shift
:getArgs
if " "%1" "==" "" " goto doneArgs
set CMD_LINE_ARGS=%CMD_LINE_ARGS% %1
shift
goto getArgs
:doneArgs

if "%JAVA%"=="" set JAVA=java
if "%TAUP_HOME%"=="" GOTO FIND
echo TAUP_HOME is no longer used and will be ignored
:FIND
PUSHD %APPDIR%
cd ..
set TAUP_HOME=%CD%
POPD

set LIB=%TAUP_HOME%\lib

set TAUP=%LIB%\TauP-2.4.5.jar
set SEISFILE=%LIB%\seisFile-1.8.0.jar
set SEEDCODEC=%LIB%\seedCodec-1.0.11.jar
set SLF4JLOG4J12=%LIB%\slf4j-log4j12-1.7.21.jar
set SLF4JAPI=%LIB%\slf4j-api-1.7.21.jar
set LOG4J=%LIB%\log4j-1.2.17.jar


if EXIST "%TAUP%" GOTO LIBEND
echo %TAUP% doesn't exist
echo TAUP requires this file to function.  It should be in the lib dir
echo parallel to the bin directory to this script in the filesystem.
echo If it seems like the lib dir is there, email sod@seis.sc.edu for help
GOTO END
:LIBEND
    
set CLASSPATH=%TAUP%;%SEISFILE%;%SEEDCODEC%;%SLF4JLOG4J12%;%SLF4JAPI%;%LOG4J%

%JAVA% -classpath %CLASSPATH%   -Xmx512m    edu.sc.seis.TauP.TauP_Path  %CMD_LINE_ARGS%
:END
