@echo off
rem Demos 3D interactive image-guided painting.
rem Author: Chris Engelsma
rem Version: 2015.10.22

setlocal

rem Where is the Java Runtime Environment (JRE)?
set JRE_HOME=C:\pro\jdk\jre

rem Where is the painting repository?
set PAINT_HOME=C:\user\box\git\3d-painting

rem Where will Java look for classes?
rem Add other jars to this list as necessary.
set CLASSPATH=^
%PAINT_HOME%\build\libs\3d-painting.jar;^
%PAINT_HOME%\libs\edu_mines_jtk.jar;^
%PAINT_HOME%\libs\gluegen-rt.jar;^
%PAINT_HOME%\libs\jogl-all.jar;^
.

rem Run a server VM with assertions enabled and a 1GB max Java heap.
rem Modify these flags as necessary for your system.
java -server -ea -Xmx1g ^
-Djava.library.path=%JAVA_LIBRARY_PATH% ^
-Djava.util.logging.config.file=C:\user\etc\java_logging_config ^
edu.mines.engelsma.paint.TensorGuidedPainting

endlocal
