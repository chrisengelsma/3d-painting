#!/bin/csh
# Demos 3D interactive image-guided painting.
# Author: Chris Engelsma
# Version: 2015.10.22

# Where is the CAE repository? (Where is the build.xml?)
setenv PAINT_HOME ${HOME}/Home/box/git/3d-painting

# Where will Java look for classes?
# Add other jars to this list as necessary.
setenv CLASSPATH ${PAINT_HOME}/build/libs/3d-painting.jar:
setenv CLASSPATH ${CLASSPATH}:${PAINT_HOME}/libs/edu_mines_jtk.jar:
setenv CLASSPATH ${CLASSPATH}:${PAINT_HOME}/libs/gluegen-rt.jar:
setenv CLASSPATH ${CLASSPATH}:${PAINT_HOME}/libs/jogl-all.jar:
setenv CLASSPATH ${CLASSPATH}:.

# Run a server 64-bit VM with assertions enabled and a 1GB max Java heap.
# Modifiy these flags as necessary for your system.
java -server -d64 -ea -Xmx1g \
-Djava.library.path=${JAVA_LIBRARY_PATH} \
-Dapple.awt.graphicsUseQuartz=true \
edu.mines.engelsma.paint.TensorGuidedPainting
