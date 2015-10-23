A 3D Image-Guided Paintbrush

Getting the source code
-----------------------

To build this from source, you must first use git to check the repository
out from GitHub (http://www.github.com/chrisengelsma/3d-painting)

bin/ - platform-dependent scripts (for running demo)
src/ - source code files (e.g., main/java/edu/mines/engelsma/paint/Painting3.java)

Tools for building
------------------

To build CAE, you need these freely available tools:
* Java SE JDK 8.0 (or later):
  http://www.oracle.com/technetwork/java/javase/downloads
* Gradle 2.6
  http://gradle.org
* Mines Java Toolkit (JTK)
  http://www.github.com/dhale/jtk
    Managed by Dave Hale and others at the Colorado School of Mines.
    CAE strongly relies on packages contained within this toolkit.
    It may be checked out and built via SVN. For further instructions on how 
    to build and run the applications contained within the JTK please refer to
    http://boole.mines.edu/jtk/trunk/readme.txt

Building 3D Painting
--------------------

Navigate to the top directory and run `gradlew build`.

Using The Paintbrush
--------------------

After you have built, you should have a JAR file 
[...]/build/libs/3d-painting.jar.
You may include this JAR file as a classpath when running Java.

To use the paintbrush, we must launch a Java virtual machine, specifying all 
of these JAR files and the locations of our JNI libraries. Provided in 
cae/trunk/bin are scripts (e.g. paintdemo.sh) that illustrate how we do this 
for different platforms. 
Similarly, provided with the Mines JTK in jtk/trunk/bin are similar
scripts. To enable CAE, we must add cae.jar to the classpath. This can be
done using a scripting language such as sh or csh as follows:

Using bash, sh or ksh (Unix/Mac):

```bash
export PAINT_HOME=/directories/to/painting/dir
export CLASSPATH=\
$PAINT_HOME/build/libs/3d-painting.jar:\
.
```

Using csh or tcsh (Unix/Mac):

```csh
setenv PAINT_HOME /directories/to/painting/dir
setenv CLASSPATH ${PAINT_HOME}/build/libs/3d-painting.jar:
# more jars here
setenv CLASSPATH ${CLASSPATH}:.
```

Using a batch file (Windows):

```bat
set PAINT_HOME=C:\path\to\painting\dir
set CLASSPATH=^
%PAINT_HOME%\build\libs\3d-painting.jar;^
.
```

By adding the 3d-painting.jar to the classpath this allows you to call 
classes that are included in this lib. 

Demos
-----

Demos are located in ```bin/``` and are platform-dependent scripts. 
Parameters in these scripts must be changed to match your system in order to 
succesfully run.

Tensor Guided Painting
* bin/paintdemo.sh  <-- sh script  (Unix/Mac)
* bin/paintdemo.csh <-- csh script (Unix/Mac)
* bin/paintdemo.bat <-- batch file (Windows)

References
----------
The code in this repo is the result of my MS reserach at the Colorado School of Mines.
