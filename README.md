# linesearcher
Find all the lines positions in a spectrum
Sergio Sousa
This is the README file for the LINESEARCHER code variation from ARES code
This package can be found in http://www.astro.up.pt/~sousasag/ares/

If you got this package from another place, please report to sousasag@astro.up.pt
or go to http://www.astro.up.pt/~sousasag/ares/ to fill the form to download 
the latest version of this code.


1-CONTENTS OF THE PACKAGE
----------
linesearcher.c	- source code for LINESEARCHER
input.search	- input parameters for ARES
README	- this file
sun_harps_ganymede.lines - example of a output file
logsearch.txt   - log file from ARES
sun_harps_ganymede.fits  - Spectrum from the Sun via Ganymede corrected from the radial velocity, taken from: http://www.ls.eso.org/lasilla/sciops/3p6/harps/monitoring/sun.html

----------

2-SYSTEM REQUIREMENTS:
-------------------
cfitsio - CFITSIO - http://heasarc.nasa.gov/fitsio/fitsio.html
gcc - GNU Compiler Colection - http://gcc.gnu.org/
gsl - GNU Scientific Library - http://www.gnu.org/software/gsl/
-------------------

3-INSTALLATION:
-------------------
You can see the webpages for each packadge enumerated before.

I have the code running in two different machines. In a PC Desktop with the Fedora Core 5, and in my laptop with the Ubuntu dapper.

In both this cases it was very easy to install the packages required to compile/run ARES C++ code.

For CFITSIO packadge i followed the instruction in the website and installations files of CFITSIO.
Careful: Remember where do you set the installation directory...

For the GCC, GSL and plotutils:

In the Fedora Core 5 you can use the YUM to find these packages, or use the graphical interface Package Manager (Add and Remove software) and use the search button to find the packages.

Same in Ubuntu, you can use the APP-GET command, or use the graphical interface Synaptic Package Manager (Add/Remove... in advance mode) and use it to search the packages.
-------------------


4-CODE COMPILATION:
------------------
In my computer u use this command to compile ARES, You can try it, but careful with the CFITSIO installation directories...

>$gcc -o runsearch linesearcher.c -L/usr/local/cfitsio/lib/ -I/usr/local/cfitsio/include/ -lcfitsio -lgsl -lgslcblas -lm

Note: -L and -I is the location of the libraries and the include files of the cfitsio package.

In case of having problems finding the libraries, you can add the path of the files missed to find (e.g. error while loading shared libraries: libgsl.so.0: cannot open shared object file: No such file or directory) to $LD_LIBRARY_PATH.


After the compilation is completed you can create a link into your favorite link directories (e.g. ~/bin)

symbolic link directory.

>$ln -s compiledfile ~/bin/runsearch (it is better to use the full paths of both files)

If you have the ~/bin/ in the $PATH system variable then you can run the program easily by typing:

>$runsearch 

make sure that you have the mine.opt file in the running directory.


------------------


5- INPUT PARAMETERS 'mine.opt' FILE:
---------------------------

specfits	: 1D fits spectrum for the analysis
fileout	: output file for the results
lambdai	: initial wavelength for the search of the lines
lambdaf	: final wavelength for the search of the lines
smoothder	: parameter for the calibration of the search of the lines. Noise smoother for the derivatives.
rejt	: parameter for the calibration of the continuum position.
lineresol	: this parameter sets the line resolution of the input spectra. If the code finds two lines closer than the value set for this parameters, 		  then we take the two lines as one line alone.



Detailed description:

To run the code it is necessary to have a file in the system running
directory named mine.opt that contains the input parameters with a 
specified format. The format of this file can be seen below. The input
parameters required to run are the following:

-specfits 	: is the location of the spectra in the FITS format. 
	  The header of this file must contain the CRVAL1 
	  and CDELT1 keywords. The spectra should be reduced
	  and calibrated in wavelength. It is supposed that 
	  the spectra should have a preliminary normalization
	  to avoid abnormal features in the 1D spectra such 
	  the ones that can appear when using reduced echelle
	  spectra.
-fileout	: the output file. The results for the identified lines
	  are prompted in this file with the following order: 
	  the central line wavelength, the depth of the line respectively.
-lambdai	: initial wavelength of the interval to search the lines.
-lambdaf	: final wavelength of the interval to search the lines.
-smoothder	: value of the smooth \textit{boxcar} to use in the 
	  numerical derivatives. Value 1 implies no smoothing.
-rejt	: parameter required to calibrate the local continuum
	  determination.
-lineresol	: minimum distance in Angstroms between lines in the spectra.






Example of a 'input.search' file: 
-----------------------------
specfits='sun_harps_ganymede.fits'
fileout='sun_harps_ganymede.lines'
lambdai=6000.
lambdaf=6500.
smoothder=4
rejt=0.995
lineresol=0.07
-----------------------------


----------------------------------------
6-OUTPUT RESULT('sun_harps_ganymede.lines'):

_________________________________________
6012.23     0.82531
6012.45     0.97689
6013.16     0.96213
6013.50     0.49826
6013.91     0.97736
6014.76     0.97985
6014.85     0.97235
6015.25     0.96125
6016.65     0.45574
6018.31     0.93221
6019.38     0.95885
6020.02     0.61637
6020.19     0.45431
________________________________________

Output report:
1st column: wavelength spectra of the line
2nd column: estimated the depth of the line


