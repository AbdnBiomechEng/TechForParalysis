Readme.txt

This zipped folder contains a bunch of functions which you 
might find useful for processing data from C3D files as well
as generating TRC files and making setup files for running 
Opensim Main programs (e.g. scale, ik, id) from the command line
in Matlab.

Once unzipped, add the folder and all subfolders to the Matlab
path using the function 'editpath' from the command window. For Matlab 
versions beyond 2012 the function 'pathtool' has replaced 'editpath'. 

Also download the example data from the Matlab_tools project 
page - https://simtk.org/home/matlab_tools
This folder contains data plus an example Matlab function which acts
like a pipeline to process data using a model that matches the data. 
Run the example pipeline (e.g. opensim_walking_pipeline.m) and use the data
from the ExampleData folder. NOTE - YOU MUST FIRST RUN THE PIPELINE ON THE
STATIC TRIAL TO CREATE A NEW SCALED MODEL AND THEN RUN THE DYNAMIC TRIAL!!

A change to the pipeline tools is the shift from C3Dserver to BTK
(Biomechanical Toolkit). This new way to access the C3D data from programs
like Matlab is a fantasic addition because it handles the data well 
(including creating relevant transforms for force data) and is also cross-
platfrom, hence it can be used in Mac and Linux as well as Windows. This 
zip file contains a folder called BTK which has the relevant files for 
running in a 64-bit version of Matlab for Windows. If this is not what you 
are using you'll need to download a version for your system from the BTK 
project site - 

http://code.google.com/p/b-tk/      (A BIG THANKYOU TO ARNAUD & STEPHANE)

You will have to add the correct BTK folder to the path of Matlab - for 
instructions see -

http://b-tk.googlecode.com/svn/doc/Matlab/0.1/index.html

Please acknowledge the BTK contribution in any scientific publications - 
e.g. All the computations were performed with MATLAB R2008b (Mathworks, USA)
and the open-source Biomechanical ToolKit package for MATLAB 
(http://code.google.com/p/b-tk).

Please also acknowledge Glen Lichtwark and any other relevant contributors
(for example the C3D functions and XML read/write functions)
for any work used in academic publications.

Please inform me of bugs / suggestions to improve as this 
will be an ongoing project.

A big thank you must go to Tim Dorn (University of Melbourne)
for the inspiration for much of these tools with his excellent
C3Dextract toolbox - https://simtk.org/home/c3dtoolbox