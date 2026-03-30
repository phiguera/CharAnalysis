# Directions for downloading and installing _CharAnalsyis_ stand-alone versions on Windows OS:

## Overview

To run _CharAnalysis_ as a stand-alone program, you need to first install the MATLAB Compiler Runtime (MCR) program on your machine. For this *you will need administrator privilages.* You only need to install the MCR program once. 

The version of the MCR program must match the version of _CharAnalysis_ you are using, and it is specific to 
the type of computer you are using (e.g., 64 bit vs. 32 bit machine, Windows 7 vs. Windows XP / Vista). Directions below are specific to two different Windows computer / OS combinations. 

## _CharAnalysis_ 1.1, for Windows 7, 64 bit machines:

### (1) Download CharAnalsyisInstaller.exe:

http://files.cfc.umt.edu/phiguera/CharAnalysis/CharAnalsyisInstaller.exe [4 mb]

### (2) Double click the .exe file and wait...

The program will direct you through the installation of the MATLAB Component Runtime (MCR)program, accessing it from a MATLAB server. The first time you install the program, this will take a few minutes. In order to be able to save output data and figure, you should save the program in a directory that has read/write privileges. Sometimes this is not the case in the default "programs" file on a Windows OS; if you run the program and it closes out before saving the data and figures, move the entire directory to a different location, outside of the "programs" directory.

## _CharAnalysis_ 1.1, for Windows XP or Vista, 32 bit machines: 

### (1) Download CharAnalysis_1_1_WinXP.exe:
http://files.cfc.umt.edu/phiguera/CharAnalysis/CharAnalysis_1_1_WinXP.exe [590 mb]

### (2) Download MATLAB Compiler Runtime (MCR): MCRInstaller_7.11.exe

http://files.cfc.umt.edu/phiguera/CharAnalysis/MCRInstaller_7.11.exe [580 mb]

### (3) Install MCR file by double-clicking on MCRInstaller_*.exe.


### (4) Once the MCR is installed, initiate CharAnalysis_*.exe by double clicking on the icon. 

## Overall notes

The first time the program runs, the startup time will be > 1 minute. After the first time, startup times will be significantly faster. 

The first task one should perform to validate the installation is to download and run the example files from Code Lake, located in the main GitHub repository. The most common problems when running the program for the first time are related to directory issues. Make sure the files you are accessing and writing to are located in the same directory as the program (or the entire path is passed to the program with the file name), and make sure any relevant files (*_charResults.csv, or the .xls file) are closed if you are to save data. 





