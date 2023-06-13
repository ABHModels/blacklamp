# blacklamp -- lamppost geometry 

## Intensity Profile of Accretion Disk in Lamppost Coronal Geometry

This C++ code calculates the intensity profile of the irradiation of the accretion disk by a point-like source located on the rotational axis of the black hole in the lamppost geometry. The generated data files can be used to tabulate and structure them into the file in .fits format using another script. The created .fits file can be used in the model relxilllp_nk, a relativistic reflection model in the lamppost geometry.  
  
When using blacklamp for your work, please cite the following the following papers:
* _Public Release of RELXILL_NK: A Relativistic Reflection Model for Testing Einstein’s Gravity_; A. B. Abdikamalov, D. Ayzenberg, C. Bambi, T. Dauser, J. A. Garcia et al.; [Astrophys.J. 878 (2019) 2, 91](https://doi.org/10.3847/1538-4357/ab1f89)
* _Reflection Spectra of Accretion Disks Illuminated by Disk-like Coronae_; S. Riaz, A. B. Abdikamalov, D. Ayzenberg, C. Bambi, H. Wang et al., [AAstrophys.J. 925 (2022) 1, 51](https://doi.org/10.3847/1538-4357/ac3827)

## How to download

        git clone https://github.com/ABHModels/blacklamp.git --branch lamppost --single-branch
  
## How to Run

In the standard from, the main.cpp file takes three input arguments:

        ./main a height mdot

where:

_a_ is the spin parameter  
_height_ is the height of the source  
_mdot_ is the accretion rate of the black hole  

The input parameters can be changed according to the medium setup.

## Output

The code generates a data file in the data directory with the name lp_a_height_mdot.dat, where a, height, and mdot are the input parameters. The file contains the following columns:  

rDisk: radial coordinate of the photons hitting the accretion disk  
III: Intensity at the given radial points on the accretion disk  
ED: emission angle of the photons at the emission point on the source  
HD: incident angle of the photons when they hit the accretion disk  

## Dependencies

The code depends on the following header files:  

eqjohannsen.h and johannsen.h: header file for calculating the metric  and other functions in Johannsen spacetime  
rk45.h: header file for implementing the Runge-Kutta method for the differential equations  

## Notes

* The code uses the RK45 method to integrate the differential equations, and it requires the definition of the cache function in the johannsen.h header file. The cache function is used to store the metric and its derivatives for each point in the trajectory of the photon.

* The output file can be used to construct a .fits file using the generatefits.py script. The .fits file can then be used in the relxilllp_nk model for relativistic reflection in the lamppost geometry.  

* Support contact: <relxill_nk@fudan.edu.cn>
