accelcomp Introduction
=========

Acceleration to broadband comparison code.  This code compares data from co-located instruments
at various GSN stations.  It is driven by a network code and event CMT.

The data pulls are written in a way specific to the ASL internal data structure.  However,
This function could be easily changed to deal with other users types of data location.

The codes basic features are as follows:

1) Pull the data from all the various instruments as defined in metadata
2) Rotate, deconvolve, filter, and trim the data to all get common units
3) Plot the comparison of the various instruments
4) Save a few residual and data comparison numbers to file

accelcomp Usage
=========

There are various user parameters that can be changed and are located in the top of accelcomp.py

[aringler@aslres01 accelcomp]$ ./accelcomp.py 
Usage: CMT ResultsName Network

[aringler@aslres01 accelcomp]$ ./accelcomp.py /SYNTHETICS/2013/C201309070013A/CMTSOLUTION check IU


accelcomp Need To Do
=========

The data reading function needs to be cleaned up and made not specific to the ASL data structure

The rotations need to be verified as there might be a rotational bug

We need to split out the parameter file so it is more easily used

The plotting feature needs to be pulled out and put in a function

This code should be merged with the synthetics to avoid duplicate functions

Need to combine both the trigger and the non-trigger versions to deal with different units


