sun-sky
=======

This is an implementation of several sky models:

* The Utah clear sky model (from "A Practical Analytic Model for Daylight")
* A table-driven version of the same for fast lookup
* A version that uses zonal harmonics to produce an approximate spec/glossy/diffuse convolution
* The CIE clear sky model
* The CIE cloudy sky model
* The CIE partly(!) cloudy sky model
    
To build and run a command-line tool exercising the code, use

    c++ SunSky.cpp SunSkyTest.cpp targa.c -o sunsky && ./sunsky

This tool is capable of producing top-down hemisphere, cube-map, and panoramic
views. See "sunsky -h" for options. 

Examples:

    sunsky -t 12  # Show noon sky for the current time of year
    sunsky -t 12 -s clearUtahBRDF -r 0.3  # Show a glossy version of the same sky
