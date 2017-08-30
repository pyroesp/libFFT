Fast Fourrier Transform
=======================

Radix-2 FFT library
	by pyroesp

This work is licensed under a Creative Commons Attribution-ShareAlike 4.0 International License.


There are 3 defines you need to change, to change the number of FFT points:

/* N */
#define FFT_POINT 512
/* N/2 */
#define FFT_POINT_2 256
/* log(N)/log(2) */
#define FFT_STAGES 9

For a 1024 point FFT you have to change the values to 1024, 512, 10, respectively.

Some functions should be used only once for initializing the FFT. 
See the header file comments for more info.


This library was first developped for an embedded system consisting of a Texas Instrument ARM microcontroller.
The library only uses standard includes, so it's cross-platform compatible.
