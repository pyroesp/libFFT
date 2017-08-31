# Fast Fourier Transformation


This is my Radix-2 FFT library, which was first developped for an embedded system consisting of a Texas Instrument ARM Cortex-M0 with a 128x96 OLED display.

It is based on the schematic representation of a Radix-2 FFT:

<img src="http://www.nicolaselectronics.be/wp-content/uploads/2013/06/FFT.gif">

## Dependencies:

* Only std includes used

The library should be cross-platform compatible.

## How to change number of FFT points

There are 3 defines you need to change, to modify the number of FFT points:

```C
	/* N */
	#define FFT_POINT 512

	/* N/2 */
	#define FFT_POINT_2 256

	/* log(N)/log(2) */
	#define FFT_STAGES 9
```

For a 1024 point FFT you have to change the values to 1024, 512, 10, respectively.

**Note:** The library (and Radix-2 fft) has been made in such a way that the FFT point value is expected to be a value of 2^X. If you try a different value, it is highly likely that it will not work as intended or not work at all.

## Library functions:

The functions are divided in two groups:

* Functions to execute only once:

```C
	void fft_BlockPerStage(uint16_t *pblocks);
	void fft_ButterfliesPerBlocks(uint16_t *pbutterflies);
	void fft_BitReversedLUT(uint16_t *pbit_reversed);
	void fft_TwiddleFactor(Complex *pW);
```

* Functions to execute continuously before and after each FFT computation:

```C
	void fft_DataToComplex(float *px, Complex *pdata_complex, uint16_t *pbit_reversed);
	void fft_Compute(Complex *pdata_complex, Complex *pW, uint16_t *pblocks, uint16_t *pbutterflies);
	void fft_ComplexToMagnPhase(Complex *pdata_complex, FFT *pspectrum, uint8_t normalize);
```

## License

Creative Commons Attribution-ShareAlike 4.0 International, see LICENSE.md.


This library was first developped for an embedded system consisting of a Texas Instrument ARM microcontroller.
The library only uses standard includes, so it's cross-platform compatible.
