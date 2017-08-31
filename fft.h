#ifndef FFT_H_INCLUDED
#define FFT_H_INCLUDED

#include <stdint.h>
/* math lib */
#include <math.h>

/* N */
#define FFT_POINT 512
/* N/2 */
#define FFT_POINT_2 256
/* log(N)/log(2) */
#define FFT_STAGES 9

/*
 * Uncomment to calculate phase
 * #define FFT_PHASE_USE
 *
 */

#ifdef FFT_PHASE_USE
	/*
	 * Convert rad to degree
	 * (180 / pi)
	 */
	#define RAD_TO_DEGREE 57.2957795130823208768
#endif

/* Complex structure */
typedef struct{
	float re;
	float im;
}Complex;

/* FFT structure */
typedef struct{
	float dB;
	float mag;
	float phase;
}FFT;

/** Use once at start of program **/
/*
 * Blocks per stage
 * N/(2^stage)
 */
void fft_BlockPerStage(uint16_t *pblocks);
/*
 * Butterflies per block per stage
 *	 2^(stage-1)
 */
void fft_ButterfliesPerBlocks(uint16_t *pbutterflies);
/* Bit Reversed LUT calculation */
void fft_BitReversedLUT(uint16_t *pbit_reversed);
/* Calculate twiddle factor */
void fft_TwiddleFactor(Complex *pW);

/** Use before each FFT **/
/* Convert x(n) array to a bit reversed complex array*/
void fft_DataToComplex(float *px, Complex *pdata_complex, uint16_t *pbit_reversed);

/** FFT funtcion **/
/* Compute FFT algorithm */
void fft_Compute(Complex *pdata_complex, Complex *pW, uint16_t *pblocks, uint16_t *pbutterflies);

/*
 * Convert real & imaginary to magnitude & phase 
 * Normalize the magnitude if normalize is not NULL
 */
void fft_ComplexToMagnPhase(Complex *pdata_complex, FFT *pspectrum, uint8_t normalize);

/* 
 * Convert real & imaginary to a dB amplitude 
 * Complex to polar + normalize + 20 * log(...)
 */
void fft_ComplexTodB(Complex *pdata_complex, FFT *pspectrum);
	
#endif
