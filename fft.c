#include "fft.h"

/*
 * Blocks per stage
 * N/(2^stage)
 */
void fft_BlockPerStage(uint16_t *pblocks){
	unsigned char uc_i;
	pblocks[0] = FFT_POINT_2;

	for (uc_i = 1; uc_i < FFT_STAGES; ++uc_i){
		pblocks[ uc_i ] = pblocks[ uc_i-1 ] >> 1;
	}
}

/*
 * Butterflies per block per stage
 *	 2^(stage-1)
 */
void fft_ButterfliesPerBlocks(uint16_t *pbutterflies){
	unsigned char uc_i;
	/* stage 0 starts with 1 butterfly per block */
	pbutterflies[0] = 1;

	for (uc_i = 1; uc_i < FFT_STAGES; ++uc_i){
		pbutterflies[ uc_i ] = pbutterflies[ uc_i-1 ] << 1;
	}
}

/* Bit Reversed LUT calculation */
void fft_BitReversedLUT(uint16_t *pbit_reversed){
	uint16_t i, j;
	uint16_t temp_var1;

	for (i = 0; i < FFT_POINT; ++i){
		temp_var1 = 0;
		for (j = 0; j < FFT_STAGES; ++j){
			temp_var1 |= !!(i & (1 << j)) << (FFT_STAGES - 1 - j);
		}

		pbit_reversed[ i ] = temp_var1;
	}
}

/* Calculate twiddle factor */
void fft_TwiddleFactor(Complex *pW){
	uint16_t i;

	/*  n
	 * W = e^(-j * 2 * pi * n / N) = cos(2 * pi * n / N) - sin(2 * pi * n / N) ; with 0 <= n <= N/2
	 *  N
	 */
	for (i = 0; i < FFT_POINT_2; ++i){
		pW[ i ].re = cosf(2 * M_PI * i / FFT_POINT);
		pW[ i ].im = sinf(2 * M_PI * i / FFT_POINT);
	}
}

/* Calculate Window */
void fft_Window(uint8_t type, float *pWin){
	uint16_t i;
	float a0, a1, a2, a3, a4;


	switch (type){
		case FFT_WIN_TRIANGLE:
			for (i = 0; i < FFT_POINT_2; ++i){
				pWin[ i ] = (float)i/(float)FFT_POINT_2;
				pWin[ FFT_POINT - 1 - i ] = pWin[ i ];
			}
			break;
		case FFT_WIN_HANNING:
			for (i = 0; i < FFT_POINT; i++){
				pWin[ i ] = 0.5 * (1.0 - cos(( 2 * M_PI * i) / (float)(FFT_POINT - 1)));
			}
			break;
		case FFT_WIN_HAMMING:
			for (i = 0; i < FFT_POINT; i++){
				pWin[ i ] = 0.54 - 0.46 * cos((2 * M_PI * i) / (float)(FFT_POINT - 1));
			}
			break;
		case FFT_WIN_BLACKMAN:
			a0 = 0.42659;
			a1 = 0.49656;
			a2 = 0.076849;
			for (i = 0; i < FFT_POINT; i++){
				pWin[ i ] = a0 - a1 * cos((2 * M_PI * i) / (float)(FFT_POINT - 1)) + a2 *  cos((4 * M_PI * i) / (float)(FFT_POINT - 1));
			}
			break;
		case FFT_WIN_NUTTALL:
			a0 = 0.355768;
			a1 = 0.487396;
			a2 = 0.144232;
			a3 = 0.012604;
			for (i = 0; i < FFT_POINT; i++){
				pWin[ i ] = a0 - a1 * cos((2 * M_PI * i) / (float)(FFT_POINT - 1)) + a2 *  cos((4 * M_PI * i) / (float)(FFT_POINT - 1)) - a3 * cos((6 * M_PI * i) / (float)(FFT_POINT - 1));
			}
			break;
		case FFT_WIN_FLAT_TOP:
			a0 = 1;
			a1 = 1.93;
			a2 = 1.29;
			a3 = 0.388;
			a4 = 0.028;
			for (i = 0; i < FFT_POINT; i++){
				pWin[ i ] = a0 - a1 * cos((2 * M_PI * i) / (float)(FFT_POINT - 1)) + a2 *  cos((4 * M_PI * i) / (float)(FFT_POINT - 1)) - a3 * cos((6 * M_PI * i) / (float)(FFT_POINT - 1)) + a4 * cos((8 * M_PI * i) / (float)(FFT_POINT - 1));
			}
			break;
		case FFT_WIN_RECTANGLE:
		default:
			for (i = 0; i < FFT_POINT; i++){
				pWin[ i ] = 0;
			}
			break;
	}
}

/* Multiply x(n) by window convert to a bit reversed complex array*/
void fft_DataToComplex(float *px, float *pWin, Complex *pdata_complex, uint16_t *pbit_reversed){
	uint16_t i;

	for (i = 0; i < FFT_POINT; ++i){
		pdata_complex[ i ].re = px[ *pbit_reversed ] * pWin[ *pbit_reversed ];
		pdata_complex[ i ].im = 0;

		++pbit_reversed;
	}
}

/* Compute FFT algorithm */
void fft_Compute(Complex *pdata_complex, Complex *pW, uint16_t *pblocks, uint16_t *pbutterflies){
	/*
	 * Offset of input of second block in stage 0
	 * Needs to change every stage * 2
	 */
	uint16_t block_offset = 2;

	/*
	 * Offset of second input in first butterfly in stage 0
	 * Needs to change every stage * 2
	 */
	uint16_t butterflies_offset = 1;

	/* for loop counters */
	uint16_t cnt_stages;
	uint16_t cnt_blocks;
	uint16_t cnt_butterflies;
	/* W offset */
	uint16_t cnt_twiddle = 0;

	/* Temporary variables */
	float temp_var1 = 0;
	float temp_var2 = 0;
	float temp_var3 = 0;
	float temp_var4 = 0;
	float temp_var1_2 = 0;
	float temp_var3_4 = 0;

	/* Index variable */
	uint16_t idx_upper = 0;
	uint16_t idx_lower = 0;
	uint16_t idx_blocks = 0;

	/* Copy of pdata_complex.re & pdata_complex.im */
	float real = 0;
	float imag = 0;

	/* Loop for each stage */
	for (cnt_stages = 0; cnt_stages < FFT_STAGES; ++cnt_stages){
		/* Loop for each block in stage */
		for (cnt_blocks = 0; cnt_blocks < *pblocks; ++cnt_blocks){
			/* Reset twiddle factor */
			cnt_twiddle = 0;
			/* Calculate index of first input of first butterfly */
			idx_blocks = cnt_blocks * block_offset;

			/* Loop for each butterfly in block */
			for (cnt_butterflies = 0; cnt_butterflies < *pbutterflies; ++cnt_butterflies){
				/* Calculate index of butterfly input */
				idx_upper = idx_blocks + cnt_butterflies;
				idx_lower = idx_upper + butterflies_offset;

				/*
				 * Temporary variables for multiplying lower butterfly input with twiddle factor
				 * lower * W = (lower.re + lower.im) * (W.re + W.im)
				 *   => var1 = lower.re * W.re  ->  real
				 *   => var2 = lower.im * W.im  ->  -real
				 *   => var3 = lower.re * W.im  ->  imaginary
				 *   => var4 = lower.im * W.re  ->  imaginary
				 *	   => var1_2 = var1 - var2
				 *	   => var3_4 = var3 + var4
				 */

				temp_var1 = pdata_complex[ idx_lower ].re * pW[ cnt_twiddle ].re;
				temp_var2 = pdata_complex[ idx_lower ].im * pW[ cnt_twiddle ].im;
				temp_var1_2 = temp_var1 - temp_var2;

				temp_var3 = pdata_complex[ idx_lower ].re * pW[ cnt_twiddle ].im;
				temp_var4 = pdata_complex[ idx_lower ].im * pW[ cnt_twiddle ].re;
				temp_var3_4 = temp_var3 + temp_var4;

				real = pdata_complex[ idx_upper ].re;
				imag = pdata_complex[ idx_upper ].im;

				/* Upper butterfly output calculation */
				pdata_complex[ idx_upper ].re = real + temp_var1_2;
				pdata_complex[ idx_upper ].im = imag + temp_var3_4;

				/* Lower butterfly output calculation */
			   	pdata_complex[ idx_lower ].re = real - temp_var1_2;
			   	pdata_complex[ idx_lower ].im = imag - temp_var3_4;

			   	/* Increase twiddle counter with factor specific to the current stage */
			   	cnt_twiddle += (FFT_POINT_2 / *pbutterflies);
			}
		}
		/* Increase butterfly pointer */
		++pbutterflies;

		/* Increase blocks pointer */
		++pblocks;

		/* Multiply offset by 2 */
		block_offset = block_offset << 1;
		butterflies_offset = butterflies_offset << 1;
	}
}

/* 
 * Convert real & imaginary to maglitude & phase 
 * Normalize magnitude if normalize is not NULL
 */
void fft_ComplexToMagnPhase(Complex *pdata_complex, FFT *pspectrum, uint8_t normalize){
	uint16_t i;
	float real_sq = 0;
	float imag_sq = 0;

	for (i = 0; i < FFT_POINT_2; ++i){
		real_sq = pdata_complex[ i ].re * pdata_complex[ i ].re;
		imag_sq = pdata_complex[ i ].im * pdata_complex[ i ].im;

		if (!normalize)
			pspectrum[ i ].mag = sqrtf(real_sq + imag_sq);
		else
			pspectrum[ i ].mag = sqrtf(real_sq + imag_sq) / (float)FFT_POINT_2;

		#ifdef FFT_PHASE_USE
			if (pdata_complex[ i ].im != 0)
				pspectrum[ i ].phase = atanf(pdata_complex[ i ].re/pdata_complex[ i ].im) * RAD_TO_DEGREE;
			else
				pspectrum[ i ].phase = 0;
		#endif
	}
}

/* 
 * Convert real & imaginary to a dB amplitude 
 * Complex to polar + normalize + 20 * log(...)
 */
void fft_ComplexTodB(Complex *pdata_complex, FFT *pspectrum){
	uint16_t i;
	float real_sq = 0;
	float imag_sq = 0;

	for (i = 0; i < FFT_POINT_2; ++i){
		real_sq = pdata_complex[ i ].re * pdata_complex[ i ].re;
		imag_sq = pdata_complex[ i ].im * pdata_complex[ i ].im;

		pspectrum[ i ].mag = sqrtf(real_sq + imag_sq) / (float)FFT_POINT_2;

		pspectrum[ i ].dB = 20.0 * log10(pspectrum[ i ].mag);

		#ifdef FFT_PHASE_USE
			if (pdata_complex[ i ].im != 0)
				pspectrum[ i ].phase = atanf(pdata_complex[ i ].re/pdata_complex[ i ].im) * RAD_TO_DEGREE;
			else
				pspectrum[ i ].phase = 0;
		#endif
	}
}

