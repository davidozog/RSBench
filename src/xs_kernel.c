#include "rsbench.h"
#include <immintrin.h>

// Reviewed
void calculate_macro_xs( double * macro_xs, int mat, double E, Input input, CalcDataPtrs data, complex double * sigTfactors ) 
{
  
//    printf("3 - nthreads:%d\n", omp_get_num_threads());
    
	// zero out macro vector
	for( int i = 0; i < 4; i++ )
		macro_xs[i] = 0;

	// for nuclide in mat
  //printf("nucs:%d\n", (data.materials).num_nucs[mat]);
	for( int i = 0; i < (data.materials).num_nucs[mat]; i++ )
	{
		double micro_xs[4];
		int nuc = (data.materials).mats[mat][i];

	//	calculate_micro_xs( micro_xs, nuc, E, input, data, sigTfactors);
	// MicroScopic XS's to Calculate
	double sigT_p, sigT;
	double sigA_p, sigA;
	double sigF_p, sigF;
	double sigE_p, sigE;

	// Calculate Window Index
	double spacing = 1.0 / data.n_windows[nuc];
	int window = (int) ( E / spacing );
	if( window == data.n_windows[nuc] )
		window--;

	// Calculate sigTfactors
	//calculate_sig_T(nuc, E, input, data, sigTfactors );
	double phi;

	for( int i = 0; i < input.numL; i++ )
	{
		phi = data.pseudo_K0RS[nuc][i] * sqrt(E);

		if( i == 1 )
			phi -= - atan( phi );
		else if( i == 2 )
			phi -= atan( 3.0 * phi / (3.0 - phi*phi));
		else if( i == 3 )
			phi -= atan(phi*(15.0-phi*phi)/(15.0-6.0*phi*phi));

		phi *= 2.0;

		sigTfactors[i] = cos(phi) - sin(phi) * _Complex_I;
	}

	// Calculate contributions from window "background" (i.e., poles outside window (pre-calculated)
	Window w = data.windows[nuc][window];
	sigT = E * w.T;
	sigA = E * w.A;
	sigF = E * w.F;

	// Loop over Poles within window, add contributions
  //printf("%d:%d, ", w.start, w.end);
//  omp_set_dynamic(244);
  double t1 = omp_get_wtime();

#pragma omp parallel shared(sigT, sigA, sigF) num_threads(30)
  {

    //printf("3 - nthreads:%d\n", omp_get_num_threads());
    sigT_p = 0;
    sigA_p = 0;
    sigF_p = 0;
//    printf("start:end => %d:%d\n", w.start, w.end);

#pragma simd reduction(+:sigT_p,sigA_p,sigF_p)
//#pragma omp for reduction(+:sigT_p,sigA_p,sigF_p)
//#pragma omp for private(sigT_p,sigA_p,sigF_p)
	  for( int i = w.start; i < w.end; i++ )
	  {
//      printf("4 - nthreads:%d\n", omp_get_num_threads());
      complex double PSIIKI;
    	complex double CDUM;
      //if (w.end > 100) {
      //   _mm_prefetch((const char *)&energy_grid[idx], _MM_HINT_T1); // vprefetch1
      //}
	  	//Pole pole = data.poles[nuc][i];
      //complex double p1 = data.poles[nuc][i].MP_EA;
      //complex double p2 = data.poles[nuc][i].MP_RT;
      //complex double p3 = data.poles[nuc][i].MP_RA;
      //complex double p4 = data.poles[nuc][i].MP_RF;
      //complex double sT = sigTfactors[data.poles[nuc][i].l_value];
	  	PSIIKI = -(0.0 - 1.0 * _Complex_I ) / ( data.poles[nuc][i].MP_EA - sqrt(E) );
	  	CDUM = PSIIKI / E;
	  	sigT_p += creal( data.poles[nuc][i].MP_RT * CDUM * sigTfactors[data.poles[nuc][i].l_value] );
	  	sigA_p += creal( data.poles[nuc][i].MP_RA * CDUM);
	  	sigF_p += creal( data.poles[nuc][i].MP_RF * CDUM);
	  }

    #pragma omp atomic
	  	sigT += sigT_p;
    #pragma omp atomic
	  	sigA += sigA_p;
    #pragma omp atomic
	  	sigF += sigF_p;

  //printf("%d - time: %f\n", omp_get_thread_num(), omp_get_wtime() - t1);
  }

	sigE = sigT - sigA;

	micro_xs[0] = sigT;
	micro_xs[1] = sigA;
	micro_xs[2] = sigF;
	micro_xs[3] = sigE;

#pragma ivdep
		for( int j = 0; j < 4; j++ )
		{
			macro_xs[j] += micro_xs[j] * data.materials.concs[mat][i];
		}
	}

	/* Debug
	printf("E = %.2lf, mat = %d, macro_xs[0] = %.2lf, macro_xs[1] = %.2lf, macro_xs[2] = %.2lf, macro_xs[3] = %.2lf\n",
	E, mat, macro_xs[0], macro_xs[1], macro_xs[2], macro_xs[3] );
	*/
	
}

void calculate_micro_xs( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data, complex double * sigTfactors)
{
	// MicroScopic XS's to Calculate
	double sigT;
	double sigA;
	double sigF;
	double sigE;

	// Calculate Window Index
	double spacing = 1.0 / data.n_windows[nuc];
	int window = (int) ( E / spacing );
	if( window == data.n_windows[nuc] )
		window--;

	// Calculate sigTfactors
	calculate_sig_T(nuc, E, input, data, sigTfactors );

	// Calculate contributions from window "background" (i.e., poles outside window (pre-calculated)
	Window w = data.windows[nuc][window];
	sigT = E * w.T;
	sigA = E * w.A;
	sigF = E * w.F;

	// Loop over Poles within window, add contributions
  //printf("%d:%d, ", w.start, w.end);
	for( int i = w.start; i < w.end; i++ )
	{
    complex double PSIIKI;
  	complex double CDUM;
		Pole pole = data.poles[nuc][i];
		PSIIKI = -(0.0 - 1.0 * _Complex_I ) / ( pole.MP_EA - sqrt(E) );
		CDUM = PSIIKI / E;
		sigT += creal( pole.MP_RT * CDUM * sigTfactors[pole.l_value] );
		sigA += creal( pole.MP_RA * CDUM);
		sigF += creal( pole.MP_RF * CDUM);
	}

	sigE = sigT - sigA;

	micro_xs[0] = sigT;
	micro_xs[1] = sigA;
	micro_xs[2] = sigF;
	micro_xs[3] = sigE;
}

void calculate_sig_T( int nuc, double E, Input input, CalcDataPtrs data, complex double * sigTfactors )
{
	double phi;

	for( int i = 0; i < input.numL; i++ )
	{
		phi = data.pseudo_K0RS[nuc][i] * sqrt(E);

		if( i == 1 )
			phi -= - atan( phi );
		else if( i == 2 )
			phi -= atan( 3.0 * phi / (3.0 - phi*phi));
		else if( i == 3 )
			phi -= atan(phi*(15.0-phi*phi)/(15.0-6.0*phi*phi));

		phi *= 2.0;

		sigTfactors[i] = cos(phi) - sin(phi) * _Complex_I;
	}
}
