#include "rsbench.h"

// JRT 
int * generate_n_poles( Input input )
{
	int total_resonances = input.avg_n_poles * input.n_nuclides;

	int * R = (int *) calloc( input.n_nuclides, sizeof(int));

	for( int i = 0; i < total_resonances; i++ )
		R[rand() % input.n_nuclides]++;

	// Ensure all nuclides have at least 1 resonance
	for( int i = 0; i < input.n_nuclides; i++ )
		if( R[i] == 0 )
			R[i] = 1;
	
	/* Debug	
	for( int i = 0; i < input.n_nuclides; i++ )
		printf("R[%d] = %d\n", i, R[i]);
	*/
	
	return R;
}

int * generate_n_windows( Input input )
{
	int total_resonances = input.avg_n_windows * input.n_nuclides;

	int * R = (int *) calloc( input.n_nuclides, sizeof(int));

	for( int i = 0; i < total_resonances; i++ )
		R[rand() % input.n_nuclides]++;

	// Ensure all nuclides have at least 1 resonance
	for( int i = 0; i < input.n_nuclides; i++ )
		if( R[i] == 0 )
			R[i] = 1;
	
	/* Debug	
	for( int i = 0; i < input.n_nuclides; i++ )
		printf("R[%d] = %d\n", i, R[i]);
	*/
	
	return R;
}

MIC_Pole generate_mypoles( Input input, int * n_poles )
{
  complex double *mp_eas = (complex double *)_mm_malloc( input.n_nuclides * input.avg_n_poles * sizeof(complex double), 64);
  complex double *mp_rts = (complex double *)_mm_malloc( input.n_nuclides * input.avg_n_poles * sizeof(complex double), 64);
  complex double *mp_ras = (complex double *)_mm_malloc( input.n_nuclides * input.avg_n_poles * sizeof(complex double), 64);
  complex double *mp_rfs = (complex double *)_mm_malloc( input.n_nuclides * input.avg_n_poles * sizeof(complex double), 64);
  short int *l_vals = (short int *)_mm_malloc( input.n_nuclides * input.avg_n_poles * sizeof(short int), 64);

  MIC_Pole mpole;
  mpole.MP_EA = mp_eas;
  mpole.MP_RT = mp_rts;
  mpole.MP_RA = mp_ras;
  mpole.MP_RF = mp_rfs;
  mpole.l_value = l_vals;

  //printf("here\n");

  return mpole;
}

Pole ** generate_poles( Input input, int * n_poles )
{
	// Allocating 2D contiguous matrix
	Pole ** R = (Pole **)_mm_malloc( input.n_nuclides * sizeof( Pole *), 64);
	Pole * contiguous = (Pole *)_mm_malloc( input.n_nuclides * input.avg_n_poles * sizeof(Pole), 64);

	int k = 0;
	for( int i = 0; i < input.n_nuclides; i++ )
	{
		R[i] = &contiguous[k];
		k += n_poles[i];
	}
	
	// fill with data
	for( int i = 0; i < input.n_nuclides; i++ )
		for( int j = 0; j < n_poles[i]; j++ )
		{
			R[i][j].MP_EA = (double) rand() / RAND_MAX + (double) rand() / RAND_MAX * _Complex_I;
			R[i][j].MP_RT = (double) rand() / RAND_MAX + (double) rand() / RAND_MAX * _Complex_I;
			R[i][j].MP_RA = (double) rand() / RAND_MAX + (double) rand() / RAND_MAX * _Complex_I;
			R[i][j].MP_RF = (double) rand() / RAND_MAX + (double) rand() / RAND_MAX * _Complex_I;
			R[i][j].l_value = rand() % input.numL;
		}
	
	/* Debug
	for( int i = 0; i < input.n_nuclides; i++ )
		for( int j = 0; j < n_poles[i]; j++ )
			printf("R[%d][%d]: Eo = %lf lambda_o = %lf Tn = %lf Tg = %lf Tf = %lf\n", i, j, R[i][j].Eo, R[i][j].lambda_o, R[i][j].Tn, R[i][j].Tg, R[i][j].Tf);
	*/

	return R;
}

Window ** generate_window_params( Input input, int * n_windows, int * n_poles )
{
	// Allocating 2D contiguous matrix
	Window ** R = (Window **) malloc( input.n_nuclides * sizeof( Window *));
	Window * contiguous = (Window *) malloc( input.n_nuclides * input.avg_n_windows * sizeof(Window));

	int k = 0;
	for( int i = 0; i < input.n_nuclides; i++ )
	{
		R[i] = &contiguous[k];
		k += n_windows[i];
	}
	
	// fill with data
	for( int i = 0; i < input.n_nuclides; i++ )
	{
		int space = n_poles[i] / n_windows[i];
		int ctr = 0;
		for( int j = 0; j < n_windows[i]; j++ )
		{
			R[i][j].T = (double) rand() / RAND_MAX;
			R[i][j].A = (double) rand() / RAND_MAX;
			R[i][j].F = (double) rand() / RAND_MAX;
			R[i][j].start = ctr*j; 
			R[i][j].end = ctr*j + space;
			if( j == n_windows[i] - 1 )
				R[i][j].end = n_poles[i] - 1;
		}
	}

	return R;
}

double ** generate_pseudo_K0RS( Input input )
{
	double ** R = (double **) malloc( input.n_nuclides * sizeof( double * ));
	double * contiguous = (double *) malloc( input.n_nuclides * input.numL * sizeof(double));

	for( int i = 0; i < input.n_nuclides; i++ )
		R[i] = &contiguous[i*input.numL];

	return R;
}

void cleanup(int *n_poles, int *n_windows, Input input, CalcDataPtrs data) {
  free(n_poles);
  free(n_windows);
  _mm_free(data.poles);
  _mm_free(data.mypoles->MP_EA);
  _mm_free(data.mypoles->MP_RT);
  _mm_free(data.mypoles->MP_RA);
  _mm_free(data.mypoles->MP_RF);
  _mm_free(data.mypoles->l_value);
  free(data.materials.num_nucs);
  free(data.materials.mats);
  free(data.materials.concs);
  free(data.windows);
  free(data.pseudo_K0RS);
}
