#include "rsbench.h"

int main(int argc, char * argv[])
{
	// =====================================================================
	// Initialization & Command Line Read-In
	// =====================================================================

	int version = 2;
	int max_procs = omp_get_num_procs();
	double start, stop;

	srand(time(NULL));
	
	// Process CLI Fields
	Input input = read_CLI( argc, argv );

	// Set number of OpenMP Threads
	omp_set_num_threads(input.nthreads); 
	
	// =====================================================================
	// Print-out of Input Summary
	// =====================================================================
	logo(version);
	center_print("INPUT SUMMARY", 79);
	border_print();
	print_input_summary(input);

	// =====================================================================
	// Prepare Pole Paremeter Grids
	// =====================================================================
	border_print();
	center_print("INITIALIZATION", 79);
	border_print();
	
	start = omp_get_wtime();

  unsigned long long vhash = 0;

	// Allocate & fill energy grids
	printf("Generating resonance distributions...\n");
	int * n_poles = generate_n_poles( input );

	// Allocate & fill Window grids
	printf("Generating window distributions...\n");
	int * n_windows = generate_n_windows( input );

	// Get material data
	printf("Loading Hoogenboom-Martin material data...\n");
	Materials materials = get_materials( input ); 

	// Prepare full resonance grid
	printf("Generating resonance parameter grid...\n");
	Pole ** poles = generate_poles( input, n_poles );

	MIC_Pole mypoles = generate_mypoles( input, n_poles );

	// Prepare full Window grid
	printf("Generating window parameter grid...\n");
	Window ** windows = generate_window_params( input, n_windows, n_poles);

	// Prepare 0K Resonances
	printf("Generating 0K l_value data...\n");
	double ** pseudo_K0RS = generate_pseudo_K0RS( input );

	CalcDataPtrs data;
	data.n_poles = n_poles;
	data.n_windows = n_windows;
	data.materials = materials;
	data.poles = poles;
	data.mypoles = &mypoles;
	data.windows = windows;
	data.pseudo_K0RS = pseudo_K0RS;

	stop = omp_get_wtime();
	printf("Initialization Complete. (%.2lf seconds)\n", stop-start);
	
	// =====================================================================
	// Cross Section (XS) Parallel Lookup Simulation Begins
	// =====================================================================
	border_print();
	center_print("SIMULATION", 79);
	border_print();
	
	printf("Beginning Simulation.\n");
	#ifndef STATUS
	printf("Calculating XS's...\n");
	#endif
	
	#ifdef PAPI
	/* initialize papi with one thread here  */
	if ( PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT){
		fprintf(stderr, "PAPI library init error!\n");
		exit(1);
	}
	#endif	

	start = omp_get_wtime();

	unsigned long seed = rand();
	int mat;
	double E;
	int i, j;

	#pragma omp parallel default(none) \
	private(seed, mat, E, i, j) \
	shared(input, data, vhash) 
	{
		double macro_xs[4];
		int thread = omp_get_thread_num();
		seed += thread;
		
		#ifdef PAPI
		int eventset = PAPI_NULL; 
		int num_papi_events;
		#pragma omp critical
		{
			counter_init(&eventset, &num_papi_events);
		}
		#endif
		complex double * sigTfactors =
			(complex double *) malloc( input.numL * sizeof(complex double) );

    #pragma omp for schedule(dynamic) private(macro_xs)
		for( i = 0; i < input.lookups; i++ )
		{
			#ifdef STATUS
			if( thread == 0 && i % 1000 == 0 )
				printf("\rCalculating XS's... (%.0lf%% completed)",
						(i / ( (double)input.lookups /
						(double) input.nthreads )) /
						(double) input.nthreads * 100.0);
			#endif
			mat = pick_mat( &seed );
			E = rn( &seed );
			//calculate_macro_xs( macro_xs, mat, E, input, data, sigTfactors ); 
	    // zero out macro vector
	    for( int i = 0; i < 4; i++ )
	    	macro_xs[i] = 0;

      //TODO: notice how "novector" sucks!
      //#pragma novector
      #pragma simd
	    for( int i = 0; i < (data.materials).num_nucs[mat]; i++ )
	    //for( int i = 0; i < 16; i++ )
	    {
	    	double micro_xs[4];
	    	int nuc = (data.materials).mats[mat][i];
        //int nuc = i;

	    	//calculate_micro_xs( micro_xs, nuc, E, input, data, sigTfactors);
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
	      //calculate_sig_T(nuc, E, input, data, sigTfactors );
	      double phi;

        //TODO: Notice that this loop is a vectorization killer!
	      //for( int i = 0; i < input.numL; i++ )
	      //{
	      //	phi = data.pseudo_K0RS[nuc][i] * sqrt(E);

	      //	if( i == 1 )
	      //		phi -= - atan( phi );
	      //	else if( i == 2 )
	      //		phi -= atan( 3.0 * phi / (3.0 - phi*phi));
	      //	else if( i == 3 )
	      //		phi -= atan(phi*(15.0-phi*phi)/(15.0-6.0*phi*phi));

	      //	phi *= 2.0;

	      //	sigTfactors[i] = cos(phi) - sin(phi) * _Complex_I;
	      //}

	      phi = data.pseudo_K0RS[nuc][0] * sqrt(E);
	      phi *= 2.0;
	      sigTfactors[0] = cos(phi) - sin(phi) * _Complex_I;

	      phi = data.pseudo_K0RS[nuc][1] * sqrt(E);
	      phi -= - atan( phi );
	      phi *= 2.0;
	      sigTfactors[1] = cos(phi) - sin(phi) * _Complex_I;

	      phi = data.pseudo_K0RS[nuc][2] * sqrt(E);
	      phi -= atan( 3.0 * phi / (3.0 - phi*phi));
	      phi *= 2.0;
	      sigTfactors[2] = cos(phi) - sin(phi) * _Complex_I;

	      phi = data.pseudo_K0RS[nuc][3] * sqrt(E);
	      phi -= atan(phi*(15.0-phi*phi)/(15.0-6.0*phi*phi));
	      phi *= 2.0;
	      sigTfactors[3] = cos(phi) - sin(phi) * _Complex_I;

	      // Calculate contributions from window "background" (i.e., poles outside window (pre-calculated)
	      Window w = data.windows[nuc][window];
	      sigT = E * w.T;
	      sigA = E * w.A;
	      sigF = E * w.F;

        //int idx = data.n_poles[nuc];
        int idx = input.avg_n_poles*data.n_poles[nuc]-1+w.start;
        //int idx = data.n_poles[nuc]-1+w.start;

        // TODO: This range must be static for astounding performance, possible...?
	      // Loop over Poles within window, add contributions
//#pragma novector
#pragma ivdep
	      for( int i = 0; i < 4; i++ )
	      //for( int i = w.start; i < w.start+4; i++ )
	      //for( int i = w.start; i < w.end; i++ )
	      {
	      	complex double PSIIKI;
	      	complex double CDUM;
          //Pole pole;
          //pole.MP_EA = 0.4;
	        //pole.MP_RT = 0.3;
	        //pole.MP_RA = 0.2;
	        //pole.MP_RF = 0.1;
	        //pole.l_value = 123;

	      	//Pole pole = data.poles[nuc][i];
	      	//Pole pole = data.mypoles[nuc];
          //Pole pole;
          //pole.MP_EA = data.poles[nuc][i].MP_EA;
          //pole.MP_RT = data.poles[nuc][i].MP_RT;
          //pole.MP_RA = data.poles[nuc][i].MP_RA;
          //pole.MP_RF = data.poles[nuc][i].MP_RF;
          //pole.l_value = data.poles[nuc][i].l_value;

	      	//PSIIKI = -(0.0 - 1.0 * _Complex_I ) / ( pole.MP_EA - sqrt(E) );
	      	//CDUM = PSIIKI / E;
	      	//sigT += creal( pole.MP_RT * CDUM * sigTfactors[pole.l_value] );
	      	//sigA += creal( pole.MP_RA * CDUM);
	      	//sigF += creal( pole.MP_RF * CDUM);

          //complex double p1 = data.poles[nuc][i].MP_EA;
          //complex double p2 = data.poles[nuc][i].MP_RT;
          //complex double p3 = data.poles[nuc][i].MP_RA;
          //complex double p4 = data.poles[nuc][i].MP_RF;
          //complex double sT = sigTfactors[data.poles[nuc][i].l_value];
	  	    //PSIIKI = -(0.0 - 1.0 * _Complex_I ) / ( p1 - sqrt(E) );
	  	    //CDUM = PSIIKI / E;
	  	    //sigT += creal( p2 * CDUM * sT );
	  	    //sigA += creal( p3 * CDUM);
	  	    //sigF += creal( p4 * CDUM);

	  	    //PSIIKI = -(0.0 - 1.0 * _Complex_I ) / ( data.poles[nuc][i].MP_EA - sqrt(E) );
	  	    //CDUM = PSIIKI / E;
	  	    //sigT += creal( data.poles[nuc][i].MP_RT * CDUM * sigTfactors[data.poles[nuc][i].l_value] );
	  	    //sigA += creal( data.poles[nuc][i].MP_RA * CDUM);
	  	    //sigF += creal( data.poles[nuc][i].MP_RF * CDUM);

          // TODO: Something like this... 1st one seg faults on host
          //int idx = input.avg_n_poles*data.n_poles[nuc] + i;
          
          idx += i;

          //int idx = data.n_poles[nuc]+i;

	  	    PSIIKI = -(0.0 - 1.0 * _Complex_I ) / ( data.mypoles->MP_EA[idx] - sqrt(E) );
	  	    CDUM = PSIIKI / E;
	  	    sigT += creal( data.mypoles->MP_RT[idx] * CDUM * sigTfactors[data.mypoles->l_value[idx]] );
	  	    sigA += creal( data.mypoles->MP_RA[idx] * CDUM);
	  	    sigF += creal( data.mypoles->MP_RF[idx] * CDUM);
	      }

	      sigE = sigT - sigA;

	      micro_xs[0] = sigT;
	      micro_xs[1] = sigA;
	      micro_xs[2] = sigF;
	      micro_xs[3] = sigE;

	    	for( int j = 0; j < 4; j++ ) {
	    		macro_xs[j] += micro_xs[j] * 0.8;
	    		//macro_xs[j] += micro_xs[j] * data.materials.concs[mat][i];
	    	}
	    }

      //#pragma omp parallel shared(sum) num_threads(inner)
      //{

        //#pragma omp for reduction(+:sum)
        //#pragma novector
        //for (j=1; j<=input.avg_n_poles; j++) {
        //  sum += j;
        //}

      //}

      //printf("sum = %f\n", sum);
		}

		free(sigTfactors);
		
		#ifdef PAPI
		if( thread == 0 )
		{
			printf("\n");
			border_print();
			center_print("PAPI COUNTER RESULTS", 79);
			border_print();
			printf("Count          \tSmybol      \tDescription\n");
		}
		{
		#pragma omp barrier
		}
		counter_stop(&eventset, num_papi_events);
		#endif
    
#ifdef VERIFY
    //printf("macro_xs= %f\n", macro_xs[0]);
    char line[256];
    sprintf(line, "%.5lf %.5lf %.5lf %.5lf",
    macro_xs[0],
    macro_xs[1],
    macro_xs[2],
    macro_xs[3]);
    unsigned long long vhash_local = hash(line, 10000);

    #pragma omp atomic
    vhash += vhash_local;
#endif

	}

	stop = omp_get_wtime();
	#ifndef PAPI
	printf("\nSimulation Complete.\n");
	#endif
  
	// =====================================================================
	// Print / Save Results and Exit
	// =====================================================================
	border_print();
	center_print("RESULTS", 79);
	border_print();

	printf("Threads:     %d\n", input.nthreads);
	printf("Runtime:     %.3lf seconds\n", stop-start);
	printf("Lookups:     "); fancy_int(input.lookups);
	printf("Lookups/s:   "); fancy_int((double) input.lookups / (stop-start));

	border_print();

#ifdef VERIFY
  printf("Verification checksum: %llu\n", vhash);
#endif


  cleanup(n_poles, n_windows, input, data);

	return 0;
}
