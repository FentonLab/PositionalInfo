/*	prob_spk_map.c 							*/

/* 21-08-2015
	created make_pdf_spk_map_hd() which performs the analogous analysis on the information about head direction.
	TSD Binary data format:
        Xposition:1byte Yposition:1byte Direction:4bytes NspksInSample:1byte TimestampOfSpk-1:4bytes...TimeStampOfSpk-n:4bytes */

/* 08-08-2009
set the fread() and fwrite() parameters to read constant 1 and 4 byte data pieces, conforming to the TS format.
Previously these were machine-specific because they were set by sizeof() calls, which could generate errors on some machines */

#include	<stdio.h>
#include	<stdlib.h>
#include	<strings.h>
#include	<math.h>
#include 	<assert.h>

#define DEB 		0	

#define	TS_1_BYTE	1
#define	TS_4_BYTES	4

typedef  char str256[256];

char *get_line (FILE *fp, char *s);

void 	make_pdf_spk_map_hd (fp_tsd, fp_pdf, fp_fr, fp_DirInfo, fp_momentary_DirInfo, abs_info)
    	FILE *fp_tsd, *fp_pdf, *fp_fr, *fp_DirInfo, *fp_momentary_DirInfo;
	int	abs_info;
{
	extern	void
		*memchr();

	extern	size_t
		strlen();

	extern unsigned int 	
		max_hdsz, max_spsz,
	 	begin_analysis, end_analysis, MinHDSample, pdf,
		own_scale_hd, own_scale_x, own_scale_y, momentary;

	extern double 		
		eps, time_step, undefined_prob, undefined_info;

	int
		hd;
	
	unsigned int
		*nmb_visits,	/*  to check that no spikes are lost  */
		**hd_spk_nmb,
		i, j, k,
		x, y, scale_x, scale_y, spks, this_max_spsz,
		nmb_header_lines,
		nmb_visited_directions,
		hd_min, hd_max;


    	static unsigned int
		*spk_nmb,
		tsd_4_bytes = TS_4_BYTES,
		tsd_1_byte = TS_1_BYTE;

		
    	unsigned long 	
		i_time_step, l_dummy,
		nmb_time_steps, begin_analysis_time_step, end_analysis_time_step,
		spk_time_curr, data_start;

	
	double 	*spk_prob,
		*spk_prob_id,	/*  spike probabilities for the idealized case when each direction is equally sampled  */	 
		*hd_rate, 
		**hd_spk_prob,
		**hd_momentary_DirInfo,
		*I_spec_k,
		mean_rate,
		prob_D, prob_D_id, 
		P_k_d, d_dummy,
	 	dir_info, dir_info_avg,
	 	H_K, H_K_hd, I_spec_d;

    	char
		s[256], *ch_ptr;

    	unsigned char
		xc, yc, spksc, ch_dummy;
	

// ALLOCATE MEMORY
    	spk_nmb = (unsigned int *) calloc ( max_spsz,  (size_t) sizeof (unsigned int) ); 
    	assert (spk_nmb != NULL);

    	spk_prob = (double *) calloc ( max_spsz,  (size_t) sizeof (double) ); 
    	assert (spk_prob != NULL);

    	nmb_visits = (unsigned int *) calloc ( (size_t) max_hdsz,  (size_t) sizeof(unsigned int )); 
    	assert (nmb_visits != NULL);
    
    	hd_spk_nmb = (unsigned int **) calloc ( (size_t) max_hdsz,  (size_t) sizeof(unsigned int *)); 
    	assert (hd_spk_nmb != NULL);

    	hd_spk_prob = (double **) calloc ( (size_t) max_hdsz,  (size_t) sizeof(double *)); 
    	assert (hd_spk_prob != NULL);

    	hd_momentary_DirInfo = (double **) calloc ( (size_t) max_hdsz,  (size_t) sizeof(double *)); 
    	assert (hd_momentary_DirInfo != NULL);

    	hd_rate = (double *) calloc ( (size_t) max_hdsz,  (size_t) sizeof(double )); assert (hd_rate != NULL); 

	for (hd=0; hd<max_hdsz; hd++){
		hd_spk_nmb[hd] = (unsigned int *) calloc ( (size_t) max_spsz,  (size_t) sizeof(unsigned int ));
        	assert (hd_spk_nmb[hd] != NULL);

        	hd_spk_prob[hd] = (double *) calloc ( (size_t) max_spsz,  (size_t) sizeof(double));
        	assert (hd_spk_prob[hd] != NULL);

        	hd_momentary_DirInfo[hd] = (double *) calloc ( (size_t) max_spsz,  (size_t) sizeof(double));
        	assert (hd_momentary_DirInfo[hd] != NULL);
	}

// READ THE HEADER OF THE TSD FILE AND DETERMINE SPATIAL SCALES
    	fscanf(fp_tsd,"%d", &nmb_header_lines);
    	for (i = 0; i < nmb_header_lines ; i++) {
		if ((ch_ptr = (char *) get_line (fp_tsd, s)) == NULL) {
	    		fprintf (stderr, "Error in reading the header");
	    		exit (-1);
		}
#if (DEB)
		printf("\nline %d: %s", i, s);
#endif
		if ( (ch_ptr = (char *) memchr (s, '%', strlen(s))) == NULL) continue; 

		if ( memcmp (ch_ptr+1, "SAMPLING_INTERVAL(samps/sec)", strlen("SAMPLING_INTERVAL(samps/sec)")) == 0 ) {
	    		time_step = (double)atof(ch_ptr + strlen("%SAMPLING_INTERVAL(samps/sec)"));
			time_step = 1.0 / time_step;
	    		printf("\nTSD file time_step = %0.2lf ms", time_step * 1000.0);
		};
		if ( memcmp (ch_ptr+1, "SCALE_Y", 7) == 0 ) {
	    		if (own_scale_y) scale_y = own_scale_y;
	    		else scale_y =  atoi(ch_ptr);
	    		// printf("\n\nscale_y = %d", scale_y);
			continue;
		};
		if ( memcmp (ch_ptr+1, "SCALE_X", 7) == 0 ) {
	    		if (own_scale_x) scale_x = own_scale_x;
	    		else scale_x =  atoi(ch_ptr);
	    		// printf("\nscale_x = %d\n", scale_x);
		};
  
  	}

	if ( own_scale_hd == 0) {
		fprintf (stderr, "Error: directional scaling factor for must be defined by option -y\n");
		exit (-1);
	}

	begin_analysis_time_step = begin_analysis ? (unsigned long) (begin_analysis / time_step) : 0;
	end_analysis_time_step = end_analysis ? (unsigned long) (end_analysis / time_step) : 0;

	// printf ("\nbegin_analysis_time_step = %d", begin_analysis_time_step);
	// printf ("\nend_analysis_time_step = %d", end_analysis_time_step);

	nmb_time_steps = 0;
	i_time_step = 0;

	data_start = ftell(fp_tsd); // the spike at position data begin here.

// LOOP FOR CALCULATING HISTOGRAMS OF SPIKES
    	while (!feof(fp_tsd)) {
		fread (&xc, tsd_1_byte, 1, fp_tsd); if (feof(fp_tsd)) break; x = (unsigned int) xc;
		fread (&yc, tsd_1_byte, 1, fp_tsd); y = (unsigned int) yc;
		fread (&hd, tsd_4_bytes, 1, fp_tsd); 
		fread (&spksc, tsd_1_byte, 1, fp_tsd); spks = (unsigned int) spksc;

		if(spks > max_spsz){
			printf("\nError: spks (%d) > max_spsz (%d)\n", spks, max_spsz);
			exit(-1);
		}

//	Skipping the rest of the line 
		for (i=0; i<spks; i++) fread (&spk_time_curr, tsd_4_bytes, 1, fp_tsd);

		i_time_step ++;
		if (i_time_step >  begin_analysis_time_step){
			if (end_analysis_time_step != 0 && i_time_step > end_analysis_time_step) break;


			if (hd < 0 ){ 	// direction of the rat for this spk is unknown

			}else{
				hd /= own_scale_hd;
				hd_spk_nmb[hd][spks]++ ; 
				nmb_visits[hd]++ ;
				nmb_time_steps++;
			}
		}
	}


	nmb_visited_directions = nmb_time_steps = 0;
	mean_rate = 0.0;

// LOOP FOR CALCULATING AVERAGES FOR ALL DIRECTIONS
	for (hd=0; hd<max_hdsz; hd++){
		if (nmb_visits[hd] < MinHDSample) continue;
		nmb_time_steps += nmb_visits[hd];
		nmb_visited_directions++;
		for (k=0; k<max_spsz; k++) {
			spk_nmb[k] += hd_spk_nmb[hd][k];
			hd_rate[hd] += k * hd_spk_nmb[hd][k];
		}
		mean_rate += hd_rate[hd];
		hd_rate[hd] /= nmb_visits[hd];
	}

    	for (k=max_spsz; --k;) 	// find the bin with the largest number of spikes
		if(spk_nmb[k] > 0) break;
	this_max_spsz = k + 1;

	mean_rate /= nmb_time_steps;

      	H_K = 0.0;
	printf("\nspk_prob:\n");
	for (k=0; k<this_max_spsz; k++) {
		spk_prob[k] = (1.0 * spk_nmb[k]) / nmb_time_steps; 
		printf("%0.6f\t", spk_prob[k]);
                if (spk_prob[k] > eps) H_K -= spk_prob[k] * log (spk_prob [k]);
	}

	if (pdf) {
	    	fprintf (fp_pdf, "#	Contentsd:\n");
		fprintf (fp_pdf, "#1st line: max_hdsz, max_spsz,  nmb_time_steps\n");
		fprintf (fp_pdf, "#2nd line: spk_prob[0] ... spk_prob[max_spsz-1]\n");
		fprintf (fp_pdf, "#from the 3rd line: for (hd = 0; hd<max_hdsz; hd++)\n");
		fprintf (fp_pdf, "#nmb_visits[hd]\thd_spk_nmb[hd][0]\t ...\t hd_spk_nmb[hd][max_spsz-1] (times dir_visit[hd])\n"); 
    
		fprintf (fp_pdf, "%d\t%d\t%d\n", max_hdsz, max_spsz, nmb_time_steps);
	}

	printf("\nspk_nmbs:\n");


	l_dummy = 0;
    	for (k=0; k<this_max_spsz; k++) {
		printf("spk_nmb[%d]=%d\t", k, spk_nmb[k]);
		if (k) l_dummy += k * spk_nmb[k];
		if (pdf) {
			if (k < max_spsz - 1) fprintf (fp_pdf, "%d\t", spk_nmb[k]); 
			else fprintf (fp_pdf, "%d\n", spk_nmb[k]);
		} else {
			if (k < max_spsz - 1) printf ("%d\t", spk_nmb[k]); 
			else  printf ("%d\n", spk_nmb[k]); 
		}
	}
	printf("\n\ntotal number of spikes =\t%d", l_dummy);
	printf  ("\nmean_rate =             \t%f", mean_rate / time_step);
	printf  ("\nnmb_time_steps =        \t%d", nmb_time_steps);
	printf  ("\nnmb_visited_directions =    \t%d", nmb_visited_directions);

    	for (hd=0; hd<max_hdsz; hd++){
		if (nmb_visits[hd] >= MinHDSample) {
			prob_D = (1.0 * nmb_visits[hd])/nmb_time_steps;
		}	

		if (pdf) {
			fprintf (fp_pdf, "%d\t%d", hd,  nmb_visits[hd]);
			for (k=0; k<max_spsz; k++) fprintf (fp_pdf, "\t%d", hd_spk_nmb[hd][k]); 
			fprintf (fp_pdf, "\n"); 			
		}
    	}	

	fprintf(fp_DirInfo, "%s\t%s\t%s\t%s\n", "dir_info", "spec_D", "nmb_visits", "prob_D");


// LOOP FOR CALCULATING DIRECTION-SPECIFIC INFORMATIONS

    	for (hd=0; hd<max_hdsz; hd++){
		if ((hd < 0) || (nmb_visits[hd] == 0)) {
			fprintf (fp_fr, "%f\n", -1.0);
			fprintf (fp_DirInfo, "%.14f\t%.14f\t%.14f\t%.14f\n", -1.0,-1.0, -1.0, -1.0);
			continue;
		}
		if ( nmb_visits[hd] < MinHDSample) {
			fprintf (fp_fr, "%f\n", -2.0);
			fprintf (fp_DirInfo, "%.14f\t%.14f\t%.14f\t%.14f\n", -2.0, -2.0, -2.0, -2.0);
			continue;
		}

		prob_D = (1.0 * nmb_visits[hd])/nmb_time_steps;

		fprintf (fp_fr, "%f\n", hd_rate[hd] /time_step);

		dir_info = H_K_hd = 0.0;
		for (k = 0; k<this_max_spsz; k++) {
			P_k_d = (1.0 * hd_spk_nmb[hd][k]) / nmb_visits[hd]; 
			if (P_k_d > eps) {
			 	dir_info += P_k_d * log (P_k_d / spk_prob[k]);
				H_K_hd -= P_k_d * log (P_k_d);
			}
		}
		
		I_spec_d = (H_K - H_K_hd) * M_LOG2E; 
		dir_info *= M_LOG2E;
		fprintf (fp_DirInfo, "%.14f\t%.14f\t%.14f\t%.14f\n", \
			dir_info, I_spec_d, (double) nmb_visits[hd], prob_D);
		
     	} // END OF THE MAIN LOOP

// LOOP FOR CALCULATING MOMENTARY DIRECTION-SPECIFIC INFORMATIONS
	if (momentary) {
		// calculate direction-specific probabilities then info
		// direction-independent probabilities already calculated and stored in spk_prob[]
		// log2[p(k|d) / p(k)]

    		for (hd=0; hd<max_hdsz; hd++) {
			for (k = 0; k<this_max_spsz; k++){ 
				if(nmb_visits[hd] >= MinHDSample){
					hd_spk_prob[hd][k] = (double)hd_spk_nmb[hd][k]/ (double)nmb_visits[hd];
					// hd_momentary_DirInfo[hd][k] = log(hd_spk_prob[hd][k] / spk_prob[hd]);
					hd_momentary_DirInfo[hd][k] = hd_spk_prob[hd][k] * log(hd_spk_prob[hd][k] / spk_prob[k]);
					if(abs_info)
						hd_momentary_DirInfo[hd][k] = fabs(hd_momentary_DirInfo[hd][k]); 
					// hd_momentary_DirInfo[hd][k] = -1.0 * hd_spk_prob[hd][k] * log(hd_spk_prob[hd][k]) + spk_prob[k] * log(spk_prob[k]);
				}else{
					hd_spk_prob[hd][k] = undefined_prob;
					hd_momentary_DirInfo[hd][k] = undefined_info;
				}
			}
		}

		// now read in all the data and at each time_step look-up the momentary pos info

		nmb_time_steps = 0;
		i_time_step = 0;
		fseek(fp_tsd, data_start, SEEK_SET); 
    		while (!feof(fp_tsd)) {
			fread (&xc, tsd_1_byte, 1, fp_tsd); if (feof(fp_tsd)) break; x = (unsigned int) xc;
			fread (&yc, tsd_1_byte, 1, fp_tsd); y = (unsigned int) yc;
                	fread (&hd, tsd_4_bytes, 1, fp_tsd);
			fread (&spksc, tsd_1_byte, 1, fp_tsd); spks = (unsigned int) spksc;

		//	Skipping the rest of the line 
			for (i=0; i<spks; i++) fread (&spk_time_curr, tsd_4_bytes, 1, fp_tsd);

			i_time_step ++;

			if (i_time_step <=  begin_analysis_time_step) continue;
			if (end_analysis_time_step != 0 && i_time_step > end_analysis_time_step) break;

			assert ( spks < max_spsz);

			if(hd < 0){
				fprintf(fp_momentary_DirInfo, "%d\t%0.6lf\t%0.6lf\n",i_time_step, undefined_info, (double)(spks/time_step));
			}else{
				x /= own_scale_x;
				y /= own_scale_y; 
				hd /= own_scale_hd;

				if((hd_spk_prob[hd][spks] == undefined_prob) || (hd_momentary_DirInfo[hd][spks] == undefined_info)){
					fprintf(fp_momentary_DirInfo, "%d\t%0.6lf\t%0.6lf\n",i_time_step, undefined_info, (double)(spks/time_step));
				}else{
					fprintf(fp_momentary_DirInfo, "%d\t%0.6lf\t%0.6lf\n",i_time_step, hd_momentary_DirInfo[hd][spks], (double)(spks/time_step));
				}
			}
		}
	}

	printf  ("\n");

//	Free memory	

    	free(spk_prob);
    	free(spk_nmb);

    	for (hd=0; hd<max_hdsz; hd++) {
		free(hd_spk_nmb[hd]);
		free(hd_spk_prob[hd]);
		free(hd_momentary_DirInfo[hd]);
    	}
  
    	free (hd_spk_nmb);
    	free (hd_spk_prob);
    	free (hd_momentary_DirInfo);
	free (hd_rate);
    	free (nmb_visits);
}
