/*	prob_spk_map.c 							*/

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

void 	make_pdf_spk_map_xy (fp_ts, fp_pdf, fp_fr, fp_pinfo, fp_momentary_pinfo, abs_info)
    	FILE *fp_ts, *fp_pdf, *fp_fr, *fp_pinfo, *fp_momentary_pinfo;
	int	abs_info;

{
	extern	void
		*memchr();

	extern	size_t
		strlen();

	extern unsigned int 	
		skaggs, max_xsz, max_ysz, max_spsz,
	 	begin_analysis, end_analysis, n_cut, pdf,
		own_scale_x, own_scale_y, momentary;

	extern double 		
		eps, time_step, undefined_prob, undefined_info;
	
	unsigned int
		**nmb_visits,	/*  to check that no spikes are lost  */
		***loc_spk_nmb,
		scale_x = 0, scale_y = 0, 
		i, j, k,
		x, y, spks, this_max_spsz,
		nmb_header_lines,
		nmb_visited_pixels,
		x_min, x_max, y_min, y_max;


    	static unsigned int
		*spk_nmb,
		ts_4_bytes = TS_4_BYTES,
		ts_1_byte = TS_1_BYTE;

		
    	unsigned long 	
		i_time_step, l_dummy,
		nmb_time_steps, begin_analysis_time_step, end_analysis_time_step,
		spk_time_curr, data_start;

	
	double 	*spk_prob,
		*spk_prob_id,	/*  spike probabilities for the idealized case when each pixel is equally sampled  */	 
		**loc_rate, 
		***loc_spk_prob,
		***loc_momentary_pinfo,
		*I_spec_k,
		mean_rate,
		prob_X, prob_X_id, 
		P_k_x, skaggs_sum, d_dummy,
	 	pos_info, pos_info_avg,
	 	H_K, H_K_loc, I_spec_x;

    	char
		s[256], *ch_ptr;

    	unsigned char
		xc, yc, spksc, ch_dummy;
	
// ALLOCATE MEMORY
    	spk_nmb = (unsigned int *) calloc ( max_spsz,  (size_t) sizeof (unsigned int) ); 
    	assert (spk_nmb != NULL);
    	spk_prob = (double *) calloc ( max_spsz,  (size_t) sizeof (double) ); 
    	assert (spk_prob != NULL);

    	nmb_visits = (unsigned int **) calloc ( (size_t) max_ysz,  (size_t) sizeof(unsigned int *)); 
    	assert (nmb_visits != NULL);
    
    	loc_spk_nmb = (unsigned int ***) calloc ( (size_t) max_ysz,  (size_t) sizeof(unsigned int **)); 
    	assert (loc_spk_nmb != NULL);

    	loc_spk_prob = (double ***) calloc ( (size_t) max_ysz,  (size_t) sizeof(double **)); 
    	assert (loc_spk_prob != NULL);

    	loc_momentary_pinfo = (double ***) calloc ( (size_t) max_ysz,  (size_t) sizeof(double **)); 
    	assert (loc_momentary_pinfo != NULL);

    	loc_rate = (double **) calloc ( (size_t) max_ysz,  (size_t) sizeof(double *)); assert (loc_rate != NULL); 

    	for (y=0; y<max_ysz; y++)  {
		nmb_visits[y] = (unsigned int *) calloc ( (size_t) max_xsz,  (size_t) sizeof(unsigned int)); 
		assert (nmb_visits[y] != NULL);

		loc_spk_nmb[y] = (unsigned int **) calloc ( (size_t) max_xsz,  (size_t) sizeof(unsigned int *)); 
		assert (loc_spk_nmb[y] != NULL);

		loc_spk_prob[y] = (double **) calloc ( (size_t) max_xsz,  (size_t) sizeof(double *)); 
		assert (loc_spk_prob[y] != NULL);

		loc_momentary_pinfo[y] = (double **) calloc ( (size_t) max_xsz,  (size_t) sizeof(double *)); 
		assert (loc_momentary_pinfo[y] != NULL);

		loc_rate[y] = (double*) calloc ( (size_t) max_xsz,  (size_t) sizeof(double)); assert (loc_rate[y] != NULL); 

		for (x=0; x<max_xsz; x++) {
	    		loc_spk_nmb[y][x] = (unsigned int *) calloc ( (size_t) max_spsz,  (size_t) sizeof(unsigned int )); 
	    		assert (loc_spk_nmb[y][x] != NULL);

	    		loc_spk_prob[y][x] = (double *) calloc ( (size_t) max_spsz,  (size_t) sizeof(double)); 
	    		assert (loc_spk_prob[y][x] != NULL);

	    		loc_momentary_pinfo[y][x] = (double *) calloc ( (size_t) max_spsz,  (size_t) sizeof(double)); 
	    		assert (loc_momentary_pinfo[y][x] != NULL);
		}
    	}	

// READ THE HEADER OF THE TS FILE AND DETERMINE SPATIAL SCALES
    	fscanf(fp_ts,"%d", &nmb_header_lines);
    	for (i = 0; i < nmb_header_lines ; i++) {
		if ((ch_ptr = (char *) get_line (fp_ts, s)) == NULL) {
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
	    		printf("\nTS file time_step = %0.2lf ms", time_step * 1000.0);
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

	if ( own_scale_y ==0 && scale_y ==0) {
			fprintf (stderr, "Error: spatial scaling factor for y must be defined by options -V\n");
			exit (-1);
	}
	if ( own_scale_x ==0 && scale_x ==0) {
			fprintf (stderr, "Error: spatial scaling factor for x must be defined by options -H\n");
			exit (-1);
	}


	begin_analysis_time_step = begin_analysis ? (unsigned long) (begin_analysis / time_step) : 0;
	end_analysis_time_step = end_analysis ? (unsigned long) (end_analysis / time_step) : 0;

	// printf ("\nbegin_analysis_time_step = %d", begin_analysis_time_step);
	// printf ("\nend_analysis_time_step = %d", end_analysis_time_step);

	nmb_time_steps = 0;
	i_time_step = 0;

	data_start = ftell(fp_ts); // the spike at position data begin here.

// LOOP FOR CALCULATING LOCAL HISTOGRAMS OF SPIKES
    	while (!feof(fp_ts)) {
		fread (&xc, ts_1_byte, 1, fp_ts); if (feof(fp_ts)) break; x = (unsigned int) xc;
		fread (&yc, ts_1_byte, 1, fp_ts); y = (unsigned int) yc;
		fread (&spksc, ts_1_byte, 1, fp_ts); spks = (unsigned int) spksc;

		if(spks > max_spsz){
			printf("\nError: spks (%d) > max_spsz (%d)\n", spks, max_spsz);
			exit(-1);
		}

//	Skipping the rest of the line 
		for (i=0; i<spks; i++) fread (&spk_time_curr, ts_4_bytes, 1, fp_ts);

		i_time_step ++;
		if (i_time_step <=  begin_analysis_time_step) goto end_while;
		if (end_analysis_time_step != 0 && i_time_step > end_analysis_time_step) break;


		if (x == 0 && y == 0 ) goto end_while;	// location of the rat for this spk is unknown

		x /= scale_x;
		y /= scale_y; 
		if ((x >= max_xsz) || (y >= max_ysz)){
			printf("\nPosition coordinates (%d, %d) out of (max_x, max_y) range (%d, %d) use larger value of scale.\n", x, y, max_xsz, max_ysz);
			exit(-1);
		}

		loc_spk_nmb [y][x][spks]++ ; 
		nmb_visits [y][x]++ ;
		nmb_time_steps++;
		end_while: ;
	}


	nmb_visited_pixels = nmb_time_steps = 0;
	mean_rate = 0.0;

// LOOP FOR CALCULATING AVERAGES FOR ALL LOCATIONS
	for (y=0; y<max_ysz; y++) for (x=0; x<max_xsz; x++) {
		if (((y == 0) && (x == 0)) || (nmb_visits[y][x] < n_cut)) continue;
		nmb_time_steps += nmb_visits[y][x];
		nmb_visited_pixels ++;
		for (k=0; k<max_spsz; k++) {
			spk_nmb[k] += loc_spk_nmb[y][x][k];
			loc_rate[y][x] += k * loc_spk_nmb[y][x][k];
		}
		mean_rate += loc_rate[y][x];
		loc_rate[y][x] /= nmb_visits[y][x];
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
	    	fprintf (fp_pdf, "#	Contents:\n");
		fprintf (fp_pdf, "#1st line: max_ysz, max_xsz, max_spsz,  nmb_time_steps\n");
		fprintf (fp_pdf, "#2nd line: spk_prob[0] ... spk_prob[max_spsz-1]\n");
		fprintf (fp_pdf, "#from the 3rd line: for (y = 0; y<max_ysz; y++) for (x=0; x<max_xsz; x++)\n");
		fprintf (fp_pdf, "#y\tx\tnmb_visits[y][x]\tloc_spk_nmb[y][x][0]\t ...\t loc_spk_nmb[y][x][max_spsz-1] (times nmb_visit[y][x])\n"); 
    
		fprintf (fp_pdf, "%d\t%d\t%d\t%d\n", max_ysz, max_xsz, max_spsz, nmb_time_steps);
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
	printf  ("\nnmb_visited_pixels =    \t%d", nmb_visited_pixels);

    	for (y=0; y<max_ysz; y++) for (x=0; x<max_xsz; x++) {
		if (nmb_visits[y][x] >= n_cut) {
			prob_X = (1.0 * nmb_visits[y][x])/nmb_time_steps;
		}	
		if (pdf) {
			fprintf (fp_pdf, "%d\t%d\t%d", y, x,  nmb_visits[y][x]);
			for (k=0; k<max_spsz; k++) fprintf (fp_pdf, "\t%d", loc_spk_nmb[y][x][k]); 
			fprintf (fp_pdf, "\n"); 			
		}
    	}	

	fprintf(fp_pinfo, "%s\t%s\t%s\t%s\n", "pos_info", "spec_x", "nmb_visits", "prob_x");

	if (skaggs) { skaggs_sum = 0.0; pos_info_avg = 0.0; }

// LOOP FOR CALCULATING LOCATION-SPECIFIC INFORMATIONS

    	for (y=0; y<max_ysz; y++) for (x=0; x<max_xsz; x++) {
		if (((y == 0) && (x == 0)) || nmb_visits[y][x] == 0 ) {
			fprintf (fp_fr, "%f\n", -1.0);
			fprintf (fp_pinfo, "%.14f\t%.14f\t%.14f\t%.14f\n", -1.0,-1.0, -1.0, -1.0);
			continue;
		}
		if ( nmb_visits[y][x] < n_cut) {
			fprintf (fp_fr, "%f\n", -2.0);
			fprintf (fp_pinfo, "%.14f\t%.14f\t%.14f\t%.14f\n", -2.0, -2.0, -2.0, -2.0);
			continue;
		}

		prob_X = (1.0 * nmb_visits[y][x])/nmb_time_steps;

		fprintf (fp_fr, "%f\n", loc_rate [y][x] /time_step);
		if (skaggs && loc_rate[y][x] > 0.0 ) {
			skaggs_sum += loc_rate[y][x] * log (loc_rate[y][x] / mean_rate) * nmb_visits[y][x];
		}
	

		pos_info = H_K_loc = 0.0;
		for (k = 0; k<this_max_spsz; k++) {
			P_k_x = (1.0 * loc_spk_nmb[y][x][k]) / nmb_visits[y][x]; 
			if (P_k_x > eps) {
			 	pos_info += P_k_x * log (P_k_x / spk_prob[k]);
				H_K_loc -= P_k_x * log (P_k_x);				}
			}
		
		I_spec_x = (H_K - H_K_loc) * M_LOG2E; 
		pos_info *= M_LOG2E;
		fprintf (fp_pinfo, "%.14f\t%.14f\t%.14f\t%.14f\n", \
			pos_info, I_spec_x, (double) nmb_visits[y][x], prob_X);
		
		if (skaggs) pos_info_avg += pos_info * nmb_visits[y][x];

     	} // END OF THE MAIN LOOP

	if (skaggs) {
		skaggs_sum = skaggs_sum * M_LOG2E / mean_rate / nmb_time_steps;
		pos_info_avg = pos_info_avg / nmb_time_steps;
		printf ("\nSkaggs' info content =  \t%f", skaggs_sum);
		printf ("\nShannon's mutual info = \t%f\n", pos_info_avg);
	}


// LOOP FOR CALCULATING MOMENTARY LOCATION-SPECIFIC INFORMATIONS
	if (momentary) {
		// calculate location-specific probabilities then info
		// location-independent probabilities already calculated and stored in spk_prob[]
		// log2[p(k|x) / p(k)]

    		for (y=0; y<max_ysz; y++) for (x=0; x<max_xsz; x++) {
			for (k = 0; k<this_max_spsz; k++){ 
				if(nmb_visits[y][x] >= n_cut){
					loc_spk_prob[y][x][k] = (double)loc_spk_nmb[y][x][k]/ (double)nmb_visits[y][x];
					// loc_momentary_pinfo[y][x][k] = log(loc_spk_prob[y][x][k] / spk_prob[k]);
					loc_momentary_pinfo[y][x][k] = loc_spk_prob[y][x][k] * log(loc_spk_prob[y][x][k] / spk_prob[k]);
					if(abs_info)
						loc_momentary_pinfo[y][x][k] = fabs(loc_momentary_pinfo[y][x][k]); 
					// loc_momentary_pinfo[y][x][k] = -1.0 * loc_spk_prob[y][x][k] * log(loc_spk_prob[y][x][k]) + spk_prob[k] * log(spk_prob[k]);
				}else{
					loc_spk_prob[y][x][k] = undefined_prob;
					loc_momentary_pinfo[y][x][k] = undefined_info;
				}
			}
		}

		// now read in all the data and at each time_step look-up the momentary pos info

		nmb_time_steps = 0;
		i_time_step = 0;
		fseek(fp_ts, data_start, SEEK_SET); 
    		while (!feof(fp_ts)) {
			fread (&xc, ts_1_byte, 1, fp_ts); if (feof(fp_ts)) break; x = (unsigned int) xc;
			fread (&yc, ts_1_byte, 1, fp_ts); y = (unsigned int) yc;
			fread (&spksc, ts_1_byte, 1, fp_ts); spks = (unsigned int) spksc;

		//	Skipping the rest of the line 
			for (i=0; i<spks; i++) fread (&spk_time_curr, ts_4_bytes, 1, fp_ts);

			i_time_step ++;

			if (i_time_step <=  begin_analysis_time_step) continue;
			if (end_analysis_time_step != 0 && i_time_step > end_analysis_time_step) break;

			assert ( spks < max_spsz);

			x /= scale_x;
			y /= scale_y; 
			
			if((loc_spk_prob[y][x][spks] == undefined_prob) || (loc_momentary_pinfo[y][x][spks] == undefined_info))
				fprintf(fp_momentary_pinfo, "%d\t%0.6lf\t%0.6lf\n",i_time_step, undefined_info, (double)(spks/time_step));
			else
				fprintf(fp_momentary_pinfo, "%d\t%0.6lf\t%0.6lf\n",i_time_step, loc_momentary_pinfo[y][x][spks], (double)(spks/time_step));
		}
	}

//	Free memory	

    	free(spk_prob);
    	free(spk_nmb);

    	for (y=0; y<max_ysz; y++) {
		for (x=0; x<max_xsz; x++){
			free(loc_spk_nmb[y][x]);
			free(loc_spk_prob[y][x]);
			free(loc_momentary_pinfo[y][x]);
		}
		free (loc_spk_nmb[y]);
		free(loc_spk_prob[y]);
		free(loc_momentary_pinfo[y]);
		free (loc_rate[y]); 
		free (nmb_visits[y]);
    	}
  
    	free (loc_spk_nmb);
    	free (loc_spk_prob);
    	free (loc_momentary_pinfo);
	free (loc_rate);
    	free (nmb_visits);
}

char *get_line (FILE *fp, char *s)
{
    int i=0;

    while (fscanf(fp, "%c", s+i) != EOF) {
        if (*(s+i) == '\n') {
            *(s+i+1) = '\0';
            return s;
        }
        i++;
    }
    return NULL;
}
