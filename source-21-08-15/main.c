#include	<stdio.h>
#include	<stdlib.h>
#include	<strings.h>
#include	<math.h>
#include 	<assert.h>

#define MAX_XSZ 		32
#define MAX_YSZ 		32
#define MAX_HDSZ 		360
#define MAX_SPSZ 		40
#define N_CUT			10
#define EPS	 		1E-10
#define UNDEFINED_PROB		-100.0
#define UNDEFINED_INFO		-100.0
#define TIME_STEP		116.666667

#define DEB	0

typedef  char str256[256];

void instruct();
void make_pdf_spk_map (FILE *, FILE *, FILE *, FILE *, FILE *, unsigned int);
void make_pdf_spk_map_hd (FILE *, FILE *, FILE *, FILE *, FILE *, unsigned int);
void    check_env (char *name);

// char *strcpy(char *, char *);

unsigned int 	skaggs 	= 1,
		pdf	= 1,
		momentary = 1,
		max_xsz = MAX_XSZ,
		max_ysz = MAX_YSZ,
		max_hdsz = MAX_HDSZ,
		max_spsz = MAX_SPSZ,
		n_cut 	 = N_CUT,
		MinHDSample = N_CUT;

unsigned int	begin_analysis = 0, end_analysis = 0,
		own_scale_x = 1, 
		own_scale_y = 1,
		own_scale_hd = 1,
		abs_info = 0,
        	dir_info = 0;

double		eps = EPS,
		time_step = TIME_STEP,  
		undefined_prob = UNDEFINED_PROB,
		undefined_info = UNDEFINED_INFO;
    
int main(int argc, char **argv)
{
	FILE	
	    *fp, *fp_ts, *fp_pdf, *fp_fr, *fp_info, *fp_momentary_info;
	

	extern unsigned int optind;
	extern char *optarg;


	str256	
	    	file_list, file_name, 
		ts_file, pdf_file, fr_file, pinfo_file, momentary_info_file, 
		ts_file_dir, pdf_dir, fr_dir, info_dir, 
		str_dummy; 
	
	char	
	    ch_dummy;    

	if (argc == 1) instruct();

	while ((ch_dummy = getopt(argc,argv,"AS:B:Dd:E:MN:PK:X:Y:T:x:y:")) != EOF) {
	    	switch(ch_dummy) {
			case 'A': abs_info = 1;  break;
			case 'B': begin_analysis = atoi (optarg);  break;
            		case 'D': dir_info = 1;  break;
			case 'd': own_scale_hd = atoi (optarg);  break;
			case 'E': end_analysis = atoi (optarg); break;
			case 'H': max_hdsz = atoi (optarg);  break;
			case 'M': momentary = 0; break;
			case 'N': n_cut = atoi (optarg);  break;
			case 'P': pdf = 0; break;
			case 'S': max_spsz = atoi (optarg);  break;
			case 'T': time_step = atof (optarg); break;
			case 'X': max_xsz = atoi (optarg);  break;
			case 'x': own_scale_x = atoi (optarg);  break;
			case 'Y': max_ysz = atoi (optarg);  break;
			case 'y': own_scale_y = atoi (optarg);  break;
		}
	}

	if(dir_info){
		check_env ("TSD_INFO_DIR");  strcpy(ts_file_dir, getenv("TSD_INFO_DIR"));
		check_env ("DINFO_DIR");  strcpy(info_dir, getenv("DINFO_DIR"));
	}else{
		check_env ("TS_INFO_DIR");  strcpy(ts_file_dir, getenv("TS_INFO_DIR"));
		check_env ("PINFO_DIR");  strcpy(info_dir, getenv("PINFO_DIR"));
	}
	check_env ("PDF_DIR"); strcpy(pdf_dir, getenv("PDF_DIR"));
	check_env ("RATE_INFO_DIR");  strcpy(fr_dir, getenv("RATE_INFO_DIR"));


    strcpy (file_list, argv[argc-1]);	
    printf ("\nFILE_LIST = %s", file_list);
    
    fp = fopen(file_list,"r"); assert (fp != NULL);
    
    while ( fscanf(fp, "%s", str_dummy) != EOF) 
    {
	printf ("\n\nfile_name = %s", str_dummy);
	strcpy (file_name, str_dummy);
	sprintf (ts_file, "%s/%s", ts_file_dir, file_name);
        fp_ts = fopen(ts_file,"r"); 	//assert (fp_ts != NULL);
	if (fp_ts == NULL) {
		fprintf (stderr, "Can not open %s for reading!\n", ts_file);
		continue;
	}
//	printf ("\nts_file = %s", ts_file);

	if ( (begin_analysis + end_analysis) == 0) { 
		sprintf (pdf_file, "%s/%s%s", pdf_dir, file_name, ".pdf");
		sprintf (fr_file, "%s/%s%s", fr_dir, file_name, ".fr");
		sprintf (pinfo_file, "%s/%s%s", info_dir, file_name, ".pinfo");
		if(momentary)
			sprintf (momentary_info_file, "%s/%s%s", info_dir, file_name, ".mpinfo");
	} else {
		sprintf (pdf_file, "%s/%s%s%d%s%d", pdf_dir, file_name, ".pdf.", begin_analysis, "_", end_analysis);
		sprintf (fr_file, "%s/%s%s%d%s%d", fr_dir, file_name, ".fr.", begin_analysis, "_", end_analysis);
		sprintf (pinfo_file, "%s/%s%s%d%s%d", info_dir, file_name, ".pinfo.", begin_analysis, "_", end_analysis);
		if(momentary)
			sprintf (momentary_info_file, "%s/%s%s%d%s%d", info_dir, file_name, ".mpinfo.", begin_analysis, "_", end_analysis);
	}

/*
	printf ("\npdf_file = %s", pdf_file);
	printf ("\nfr_file = %s", fr_file);
	printf ("\npinfo_file = %s", pinfo_file);
	if(momentary)
		printf ("\nmomentary_info_file = %s", momentary_info_file);
*/

        fp_pdf = fopen(pdf_file,"w"); assert (fp_pdf != NULL);
        fp_fr = fopen(fr_file,"w"); assert (fp_fr != NULL);
        fp_info = fopen(pinfo_file,"w"); assert (fp_info != NULL);
        fp_momentary_info = fopen(momentary_info_file,"w"); assert (fp_momentary_info != NULL);
        
        // condition for head direction or location
        if(dir_info){
            make_pdf_spk_map_hd (fp_ts, fp_pdf, fp_fr, fp_info, fp_momentary_info, abs_info);
        }else{
            make_pdf_spk_map_xy (fp_ts, fp_pdf, fp_fr, fp_info, fp_momentary_info, abs_info);
        }
        
        
	fclose (fp_ts); 
	fclose (fp_pdf); 
	fclose (fp_fr); 
	fclose (fp_info);
	fclose (fp_momentary_info);
    
    }
    fclose (fp); 
    return 0;
}


void instruct()
{
	(void)fprintf(stderr,"\n Four directories must be specified:");
	(void)fprintf(stderr,"\n TS_INFO_DIR  - for the TS input file with pixels (or optionally TSD files with direction) and spikes in them for each time step in time-series (.ts) or time-series direction (.tsd) format");
	(void)fprintf(stderr,"\n PDF_DIR - for the output *.pdf file with pixels and distributions of spikes in them,");
	(void)fprintf(stderr,"\n RATE_INFO_DIR  - for the output *.fr file with pixels and firing rates in them.");
	(void)fprintf(stderr,"\n PINFO_DIR  - for the output *.pinfo file with pixels and positional (total, total-position, etc.) info  in them.");
	(void)fprintf(stderr,"\n            - for the output *.mpinfo file (option M).\n");

	(void)fprintf(stderr,"\n Call with (option and)  a file that contains list of files from TS\n");
	(void)fprintf(stderr,"\n Options:");
	(void)fprintf(stderr,"\n A   - Output the absolute value of the momentary positional info; by default the sign is preserved;");
	(void)fprintf(stderr,"\n B int - beginning of the time interval to be analyzed (in seconds); by default analysis is performed from the very beginning of the session;");	
    	(void)fprintf(stderr,"\n D   - Do directional information on a HD TSD file;");
	(void)fprintf(stderr,"\n d int  - head direction scale. Directions will be divided by this factor.");
	(void)fprintf(stderr,"\n E int - end of the time interval to be analyzed (in seconds); by default analysis is performed up to the very end of the session;");
	(void)fprintf(stderr,"\n H int  - maximal number of directions to consider; default is %d", MAX_HDSZ);
	(void)fprintf(stderr,"\n M   - DO NOT calculate momentary positional information. Store in PINFO_DIR/*.mpinfo.");
	(void)fprintf(stderr,"\n       .mpinfo Format: time_step m_pinfo instantaneous_rate");
	(void)fprintf(stderr,"\n                       time_step m_pinfo instantaneous_rate");
	(void)fprintf(stderr,"\n                       ...");
	(void)fprintf(stderr,"\n                       if m_pinfo is undefined m_pinfo = %0.4lf\n", UNDEFINED_INFO);
	(void)fprintf(stderr,"\n N int  - minimal sampling (pixels visited for not less than N time-steps are considered only); default is 10");
	(void)fprintf(stderr,"\n P   - put minimal data to *.pdf files; by default everything is written there");
	(void)fprintf(stderr,"\n S int  - maximal number of spikes during a time step; default is %d", MAX_SPSZ);
	(void)fprintf(stderr,"\n T float  - original time-step in the TS/TSD file in milliseconds; default is 116.6667 ms");
	(void)fprintf(stderr,"\n X int  - maximal number of pixels in x direction; default is %d", MAX_XSZ);
	(void)fprintf(stderr,"\n x int  - horizontal scale (for x), if it is different from the scale_x in TS file");
	(void)fprintf(stderr,"\n Y int  - maximal number of pixels in y direction; default is %d", MAX_YSZ);
	(void)fprintf(stderr,"\n y int  - vertical scale (for y), if it is different from the scale_y in TS file\n\n");
	exit(2);
}

void    check_env (char *name)
{
        if(getenv(name) == NULL){
                (void)fprintf(stderr,"%s is not defined\n", name);
		printf("\n\nHERE IS YOUR 'DIR' ENVIRONMENT\n\n"); system("env | grep DIR\n");
                instruct ();
        }
}
