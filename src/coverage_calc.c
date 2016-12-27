#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "subread.h" 
#include "core.h" 
#include "HelperFunctions.h"
#include "sambam-file.h" 
#include "gene-algorithms.h" 
#include "input-files.h" 
#include "hashtable.h"

#define COVERAGE_MAX_INT 0x7ffffff0 
unsigned long long all_counted;
typedef unsigned int coverage_bin_entry_t;
int is_BAM_input = 0;
int max_M = 10;
char input_file_name[300];
char output_file_name[300];
HashTable * cov_bin_table;


static struct option cov_calc_long_options[] =
{
	{"maxMOp",required_argument, 0, 'M'},
	{"primary",no_argument, 0, 0},
	{0, 0, 0, 0}
};


void calcCount_usage()
{
	SUBREADprintf("\ncoverageCount Version %s\n\n", SUBREAD_VERSION);
	SUBREADputs("  This program calculates the coverage of mapped reads at each location on");
	SUBREADputs("the reference genome. It generates a binary file for each chromosome by concate-");
	SUBREADputs("nating the coverage levels as 4-bytes integer numbers.");
	SUBREADputs("");
	SUBREADputs("Usage");
	SUBREADputs("");
	SUBREADputs("  ./coverageCount [options] -i <input_file> -o <output_prefix>");
	SUBREADputs("");
	SUBREADputs("Required arguments:");
	SUBREADputs("");
	SUBREADputs("  -i <string>  Name of input file in SAM or BAM format.");
	SUBREADputs("");
	SUBREADputs("  -o <string>  Prefix of the output files. Each output file contains Four-byte");
	SUBREADputs("               integer numbers");
	SUBREADputs("");
	SUBREADputs("Optional arguments:");
	SUBREADputs("");
	SUBREADputs("  --maxMOp <int> Maximum number of 'M' operations allowed in a CIGAR string.");
	SUBREADputs("               10 by default. Both 'X' and '=' are treated as 'M' and adjacent");
	SUBREADputs("               'M' operations are merged in the CIGAR string.");
	SUBREADputs("");
}

void add_chro(char *sam_h)
{
	char *chro_name = malloc(200);
	unsigned int chro_len = 0;

	char nch;
	int cur = 0, tabs=0, state = 0, txtcur = 0;

	chro_name[0]=0;
	while(1){
		nch = sam_h[cur++];
		if(!nch || nch == '\n')break;

		if(nch == '\t')
		{
			txtcur = 0;
			tabs ++;
			state = 0;
		}
		else if(nch == ':')
		{
			state ++;
		}
		else
		{
			if(state == 1 && tabs == 1)
			{
				chro_name[txtcur++]=nch;
				chro_name[txtcur]=0;
			}
			else if(state == 1 && tabs ==2)
				chro_len = chro_len*10 + (nch-'0');
		}
	}

	if(chro_name[0]==0 || chro_len<1)
	{
		SUBREADprintf("ERROR: incorrect SAM format: %s\n", sam_h);
	}

	void ** bin_entry = malloc(sizeof(void *)*2);
	coverage_bin_entry_t * chro_bin = calloc(sizeof(coverage_bin_entry_t) , chro_len);
	if(!chro_bin)
	{
		SUBREADprintf("ERROR: cannot allocate the memory block. You need at least 4GB of memory,\n");
	}
	bin_entry[0] = (void *)chro_bin;
	bin_entry[1] = (void *)(NULL + chro_len);
	HashTablePut(cov_bin_table, chro_name, bin_entry);
	SUBREADprintf("Added a new chromosome : %s [%u]\n", chro_name, chro_len);
}

void get_read_info(char * fl, char * chro, unsigned int * pos , char * cigar, int *flags){
	char * tmp_tok = NULL;

	char * tmp_res = strtok_r(fl, "\t", &tmp_tok);

	tmp_res = strtok_r(NULL, "\t", &tmp_tok);
	(*flags) = atoi(tmp_res);

	tmp_res = strtok_r(NULL, "\t", &tmp_tok); // chro
	strcpy(chro, tmp_res);

	tmp_res = strtok_r(NULL, "\t", &tmp_tok); // pos
	(*pos) = atoi(tmp_res);

	tmp_res = strtok_r(NULL, "\t", &tmp_tok); // qual
	tmp_res = strtok_r(NULL, "\t", &tmp_tok); // cigar
	strcpy(cigar, tmp_res);
}

int covCalc()
{

	cov_bin_table = HashTableCreate(200);
	HashTableSetHashFunction(cov_bin_table, fc_chro_hash);
	HashTableSetKeyComparisonFunction(cov_bin_table , fc_strcmp_chro);
	HashTableSetDeallocationFunctions(cov_bin_table , free, free);


	SamBam_FILE * in_fp = SamBam_fopen(input_file_name, is_BAM_input?SAMBAM_FILE_BAM:SAMBAM_FILE_SAM);
	char * line_buffer = malloc(3000);
	
	while(1)
	{
		char * is_ret = SamBam_fgets(in_fp, line_buffer, 2999, 1);
		if(!is_ret) break;
		if(line_buffer[0]=='@'){
			if(strstr(line_buffer,"@SQ\t"))
				add_chro(line_buffer);
		}
		else
		{
			char * Chros[FC_CIGAR_PARSER_ITEMS];
			unsigned int Staring_Points[FC_CIGAR_PARSER_ITEMS];
			unsigned short Staring_Read_Points[FC_CIGAR_PARSER_ITEMS];
			unsigned short Section_Lengths[FC_CIGAR_PARSER_ITEMS];
	
			int flags=0, x1, is_junc = 0;
			char cigar_str[200];
			char chro[200];
			unsigned int pos = 0;
			cigar_str[0]=0;
			chro[0]=0;

			get_read_info(line_buffer, chro, &pos, cigar_str, &flags);

			if(flags & 4) continue;

			void ** bin_entry = HashTableGet(cov_bin_table, chro);
			if(NULL == bin_entry)
			{
				SUBREADprintf("ERROR: The chromosome name is not in header:%s\n", chro);
			}

			coverage_bin_entry_t * chrbin = (coverage_bin_entry_t*) bin_entry[0];
			unsigned int chrlen = (void *)( bin_entry[1]) - NULL;
			int cigar_sections = RSubread_parse_CIGAR_string(chro, pos, cigar_str, max_M, Chros, Staring_Points, Staring_Read_Points, Section_Lengths, &is_junc);
			for(x1 = 0; x1 < cigar_sections; x1++)
			{
				unsigned int x2;
				for(x2 = Staring_Points[x1]; x2<Staring_Points[x1]+Section_Lengths[x1]; x2++)
				{
					if(x2 < chrlen) {
						if(chrbin[x2] <= COVERAGE_MAX_INT)chrbin[x2] ++;
						all_counted ++;
						if(all_counted % 10000000 == 0)
							SUBREADprintf("Processed %llu bases.\n", all_counted);
					} else {
						SUBREADprintf("%s:%s %u [%s] :: %u-%u\n", line_buffer , chro, pos, cigar_str, Staring_Points[x1], Staring_Points[x1]+Section_Lengths[x1]);
						SUBREADprintf("Read %s overhangs the boundary of chromosome %s (%u >= %u)\n", line_buffer, chro, x2, chrlen);
						//exit(-1);
					}
				}
			}
		}
	}

	free(line_buffer);

	SamBam_fclose(in_fp);


	SUBREADprintf("Processed totally %llu bases.\nNow write results.\n", all_counted);

	int bucket;
	KeyValuePair *cursor;
	for(bucket=0; bucket < cov_bin_table  -> numOfBuckets; bucket++)
	{
		cursor = cov_bin_table -> bucketArray[bucket];
		while (1)
		{
			if (!cursor) break;
			coverage_bin_entry_t * chrbin = (coverage_bin_entry_t*)(((void **) cursor -> value)[0]);
			unsigned int chrlen = (((void **) cursor -> value)[1]) - NULL; 
			char * chro = (char *)(cursor -> key);
			char out_name[340];
			sprintf(out_name,"%s-%s.bin", output_file_name, chro);

			FILE * fpo = fopen(out_name,"w");
			fwrite(chrbin, sizeof(coverage_bin_entry_t), chrlen, fpo);
			fclose(fpo);
			free(chrbin);

			SUBREADprintf("Wrote bin for %s\n", chro);
			cursor = cursor->next;
		}	
	}

	HashTableDestroy(cov_bin_table);

	SUBREADprintf("Calculation finished.\n");
	return 0;
}



#ifdef MAKE_STANDALONE
int main(int argc, char ** argv)
#else
int cov_calc_main(int argc, char ** argv)
#endif
{
	int ret = 0;
	int c=0;
	int option_index=0;
	input_file_name[0]=0;
	output_file_name[0]=0;
	max_M = 10;
	is_BAM_input=0;
	all_counted = 0;

	optind=0;
	opterr=1;
	optopt=63;

	while ((c = getopt_long (argc, argv, "BM:i:o:?", cov_calc_long_options, &option_index)) != -1)
		switch(c)
		{
			case 'M':
				if(!is_valid_digit_range(optarg, "maxMOp", 1 , 64))
					exit(-1);
				max_M = atoi(optarg);
			break;
			case 'i':
				strcpy(input_file_name, optarg);
			break;
			case 'o':
				strcpy(output_file_name, optarg);
			break;

			case '?':
			default :
				calcCount_usage();
				return -1;
	
		}
	

	if((!output_file_name[0])||(!input_file_name[0]))
	{
		calcCount_usage();
		return 0;
	}

	int is_bam = is_certainly_bam_file(input_file_name, NULL, NULL);

	if(1==is_bam) is_BAM_input = 1;
	else if(is_bam < 0)
	{
		ret = -1;
		SUBREADprintf("Unable to open input file '%s' or the input file is empty!\n", input_file_name);
	}

	ret = ret || covCalc();

	return ret;
}
