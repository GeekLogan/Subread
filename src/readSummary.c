/***************************************************************

   The Subread and Rsubread software packages are free
   software packages:
 
   you can redistribute it and/or modify it under the terms
   of the GNU General Public License as published by the 
   Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   Subread is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty
   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
   
   See the GNU General Public License for more details.

   Authors: Drs Yang Liao and Wei Shi

  ***************************************************************/
  
  
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>


#ifndef MAKE_STANDALONE
  #include <R.h>
#endif

#include <zlib.h>
#include <math.h>
#include <pthread.h>
#include <getopt.h>
#include "subread.h"
#include "interval_merge.h"
#include "core.h"
#include "gene-algorithms.h"
#include "sambam-file.h"
#include "input-files.h"
#include "hashtable.h"
#include "HelperFunctions.h"

/********************************************************************/
/********************************************************************/
/********************************************************************/
//  NEW FUNCTION FOR MULTI-THREADING
/********************************************************************/
/********************************************************************/
/********************************************************************/
#define FEATURE_NAME_LENGTH  256 
#define CHROMOSOME_NAME_LENGTH 256 
#define MAX_LINE_LENGTH 3000
#define FILE_TYPE_RSUBREAD 10
#define FILE_TYPE_GTF 100

#define ALLOW_ALL_MULTI_MAPPING 1
#define ALLOW_PRIMARY_MAPPING 2
#define MAX_FC_READ_LENGTH 10001

#define MAX_HIT_NUMBER 3000

typedef struct{
	char gene_name[FEATURE_NAME_LENGTH];
	unsigned int pos_first_base;
	unsigned int pos_last_base;
} fc_junction_gene_t;



typedef struct {
	int space;
	int used;
	fc_junction_gene_t ** genes;
} gene_info_list_t;

typedef struct {
	char chromosome_name_left[CHROMOSOME_NAME_LENGTH + 1];
	char chromosome_name_right[CHROMOSOME_NAME_LENGTH + 1];
	unsigned int last_exon_base_left;
	unsigned int first_exon_base_right;
} fc_junction_info_t;

typedef struct {
	unsigned int feature_name_pos;
	unsigned int start;
	unsigned int end;
	unsigned int sorted_order;

	unsigned short chro_name_pos_delta;
	char is_negative_strand;
} fc_feature_info_t;

typedef struct {
	unsigned long long assigned_reads;
	unsigned long long unassigned_ambiguous;
	unsigned long long unassigned_multimapping;
	unsigned long long unassigned_nofeatures;
	unsigned long long unassigned_unmapped;
	unsigned long long unassigned_mappingquality;
	unsigned long long unassigned_fragmentlength;
	unsigned long long unassigned_chimericreads;
	unsigned long long unassigned_secondary;
	unsigned long long unassigned_junction_condition;
	unsigned long long unassigned_duplicate;
} fc_read_counters;

typedef unsigned long long read_count_type_t;

typedef struct {
	unsigned short thread_id;
	unsigned long long int nreads_mapped_to_exon;
	unsigned long long int all_reads;
	//unsigned short current_read_length1;
	//unsigned short current_read_length2;
	read_count_type_t * count_table;
	read_count_type_t unpaired_fragment_no;
	unsigned int chunk_read_ptr;
	pthread_t thread_object;

	unsigned int hits_start_pos1[MAX_HIT_NUMBER];
	unsigned int hits_start_pos2[MAX_HIT_NUMBER];

	unsigned short hits_length1[MAX_HIT_NUMBER];
	unsigned short hits_length2[MAX_HIT_NUMBER];

	char * hits_chro1[MAX_HIT_NUMBER];
	char * hits_chro2[MAX_HIT_NUMBER];

	long hits_indices1 [MAX_HIT_NUMBER];
	long hits_indices2 [MAX_HIT_NUMBER];

	char ** scoring_buff_gap_chros;
	unsigned int * scoring_buff_gap_starts;
	unsigned short * scoring_buff_gap_lengths;

	unsigned int scoring_buff_numbers[MAX_HIT_NUMBER * 2];
	unsigned int scoring_buff_flags[MAX_HIT_NUMBER * 2];
	unsigned short scoring_buff_overlappings[MAX_HIT_NUMBER * 2];
	long scoring_buff_exon_ids[MAX_HIT_NUMBER * 2];

	char * chro_name_buff;
	z_stream * strm_buffer;

	HashTable * junction_counting_table;   // key: string chro_name \t last_base_previous_exont \t first_base_next_exon
	HashTable * splicing_point_table;
	fc_read_counters read_counters;

	SamBam_Alignment aln_buffer;
} fc_thread_thread_context_t;

#define REVERSE_TABLE_BUCKET_LENGTH 131072
#define REDUCE_TO_5_PRIME_END 5
#define REDUCE_TO_3_PRIME_END 3

typedef struct {
	unsigned int chro_number;
	unsigned int chro_features;
	unsigned int chro_feature_table_start;
	unsigned int chro_block_table_start;
	unsigned int chro_block_table_end;
	unsigned int chro_possible_length;

	unsigned short chro_reverse_table_current_size;
	unsigned int * reverse_table_start_index;
	//unsigned int * reverse_table_end_index;
} fc_chromosome_index_info;

typedef struct {
	int is_gene_level;
	int is_paired_end_input_file;
	int is_paired_end_mode_assign;
	int is_multi_overlap_allowed;
	int restricted_no_multi_overlap;
	int is_strand_checked;
	int is_both_end_required;
	int is_chimertc_disallowed;
	int is_PE_distance_checked;
	int is_multi_mapping_allowed;
	int is_SAM_file;
	int is_read_details_out;
	int is_junction_no_chro_shown;
	int is_SEPEmix_warning_shown;
	int is_unpaired_warning_shown;
	int is_stake_warning_shown;
	int is_split_or_exonic_only;
	int is_duplicate_ignored;
	int is_first_read_reversed;
	int is_second_read_straight;
	int do_not_sort;
	int reduce_5_3_ends_to_one;
	int isCVersion;
	int use_fraction_multi_mapping;
	int do_junction_counting;

	int need_calculate_overlap_len;
	int need_calculate_fragment_len;

	int min_mapping_quality_score;
	int min_paired_end_distance;
	int max_paired_end_distance;
	int max_M;
	int feature_block_size;
	int read_length;
	int line_length;
	int longest_chro_name;
	int five_end_extension;
	int three_end_extension;
	int fragment_minimum_overlapping;
	float fractional_minimum_overlapping; 
	int use_overlapping_break_tie;

	unsigned long long int all_reads;

	unsigned short thread_number;
	fc_thread_thread_context_t * thread_contexts;
	int is_all_finished;
	int sambam_chro_table_items;
	int is_input_bad_format;
	SamBam_Reference_Info * sambam_chro_table;
	pthread_spinlock_t sambam_chro_table_lock;

	SAM_pairer_context_t read_pairer;

	char * debug_command;
	char * unistr_buffer_space;
	long long max_BAM_header_size;
	unsigned int unistr_buffer_size;
	unsigned int unistr_buffer_used;
	HashTable * junction_features_table;
	HashTable * junction_bucket_table;
	fasta_contigs_t * fasta_contigs;
	HashTable * gene_name_table;	// gene_name -> gene_number
	HashTable * annot_chro_name_alias_table;	// name in annotation file -> alias name
	char alias_file_name[300];
	char input_file_name[300];
	char * input_file_short_name;
	char raw_input_file_name[300];
	char output_file_name[300];
	char output_file_path[300];
	char temp_file_dir[300];
	unsigned char ** gene_name_array;	// gene_internal_number -> gene_name 
	int input_file_unique;

	HashTable * exontable_chro_table;	// gene_name -> fc_chromosome_index_info structure (contains chro_number, feature_number, block_start, block_end, etc) 
	int exontable_nchrs;
	int exontable_exons;
	int * exontable_geneid;
	char * exontable_strand;
	char ** exontable_chr;
	long * exontable_start;
	long * exontable_stop;
	char feature_name_column[100];
	char gene_id_column[100];

	long * exontable_block_end_index;
	long * exontable_block_max_end;
	long * exontable_block_min_start;

	char ** exontable_anno_chrs;
	char * exontable_anno_chr_2ch;
	long * exontable_anno_chr_heads;

	FILE * SAM_output_fp;
	double start_time;

	char * cmd_rebuilt;

	char   redo;

	fc_read_counters read_counters;
	
} fc_thread_global_context_t;

unsigned int tick_time = 1000;


int fetch_boundaries(char * chroname,char * cigar, unsigned int pos, char strand, int *has_left, unsigned short *left_on_read, unsigned int *left_pos, int *has_right, unsigned short *right_on_read, unsigned int *right_pos, fc_junction_info_t *  result_junctions, int junction_space){

	int cigar_cursor = 0, nch, read_len = 0, ret = 0;
	unsigned int chro_cursor = pos, tmpi = 0;
	unsigned int right_boundary = 0;
	unsigned short left_clipped = 0;
	unsigned short right_clipped = 0;
	*has_right = 0;
	*has_left = 0;

	for(; (nch = cigar[cigar_cursor])!=0 ; cigar_cursor++){
		if(isdigit(nch)){
			tmpi = tmpi*10 + (nch - '0');
		} else {
			if (nch == 'S'){
				if(chro_cursor == pos) left_clipped = tmpi;else right_clipped=tmpi;
				read_len += tmpi;
			} else if(nch == 'M' || nch == 'D'){
				if(nch == 'M')read_len += tmpi;

				chro_cursor += tmpi;
				right_boundary = chro_cursor -1;
			} else if(nch == 'N'){
				unsigned int last_exon_last_base = chro_cursor - 1;
				unsigned int next_exon_first_base = chro_cursor + tmpi;
				chro_cursor += tmpi;

				if(ret < junction_space){
					result_junctions[ret].last_exon_base_left = last_exon_last_base;
					result_junctions[ret].first_exon_base_right = next_exon_first_base;
					strcpy(result_junctions[ret].chromosome_name_left, chroname);
					strcpy(result_junctions[ret].chromosome_name_right, chroname);

					ret ++;
				}


			} else if(nch == 'I') read_len += tmpi;
			tmpi = 0;
		}
	}
	if(left_clipped){
		*has_left = 1;
		*left_on_read = left_clipped;
		*left_pos = pos;
	}
	if(right_clipped){
		*has_right = 1;
		*right_on_read = read_len - right_clipped - 1;
		*right_pos = right_boundary;
	}
	return ret;
}

// This function parses the cigar string and returns the number of exon-exon junctions found in the cigar.
// It returns 0 if no junctions are found.
int calc_junctions_from_cigar(fc_thread_global_context_t * global_context, int flag, char * chroname, unsigned int pos, char * cigar , char * extra_tags, fc_junction_info_t * result_junctions){
	unsigned short boundaries_inclusive_base_on_read[global_context -> max_M];
	unsigned int boundaries_inclusive_base_pos[global_context -> max_M];
	char boundaries_chromosomes[global_context -> max_M][MAX_CHROMOSOME_NAME_LEN];
	char boundaries_extend_to_left_on_read[global_context -> max_M];
	int boundaries = 0;

	int cigar_cursor = 0, nch, ret = 0, read_len = 0, x1, x2;
	unsigned int chro_cursor = pos, tmpi = 0;
	unsigned int right_boundary = 0;
	unsigned short left_clipped = 0;
	unsigned short right_clipped = 0;

	for(; (nch = cigar[cigar_cursor])!=0 ; cigar_cursor++){
		if(isdigit(nch)){
			tmpi = tmpi*10 + (nch - '0');
		} else {
			if (nch == 'S'){
				if(chro_cursor == pos) left_clipped = tmpi;else right_clipped=tmpi;
				read_len += tmpi;
			} else if(nch == 'M' || nch == 'D'){
				if(nch == 'M')read_len += tmpi;

				chro_cursor += tmpi;
				right_boundary = chro_cursor -1;
			} else if(nch == 'N'){
				unsigned int last_exon_last_base = chro_cursor - 1;
				unsigned int next_exon_first_base = chro_cursor + tmpi;
				if(ret <= global_context -> max_M - 1){
					result_junctions[ret].last_exon_base_left = last_exon_last_base;
					result_junctions[ret].first_exon_base_right = next_exon_first_base;
					strcpy(result_junctions[ret].chromosome_name_left, chroname);
					strcpy(result_junctions[ret].chromosome_name_right, chroname);

					ret ++;
				}
				chro_cursor += tmpi;
			} else if(nch == 'I') read_len += tmpi;
			tmpi = 0;
		}
	}
	if(left_clipped){
		strcpy(boundaries_chromosomes[boundaries] , chroname);
		boundaries_extend_to_left_on_read[boundaries] = 0;
		boundaries_inclusive_base_pos[boundaries] = pos;
		boundaries_inclusive_base_on_read[boundaries++] = left_clipped;
	}
	if(right_clipped){
		strcpy(boundaries_chromosomes[boundaries] , chroname);
		boundaries_extend_to_left_on_read[boundaries] = 1;
		boundaries_inclusive_base_pos[boundaries] = chro_cursor - 1;
		boundaries_inclusive_base_on_read[boundaries++] = read_len - right_clipped - 1;
	}

	int tag_cursor=0;

	//if(strstr(extra_tags, "CG:Z")) {
	//	SUBREADprintf("CIGAR=%s, EXTRA=%s\n", cigar, extra_tags);
	//}
	int status = PARSE_STATUS_TAGNAME;
	char tag_name[2], typechar=0;
	int tag_inner_cursor=0;

	char read_main_strand = (((flag & 0x10) == 0x10) == ((flag & 0x40)==0x40))?'-':'+';
	char current_fusion_char[MAX_CHROMOSOME_NAME_LEN];
	unsigned int current_fusion_pos = 0;
	char current_fusion_strand = 0;
	char current_fusion_cigar[global_context -> max_M * 15];
	current_fusion_cigar [0] =0;
	current_fusion_char [0]=0;

	while(1){
		int nch = extra_tags[tag_cursor];
		if(status == PARSE_STATUS_TAGNAME){
			tag_name[tag_inner_cursor++] = nch;
			if(tag_inner_cursor == 2){
				status = PARSE_STATUS_TAGTYPE;
				tag_cursor += 1;
				assert(extra_tags[tag_cursor] == ':');
			}
		}else if(status == PARSE_STATUS_TAGTYPE){
			typechar = nch;
			tag_cursor +=1;
			assert(extra_tags[tag_cursor] == ':');
			tag_inner_cursor = 0;
			status = PARSE_STATUS_TAGVALUE;
		}else if(status == PARSE_STATUS_TAGVALUE){
			if(nch == '\t' || nch == 0){
				if(current_fusion_cigar[0] && current_fusion_char[0] && current_fusion_pos && current_fusion_strand){

					unsigned int left_pos = 0, right_pos = 0;
					unsigned short left_on_read = 0, right_on_read = 0;
					int has_left = 0, has_right = 0;

					unsigned int start_pos = current_fusion_pos;
					if(current_fusion_strand!=read_main_strand)
						start_pos = find_left_end_cigar(current_fusion_pos, current_fusion_cigar);

					ret += fetch_boundaries(current_fusion_char, current_fusion_cigar, start_pos, current_fusion_strand, &has_left, &left_on_read, &left_pos, &has_right, &right_on_read, &right_pos, result_junctions + ret, global_context -> max_M - ret );

					if(has_left){
						strcpy(boundaries_chromosomes[boundaries] , current_fusion_char);
						boundaries_extend_to_left_on_read[boundaries] = 0;
						boundaries_inclusive_base_pos[boundaries] = left_pos;
						boundaries_inclusive_base_on_read[boundaries++] = left_on_read;
					}
					if(has_right){
						strcpy(boundaries_chromosomes[boundaries] , current_fusion_char);
						boundaries_extend_to_left_on_read[boundaries] = 1;
						boundaries_inclusive_base_pos[boundaries] = right_pos;
						boundaries_inclusive_base_on_read[boundaries++] = right_on_read;
					}
	

			//		SUBREADprintf("BOUND_EXT: %s:%u (at %u) (%c)  ~  %s:%u (at %u) (%c)\n", current_fusion_char, left_pos, left_on_read, has_left?'Y':'X' , current_fusion_char, right_pos, right_on_read,  has_right?'Y':'X');

					current_fusion_pos = 0;
					current_fusion_strand = 0;
					current_fusion_cigar [0] =0;
					current_fusion_char [0]=0;
				}

				tag_inner_cursor = 0;
				status = PARSE_STATUS_TAGNAME;
			}else{
				if(tag_name[0]=='C' && tag_name[1]=='C' && typechar == 'Z'){
					current_fusion_char[tag_inner_cursor++]=nch;
					current_fusion_char[tag_inner_cursor]=0;
				}else if(tag_name[0]=='C' && tag_name[1]=='G' && typechar == 'Z'){
					current_fusion_cigar[tag_inner_cursor++]=nch;
					current_fusion_cigar[tag_inner_cursor]=0;
				}else if(tag_name[0]=='C' && tag_name[1]=='P' && typechar == 'i'){
					current_fusion_pos = current_fusion_pos * 10 + (nch - '0');
				}else if(tag_name[0]=='C' && tag_name[1]=='T' && typechar == 'Z'){
					current_fusion_strand = nch;
				}
			}
		}

		if(nch == 0){
			assert(status == PARSE_STATUS_TAGNAME);
			break;
		}

		tag_cursor++;
	}


	//for(x1 = 0; x1 < boundaries; x1++)
	//	SUBREADprintf("HAS: LR:%d, READ:%d\n", boundaries_extend_to_left_on_read[x1], boundaries_inclusive_base_on_read[x1]);

	for(x1 = 0; x1 < boundaries; x1++)
		for(x2 = 0; x2 < boundaries; x2++){
			if(x1==x2) continue;
			if(boundaries_chromosomes[x1][0]==0 || boundaries_chromosomes[x2][0]==0) continue;
			if(boundaries_extend_to_left_on_read[x1] == 1 && boundaries_extend_to_left_on_read[x2] == 0){
				if( boundaries_inclusive_base_on_read[x1] == boundaries_inclusive_base_on_read[x2]-1 ){

					if(ret <= global_context -> max_M - 1){
						result_junctions[ret].last_exon_base_left = boundaries_inclusive_base_pos[x1];
						result_junctions[ret].first_exon_base_right = boundaries_inclusive_base_pos[x2];
						strcpy(result_junctions[ret].chromosome_name_left, boundaries_chromosomes[x1]);
						strcpy(result_junctions[ret].chromosome_name_right, boundaries_chromosomes[x2]);
						ret++;
					}


	//				SUBREADprintf("MATCH: %d ~ %d\n", boundaries_inclusive_base_on_read[x1], boundaries_inclusive_base_on_read[x2]);
					boundaries_chromosomes[x1][0]=0;
					boundaries_chromosomes[x2][0]=0;
				}
			}
		}

	//for(x1 = 0; x1 < boundaries; x1++)
	//	if(boundaries_chromosomes[x1][0])
	//		SUBREADprintf("LEFT: LR:%d, READ:%d\n", boundaries_extend_to_left_on_read[x1], boundaries_inclusive_base_on_read[x1]);
	return ret;
}


unsigned int unistr_cpy(fc_thread_global_context_t * global_context, char * str, int strl)
{
	unsigned int ret;
	if(global_context->unistr_buffer_used + strl >= global_context->unistr_buffer_size-1)
	{
		if( global_context->unistr_buffer_size < 3435973835u) // 4G / 5 * 4 - 5
		{
			global_context -> unistr_buffer_size = global_context->unistr_buffer_size /4 *5;
			global_context -> unistr_buffer_space = realloc(global_context -> unistr_buffer_space, global_context->unistr_buffer_size);
		}
		else
		{
			SUBREADprintf("Error: exceed memory limit (4GB) for storing annotation data.\n");
			return 0xffffffffu;
		}
	}

	strcpy(global_context -> unistr_buffer_space + global_context->unistr_buffer_used, str);
	ret = global_context->unistr_buffer_used;

	global_context->unistr_buffer_used += strl +1;

	return ret;
}

int print_FC_configuration(fc_thread_global_context_t * global_context, char * annot, char * sam, char * out, int is_sam, int is_GTF, int *n_input_files, int isReadSummaryReport)
{
	char * tmp_ptr1 = NULL , * next_fn, *sam_used = malloc(strlen(sam)+300), sam_ntxt[30],bam_ntxt[30], next_ntxt[50];
	int nfiles=1, nBAMfiles = 0, nNonExistFiles = 0;

	sprintf(sam_used, "%s/featureCounts_test_file_writable.tmp", global_context -> temp_file_dir);
	FILE * fp = fopen(sam_used,"w");
	if(fp){
		fclose(fp);
		unlink(sam_used);
	}else{
		SUBREADprintf("\nERROR: temporary directory is not writable: '%s'\n\n", global_context -> temp_file_dir);
		return 1;
	}

	strcpy(sam_used, sam);
	nfiles = 0;
	while(1)
	{
		next_fn = strtok_r(nfiles==0?sam_used:NULL, ";", &tmp_ptr1);
		if(next_fn == NULL || strlen(next_fn)<1) break;
		nfiles++;

		long long BAM_header_size = -1;
		int file_probe = is_certainly_bam_file(next_fn, NULL, &BAM_header_size);
		//SUBREADprintf(" >>> %s : header=%lld , BSIZE=%lld\n", next_fn,BAM_header_size, global_context -> max_BAM_header_size );
		if(BAM_header_size>0) global_context -> max_BAM_header_size = max( global_context -> max_BAM_header_size , BAM_header_size + 180000);
		if(file_probe==-1){
			nNonExistFiles++;
			SUBREADprintf("\nERROR: invalid parameter: '%s'\n\n", next_fn);
			return 1;
		}
		if(file_probe == 1) nBAMfiles++;		
	}

	SUBREADputs("");
	print_subread_logo();
	SUBREADputs("");
	print_in_box(80,1,1,"featureCounts setting");
	print_in_box(80,0,0,"");

	sam_ntxt[0]=0;
	bam_ntxt[0]=0;
	next_ntxt[0]=0;

	if(nNonExistFiles)
		sprintf(next_ntxt, "%d unknown file%s", nNonExistFiles, nNonExistFiles>1?"s":"");
	if(nBAMfiles)
		sprintf(bam_ntxt, "%d BAM file%s  ", nBAMfiles, nBAMfiles>1?"s":"");
	if(nfiles-nNonExistFiles-nBAMfiles)
		sprintf(sam_ntxt, "%d SAM file%s  ", nfiles-nNonExistFiles-nBAMfiles , (nfiles-nNonExistFiles-nBAMfiles)>1?"s":"");


	strcpy(sam_used, sam);

	print_in_box(80,0,0,"            Input files : %s%s%s", sam_ntxt, bam_ntxt, next_ntxt);
	nfiles=0;

	while(1)
	{
		next_fn = strtok_r(nfiles==0?sam_used:NULL, ";", &tmp_ptr1);
		if(next_fn == NULL || strlen(next_fn)<1) break;
		int is_first_read_PE = 0 , file_probe = is_certainly_bam_file(next_fn, &is_first_read_PE, NULL);

		char file_chr = 'S';
		if(file_probe == -1) file_chr = '?';
		else if(is_first_read_PE == 1) file_chr = 'P';
		//file_chr = 'o';

		print_in_box(94,0,0,"                          %c[32m%c%c[36m %s%c[0m",CHAR_ESC, file_chr,CHAR_ESC, next_fn,CHAR_ESC);
		nfiles++;
	}

	(*n_input_files) = nfiles;
	print_in_box(80,0,0,"");
	print_in_box(80,0,0,"            Output file : %s", out);
	print_in_box(80,0,0,"                Summary : %s.summary", out);
	print_in_box(80,0,0,"             Annotation : %s (%s)", annot, is_GTF?"GTF":"SAF");
	if(isReadSummaryReport){
		print_in_box(80,0,0,"     Assignment details : <input_file>.featureCounts");
		print_in_box(80,0,0,"                     (Note that files are saved to the output directory)");
		print_in_box(80,0,0,"");
	}
	if(global_context -> do_junction_counting)
		print_in_box(80,0,0,"      Junction Counting : <output_file>.jcounts");

	if(global_context -> alias_file_name[0])
		print_in_box(80,0,0,"  Chromosome alias file : %s", global_context -> alias_file_name);
	print_in_box(80,0,0,"     Dir for temp files : %s", global_context->temp_file_dir);

	print_in_box(80,0,0,"");
	print_in_box(80,0,0,"                Threads : %d", global_context->thread_number);
	print_in_box(80,0,0,"                  Level : %s level", global_context->is_gene_level?"meta-feature":"feature");
	print_in_box(80,0,0,"             Paired-end : %s", global_context->is_paired_end_mode_assign?"yes":"no");
	if(global_context -> do_not_sort && global_context->is_paired_end_mode_assign)
	{
		print_in_box(80,0,0,"       Sorting PE Reads : never");
		print_in_box(80,0,0,"");
	}

	print_in_box(80,0,0,"        Strand specific : %s", global_context->is_strand_checked?(global_context->is_strand_checked==1?"stranded":"reversely stranded"):"no");
	char * multi_mapping_allow_mode = "not counted";
	if(global_context->is_multi_mapping_allowed == ALLOW_PRIMARY_MAPPING)
		multi_mapping_allow_mode = "primary only";
	else if(global_context->is_multi_mapping_allowed == ALLOW_ALL_MULTI_MAPPING)
		multi_mapping_allow_mode = global_context -> use_fraction_multi_mapping?"counted (as fractions)": "counted (as integer)";

	print_in_box(80,0,0,"     Multimapping reads : %s", multi_mapping_allow_mode);
	print_in_box(80,0,0,"Multi-overlapping reads : %s", global_context->is_multi_overlap_allowed?"counted":"not counted");
	if(global_context -> is_split_or_exonic_only)
		print_in_box(80,0,0,"       Split alignments : %s", (1 == global_context -> is_split_or_exonic_only)?"only split alignments":"only exonic alignments");
	if(global_context -> fragment_minimum_overlapping !=1)
		print_in_box(80,0,0,"      Overlapping bases : %d", global_context -> fragment_minimum_overlapping);
	if(global_context -> fractional_minimum_overlapping !=1)
		print_in_box(81,0,0,"      Overlapping bases : %0.1f%%%%", global_context -> fractional_minimum_overlapping*100);
	if(global_context -> five_end_extension || global_context -> three_end_extension)
		print_in_box(80,0,0,"        Read extensions : %d on 5' and %d on 3' ends", global_context -> five_end_extension , global_context -> three_end_extension);
	if(global_context -> reduce_5_3_ends_to_one)
		print_in_box(80,0,0,"      Read reduction to : %d' end" , global_context -> reduce_5_3_ends_to_one == REDUCE_TO_5_PRIME_END ?5:3);
	if(global_context -> is_duplicate_ignored)
		print_in_box(80,0,0,"       Duplicated Reads : ignored");
	//print_in_box(80,0,0,"      Read orientations : %c%c", global_context->is_first_read_reversed?'r':'f', global_context->is_second_read_straight?'f':'r' );

	if(global_context->is_paired_end_mode_assign)
	{
		print_in_box(80,0,0,"");
		print_in_box(80,0,0,"         Chimeric reads : %s", global_context->is_chimertc_disallowed?"not counted":"counted");
		print_in_box(80,0,0,"       Both ends mapped : %s", global_context->is_both_end_required?"required":"not required");

		if(global_context->is_PE_distance_checked)
			print_in_box(80,0,0,"        Fragment length : %d - %d", global_context -> min_paired_end_distance, global_context -> max_paired_end_distance);
	}

	print_in_box(80,0,0,"");
	print_in_box(80,2,1,"http://subread.sourceforge.net/");
	SUBREADputs("");
	print_in_box(80,1,1,"Running");
	print_in_box(80,0,0,"");
	if( global_context -> max_BAM_header_size > 32 * 1024 * 1024 ){
	}
	if(global_context->annot_chro_name_alias_table)
		print_in_box(80,0,0,"%ld chromosome name aliases are loaded.", global_context -> annot_chro_name_alias_table ->numOfElements);

	free(sam_used);
	return 0;
}

void print_FC_results(fc_thread_global_context_t * global_context)
{
	print_in_box(89,0,1,"%c[36mRead assignment finished.%c[0m", CHAR_ESC, CHAR_ESC);
	print_in_box(80,0,0,"");
	print_in_box(80,2,1,"http://subread.sourceforge.net/");
	SUBREADputs("");
	return;


	if(0){
		print_in_box(80,1,1,"Summary");
		print_in_box(80,0,0,"");
		if(global_context->is_paired_end_mode_assign)
			print_in_box(80,0,0,"        All fragments : %llu", global_context -> all_reads);
		else
			print_in_box(80,0,0,"            All reads : %llu", global_context -> all_reads);

		if(global_context->is_gene_level)
			print_in_box(80,0,0,"        Meta-features : %lu", global_context -> gene_name_table -> numOfElements);
		else
			print_in_box(80,0,0,"             Features : %u", global_context -> exontable_exons);

		if(global_context->is_paired_end_mode_assign)
			print_in_box(80,0,0,"   Assigned fragments : %llu", global_context -> read_counters.assigned_reads);
		else
			print_in_box(80,0,0,"       Assigned reads : %llu", global_context -> read_counters.assigned_reads);

		print_in_box(80,0,0,"            Time cost : %.3f minutes", (miltime() - global_context -> start_time)/60);
		print_in_box(80,0,0,"");
		print_in_box(80,2,1,"http://subread.sourceforge.net/");
	}
	SUBREADputs("");
}

int fc_strcmp(const void * s1, const void * s2)
{
	return strcmp((char*)s1, (char*)s2);
}


int is_comment_line(const char * l, int file_type, unsigned int lineno)
{
	int tabs = 0, xk1 = 0;
	if(l[0]=='#') return 1;

	if(isalpha(l[0]) && file_type == FILE_TYPE_RSUBREAD)
	{
		char target_chr[16];
		memcpy(target_chr, l, 16);
		for(xk1=0; xk1<16; xk1++)
			target_chr[xk1] = tolower(target_chr[xk1]);

		if(memcmp(target_chr, "geneid\tchr\tstart",16)==0) return 1;
	}

	xk1=0;
	while(l[xk1]) tabs += (l[xk1++] == '\t');

	return tabs < ((file_type == FILE_TYPE_GTF)?8:4);
}

void register_junc_feature(fc_thread_global_context_t *global_context, char * feature_name, char * chro, unsigned int start, unsigned int stop){
	HashTable * gene_table = HashTableGet(global_context -> junction_features_table, chro);
	//SUBREADprintf("REG %s : %p\n", chro, gene_table);
	if(NULL == gene_table){
		gene_table = HashTableCreate(48367);
		HashTableSetDeallocationFunctions(gene_table, NULL, free);
		HashTableSetKeyComparisonFunction(gene_table, fc_strcmp);
		HashTableSetHashFunction(gene_table, fc_chro_hash);

		char * new_name = malloc(strlen(chro)+1);
		strcpy(new_name, chro);
		HashTablePut(global_context -> junction_features_table, new_name, gene_table);
	}
	fc_junction_gene_t * gene_info = HashTableGet(gene_table, feature_name);
	if(NULL == gene_info){
		gene_info = malloc(sizeof(fc_junction_gene_t));
		strcpy(gene_info -> gene_name, feature_name);
		gene_info -> pos_first_base = start;
		gene_info -> pos_last_base = stop;

		HashTablePut(gene_table, gene_info -> gene_name, gene_info);
	}else{
		gene_info -> pos_first_base = min(start, gene_info -> pos_first_base);
		gene_info -> pos_last_base = max(stop, gene_info -> pos_last_base);
	}
}

void free_bucket_table_list(void * pv){
	gene_info_list_t * list = (gene_info_list_t*) pv;
	free(list -> genes);
	free(list);
}

#define JUNCTION_BUCKET_STEP (128*1024)

int locate_junc_features(fc_thread_global_context_t *global_context, char * chro, unsigned int pos, fc_junction_gene_t ** ret_info, int max_ret_info_size){
	gene_info_list_t * list = NULL;
	char bucket_key[CHROMOSOME_NAME_LENGTH + 20];

	if(global_context -> annot_chro_name_alias_table) {
		char * anno_chro_name = HashTableGet( global_context -> annot_chro_name_alias_table , chro);
		if(anno_chro_name){
			sprintf(bucket_key, "%s:%u", anno_chro_name, pos - pos % JUNCTION_BUCKET_STEP);
			list = HashTableGet(global_context -> junction_bucket_table, bucket_key);
		}
	}

	if(list == NULL){
		sprintf(bucket_key, "%s:%u", chro, pos - pos % JUNCTION_BUCKET_STEP);
		list = HashTableGet(global_context -> junction_bucket_table, bucket_key);
	}

	if(list == NULL && strlen(chro)>3 && memcmp(chro, "chr", 3)==0){
		sprintf(bucket_key, "%s:%u", chro+3, pos - pos % JUNCTION_BUCKET_STEP);
		list = HashTableGet(global_context -> junction_bucket_table, bucket_key);
	}

	if(list == NULL){
		sprintf(bucket_key, "chr%s:%u", chro, pos - pos % JUNCTION_BUCKET_STEP);
		list = HashTableGet(global_context -> junction_bucket_table, bucket_key);
	}

	int ret = 0;

	if(list){
		int x1; 
		for(x1 = 0; x1 < list -> used; x1++){
			fc_junction_gene_t * gene_info = list -> genes[x1];
			if(gene_info -> pos_first_base <= pos && gene_info -> pos_last_base >= pos){
				if(ret < max_ret_info_size)
					ret_info [ret ++] = gene_info;
			}
		}
	}
	
	return ret;
}

// This function loads annotations from the file.
// It returns the number of featres loaded, or -1 if something is wrong. 
// Memory will be allowcated in this function. The pointer is saved in *loaded_features.
// The invoker must release the memory itself.

int load_feature_info(fc_thread_global_context_t *global_context, const char * annotation_file, int file_type, fc_feature_info_t ** loaded_features)
{
	unsigned int features = 0, xk1 = 0, lineno=0;
	char * file_line = malloc(MAX_LINE_LENGTH+1);
	FILE * fp = f_subr_open(annotation_file,"r"); 
	int is_GFF_warned = 0;
	if(!fp) return -1;

	HashTable * chro_name_table = HashTableCreate(1603);
	HashTableSetHashFunction(chro_name_table, fc_chro_hash);
	HashTableSetKeyComparisonFunction(chro_name_table, fc_strcmp_chro);
	global_context -> longest_chro_name = 0;

	if(global_context -> do_junction_counting){
		global_context -> junction_bucket_table = HashTableCreate(76037);
		HashTableSetDeallocationFunctions(global_context -> junction_bucket_table, free, free_bucket_table_list);
		HashTableSetKeyComparisonFunction(global_context -> junction_bucket_table, fc_strcmp);
		HashTableSetHashFunction(global_context -> junction_bucket_table, fc_chro_hash);

		global_context -> junction_features_table = HashTableCreate(1603);
		HashTableSetDeallocationFunctions(global_context -> junction_features_table, free, (void (*)(void *))HashTableDestroy);
		HashTableSetKeyComparisonFunction(global_context -> junction_features_table, fc_strcmp);
		HashTableSetHashFunction(global_context -> junction_features_table, fc_chro_hash);
	}


	// first scan: get the chromosome size, etc
	while(1)
	{
		char * fgets_ret = fgets(file_line, MAX_LINE_LENGTH, fp);
		char * token_temp, *chro_name;
		fc_chromosome_index_info * chro_stab;
		unsigned int feature_pos = 0;
		if(!fgets_ret) break;

		lineno++;
		if(is_comment_line(file_line, file_type, lineno-1))continue;
		if(file_type == FILE_TYPE_GTF)
		{
			chro_name = strtok_r(file_line,"\t",&token_temp);
			strtok_r(NULL,"\t", &token_temp); // lib_name (not needed)
			char * feature_type = strtok_r(NULL,"\t", &token_temp);
			if(strcmp(feature_type, global_context -> feature_name_column)==0)
			{
				strtok_r(NULL,"\t", &token_temp); // feature_start
				feature_pos = atoi(strtok_r(NULL,"\t", &token_temp));// feature_end
				features++;
			}
			else chro_name = NULL;
		}
		else
		{
			strtok_r(file_line,"\t", &token_temp);
			chro_name = strtok_r(NULL,"\t",&token_temp);
			strtok_r(NULL,"\t",&token_temp);	// feature_start
			feature_pos = atoi(strtok_r(NULL,"\t", &token_temp));// feature_end

			features++;
		}

		if(chro_name)
		{
			if(strlen(chro_name)>=CHROMOSOME_NAME_LENGTH) 
				chro_name[CHROMOSOME_NAME_LENGTH-1]=0;
			chro_stab = HashTableGet(chro_name_table, chro_name);

			if(chro_stab)
			{
				chro_stab -> chro_possible_length = max(chro_stab -> chro_possible_length , feature_pos+1);
			}else
			{
				char * tmp_chro_name = malloc(CHROMOSOME_NAME_LENGTH);
				term_strncpy(tmp_chro_name, chro_name, CHROMOSOME_NAME_LENGTH);
				chro_stab = calloc(sizeof(fc_chromosome_index_info),1);
				chro_stab -> chro_number = chro_name_table->numOfElements;
				chro_stab -> chro_possible_length = feature_pos+1;
				chro_stab -> reverse_table_start_index = NULL;
				HashTablePut(chro_name_table, tmp_chro_name, chro_stab);
			}

			chro_stab -> chro_features ++;
		}
	}

	fseek(fp,0,SEEK_SET);

	fc_feature_info_t * ret_features = malloc(sizeof(fc_feature_info_t) * features);

	lineno = 0;
	while(xk1 < features)
	{
		int is_gene_id_found = 0;
		fgets(file_line, MAX_LINE_LENGTH, fp);
		lineno++;
		char * token_temp;
		if(is_comment_line(file_line, file_type, lineno-1))continue;

		if(file_type == FILE_TYPE_RSUBREAD)
		{
			char * feature_name = strtok_r(file_line,"\t",&token_temp);
			int feature_name_len = strlen(feature_name);
			if(feature_name_len > FEATURE_NAME_LENGTH) feature_name[FEATURE_NAME_LENGTH -1 ] = 0;
			ret_features[xk1].feature_name_pos = unistr_cpy(global_context, (char *)feature_name, feature_name_len);

			char * seq_name = strtok_r(NULL,"\t", &token_temp);
			int chro_name_len = strlen(seq_name);
			if(chro_name_len > CHROMOSOME_NAME_LENGTH) seq_name[CHROMOSOME_NAME_LENGTH -1 ] = 0;
			unsigned int chro_name_pos = unistr_cpy(global_context, (char *)seq_name, chro_name_len);
			global_context -> longest_chro_name = max(chro_name_len, global_context -> longest_chro_name);


			char * start_ptr = strtok_r(NULL,"\t", &token_temp);
			char * end_ptr = strtok_r(NULL,"\t", &token_temp);

			if(start_ptr == NULL || end_ptr == NULL){
				SUBREADprintf("\nWarning: the format on the %d-th line is wrong.\n", lineno);
			}
			long long int tv1 = atoll(start_ptr);
			long long int tv2 = atoll(end_ptr);

			if( isdigit(start_ptr[0]) && isdigit(end_ptr[0]) ){
				if(strlen(start_ptr) > 10 || strlen(end_ptr) > 10 || tv1 > 0x7fffffff || tv2> 0x7fffffff){
					SUBREADprintf("\nError: Line %d contains a coordinate greater than 2^31!\n", lineno);
					return -2;
				}
			}else{
				SUBREADprintf("\nError: Line %d contains a format error. The expected annotation format is SAF.\n", lineno);
				return -2;
			}

			ret_features[xk1].chro_name_pos_delta = chro_name_pos - ret_features[xk1].feature_name_pos;
			ret_features[xk1].start = atoi( start_ptr );// start 
			if(ret_features[xk1].start>0x7fffffff)
			{
				ret_features[xk1].start = 0;
				print_in_box(80,0,0,"WARNING the %d-th line has a negative chro coordinate.", lineno);
			}

			ret_features[xk1].end = atoi( end_ptr );//end 
			if(ret_features[xk1].end>0x7fffffff)
			{
				ret_features[xk1].end = 0;
				print_in_box(80,0,0,"WARNING the %d-th line has a negative chro coordinate.", lineno);
			}




			char * strand_str = strtok_r(NULL,"\t", &token_temp); 
			if(strand_str == NULL)
				ret_features[xk1].is_negative_strand = 0;
			else
				ret_features[xk1].is_negative_strand = ('-' ==strand_str[0]);
			ret_features[xk1].sorted_order = xk1;

			int bin_location = ret_features[xk1].start / REVERSE_TABLE_BUCKET_LENGTH;
			
			fc_chromosome_index_info * chro_stab = HashTableGet(chro_name_table, seq_name);
			if(!chro_stab -> reverse_table_start_index)
			{
				chro_stab -> reverse_table_start_index = malloc(sizeof(int) *( chro_stab->chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH +2));
				memset(chro_stab -> reverse_table_start_index, 0 , sizeof(int) *( chro_stab->chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH +2));
			}
			chro_stab -> reverse_table_start_index[bin_location]++;

			is_gene_id_found = 1;

			assert(feature_name);
			if(global_context -> do_junction_counting)
				register_junc_feature(global_context , feature_name, seq_name, ret_features[xk1].start, ret_features[xk1].end);

			xk1++;
		}
		else if(file_type == FILE_TYPE_GTF)
		{
			char feature_name_tmp[FEATURE_NAME_LENGTH];
			sprintf(feature_name_tmp, "LINE_%07u", xk1 + 1);
			char * seq_name = strtok_r(file_line,"\t",&token_temp);
			strtok_r(NULL,"\t", &token_temp);// source
			char * feature_type = strtok_r(NULL,"\t", &token_temp);// feature_type
			if(strcmp(feature_type, global_context -> feature_name_column)==0)
			{
				char * start_ptr = strtok_r(NULL,"\t", &token_temp);
				char * end_ptr = strtok_r(NULL,"\t", &token_temp);

				if(start_ptr == NULL || end_ptr == NULL){
					SUBREADprintf("\nWarning: the format on the %d-th line is wrong.\n", lineno);
				}
				long long int tv1 = atoll(start_ptr);
				long long int tv2 = atoll(end_ptr);

				if( isdigit(start_ptr[0]) && isdigit(end_ptr[0]) ){
					if(strlen(start_ptr) > 10 || strlen(end_ptr) > 10 || tv1 > 0x7fffffff || tv2> 0x7fffffff){
						SUBREADprintf("\nError: Line %d contains a coordinate greater than 2^31!\n", lineno);
						return -2;
					}
				}else{
					SUBREADprintf("\nError: Line %d contains a format error. The expected annotation format is GTF/GFF.\n", lineno);
					return -2;
				}
				ret_features[xk1].start = atoi(start_ptr);// start 
				ret_features[xk1].end = atoi(end_ptr);//end 

				if(ret_features[xk1].start < 1 || ret_features[xk1].end<1 ||  ret_features[xk1].start > 0x7fffffff ||  ret_features[xk1].end > 0x7fffffff || ret_features[xk1].start > ret_features[xk1].end)
					SUBREADprintf("\nWarning: the feature on the %d-th line has zero coordinate or zero lengths\n\n", lineno);


				strtok_r(NULL,"\t", &token_temp);// score 
				ret_features[xk1].is_negative_strand = ('-' == (strtok_r(NULL,"\t", &token_temp)[0]));//strand 
				ret_features[xk1].sorted_order = xk1;
				strtok_r(NULL,"\t",&token_temp);	// "frame"
				char * extra_attrs = strtok_r(NULL,"\t",&token_temp);	// name_1 "val1"; name_2 "val2"; ... 
				if(extra_attrs && (strlen(extra_attrs)>2))
				{
					int attr_val_len = GTF_extra_column_value(extra_attrs , global_context -> gene_id_column , feature_name_tmp, FEATURE_NAME_LENGTH);
					if(attr_val_len>0) is_gene_id_found=1;
			//		printf("V=%s\tR=%d\n", extra_attrs , attr_val_len);
				}

				if(is_gene_id_found)
				{
				}
				else
				{
					if(!is_GFF_warned)
					{
						int ext_att_len = strlen(extra_attrs);
						if(extra_attrs[ext_att_len-1] == '\n') extra_attrs[ext_att_len-1] =0;
						SUBREADprintf("\nWarning: failed to find the gene identifier attribute in the 9th column of the provided GTF file.\nThe specified gene identifier attribute is '%s' \nThe attributes included in your GTF annotation are '%s' \n\n",  global_context -> gene_id_column, extra_attrs);
					}
					is_GFF_warned++;
				}

				int feature_name_len = strlen(feature_name_tmp);
				if(feature_name_len > FEATURE_NAME_LENGTH) feature_name_tmp[FEATURE_NAME_LENGTH -1 ] = 0;
				ret_features[xk1].feature_name_pos = unistr_cpy(global_context, (char *)feature_name_tmp, feature_name_len);

				int chro_name_len = strlen(seq_name);
				if(chro_name_len > CHROMOSOME_NAME_LENGTH) seq_name[CHROMOSOME_NAME_LENGTH -1 ] = 0;
				unsigned int chro_name_pos = unistr_cpy(global_context, (char *)seq_name, chro_name_len);
				global_context -> longest_chro_name = max(chro_name_len, global_context -> longest_chro_name);

				ret_features[xk1].chro_name_pos_delta = chro_name_pos - ret_features[xk1].feature_name_pos;

				int bin_location = ret_features[xk1].start / REVERSE_TABLE_BUCKET_LENGTH;
				fc_chromosome_index_info * chro_stab = HashTableGet(chro_name_table, seq_name);
				if(!chro_stab -> reverse_table_start_index)
				{
					chro_stab -> reverse_table_start_index = malloc(sizeof(int) *( chro_stab->chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH +2));
					memset(chro_stab -> reverse_table_start_index, 0 , sizeof(int) *( chro_stab->chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH +2));
				}
				chro_stab -> reverse_table_start_index[bin_location]++;

				if(global_context -> do_junction_counting)
					register_junc_feature(global_context , feature_name_tmp, seq_name, ret_features[xk1].start, ret_features[xk1].end);

				xk1++;
			}
		}
	}
	fclose(fp);
	free(file_line);

	(*loaded_features) = ret_features;
	global_context -> exontable_nchrs = (int)chro_name_table-> numOfElements;
	global_context -> exontable_chro_table = chro_name_table;

	print_in_box(80,0,0,"   Features : %d\n", features);
	if(features < 1)
	{
		print_in_box(80,0,0,"WARNING no features were loaded in format %s.", file_type == FILE_TYPE_GTF?"GTF":"SAF");
		print_in_box(80,0,0,"        annotation format can be specified using '-F'.");
	}
	return features;
}

int find_or_insert_gene_name(fc_thread_global_context_t * global_context, unsigned char * feature_name)
{
	HashTable * genetable = global_context -> gene_name_table;

	long long int gene_number = HashTableGet(genetable, feature_name) - NULL;
	if(gene_number>0)
		return gene_number-1;
	else
	{
		gene_number = genetable -> numOfElements; 
		HashTablePut(genetable, feature_name, NULL+gene_number+1);
		global_context -> gene_name_array[gene_number] = feature_name;
			// real memory space of feature_name is in the "loaded_features" data structure.
			// now we only save its pointer.

		return gene_number;
	}
}

void register_reverse_table(int block_no, long this_block_min_start, long this_block_max_end, fc_chromosome_index_info * chro_inf)
{

	unsigned int reversed_bucket_start = this_block_min_start /  REVERSE_TABLE_BUCKET_LENGTH;
	unsigned int reversed_bucket_end = this_block_max_end / REVERSE_TABLE_BUCKET_LENGTH;
	assert(this_block_min_start <= this_block_max_end);
	assert(reversed_bucket_end < chro_inf -> chro_possible_length);
	int x1;
	for(x1 = reversed_bucket_start; x1 <= reversed_bucket_end; x1++)
	{
		chro_inf->reverse_table_start_index[x1] = min(chro_inf->reverse_table_start_index[x1], block_no);
		//chro_inf->reverse_table_end_index[x1] = max(chro_inf->reverse_table_end_index[x1], block_no+1);
	}

}

void feature_merge(void * arrv, int start, int items, int items2)
{

	void ** arr = (void **) arrv;

	long * ret_start = (long *) arr[0];
	long * ret_end = (long *) arr[1];
	unsigned char * ret_strand = (unsigned char *) arr[2];
	int * ret_entyrez = (int *) arr[3];
	fc_feature_info_t ** old_info_ptr = (fc_feature_info_t **) arr[4];

	int total_items = items+items2;
	long * tmp_start = malloc(sizeof(long) * total_items);
	long * tmp_end = malloc(sizeof(long) * total_items);
	unsigned char * tmp_strand = malloc(sizeof(char) * total_items);
	int * tmp_entyrez = malloc(sizeof(int) * total_items);
	fc_feature_info_t ** tmp_info_ptr = malloc(sizeof(fc_feature_info_t*) * total_items);

	int read_1_ptr = start;
	int read_2_ptr = start+items;
	int write_ptr;

	for(write_ptr=0; write_ptr<total_items; write_ptr++)
	{
		if((read_1_ptr >= start+items)||(read_2_ptr < start+total_items && ret_start[read_1_ptr] >= ret_start[read_2_ptr]))
		{
			tmp_start[write_ptr] = ret_start[read_2_ptr];
			tmp_end[write_ptr] = ret_end[read_2_ptr];
			tmp_strand[write_ptr] = ret_strand[read_2_ptr];
			tmp_entyrez[write_ptr] = ret_entyrez[read_2_ptr];
			tmp_info_ptr[write_ptr] = old_info_ptr[read_2_ptr];
			read_2_ptr++;
		}
		else
		{
			tmp_start[write_ptr] = ret_start[read_1_ptr];
			tmp_end[write_ptr] = ret_end[read_1_ptr];
			tmp_strand[write_ptr] = ret_strand[read_1_ptr];
			tmp_entyrez[write_ptr] = ret_entyrez[read_1_ptr];
			tmp_info_ptr[write_ptr] = old_info_ptr[read_1_ptr];
			read_1_ptr++;
		}
	}

	memcpy(ret_start+ start, tmp_start, sizeof(long) * total_items);
	memcpy(ret_end+ start, tmp_end, sizeof(long) * total_items);
	memcpy(ret_strand+ start, tmp_strand, sizeof(char) * total_items);
	memcpy(ret_entyrez+ start, tmp_entyrez, sizeof(int) * total_items);
	memcpy(old_info_ptr+ start, tmp_info_ptr, sizeof(fc_feature_info_t*) * total_items);

	free(tmp_start);
	free(tmp_end);
	free(tmp_strand);
	free(tmp_entyrez);
	free(tmp_info_ptr);
}


int feature_sort_compare(void * arrv, int l, int r)
{
	void ** arr = (void **) arrv;
	long * ret_start = (long *)arr[0];
	long ll = ret_start[l];
	long rl = ret_start[r];

	if(ll==rl) return 0;
	else if(ll>rl) return 1;
	else return -1;
}

void feature_sort_exchange(void * arrv, int l, int r)
{
	void ** arr = (void **) arrv;
	long tmp;
	fc_feature_info_t * tmpptr;

	long * ret_start = (long *) arr[0];
	long * ret_end = (long *) arr[1];
	unsigned char * ret_strand = (unsigned char *) arr[2];
	int * ret_entyrez = (int *) arr[3];
	fc_feature_info_t ** old_info_ptr = (fc_feature_info_t **) arr[4];

	
	tmp = ret_start[r];
	ret_start[r]=ret_start[l];
	ret_start[l]=tmp;

	tmp = ret_end[r];
	ret_end[r]=ret_end[l];
	ret_end[l]=tmp;

	tmp = ret_strand[r];
	ret_strand[r]=ret_strand[l];
	ret_strand[l]=tmp;

	tmp = ret_entyrez[r];
	ret_entyrez[r]=ret_entyrez[l];
	ret_entyrez[l]=tmp;

	tmpptr = old_info_ptr[r];
	old_info_ptr[r]=old_info_ptr[l];
	old_info_ptr[l]=tmpptr;

}



void sort_feature_info(fc_thread_global_context_t * global_context, unsigned int features, fc_feature_info_t * loaded_features, char *** sorted_chr_names, int ** sorted_entrezid, long ** sorted_start, long ** sorted_end, unsigned char ** sorted_strand, char ** anno_chr_2ch, char *** anno_chrs, long ** anno_chr_head, long ** block_end_index, long ** block_min_start_pos, long ** block_max_end_pos)
{
	unsigned int chro_pnt;
	unsigned int xk1,xk2;
	int * ret_entrez = malloc(sizeof(int) * features);
	long * ret_start = malloc(sizeof(long) * features);
	long * ret_end = malloc(sizeof(long) * features);
	int current_block_buffer_size = 2000;

	long * ret_block_end_index = malloc(sizeof(long) * current_block_buffer_size);
	long * ret_block_min_start = malloc(sizeof(long) * current_block_buffer_size);
	long * ret_block_max_end = malloc(sizeof(long) * current_block_buffer_size);
	unsigned char * ret_strand = malloc(features);
	char ** ret_char_name = malloc(sizeof(void *) * features);
	fc_feature_info_t ** old_info_ptr = malloc(sizeof(void *) * features);
	(*anno_chrs) = malloc(sizeof(void *) * global_context -> exontable_nchrs);
	(*anno_chr_head) = malloc(sizeof(long) * global_context -> exontable_nchrs);
	(*anno_chr_2ch) = malloc(sizeof(char) * global_context -> exontable_nchrs*2); 
	unsigned int * chro_feature_ptr = calloc(sizeof(int) * global_context -> exontable_nchrs,1);
	fc_chromosome_index_info ** tmp_chro_info_ptrs = malloc(global_context -> exontable_nchrs * sizeof(fc_chromosome_index_info *));

	global_context -> gene_name_array = malloc(sizeof(char *) * features);	// there should be much less identical names.
	global_context -> gene_name_table = HashTableCreate(5000);
	HashTableSetHashFunction(global_context -> gene_name_table, HashTableStringHashFunction);
	HashTableSetKeyComparisonFunction(global_context -> gene_name_table, fc_strcmp);

	// init start positions of each chromosome block.
	if(1)
	{
		KeyValuePair * cursor;
		int bucket;
		unsigned int sum_ptr = 0;
		for(bucket=0; bucket < global_context -> exontable_chro_table  -> numOfBuckets; bucket++)
		{
			cursor = global_context -> exontable_chro_table -> bucketArray[bucket];
			while (1)
			{
				if (!cursor) break;
				fc_chromosome_index_info * tmp_chro_inf = cursor -> value;
				cursor = cursor->next;
				//tmp_chro_inf -> reverse_table_end_index = calloc(sizeof(int), tmp_chro_inf->chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH +2);
				chro_feature_ptr [tmp_chro_inf -> chro_number] = tmp_chro_inf -> chro_features;
				tmp_chro_info_ptrs[tmp_chro_inf -> chro_number] = tmp_chro_inf;
			}
		}

		for(xk1 = 0; xk1 < global_context -> exontable_nchrs; xk1++)
		{
			unsigned int tmpv = chro_feature_ptr[xk1];
			chro_feature_ptr[xk1] = sum_ptr;
			tmp_chro_info_ptrs[xk1] -> chro_feature_table_start = sum_ptr;
		//		printf("SII=%u  +  %u\n", sum_ptr, tmpv);
			sum_ptr += tmpv;
		}

	}
	int current_block_id = 0, sort_i = 0;

	(*sorted_chr_names) = ret_char_name;
	(*sorted_entrezid) = ret_entrez;
	(*sorted_start) = ret_start;
	(*sorted_end) = ret_end;
	(*sorted_strand) = ret_strand;
	int curr_chro_number = 0;

	for(chro_pnt=0; chro_pnt < features; chro_pnt++)
	{
		char * this_chro_name = global_context -> unistr_buffer_space + loaded_features[chro_pnt].feature_name_pos + loaded_features[chro_pnt].chro_name_pos_delta;
		fc_chromosome_index_info * this_chro_info = HashTableGet(global_context -> exontable_chro_table , this_chro_name);
		assert(this_chro_info);
		unsigned int this_chro_number = this_chro_info -> chro_number;
		unsigned int this_chro_table_ptr = chro_feature_ptr[this_chro_number];

		ret_char_name[this_chro_table_ptr] = this_chro_name;// (char *)loaded_features[chro_pnt].chro;
		ret_entrez[this_chro_table_ptr] = find_or_insert_gene_name(global_context, (unsigned char *)(global_context -> unistr_buffer_space + loaded_features[chro_pnt].feature_name_pos));
		ret_start[this_chro_table_ptr] = loaded_features[chro_pnt].start;
		ret_end[this_chro_table_ptr] = loaded_features[chro_pnt].end;
		ret_strand[this_chro_table_ptr] = loaded_features[chro_pnt].is_negative_strand;
		old_info_ptr[this_chro_table_ptr] = &loaded_features[chro_pnt];

		chro_feature_ptr[this_chro_number]++;
	}

	for(xk1 = 0; xk1 < global_context -> exontable_nchrs; xk1++)
	{
		fc_chromosome_index_info * tmp_chro_inf = tmp_chro_info_ptrs[xk1];
		int bins_in_chr = ( tmp_chro_inf->chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH +2);
		short * features_per_block_bins = malloc(sizeof(short)*bins_in_chr);
		for(xk2=0; xk2<bins_in_chr; xk2++)
		{
			features_per_block_bins[xk2] = max(1,min(1000,(int)(0.9999999+sqrt(tmp_chro_inf -> reverse_table_start_index[xk2]))));
			//printf("CHR%d : SQR[%d]=%d (%d)\n",  tmp_chro_inf -> chro_number,xk2, features_per_block_bins[xk2], tmp_chro_inf -> reverse_table_start_index[xk2] );
		}

		memset(tmp_chro_inf -> reverse_table_start_index, 0xff, sizeof(int) *bins_in_chr);

		tmp_chro_inf -> chro_block_table_start = current_block_id; 
		unsigned int this_block_items = 0;
		long this_block_min_start = 0x7fffffff, this_block_max_end = 0;
		unsigned int this_chro_tab_end =  tmp_chro_inf -> chro_features + tmp_chro_inf -> chro_feature_table_start;

		void * in_array[5];
		in_array[0] = ret_start + tmp_chro_inf -> chro_feature_table_start; 
		in_array[1] = ret_end + tmp_chro_inf -> chro_feature_table_start; 
		in_array[2] = ret_strand + tmp_chro_inf -> chro_feature_table_start; 
		in_array[3] = ret_entrez + tmp_chro_inf -> chro_feature_table_start; 
		in_array[4] = old_info_ptr + tmp_chro_inf -> chro_feature_table_start; 

		merge_sort(in_array, this_chro_tab_end - tmp_chro_inf -> chro_feature_table_start, feature_sort_compare, feature_sort_exchange, feature_merge);

		for(sort_i = tmp_chro_inf -> chro_feature_table_start; sort_i< this_chro_tab_end ; sort_i++)
		{
			// NOW THE FEATURES (ret_start, ret_end, ret_strand, ret_entrez, old_info_ptr) ARE ALL SORTED!
			//printf("NT=%lu\tCHRO=%d\n", ret_start[sort_i], tmp_chro_inf->chro_number);
			old_info_ptr[sort_i]->sorted_order = sort_i;

			int feature_bin_location = ret_start[sort_i] / REVERSE_TABLE_BUCKET_LENGTH;
			int block_bin_location = this_block_min_start / REVERSE_TABLE_BUCKET_LENGTH;

			if(this_block_items && (this_block_items > features_per_block_bins[block_bin_location] || feature_bin_location != block_bin_location))//global_context -> feature_block_size)
			{

				if(current_block_id >= current_block_buffer_size - 1)
				{
					current_block_buffer_size *= 1.3;
					ret_block_min_start = realloc(ret_block_min_start, sizeof(long)*current_block_buffer_size);
					ret_block_max_end = realloc(ret_block_max_end, sizeof(long)*current_block_buffer_size);
					ret_block_end_index = realloc(ret_block_end_index, sizeof(long)*current_block_buffer_size);
				}


				ret_block_end_index[current_block_id] = sort_i;	// FIRST UNWANTED ID
				ret_block_min_start[current_block_id] = this_block_min_start;
				ret_block_max_end[current_block_id] = this_block_max_end;
				register_reverse_table(current_block_id, this_block_min_start, this_block_max_end, tmp_chro_inf);
				//printf("B=%d; ST=%ld, END=%ld, ITM=%d\n", current_block_id, this_block_min_start, this_block_max_end, this_block_items);
				current_block_id++;
				this_block_max_end = 0;
				this_block_items = 0;
				this_block_min_start = 0x7fffffff;
			}

			this_block_max_end = max(this_block_max_end, ret_end[sort_i]);
			this_block_min_start = min(this_block_min_start, ret_start[sort_i]);
			this_block_items ++;
		
		}
		if(this_block_items)
		{
			if(current_block_id >= current_block_buffer_size)
			{
				current_block_buffer_size *= 1.3;
				ret_block_min_start = realloc(ret_block_min_start, sizeof(long)*current_block_buffer_size);
				ret_block_max_end = realloc(ret_block_max_end, sizeof(long)*current_block_buffer_size);
				ret_block_end_index = realloc(ret_block_end_index, sizeof(long)*current_block_buffer_size);
			}

			ret_block_end_index[current_block_id] = this_chro_tab_end;	// FIRST UNWANTED ID
			ret_block_min_start[current_block_id] = this_block_min_start;
			ret_block_max_end[current_block_id] = this_block_max_end;
			register_reverse_table(current_block_id, this_block_min_start, this_block_max_end, tmp_chro_inf);
			current_block_id++;
		}

		(*anno_chr_head) [curr_chro_number] = current_block_id; 
		tmp_chro_inf -> chro_block_table_end = current_block_id; 
		free(features_per_block_bins);
	}

	(*block_end_index) = ret_block_end_index;
	(*block_min_start_pos) = ret_block_min_start;
	(*block_max_end_pos) = ret_block_max_end;

	//print_in_box(80, 0,0,"The %u features are sorted.\n", sort_i);
	free(old_info_ptr);
	free(tmp_chro_info_ptrs);
	free(chro_feature_ptr);
}

int strcmp_slash(char * s1, char * s2)
{
	char nch;
	while(0!=(nch = *(s1++))){
		if(nch == '/') break;
		if(nch != (*s2)) return 1;
		s2++;
	}
	return nch != *s2;
}

#define NH_FRACTION_INT 65536

unsigned int calculate_multi_overlap_fraction(fc_thread_global_context_t * global_context, unsigned int fixed_fractional_count, int maximum_total_count){
	if(global_context -> use_fraction_multi_mapping) return fixed_fractional_count / maximum_total_count;
	else return fixed_fractional_count;
}

unsigned int calc_fixed_fraction(int nh){
	if(nh==1) return NH_FRACTION_INT;
	else if(nh == 2) return NH_FRACTION_INT>>1;
	else return NH_FRACTION_INT / nh; 
}


int calc_float_fraction(read_count_type_t score, read_count_type_t * integer_count, double * float_count){
	if(score % NH_FRACTION_INT == 0){
		(*integer_count) = score / NH_FRACTION_INT;
		return 0;
	}else{
		(*float_count) = score * 1./NH_FRACTION_INT;
		return 1;
	}
}


void print_read_wrapping(char * rl, int is_second){
	int refill_spaces = 3;

	int read_length = 0, x1 = 0, spaces=0;

	for(x1 = 0; x1 < 3100; x1++){
		if(rl[x1]==0 && rl[x1+1]==0)break;
		if(rl[x1]=='0' || rl[x1]=='\t') spaces++;
		read_length ++;
	}

	char *out_buf1 = malloc(read_length + spaces * refill_spaces + 1), out_buf2[100];
	int ox=0;

	for(x1 = 0; x1 < 3000; x1++){
		if(rl[x1]=='\n' || (rl[x1]==0 && rl[x1+1]==0)){
			out_buf1[ox]=0;
			break;
		} else if((rl[x1]==0 && rl[x1+1]!=0) || rl[x1] == '\t'){
			int x2;
			for(x2 = 0; x2 < refill_spaces ; x2++){
				out_buf1[ox]=' ';
				ox++;
			}
		} else {
			out_buf1[ox]=rl[x1];
			ox++;
		}
	}
	out_buf1[ox] = 0;

	x1=0;

	while(1){
		int x2;
		for(x2 = 0; x2 < 67 ; x2 ++){
			char nch = out_buf1[x1];
			if(nch == 0) break;
			out_buf2[x2] = nch;
			x1++;
		}
		out_buf2[x2] = 0;

		print_in_box(80,0,PRINT_BOX_NOCOLOR_FOR_COLON,"      %s", out_buf2);
		if(out_buf1[x1] == 0)break;
	}

	free(out_buf1);

}


void process_pairer_reset(void * pairer_vp){
	SAM_pairer_context_t * pairer = (SAM_pairer_context_t *) pairer_vp;
	fc_thread_global_context_t * global_context = (fc_thread_global_context_t * )pairer -> appendix1;
	if(global_context -> sambam_chro_table) free(global_context -> sambam_chro_table);
	global_context -> sambam_chro_table = NULL;
	global_context -> sambam_chro_table_items = 0;

	int xk1, xk2;
	for(xk1=0; xk1<global_context-> thread_number; xk1++)
	{
		for(xk2=0; xk2<global_context -> exontable_exons; xk2++)
		{
			global_context -> thread_contexts[xk1].count_table[xk2] = 0;
		}

		global_context -> thread_contexts[xk1].all_reads = 0;
		global_context -> thread_contexts[xk1].nreads_mapped_to_exon = 0;
		global_context -> thread_contexts[xk1].unpaired_fragment_no = 0;



		global_context -> thread_contexts[xk1].read_counters.unassigned_ambiguous = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_nofeatures = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_unmapped = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_mappingquality = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_fragmentlength = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_chimericreads = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_multimapping = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_secondary = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_junction_condition = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_duplicate = 0;
		global_context -> thread_contexts[xk1].read_counters.assigned_reads = 0;
	}

	if(global_context -> SAM_output_fp){
		ftruncate(fileno(global_context -> SAM_output_fp), 0);
		fseek(global_context -> SAM_output_fp, 0 , SEEK_SET);
	}
}

int process_pairer_header (void * pairer_vp, int thread_no, int is_text, unsigned int items, char * bin, unsigned int bin_len){


	SAM_pairer_context_t * pairer = (SAM_pairer_context_t *) pairer_vp;
	fc_thread_global_context_t * global_context = (fc_thread_global_context_t * )pairer -> appendix1;

	//SUBREADprintf("ENTER PROCESS (THRD %d): IS_TXT=%d,  ITEMS = %d, CURRENT_ITEMS=%d\n", thread_no, is_text, items, global_context -> sambam_chro_table_items);
	pthread_spin_lock(&global_context -> sambam_chro_table_lock);

	if( !is_text ){
		if(global_context -> sambam_chro_table)
			global_context -> sambam_chro_table = delay_realloc(global_context -> sambam_chro_table, global_context -> sambam_chro_table_items * sizeof(SamBam_Reference_Info), (items + global_context -> sambam_chro_table_items) * sizeof(SamBam_Reference_Info));
		else global_context -> sambam_chro_table = malloc(items * sizeof(SamBam_Reference_Info));

		int x1, bin_ptr = 0;
		for(x1 =  global_context -> sambam_chro_table_items; x1 <  global_context -> sambam_chro_table_items+items; x1++){
			int l_name;
			memcpy(&l_name, bin + bin_ptr, 4);
			assert(l_name < MAX_CHROMOSOME_NAME_LEN);
			bin_ptr += 4;
			memcpy(global_context -> sambam_chro_table[x1].chro_name ,  bin + bin_ptr, l_name);
			//SUBREADprintf("The %d-th is '%s'\n", x1, global_context -> sambam_chro_table[x1].chro_name);
			bin_ptr += l_name;
			memcpy(&global_context -> sambam_chro_table[x1].chro_length ,  bin + bin_ptr, 4);
			bin_ptr += 4;
		}
		global_context -> sambam_chro_table_items += items;
	}
	pthread_spin_unlock(&global_context -> sambam_chro_table_lock);
	return 0;
}

void process_line_buffer(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context, char * bin1, char * bin2);

void make_dummy(char * rname, char * bin1, char * out_txt2,  SamBam_Reference_Info * sambam_chro_table){
	char * tmptr = NULL;

	//SUBREADprintf("S=%s  ", rname);
	char * realname = strtok_r(rname, "\027", &tmptr);
	//int len_name = strlen(realname);
	int r1_chro = atoi(strtok_r(NULL, "\027", &tmptr));
	int r1_pos = atoi(strtok_r(NULL, "\027", &tmptr));
	int r2_chro = atoi(strtok_r(NULL, "\027", &tmptr));
	int r2_pos = atoi(strtok_r(NULL, "\027", &tmptr));
	int HItag = atoi(strtok_r(NULL, "\027", &tmptr));
	int mate_FLAG = 0;
	memcpy(&mate_FLAG, bin1 + 16, 4);
	mate_FLAG = 0xffff&(mate_FLAG >>16);
	int mate_tlen = 0;
	memcpy(&mate_tlen, bin1 + 32, 4);

	if(r1_chro<0) r1_pos=-1;
	if(r2_chro<0) r2_pos=-1;

	int my_chro = (mate_FLAG&0x40)? r2_chro : r1_chro;
	int my_pos = (mate_FLAG&0x40)? r2_pos : r1_pos;
	int mate_chro = (mate_FLAG&0x40)? r1_chro : r2_chro;
	int mate_pos = (mate_FLAG&0x40)? r1_pos : r2_pos;

	//int bin_mq_nl = (len_name+1);
	int my_flag = (mate_FLAG&0x40)? 0x80:0x40;
	my_flag |= 1;
	if(mate_FLAG & 8)my_flag |=4;
	if(mate_FLAG & 4)my_flag |=8;
	if(mate_FLAG & 0x10) my_flag |= 0x20;
	if(mate_FLAG & 0x20) my_flag |= 0x10;

	char HItagStr[20];
	if(HItag>=0){
		sprintf(HItagStr, "\tHI:i:%d", HItag);
	}else{
		HItagStr[0]=0;
	}

	char * my_chro_str = "*";
	if(my_chro >= 0) my_chro_str = sambam_chro_table[my_chro].chro_name;

	char * mate_chro_str = "*";
	if(mate_chro >= 0) mate_chro_str = sambam_chro_table[mate_chro].chro_name;

	sprintf(out_txt2, "%s\t%d\t%s\t%d\t0\t*\t%s\t%d\t0\tN\tI\t%s", realname, my_flag, my_chro_str, max(0, my_pos),
		mate_chro_str, max(0,mate_pos), HItagStr);
}


void convert_bin_to_read(char * bin, char * txt, SamBam_Reference_Info * sambam_chro_table){
	unsigned int block_len;
	memcpy(&block_len, bin, 4);
	int ref_id;
	memcpy(&ref_id, bin + 4, 4);
	int pos;
	memcpy(&pos, bin + 8, 4);
	unsigned int bin_mq_nl;
	memcpy(&bin_mq_nl, bin + 12, 4);
	unsigned int flag_nc;
	memcpy(&flag_nc, bin + 16, 4);
	int l_seq;
	memcpy(&l_seq, bin + 20, 4);
	int next_refID;
	memcpy(&next_refID, bin + 24, 4);
	int next_pos;
	memcpy(&next_pos, bin + 28, 4);
	int tlen;
	memcpy(&tlen, bin + 32, 4);

	int txt_ptr = 0;
	int l_read_name = bin_mq_nl & 0xff;
	memcpy(txt , bin + 36, l_read_name);
	txt_ptr += l_read_name - 1;
	txt_ptr += sprintf(txt+txt_ptr, "\t%d", flag_nc >> 16);
	if(ref_id < 0){
		strcpy(txt+txt_ptr, "\t*\t0\t0");
		txt_ptr += 6;
	}else 	txt_ptr += sprintf(txt+txt_ptr, "\t%s\t%d\t%d", sambam_chro_table[ref_id].chro_name, pos + 1, (bin_mq_nl >> 8 & 0xff));

	int cigar_ops = flag_nc & 0xffff;
	if(cigar_ops < 1){
		strcpy(txt+txt_ptr, "\t*");
		txt_ptr += 2;
	}else{
		int x1;
		strcpy(txt+txt_ptr, "\t");
		txt_ptr++;
		for(x1=0; x1 < cigar_ops; x1++){
			unsigned int cigar_sec;
			memcpy(&cigar_sec, bin + 36 + l_read_name + 4 * x1 , 4);
			txt_ptr += sprintf(txt+txt_ptr, "%u%c", cigar_sec >> 4 , cigar_op_char( cigar_sec & 15 ));
		}
	}

	if(next_refID < 0)
		txt_ptr += sprintf(txt+txt_ptr, "\t*\t0\t%d", tlen);
	else 	txt_ptr += sprintf(txt+txt_ptr, "\t%s\t%d\t%d", sambam_chro_table[next_refID].chro_name, next_pos + 1, tlen);
	strcpy(txt+txt_ptr, "\tN\tI");
	txt_ptr += 4;

	int bin_ptr = 36 + l_read_name + 4 * cigar_ops + l_seq + (l_seq+1)/2;

	while(bin_ptr < block_len + 4){
		char tag_name[3];
		tag_name[0]=bin[bin_ptr];
		tag_name[1]=bin[bin_ptr+1];
		tag_name[2]=0;

		char tagtype = bin[bin_ptr+2];
		int delta = 0;
		int tmpi = 0;
		if(tagtype == 'i' || tagtype == 'I'){
			delta = 4;
			memcpy(&tmpi, bin + bin_ptr + 3, 4);
			txt_ptr += sprintf(txt+txt_ptr, "\t%s:i:%d", tag_name,tmpi);
		}else if(tagtype == 's' || tagtype == 'S'){
			delta = 2;
			memcpy(&tmpi, bin + bin_ptr + 3, 2);
			txt_ptr += sprintf(txt+txt_ptr, "\t%s:i:%d", tag_name,tmpi);
		}else if(tagtype == 'c' || tagtype == 'C'){
			delta = 1;
			memcpy(&tmpi, bin + bin_ptr + 3, 1);
			txt_ptr += sprintf(txt+txt_ptr, "\t%s:i:%d", tag_name,tmpi);
		}else if(tagtype == 'A'){
			delta = 1;
			txt_ptr += sprintf(txt+txt_ptr, "\t%s:%c:%c", tag_name, tagtype, *(bin + bin_ptr + 3));
		}else if(tagtype == 'f')
			delta = 4;
		else if(tagtype == 'Z' ||tagtype == 'H'){
			txt_ptr += sprintf(txt+txt_ptr, "\t%s:%c", tag_name, tagtype);
			while(bin[bin_ptr + 3+delta]){
				*(txt+txt_ptr) = bin[bin_ptr + 3+delta];
				txt_ptr ++;
				delta ++;
			}
			*(txt+txt_ptr) = 0;
		}else if(tagtype == 'B'){
			char celltype = bin[bin_ptr + 4];
			int cellitems ;
			memcpy(&cellitems, bin + bin_ptr + 5, 4);
			int celldelta = 1;
			if(celltype == 's' || celltype == 'S') celldelta = 2;
			else if(celltype == 'i' || celltype == 'I' || celltype == 'f') celldelta = 4;
			delta = cellitems * celldelta;
		}
		bin_ptr += 3 + delta;
	}
}
int reverse_flag(int mf){
	int ret = mf & 3;
	if(mf & 4) ret |= 8;
	if(mf & 8) ret |= 4;

	if(mf & 0x10) ret |= 0x20;
	if(mf & 0x20) ret |= 0x10;

	if(mf & 0x40) ret |= 0x80;
	if(mf & 0x80) ret |= 0x40;
	return ret;
}

#define MAXIMUM_INSERTION_IN_SECTION 8

typedef struct {
	char * chro;
	unsigned int start_pos;
	unsigned int chromosomal_length;
	short insertions;
	unsigned int insertion_start_pos[ MAXIMUM_INSERTION_IN_SECTION ];
	unsigned short insertion_lengths[ MAXIMUM_INSERTION_IN_SECTION ];
} CIGAR_interval_t;


int calc_total_frag_one_len(CIGAR_interval_t * intvs, int intvn){
	int ret = 0, x1;
	for(x1 = 0; x1 < intvn; x1++){
		int x2;
		for(x2 = 0; x2 < intvs[x1].insertions; x2++) ret += intvs[x1].insertion_lengths[x2];
		ret += intvs[x1].chromosomal_length;
	}
	return ret;
}

int calc_total_has_overlap(unsigned int r1_start, unsigned int r1_end, unsigned int r2_start, unsigned int r2_end, unsigned int * overlap_start, unsigned int * overlap_end){
	if((r1_start <= r2_start && r1_end > r2_start) || (r2_start <= r1_start && r2_end > r1_start) ){
		(*overlap_start) = max( r1_start, r2_start );
		(*overlap_end) = min( r1_end, r2_end );
		return 1;
	}
	return 0;
}

int calc_total_frag_len( fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context, CIGAR_interval_t * CIGAR_intervals_R1, int CIGAR_intervals_R1_sections, CIGAR_interval_t * CIGAR_intervals_R2, int CIGAR_intervals_R2_sections, char * read_name){
	if     ( CIGAR_intervals_R1_sections == 0 && CIGAR_intervals_R2_sections > 0) return calc_total_frag_one_len( CIGAR_intervals_R2,CIGAR_intervals_R2_sections );
	else if( CIGAR_intervals_R1_sections  > 0 && CIGAR_intervals_R2_sections== 0) return calc_total_frag_one_len( CIGAR_intervals_R1,CIGAR_intervals_R1_sections );
	else if( CIGAR_intervals_R1_sections == 0 && CIGAR_intervals_R2_sections== 0) return 0;

	if(CIGAR_intervals_R1_sections > 0 && CIGAR_intervals_R2_sections > 0 && strcmp(CIGAR_intervals_R1[0].chro, CIGAR_intervals_R2[0].chro )!=0 )
		// two reads are from different chromosomes
		return calc_total_frag_one_len( CIGAR_intervals_R2,CIGAR_intervals_R2_sections ) + calc_total_frag_one_len( CIGAR_intervals_R1,CIGAR_intervals_R1_sections );

	
	unsigned int merged_section_count = 0;
	unsigned int merged_section_start_positions[ MAXIMUM_INSERTION_IN_SECTION * 3 ];
	unsigned short merged_section_lengths[ MAXIMUM_INSERTION_IN_SECTION * 3 ];
	unsigned char merged_section_belong_R1[ MAXIMUM_INSERTION_IN_SECTION * 3 ];
	unsigned char merged_section_belong_R2[ MAXIMUM_INSERTION_IN_SECTION * 3 ];
	unsigned int merged_section_indel_counts[ MAXIMUM_INSERTION_IN_SECTION * 3 ];
	unsigned int merged_section_indel_positions[ MAXIMUM_INSERTION_IN_SECTION * 3 ][ MAXIMUM_INSERTION_IN_SECTION ];
	unsigned short merged_section_indel_lengths[ MAXIMUM_INSERTION_IN_SECTION * 3 ][ MAXIMUM_INSERTION_IN_SECTION ];

	int R1_i = 0 , R2_i = 0;
	while (1){
		//SUBREADprintf("FRAGDEBUG %s : %d < %d & %d < %d; MC=%d; INS1=%d INS2=%d\n", read_name, R1_i,CIGAR_intervals_R1_sections,R2_i,CIGAR_intervals_R2_sections, merged_section_count, CIGAR_intervals_R1[R1_i].insertions, CIGAR_intervals_R2[R2_i].insertions);
		if( R1_i >= CIGAR_intervals_R1_sections &&  R2_i >= CIGAR_intervals_R2_sections ) break;

		if( R1_i < CIGAR_intervals_R1_sections && R2_i < CIGAR_intervals_R2_sections){
			// see if R1 and R2 overlap
			// if not: add R2 to specific sction; R2_i ++
			// elif overlap: add the R1 first_half and/or R2 first_half or zero to specific section, and add overlapping part to overlapping section; DO NOT add the second specific half!
			// 	if R1_end > R2_end: R1_section_start = overlapping_end; R2_i ++
			// 	elif R2_end > R1_end: R2_section_start = overlapping_end; R1_i ++ 
			// 	elif R2_end == R1_end: R1_i++; R2_i++

			unsigned int overlapping_start= 0 ,  overlapping_end = 0;

			int is_r1r2_overlap = 0;

			is_r1r2_overlap = calc_total_has_overlap( CIGAR_intervals_R1[R1_i].start_pos, CIGAR_intervals_R1[R1_i].start_pos + CIGAR_intervals_R1[R1_i].chromosomal_length , CIGAR_intervals_R2[R2_i].start_pos, CIGAR_intervals_R2[R2_i].start_pos + CIGAR_intervals_R2[R2_i].chromosomal_length , & overlapping_start , & overlapping_end);

			if( is_r1r2_overlap ){
				if (CIGAR_intervals_R1[R1_i].start_pos > CIGAR_intervals_R2[R2_i].start_pos ){
					//first half_R2 add special
					merged_section_start_positions[merged_section_count] = CIGAR_intervals_R2[R2_i].start_pos;
					merged_section_lengths[merged_section_count] = overlapping_start - CIGAR_intervals_R2[R2_i].start_pos;
					merged_section_belong_R1[merged_section_count]=0;
					merged_section_belong_R2[merged_section_count]=1;

					int indel_i;
					for(indel_i = 0; indel_i < CIGAR_intervals_R2[R2_i].insertions; indel_i++){
						if( CIGAR_intervals_R2[R2_i].insertion_start_pos[indel_i] >= overlapping_start ){
							if(indel_i>0){
								int insmov_i, ins_dist_i = 0;
								for(insmov_i = indel_i ; insmov_i < CIGAR_intervals_R2[R2_i].insertions; insmov_i++){
									CIGAR_intervals_R2[R2_i].insertion_start_pos[ins_dist_i] = CIGAR_intervals_R2[R2_i].insertion_start_pos[insmov_i];
									CIGAR_intervals_R2[R2_i].insertion_lengths[ins_dist_i] = CIGAR_intervals_R2[R2_i].insertion_lengths[insmov_i];
									ins_dist_i++;
								}
								CIGAR_intervals_R2[R2_i].insertions = ins_dist_i;
							}
							break;
						}
						merged_section_indel_positions[merged_section_count][indel_i] = CIGAR_intervals_R2[R2_i].insertion_start_pos[indel_i];
						merged_section_indel_lengths[merged_section_count][indel_i] = CIGAR_intervals_R2[R2_i].insertion_lengths[indel_i];
					}
					merged_section_indel_counts[merged_section_count] = indel_i;

					merged_section_count ++;
	
				}else if( CIGAR_intervals_R1[R1_i].start_pos < CIGAR_intervals_R2[R2_i].start_pos ){
					//first half_R1 add special
					merged_section_start_positions[merged_section_count] = CIGAR_intervals_R1[R1_i].start_pos;
					merged_section_lengths[merged_section_count] = overlapping_start - CIGAR_intervals_R1[R1_i].start_pos;
					merged_section_belong_R1[merged_section_count]=1;
					merged_section_belong_R2[merged_section_count]=0;

					int indel_i;
					for(indel_i = 0; indel_i < CIGAR_intervals_R1[R1_i].insertions; indel_i++){
						if( CIGAR_intervals_R1[R1_i].insertion_start_pos[indel_i] >= overlapping_start ){
							if(indel_i>0){
								int insmov_i, ins_dist_i = 0;
								for(insmov_i = indel_i ; insmov_i < CIGAR_intervals_R1[R1_i].insertions; insmov_i++){
									CIGAR_intervals_R1[R1_i].insertion_start_pos[ins_dist_i] = CIGAR_intervals_R1[R1_i].insertion_start_pos[insmov_i];
									CIGAR_intervals_R1[R1_i].insertion_lengths[ins_dist_i] = CIGAR_intervals_R1[R1_i].insertion_lengths[insmov_i];
									ins_dist_i++;
								}
								CIGAR_intervals_R1[R1_i].insertions = ins_dist_i;
							}
							break;
						}

						merged_section_indel_positions[merged_section_count][indel_i] = CIGAR_intervals_R1[R1_i].insertion_start_pos[indel_i];
						merged_section_indel_lengths[merged_section_count][indel_i] = CIGAR_intervals_R1[R1_i].insertion_lengths[indel_i];
					}
					merged_section_indel_counts[merged_section_count] = indel_i;

					merged_section_count ++;
				}

				merged_section_start_positions[merged_section_count] = overlapping_start;
				merged_section_lengths[merged_section_count] = overlapping_end - overlapping_start;
				merged_section_belong_R1[merged_section_count]=1;
				merged_section_belong_R2[merged_section_count]=1;
				merged_section_indel_counts[merged_section_count] = 0;


				int indel_i_R1 = 0, indel_i_R2 = 0;
				while(1){
					//SUBREADprintf("FRAGDEBUG: CC[%d] = %d ; II1=%d < %d; II2=%d < %d\n", merged_section_count,  merged_section_indel_counts[merged_section_count], indel_i_R1, CIGAR_intervals_R1[R1_i].insertions , indel_i_R2, CIGAR_intervals_R2[R2_i].insertions);

					if( indel_i_R1 >= CIGAR_intervals_R1[R1_i].insertions ||  indel_i_R2 >= CIGAR_intervals_R2[R2_i].insertions ){
						if(indel_i_R1 > 0){
							int insmov_i, ins_dist_i = 0;
							for(insmov_i = indel_i_R1 ; insmov_i < CIGAR_intervals_R1[R1_i].insertions; insmov_i++){
								CIGAR_intervals_R1[R1_i].insertion_start_pos[ins_dist_i] = CIGAR_intervals_R1[R1_i].insertion_start_pos[insmov_i];
								CIGAR_intervals_R1[R1_i].insertion_lengths[ins_dist_i] = CIGAR_intervals_R1[R1_i].insertion_lengths[insmov_i];
								ins_dist_i++;
							}
							CIGAR_intervals_R1[R1_i].insertions = ins_dist_i;
						}
						if(indel_i_R2 > 0){
							int insmov_i, ins_dist_i = 0;
							for(insmov_i = indel_i_R2 ; insmov_i < CIGAR_intervals_R2[R2_i].insertions; insmov_i++){
								CIGAR_intervals_R2[R2_i].insertion_start_pos[ins_dist_i] = CIGAR_intervals_R2[R2_i].insertion_start_pos[insmov_i];
								CIGAR_intervals_R2[R2_i].insertion_lengths[ins_dist_i] = CIGAR_intervals_R2[R2_i].insertion_lengths[insmov_i];
								ins_dist_i++;
							}
							CIGAR_intervals_R2[R2_i].insertions = ins_dist_i;
						}
						break;
					}

					if( CIGAR_intervals_R1[R1_i].insertion_start_pos[indel_i_R1] > CIGAR_intervals_R2[R2_i].insertion_start_pos[indel_i_R2] ) indel_i_R2 ++;
					else if( CIGAR_intervals_R1[R1_i].insertion_start_pos[indel_i_R1] < CIGAR_intervals_R2[R2_i].insertion_start_pos[indel_i_R2]  ) indel_i_R1 ++;
					else{
						if( CIGAR_intervals_R1[R1_i].insertion_lengths[ indel_i_R1 ] == CIGAR_intervals_R2[R2_i].insertion_lengths[ indel_i_R2 ] ){
							merged_section_indel_positions[merged_section_count][ merged_section_indel_counts[merged_section_count] ] = CIGAR_intervals_R1[R1_i].insertion_start_pos[indel_i_R1];
							merged_section_indel_lengths[merged_section_count][ merged_section_indel_counts[merged_section_count] ] = CIGAR_intervals_R1[R1_i].insertion_lengths[indel_i_R1];
							merged_section_indel_counts[merged_section_count] ++;
						}
						indel_i_R2++;
						indel_i_R1++;
					}
				}

				merged_section_count ++;

				// add common

				if(CIGAR_intervals_R1[R1_i].start_pos + CIGAR_intervals_R1[R1_i].chromosomal_length > CIGAR_intervals_R2[R2_i].start_pos + CIGAR_intervals_R2[R2_i].chromosomal_length){
					CIGAR_intervals_R1[R1_i].chromosomal_length -= ( overlapping_end - CIGAR_intervals_R1[R1_i].start_pos );
					CIGAR_intervals_R1[R1_i].start_pos = overlapping_end;
					R2_i ++;
				}else if(CIGAR_intervals_R1[R1_i].start_pos + CIGAR_intervals_R1[R1_i].chromosomal_length < CIGAR_intervals_R2[R2_i].start_pos + CIGAR_intervals_R2[R2_i].chromosomal_length){
					CIGAR_intervals_R2[R2_i].chromosomal_length -= ( overlapping_end - CIGAR_intervals_R2[R2_i].start_pos );
					CIGAR_intervals_R2[R2_i].start_pos = overlapping_end;
					R1_i ++;
				}else{
					R1_i ++;
					R2_i ++;
				}

			}else if(CIGAR_intervals_R1[R1_i].start_pos >  CIGAR_intervals_R2[R2_i].start_pos){
				merged_section_start_positions[merged_section_count] = CIGAR_intervals_R2[R2_i].start_pos;
				merged_section_lengths[merged_section_count] = CIGAR_intervals_R2[R2_i].chromosomal_length;
				merged_section_belong_R1[merged_section_count]=0;
				merged_section_belong_R2[merged_section_count]=1;

				int indel_i;
				for(indel_i = 0; indel_i < CIGAR_intervals_R2[R2_i].insertions; indel_i++){
					merged_section_indel_positions[merged_section_count][indel_i] = CIGAR_intervals_R2[R2_i].insertion_start_pos[indel_i];
					merged_section_indel_lengths[merged_section_count][indel_i] = CIGAR_intervals_R2[R2_i].insertion_lengths[indel_i];
				}
				merged_section_indel_counts[merged_section_count] = CIGAR_intervals_R2[R2_i].insertions;

				merged_section_count ++;
				R2_i ++;
			}else{
				merged_section_start_positions[merged_section_count] = CIGAR_intervals_R1[R1_i].start_pos;
				merged_section_lengths[merged_section_count] = CIGAR_intervals_R1[R1_i].chromosomal_length;
				merged_section_belong_R1[merged_section_count]=1;
				merged_section_belong_R2[merged_section_count]=0;

				int indel_i;
				for(indel_i = 0; indel_i < CIGAR_intervals_R1[R1_i].insertions; indel_i++){
					merged_section_indel_positions[merged_section_count][indel_i] = CIGAR_intervals_R1[R1_i].insertion_start_pos[indel_i];
					merged_section_indel_lengths[merged_section_count][indel_i] = CIGAR_intervals_R1[R1_i].insertion_lengths[indel_i];
				}
				merged_section_indel_counts[merged_section_count] = CIGAR_intervals_R1[R1_i].insertions;

				merged_section_count ++;
				R1_i ++;
			}
		}else if(R1_i < CIGAR_intervals_R1_sections){ 
			// add R1 section to specific section
			// R1_i ++
			merged_section_start_positions[merged_section_count] = CIGAR_intervals_R1[R1_i].start_pos;
			merged_section_lengths[merged_section_count] = CIGAR_intervals_R1[R1_i].chromosomal_length;
			merged_section_belong_R1[merged_section_count]=1;
			merged_section_belong_R2[merged_section_count]=0;

			int indel_i;
			for(indel_i = 0; indel_i < CIGAR_intervals_R1[R1_i].insertions; indel_i++){
				merged_section_indel_positions[merged_section_count][indel_i] = CIGAR_intervals_R1[R1_i].insertion_start_pos[indel_i];
				merged_section_indel_lengths[merged_section_count][indel_i] = CIGAR_intervals_R1[R1_i].insertion_lengths[indel_i];
			}
			merged_section_indel_counts[merged_section_count] = CIGAR_intervals_R1[R1_i].insertions;

			merged_section_count ++;
			R1_i ++;
		}else if(R2_i < CIGAR_intervals_R2_sections){
			merged_section_start_positions[merged_section_count] = CIGAR_intervals_R2[R2_i].start_pos;
			merged_section_lengths[merged_section_count] = CIGAR_intervals_R2[R2_i].chromosomal_length;
			merged_section_belong_R1[merged_section_count]=0;
			merged_section_belong_R2[merged_section_count]=1;

			int indel_i;
			for(indel_i = 0; indel_i < CIGAR_intervals_R2[R2_i].insertions; indel_i++){
				merged_section_indel_positions[merged_section_count][indel_i] = CIGAR_intervals_R2[R2_i].insertion_start_pos[indel_i];
				merged_section_indel_lengths[merged_section_count][indel_i] = CIGAR_intervals_R2[R2_i].insertion_lengths[indel_i];
			}
			merged_section_indel_counts[merged_section_count] = CIGAR_intervals_R2[R2_i].insertions;

			merged_section_count ++;
			R2_i ++;
		}
	}

	int ret = 0, x1, x2;
	for(x1 = 0; x1 < merged_section_count ; x1++){
		ret += merged_section_lengths[x1];
		for(x2 = 0; x2 < merged_section_indel_counts[x1]; x2++)
			ret += merged_section_indel_lengths[x1][x2];
//		SUBREADprintf("FRAGDEBUG %s [%d] : len = %d , indels = %d\n" , read_name, x1, merged_section_lengths[x1] , merged_section_indel_counts[x1]);
	}

	return ret; 
}

void parse_bin(SamBam_Reference_Info * sambam_chro_table, char * bin, char * bin2, char ** read_name, int * flag, char ** chro, long * pos, int * mapq, char ** mate_chro, long * mate_pos, long * tlen, int * is_junction_read, int * cigar_sect, unsigned int * Starting_Chro_Points, unsigned short * Starting_Read_Points, unsigned short * Section_Read_Lengths, char ** ChroNames, char * Event_After_Section, int * NH_value, int max_M, CIGAR_interval_t * intervals_buffer, int * intervals_i){
	int x1, len_of_S1 = 0;
	*cigar_sect = 0;
	*NH_value = 1;
	*flag = 0;
	*is_junction_read = 0;
	assert(bin||bin2);

	if(bin){
		(*read_name) = bin + 36;
		memcpy(flag, bin + 16, 4);
		int cigar_opts = (*flag) & 0xffff;
		(*flag) = (*flag) >> 16;
		int refID, mate_refID;
		memcpy(&refID, bin + 4, 4);
		if(refID >= 0){
			/*if(sambam_chro_table[refID].chro_name < NULL + 0xfffff){
				SUBREADprintf("DANGEROUS: PARSE: chro[%d] = %p, TABLE_PTR=%p\n", refID , sambam_chro_table[refID].chro_name, sambam_chro_table);
			}*/
			(*chro) = sambam_chro_table[refID].chro_name;
		}
		else (*chro) = NULL;
		(*pos) = 0;
		memcpy(pos, bin+8, 4);
		(*pos) ++;
		
		memcpy(mapq, bin+12, 4);
		int l_read_name = (*mapq)& 0xff;
		(*mapq) = ((*mapq)>>8)&0xff;

		int seq_len;
		memcpy(&seq_len, bin + 20,4);
		memcpy(&mate_refID, bin+24, 4);
		if(mate_refID>=0){
			/*if(sambam_chro_table[mate_refID].chro_name < NULL + 0xfffff){
				SUBREADprintf("DANGEROUS: PARSE: matechro[%d] = %p, TABLE_PTR=%p\n", mate_refID , sambam_chro_table[mate_refID].chro_name, sambam_chro_table);
			}*/
			(*mate_chro) = sambam_chro_table[mate_refID].chro_name;
		}
		else	(*mate_chro) = NULL;

		(*mate_pos)=0;
		memcpy(mate_pos, bin+28, 4);
		(*mate_pos)++;

		int tlen_int;
		memcpy(&tlen_int, bin+32, 4);
		(*tlen) = tlen_int;

		int * cigar_opt_ints = (int *)(bin + 36 + l_read_name);
		unsigned int chro_cursor = (*pos), section_start_chro = (*pos);
		unsigned short read_cursor = 0, this_section_length = 0, section_start_read = 0;

		if(intervals_buffer){
			intervals_buffer[ *intervals_i ].start_pos = chro_cursor;
			intervals_buffer[ *intervals_i ].chro = *chro;
		}

		for(x1 = 0 ; x1 < cigar_opts; x1++){
			int optype = cigar_opt_ints[x1]&0xf;
			int optval = (cigar_opt_ints[x1]>>4)& 0xfffffff;
			if(optype == 0 || optype == 7 || optype == 8){ // 'M' , '=', 'X'
				chro_cursor += optval;
				read_cursor += optval;
				this_section_length += optval;
/*			}else if(optype == 1){ // 'I'
				read_cursor += optval;
			}else if(optype == 2){ // 'D'
				chro_cursor += optval;
*/			}else if(optype == 1 || optype == 2 || optype == 3){ // 'I', 'D' or 'N'
				if(3 == optype)
					(*is_junction_read) = 1;
				char event_char=0;
				if(optype == 3) event_char = 'N';
				if(optype == 2) event_char = 'D';
				else if(optype == 1){
					if(intervals_buffer && intervals_buffer[ *intervals_i ].insertions < MAXIMUM_INSERTION_IN_SECTION){
						intervals_buffer[ *intervals_i ].insertion_start_pos[  intervals_buffer[ *intervals_i ].insertions  ] = chro_cursor;
						intervals_buffer[ *intervals_i ].insertion_lengths[ intervals_buffer[ *intervals_i ].insertions ] = optval;
						intervals_buffer[ *intervals_i ].insertions ++;
					}
					event_char = 'I';
				}

				if( (*cigar_sect) < max_M){
					Event_After_Section[*cigar_sect] = event_char;
					Starting_Chro_Points[*cigar_sect] = section_start_chro; 
					Starting_Read_Points[*cigar_sect] = section_start_read;
					Section_Read_Lengths[*cigar_sect] = this_section_length;
					ChroNames[*cigar_sect] = (*chro);
					(*cigar_sect)++;

					if(intervals_buffer){
						intervals_buffer[ *intervals_i ].chromosomal_length = chro_cursor - intervals_buffer[ *intervals_i ].start_pos;
						(*intervals_i) ++;
					}
				}

				if(optype == 2 || optype == 3)// N or D
					chro_cursor += optval;
				else
					read_cursor += optval;

				if(intervals_buffer && (*cigar_sect) < max_M){
					intervals_buffer[ *intervals_i ].start_pos = chro_cursor;
					intervals_buffer[ *intervals_i ].chro = *chro;
				}

				section_start_chro = chro_cursor;
				section_start_read = read_cursor;
				this_section_length = 0;
			}else if(optype == 4){ // 'S'
				if(read_cursor==0)
				{
					read_cursor += optval;
					section_start_read = read_cursor;

					if(intervals_buffer){
						if(intervals_buffer[ *intervals_i ].start_pos > optval) intervals_buffer[ *intervals_i ].start_pos -= optval;
						else intervals_buffer[ *intervals_i ].start_pos = 0;
					}
				}else	len_of_S1 = optval;
			}	// H and P do not have effect on cigar parsing.
		}
		if(this_section_length>0){
			// add new section
			if( (*cigar_sect) < max_M){
				if(intervals_buffer){
					intervals_buffer[ *intervals_i ].chromosomal_length = chro_cursor - intervals_buffer[ *intervals_i ].start_pos + len_of_S1;
					(*intervals_i)++;
				}
				Starting_Chro_Points[*cigar_sect] = section_start_chro; 
				Starting_Read_Points[*cigar_sect] = section_start_read;
				Section_Read_Lengths[*cigar_sect] = this_section_length ;
				ChroNames[*cigar_sect] = (*chro);
				(*cigar_sect)++;
			}
		}
		
		int bin_ptr = 36 + l_read_name + seq_len + (seq_len+1)/2 + 4 * cigar_opts;
		int block_len;
		memcpy(&block_len, bin, 4);
		int found_NH = SAM_pairer_iterate_int_tags((unsigned char *)bin+bin_ptr, block_len + 4 - bin_ptr, "NH", NH_value);
		if(!found_NH) *(NH_value) = 1;
		//SUBREADprintf("FOUND=%d, NH=%d, TAG=%.*s\n", found_NH, *(NH_value), 3 , bin+bin_ptr);
	}else{
		(*read_name) = bin2 + 36;
		int mate_flag;
		memcpy(&mate_flag, bin2 + 16, 4);
		mate_flag = mate_flag >> 16;
		(*flag) = reverse_flag(mate_flag);

		int refID, mate_refID;
		memcpy(&refID, bin2 + 24, 4);
		memcpy(&mate_refID, bin2 + 4, 4);
		if(refID < 0) *chro = NULL; 
		else (*chro) = sambam_chro_table[refID].chro_name;

		if(mate_refID < 0) *mate_chro = NULL;
		else (*mate_chro) = sambam_chro_table[mate_refID].chro_name; 

		*pos=0;
		memcpy(pos, bin2+28, 4);
		(*pos)++;

		*mate_pos=0;
		memcpy(mate_pos, bin2+8, 4);
		(*mate_pos)++;
	
		(*tlen) = 0;
		memcpy(tlen, bin2+32, 4);
		(*tlen) = -(*tlen);
	}
}

/*
typedef struct {
        char chromosome_name_left[CHROMOSOME_NAME_LENGTH + 1];
        char chromosome_name_right[CHROMOSOME_NAME_LENGTH + 1];
        unsigned int last_exon_base_left;
        unsigned int first_exon_base_right;
} fc_junction_info_t;

*/
int calc_junctions_from_cigarInts(fc_thread_global_context_t * global_context, int alignment_masks , int cigar_sections, unsigned int * Starting_Chro_Points, unsigned short * Starting_Read_Points, unsigned short * Section_Lengths, char ** ChroNames, char * Event_After_Section, fc_junction_info_t * junctions_current){
	int x1, ret = 0;
	unsigned int last_base_pos = Starting_Chro_Points[0] + Section_Lengths[0] - 1;
	for(x1 = 1; x1 < cigar_sections; x1++){
		if(Event_After_Section[x1-1] == 'N'){
			unsigned int first_base_pos = Starting_Chro_Points[x1];
			junctions_current[ret].last_exon_base_left = last_base_pos;
			junctions_current[ret].first_exon_base_right = first_base_pos;
			strcpy(junctions_current[ret].chromosome_name_left, ChroNames[x1]);
			strcpy(junctions_current[ret].chromosome_name_right, ChroNames[x1]);
			ret ++;
		}

		last_base_pos = Starting_Chro_Points[x1] + Section_Lengths[x1] - 1;
	}
	return ret;
}

void add_fragment_supported_junction(	fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context, fc_junction_info_t * supported_junctions1,
					int njunc1, fc_junction_info_t * supported_junctions2, int njunc2);

void process_line_junctions(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context, char * bin1, char * bin2) {
	fc_junction_info_t supported_junctions1[global_context -> max_M], supported_junctions2[global_context -> max_M];
	int is_second_read, njunc1=0, njunc2=0, is_junction_read, cigar_sections;
	int alignment_masks, mapping_qual, NH_value;

	for(is_second_read = 0 ; is_second_read < 2; is_second_read++){
		char * read_chr, *read_name, *mate_chr;
		long read_pos, fragment_length = 0, mate_pos;
		unsigned int Starting_Chro_Points[global_context -> max_M];
		unsigned short Starting_Read_Points[global_context -> max_M];
		unsigned short Section_Read_Lengths[global_context -> max_M];
		char * ChroNames[global_context -> max_M];
		char Event_After_Section[global_context -> max_M];
		if(is_second_read && !global_context -> is_paired_end_mode_assign) break;

		parse_bin(global_context -> sambam_chro_table, is_second_read?bin2:bin1, is_second_read?bin1:bin2 , &read_name,  &alignment_masks , &read_chr, &read_pos, &mapping_qual, &mate_chr, &mate_pos, &fragment_length, &is_junction_read, &cigar_sections, Starting_Chro_Points, Starting_Read_Points, Section_Read_Lengths, ChroNames, Event_After_Section, &NH_value, global_context -> max_M, NULL, NULL);
		assert(cigar_sections <= global_context -> max_M);

		int * njunc_current = is_second_read?&njunc2:&njunc1;
		fc_junction_info_t * junctions_current = is_second_read?supported_junctions2:supported_junctions1;
		(*njunc_current) = calc_junctions_from_cigarInts(global_context, alignment_masks , cigar_sections, Starting_Chro_Points, Starting_Read_Points, Section_Read_Lengths, ChroNames, Event_After_Section, junctions_current);

		//if(0 && FIXLENstrcmp("HWI-ST212:219:C0C1TACXX:1:1101:13391:171460", read_name)==0){
		//	SUBREADprintf("JUNC_FOUND_IN_READ OF %s : %d\n", read_name , *njunc_current);
		//}
	}
	if(njunc1 >0 || njunc2>0)
		add_fragment_supported_junction(global_context, thread_context, supported_junctions1, njunc1, supported_junctions2, njunc2);

}

int process_pairer_output(void * pairer_vp, int thread_no, char * rname, char * bin1, char * bin2){
	SAM_pairer_context_t * pairer = (SAM_pairer_context_t *) pairer_vp;
	fc_thread_global_context_t * global_context = (fc_thread_global_context_t * )pairer -> appendix1;
	fc_thread_thread_context_t * thread_context = global_context -> thread_contexts + thread_no;

	//#warning "++++++ REMOVE THIS RETURN ++++++"
	//return 0;

	/*if(bin1) convert_bin_to_read( bin1, thread_context -> line_buffer1 , global_context -> sambam_chro_table);
	else    make_dummy(rname, bin2, thread_context -> line_buffer1,  global_context -> sambam_chro_table);
	if(bin2) convert_bin_to_read( bin2, thread_context -> line_buffer2 , global_context -> sambam_chro_table );
	else	make_dummy(rname, bin1, thread_context -> line_buffer2,  global_context -> sambam_chro_table);*/
	process_line_buffer(global_context, thread_context, bin1, bin2);
	if(global_context -> do_junction_counting){
		process_line_junctions(global_context, thread_context, bin1, bin2);
	}
	return 0;
}

void sort_bucket_table(fc_thread_global_context_t * global_context);
void vote_and_add_count(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context,
			long * hits_indices1, int nhits1, long * hits_indices2, int nhits2, unsigned int total_frag_len,
			char ** hits_chro1, char ** hits_chro2, unsigned int * hits_start_pos1, unsigned int * hits_start_pos2, unsigned short * hits_length1, unsigned short * hits_length2,
			int fixed_fractional_count, char * read_name);

void process_line_buffer(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context, char * bin1, char * bin2)
{

	char * read_chr, *read_name, *mate_chr;
	long read_pos, fragment_length = 0, mate_pos;
	unsigned int search_start = 0, search_end;
	int nhits1 = 0, nhits2 = 0, alignment_masks, search_block_id, search_item_id, mapping_qual;
	long * hits_indices1 = thread_context -> hits_indices1, * hits_indices2 = thread_context -> hits_indices2;
	unsigned int * hits_start_pos1 = thread_context -> hits_start_pos1 ,  * hits_start_pos2 = thread_context -> hits_start_pos2;
	unsigned short * hits_length1 = thread_context -> hits_length1 ,  * hits_length2 = thread_context -> hits_length2;
	char ** hits_chro1 = thread_context -> hits_chro1 , **hits_chro2 = thread_context -> hits_chro2;
	unsigned int  total_frag_len =0;

	int cigar_sections, is_junction_read;
	unsigned int Starting_Chro_Points[global_context -> max_M];
	unsigned short Starting_Read_Points[global_context -> max_M];
	unsigned short Section_Read_Lengths[global_context -> max_M];
	char * ChroNames[global_context -> max_M];
	char Event_After_Section[global_context -> max_M];

	int is_second_read;
	int maximum_NH_value = 1, NH_value;
	int skipped_for_exonic = 0;
	int first_read_quality_score = 0, CIGAR_intervals_R1_sections = 0, CIGAR_intervals_R2_sections = 0;

	CIGAR_interval_t CIGAR_intervals_R1 [ global_context -> max_M ];
	CIGAR_interval_t CIGAR_intervals_R2 [ global_context -> max_M ];

	if(global_context -> need_calculate_overlap_len ){
		memset( CIGAR_intervals_R1, 0, sizeof(CIGAR_interval_t) *  global_context -> max_M  );
		memset( CIGAR_intervals_R2, 0, sizeof(CIGAR_interval_t) *  global_context -> max_M  );
	}

	thread_context->all_reads++;
	//if(thread_context->all_reads>1000000) printf("TA=%llu\n%s\n",thread_context->all_reads, thread_context -> line_buffer1);


	for(is_second_read = 0 ; is_second_read < 2; is_second_read++)
	{
		if(is_second_read && !global_context -> is_paired_end_mode_assign) break;

		parse_bin(global_context -> sambam_chro_table, is_second_read?bin2:bin1, is_second_read?bin1:bin2 , &read_name,  &alignment_masks , &read_chr, &read_pos, &mapping_qual, &mate_chr, &mate_pos, &fragment_length, &is_junction_read, &cigar_sections, Starting_Chro_Points, Starting_Read_Points, Section_Read_Lengths, ChroNames, Event_After_Section, &NH_value, global_context -> max_M , global_context -> need_calculate_overlap_len?(is_second_read?CIGAR_intervals_R2:CIGAR_intervals_R1):NULL, is_second_read?&CIGAR_intervals_R2_sections:&CIGAR_intervals_R1_sections );
	//	SUBREADprintf("  RNAME=%s\n", read_name);

		//#warning "==================== REMOVE WHEN RELEASE ========================"
		//if(global_context -> SAM_output_fp)
		//	fprintf(global_context -> SAM_output_fp, "SAMDEBUG: %s\t\t%s, %ld\n", read_name, read_chr, read_pos);
		if(is_second_read == 0)
		{
			//skip the read if unmapped (its mate will be skipped as well if paired-end)
			if( ((!global_context -> is_paired_end_mode_assign) &&  (alignment_masks & SAM_FLAG_UNMAPPED) ) ||
			    ((alignment_masks & SAM_FLAG_UNMAPPED)   &&  (alignment_masks & SAM_FLAG_MATE_UNMATCHED) && global_context -> is_paired_end_mode_assign) ||
			    (((alignment_masks & SAM_FLAG_UNMAPPED) || (alignment_masks & SAM_FLAG_MATE_UNMATCHED)) && global_context -> is_paired_end_mode_assign && global_context -> is_both_end_required)
			  ){
				thread_context->read_counters.unassigned_unmapped ++;

				if(global_context -> SAM_output_fp)
					fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_Unmapped\t*\t*\n", read_name);

				return;	// do nothing if a read is unmapped, or the first read in a pair of reads is unmapped.
			}
		}

		if(global_context -> is_paired_end_mode_assign && (!global_context ->is_SEPEmix_warning_shown)){
			if(((!global_context -> is_paired_end_input_file)  && ( alignment_masks & SAM_FLAG_PAIRED_TASK )) || ((global_context -> is_paired_end_input_file)  && 0 == ( alignment_masks & SAM_FLAG_PAIRED_TASK ))){
				print_in_box(85,0,0,"   %c[31mBoth single-end and paired-end reads were found.", 27);
				//SUBREADprintf("BAD READ:%s, FLAG=%d\n", read_name, alignment_masks);
				global_context ->is_SEPEmix_warning_shown = 1;
			}
		}

		if(global_context -> min_mapping_quality_score>0)
		{
			//printf("SECOND=%d; FIRST=%d; THIS=%d; Q=%d\n", is_second_read, first_read_quality_score, mapping_qual, );
			if(( mapping_qual < global_context -> min_mapping_quality_score  && ! global_context -> is_paired_end_mode_assign)||( is_second_read  && max( first_read_quality_score, mapping_qual ) < global_context -> min_mapping_quality_score))
			{
				thread_context->read_counters.unassigned_mappingquality ++;

				if(global_context -> SAM_output_fp)
				{
					fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_MappingQuality\t*\tMapping_Quality=%d,%d\n", read_name, first_read_quality_score, mapping_qual);
				}
				return;
			}
			if(is_second_read==0 && global_context -> is_paired_end_mode_assign)
			{
				first_read_quality_score = mapping_qual;
			}
		}

		if(is_second_read == 0 && global_context -> is_paired_end_mode_assign && 
	   	  (global_context -> is_PE_distance_checked || global_context -> is_chimertc_disallowed)
		  )
		{
			int is_half_mapped = (alignment_masks & SAM_FLAG_UNMAPPED) || (alignment_masks & SAM_FLAG_MATE_UNMATCHED);

			if(!is_half_mapped)
			{
				fragment_length = abs( fragment_length ); //get the fragment length

				int is_first_read_negative_strand = (alignment_masks & SAM_FLAG_REVERSE_STRAND_MATCHED)?1:0; 
				int is_second_read_negative_strand = (alignment_masks & SAM_FLAG_MATE_REVERSE_STRAND_MATCHED)?1:0; 

				if(mate_chr == read_chr && is_first_read_negative_strand!=is_second_read_negative_strand) {
				 //^^^^^^^^^^^^^^^^^^^^ They are directly compared because they are both pointers in the same contig name table.
				 //
					if(global_context -> is_PE_distance_checked && ((fragment_length > global_context -> max_paired_end_distance) || (fragment_length < global_context -> min_paired_end_distance))) {
						thread_context->read_counters.unassigned_fragmentlength ++;

						if(global_context -> SAM_output_fp)
							fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_FragmentLength\t*\tLength=%ld\n", read_name, fragment_length);
						return;
					}
				} else {
					if(global_context -> is_chimertc_disallowed) {
						thread_context->read_counters.unassigned_chimericreads ++;

						if(global_context -> SAM_output_fp)
							fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_Chimera\t*\t*\n", read_name);
						return;
					}
				}
			}
		}

		// This filter has to be put here because the 0x400 FLAG is not about mapping but about sequencing.
		// A unmapped read with 0x400 FLAG should be able to kill the mapped mate which may have no 0x400 FLAG. 
		if(global_context -> is_duplicate_ignored)
		{
			if(alignment_masks & SAM_FLAG_DUPLICATE)
			{
				thread_context->read_counters.unassigned_duplicate ++;
				if(global_context -> SAM_output_fp)
					fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_Duplicate\t*\t*\n", read_name);

				return;
			}

		}

		if(SAM_FLAG_UNMAPPED & alignment_masks) continue;

		if( NH_value > 1 ) {
			if(global_context -> is_multi_mapping_allowed == 0)
			{
				// now it is a NH>1 read!
				// not allow multimapping -> discard!
				thread_context->read_counters.unassigned_multimapping ++;

				if(global_context -> SAM_output_fp)
					fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_MultiMapping\t*\t*\n", read_name);

				return;
			}
		}

		maximum_NH_value = max(maximum_NH_value, NH_value);

		// if a pair of reads have one secondary, the entire fragment is seen as secondary.
		if((alignment_masks & SAM_FLAG_SECONDARY_MAPPING) && (global_context -> is_multi_mapping_allowed == ALLOW_PRIMARY_MAPPING))
		{
			thread_context->read_counters.unassigned_secondary ++;

			if(global_context -> SAM_output_fp)
				fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_Secondary\t*\t*\n", read_name);
			return;
		}

		int is_this_negative_strand = (alignment_masks & SAM_FLAG_REVERSE_STRAND_MATCHED)?1:0; 
		int is_fragment_negative_strand = is_this_negative_strand;

		if( global_context -> is_paired_end_mode_assign ){
			int is_second_read_in_pair = alignment_masks & SAM_FLAG_SECOND_READ_IN_PAIR;
			//is_fragment_negative_strand = is_second_read_in_pair?(!is_this_negative_strand):is_this_negative_strand;
			if(is_second_read_in_pair)
				is_fragment_negative_strand = global_context -> is_second_read_straight?is_this_negative_strand:(!is_this_negative_strand);
			else
				is_fragment_negative_strand = global_context -> is_first_read_reversed?(!is_this_negative_strand):is_this_negative_strand;
		}

		int nhits = 0;

		int cigar_section_id;
		long * hits_indices = (is_second_read?hits_indices2:hits_indices1);
		unsigned int * hits_start_pos = is_second_read?hits_start_pos2:hits_start_pos1;
		unsigned short * hits_length = is_second_read?hits_length2:hits_length1;
		char ** hits_chro = is_second_read?hits_chro2:hits_chro1;

		if(global_context->is_split_or_exonic_only == 1 && !is_junction_read) {
			skipped_for_exonic ++;

			if(skipped_for_exonic == 1 + global_context -> is_paired_end_mode_assign){
				if(global_context -> SAM_output_fp)
					fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_%s\t*\t*\n", read_name, (global_context->is_split_or_exonic_only == 2)?"Hasjunction":"Nonjunction");

				thread_context->read_counters.unassigned_junction_condition ++;
				return;
			}
		}


		if(global_context->is_split_or_exonic_only == 2 && is_junction_read) {
			if(global_context -> SAM_output_fp)
				fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_%s\t*\t*\n", read_name, (global_context->is_split_or_exonic_only == 2)?"Hasjunction":"Nonjunction");
			thread_context->read_counters.unassigned_junction_condition ++;
			return;
		}

		if(1) {
		//#warning "=================== COMMENT THESE 2 LINES ================================"
		//for(cigar_section_id = 0; cigar_section_id<cigar_sections; cigar_section_id++)
		//	SUBREADprintf("BCCC: %llu , sec[%d] %s: %u ~ %u ; secs=%d ; flags=%d ; second=%d\n", read_pos, cigar_section_id , ChroNames[cigar_section_id] , Starting_Chro_Points[cigar_section_id], Section_Lengths[cigar_section_id], cigar_sections, alignment_masks, is_second_read);

			if(global_context -> reduce_5_3_ends_to_one)
			{
				if((REDUCE_TO_5_PRIME_END == global_context -> reduce_5_3_ends_to_one) + is_this_negative_strand == 1) // reduce to 5' end (small coordinate if positive strand / large coordinate if negative strand)
				{
					Section_Read_Lengths[0]=1;
				}
				else
				{
					Starting_Chro_Points[0] = Starting_Chro_Points[cigar_sections-1] + Section_Read_Lengths[cigar_sections-1] - 1;
					Section_Read_Lengths[0]=1;
				}

				cigar_sections = 1;
			}

			// Extending the reads to the 3' and 5' ends. (from the read point of view) 
			if(global_context -> five_end_extension)
			{
				if(is_this_negative_strand){
					Section_Read_Lengths [cigar_sections - 1] += global_context -> five_end_extension;
				}else{
					//SUBREADprintf("5-end extension: %d [%d]\n", Starting_Chro_Points[0], Section_Lengths[0]);
					if( read_pos > global_context -> five_end_extension)
					{
						Section_Read_Lengths [0] += global_context -> five_end_extension;
						Starting_Chro_Points [0] -= global_context -> five_end_extension;
					}
					else
					{
						Section_Read_Lengths [0] += read_pos-1;
						Starting_Chro_Points [0] -= read_pos-1;
					}
				}
			}

			if(global_context -> three_end_extension)
			{

				if(is_this_negative_strand){
					if( read_pos > global_context -> three_end_extension)
					{
						Section_Read_Lengths [0] += global_context -> three_end_extension;
						Starting_Chro_Points [0] -= global_context -> three_end_extension;
					}
					else
					{
						Section_Read_Lengths [0] += read_pos - 1;
						Starting_Chro_Points [0] -= read_pos - 1;
					}
				}
				else	Section_Read_Lengths [cigar_sections - 1] += global_context -> three_end_extension;

			}

			for(cigar_section_id = 0; cigar_section_id<cigar_sections; cigar_section_id++)
			{
				long section_begin_pos = Starting_Chro_Points[cigar_section_id];
				long section_end_pos = Section_Read_Lengths[cigar_section_id] + section_begin_pos - 1;

				
				int start_reverse_table_index = section_begin_pos / REVERSE_TABLE_BUCKET_LENGTH;
				int end_reverse_table_index = (1+section_end_pos) / REVERSE_TABLE_BUCKET_LENGTH;

				/*if(ChroNames[cigar_section_id] < (char *)NULL + 0xfffff)
					SUBREADprintf("DANGEROUS! RNAME=%s,  CNAME=%p,  LEN_P=%d,  SECID=%d\n", read_name, ChroNames[cigar_section_id], Section_Read_Lengths[cigar_section_id], cigar_section_id);*/

				fc_chromosome_index_info * this_chro_info = HashTableGet(global_context -> exontable_chro_table, ChroNames[cigar_section_id]);
				if(this_chro_info == NULL)
				{
					if(global_context -> annot_chro_name_alias_table)
					{
						char * anno_chro_name = HashTableGet( global_context -> annot_chro_name_alias_table , ChroNames[cigar_section_id]);
						if(anno_chro_name)
							this_chro_info = HashTableGet(global_context -> exontable_chro_table, anno_chro_name);
					}
					if(this_chro_info == NULL && memcmp(ChroNames[cigar_section_id], "chr", 3)==0)
					{
						this_chro_info = HashTableGet(global_context -> exontable_chro_table, ChroNames[cigar_section_id]+3);
					//	SUBREADprintf("INQ: %p : '%s'\n", this_chro_info , ChroNames[cigar_section_id]+3);
					}
					if(this_chro_info == NULL && strlen(ChroNames[cigar_section_id])<=2)
					{
						strcpy(thread_context -> chro_name_buff, "chr");
						strcpy(thread_context -> chro_name_buff+3, ChroNames[cigar_section_id]);
						this_chro_info = HashTableGet(global_context -> exontable_chro_table, thread_context -> chro_name_buff);
					}
				}

				//SUBREADprintf("INF: %p : %s\n", this_chro_info , ChroNames[cigar_section_id]);

				if(this_chro_info)
				{
					start_reverse_table_index = min(start_reverse_table_index, this_chro_info-> chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH);
					end_reverse_table_index = min(end_reverse_table_index, this_chro_info-> chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH+ 1);

					while(start_reverse_table_index<=end_reverse_table_index)
					{
						search_start = this_chro_info -> reverse_table_start_index [start_reverse_table_index];
						if(search_start<0xffffff00)break;
						start_reverse_table_index++;
					}
					if(search_start>0xffffff00) continue;

					//search_start = this_chro_info -> chro_block_table_start;

					search_end = this_chro_info -> chro_block_table_end;//reverse_table_end_index [end_reverse_table_index];
		
					for(search_block_id=search_start;search_block_id<search_end;search_block_id++){
						if (global_context -> exontable_block_min_start[search_block_id] > section_end_pos) break;
						if (global_context -> exontable_block_max_end[search_block_id] < section_begin_pos) continue;

						int search_item_start = 0, search_item_end = global_context -> exontable_block_end_index[search_block_id];
						if(search_block_id>0)search_item_start = global_context -> exontable_block_end_index[search_block_id-1];

						// search_item_id is the inner number of the exons.
						// the exontables in global_index has search_item_id as the index.

						for(search_item_id = search_item_start ; search_item_id < search_item_end; search_item_id++)
						{
							if (global_context -> exontable_stop[search_item_id] >= section_begin_pos)
							{
								if (global_context -> exontable_start[search_item_id] > section_end_pos) break;
								// there is an overlap >=1 between read and feature.
								// the overlap length is min(end_r, end_F) - max(start_r, start_F) + 1
								
								int is_strand_ok =1;

								if(global_context->is_strand_checked){
									if(global_context->is_strand_checked == 1)
										is_strand_ok = (is_fragment_negative_strand == global_context -> exontable_strand[search_item_id]);
									else// if(global_context->is_strand_checked == 2)
										is_strand_ok = (is_fragment_negative_strand != global_context -> exontable_strand[search_item_id]);
									//SUBREADprintf("%d = %d == %d\n", is_strand_ok, is_fragment_negative_strand, global_context -> exontable_strand[search_item_id]);
								}

								if(is_strand_ok){
									if(nhits<=MAX_HIT_NUMBER - 1)
									{
										hits_indices[nhits] = search_item_id;

										if(global_context -> need_calculate_overlap_len) {
											hits_start_pos[nhits] = max(Starting_Chro_Points[cigar_section_id], global_context -> exontable_start[search_item_id]);
											hits_length[nhits] =  min(global_context -> exontable_stop[search_item_id] , section_end_pos)+1 - hits_start_pos[nhits] ;
											hits_chro[nhits] = ChroNames[cigar_section_id];
											if(0 && FIXLENstrcmp("V0112_0155:7:1101:10214:3701", read_name)==0)
												SUBREADprintf("QNAME: [%d] %s %d ~ %d\n", nhits, hits_chro[nhits],  hits_start_pos[nhits],  hits_start_pos[nhits]+hits_length[nhits]);
										}

										nhits++;
									}
									else break;
								}
							} 
						}
					}
				}
			}
		}


		if(is_second_read) nhits2 = nhits;
		else	nhits1 = nhits;
	}	// loop for is_second_read


	if(global_context -> need_calculate_fragment_len )
		total_frag_len = calc_total_frag_len( global_context, thread_context, CIGAR_intervals_R1, CIGAR_intervals_R1_sections, CIGAR_intervals_R2, CIGAR_intervals_R2_sections , read_name);

	//SUBREADprintf("FRAGLEN: %s %d\n", read_name, total_frag_len);

	int fixed_fractional_count = global_context -> use_fraction_multi_mapping ?calc_fixed_fraction(maximum_NH_value): NH_FRACTION_INT;

	// we have hits_indices1 and hits_indices2 and nhits1 and nhits2 here
	// we also have fixed_fractional_count which is the value to add

	vote_and_add_count(global_context, thread_context,
			   hits_indices1,  nhits1, hits_indices2,  nhits2, total_frag_len,
			   hits_chro1, hits_chro2, hits_start_pos1, hits_start_pos2, hits_length1, hits_length2,
			   fixed_fractional_count, read_name);
	return;
}

void add_bitmap_overlapping(char * x1_bitmap, short start_base, short len){
	int x1;
	int rl16 = start_base+len-16;
	for(x1 = start_base; x1 < start_base+len; x1++){
		int bit = x1 % 8;
		int byte = x1 / 8;
		if(bit == 0 && x1 < rl16){
			x1_bitmap[byte]=-1;
			x1_bitmap[byte+1]=-1;
			x1+=15;
		}else{
			x1_bitmap[byte] |= (1<<bit);
		}
	}
}

int count_bitmap_overlapping(char * x1_bitmap, unsigned short rl){

	int x1;
	int ret = 0;
	for(x1 = 0; x1 < rl; x1++){
		int byte = x1 / 8;
		int bit = x1 % 8;

		if(bit == 0 && x1_bitmap[byte]==-1){
			x1 += 7;
			ret += 8;
		}else if(x1_bitmap[byte] &  (1<<bit)) ret ++;
	}
	return ret;
}

void add_fragment_supported_junction(	fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context, fc_junction_info_t * supported_junctions1,
					int njunc1, fc_junction_info_t * supported_junctions2, int njunc2){
	assert(njunc1 >= 0 && njunc1 <= global_context -> max_M -1 );
	assert(njunc2 >= 0 && njunc2 <= global_context -> max_M -1 );
	int x1,x2, in_total_junctions = njunc2 + njunc1;
	for(x1 = 0; x1 < in_total_junctions; x1 ++){
		fc_junction_info_t * j_one = (x1 >= njunc1)?supported_junctions2+(x1-njunc1):(supported_junctions1+x1);
		if(j_one->chromosome_name_left[0]==0) continue;

		for(x2 = x1+1; x2 < in_total_junctions ; x2 ++){
			fc_junction_info_t * j_two = (x2 >= njunc1)?supported_junctions2+(x2-njunc1):(supported_junctions1+x2);
			if(j_two->chromosome_name_left[0]==0) continue;
			if(
				j_one -> last_exon_base_left == j_two -> last_exon_base_left &&
				j_one -> first_exon_base_right == j_two -> first_exon_base_right &&
				strcmp(j_one -> chromosome_name_left, j_two -> chromosome_name_left) == 0 &&
				strcmp(j_one -> chromosome_name_right, j_two -> chromosome_name_right) == 0
			) j_two -> chromosome_name_left[0]=0;
		}

		char * this_key = malloc(strlen(j_one->chromosome_name_left) + strlen(j_one->chromosome_name_right)  + 36);
		sprintf(this_key, "%s\t%u\t%s\t%u", j_one->chromosome_name_left, j_one -> last_exon_base_left, j_one->chromosome_name_right, j_one -> first_exon_base_right);
		void * count_ptr = HashTableGet(thread_context -> junction_counting_table, this_key);
		unsigned long long count_junc = count_ptr - NULL;
		HashTablePut(thread_context -> junction_counting_table, this_key, NULL+count_junc + 1);

//		#warning "CONTINUE SHOULD BE REMOVED!!!!"
//			continue;

		char * left_key = malloc(strlen(j_one->chromosome_name_left) + 16);
		char * right_key = malloc(strlen(j_one->chromosome_name_right) + 16);
		sprintf(left_key, "%s\t%u", j_one->chromosome_name_left, j_one -> last_exon_base_left);
		sprintf(right_key, "%s\t%u", j_one->chromosome_name_right, j_one -> first_exon_base_right);

		for( x2 = 0 ; x2 < 2 ; x2++ ){
			char * lr_key = x2?right_key:left_key;
			count_ptr = HashTableGet(thread_context -> splicing_point_table, lr_key);
			count_junc = count_ptr - NULL;
			HashTablePut(thread_context -> splicing_point_table, lr_key, NULL + count_junc + 1);
		}
	}
}

int overlap_compare(void * arr, int L, int R){
	unsigned int * pos = (unsigned int *)arr;
	return pos[ L*2 ] -  pos[R*2];
}

void overlap_exchange(void * arr, int L, int R){
	unsigned int * pos = (unsigned int *)arr, tt;
	tt=pos[L*2];
	pos[L*2] = pos[R*2];
	pos[R*2] = tt;

	tt=pos[L*2+1];
	pos[L*2+1] = pos[R*2+1];
	pos[R*2+1] = tt;
}

unsigned short calc_score_overlaps(fc_thread_global_context_t * global_context,  fc_thread_thread_context_t * thread_context, char ** chros, unsigned int * start_poses, unsigned short * lens, int sections){
	unsigned int in_intervals[ 2*sections ];
	unsigned int out_intervals[ 2*sections ], x1;
	char used_interval[ sections ];

	memset(used_interval, 0 , sections);
	unsigned int ret = 0;

	for(x1 = 0  ; x1 < sections ; x1++){
		if( used_interval [x1] )continue;

		in_intervals[x1*2] = start_poses[x1];
		in_intervals[x1*2 + 1] = start_poses[x1] + lens[x1];
		used_interval[x1]=1;
	
		int x2, this_sections = 1;
		for(x2 = x1 + 1; x2 < sections; x2++){
			if(strcmp( chros[x2], chros[x1] ) == 0){
				in_intervals[this_sections*2] = start_poses[x2];
				in_intervals[this_sections*2 + 1] = start_poses[x2] + lens[x2];
				used_interval[x2]=1;
				this_sections++;
			}
		}

		basic_sort( in_intervals, this_sections, overlap_compare, overlap_exchange );

		int merged_secs = mergeIntervals( in_intervals, out_intervals, this_sections );
		for(x2 = 0; x2 < merged_secs; x2++){
			ret += ( out_intervals[x2*2+1] - out_intervals[x2*2] );
		}
	}
	return ret;
}


void vote_and_add_count(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context,
			long * hits_indices1, int nhits1, long * hits_indices2, int nhits2, unsigned int total_frag_len,
			char ** hits_chro1, char ** hits_chro2, unsigned int * hits_start_pos1, unsigned int * hits_start_pos2, unsigned short * hits_length1, unsigned short * hits_length2,
			int fixed_fractional_count, char * read_name){
	if(global_context -> need_calculate_overlap_len == 0 && nhits2+nhits1==1) {
		long hit_exon_id = nhits2?hits_indices2[0]:hits_indices1[0];
		thread_context->count_table[hit_exon_id] += fixed_fractional_count;
		thread_context->nreads_mapped_to_exon++;
		if(global_context -> SAM_output_fp)
		{
			int final_gene_number = global_context -> exontable_geneid[hit_exon_id];
			unsigned char * final_feture_name = global_context -> gene_name_array[final_gene_number];
			fprintf(global_context -> SAM_output_fp,"%s\tAssigned\t%s\tTotal=1\n", read_name, final_feture_name);
		}
		thread_context->read_counters.assigned_reads ++;
	} else if(global_context -> need_calculate_overlap_len == 0 && nhits2 == 1 && nhits1 == 1 && hits_indices2[0]==hits_indices1[0]) {
		long hit_exon_id = hits_indices1[0];
		thread_context->count_table[hit_exon_id] += fixed_fractional_count;
		thread_context->nreads_mapped_to_exon++;
		if(global_context -> SAM_output_fp)
		{
			int final_gene_number = global_context -> exontable_geneid[hit_exon_id];
			unsigned char * final_feture_name = global_context -> gene_name_array[final_gene_number];
			fprintf(global_context -> SAM_output_fp,"%s\tAssigned\t%s\tTotal=1\n", read_name, final_feture_name);
		}
		thread_context->read_counters.assigned_reads ++;
	} else {
		// Build a voting table.
		// The voting table should be:
		//      total_length [nhit_final] = total_length_overlapping
		//      final_id [nhit_final] = final_exon_id

		// if is_gene_leven, then decision_table_exon_ids[nhit_final] is the exon id where the count is added.

		// After all, the count is added to all hits where total_length has the maximum value.
		// If there are more than one locations having the same total_length, then the fragment is ambiguous. 
		// Count is added when "-O" is specified.

		// merge feature : if a read overlaps with an EXON twice or more times (by >=2 segments in cigar),
		//                 then the total length of the overlapped bases is calculated.
		// 
		// two ends in a fragment is considered individually; the overlapping bases are not added up.
		//

		
		unsigned int * scoring_numbers = thread_context -> scoring_buff_numbers;	// size is : MAX_HIT_NUMBER *2
		unsigned int * scoring_flags = thread_context -> scoring_buff_flags;		// size is : MAX_HIT_NUMBER *2
		unsigned short * scoring_overlappings = thread_context -> scoring_buff_overlappings;		// size is : MAX_HIT_NUMBER *2
		long * scoring_exon_ids = thread_context -> scoring_buff_exon_ids;		// size is : MAX_HIT_NUMBER *2
		int scoring_count = 0,  score_x1;


		if( global_context -> need_calculate_overlap_len ){

			char ** scoring_gap_chros = thread_context -> scoring_buff_gap_chros;
			unsigned int * scoring_gap_starts = thread_context -> scoring_buff_gap_starts; // size is : MAX_HIT_NUMBER *2;
			unsigned short * scoring_gap_lengths = thread_context -> scoring_buff_gap_lengths; 	// size is : MAX_HIT_NUMBER *2*  global_context -> max_M*2

			int end1, end2, hit_x1, hit_x2;
			char used_hit1 [nhits1];
			char used_hit2 [nhits2];
			memset(used_hit1 , 0 , nhits1);
			memset(used_hit2 , 0 , nhits2);

			for(end1 = 0; end1 < global_context -> is_paired_end_mode_assign + 1 ; end1++){
				long * hits_indices_X1 = end1?hits_indices2:hits_indices1;
				char * used_hit_X1 = end1?used_hit2:used_hit1;
				int nhit_X1 = end1?nhits2:nhits1;

				for( hit_x1 = 0 ; hit_x1 < nhit_X1; hit_x1 ++ ){
					if(used_hit_X1[hit_x1])continue;

					int gaps = 0;
					long tmp_exon_id = hits_indices_X1[hit_x1];
					long score_merge_key;
					if (global_context -> is_gene_level )
						score_merge_key = global_context -> exontable_geneid[tmp_exon_id];
					else	score_merge_key = tmp_exon_id;


					scoring_gap_chros[0 ] = (end1?hits_chro2:hits_chro1)[hit_x1]; 
					scoring_gap_starts[0 ] = (end1?hits_start_pos2:hits_start_pos1)[hit_x1]; 
					scoring_gap_lengths[0 ] = (end1?hits_length2:hits_length1)[hit_x1]; 

					gaps=1;

					scoring_flags[scoring_count] = end1?2:1;
					scoring_numbers[scoring_count] =1;
					scoring_exon_ids[scoring_count] = tmp_exon_id;

					used_hit_X1[ hit_x1 ]=1;

					for(end2 = 0; end2 < global_context -> is_paired_end_mode_assign + 1 ; end2++){
						long * hits_indices_X2 = end2?hits_indices2:hits_indices1;
						char * used_hit_X2 = end2?used_hit2:used_hit1;
						int nhit_X2 = end2?nhits2:nhits1;

						for( hit_x2 = 0 ; hit_x2 < nhit_X2; hit_x2 ++ ){
							if(used_hit_X2[hit_x2])continue;

							long X2_merge_key;
							if (global_context -> is_gene_level )
								X2_merge_key = global_context -> exontable_geneid[ hits_indices_X2[hit_x2] ];
							else	X2_merge_key = hits_indices_X2[hit_x2];

							if( X2_merge_key == score_merge_key ){
								used_hit_X2[ hit_x2 ]=1;
								scoring_gap_chros[ gaps ] = (end2?hits_chro2:hits_chro1)[hit_x2]; 
								scoring_gap_starts[ gaps ] = (end2?hits_start_pos2:hits_start_pos1)[hit_x2]; 
								scoring_gap_lengths[ gaps ] = (end2?hits_length2:hits_length1)[hit_x2]; 

								if((scoring_flags[scoring_count] & (end2?2:1))== 0 ){
									scoring_flags[scoring_count] |= end2?2:1;
									scoring_numbers[scoring_count] ++;
								}
								gaps ++;
							}
						}
					}

					scoring_overlappings [scoring_count] = calc_score_overlaps(global_context, thread_context, scoring_gap_chros, scoring_gap_starts, scoring_gap_lengths, gaps);
					if( global_context -> use_overlapping_break_tie )
						scoring_numbers[scoring_count] = scoring_overlappings [scoring_count];
					scoring_count++;
				}
			}
		}else{
			int ends;
			for(ends =0 ; ends < global_context -> is_paired_end_mode_assign + 1 ; ends++){
				int nhits = ends?nhits2:nhits1;
				long * hits_indices = ends?hits_indices2:hits_indices1;

				int hit_x1;
				for(hit_x1 = 0; hit_x1 < nhits; hit_x1++){
					long tmp_exon_id = hits_indices[hit_x1], score_merge_key;
					int found = 0;
					if (global_context -> is_gene_level )
						score_merge_key = global_context -> exontable_geneid[tmp_exon_id];
					else	score_merge_key = tmp_exon_id;

					for(score_x1 = 0; score_x1 < scoring_count; score_x1 ++){
						long score_x1_key ; 
						if (global_context -> is_gene_level )
							score_x1_key = global_context -> exontable_geneid[ scoring_exon_ids[score_x1] ];
						else	score_x1_key = scoring_exon_ids[score_x1] ;

						//fprintf(stderr, "Q222KEY: exon=%ld, gene=%ld\n", scoring_exon_ids[score_x1] , score_x1_key  );
						if( score_x1_key == score_merge_key ){
							if((scoring_flags[score_x1] & ( ends?2:1 )) == 0) {
								scoring_flags[score_x1] |= (ends?2:1);
								scoring_numbers[score_x1] ++;
							}

							found = 1;
							break;
						}
					}

					if(0 == found){
						scoring_exon_ids[scoring_count] = tmp_exon_id;
						scoring_flags[scoring_count] = ends?2:1;
						scoring_numbers[scoring_count] = 1;

						scoring_count++;
					}
				}
			}
		}


		int maximum_score = 0;
		int maximum_total_count = 0;
		int maximum_score_x1 = 0;
		int applied_fragment_minimum_overlapping = 1;

		if( global_context -> fragment_minimum_overlapping > 1 ||  global_context -> need_calculate_fragment_len){
			applied_fragment_minimum_overlapping = max( global_context -> fragment_minimum_overlapping, global_context -> fractional_minimum_overlapping * ( total_frag_len) );
		}


		for(score_x1 = 0; score_x1 < scoring_count ; score_x1++){

			if(0 && FIXLENstrcmp("V0112_0155:7:1101:5387:6362", read_name)==0) SUBREADprintf("Scoring Overlap %s = %d >=%d, score=%d, exonid=%ld\n", read_name, scoring_overlappings[score_x1], applied_fragment_minimum_overlapping, scoring_numbers[score_x1], scoring_exon_ids[score_x1]);
			//SUBREADprintf("RLTEST: %s %d\n", read_name, scoring_overlappings[score_x1]);
			if( applied_fragment_minimum_overlapping > 1 )
				if( applied_fragment_minimum_overlapping > scoring_overlappings[score_x1] ){
					scoring_numbers[score_x1] = 0;
					continue;
				}
			
			if( maximum_score < scoring_numbers[score_x1] ){
				maximum_total_count = 1;
				maximum_score = scoring_numbers[score_x1];
				maximum_score_x1 = score_x1;
			}else if( maximum_score == scoring_numbers[score_x1] )
				maximum_total_count++;
		}

		if(maximum_total_count == 0){
			if(global_context -> SAM_output_fp)
				fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_NoFeatures\t*\t*\n", read_name);

			thread_context->read_counters.unassigned_nofeatures ++;
		}else{

			// final adding votes.
			if(1 == maximum_total_count && !global_context -> is_multi_overlap_allowed) {
				// simple add to the exon ( EXON_ID = decision_table_exon_ids[maximum_decision_no])
				long max_exon_id = scoring_exon_ids[maximum_score_x1];
				thread_context->count_table[max_exon_id] += fixed_fractional_count;
				thread_context->nreads_mapped_to_exon++;
				if(global_context -> SAM_output_fp)
				{
					int final_gene_number = global_context -> exontable_geneid[max_exon_id];
					unsigned char * final_feture_name = global_context -> gene_name_array[final_gene_number];
					if(scoring_count>1)
						fprintf(global_context -> SAM_output_fp,"%s\tAssigned\t%s\tTotal=1;%s/Targets=%d/%d\n", read_name, final_feture_name, global_context -> use_overlapping_break_tie? "MaximumOverlapping":"Votes", maximum_score, scoring_count);
					else
						fprintf(global_context -> SAM_output_fp,"%s\tAssigned\t%s\tTotal=1\n", read_name, final_feture_name);
				}
				thread_context->read_counters.assigned_reads ++;
			}else if(global_context -> is_multi_overlap_allowed) {
				char final_feture_names[1000];
				int assigned_no = 0, xk1;
				final_feture_names[0]=0;
				for(xk1 = 0; xk1 < scoring_count; xk1++)
				{

					// This change was made on 31/MAR/2016
					if( scoring_numbers[xk1] < maximum_score ) continue ; 

					long tmp_voter_id = scoring_exon_ids[xk1];
					thread_context->count_table[tmp_voter_id] += calculate_multi_overlap_fraction(global_context, fixed_fractional_count, maximum_total_count);

					if(global_context -> SAM_output_fp)
					{
						if(strlen(final_feture_names)<700)
						{
							int final_gene_number = global_context -> exontable_geneid[tmp_voter_id];
							unsigned char * final_feture_name = global_context -> gene_name_array[final_gene_number];
							strncat(final_feture_names, (char *)final_feture_name, 999);
							strncat(final_feture_names, ",", 999);
							assigned_no++;
						}
					}
				}
				final_feture_names[999]=0;
				thread_context->nreads_mapped_to_exon++;
				if(global_context -> SAM_output_fp)
				{
					int ffnn = strlen(final_feture_names);
					if(ffnn>0) final_feture_names[ffnn-1]=0;
					// overlapped but still assigned 
					fprintf(global_context -> SAM_output_fp,"%s\tAssigned\t%s\tTotal=%d\n", read_name, final_feture_names, assigned_no);
				}
				thread_context->read_counters.assigned_reads ++;
			} else {
				if(global_context -> SAM_output_fp)
					fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_Ambiguity\t*\tNumber_Of_Overlapped_Genes=%d\n", read_name, maximum_total_count);

				thread_context->read_counters.unassigned_ambiguous ++;
			}
		}
	}
}

void fc_thread_merge_results(fc_thread_global_context_t * global_context, read_count_type_t * nreads , unsigned long long int *nreads_mapped_to_exon, fc_read_counters * my_read_counter, HashTable * junction_global_table, HashTable * splicing_global_table)
{
	int xk1, xk2;

	long long int total_input_reads = 0 ;
	read_count_type_t unpaired_fragment_no = 0;

	(*nreads_mapped_to_exon)=0;

	SAM_pairer_destroy(&global_context -> read_pairer);

	for(xk1=0; xk1<global_context-> thread_number; xk1++)
	{
		for(xk2=0; xk2<global_context -> exontable_exons; xk2++)
		{
			nreads[xk2]+=global_context -> thread_contexts[xk1].count_table[xk2];
		}
		total_input_reads += global_context -> thread_contexts[xk1].all_reads;
		(*nreads_mapped_to_exon) += global_context -> thread_contexts[xk1].nreads_mapped_to_exon;
		unpaired_fragment_no += global_context -> thread_contexts[xk1].unpaired_fragment_no;

		global_context -> read_counters.unassigned_ambiguous += global_context -> thread_contexts[xk1].read_counters.unassigned_ambiguous;
		global_context -> read_counters.unassigned_nofeatures += global_context -> thread_contexts[xk1].read_counters.unassigned_nofeatures;
		global_context -> read_counters.unassigned_unmapped += global_context -> thread_contexts[xk1].read_counters.unassigned_unmapped;
		global_context -> read_counters.unassigned_mappingquality += global_context -> thread_contexts[xk1].read_counters.unassigned_mappingquality;
		global_context -> read_counters.unassigned_fragmentlength += global_context -> thread_contexts[xk1].read_counters.unassigned_fragmentlength;
		global_context -> read_counters.unassigned_chimericreads += global_context -> thread_contexts[xk1].read_counters.unassigned_chimericreads;
		global_context -> read_counters.unassigned_multimapping += global_context -> thread_contexts[xk1].read_counters.unassigned_multimapping;
		global_context -> read_counters.unassigned_secondary += global_context -> thread_contexts[xk1].read_counters.unassigned_secondary;
		global_context -> read_counters.unassigned_junction_condition += global_context -> thread_contexts[xk1].read_counters.unassigned_junction_condition;
		global_context -> read_counters.unassigned_duplicate += global_context -> thread_contexts[xk1].read_counters.unassigned_duplicate;
		global_context -> read_counters.assigned_reads += global_context -> thread_contexts[xk1].read_counters.assigned_reads;

		my_read_counter->unassigned_ambiguous += global_context -> thread_contexts[xk1].read_counters.unassigned_ambiguous;
		my_read_counter->unassigned_nofeatures += global_context -> thread_contexts[xk1].read_counters.unassigned_nofeatures;
		my_read_counter->unassigned_unmapped += global_context -> thread_contexts[xk1].read_counters.unassigned_unmapped;
		my_read_counter->unassigned_mappingquality += global_context -> thread_contexts[xk1].read_counters.unassigned_mappingquality;
		my_read_counter->unassigned_fragmentlength += global_context -> thread_contexts[xk1].read_counters.unassigned_fragmentlength;
		my_read_counter->unassigned_chimericreads += global_context -> thread_contexts[xk1].read_counters.unassigned_chimericreads;
		my_read_counter->unassigned_multimapping += global_context -> thread_contexts[xk1].read_counters.unassigned_multimapping;
		my_read_counter->unassigned_secondary += global_context -> thread_contexts[xk1].read_counters.unassigned_secondary;
		my_read_counter->unassigned_junction_condition += global_context -> thread_contexts[xk1].read_counters.unassigned_junction_condition;
		my_read_counter->unassigned_duplicate += global_context -> thread_contexts[xk1].read_counters.unassigned_duplicate;
		my_read_counter->assigned_reads += global_context -> thread_contexts[xk1].read_counters.assigned_reads;

		if(global_context -> do_junction_counting){
			int bucket_i;
			for(bucket_i = 0 ; bucket_i < global_context -> thread_contexts[xk1].junction_counting_table -> numOfBuckets; bucket_i++){
				KeyValuePair * cursor;
				cursor = global_context -> thread_contexts[xk1].junction_counting_table -> bucketArray[bucket_i];
				while(cursor){
					char * junckey = (char *) cursor -> key;

					void * globval = HashTableGet(junction_global_table, junckey);
					char * new_key = malloc(strlen(junckey)+1);
					strcpy(new_key, junckey);
					globval += (cursor -> value - NULL);
					HashTablePut(junction_global_table, new_key, globval);
						// new_key will be freed when it is replaced next time or when the global table is destroyed.

					cursor = cursor->next;
				}
			}

			for(bucket_i = 0 ; bucket_i < global_context -> thread_contexts[xk1].splicing_point_table -> numOfBuckets; bucket_i++){
				KeyValuePair * cursor;
				cursor = global_context -> thread_contexts[xk1].splicing_point_table -> bucketArray[bucket_i];
				while(cursor){
					char * junckey = (char *) cursor -> key;
					void * globval = HashTableGet(splicing_global_table, junckey);
					char * new_key = malloc(strlen(junckey)+1);
					strcpy(new_key, junckey);

					//if(xk1>0)
					//SUBREADprintf("MERGE THREAD-%d : %s    VAL=%u, ADD=%u\n", xk1, junckey, globval - NULL, cursor -> value - NULL);

					globval += (cursor -> value - NULL);
					HashTablePut(splicing_global_table, new_key, globval);
					cursor = cursor->next;
				}
			}
		}
	}

	char pct_str[10];
	if(total_input_reads>0)
		sprintf(pct_str,"(%.1f%%%%)", (*nreads_mapped_to_exon)*100./total_input_reads);
	else	pct_str[0]=0;

	if(unpaired_fragment_no){
		print_in_box(80,0,0,"   Not properly paired fragments : %llu", unpaired_fragment_no);
	}
	print_in_box(80,0,0,"   Total %s : %llu", global_context -> is_paired_end_mode_assign?"fragments":"reads", total_input_reads); 
	print_in_box(pct_str[0]?81:80,0,0,"   Successfully assigned %s : %llu %s", global_context -> is_paired_end_mode_assign?"fragments":"reads", *nreads_mapped_to_exon,pct_str); 
	print_in_box(80,0,0,"   Running time : %.2f minutes", (miltime() - global_context -> start_time)/60);
	print_in_box(80,0,0,"");
}

HashTable * load_alias_table(char * fname)
{
	FILE * fp = f_subr_open(fname, "r");
	if(!fp)
	{
		print_in_box(80,0,0,"WARNING unable to open alias file '%s'", fname);
		return NULL;
	}

	char * fl = malloc(2000);

	HashTable * ret = HashTableCreate(1013);
	HashTableSetDeallocationFunctions(ret, free, free);
	HashTableSetKeyComparisonFunction(ret, fc_strcmp);
	HashTableSetHashFunction(ret, fc_chro_hash);
	
	while (1)
	{
		char *ret_fl = fgets(fl, 1999, fp);
		if(!ret_fl) break;
		if(fl[0]=='#') continue;
		char * sam_chr = NULL;
		char * anno_chr = strtok_r(fl, ",", &sam_chr);
		if((!sam_chr)||(!anno_chr)) continue;

		sam_chr[strlen(sam_chr)-1]=0;
		char * anno_chr_buf = malloc(strlen(anno_chr)+1);
		strcpy(anno_chr_buf, anno_chr);
		char * sam_chr_buf = malloc(strlen(sam_chr)+1);
		strcpy(sam_chr_buf, sam_chr);
		
		//printf("ALIAS: %s -> %s\n", sam_chr, anno_chr);
		HashTablePut(ret, sam_chr_buf, anno_chr_buf);
	}

	fclose(fp);

	free(fl);
	return ret;
}

void get_temp_dir_from_out(char * tmp, char * out){
	char * slash = strrchr(out,'/');
	if(NULL == slash){
		strcpy(tmp, "./");
	}else{
		memcpy(tmp, out, slash - out);
		tmp[slash - out]=0;
	}
}

void fc_thread_init_global_context(fc_thread_global_context_t * global_context, unsigned int buffer_size, unsigned short threads, int line_length , int is_PE_data, int min_pe_dist, int max_pe_dist, int is_gene_level, int is_overlap_allowed, int is_strand_checked, char * output_fname, int is_sam_out, int is_both_end_required, int is_chimertc_disallowed, int is_PE_distance_checked, char *feature_name_column, char * gene_id_column, int min_map_qual_score, int is_multi_mapping_allowed, int is_SAM, char * alias_file_name, char * cmd_rebuilt, int is_input_file_resort_needed, int feature_block_size, int isCVersion, int fiveEndExtension,  int threeEndExtension, int minFragmentOverlap, int is_split_or_exonic_only, int reduce_5_3_ends_to_one, char * debug_command, int is_duplicate_ignored, int is_not_sort, int use_fraction_multimapping, int useOverlappingBreakTie, char * pair_orientations, int do_junction_cnt, int max_M, int isRestrictlyNoOvelrapping, float fracOverlap, char * temp_dir)
{
	int x1;

	memset(global_context, 0, sizeof(fc_thread_global_context_t));
	global_context -> max_BAM_header_size = buffer_size;
	global_context -> all_reads = 0;
	global_context -> redo = 0;
	global_context -> SAM_output_fp = NULL;

	global_context -> isCVersion = isCVersion;
	global_context -> is_read_details_out = is_sam_out;
	global_context -> is_multi_overlap_allowed = is_overlap_allowed;
	global_context -> restricted_no_multi_overlap = isRestrictlyNoOvelrapping;
	global_context -> is_paired_end_mode_assign = is_PE_data;
	global_context -> is_gene_level = is_gene_level;
	global_context -> is_strand_checked = is_strand_checked;
	global_context -> is_both_end_required = is_both_end_required;
	global_context -> is_chimertc_disallowed = is_chimertc_disallowed;
	global_context -> is_PE_distance_checked = is_PE_distance_checked;
	global_context -> is_multi_mapping_allowed = is_multi_mapping_allowed;
	global_context -> is_split_or_exonic_only = is_split_or_exonic_only;
	global_context -> is_duplicate_ignored = is_duplicate_ignored;
	//global_context -> is_first_read_reversed = (pair_orientations[0]=='r');
	//global_context -> is_second_read_straight = (pair_orientations[1]=='f');

	global_context -> reduce_5_3_ends_to_one = reduce_5_3_ends_to_one;
	global_context -> do_not_sort = is_not_sort;
	global_context -> is_SAM_file = is_SAM;
	global_context -> use_fraction_multi_mapping = use_fraction_multimapping;
	global_context -> do_junction_counting = do_junction_cnt;

	global_context -> thread_number = threads;
	global_context -> min_mapping_quality_score = min_map_qual_score;
	global_context -> unistr_buffer_size = 1024*1024*2;
	global_context -> unistr_buffer_used = 0;
	global_context -> unistr_buffer_space = malloc(global_context -> unistr_buffer_size);
	global_context -> annot_chro_name_alias_table = NULL;
	global_context -> cmd_rebuilt = cmd_rebuilt;
	global_context -> feature_block_size = feature_block_size;
	global_context -> five_end_extension = fiveEndExtension;
	global_context -> three_end_extension = threeEndExtension;
	global_context -> fragment_minimum_overlapping = minFragmentOverlap;
	global_context -> fractional_minimum_overlapping = fracOverlap;
	global_context -> use_overlapping_break_tie = useOverlappingBreakTie;
	global_context -> need_calculate_fragment_len = ( global_context -> fractional_minimum_overlapping > 1E-10 );
	global_context -> need_calculate_overlap_len = global_context -> fractional_minimum_overlapping > 1E-10 || (global_context -> fragment_minimum_overlapping > 1) || global_context -> use_overlapping_break_tie;
	global_context -> debug_command = debug_command;
	global_context -> max_M = max_M;
	global_context -> max_BAM_header_size = buffer_size;

	global_context -> read_counters.unassigned_ambiguous=0;
	global_context -> read_counters.unassigned_nofeatures=0;
	global_context -> read_counters.unassigned_unmapped=0;
	global_context -> read_counters.unassigned_mappingquality=0;
	global_context -> read_counters.unassigned_fragmentlength=0;
	global_context -> read_counters.unassigned_chimericreads=0;
	global_context -> read_counters.unassigned_multimapping=0;
	global_context -> read_counters.unassigned_secondary=0;
	global_context -> read_counters.unassigned_junction_condition=0;
	global_context -> read_counters.unassigned_duplicate=0;
	global_context -> read_counters.assigned_reads=0;
	
	if(alias_file_name && alias_file_name[0])
	{
		strcpy(global_context -> alias_file_name,alias_file_name);
		global_context -> annot_chro_name_alias_table = load_alias_table(alias_file_name);
	}
	else	global_context -> alias_file_name[0]=0;

	strcpy(global_context -> feature_name_column,feature_name_column);
	strcpy(global_context -> gene_id_column,gene_id_column);
	strcpy(global_context -> output_file_name, output_fname);
	global_context -> output_file_path[0]=0;
	for( x1 = strlen(output_fname)-1; x1 >= 0; x1 --){
		if(output_fname[x1]=='/'){
			memcpy(global_context -> output_file_path, output_fname, x1);
			global_context -> output_file_path[x1]=0;
			break;
		}
	}
	if(0 == global_context -> output_file_path[0]){
		strcpy(global_context -> output_file_path, ".");
	}

	if(temp_dir == NULL)get_temp_dir_from_out(global_context -> temp_file_dir, output_fname);
	else strcpy(global_context -> temp_file_dir, temp_dir);
	//SUBREADprintf("OFPP:%s, OFNN:%s\n", global_context -> output_file_path, global_context -> output_file_name);

	global_context -> min_paired_end_distance = min_pe_dist;
	global_context -> max_paired_end_distance = max_pe_dist;
	global_context -> thread_number = threads;
	global_context -> line_length = line_length;
}


void pairer_unsorted_notification(void * pairer_vp, char * bin1, char * bin2){
	print_in_box(80,0,0,"");
	print_in_box(80,0,PRINT_BOX_NOCOLOR_FOR_COLON,"   WARNING: reads from the same pair were found not adjacent to each");
	print_in_box(80,0,0,"            other in the input (due to read sorting by location or");
	print_in_box(80,0,0,"            reporting of multi-mapping read pairs).");
	print_in_box(80,0,0,"");
	print_in_box(80,0,0,"   Read re-ordering is performed.");
	print_in_box(80,0,0,"");
}



int fc_thread_start_threads(fc_thread_global_context_t * global_context, int et_exons, int * et_geneid, char ** et_chr, long * et_start, long * et_stop, unsigned char * et_strand, char * et_anno_chr_2ch, char ** et_anno_chrs, long * et_anno_chr_heads, long * et_bk_end_index, long * et_bk_min_start, long * et_bk_max_end, int read_length)
{
	int xk1;

	global_context -> read_length = read_length;
	global_context -> is_unpaired_warning_shown = 0;
	global_context -> is_stake_warning_shown = 0;

	if(global_context -> is_read_details_out)
	{
		char tmp_fname[350], *modified_fname;
		int i=0;
		if( global_context -> input_file_unique ){
			sprintf(tmp_fname, "%s/%s.featureCounts", global_context -> output_file_path, global_context -> input_file_short_name);
			global_context -> SAM_output_fp = f_subr_open(tmp_fname, "w");
			//SUBREADprintf("FCSSF=%s\n", tmp_fname);
		} else {
			sprintf(tmp_fname, "%s.featureCounts", global_context -> raw_input_file_name);
			modified_fname = tmp_fname;
			while(modified_fname[0]=='/' || modified_fname[0]=='.' || modified_fname[0]=='\\'){
				modified_fname ++;
			}
			while(modified_fname[i]){
				if(modified_fname[i]=='\\' || modified_fname[i]=='/'||modified_fname[i]==' ')modified_fname[i]='.';
				i++;
			}
			char tmp_fname2[350];
			sprintf(tmp_fname2, "%s/%s",  global_context -> output_file_path, modified_fname);
			global_context -> SAM_output_fp = f_subr_open(tmp_fname2, "w");
			//SUBREADprintf("FCSSF=%s\n", tmp_fname2);
		}
		if(!global_context -> SAM_output_fp)
		{
			SUBREADprintf("Unable to create file '%s'; the read assignment details are not written.\n", tmp_fname);
		}
	}
	else
		global_context -> SAM_output_fp = NULL;

	global_context -> redo = 0;
	global_context -> exontable_geneid = et_geneid;
	global_context -> exontable_chr = et_chr;
	global_context -> exontable_start = et_start;
	global_context -> exontable_stop = et_stop;
	global_context -> exontable_strand = (char *)et_strand;
	global_context -> exontable_anno_chr_2ch = et_anno_chr_2ch;
	global_context -> exontable_anno_chrs = et_anno_chrs;
	global_context -> exontable_anno_chr_heads = et_anno_chr_heads;
	global_context -> exontable_block_end_index = et_bk_end_index;
	global_context -> exontable_block_max_end = et_bk_max_end;
	global_context -> exontable_block_min_start = et_bk_min_start;
	global_context -> sambam_chro_table_items = 0;
	global_context -> sambam_chro_table = NULL;
	pthread_spin_init(&global_context->sambam_chro_table_lock, PTHREAD_PROCESS_PRIVATE);

	global_context -> is_all_finished = 0;
	global_context -> thread_contexts = malloc(sizeof(fc_thread_thread_context_t) * global_context -> thread_number);
	for(xk1=0; xk1<global_context -> thread_number; xk1++)
	{
	//	printf("CHRR_MALLOC\n");
		global_context -> thread_contexts[xk1].thread_id = xk1;
		global_context -> thread_contexts[xk1].chunk_read_ptr = 0;
		global_context -> thread_contexts[xk1].count_table = calloc(sizeof(read_count_type_t), et_exons);
		global_context -> thread_contexts[xk1].nreads_mapped_to_exon = 0;
		global_context -> thread_contexts[xk1].all_reads = 0;
		global_context -> thread_contexts[xk1].chro_name_buff = malloc(CHROMOSOME_NAME_LENGTH);
		global_context -> thread_contexts[xk1].strm_buffer = malloc(sizeof(z_stream));

		global_context -> thread_contexts[xk1].unpaired_fragment_no = 0;
		global_context -> thread_contexts[xk1].read_counters.assigned_reads = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_ambiguous = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_nofeatures = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_unmapped = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_mappingquality = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_fragmentlength = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_chimericreads = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_multimapping = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_secondary = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_junction_condition = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_duplicate = 0;

		if(global_context -> need_calculate_overlap_len){
			global_context -> thread_contexts[xk1].scoring_buff_gap_chros = malloc( sizeof(char *) * MAX_HIT_NUMBER * 2 * global_context -> max_M *2);
			global_context -> thread_contexts[xk1].scoring_buff_gap_starts = malloc( sizeof(unsigned int ) * MAX_HIT_NUMBER * 2 * global_context -> max_M *2);
			global_context -> thread_contexts[xk1].scoring_buff_gap_lengths = malloc( sizeof(unsigned short) * MAX_HIT_NUMBER * 2 * global_context -> max_M *2);
		} else global_context -> thread_contexts[xk1].scoring_buff_gap_chros = NULL;

		if(global_context -> do_junction_counting)
		{
			global_context -> thread_contexts[xk1].junction_counting_table = HashTableCreate(131317);
			HashTableSetHashFunction(global_context -> thread_contexts[xk1].junction_counting_table,HashTableStringHashFunction);
			HashTableSetDeallocationFunctions(global_context -> thread_contexts[xk1].junction_counting_table, free, NULL);
			HashTableSetKeyComparisonFunction(global_context -> thread_contexts[xk1].junction_counting_table, fc_strcmp_chro);
			
			global_context -> thread_contexts[xk1].splicing_point_table = HashTableCreate(131317);
			HashTableSetHashFunction(global_context -> thread_contexts[xk1].splicing_point_table,HashTableStringHashFunction);
			HashTableSetDeallocationFunctions(global_context -> thread_contexts[xk1].splicing_point_table, free, NULL);
			HashTableSetKeyComparisonFunction(global_context -> thread_contexts[xk1].splicing_point_table, fc_strcmp_chro);
		}

		if(!global_context ->  thread_contexts[xk1].count_table) return 1;
		void ** thread_args = malloc(sizeof(void *)*2);
		thread_args[0] = global_context;
		thread_args[1] = & global_context -> thread_contexts[xk1];
	}

	char rand_prefix[500];
	char MAC_or_random[13];
	mac_or_rand_str(MAC_or_random);
	sprintf(rand_prefix, "%s/temp-core-%06u-%s.sam", global_context -> temp_file_dir, getpid(), MAC_or_random);

	//#warning "REMOVE ' * 2 ' FROM NEXT LINE !!!!!!"
	SAM_pairer_create(&global_context -> read_pairer, global_context -> thread_number , global_context -> max_BAM_header_size/1024/1024+2, !global_context-> is_SAM_file, 1, !global_context -> is_paired_end_mode_assign, global_context ->is_paired_end_mode_assign && global_context -> do_not_sort ,0, global_context -> input_file_name, process_pairer_reset, process_pairer_header, process_pairer_output, rand_prefix, global_context);
	SAM_pairer_set_unsorted_notification(&global_context -> read_pairer, pairer_unsorted_notification);

	return 0;
}

void fc_thread_destroy_thread_context(fc_thread_global_context_t * global_context)
{
	int xk1;
	if(global_context -> is_read_details_out)
	{
		fclose(global_context -> SAM_output_fp);
		global_context -> SAM_output_fp = NULL;
	}

	for(xk1=0; xk1<global_context-> thread_number; xk1++)
	{
		//printf("CHRR_FREE\n");
		free(global_context -> thread_contexts[xk1].count_table);	
		free(global_context -> thread_contexts[xk1].chro_name_buff);
		free(global_context -> thread_contexts[xk1].strm_buffer);
		if(global_context -> thread_contexts[xk1].scoring_buff_gap_chros){
			free(global_context -> thread_contexts[xk1].scoring_buff_gap_chros);
			free(global_context -> thread_contexts[xk1].scoring_buff_gap_starts);
			free(global_context -> thread_contexts[xk1].scoring_buff_gap_lengths);
		}
		if(global_context -> do_junction_counting){
			HashTableDestroy(global_context -> thread_contexts[xk1].junction_counting_table);
			HashTableDestroy(global_context -> thread_contexts[xk1].splicing_point_table);
		}
	}

	pthread_spin_destroy(&global_context->sambam_chro_table_lock);
	free(global_context -> thread_contexts);
}
void fc_thread_wait_threads(fc_thread_global_context_t * global_context)
{
	global_context -> is_input_bad_format |= SAM_pairer_run(&global_context -> read_pairer); 
}

void BUFstrcat(char * targ, char * src, char ** buf){
	int srclen = strlen(src);
	if( (*buf) == NULL){
		(*buf) = targ;
	}
	memcpy((*buf), src, srclen);
	(*buf) += srclen;
	(**buf) = 0;
}

void fc_write_final_gene_results(fc_thread_global_context_t * global_context, int * et_geneid, char ** et_chr, long * et_start, long * et_stop, unsigned char * et_strand, const char * out_file, int features, read_count_type_t ** column_numbers, char * file_list, int n_input_files, fc_feature_info_t * loaded_features, int header_out)
{
	int xk1;
	int genes = global_context -> gene_name_table -> numOfElements;
	read_count_type_t *gene_columns;

	FILE * fp_out = f_subr_open(out_file,"w");
	if(!fp_out){
		SUBREADprintf("Failed to create file %s\n", out_file);
		return;
	}

	if(header_out)
	{
		fprintf(fp_out, "# Program:featureCounts v%s", SUBREAD_VERSION);
		if(global_context->cmd_rebuilt)
			fprintf(fp_out, "; Command:%s", global_context->cmd_rebuilt);
		fprintf(fp_out, "\n");
	}

	char * tmp_ptr = NULL, * next_fn;
	int non_empty_files = 0, i_files=0;
	fprintf(fp_out,"Geneid\tChr\tStart\tEnd\tStrand\tLength");
	next_fn = strtok_r(file_list, ";", &tmp_ptr);
	while(1){
		if(!next_fn||strlen(next_fn)<1) break;
		if(column_numbers[i_files])
		{
			fprintf(fp_out,"\t%s", next_fn);
			non_empty_files ++;
		}
		next_fn = strtok_r(NULL, ";", &tmp_ptr);
		i_files++;
	}
	fprintf(fp_out,"\n");

	gene_columns = calloc(sizeof(read_count_type_t) , genes * non_empty_files);
	unsigned int * gene_exons_number = calloc(sizeof(unsigned int) , genes);
	unsigned int * gene_exons_pointer = calloc(sizeof(unsigned int) , genes);
	unsigned int * gene_exons_start = malloc(sizeof(unsigned int) * features);
	unsigned int * gene_exons_end = malloc(sizeof(unsigned int) * features);
	char ** gene_exons_chr = malloc(sizeof(char *) * features);
	char * gene_exons_strand = malloc(features);

	for(xk1 = 0; xk1 < features; xk1++)
	{
		int gene_id = et_geneid[xk1];
		gene_exons_number[gene_id]++;
	}

	unsigned int accumulative_no = 0;
	unsigned longest_gene_exons = 0;
	for(xk1 = 0 ; xk1 < genes; xk1++)
	{
		unsigned int tmpv = gene_exons_number[xk1];
		longest_gene_exons = max(longest_gene_exons, tmpv);
		accumulative_no += gene_exons_number[xk1];
		gene_exons_number[xk1] = accumulative_no - tmpv;
	}

	for(xk1 = 0; xk1 < features; xk1++)
	{
		int gene_id = et_geneid[xk1];
		int gene_write_ptr = gene_exons_number[gene_id] + gene_exons_pointer[gene_id];

		gene_exons_chr[gene_write_ptr] = et_chr[xk1];
		gene_exons_start[gene_write_ptr] = et_start[xk1]; 
		gene_exons_end[gene_write_ptr] = et_stop[xk1]; 
		gene_exons_strand[gene_write_ptr] = et_strand[xk1]; 

		gene_exons_pointer[gene_id]++;
	}

	for(xk1 = 0; xk1 < features; xk1++)
	{
		int gene_id = et_geneid[xk1], k_noempty = 0;
		for(i_files=0;i_files < n_input_files; i_files++)
		{
			if(column_numbers[i_files]==NULL) continue;
			gene_columns[gene_id * non_empty_files + k_noempty ] += column_numbers[i_files][xk1];
			k_noempty++;
		}
	}


	char *is_occupied = malloc(longest_gene_exons);
	unsigned int * input_start_stop_list = malloc(longest_gene_exons * sizeof(int) * 2);
	unsigned int * output_start_stop_list = malloc(longest_gene_exons * sizeof(int) * 2);

	char * out_chr_list = malloc(longest_gene_exons * (1+global_context -> longest_chro_name) + 1), * tmp_chr_list = NULL;
	char * out_start_list = malloc(11 * longest_gene_exons + 1), * tmp_start_list = NULL;
	char * out_end_list = malloc(11 * longest_gene_exons + 1), * tmp_end_list = NULL;
	char * out_strand_list = malloc(2 * longest_gene_exons + 1), * tmp_strand_list = NULL;

	for(xk1 = 0 ; xk1 < genes; xk1++)
	{
		int xk2;
		
		memset(is_occupied,0,gene_exons_pointer[xk1]);
		tmp_chr_list = NULL;
		tmp_start_list = NULL;
		tmp_end_list = NULL;
		tmp_strand_list = NULL;
		out_chr_list[0]=0;
		out_start_list[0]=0;
		out_end_list[0]=0;
		out_strand_list[0]=0;
		int gene_nonoverlap_len =0;

		for(xk2=0; xk2<gene_exons_pointer[xk1]; xk2++)
		{
			if(!is_occupied[xk2])
			{
				int xk3;
				char * matched_chr = gene_exons_chr[xk2 + gene_exons_number[xk1]];
				char matched_strand = gene_exons_strand[xk2 + gene_exons_number[xk1]];

				memset(input_start_stop_list, 0, gene_exons_pointer[xk1] * sizeof(int) * 2);
				int gap_merge_ptr = 1;
				input_start_stop_list[0] = gene_exons_start[xk2 + gene_exons_number[xk1]];
				input_start_stop_list[1] = gene_exons_end[xk2 + gene_exons_number[xk1]] + 1;

				for(xk3 = xk2+1; xk3 < gene_exons_pointer[xk1]; xk3++)
				{
					if((!is_occupied[xk3]) && strcmp(matched_chr, gene_exons_chr[xk3+gene_exons_number[xk1]])==0 && matched_strand == gene_exons_strand[xk3 + gene_exons_number[xk1]])
					{
						is_occupied[xk3]=1;
						input_start_stop_list[gap_merge_ptr*2] = gene_exons_start[xk3+gene_exons_number[xk1]]; 
						input_start_stop_list[gap_merge_ptr*2+1] = gene_exons_end[xk3+gene_exons_number[xk1]]+1;
						gap_merge_ptr++;
					}
				}

				{
						int merged_gaps = mergeIntervals(input_start_stop_list, output_start_stop_list, gap_merge_ptr);

						for(xk3=0; xk3<gap_merge_ptr; xk3++)
						{
							char numbbuf[12];
							BUFstrcat(out_chr_list, matched_chr, &tmp_chr_list);
							BUFstrcat(out_chr_list, ";", &tmp_chr_list);

							sprintf(numbbuf,"%u;", input_start_stop_list[xk3 * 2]);
							BUFstrcat(out_start_list, numbbuf, &tmp_start_list);
							sprintf(numbbuf,"%u;", input_start_stop_list[xk3 * 2 + 1] - 1);
							BUFstrcat(out_end_list, numbbuf, &tmp_end_list);
							sprintf(numbbuf,"%c;", matched_strand?'-':'+');
							BUFstrcat(out_strand_list, numbbuf, &tmp_strand_list);

						}
						for(xk3=0; xk3<merged_gaps; xk3++)
							gene_nonoverlap_len += output_start_stop_list[xk3 * 2 + 1] - output_start_stop_list[xk3 * 2];
				}	
			}
		}

		unsigned char * gene_symbol = global_context -> gene_name_array [xk1];

		#define _cut_tail(x) (x)[strlen(x)-1]=0

		_cut_tail(out_chr_list);
		_cut_tail(out_start_list);
		_cut_tail(out_end_list);
		_cut_tail(out_strand_list);

		fprintf(fp_out, "%s\t%s\t%s\t%s\t%s\t%d"    , gene_symbol, out_chr_list, out_start_list, out_end_list, out_strand_list, gene_nonoverlap_len);

		// all exons: gene_exons_number[xk1] : gene_exons_pointer[xk1]
		int non_empty_file_index = 0;
		for(i_files=0; i_files< n_input_files; i_files++)
		{
			if(column_numbers[i_files])
			{
				read_count_type_t longlong_res = 0;
				double double_res = 0;
				int is_double_number = calc_float_fraction(gene_columns[non_empty_file_index+non_empty_files*xk1], &longlong_res, &double_res);
				if(is_double_number){
					fprintf(fp_out,"\t%.2f", double_res);
				}else{
					fprintf(fp_out,"\t%llu", longlong_res);
				}
				non_empty_file_index ++;
			}
		}
		fprintf(fp_out,"\n");

	}
	free(is_occupied);
	free(input_start_stop_list);
	free(output_start_stop_list);
	free(out_chr_list);
	free(out_strand_list);
	free(out_start_list);
	free(out_end_list);

	free(gene_exons_number);
	free(gene_exons_pointer);
	free(gene_columns);
	free(gene_exons_chr);
	free(gene_exons_start);
	free(gene_exons_end);
	free(gene_exons_strand);
	fclose(fp_out);
}

void fc_write_final_counts(fc_thread_global_context_t * global_context, const char * out_file, int nfiles, char * file_list, read_count_type_t ** column_numbers, fc_read_counters *read_counters, int isCVersion)
{
	char fname[300];
	int i_files, xk1;

	sprintf(fname, "%s.summary", out_file);
	FILE * fp_out = f_subr_open(fname,"w");

	if(!fp_out){
		SUBREADprintf("Unable to create summary file '%s'\n", fname);
		return;
	}

	fprintf(fp_out,"Status");
	char * next_fn = file_list;
	
	for(i_files=0; i_files<nfiles; i_files++)
	{
		if(!next_fn||strlen(next_fn)<1) break;
		if(column_numbers[i_files])
			fprintf(fp_out,"\t%s", next_fn);

		next_fn += strlen(next_fn)+1;
	}

	fprintf(fp_out,"\n");
	char * keys [] ={ "Assigned" , "Unassigned_Ambiguity", "Unassigned_MultiMapping" ,"Unassigned_NoFeatures", "Unassigned_Unmapped", "Unassigned_MappingQuality", "Unassigned_FragmentLength", "Unassigned_Chimera", "Unassigned_Secondary", (global_context->is_split_or_exonic_only == 2)?"Unassigned_Hasjunction":"Unassigned_Nonjunction", "Unassigned_Duplicate"};

	for(xk1=0; xk1<11; xk1++)
	{
		fprintf(fp_out,"%s", keys[xk1]);
		for(i_files = 0; i_files < nfiles; i_files ++)
		{
			unsigned long long * array_0 = (unsigned long long *)&(read_counters[i_files]);
			unsigned long long * cntr = array_0 + xk1;
			if(column_numbers[i_files])
				fprintf(fp_out,"\t%llu", *cntr);
		}
		fprintf(fp_out,"\n");
	}


	fclose(fp_out);
}
void fc_write_final_results(fc_thread_global_context_t * global_context, const char * out_file, int features, read_count_type_t ** column_numbers, char * file_list, int n_input_files, fc_feature_info_t * loaded_features, int header_out)
{
	/* save the results */
	FILE * fp_out;
	int i, i_files = 0;
	fp_out = f_subr_open(out_file,"w");
	if(!fp_out){
		SUBREADprintf("Failed to create file %s\n", out_file);
			return;
		}

	if(header_out)
	{
		fprintf(fp_out, "# Program:featureCounts v%s", SUBREAD_VERSION);
		if(global_context->cmd_rebuilt)
			fprintf(fp_out, "; Command:%s", global_context->cmd_rebuilt);
		fprintf(fp_out, "\n");
	}



	char * tmp_ptr = NULL, * next_fn;
	fprintf(fp_out,"Geneid\tChr\tStart\tEnd\tStrand\tLength");
	next_fn = strtok_r(file_list, ";", &tmp_ptr);
	while(1){
		if(!next_fn||strlen(next_fn)<1) break;
		if(column_numbers[i_files])
			fprintf(fp_out,"\t%s", next_fn);
		next_fn = strtok_r(NULL, ";", &tmp_ptr);
		i_files++;
	}
	fprintf(fp_out,"\n");
	for(i=0;i<features;i++)
	{
		fprintf(fp_out,"%s\t%s\t%u\t%u\t%c\t%d", global_context -> unistr_buffer_space + loaded_features[i].feature_name_pos,
 							   global_context -> unistr_buffer_space + loaded_features[i].feature_name_pos + loaded_features[i].chro_name_pos_delta,
						   	   loaded_features[i].start, loaded_features[i].end, loaded_features[i].is_negative_strand?'-':'+',loaded_features[i].end-loaded_features[i].start+1);
		for(i_files=0; i_files<n_input_files; i_files++)
		{
			if(column_numbers[i_files])
			{
				int sorted_exon_no = loaded_features[i].sorted_order;
				unsigned long long count_frac_raw =  column_numbers[i_files][sorted_exon_no], longlong_res = 0;

				double double_res = 0;
				int is_double_number = calc_float_fraction(count_frac_raw, &longlong_res, &double_res);
				if(is_double_number){
					fprintf(fp_out,"\t%.2f", double_res);
				}else{
					fprintf(fp_out,"\t%llu", longlong_res);
				}


			}
		}
		fprintf(fp_out,"\n");
	}

	fclose(fp_out);
}

static struct option long_options[] =
{
	{"primary",no_argument, 0, 0},
	{"readExtension5", required_argument, 0, 0},
	{"readExtension3", required_argument, 0, 0},
	{"read2pos", required_argument, 0, 0},
	{"minOverlap", required_argument, 0, 0},
	{"fracOverlap", required_argument, 0, 0},
	{"splitOnly", no_argument, 0, 0},
	{"nonSplitOnly", no_argument, 0, 0},
	{"debugCommand", required_argument, 0, 0},
	{"ignoreDup", no_argument, 0, 0},
	{"donotsort", no_argument, 0, 0},
	{"restrictedlyNoOverlap", no_argument, 0, 0},
	{"fraction", no_argument, 0, 0},
	{"order", required_argument, 0, 'S'},
	{"genome", required_argument, 0, 'G'},
	{"maxMOp", required_argument, 0, 0},
	{"tmpDir", required_argument, 0, 0},
	{"largestOverlap", no_argument, 0,0},
	{0, 0, 0, 0}
};

void print_usage()
{
	SUBREADprintf("\nVersion %s\n\n", SUBREAD_VERSION);

	SUBREADputs("Usage: featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2] ... \n");
	SUBREADputs("## Required arguments:");
	SUBREADputs("");
	SUBREADputs("  -a <string>         Name of an annotation file. GTF/GFF format by default.");
	SUBREADputs("                      See -F option for more formats.");
	SUBREADputs("");
	SUBREADputs("  -o <string>         Name of the output file including read counts. A separate"); 
	SUBREADputs("                      file including summary statistics of counting results is");
	SUBREADputs("                      also included in the output (`<string>.summary')");
	SUBREADputs("");
	SUBREADputs("  input_file1 [input_file2] ...   A list of SAM or BAM format files.");
	SUBREADputs("");
	SUBREADputs("## Options:");
	SUBREADputs("# Annotation");
	SUBREADputs("");
	SUBREADputs("  -F <string>         Specify format of provided annotation file. Acceptable");
	SUBREADputs("                      formats include `GTF/GFF' and `SAF'. `GTF/GFF' by default.");
	SUBREADputs("                      See Users Guide for description of SAF format.");
	SUBREADputs("");
	SUBREADputs("  -t <string>         Specify feature type in GTF annotation. `exon' by ");
	SUBREADputs("                      default. Features used for read counting will be ");
	SUBREADputs("                      extracted from annotation using the provided value.");
	SUBREADputs("");
	SUBREADputs("  -g <string>         Specify attribute type in GTF annotation. `gene_id' by ");
	SUBREADputs("                      default. Meta-features used for read counting will be ");
	SUBREADputs("                      extracted from annotation using the provided value.");
	SUBREADputs("");
	SUBREADputs("  -A <string>         Provide a chromosome name alias file to match chr names in");
	SUBREADputs("                      annotation with those in the reads. This should be a two-");
	SUBREADputs("                      column comma-delimited text file. Its first column should");
	SUBREADputs("                      include chr names in the annotation and its second column");
	SUBREADputs("                      should include chr names in the reads. Chr names are case");
	SUBREADputs("                      sensitive. No column header should be included in the");
	SUBREADputs("                      file.");
	SUBREADputs("");

	SUBREADputs("# Level of summarization");
	SUBREADputs("");
	SUBREADputs("  -f                  Perform read counting at feature level (eg. counting ");
	SUBREADputs("                      reads for exons rather than genes).");
	SUBREADputs("");

	SUBREADputs("# Overlap between reads and features");
	SUBREADputs("");
	SUBREADputs("  -O                  Assign reads to all their overlapping meta-features (or ");
	SUBREADputs("                      features if -f is specified).");
	SUBREADputs("");
	SUBREADputs("  --minOverlap <int>  Minimum number of overlapping bases in a read that is");
	SUBREADputs("                      required for read assignment. 1 by default. Number of");
	SUBREADputs("                      overlapping bases is counted from both reads if paired");
	SUBREADputs("                      end. If a negative value is provided, then a gap of up");
	SUBREADputs("                      to specified size will be allowed between read and the");
	SUBREADputs("                      feature that the read is assigned to.");
        SUBREADputs("");
	SUBREADputs("  --fracOverlap <value> Minimum fraction of overlapping bases in a read that is");
	SUBREADputs("                      required for read assignment. Value should be within range");
	SUBREADputs("                      [0,1]. 0 by default. Number of overlapping bases is");
	SUBREADputs("                      counted from both reads if paired end. Both this option");
	SUBREADputs("                      and '--minOverlap' option need to be satisfied for read");
	SUBREADputs("                      assignment.");
	SUBREADputs("");
	SUBREADputs("  --largestOverlap    Assign reads to a meta-feature/feature that has the ");
	SUBREADputs("                      largest number of overlapping bases.");
	SUBREADputs("");
	SUBREADputs("  --readExtension5 <int> Reads are extended upstream by <int> bases from their");
	SUBREADputs("                      5' end.");
	SUBREADputs("");
	SUBREADputs("  --readExtension3 <int> Reads are extended upstream by <int> bases from their");
	SUBREADputs("                      3' end.");
	SUBREADputs("");
	SUBREADputs("  --read2pos <5:3>    Reduce reads to their 5' most base or 3' most base. Read");
	SUBREADputs("                      counting is then performed based on the single base the ");
	SUBREADputs("                      read is reduced to.");
	SUBREADputs("");

	SUBREADputs("# Multi-mapping reads");
	SUBREADputs("");
	SUBREADputs("  -M                  Multi-mapping reads will also be counted. For a multi-");
	SUBREADputs("                      mapping read, all its reported alignments will be ");
	SUBREADputs("                      counted. The `NH' tag in BAM/SAM input is used to detect ");
	SUBREADputs("                      multi-mapping reads.");
	SUBREADputs("");
	SUBREADputs("# Fractional counting");
	SUBREADputs("");
	SUBREADputs("  --fraction          Assign fractional counts to features. This option must");
	SUBREADputs("                      be used together with '-M' or '-O' or both. When '-M' is");
	SUBREADputs("                      specified, each reported alignment from a multi-mapping");
	SUBREADputs("                      read (identified via 'NH' tag) will carry a fractional");
	SUBREADputs("                      count of 1/x, instead of 1 (one), where x is the total");
	SUBREADputs("                      number of alignments reported for the same read. When '-O'");
	SUBREADputs("                      is specified, each overlapping feature will receive a");
	SUBREADputs("                      fractional count of 1/y, where y is the total number of");
	SUBREADputs("                      features overlapping with the read. When both '-M' and");
	SUBREADputs("                      '-O' are specified, each alignment will carry a fraction");
	SUBREADputs("                      count of 1/(x*y).");
	SUBREADputs("");
	

	SUBREADputs("# Read filtering");
	SUBREADputs("");
	SUBREADputs("  -Q <int>            The minimum mapping quality score a read must satisfy in");
	SUBREADputs("                      order to be counted. For paired-end reads, at least one");
	SUBREADputs("                      end should satisfy this criteria. 0 by default.");
	SUBREADputs("");
	SUBREADputs("  --splitOnly         Count split alignments only (ie. alignments with CIGAR");
	SUBREADputs("                      string containing 'N'). An example of split alignments is");
	SUBREADputs("                      exon-spanning reads in RNA-seq data.");
	SUBREADputs("");
	SUBREADputs("  --nonSplitOnly      If specified, only non-split alignments (CIGAR strings do");
	SUBREADputs("                      not contain letter 'N') will be counted. All the other");
	SUBREADputs("                      alignments will be ignored.");
	SUBREADputs("");
	SUBREADputs("  --primary           Count primary alignments only. Primary alignments are ");
	SUBREADputs("                      identified using bit 0x100 in SAM/BAM FLAG field.");
	SUBREADputs("");
	SUBREADputs("  --ignoreDup         Ignore duplicate reads in read counting. Duplicate reads ");
	SUBREADputs("                      are identified using bit Ox400 in BAM/SAM FLAG field. The ");
	SUBREADputs("                      whole read pair is ignored if one of the reads is a ");
	SUBREADputs("                      duplicate read for paired end data.");
	SUBREADputs("");

	SUBREADputs("# Strandness");
	SUBREADputs("");
	SUBREADputs("  -s <int>            Perform strand-specific read counting. Acceptable values:");
	SUBREADputs("                      0 (unstranded), 1 (stranded) and 2 (reversely stranded).");
	SUBREADputs("                      0 by default.");
	SUBREADputs("");

	SUBREADputs("# Exon-exon junctions");
	SUBREADputs("");
	SUBREADputs("  -J                  Count number of reads supporting each exon-exon junction.");
	SUBREADputs("                      Junctions were identified from those exon-spanning reads");
	SUBREADputs("                      in the input (containing 'N' in CIGAR string). Counting");
	SUBREADputs("                      results are saved to a file named '<output_file>.jcounts'");
	SUBREADputs("");
	SUBREADputs("  -G <string>         Provide the name of a FASTA-format file that contains the");
	SUBREADputs("                      reference sequences used in read mapping that produced the");
	SUBREADputs("                      provided SAM/BAM files. This optional argument can be used");
	SUBREADputs("                      with '-J' option to improve read counting for junctions.");
	SUBREADputs("");

	SUBREADputs("# Parameters specific to paired end reads");
	SUBREADputs("");
	SUBREADputs("  -p                  If specified, fragments (or templates) will be counted");
	SUBREADputs("                      instead of reads. This option is only applicable for");
	SUBREADputs("                      paired-end reads.");
	SUBREADputs("");
	SUBREADputs("  -B                  Count read pairs that have both ends successfully aligned ");
	SUBREADputs("                      only.");
	SUBREADputs("");
	SUBREADputs("  -P                  Check validity of paired-end distance when counting read ");
	SUBREADputs("                      pairs. Use -d and -D to set thresholds.");
	SUBREADputs("");
	SUBREADputs("  -d <int>            Minimum fragment/template length, 50 by default.");
	SUBREADputs("");
	SUBREADputs("  -D <int>            Maximum fragment/template length, 600 by default.");
	SUBREADputs("");
	SUBREADputs("  -C                  Do not count read pairs that have their two ends mapping ");
	SUBREADputs("                      to different chromosomes or mapping to same chromosome ");
	SUBREADputs("                      but on different strands.");
	SUBREADputs("");
	SUBREADputs("  --donotsort         Do not sort reads in BAM/SAM input. Note that reads from ");
	SUBREADputs("                      the same pair are required to be located next to each ");
	SUBREADputs("                      other in the input.");
	SUBREADputs("");

	SUBREADputs("# Number of CPU threads");
	SUBREADputs("");
	SUBREADputs("  -T <int>            Number of the threads. 1 by default.");
	SUBREADputs("");
	SUBREADputs("# Miscellaneous");
	SUBREADputs("");
	SUBREADputs("  -R                  Output detailed assignment result for each read. A text ");
	SUBREADputs("                      file will be generated for each input file, including ");
	SUBREADputs("                      names of reads and meta-features/features reads were ");
	SUBREADputs("                      assigned to. See Users Guide for more details.");
	SUBREADputs("");
	SUBREADputs("  --tmpDir <string>   Directory under which intermediate files are saved (later");
	SUBREADputs("                      removed). By default, intermediate files will be saved to");
	SUBREADputs("                      the directory specified in '-o' argument.");
	SUBREADputs("");
	SUBREADputs("  --maxMOp <int>      Maximum number of 'M' operations allowed in a CIGAR");
	SUBREADputs("                      string. 10 by default. Both 'X' and '=' are treated as 'M'");
	SUBREADputs("                      and adjacent 'M' operations are merged in the CIGAR");
	SUBREADputs("                      string.");
	SUBREADputs("");
	SUBREADputs("  -v                  Output version of the program.");
	SUBREADputs("");

}

int junckey_sort_compare(void * inptr, int i, int j){
	char ** inp = (char **) inptr;
	int x1;

	int chrI=-1, chrJ=-1;

	if(atoi(inp[i])>0) chrI = atoi(inp[i]);
	if(atoi(inp[j])>0) chrJ = atoi(inp[j]);

	if(inp[i][0]=='X' && !isdigit(inp[i][1])&& !isalpha(inp[i][1])) chrI = 90;
	if(inp[i][0]=='Y' && !isdigit(inp[i][1])&& !isalpha(inp[i][1])) chrI = 91;
	if(inp[i][0]=='M' && !isdigit(inp[i][1])&& !isalpha(inp[i][1])) chrI = 99;
	if(inp[j][0]=='X' && !isdigit(inp[j][1])&& !isalpha(inp[j][1])) chrJ = 90;
	if(inp[j][0]=='Y' && !isdigit(inp[j][1])&& !isalpha(inp[j][1])) chrJ = 91;
	if(inp[j][0]=='M' && !isdigit(inp[j][1])&& !isalpha(inp[j][1])) chrJ = 99;



	if(memcmp(inp[i], "chr", 3)==0){
		chrI=atoi(inp[i]+3);
		if(0 == chrI && inp[i][3] == 'X') chrI = 90;
		if(0 == chrI && inp[i][3] == 'Y') chrI = 91;
		if(0 == chrI && inp[i][3] == 'M') chrI = 99;
	}
	if(memcmp(inp[j], "chr", 3)==0){
		chrJ=atoi(inp[j]+3);
		if(0 == chrJ && inp[j][3] == 'X') chrJ = 90;
		if(0 == chrJ && inp[j][3] == 'Y') chrJ = 91;
		if(0 == chrJ && inp[j][3] == 'M') chrJ = 99;
	}

	int len_I_long = 9;
	for(x1 = 0 ; x1 < FEATURE_NAME_LENGTH + 15 ; x1++){
		int c1 = inp[i][x1];
		int c2 = inp[j][x1];
		if(c1 == '\t' && c2 != '\t')
			len_I_long = -1;
		else if(c1 != '\t' && c2 == '\t')
			len_I_long = 1;
		else if(c1 == '\t' && c2 == '\t')
			len_I_long = 0;

		if(len_I_long != 9) break;
	}

	if(chrI != chrJ || len_I_long != 0){
		return (chrI * 100 + len_I_long) - (chrJ * 100);
	}

	for(x1 = 0 ; x1 < FEATURE_NAME_LENGTH + 15 ; x1++){
		int c1 = inp[i][x1];
		int c2 = inp[j][x1];
		if(c1 != c2){
			return c1 - c2;
		}else if(c1 == '\t' && c1 == c2){
			int pos1 = atoi(inp[i]+x1+1);
			int pos2 = atoi(inp[j]+x1+1);
			if( pos1 == pos2)
				return strcmp(inp[i], inp[j]);
			else
				return pos1 - pos2;
		}

		if(c1 == 0 || c2 == 0)return c1 - c2;
	}
	return 0;
}

void junckey_sort_exchange(void * inptr, int i, int j){

	char ** inp = (char **) inptr;
	char * tmpp = inp[j];
	inp[j]=inp[i];
	inp[i]=tmpp;
}

void junckey_sort_merge(void * inptr, int start, int items1, int items2){
	char ** inp = (char **) inptr;
	char ** tmpp = malloc(sizeof(char *) * (items1+items2));
	int read_1_ptr = start, read_2_ptr = start+items1, outptr = 0;
	while(1){
		if(read_1_ptr == start+items1 && read_2_ptr == start+items1+items2) break;
		if((read_1_ptr == start+items1)||(read_2_ptr < start+items1+items2 &&  junckey_sort_compare(inptr, read_1_ptr, read_2_ptr) > 0 )) {
			// select 2
			tmpp[outptr++]=inp[read_2_ptr++];
		} else {
			// select 1
			tmpp[outptr++]=inp[read_1_ptr++];
		}
	}
	memcpy(inp + start, tmpp, sizeof(char *)*(items1+items2));
	free(tmpp);
}

int junccmp(fc_junction_gene_t * j1, fc_junction_gene_t * j2){
	if(strcmp( j1 -> gene_name, j2 -> gene_name ) == 0)
		return 0;
	return 1;
}


void fc_write_final_junctions(fc_thread_global_context_t * global_context,  char * output_file_name,  read_count_type_t ** table_columns, char * input_file_names, int n_input_files, HashTable ** junction_global_table_list, HashTable ** splicing_global_table_list){
	int infile_i;

	HashTable * merged_junction_table = HashTableCreate(156679);

	HashTableSetHashFunction(merged_junction_table,HashTableStringHashFunction);
	HashTableSetDeallocationFunctions(merged_junction_table, NULL, NULL);
	HashTableSetKeyComparisonFunction(merged_junction_table, fc_strcmp_chro);

	HashTable * merged_splicing_table = HashTableCreate(156679);

	HashTableSetHashFunction(merged_splicing_table,HashTableStringHashFunction);
	HashTableSetDeallocationFunctions(merged_splicing_table, NULL, NULL);
	HashTableSetKeyComparisonFunction(merged_splicing_table, fc_strcmp_chro);


	for(infile_i = 0 ; infile_i < n_input_files ; infile_i ++){
		if(!table_columns[infile_i]) continue;	// bad input file
		KeyValuePair * cursor;
		int bucket;
		for(bucket=0; bucket < splicing_global_table_list[infile_i]  -> numOfBuckets; bucket++)
		{
			cursor = splicing_global_table_list[infile_i] -> bucketArray[bucket];
			while (cursor)
			{
				char * ky = (char *)cursor -> key;
				unsigned int old_supp = HashTableGet(merged_splicing_table, ky) - NULL;
				old_supp += (cursor -> value - NULL);
				HashTablePut(merged_splicing_table, ky, NULL+old_supp);
				cursor = cursor -> next;
			}	
		}
	}

	for(infile_i = 0 ; infile_i < n_input_files ; infile_i ++){
		if(!table_columns[infile_i]) continue;	// bad input file
		KeyValuePair * cursor;
		int bucket;
		for(bucket=0; bucket < junction_global_table_list[infile_i]  -> numOfBuckets; bucket++)
		{
			cursor = junction_global_table_list[infile_i] -> bucketArray[bucket];
			while (cursor)
			{
				char * ky = (char *)cursor -> key;

				if(HashTableGet(merged_junction_table, ky)==NULL)
					HashTablePut(merged_junction_table, ky, NULL+1);
				cursor = cursor -> next;
			}	
		}
	}

	char ** key_list;
	key_list = malloc(sizeof(char *) * merged_junction_table -> numOfElements);

	KeyValuePair * cursor;
	int bucket, ky_i = 0;
	for(bucket=0; bucket < merged_junction_table -> numOfBuckets; bucket++){
		cursor = merged_junction_table -> bucketArray[bucket];
		while (cursor){
			char * ky = (char *)cursor -> key;

			key_list[ky_i ++] = ky;
			cursor = cursor -> next;
		}
	}

	merge_sort(key_list,  merged_junction_table -> numOfElements , junckey_sort_compare, junckey_sort_exchange, junckey_sort_merge);

	char outfname[300];
	sprintf(outfname, "%s.jcounts", output_file_name);

	int max_junction_genes = 3000;
	char * gene_names = malloc(max_junction_genes * FEATURE_NAME_LENGTH), * gene_name_tail;
	fc_junction_gene_t ** ret_juncs_small = malloc(sizeof(fc_junction_gene_t *) * max_junction_genes);
	fc_junction_gene_t ** ret_juncs_large = malloc(sizeof(fc_junction_gene_t *) * max_junction_genes);
	fc_junction_gene_t ** junction_key_list = malloc(sizeof(fc_junction_gene_t *)* max_junction_genes * 2);
	unsigned int * junction_support_list = malloc(sizeof(int)* max_junction_genes * 2);
	unsigned char * junction_source_list = malloc(sizeof(char)* max_junction_genes * 2 );

	int ky_i1, ky_i2;
	FILE * ofp = fopen(outfname, "w");
	char * tmpp = NULL;
	char * next_fn = input_file_names;

	fprintf(ofp, "PrimaryGene\tSecondaryGenes\tSite1_chr\tSite1_location\tSite1_strand\tSite2_chr\tSite2_location\tSite2_strand");

	for(infile_i=0; infile_i < n_input_files; infile_i++)
	{
		if(!next_fn||strlen(next_fn)<1) break;
		if(table_columns[infile_i])
			fprintf(ofp,"\t%s", next_fn);

		next_fn += strlen(next_fn)+1;
	}
	fprintf(ofp, "\n");

	for(ky_i = 0; ky_i < merged_junction_table -> numOfElements ; ky_i ++){

		//SUBREADprintf("KY=%s\n", key_list[ky_i]);

		int unique_junctions = 0;
		char * chro_small = strtok_r( key_list[ky_i] , "\t", &tmpp);
		char * pos_small_str = strtok_r( NULL, "\t", &tmpp);
		char * chro_large = strtok_r( NULL, "\t", &tmpp);
		char * pos_large_str = strtok_r( NULL, "\t", &tmpp);

		unsigned int pos_small = atoi(pos_small_str);
		unsigned int pos_large = atoi(pos_large_str);

		int found_features_small = locate_junc_features(global_context, chro_small, pos_small, ret_juncs_small , max_junction_genes); 
		int found_features_large = locate_junc_features(global_context, chro_large, pos_large, ret_juncs_large , max_junction_genes);

		char * strand = "NA";
		if(global_context -> fasta_contigs){
			char donor[3], receptor[3];
			donor[2]=receptor[2]=0;
			int has = !get_contig_fasta(global_context -> fasta_contigs, chro_small, pos_small, 2, donor);
			has = has && !get_contig_fasta(global_context -> fasta_contigs, chro_large, pos_large-3, 2, receptor);
			if(has){
				if(donor[0]=='G' && donor[1]=='T' && receptor[0]=='A' && receptor[1]=='G') strand = "+";
				else if(donor[0]=='C' && donor[1]=='T' && receptor[0]=='A' && receptor[1]=='C') strand = "-";
			}else if(!global_context ->is_junction_no_chro_shown){
				global_context ->is_junction_no_chro_shown = 1;
				print_in_box(80,0,0, "   WARNING contig `%s' is not found in the", chro_small);
				print_in_box(80,0,0, "   provided genome file!");
				print_in_box(80,0,0,"");

			}
		}

		//SUBREADprintf("FOUND=%d, %d\n", found_features_small, found_features_large);

		gene_name_tail = gene_names;
		gene_names[0]=0;

		// rules to choose the primary gene:
		// (1) if some genes have one support but the other have multiple supporting reads: remove the lowly supported genes
		// (2) if all genes have only one support but from different ends of the fragment, then remove the genes that are assigned to the end having lower supporting fragments
		// (3) choose the gene that have the smallest coordinate.

		int max_supp = 0;
		for(ky_i1 = 0; ky_i1 < found_features_small + found_features_large; ky_i1++){
			int is_duplicate = 0;
			fc_junction_gene_t * tested_key = (ky_i1 < found_features_small)?ret_juncs_small[ky_i1] :ret_juncs_large[ky_i1 - found_features_small];
			for(ky_i2 = 0; ky_i2 < unique_junctions; ky_i2 ++){
				if(junccmp( tested_key, junction_key_list[ky_i2]  )==0){
					junction_support_list[ ky_i2 ] ++;
					junction_source_list[ky_i2] |= ( (ky_i1 < found_features_small)? 1 : 2 );
					is_duplicate = 1;
					break;
				}
			}

			if(!is_duplicate){
				junction_key_list[unique_junctions] = tested_key;
				junction_support_list[unique_junctions] = 1;
				junction_source_list[unique_junctions] = ( (ky_i1 < found_features_small)? 1 : 2 );
				max_supp = max(junction_support_list[unique_junctions], max_supp);
				unique_junctions++;
			}
		}

		if(1 == max_supp){
			if(found_features_small > 0 && found_features_large > 0){
				char junc_key [FEATURE_NAME_LENGTH + 15]; 
				sprintf(junc_key, "%s\t%u", chro_small, pos_small);
				unsigned int supp_small = HashTableGet(merged_splicing_table, junc_key) - NULL;
				sprintf(junc_key, "%s\t%u", chro_large, pos_large);
				unsigned int supp_large = HashTableGet(merged_splicing_table, junc_key) - NULL;

				if(supp_small !=supp_large){
					for(ky_i2 = 0; ky_i2 < unique_junctions; ky_i2 ++){
						if(supp_small > supp_large && junction_source_list[ky_i2] == 1) junction_key_list[ky_i2] = NULL;
						else if(supp_small < supp_large && junction_source_list[ky_i2] == 2) junction_key_list[ky_i2] = NULL;
					}
				} 
			}
		}

		int smallest_coordinate_gene = 0x7fffffff;
		fc_junction_gene_t * primary_gene = NULL;
		
		for(ky_i2 = 0; ky_i2 < unique_junctions; ky_i2 ++){
			fc_junction_gene_t * tested_key = junction_key_list[ky_i2];
			if(tested_key != NULL && tested_key -> pos_first_base < smallest_coordinate_gene){
				primary_gene = tested_key;
				smallest_coordinate_gene = tested_key -> pos_first_base;
			}
		}

		if(primary_gene == NULL){
			strcpy(gene_names, "NA");
		}else{
			strcpy(gene_names, primary_gene -> gene_name);
		}

		*(pos_small_str-1)='\t';
		*(pos_large_str-1)='\t';

		fprintf(ofp, "%s", gene_names);
	
		gene_name_tail = gene_names;
		gene_names[0]=0;
		for(ky_i2 = 0; ky_i2 < unique_junctions; ky_i2 ++){
			fc_junction_gene_t * tested_key = junction_key_list[ky_i2];
			if(tested_key && tested_key != primary_gene)
				gene_name_tail += sprintf(gene_name_tail, "%s,", tested_key -> gene_name);
		}
		if( gene_names[0] ) gene_name_tail[-1]=0;
		else strcpy(gene_names, "NA");
		fprintf(ofp, "\t%s", gene_names);

		fprintf(ofp, "\t%s\t%s\t%s\t%s", chro_small, strand, chro_large, strand);

		chro_large[-1]='\t';

		for(infile_i = 0 ; infile_i < n_input_files ; infile_i ++){
			if(!table_columns[infile_i]) continue;
			unsigned long count = HashTableGet(junction_global_table_list[infile_i]  , key_list[ky_i]) - NULL;
			fprintf(ofp,"\t%lu", count);
		}
		fprintf(ofp, "\n");
	}
	fclose(ofp);
	free(junction_key_list);
	free(gene_names);
	free(ret_juncs_small);
	free(ret_juncs_large);
	free(junction_support_list);
	free(key_list);
	free(junction_source_list);

	print_in_box(80,0,PRINT_BOX_CENTER,"Found %llu junctions in all the input files.", merged_junction_table -> numOfElements);
	print_in_box(80,0,0,"");

	HashTableDestroy(merged_junction_table);
	HashTableDestroy(merged_splicing_table);
}

char * get_short_fname(char * lname){
	char * ret = lname;

	int x1;
	for(x1 = strlen(lname)-1; x1>=0; x1--){
		if(lname [x1] == '/'){
			ret = lname + x1 + 1;
			break;
		}
	}
	return ret;
}

int readSummary_single_file(fc_thread_global_context_t * global_context, read_count_type_t * column_numbers, int nexons,  int * geneid, char ** chr, long * start, long * stop, unsigned char * sorted_strand, char * anno_chr_2ch, char ** anno_chrs, long * anno_chr_head, long * block_end_index, long * block_min_start , long * block_max_end, fc_read_counters * my_read_counter, HashTable * junc_glob_tab, HashTable * splicing_glob_tab);

int readSummary(int argc,char *argv[]){

	/*
	   This function counts the number of reads falling into each exon region.
	   The order of exons in the output is the same as that of exons included in the annotation.
	   The annotation, if provided as a file, should be sorted by chromosome name.

	   Parameters passed from the featureCounts R function:
	0: "readSummary"
	1: ann
	2: files[i]
	3: fout
	4: as.numeric(isPairedEnd)
	5: min.distance
	6: max.distance
	7: as.numeric(tolower(file.type)=="sam")
	8: as.numeric(allowMultiOverlap)
	9: as.numeric(isGeneLevel)
	10: as.numeric(nthreads)
	11: as.numeric(isGTFannotation)
	12: as.numeric(isStrandChecked)
	13: as.numeric(isReadSummaryReported)
	14: as.numeric(isBothEndMapped)
	15: as.numeric(isChimericDisallowed)
	16: as.numeric(isPEDistChecked)
	17: nameFeatureTypeColumn 
	18: nameGeneIDColumn
	19: min.MappingQualityScore
	20: as.numeric(isMultiMappingAllowed)
	21: Annotation Chromosome Alias Name File. If the file is not specified, set this value to NULL or a zero-length string.
	22: Command line for CfeatureCounts header output; RfeatureCounts should set this value to NULL or a zero-length string or a space (' ').
	23: as.numeric(isInputFileResortNeeded)
	24: NOT IN USE: as.numeric(feature_block_size) # This parameter is no longer used. Give "14" for safe. 
	25: as.numeric(Five_End_Extension_Length)  # 5' end extension
	26: as.numeric(Three_End_Extension_Length)  # 3' end extension
	27: as.numeric(Minimum_Overlap_Between_Read_And_Feature) # 1 by default
	28: as.numeric(is_Split_or_Exonic_Only) # 0 by default; 0: all reads are counted ; 1: only split (Cigar has "N") reads are counted ; 2: only exonic (no "N" in Cigar) are counted.
	29: as.numeric(reduce_5_3_ends_to_one) # 0= no reduction; 1= reduce to 5' end; 2= reduce to 3' end
	30: debug_command # This is for debug only; RfeatureCounts should pass a space (" ") to this parameter, disabling the debug command.
	31: as.numeric(is_duplicate_ignored) # 0 = INCLUDE DUPLICATE READS; 1 = IGNORE DUPLICATE READS (0x400 FLAG IS SET) ; "0" by default.
	32: as.numeric(do_not_sort)   # 1 = NEVER SORT THE PE BAM/SAM FILES; 0 = SORT THE BAM/SAM FILE IF IT IS FOUND NOT SORTED.
	33: as.numeric(fractionMultiMapping) # 1 = calculate fraction numbers if a read overlaps with multiple features or meta-features. "-M" must be specified when fractions are caculated.
	34: as.numeric(useOverlappingBreakTie) # 1 = Select features or meta-features with a longer overlapping length; 0 = just use read-voting strategy: one overlapping read = 1 vote
	35: Pair_Orientations # FF, FR, RF or RR. This parameter matters only if "-s" option is 1 or 2.
	36: as.numeric(doJunctionCounting)  # 1 = count the number of junction reads spaining each exon-exon pairs;  0 = do not.
	37: file name of genome fasta (for determine the strandness of junctions by looking for GT/AG or CT/AC).
	38: as.numeric(max_M_Ops) # maximum "M" sections allowed in the CIGAR string. This parameter is needed in parse_BIN()
	39: as.numeric(is_Restrictly_No_Overlapping) # when "1", disable the voting-based tie breaking (e.g., when the reads are paired-end and one gene receives two votes but the other gene only has one.). "0" by default.
	40: as.numeric(min_Fractional_Overlap) # A fractioal number.  0.00 : at least 1 bp overlapping
	41: temp_directory # the directory to put temp files. "<use output directory>" by default, namely find it from the output file dir.
	 */

	int isStrandChecked, isCVersion, isChimericDisallowed, isPEDistChecked, minMappingQualityScore=0, isInputFileResortNeeded, feature_block_size = 20, reduce_5_3_ends_to_one;
	float fracOverlap;
	char **chr;
	long *start, *stop;
	int *geneid;

	char *nameFeatureTypeColumn, *nameGeneIDColumn,*debug_command, *pair_orientations="fr", *temp_dir;
	long nexons;


	long * anno_chr_head, * block_min_start, *block_max_end, *block_end_index;
	char ** anno_chrs, * anno_chr_2ch;
	long curchr, curpos;
	char * curchr_name, * fasta_contigs_fname;
	unsigned char * sorted_strand;
	curchr = 0;
	curpos = 0;
	curchr_name = "";

	int isPE, minPEDistance, maxPEDistance, isReadSummaryReport, isBothEndRequired, isMultiMappingAllowed, fiveEndExtension, threeEndExtension, minFragmentOverlap, isSplitOrExonicOnly, is_duplicate_ignored, doNotSort, fractionMultiMapping, useOverlappingBreakTie, doJuncCounting, max_M, isRestrictlyNoOvelrapping;

	int  isGTF, n_input_files=0;
	char *  alias_file_name = NULL, * cmd_rebuilt = NULL;

	int isMultiOverlapAllowed, isGeneLevel;

	isCVersion = ((argv[0][0]=='C')?1:0);

	isPE = atoi(argv[4]);
	minPEDistance = atoi(argv[5]);
	maxPEDistance = atoi(argv[6]);

	//  isSAM = atoi(argv[7]);
	isMultiOverlapAllowed = atoi(argv[8]);
	isGeneLevel = atoi(argv[9]);
	unsigned short thread_number;
	if(argc > 10)
		thread_number = atoi(argv[10]);
	else	thread_number = 4;
	if(argc > 11)
		isGTF = atoi(argv[11]);
	else	isGTF = 0;
	if(argc > 12)
		isStrandChecked = atoi(argv[12]);
	else	isStrandChecked = 0;
	if(argc > 13)
		isReadSummaryReport = atoi(argv[13]);
	else	isReadSummaryReport = 0;
	if(argc > 14)
		isBothEndRequired = atoi(argv[14]);
	else	isBothEndRequired = 0;
	if(argc > 15)
		isChimericDisallowed = atoi(argv[15]);
	else	isChimericDisallowed = 0;
	if(argc > 16)
		isPEDistChecked = atoi(argv[16]);
	else	isPEDistChecked = 0;
	if(argc > 17)
		nameFeatureTypeColumn = argv[17];
	else	nameFeatureTypeColumn = "exon";
	if(argc > 18)
		nameGeneIDColumn = argv[18];
	else	nameGeneIDColumn = "gene_id";
	if(argc > 19)
		minMappingQualityScore = atoi(argv[19]);
	else	minMappingQualityScore = 0;
	if(argc > 20)
		isMultiMappingAllowed = atoi(argv[20]);
	else	isMultiMappingAllowed = 1;
	if(argc > 21)
	{
		alias_file_name = argv[21];
		if(alias_file_name == NULL || alias_file_name[0]==' ' || alias_file_name[0]==0)
			alias_file_name = NULL;
	}
	else	alias_file_name = NULL;
	if(argc > 22)
	{
		cmd_rebuilt = argv[22];
		if(cmd_rebuilt == NULL || cmd_rebuilt[0]==' '||cmd_rebuilt[0]==0)
			cmd_rebuilt=NULL;
	}
	else	cmd_rebuilt = NULL;
	if(argc>23)
		isInputFileResortNeeded = atoi(argv[23]);
	else	isInputFileResortNeeded = 0;
	if(thread_number<1) thread_number=1;
	if(thread_number>16)thread_number=16;

	int Param_fiveEndExtension, Param_threeEndExtension;
	if(argc>25)
		Param_fiveEndExtension = atoi(argv[25]);
	else    Param_fiveEndExtension = 0;

	if(argc>26)
		Param_threeEndExtension = atoi(argv[26]);
	else    Param_threeEndExtension = 0;

	if(argc>27)
		minFragmentOverlap = atoi(argv[27]);
	else    minFragmentOverlap = 1;

	if(minFragmentOverlap <1){
		fiveEndExtension = 1 - minFragmentOverlap;
		threeEndExtension = 1 - minFragmentOverlap;
		minFragmentOverlap = 1;
	}else{
		fiveEndExtension = Param_fiveEndExtension;
		threeEndExtension = Param_threeEndExtension;
	}

	if(argc>28)
		isSplitOrExonicOnly = atoi(argv[28]);
	else	isSplitOrExonicOnly = 0;

	if(argc>29)
		reduce_5_3_ends_to_one = atoi(argv[29]);	// 0 : no reduce; 1: reduce to 5' end; 2: reduce to 3' end.
	else	reduce_5_3_ends_to_one = 0;


	if(argc>30 && strlen(argv[30])>0 && argv[30][0]!=' ')
		debug_command = argv[30];
	else
		debug_command = " ";

	if(argc>31)
		is_duplicate_ignored = atoi(argv[31]);
	else
		is_duplicate_ignored = 0;

	if(argc>32)
		doNotSort = atoi(argv[32]);
	else
		doNotSort = 0;

	if(argc>33)
		fractionMultiMapping = atoi(argv[33]);
	else
		fractionMultiMapping = 0;

	if(argc>34)
		useOverlappingBreakTie = atoi(argv[34]);
	else	useOverlappingBreakTie = 0;


	/*if(argc>35) "-S" is depreciated.
		pair_orientations = argv[35];
	else	pair_orientations = "FR";
	*/

	if(argc>36)
		doJuncCounting = atoi(argv[36]);
	else	doJuncCounting = 0;

	fasta_contigs_fname = NULL;
	if(argc>37)
		if(argv[37][0] != 0 && argv[37][0]!=' ')
			fasta_contigs_fname = argv[37];

	if(argc>38)
		max_M = atoi(argv[38]);
	else	max_M = 10;

	if(argc>39)
		isRestrictlyNoOvelrapping = atoi(argv[39]);
	else	isRestrictlyNoOvelrapping = 0;

	if(argc>40)
		fracOverlap = atof(argv[40]);
	else	fracOverlap= 0.0;
	
	if(argc>41){
		if(strcmp("<use output directory>", argv[41])!=0)temp_dir = argv[41];
		else temp_dir = NULL;
	}
	else	temp_dir = NULL;//	get_temp_dir_from_out(temp_dir, (char *)argv[3]);


	if(SAM_pairer_warning_file_open_limit()) return -1;

	fc_thread_global_context_t global_context;

	fc_thread_init_global_context(& global_context, FEATURECOUNTS_BUFFER_SIZE, thread_number, MAX_LINE_LENGTH, isPE, minPEDistance, maxPEDistance,isGeneLevel, isMultiOverlapAllowed, isStrandChecked, (char *)argv[3] , isReadSummaryReport, isBothEndRequired, isChimericDisallowed, isPEDistChecked, nameFeatureTypeColumn, nameGeneIDColumn, minMappingQualityScore,isMultiMappingAllowed, 0, alias_file_name, cmd_rebuilt, isInputFileResortNeeded, feature_block_size, isCVersion, fiveEndExtension, threeEndExtension , minFragmentOverlap, isSplitOrExonicOnly, reduce_5_3_ends_to_one, debug_command, is_duplicate_ignored, doNotSort, fractionMultiMapping, useOverlappingBreakTie, pair_orientations, doJuncCounting, max_M, isRestrictlyNoOvelrapping, fracOverlap, temp_dir);


	if( global_context.is_multi_mapping_allowed != ALLOW_ALL_MULTI_MAPPING && (!isMultiOverlapAllowed) && global_context.use_fraction_multi_mapping)
	{
		SUBREADprintf("ERROR: '--fraction' option should be used together with '-M' or '-O'. Please change the parameters to allow multi-mapping reads and/or multi-overlapping features.\n");
		return -1;
	}
	if( print_FC_configuration(&global_context, argv[1], argv[2], argv[3], global_context.is_SAM_file, isGTF, & n_input_files, isReadSummaryReport) )
		return -1;


	//print_in_box(80,0,0,"IG=%d, IS=%d", isGeneLevel, isSplitOrExonicOnly);
	if(0)if(isSplitOrExonicOnly && ( isGeneLevel || !isMultiOverlapAllowed) )
	{
		print_in_box(80,0,0,"NOTICE --splitOnly is specified, but '-O' and '-f' are not");
		print_in_box(80,0,0,"       both specified. Please read the manual for details.");
		print_in_box(80,0,0,"");
	}

	// Loading the annotations.
	// Nothing is done if the annotation does not exist.
	fc_feature_info_t * loaded_features;
	print_in_box(84,0,0,"Load annotation file %s %c[0m...", argv[1], CHAR_ESC);
	nexons = load_feature_info(&global_context,argv[1], isGTF?FILE_TYPE_GTF:FILE_TYPE_RSUBREAD, &loaded_features);
	if(nexons<1){
		if(nexons >= -1) SUBREADprintf("Failed to open the annotation file %s, or its format is incorrect, or it contains no '%s' features.\n",argv[1], nameFeatureTypeColumn);
		return -1;
	}

	sort_feature_info(&global_context, nexons, loaded_features, &chr, &geneid, &start, &stop, &sorted_strand, &anno_chr_2ch, &anno_chrs, &anno_chr_head, & block_end_index, & block_min_start, & block_max_end);
	if(global_context.do_junction_counting){
		sort_bucket_table(&global_context);
	}
	print_in_box(80,0,0,"   Meta-features : %d", global_context . gene_name_table -> numOfElements);
	print_in_box(80,0,0,"   Chromosomes/contigs : %d", global_context . exontable_nchrs);

	print_in_box(80,0,0,"");


	if(fasta_contigs_fname){
		print_in_box(80,0,0,"Loading FASTA contigs : %s", fasta_contigs_fname);
		global_context.fasta_contigs = malloc(sizeof(fasta_contigs_t));
		int ret_fq = read_contig_fasta(global_context.fasta_contigs, fasta_contigs_fname);
		if(ret_fq){
			print_in_box(80,0,0,"   WARNING unable to open the FASTA file.");
			print_in_box(80,0,0,"");
			free(global_context.fasta_contigs);
			global_context.fasta_contigs = NULL;
		}else{
			print_in_box(80,0,0,"   %lu contigs were loaded", global_context.fasta_contigs -> contig_table -> numOfElements);
			print_in_box(80,0,0,"");
		}
	}else	global_context.fasta_contigs = NULL;
	

	global_context.exontable_exons = nexons;
	unsigned int x1, * nreads = (unsigned int *) calloc(nexons,sizeof(int));




	char * tmp_pntr = NULL;
	char * file_list_used = malloc(strlen(argv[2])+1);
	char * file_list_used2 = malloc(strlen(argv[2])+1);
	char * is_unique = malloc(strlen(argv[2])+1);
	strcpy(file_list_used, argv[2]);
	for(x1 = 0;;x1++){
		char * test_fn = strtok_r(x1?NULL:file_list_used,";", &tmp_pntr);
		if(NULL == test_fn) break; 
		char * short_fname = get_short_fname(test_fn);
		strcpy(file_list_used2, argv[2]);

		is_unique[x1]=1;
		char * loop_ptr = NULL;
		int x2;
		for(x2 = 0;;x2++){
			char * test_loopfn = strtok_r(x2?NULL:file_list_used2, ";", &loop_ptr);
			if(NULL == test_loopfn) break;
			if(x1==x2)continue;

			char * short_loop_fname = get_short_fname(test_loopfn);

			if(strcmp(short_loop_fname, short_fname)==0) {
				is_unique[x1] = 0;
				break;
			}
		}
	}
	free(file_list_used2);

	tmp_pntr = NULL;
	strcpy(file_list_used, argv[2]);
	char * next_fn = strtok_r(file_list_used,";", &tmp_pntr);
	read_count_type_t ** table_columns = calloc( n_input_files , sizeof(read_count_type_t *)), i_files=0;
	fc_read_counters * read_counters = calloc(n_input_files , sizeof(fc_read_counters)); 
	HashTable ** junction_global_table_list = NULL;
	HashTable ** splicing_global_table_list = NULL;

	if(global_context.do_junction_counting){
		junction_global_table_list = calloc(n_input_files, sizeof(HashTable *));
		splicing_global_table_list = calloc(n_input_files, sizeof(HashTable *));
	}

	for(x1 = 0;;x1++){
		int orininal_isPE = global_context.is_paired_end_mode_assign;
		if(next_fn==NULL || strlen(next_fn)<1) break;

		read_count_type_t * column_numbers = calloc(nexons, sizeof(read_count_type_t));
		HashTable * junction_global_table = NULL;
		HashTable * splicing_global_table = NULL;

		strcpy(global_context.input_file_name, next_fn);
		strcpy(global_context.raw_input_file_name, next_fn);
		global_context.input_file_unique = is_unique[x1];
		global_context.input_file_short_name = get_short_fname(next_fn);
		//SUBREADprintf("UNQQ=%d ; FNAME=%s\n", global_context.input_file_unique, global_context.input_file_name);
		global_context.redo=0;
		

		if(global_context.do_junction_counting){
			junction_global_table = HashTableCreate(156679);
			splicing_global_table = HashTableCreate(156679);

			HashTableSetHashFunction(junction_global_table,HashTableStringHashFunction);
			HashTableSetDeallocationFunctions(junction_global_table, free, NULL);
			HashTableSetKeyComparisonFunction(junction_global_table, fc_strcmp_chro);

			HashTableSetHashFunction(splicing_global_table,HashTableStringHashFunction);
			HashTableSetDeallocationFunctions(splicing_global_table, free, NULL);
			HashTableSetKeyComparisonFunction(splicing_global_table, fc_strcmp_chro);
		}

		fc_read_counters * my_read_counter = &(read_counters[i_files]);
		memset(my_read_counter, 0, sizeof(fc_read_counters));

		int ret_int = readSummary_single_file(& global_context, column_numbers, nexons, geneid, chr, start, stop, sorted_strand, anno_chr_2ch, anno_chrs, anno_chr_head, block_end_index, block_min_start, block_max_end, my_read_counter, junction_global_table, splicing_global_table);
		if(ret_int!=0){
			// give up this file.

			table_columns[i_files] = NULL;
			if(global_context.do_junction_counting){
				HashTableDestroy(junction_global_table);
				HashTableDestroy(splicing_global_table);
			}
			free(column_numbers);
		} else {
			// finished
			table_columns[i_files] = column_numbers;
			if(global_context.do_junction_counting){
				junction_global_table_list[ i_files ] = junction_global_table;
				splicing_global_table_list[ i_files ] = splicing_global_table;
			}
		}
		global_context.is_paired_end_mode_assign = orininal_isPE;

		i_files++;
		next_fn = strtok_r(NULL, ";", &tmp_pntr);
	}

	free(file_list_used);

	if(global_context.is_input_bad_format){
		SUBREADprintf("\nFATAL Error: an input file has wrong format! The program has to terminate and no counting file is generated.\n\n");
	}else{
		if(isGeneLevel)
			fc_write_final_gene_results(&global_context, geneid, chr, start, stop, sorted_strand, argv[3], nexons,  table_columns, argv[2], n_input_files , loaded_features, isCVersion);
		else
			fc_write_final_results(&global_context, argv[3], nexons, table_columns, argv[2], n_input_files ,loaded_features, isCVersion);
	}
	if(global_context.do_junction_counting)
		fc_write_final_junctions(&global_context, argv[3], table_columns, argv[2], n_input_files , junction_global_table_list, splicing_global_table_list);

	fc_write_final_counts(&global_context, argv[3], n_input_files, argv[2], table_columns, read_counters, isCVersion);

	int total_written_coulmns = 0;
	for(i_files=0; i_files<n_input_files; i_files++)
		if(table_columns[i_files]){
			free(table_columns[i_files]);
			if(global_context.do_junction_counting){
				HashTableDestroy(junction_global_table_list[i_files]);
				HashTableDestroy(splicing_global_table_list[i_files]);
			}

			total_written_coulmns++;

		}
	free(table_columns);


	if(global_context.is_input_bad_format == 0) print_FC_results(&global_context);
	KeyValuePair * cursor;
	int bucket;
	for(bucket=0; bucket < global_context.exontable_chro_table  -> numOfBuckets; bucket++)
	{
		cursor = global_context.exontable_chro_table -> bucketArray[bucket];
		while (1)
		{
			if (!cursor) break;
			fc_chromosome_index_info * del_chro_info = cursor->value;
			free(del_chro_info->reverse_table_start_index);
			//free(del_chro_info->reverse_table_end_index);
			free((void *)cursor -> key);
			free(del_chro_info);
			cursor = cursor->next;
		}
	}

	if(global_context.SAM_output_fp) fclose(global_context. SAM_output_fp);
	HashTableDestroy(global_context.gene_name_table);
	free(global_context.gene_name_array);

	HashTableDestroy(global_context.exontable_chro_table);
	if(global_context.fasta_contigs){
		destroy_contig_fasta(global_context.fasta_contigs);
		free(global_context.fasta_contigs);
	}
	if(global_context.annot_chro_name_alias_table)
		HashTableDestroy(global_context.annot_chro_name_alias_table);
	if(global_context.do_junction_counting){
		HashTableDestroy(global_context.junction_bucket_table);
		HashTableDestroy(global_context.junction_features_table);
		free(junction_global_table_list);
		free(splicing_global_table_list);
	}

	free(global_context.unistr_buffer_space);
	free(loaded_features);
	free(geneid);
	free(chr);
	free(start);
	free(sorted_strand);
	free(anno_chr_2ch);
	free(anno_chrs);
	free(anno_chr_head);
	free(block_min_start);
	free(block_max_end);
	free(block_end_index);
	free(stop);
	free(nreads);


	return total_written_coulmns?0:-1;
}

void register_buckets(fc_thread_global_context_t * global_context , HashTable * gene_feature_table, char * chro_name){
	KeyValuePair * cursor;
	int bucket;
	for(bucket=0; bucket < gene_feature_table -> numOfBuckets; bucket++){
		cursor = gene_feature_table -> bucketArray[bucket];
		while(1){
			if (!cursor) break;
			fc_junction_gene_t * gene = (fc_junction_gene_t *) cursor -> value;
			unsigned int x1;
			
			for(x1 = gene -> pos_first_base - gene -> pos_first_base % JUNCTION_BUCKET_STEP; x1 <= gene -> pos_last_base ; x1 += JUNCTION_BUCKET_STEP){
				char bucket_key[CHROMOSOME_NAME_LENGTH + 20];
				sprintf(bucket_key, "%s:%u", chro_name, x1);
				gene_info_list_t * list = HashTableGet(global_context -> junction_bucket_table, bucket_key);
				if(list == NULL){
					list = malloc(sizeof(gene_info_list_t));
					list -> space = 3;
					list -> used = 0;
					list -> genes = malloc(sizeof(void *) * list -> space);
					char * mem_bucket_key = malloc(strlen(bucket_key) + 1);
					strcpy(mem_bucket_key , bucket_key);
					HashTablePut(global_context -> junction_bucket_table, mem_bucket_key , list);
				}

				if(list -> used  ==  list -> space){
					list -> space = max(list -> space + 3, list -> space * 1.3);
					list -> genes = realloc(list -> genes , list -> space * sizeof(void *));
				}
				list -> genes[list -> used++] = gene;
			}
			cursor = cursor -> next;
		}
	}
}

void sort_bucket_table(fc_thread_global_context_t * global_context){
	KeyValuePair * cursor;
	int bucket;
	for(bucket=0; bucket < global_context -> junction_features_table -> numOfBuckets; bucket++){
		cursor = global_context -> junction_features_table -> bucketArray[bucket];
		while(1){
			if (!cursor) break;
			HashTable * gene_feature_table = cursor -> value;
			char * chro_name = (char *)cursor -> key;
			register_buckets(global_context , gene_feature_table, chro_name);
			cursor = cursor -> next;
		}
	}
}




int readSummary_single_file(fc_thread_global_context_t * global_context, read_count_type_t * column_numbers, int nexons,  int * geneid, char ** chr, long * start, long * stop, unsigned char * sorted_strand, char * anno_chr_2ch, char ** anno_chrs, long * anno_chr_head, long * block_end_index, long * block_min_start , long * block_max_end, fc_read_counters * my_read_counter, HashTable * junction_global_table, HashTable * splicing_global_table)
{
	FILE *fp_in = NULL;
	int read_length = 0;
	int is_first_read_PE=0;
	char * line = (char*)calloc(MAX_LINE_LENGTH, 1);
	char * file_str = "";

	if(strcmp( global_context->input_file_name,"STDIN")!=0)
	{
		int file_probe = is_certainly_bam_file(global_context->input_file_name, &is_first_read_PE, NULL);
		
		global_context -> is_paired_end_input_file = is_first_read_PE;
		// a Singel-end SAM/BAM file cannot be assigned as a PE SAM/BAM file;
		// but a PE SAM/BAM file may be assigned as a SE file if the user wishes to do so.
		if(is_first_read_PE==0)
				global_context -> is_paired_end_mode_assign = 0;

		global_context->is_SAM_file = 1;
		if(file_probe == 1) global_context->is_SAM_file = 0;

		global_context -> start_time = miltime();

		file_str = "SAM";
		if(file_probe == 1) file_str = "BAM" ;
		if(file_probe == -1) file_str = "Unknown";

		if(!global_context->redo)
		{
			print_in_box(80,0,0,"Process %s file %s...", file_str, global_context->input_file_name);
			if(is_first_read_PE)
				print_in_box(80,0,0,"   Paired-end reads are included.");
			else
				print_in_box(80,0,0,"   Single-end reads are included.");
		}
		
	}

	if(strcmp( global_context->input_file_name,"STDIN")!=0)
	{
		FILE * exist_fp = f_subr_open( global_context->input_file_name,"r");
		if(!exist_fp)
		{
			print_in_box(80,0,0,"Failed to open file %s",  global_context->input_file_name);
			print_in_box(80,0,0,"No counts were generated for this file.");
			print_in_box(80,0,0,"");
			return -1;
		}
		fclose(exist_fp);
	}

	/*
	if(strcmp(global_context->input_file_name,"STDIN")!=0)
		if(warning_file_type(global_context->input_file_name, global_context->is_SAM_file?FILE_TYPE_SAM:FILE_TYPE_BAM))
			global_context->is_unpaired_warning_shown=1;
	*/
	// Open the SAM/BAM file
	// Nothing is done if the file does not exist.

	#ifdef MAKE_STANDALONE
	if(strcmp("STDIN",global_context->input_file_name)==0)
		fp_in = stdin;
	else
		fp_in = f_subr_open(global_context->input_file_name,"r");
	#else
		fp_in = f_subr_open(global_context->input_file_name,"r");
	#endif


	// begin to load-in the data.
	if(!global_context->redo)
	{
		if( global_context->is_paired_end_mode_assign)
		{
			print_in_box(80,0,0,"   Assign fragments (read pairs) to features...");
//				print_in_box(80,0,0,"   Each fragment is counted no more than once.");
		}
		else
			print_in_box(80,0,0,"   Assign reads to features...");
	}

	fc_thread_start_threads(global_context, nexons, geneid, chr, start, stop, sorted_strand, anno_chr_2ch, anno_chrs, anno_chr_head, block_end_index, block_min_start , block_max_end, read_length);

	global_context->is_all_finished = 1;
	fc_thread_wait_threads(global_context);

	unsigned long long int nreads_mapped_to_exon = 0;
	fc_thread_merge_results(global_context, column_numbers , &nreads_mapped_to_exon, my_read_counter, junction_global_table, splicing_global_table);
	fc_thread_destroy_thread_context(global_context);

	//global_context .read_counters.assigned_reads = nreads_mapped_to_exon;

	#ifdef MAKE_STANDALONE
	if(strcmp("STDIN",global_context->input_file_name)!=0)
	#endif
		fclose(fp_in);

	if(global_context -> sambam_chro_table) free(global_context -> sambam_chro_table);
	global_context -> sambam_chro_table = NULL;

	free(line);
	if(global_context -> is_input_bad_format) return -1;
	return 0;
}


#ifdef MAKE_STANDALONE
int main(int argc, char ** argv)
#else
int feature_count_main(int argc, char ** argv)
#endif
{
	char * Rargv[42];
	char annot_name[300];
	char temp_dir[300];
	char * out_name = malloc(300);
	char * fasta_contigs_name = malloc(300);
	char * alias_file_name = malloc(300);
	int cmd_rebuilt_size = 200;
	char * cmd_rebuilt = malloc(cmd_rebuilt_size);
	char max_M_str[8];
	char nameFeatureTypeColumn[66];
	char nameGeneIDColumn[66];
	int min_qual_score = 0;
	int min_dist = 50;
	int max_dist = 600;
	char debug_command[10];
	char min_dist_str[11];
	char max_dist_str[11];
	char min_qual_score_str[11];
	char feature_block_size_str[11];
	char Strand_Sensitive_Str[11];
	char Pair_Orientations[3];
	char * very_long_file_names;
	int is_Input_Need_Reorder = 0;
	int is_PE = 0;
	int is_SAM = 1;
	int is_GeneLevel = 1;
	int is_Overlap = 0;
	int is_Both_End_Mapped = 0;
	int is_Restrictedly_No_Overlap = 0;
	int feature_block_size = 14;
	int Strand_Sensitive_Mode = 0;
	int is_ReadSummary_Report = 0;
	int is_Chimeric_Disallowed = 0;
	int is_PE_Dist_Checked = 0;
	int is_Multi_Mapping_Allowed = 0;
	int is_Split_or_Exonic_Only = 0;
	int is_duplicate_ignored = 0;
	int do_not_sort = 0;
	int do_junction_cnt = 0;
	int reduce_5_3_ends_to_one = 0;
	int use_fraction_multimapping = 0;
	int threads = 1;
	int isGTF = 1;
	int use_overlapping_length_break_tie = 0;
	char nthread_str[4];
	int option_index = 0;
	int c;
	int very_long_file_names_size = 200;
	int fiveEndExtension = 0, threeEndExtension = 0, minFragmentOverlap = 1;
	float fracOverlap = 0.0;
	char strFiveEndExtension[11], strThreeEndExtension[11], strMinFragmentOverlap[11], fracOverlapStr[20];
	very_long_file_names = malloc(very_long_file_names_size);
	very_long_file_names [0] = 0;
	fasta_contigs_name[0]=0;

	alias_file_name[0]=0;
	debug_command[0] = 0;

	strcpy(nameFeatureTypeColumn,"exon");
	strcpy(nameGeneIDColumn,"gene_id");
	strcpy(temp_dir, "<use output directory>");
	annot_name[0]=0;out_name[0]=0;
	

	cmd_rebuilt[0]=0;
	for(c = 0; c<argc;c++)
	{
		if(strlen(cmd_rebuilt) + 300 > cmd_rebuilt_size)
		{
			cmd_rebuilt_size*=2;
			cmd_rebuilt = realloc(cmd_rebuilt, cmd_rebuilt_size);
		}
		sprintf(cmd_rebuilt+strlen(cmd_rebuilt), "\"%s\" ", argv[c]);
	}

	optind=0;
	opterr=1;
	optopt=63;
	strcpy(max_M_str, "10");
	strcpy(Pair_Orientations,"fr");

	while ((c = getopt_long (argc, argv, "G:A:g:t:T:o:a:d:D:L:Q:pbF:fs:S:CBJPMORv?", long_options, &option_index)) != -1)
		switch(c)
		{
			case 'S':
				/*
				if(strlen(optarg)!=2 || (strcmp(optarg, "ff")!=0 && strcmp(optarg, "rf")!=0 && strcmp(optarg, "fr")!=0)){
					SUBREADprintf("The order parameter can only be ff, fr or rf.\n");
					print_usage();
					return -1;
				}
				Pair_Orientations[0]=(optarg[0]=='r'?'r':'f');
				Pair_Orientations[1]=(optarg[1]=='f'?'f':'r');
				Pair_Orientations[2]=0;
				*/
				SUBREADprintf("The \"-S\" option has been depreciated.\n");

				break;
			case 'G':
				strcpy(fasta_contigs_name , optarg);
				break;
			case 'J':
				do_junction_cnt = 1;
				break;
			case 'A':
				strcpy(alias_file_name, optarg);
				break;
			case 'M':
				if(0 == is_Multi_Mapping_Allowed)
					is_Multi_Mapping_Allowed = ALLOW_ALL_MULTI_MAPPING;
				break;
			case 'v':
				core_version_number("featureCounts");
				return 0;
			case 'Q':
				if(!is_valid_digit_range(optarg, "Q", 0 , 255))
					STANDALONE_exit(-1);

				min_qual_score = atoi(optarg);
				break;
			case 't':
				strcpy(nameFeatureTypeColumn, optarg);
				break;
			case 'g':
				while((*optarg) == ' ') optarg++;
				strcpy(nameGeneIDColumn, optarg);
				break;
			case 'T':
				if(!is_valid_digit_range(optarg, "T", 1, 64))
					STANDALONE_exit(-1);

				threads = atoi(optarg);
				break;
			case 'd':
				if(!is_valid_digit(optarg, "d"))
					STANDALONE_exit(-1);

				min_dist = atoi(optarg);
				break;
			case 'D':
				if(!is_valid_digit(optarg, "D"))
					STANDALONE_exit(-1);

				max_dist = atoi(optarg);
				break;
			case 'p':
				is_PE = 1;
				break;
			case 'b':
				SUBREADprintf("The '-b' option has been deprecated.\n FeatureCounts will automatically examine the file format.\n");
				is_SAM = 0;
				break;
			case 'C':
				is_Chimeric_Disallowed = 1;
				break;
			case 'P':
				is_PE_Dist_Checked = 1;
				break;
			case 'B':
				is_Both_End_Mapped = 1;
				break;
			case 'f':
				is_GeneLevel = 0;
				break;
			case 'F':
				isGTF = 1;
				if(strcmp("SAF", optarg)==0) isGTF=0;
				else if(strcmp("GTF", optarg)==0) isGTF=1;
				else SUBREADprintf("\nWarning: Unknown annotation format: %s. GTF format is used.\n\n", optarg); 
				break;
			case 'O':
				is_Overlap = 1;
				break;
			case 'R':
				is_ReadSummary_Report = 1;
				break;
			case 's':
				if(!is_valid_digit_range(optarg, "s", 0 , 2))
					STANDALONE_exit(-1);

				Strand_Sensitive_Mode = atoi(optarg);
					
				break;
//			case 'i':
//				term_strncpy(sam_name, optarg,299);
//				break;
			case 'o':
				term_strncpy(out_name, optarg,299);
				break;
			case 'a':
				term_strncpy(annot_name, optarg,299);
				break;
			case 'L':
				feature_block_size = atoi(optarg);
				break;
			case 0 :	// long options

				if(strcmp("primary", long_options[option_index].name)==0)
				{
					is_Multi_Mapping_Allowed = ALLOW_PRIMARY_MAPPING;
				}

				if(strcmp("readExtension5", long_options[option_index].name)==0)
				{
					if(!is_valid_digit(optarg, "readExtension5"))
						STANDALONE_exit(-1);
					fiveEndExtension = atoi(optarg);
					fiveEndExtension = max(0, fiveEndExtension);
				}

				if(strcmp("readExtension3", long_options[option_index].name)==0)
				{
					if(!is_valid_digit(optarg, "readExtension3"))
						STANDALONE_exit(-1);
					threeEndExtension = atoi(optarg);
					threeEndExtension = max(0, threeEndExtension);
				}

				if(strcmp("fracOverlap", long_options[option_index].name)==0)
				{
					if(!is_valid_float(optarg, "fracOverlap"))
						STANDALONE_exit(-1);
					fracOverlap = atof(optarg);
				}

				if(strcmp("minOverlap", long_options[option_index].name)==0)
				{
					if(!is_valid_digit(optarg, "minOverlap"))
						STANDALONE_exit(-1);
					minFragmentOverlap = atoi(optarg);
				}

				if(strcmp("debugCommand", long_options[option_index].name)==0)
				{
					strcpy(debug_command, optarg);
				}


				if(strcmp("ignoreDup", long_options[option_index].name)==0)
				{
					is_duplicate_ignored = 1 ;
				}

				if(strcmp("fraction", long_options[option_index].name)==0)
				{
					use_fraction_multimapping = 1;
				}
				if(strcmp("tmpDir", long_options[option_index].name)==0){
					strcpy(temp_dir, optarg);
				}
				if(strcmp("maxMOp", long_options[option_index].name)==0){
					if(!is_valid_digit_range(optarg, "maxMOp", 1 , 64))
						STANDALONE_exit(-1);
					strcpy(max_M_str, optarg);
				}
				if(strcmp("read2pos", long_options[option_index].name)==0)
				{
					if(optarg[0]=='3')
						reduce_5_3_ends_to_one = REDUCE_TO_3_PRIME_END;
					else if(optarg[0]=='5')
						reduce_5_3_ends_to_one = REDUCE_TO_5_PRIME_END;
						
				}				

				if(strcmp("largestOverlap", long_options[option_index].name)==0)
				{
					use_overlapping_length_break_tie = 1;
				}

				if(strcmp("donotsort", long_options[option_index].name)==0)
				{
					do_not_sort = 1;
				}

				if(strcmp("splitOnly", long_options[option_index].name)==0)
				{
					is_Split_or_Exonic_Only = 1;
				}

				if(strcmp("restrictedlyNoOverlap", long_options[option_index].name)==0)
				{
					is_Restrictedly_No_Overlap = 1;
				}
				if(strcmp("nonSplitOnly", long_options[option_index].name)==0)
				{
					is_Split_or_Exonic_Only = 2;
				}
				break;
			case '?':
			default :
				print_usage();
				return -1;
				break;
		}


	if(minFragmentOverlap<1)
	{
		fiveEndExtension = - minFragmentOverlap + 1;
		threeEndExtension =  - minFragmentOverlap + 1;
		minFragmentOverlap = 1;
	}

	if(out_name[0]==0 || annot_name[0]==0||argc == optind)
	{
		print_usage();
		return -1;
	}

	for(; optind < argc; optind++)
	{
		int curr_strlen = strlen(very_long_file_names);
		if( very_long_file_names_size - curr_strlen <300)
		{
			very_long_file_names_size *=2;
			//printf("CL=%d ; NS=%d\n", curr_strlen , very_long_file_names_size);
			very_long_file_names=realloc(very_long_file_names , very_long_file_names_size);
		}

		strcat(very_long_file_names, argv[optind]);
		strcat(very_long_file_names, ";");
	}

	very_long_file_names[strlen(very_long_file_names)-1]=0;

	sprintf(strFiveEndExtension, "%d", fiveEndExtension);
	sprintf(strThreeEndExtension, "%d", threeEndExtension);
	sprintf(strMinFragmentOverlap, "%d", minFragmentOverlap);
	sprintf(nthread_str,"%d", threads);
	sprintf(min_dist_str,"%d",min_dist);
	sprintf(max_dist_str,"%d",max_dist);
	sprintf(min_qual_score_str,"%d", min_qual_score);
	sprintf(feature_block_size_str,"%d", feature_block_size);
	sprintf(Strand_Sensitive_Str,"%d", Strand_Sensitive_Mode);
	sprintf(fracOverlapStr, "%g", fracOverlap);
	Rargv[0] = "CreadSummary";
	Rargv[1] = annot_name;
	Rargv[2] = very_long_file_names;
	Rargv[3] = out_name;
	Rargv[4] = is_PE?"1":"0";
	Rargv[5] = min_dist_str;
	Rargv[6] = max_dist_str;
	Rargv[7] = is_SAM?"1":"0";
	Rargv[8] = is_Overlap?"1":"0";
	Rargv[9] = is_GeneLevel?"1":"0";
	Rargv[10] = nthread_str;
	Rargv[11] = isGTF?"1":"0";
	Rargv[12] = Strand_Sensitive_Str;
	Rargv[13] = is_ReadSummary_Report?"1":"0";
	Rargv[14] = is_Both_End_Mapped?"1":"0";
	Rargv[15] = is_Chimeric_Disallowed?"1":"0";
	Rargv[16] = is_PE_Dist_Checked?"1":"0";
	Rargv[17] = nameFeatureTypeColumn;
	Rargv[18] = nameGeneIDColumn;
	Rargv[19] = min_qual_score_str;
	Rargv[20] = is_Multi_Mapping_Allowed == ALLOW_PRIMARY_MAPPING?"2":(is_Multi_Mapping_Allowed == ALLOW_ALL_MULTI_MAPPING?"1":"0");
	Rargv[21] = alias_file_name;
	Rargv[22] = cmd_rebuilt;
	Rargv[23] = is_Input_Need_Reorder?"1":"0";
	Rargv[24] = feature_block_size_str;
	Rargv[25] = strFiveEndExtension;
	Rargv[26] = strThreeEndExtension;
	Rargv[27] = strMinFragmentOverlap;
	Rargv[28] = is_Split_or_Exonic_Only == 1?"1":(is_Split_or_Exonic_Only ==  2 ? "2":"0");
	Rargv[29] = (reduce_5_3_ends_to_one == 0?"0":(reduce_5_3_ends_to_one==REDUCE_TO_3_PRIME_END?"3":"5"));
	Rargv[30] = debug_command;
	Rargv[31] = is_duplicate_ignored?"1":"0";
	Rargv[32] = do_not_sort?"1":"0";
	Rargv[33] = use_fraction_multimapping?"1":"0";
	Rargv[34] = use_overlapping_length_break_tie?"1":"0";
	Rargv[35] = Pair_Orientations;
	Rargv[36] = do_junction_cnt?"1":"0";
	Rargv[37] = fasta_contigs_name;
	Rargv[38] = max_M_str;
	Rargv[39] = is_Restrictedly_No_Overlap?"1":"0"; 
	Rargv[40] = fracOverlapStr;
	Rargv[41] = temp_dir;
	int retvalue = readSummary(42, Rargv);

	free(very_long_file_names);
	free(out_name);
	free(alias_file_name);
	free(fasta_contigs_name);
	free(cmd_rebuilt);

	return retvalue;

}


