#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>

#include "subread.h"
#include "input-files.h"
#include "core.h"


static struct option long_options[] =
{
	{"index",  required_argument, 0, 'i'},
	{"read",  required_argument, 0, 'r'},
	{"read2",  required_argument, 0, 'R'},
	{"output",  required_argument, 0, 'o'},
	{"subreads",  required_argument, 0, 'n'},
	{"minvotes",  required_argument, 0, 'm'},
	{"minvotes2",  required_argument, 0, 'p'},
	{"singleSAM",  required_argument, 0, '1'},
	{"pairedSAM",  required_argument, 0, '2'},
	{"threads",  required_argument, 0, 'T'},
	{"indel",  required_argument, 0, 'I'},
	{"phred",  required_argument, 0, 'P'},
	{"mindist",  required_argument, 0, 'd'},
	{"maxdist",  required_argument, 0, 'D'},
	{"order",  required_argument, 0, 'S'},
	{"trim5", required_argument, 0, '5'},
	{"trim3", required_argument, 0, '3'},
	{"color-convert",  no_argument, 0, 'b'},
	{"junctionIns", required_argument, 0, 0},
	{"multi",  required_argument, 0, 'B'},
	{"rg",  required_argument, 0, 0},
	{"rg-id",  required_argument, 0, 0},
	{"gzFASTQinput",  no_argument, 0, 0},
	{"unique",  no_argument, 0, 'u'},
	{"SAMoutput", no_argument, 0, 0},
	{"BAMinput", no_argument, 0, 0},
	{"SAMinput", no_argument, 0, 0},
	{"hamming",  no_argument, 0, 'H'},
	{"quality",  no_argument, 0, 'Q'},
	{"fast",  no_argument, 0, 0},
	{"DPMismatch",  required_argument, 0, 'X'},
	{"DPMatch",  required_argument, 0, 'Y'},
	{"DPGapOpen",  required_argument, 0, 'G'},
	{"DPGapExt",  required_argument, 0, 'E'},
	{"extendIndelDetection", no_argument, 0, 0},
	{"allJunctions",  no_argument, 0, 0},
	{"memoryMultiplex",  required_argument, 0, 0},
	{"ignoreUnmapped",  no_argument, 0, 0},
	{"extraColumns",  no_argument, 0, 0},
	{"disableBigMargin",  no_argument, 0, 0},
	{"relaxMismatchedBases",  no_argument, 0, 0},
	{"reportPairedMultiBest",  no_argument, 0, 0},
	{"maxMismatches",  required_argument, 0, 'M'},
	{"exonicSubreadFrac",  required_argument, 0, 0},
	{"SVdetection", no_argument, 0, 0},
	{"maxVoteSimples",  required_argument, 0, 0},
	{"maxRealignLocations",  required_argument, 0, 0},
	{"minVoteCutoff",  required_argument, 0, 0},
	{"minMappedFraction",  required_argument, 0, 0},
	{"complexIndels", no_argument, 0, 0},
	{0, 0, 0, 0}
};

void print_usage_core_subjunc()
{

	SUBREADprintf("\nVersion %s\n\n", SUBREAD_VERSION);
	SUBREADputs("Usage:");
	SUBREADputs("");
	SUBREADputs(" ./subjunc [options] -i <index_name> -r <input> -o <output>");
	SUBREADputs("");
	SUBREADputs("Required arguments:");
	SUBREADputs("");
	SUBREADputs("  -i <index>        Base name of the index.");
	SUBREADputs("");
	SUBREADputs("  -r <string>       Name of the input file. Input formats including gzipped");
	SUBREADputs("                    fastq, fastq, and fasta can be automatically detected. If");
	SUBREADputs("                    paired-end, this should give the name of file including");
	SUBREADputs("                    first reads.");
	SUBREADputs("");
	SUBREADputs("Optional arguments:");
	SUBREADputs("");
	SUBREADputs("  -o <string>       Name of the output file. By default, the output is in BAM");
	SUBREADputs("                    format.");
	SUBREADputs("");
	SUBREADputs("  -n <int>          Number of selected subreads, 14 by default.");
	SUBREADputs("");
	SUBREADputs("  -m <int>          Consensus threshold for reporting a hit (minimal number of");
	SUBREADputs("                    subreads that map in consensus) . If paired-end, this gives");
	SUBREADputs("                    the consensus threshold for the anchor read. 1 by default");
	SUBREADputs("");
	SUBREADputs("  -M <int>          Specify the maximum number of mis-matched bases allowed in");
	SUBREADputs("                    the alignment. 3 by default. Mis-matches found in soft-");
	SUBREADputs("                    clipped bases are not counted.");
	SUBREADputs("");
	SUBREADputs("  -T <int>          Number of CPU threads used, 1 by default.");
	SUBREADputs("");
	SUBREADputs("  -I <int>          Maximum length (in bp) of indels that can be detected. 5 by");
	SUBREADputs("                    default. The program can detect indels of up to 200bp long.");
	SUBREADputs("");
	SUBREADputs("  -B <int>          Maximal number of equally-best mapping locations to be");
	SUBREADputs("                    reported. 1 by default. Note that -u option takes precedence");
	SUBREADputs("                    over -B.");
	SUBREADputs("");
	SUBREADputs("  -P <3:6>          Format of Phred scores used in input files, '3' for phred+33");
	SUBREADputs("                    and '6' for phred+64. '3' by default.");
	SUBREADputs("");
	SUBREADputs("  -u                Report uniquely mapped reads only. Number of mis-matched");
	SUBREADputs("                    bases is used to break the tie.");
	SUBREADputs("");
	SUBREADputs("  -b                Convert color-space read bases to base-space read bases in");
	SUBREADputs("                    the mapping output. Note that read mapping is performed at");
	SUBREADputs("                    color-space.");
	SUBREADputs("");
	SUBREADputs("  --SAMinput        Input reads are in SAM format.");
	SUBREADputs("");
	SUBREADputs("  --BAMinput        Input reads are in BAM format.");
	SUBREADputs("");
	SUBREADputs("  --SAMoutput       Save mapping result in SAM format.");
	SUBREADputs("");
	SUBREADputs("  --trim5 <int>     Trim off <int> number of bases from 5' end of each read. 0");
	SUBREADputs("                    by default.");
	SUBREADputs("");
	SUBREADputs("  --trim3 <int>     Trim off <int> number of bases from 3' end of each read. 0");
	SUBREADputs("                    by default.");
	SUBREADputs("");
	SUBREADputs("  --rg-id <string>  Add read group ID to the output.");
	SUBREADputs("");
	SUBREADputs("  --rg <string>     Add <tag:value> to the read group (RG) header in the output.");
	SUBREADputs("");
	SUBREADputs("  --DPGapOpen <int> Penalty for gap opening in short indel detection. -1 by");
	SUBREADputs("                    default.");
	SUBREADputs("");
	SUBREADputs("  --DPGapExt <int>  Penalty for gap extension in short indel detection. 0 by");
	SUBREADputs("                    default.");
	SUBREADputs("");
	SUBREADputs("  --DPMismatch <int> Penalty for mismatches in short indel detection. 0 by");
	SUBREADputs("                    default.");
	SUBREADputs("");
	SUBREADputs("  --DPMatch <int>   Score for matched bases in short indel detection. 2 by");
	SUBREADputs("                    default.");
	SUBREADputs("");
	SUBREADputs("  --allJunctions    Detect exon-exon junctions (both canonical and non-canonical");
	SUBREADputs("                    junctions) and structural variants in RNA-seq data. Refer to");
	SUBREADputs("                    Users Guide for reporting of junctions and fusions.");
	SUBREADputs("");
	SUBREADputs("  --complexIndels   Detect multiple short indels that occur concurrently in a");
	SUBREADputs("                    small genomic region (these indels could be as close as 1bp");
	SUBREADputs("                    apart).");
	SUBREADputs("");
	SUBREADputs("  -v                Output version of the program.");
	SUBREADputs("");
	SUBREADputs("Optional arguments for paired-end reads:");
	SUBREADputs("");
	SUBREADputs("  -R <string>       Name of the file including second reads.");
	SUBREADputs("");
	SUBREADputs("  -p <int>          Consensus threshold for the non-anchor read (receiving less");
	SUBREADputs("                    votes than the anchor read from the same pair). 1 by");
	SUBREADputs("                    default.");
	SUBREADputs("");
	SUBREADputs("  -d <int>          Minimum fragment/insert length, 50bp by default.");
	SUBREADputs("");
	SUBREADputs("  -D <int>          Maximum fragment/insert length, 600bp by default.");
	SUBREADputs("");
	SUBREADputs("  -S <ff:fr:rf>     Orientation of first and second reads, 'fr' by default (");
	SUBREADputs("                    forward/reverse).");
	SUBREADputs("");
	SUBREADputs("Refer to Users Manual for detailed description to the arguments.");
	SUBREADputs("");


}

int parse_opts_subjunc(int argc , char ** argv, global_context_t * global_context)
{
	int c;
	int option_index = 0;	
	int is_64_bit_computer = sizeof(char *)>4;

	optind = 0;
	opterr = 1;
	optopt = 63;

	global_context->config.entry_program_name = CORE_PROGRAM_SUBJUNC;
	global_context->config.max_mismatch_exonic_reads = 3;
	global_context->config.max_mismatch_junction_reads = 3;
	global_context->config.ambiguous_mapping_tolerance = 39 - 20 ;
	global_context->config.extending_search_indels = 0;
	global_context->config.do_fusion_detection =0;
	global_context->config.use_dynamic_programming_indel = 1;

	global_context->config.do_breakpoint_detection = 1;
	global_context->config.total_subreads = 14;
	global_context->config.minimum_subread_for_first_read =1;
	global_context->config.minimum_subread_for_second_read = 1;
	global_context->config.minimum_exonic_subread_fraction = 0.3;
	global_context->config.high_quality_base_threshold = 990000;
	global_context->config.do_big_margin_filtering_for_junctions = 1;
	global_context->config.report_no_unpaired_reads = 0;
	global_context->config.experiment_type = CORE_EXPERIMENT_RNASEQ;

	//#warning " ========================= REMOVE ' + 1 ' FROM THE NEXT LINE !! =========================="
	global_context->config.limited_tree_scan = 0;
	global_context->config.use_hamming_distance_in_exon = 0;
	//global_context->config.big_margin_record_size = 24;

	if(argc<2)
	{
		print_usage_core_subjunc();
		return -1;
	}
	while ((c = getopt_long (argc, argv, "vxsJ1:2:S:L:AHd:D:n:m:p:P:R:r:i:l:o:G:Y:E:X:T:I:B:bQFcuUfM:3:5:9:?", long_options, &option_index)) != -1)
	{
		switch(c)
		{
			case 'v':
				core_version_number("Subjunc");
				return -1;
			case 'G':
				if(!is_valid_digit(optarg, "G"))
					STANDALONE_exit(-1);

				global_context->config.DP_penalty_create_gap = atoi(optarg);
				break;
			case 'Y':
				if(!is_valid_digit(optarg, "Y"))
					STANDALONE_exit(-1);

				global_context->config.DP_match_score = atoi(optarg);
				break;
			case 'E':
				if(!is_valid_digit(optarg, "E"))
					STANDALONE_exit(-1);

				global_context->config.DP_penalty_extend_gap = atoi(optarg);
				break;
			case 'X':
				if(!is_valid_digit(optarg, "X"))
					STANDALONE_exit(-1);

				global_context->config.DP_mismatch_penalty = atoi(optarg);
				break;
			case '3':
				if(!is_valid_digit(optarg, "3"))
					STANDALONE_exit(-1);

				global_context->config.read_trim_3 = atoi(optarg); 
				break;
			case '5':
				if(!is_valid_digit(optarg, "5"))
					STANDALONE_exit(-1);

				global_context->config.read_trim_5 = atoi(optarg); 
				break;

			case '1':
			case '2':
			case 'r':
				strncpy(global_context->config.first_read_file, optarg, MAX_FILE_NAME_LENGTH-1);
				break;
			case 'J':
				global_context->config.show_soft_cliping = 0;
				break;
			case 'Q':
				global_context->config.use_quality_score_break_ties = 1;
				break;
			case 'H':
				global_context->config.use_hamming_distance_break_ties = 1;
				break;
			case 's':
				global_context->config.downscale_mapping_quality = 1;
				break;
			case 'A':
				global_context->config.report_sam_file = 0;
				break;
			case 'S':
				global_context->config.is_first_read_reversed = optarg[0]=='r'?1:0;
				global_context->config.is_second_read_reversed = optarg[1]=='f'?0:1;
				break;
			case 'U':
				global_context->config.report_no_unpaired_reads = 1;
				break;
			case 'u':
				global_context->config.report_multi_mapping_reads = 0;
				global_context->config.use_hamming_distance_break_ties = 1;
				break;
			case 'b':
				global_context->config.convert_color_to_base = 1;
				break;
			case 'D':
				if(!is_valid_digit(optarg, "D"))
					STANDALONE_exit(-1);

				global_context->config.maximum_pair_distance = atoi(optarg);
				break;
			case 'd':
				if(!is_valid_digit(optarg, "d"))
					STANDALONE_exit(-1);

				global_context->config.minimum_pair_distance = atoi(optarg);
				if(global_context->config.minimum_pair_distance <0)
					global_context->config.restrected_read_order = 0;
				break;
			case 'n':
				if(!is_valid_digit(optarg, "n"))
					STANDALONE_exit(-1);

				global_context->config.total_subreads = atoi(optarg);
				//global_context->config.total_subreads = min(31,global_context->config.total_subreads);
				break;
			case 'm':
				if(!is_valid_digit(optarg, "m"))
					STANDALONE_exit(-1);

				global_context->config.minimum_subread_for_first_read = atoi(optarg);
				break;
			case 'T':
				if(!is_valid_digit_range(optarg, "T", 1, 64))
					STANDALONE_exit(-1);

				global_context->config.all_threads = atoi(optarg);
				if(global_context->config.all_threads <1) global_context->config.all_threads = 1;
				if(global_context->config.all_threads > MAX_THREADS) global_context->config.all_threads = MAX_THREADS;

				break;
			case 'M':
				if(!is_valid_digit(optarg, "M"))
					STANDALONE_exit(-1);

				global_context->config.max_mismatch_exonic_reads = atoi(optarg);
				global_context->config.max_mismatch_junction_reads = atoi(optarg);
				break;
			case 'R':
				global_context->input_reads.is_paired_end_reads = 1;
				strncpy(global_context->config.second_read_file, optarg, MAX_FILE_NAME_LENGTH-1);
				break;
			case 'i':
				strncpy(global_context->config.index_prefix, optarg, MAX_FILE_NAME_LENGTH-1);
				break;
			case 'o':
				strncpy(global_context->config.output_prefix, optarg, MAX_FILE_NAME_LENGTH-1);
				break;
			case 'I':
				if(!is_valid_digit(optarg, "I"))
					STANDALONE_exit(-1);

				global_context->config.max_indel_length = atoi(optarg);

				if(!is_64_bit_computer) global_context->config.max_indel_length = min(global_context->config.max_indel_length , 16); 
				if(global_context->config.max_indel_length <0)global_context->config.max_indel_length =0;
				if(global_context->config.max_indel_length > MAX_INSERTION_LENGTH)global_context->config.max_indel_length = MAX_INSERTION_LENGTH;
				if(global_context->config.max_indel_length > 16)
				{
					global_context->config.reassembly_subread_length = 12;
					global_context->config.reassembly_window_multiplex = 3;
					global_context->config.reassembly_start_read_number = 5;
					global_context->config.reassembly_tolerable_voting = 0;
					global_context->config.reassembly_window_alleles = 2;
					global_context->config.reassembly_key_length = 28;
					global_context->config.flanking_subread_indel_mismatch = 0;

					global_context->config.is_third_iteration_running = 1;
					//global_context->config.max_mismatch_exonic_reads = 2;
					//global_context->config.max_mismatch_junction_reads = 2;
					//global_context->config.total_subreads = 28;
					//global_context->config.minimum_subread_for_first_read = 3;
					//global_context->config.minimum_subread_for_second_read = 1;
					//global_context->config.do_big_margin_filtering_for_reads = 1;

					global_context->config.do_superlong_indel_detection = 0;
				}
				break;
			case 'P':
				if (optarg[0]=='3')
					global_context->config.phred_score_format = FASTQ_PHRED33;
				else
					global_context->config.phred_score_format = FASTQ_PHRED64;
				break;
			case 'p':
				if(!is_valid_digit(optarg, "p"))
					STANDALONE_exit(-1);

				global_context->config.minimum_subread_for_second_read = atoi(optarg);
				break;
			case 'F':
				global_context->config.is_second_iteration_running = 0;
				global_context->config.report_sam_file = 0;
				break;
			case 'B':
				if(!is_valid_digit(optarg, "B"))
					STANDALONE_exit(-1);

				global_context->config.multi_best_reads = atoi(optarg); 

				if(global_context->config.multi_best_reads<1)
					global_context->config.multi_best_reads=1;

				global_context->config.reported_multi_best_reads = global_context->config.multi_best_reads;

				global_context->config.max_vote_combinations = max(global_context->config.max_vote_combinations, global_context->config.reported_multi_best_reads + 1);
				global_context->config.max_vote_simples = max(global_context->config.max_vote_simples, global_context->config.reported_multi_best_reads + 1);

				break;
			case 'c':
				global_context->config.space_type = GENE_SPACE_COLOR; 
				break;
				
			case 0:
				if(strcmp("memoryMultiplex", long_options[option_index].name)==0) 
				{
					global_context->config.memory_use_multiplex = atof(optarg);
				}
				else if(strcmp("ignoreUnmapped", long_options[option_index].name)==0) 
				{
					global_context->config.ignore_unmapped_reads = 1;
				}
				else if(strcmp("rg-id", long_options[option_index].name)==0) 
				{
					strcpy(global_context->config.read_group_id, optarg);
				}
				else if(strcmp("rg", long_options[option_index].name)==0) 
				{
					strcat(global_context->config.read_group_txt, "\t");
					strcat(global_context->config.read_group_txt, optarg);
				}
				else if(strcmp("SAMoutput", long_options[option_index].name)==0) 
				{
					global_context->config.is_BAM_output = 0;
				}
				else if(strcmp("BAMinput", long_options[option_index].name)==0) 
				{
					global_context->config.is_BAM_input = 1;
					global_context->config.is_SAM_file_input = 1;
				}
				else if(strcmp("SAMinput", long_options[option_index].name)==0) 
				{
					global_context->config.is_BAM_input = 0;
					global_context->config.is_SAM_file_input = 1;
				}
				else if(strcmp("junctionIns", long_options[option_index].name)==0)
				{
					if(!is_valid_digit(optarg, "junctionIns"))
						STANDALONE_exit(-1);
					global_context->config.check_donor_at_junctions=0;
					global_context->config.limited_tree_scan = 0;
					global_context->config.max_insertion_at_junctions = atoi(optarg);
				}
				else if(strcmp("extraColumns", long_options[option_index].name)==0) 
				{
					global_context->config.SAM_extra_columns=1;
				}
				else if(strcmp("minMappedFraction", long_options[option_index].name)==0) 
				{
					if(!is_valid_digit(optarg, "minMappedFraction"))
						STANDALONE_exit(-1);
					global_context->config.min_mapped_fraction = atoi(optarg);
				}
				else if(strcmp("relaxMismatchedBases", long_options[option_index].name)==0) 
				{
					global_context->config.max_mismatch_junction_reads = 999;
					global_context->config.max_mismatch_exonic_reads = 999;
					global_context->config.min_mapped_fraction = 61;
				}
				else if(strcmp("reportPairedMultiBest", long_options[option_index].name)==0) 
				{
					global_context->config.report_multiple_best_in_pairs = 1;
				}
				else if(strcmp("exonicSubreadFrac", long_options[option_index].name)==0) 
				{
					if(atof(optarg)>0)
						global_context->config.minimum_exonic_subread_fraction = atof(optarg);
					else	SUBREADprintf("WARNING: unknown parameter: --exonicSubreadFrac '%s'\n", optarg);
				}
				else if(strcmp("fast", long_options[option_index].name)==0) 
				{
					global_context -> config.fast_run = 1;
				}
				else if(strcmp("SVdetection", long_options[option_index].name)==0) 
				{
					global_context -> config.do_structural_variance_detection = 1;
					global_context -> config.use_memory_buffer = 1;
					global_context -> config.reads_per_chunk = 300llu*1024*1024;
				}
				else if(strcmp("minDistanceBetweenVariants", long_options[option_index].name)==0)
				{
					int newdist = atoi(optarg);
					newdist = max(newdist, 1);
					newdist = min(newdist, MAX_READ_LENGTH);
					global_context->config.realignment_minimum_variant_distance = newdist;
				}
				else if(strcmp("disableBigMargin", long_options[option_index].name)==0) 
				{
					global_context->config.do_big_margin_filtering_for_junctions = 0;
					global_context->config.limited_tree_scan = 0;
				}
				else if(strcmp("maxVoteSimples", long_options[option_index].name)==0)
				{
					global_context->config.max_vote_simples = atoi(optarg);
				}
				else if(strcmp("maxRealignLocations", long_options[option_index].name)==0)
				{
					global_context->config.max_vote_combinations = atoi(optarg);
					global_context->config.multi_best_reads = atoi(optarg);
				}
				else if(strcmp("complexIndels", long_options[option_index].name)==0)
				{
					global_context->config.maximise_sensitivity_indel = 1;
					global_context->config.realignment_minimum_variant_distance = 1;
					global_context->config.max_indel_length = 16;
				}
				else if(strcmp("disableBigMargin", long_options[option_index].name)==0)
				{
					global_context->config.big_margin_record_size = 0;
				}
				else if(strcmp("extendIndelDetection", long_options[option_index].name)==0)
				{
					global_context->config.extending_search_indels = 1;
				}
				else if(strcmp("minVoteCutoff", long_options[option_index].name)==0)
				{
					global_context->config.max_vote_number_cutoff  = atoi(optarg);
				}
				else if(strcmp("allJunctions", long_options[option_index].name)==0)
				{
					global_context->config.do_fusion_detection = 1;
				}

				break;


			case '?':
			default:
				print_usage_core_subjunc();
				return -1 ;
		}
	}

	if(argc > optind){
		SUBREADprintf("Invalid parameter '%s'\n", argv[optind]);
		return -1;
	}


	if(global_context->config.is_SAM_file_input) global_context->config.phred_score_format = FASTQ_PHRED33;
	global_context->config.more_accurate_fusions = global_context->config.more_accurate_fusions && global_context->config.do_fusion_detection;

	if(global_context->config.more_accurate_fusions)
	{
		global_context->config.high_quality_base_threshold = 999999;
		//#warning "============ REMOVE THE NEXT LINE ======================"
		global_context->config.show_soft_cliping = 1;
		//#warning "============ REMOVE ' + 3' FROM NEXT LINE =============="
		global_context->config.max_mismatch_junction_reads = 0 + 3;

		//#warning "============ REMOVE ' - 1' FROM NEXT LINE =============="
		global_context->config.do_big_margin_filtering_for_junctions = 1 - 1;
		global_context->config.total_subreads = 28;

		//#warning "============ REMOVE THE NEXT LINE BEFORE RELEASE ==============="
		//global_context->config.multi_best_reads = 1;
	}



	return 0;
}





#if defined MAKE_STANDALONE
int main(int argc , char ** argv)
{
#elif defined RUNNING_ENV_JAVA
int subread_subjunc_main(int argc , char ** argv)
{
#else
int main_junction(int argc , char ** argv)
{
#endif
	return core_main(argc, argv, parse_opts_subjunc);
}

