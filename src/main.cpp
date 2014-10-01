#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

#include <reference.hpp>
#include <alignment.hpp>

#include <getopt.h>
static struct option long_options[] =
{
	{"r3223", required_argument, 0, 1000},
	{"r3236", required_argument, 0, 1010},
	{0, 0, 0, 0}
};

std::string referenceFile_3223;
std::string referenceFile_3236;

int minCoverage = 10;
double nonAmbigFrac = 0.9;
std::vector<std::string> inputFiles;
bool help_flag = false;

void parse_arguments(int argc, char *argv[])
{
	int c, option_index = 0;
	
	while ((c = getopt_long(argc, argv, "m:d:h", long_options, &option_index)) != -1)
	{
		switch (c)
		{
			case 'm':
				minCoverage = atoi(optarg);
				break;
				
			case 'd':
				nonAmbigFrac = atof(optarg);
				break;
				
			case 'h':
				help_flag = true;
				break;
				
			case 1000:
				referenceFile_3223 = optarg;
				break;
			
			case 1010:
				referenceFile_3236 = optarg;
				break;
				
			case 0:
				/* If this option set a flag, do nothing else now. */
				if (long_options[option_index].flag != 0)
					break;
				printf ("option %s", long_options[option_index].name);
				if (optarg)
					printf (" with arg %s", optarg);
				printf ("\n");
				break;
				
			default:
				abort ();
		}
	}
	
	if ((help_flag) || (argc == 1))
	{
		std::cout << "Help!\n";
		exit(EXIT_SUCCESS);
	}
	
	if (optind < argc)
	{
		for (int i = optind; i < argc; ++i)
			inputFiles.push_back(argv[i]);
	}
	else
	{
		std::cout << "Missing input file!\n";
		exit(EXIT_FAILURE);
	}
	
	for(const std::string& i : inputFiles)
	{
		if ((i.find("3223") == std::string::npos) && (i.find("3236") == std::string::npos))
		{
			std::cout << "Input filename " << i << " does not contain either '3223' or '3236', which is required for reference assignment purposes\n";
			exit(EXIT_FAILURE);
		}
	}
	
	std::cout << "Minimum calling coverage:  " << minCoverage << '\n';
	std::cout << "Minimum majority fraction: " << nonAmbigFrac << '\n';
}

reference ref3223;
reference ref3236;

int main(int argc, char *argv[])
{
	parse_arguments(argc, argv);
	
	// 1.) load reference
	std::cout << "1. Loading references\n";
	ref3223 = reference(referenceFile_3223, START_3223);
	ref3236 = reference(referenceFile_3236, START_3236);
	
	ref3223.display_strains_verbose();
	ref3236.display_strains_verbose();
	
	//ref3223.display_strains_hetero();
	//ref3236.display_strains_hetero();
	
	// 2.) load alignments
	std::cout << "2. Loading alignments\n";
	alignments DATA(inputFiles);
	
	// 3.) process data
	std::cout << "3. Processing data\n";
	DATA.remove_primerID_collisions(minCoverage, nonAmbigFrac, true);
	
	/*
	// 4.) calculate statistics
	std::cout << "4. Performing statistics\n";
	DATA.estimate_RT_substitution_rate(true);
	DATA.estimate_RT_recombination_rate(true);
	DATA.estimate_PCR_substitution_rate(true);
	
	DATA.estimate_effective_RNA_number(true);
	
	DATA.display_raw_and_primerID_counts();
	*/
	
	return EXIT_SUCCESS;
}