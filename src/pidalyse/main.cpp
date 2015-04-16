#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

#include "reference.hpp"
#include "alignment.hpp"

#include <boost/program_options.hpp>

reference ref3223;
reference ref3236;

int main(int argc, char* argv[])
{
	// program options
	std::string referenceFile_3223;
	std::string referenceFile_3236;

	int minCoverage = 10;
	double nonAmbigFrac = 0.8;
	std::vector<std::string> inputFiles;
	bool help_flag = false;
	
  // generic options
  boost::program_options::options_description generic("Generic options");
  generic.add_options()("help,h", "Print this help");

  // configuration options
  boost::program_options::options_description config("Configuration");
  config.add_options()
		("r3223", boost::program_options::value<std::string>(&referenceFile_3223)->required(), "File containing the 3223 references (REQUIRED)")
		("r3236", boost::program_options::value<std::string>(&referenceFile_3236)->required(), "File containing the 3236 references (REQUIRED)")
		("minCov", boost::program_options::value<int>(&minCoverage)->default_value(10), "Minimum number of reads per pID")
		(",f", boost::program_options::value<double>(&nonAmbigFrac)->default_value(0.8, "0.8"), "Minimum fraction of non-ambiguous bases A,C,G,T to call majority base");

  // hidden options, i.e., input files
  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
		("input-file", boost::program_options::value<std::vector<std::string>>()->required(), "input file");

  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(config).add(hidden);

  boost::program_options::options_description visible("Allowed options");
  visible.add(generic).add(config);

  boost::program_options::positional_options_description p;
  p.add("input-file", -1);

  boost::program_options::variables_map global_options;

  try
  {
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), global_options);
		
		// show help options
	  if (global_options.count("help"))
	  {
	    std::cout << visible << "\n";
	    return EXIT_SUCCESS;
	  }
		
    boost::program_options::notify(global_options);
  }
  catch (boost::program_options::required_option& e)
  {
    if (e.get_option_name() == "--input-file")
			std::cerr << "ERROR: no input files provided\n";
		else
			std::cerr << "ERROR: " << e.what() << "\n";
    return EXIT_FAILURE;
  }
  catch (boost::program_options::error& e)
  {
    std::cerr << "ERROR: " << e.what() << "\n";
    return EXIT_FAILURE;
  }

	// load input files
  inputFiles = global_options["input-file"].as<std::vector<std::string>>();

  for (const auto& i : inputFiles)
  {
    std::cout << i << '\n';
  }

  exit(0);

  // 1.) load reference
  std::cout << "1. Loading references\n";
  ref3223 = reference(referenceFile_3223, START_3223);
  ref3236 = reference(referenceFile_3236, START_3236);

  // ref3223.display_strains_hetero();
  // ref3236.display_strains_hetero();

  // 2.) load alignments
  std::cout << "2. Loading alignments\n";
  alignments DATA(inputFiles);

  //DATA.write_raw_to_fasta();

  // 3.) process data
  std::cout << "3. Processing data\n";
  DATA.filtering_QA();

  DATA.remove_primerID_collisions(minCoverage, nonAmbigFrac, true);
  //DATA.show_clone_frequencies();

  // 4.) calculate statistics
  std::cout << "4. Performing statistics\n";

  DATA.estimate_RT_substitution_rate(true);
  DATA.estimate_RT_recombination_rate(true);
  DATA.estimate_PCR_substitution_rate(true);

  DATA.write_prob_to_csv();

  DATA.display_raw_and_primerID_counts();

  /*
	// DATA.calculate_RT_bias_pvalue();
  // DATA.estimate_effective_RNA_number(true);
  // DATA.show_recombination_patterns();
  // DATA.plot_RT_recombination_LogLik(1E-6, 2000);
	*/

  return EXIT_SUCCESS;
}