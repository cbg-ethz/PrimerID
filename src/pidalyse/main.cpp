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
    std::string referenceFile_3223;
    std::string referenceFile_3236;

    bool removeCollisions = true;

    int minCoverage;
    double nonAmbigFrac;
    std::vector<std::string> inputFiles;

    // program options
    boost::program_options::options_description generic("Generic options");
    generic.add_options()("help,h", "Print this help");

    // configuration options
    boost::program_options::options_description config("Configuration");
    config.add_options()("r3223", boost::program_options::value<std::string>(&referenceFile_3223)->required(), "File containing the 3223 references (REQUIRED)")("r3236", boost::program_options::value<std::string>(&referenceFile_3236)->required(), "File containing the 3236 references (REQUIRED)")("leave-col", "Leave collisions (i.e., do not attempt to remove them)")("minCov", boost::program_options::value<int>(&minCoverage)->default_value(10), "Minimum number of reads per pID")(",f", boost::program_options::value<double>(&nonAmbigFrac)->default_value(0.8, "0.8"), "Minimum fraction of non-ambiguous bases A,C,G,T to call majority base");

    // hidden options, i.e., input files
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()("input-file", boost::program_options::value<std::vector<std::string>>()->required(), "input file");

    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(config).add(hidden);
    boost::program_options::options_description visible("Allowed options");
    visible.add(generic).add(config);
    boost::program_options::positional_options_description p;
    p.add("input-file", -1);
    boost::program_options::variables_map global_options;

    // parse command-line arguments
    try {
        boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), global_options);

        // show help options
        if (global_options.count("help")) {
            std::cout << visible << "\n";
            return EXIT_SUCCESS;
        }

        boost::program_options::notify(global_options);
    }
    catch (boost::program_options::required_option& e) {
        if (e.get_option_name() == "--input-file")
            std::cerr << "ERROR: no input files provided\n";
        else
            std::cerr << "ERROR: " << e.what() << "\n";
        return EXIT_FAILURE;
    }
    catch (boost::program_options::error& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    // 1.) parse command line options
    removeCollisions = !(global_options.count("leave-col"));
    inputFiles = global_options["input-file"].as<std::vector<std::string>>();

    // 2.) show parameters
    std::cout << "Parameters:\n";
    std::cout << "            Min coverage for calling consensus: " << minCoverage << '\n';
    std::cout << "  Minimum plurality required for majority base: " << nonAmbigFrac << '\n';
    std::cout << "                         Number of input files: " << inputFiles.size() << '\n';
    std::cout << (removeCollisions ? "Collisions should be removed" : "Collisions will not be removed") << "\n\n";

    // 3.) load reference
    std::cout << "1. Loading references\n";
    ::ref3223 = reference(referenceFile_3223, referenceType::START_3223);
    ::ref3236 = reference(referenceFile_3236, referenceType::START_3236);
    set_reference_intersection(::ref3223, ::ref3236);

    ref3223.display_strains_hetero();
    ref3236.display_strains_hetero();

    ref3223.display_hamming_distance();
    ref3236.display_hamming_distance();

    //ref3223.display_strains_verbose();
    //ref3236.display_strains_verbose();

    // 4.) load alignments
    std::cout << "2. Loading alignments\n";
    alignments DATA(inputFiles);

    // 5.) perform QA and filtering
    std::cout << "3. Filtering and performing QA\n";
    DATA.filtering_QA();

    // 6.) Remove collisions
    std::cout << "4. Removing collisions\n";
    DATA.remove_pID_collisions(minCoverage, nonAmbigFrac, removeCollisions);

    // 7.) Calculate estimators
    std::cout << "5. Calculating estimators\n";
    DATA.estimate_PCR_substitution_rate();

    DATA.estimate_RT_substitution_rate();

    DATA.estimate_RTPCR_recombination_rate();
    DATA.plot_RTPCR_recombination_LogLik(1E-6);

    // 8.) Perform remaining statistics
    std::cout << "6. Writing statistics to file\n";
    DATA.show_clone_frequencies();
    DATA.write_frequencies_to_MATLAB();
    DATA.write_all_statistics();

    return EXIT_SUCCESS;
}