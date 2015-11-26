#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <getopt.h>
static struct option long_options[] = {
    { 0, 0, 0, 0 }
};

template <typename Store>
std::vector<Store> readFasta(const std::string& fileName, std::function<bool(const std::string& header, const std::string& Seq)> acceptFunc = [](const std::string& header, const std::string& Seq) {return true; })
{
    // 1.) load FASTA file
    std::ifstream input;
    input.open(fileName.c_str());

    std::vector<Store> result;

    if (input.is_open()) {
        std::string temp, seq, id;

        while (input.good()) {
            getline(input, temp);
            temp.erase(std::remove(temp.begin(), temp.end(), '\n'), temp.end());

            if (temp.length() == 0)
                continue;

            if (temp[0] == '>') {
                // header
                if (!(seq.empty())) {
                    if (acceptFunc(id, seq))
                        result.emplace_back(id, std::move(seq));
                    seq.clear();
                }
                id = temp.substr(1, 1000);
            }
            else {
                // sequence
                seq.append(temp);
            }
        }

        if (acceptFunc(id, seq))
            result.emplace_back(id, std::move(seq));
    }
    else {
        std::cout << fileName << " is not readable!\n";
        exit(EXIT_FAILURE);
    }
    input.close();

    return result;
}

// >Clique_xaaaaaaaaaa-dir_14/1_0.150535

struct fasta_record {
    std::string _name;
    std::string _seq;

    fasta_record() = delete;
    fasta_record(std::string&& seq)
        : _seq(std::move(seq)) {}
    fasta_record(const std::string& name, std::string&& seq)
        : _name(name), _seq(std::move(seq)) {}
};

auto keep_HC_variant = [](const std::string& header, const std::string& Seq) {
	int i = header.find_last_of('_')+1;
	double freq = stod(header.substr(i, std::string::npos));
	
	return ((header.find("/1") != std::string::npos) && (freq > 0.01));
};

struct HC_variant : public fasta_record {
    double _freq;

    HC_variant() = delete;
    HC_variant(const std::string& name, std::string&& seq)
        : fasta_record(std::move(seq))
    {
        int i = name.find_last_of('_') + 1;
        _freq = stod(name.substr(i, std::string::npos));
        _name = name.substr(1, i - 3);
    }
};

std::string refFile;
std::string outputFile;
std::vector<std::string> inputFiles;

void parse_arguments(int argc, char* argv[])
{
    int c, option_index = 0;

    while ((c = getopt_long(argc, argv, "r:o:", long_options, &option_index)) != -1) {
        switch (c) {
        case 'r':
            refFile = optarg;
            break;

        case 'o':
            outputFile = optarg;
            break;

        case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
                break;

            printf("option %s", long_options[option_index].name);
            if (optarg)
                printf(" with arg %s", optarg);

            printf("\n");
            break;

        default:
            abort();
        }
    }

    if (optind < argc) {
        for (int i = optind; i < argc; ++i)
            inputFiles.push_back(argv[i]);
    }
    else {
        std::cout << "Missing input file!\n";
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char* argv[])
{
    parse_arguments(argc, argv);
    std::vector<fasta_record> refs = readFasta<fasta_record>(refFile);
    std::vector<std::vector<HC_variant>> HCs;
    HCs.resize(inputFiles.size());

    // normalize
    for (int i = 0; i < inputFiles.size(); ++i) {
        std::vector<HC_variant>& HC_dataset = HCs[i];
        HC_dataset = readFasta<HC_variant>(inputFiles[i], keep_HC_variant);

        double sum = 0;
        for (const auto& j : HC_dataset)
            sum += j._freq;
        for (auto& j : HC_dataset)
            j._freq /= sum;
    }

    // output
    bool writeFile = !(outputFile.empty());
    std::ofstream output;
    if (writeFile)
        output.open(outputFile.c_str(), std::ios_base::app);

    for (const auto& i : refs)
        std::cout << "    \t" << i._name;
    std::cout << '\n';

    for (int i = 0; i < HCs.size(); ++i) {
        std::cout << inputFiles[i].substr(inputFiles[i].find("32"), 5);
        bool firstHap = true;

        for (const auto& j : refs) {
            for (const auto& k : HCs[i]) {
                if (j._seq.find(k._seq) != std::string::npos) {
                    // stdout
                    std::cout << std::fixed << std::setprecision(6) << '\t' << k._freq;

                    // to file
                    if (writeFile)
                        output << std::fixed << std::setprecision(6) << (firstHap ? firstHap = false, "" : ",") << k._freq;
                }
            }
        }
        std::cout << '\n';
        if (writeFile)
            output << ";\n";
    }

    if (writeFile)
        output.close();

    return 0;
}