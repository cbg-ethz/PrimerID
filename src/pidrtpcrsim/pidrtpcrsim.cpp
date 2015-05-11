#include <iostream>
#include <fstream>

#include <vector>
#include <bitset>
#include <set>
#include <utility>
#include <numeric>
#include <cstdint>

#include <random>
#include <gsl/gsl_randist.h>
#include <boost/program_options.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/thread.hpp>
#include <boost/format.hpp>

typedef unsigned int intType;
typedef std::vector<intType> intVector;

intType bins;
intType N_RT;
intType cycles;
double p;
intType N_trials;
intType N_reads;
intType rank_cutoff;
intType min_coverage;
intType number_threads;
bool silent;

intType used_with_replacement = 0;
intType used_without_replacement = 0;
double average_PCR_molecules = 0;

std::vector<intType> ran_mult_multinomial(const std::vector<uint64_t>& pop, intType sample_size, const gsl_rng* r, intType min_coverage = 0)
{
    std::vector<intType> sample;
    sample.resize(bins);

    // 1.) first generate N_reads samples, not all of which will satisfy >= min_coverage
    std::vector<double> p_PCR(pop.begin(), pop.end());
    gsl_ran_multinomial(r, bins, sample_size, p_PCR.data(), sample.data());

    uint64_t valid_samples = 0;
    for (const auto& i : sample) {
        if (i >= min_coverage) {
            valid_samples += i;
        }
    }

    // 2.) proceed to add remaining samples
    if (valid_samples < sample_size) {
        gsl_ran_discrete_t* categorical_distrib = gsl_ran_discrete_preproc(bins, p_PCR.data());
        intType ind;

        while (valid_samples < sample_size) {
            ind = gsl_ran_discrete(r, categorical_distrib);
            ++sample[ind];

            if (sample[ind] > min_coverage) {
                ++valid_samples;
            }
            else {
                if (sample[ind] == min_coverage)
                    valid_samples += min_coverage;
            }
        }

        gsl_ran_discrete_free(categorical_distrib);
    }

    return sample;
}

/*
// bitset version for bookkeeping of sampled items
std::vector<intType> ran_mult_hypergeometric(const std::vector<uint64_t>& pop, intType sample_size, const gsl_rng* r, std::default_random_engine& generator, intType min_coverage = 0)
{
  std::vector<intType> sample(bins, 0);
  std::vector<uint64_t> part_sum;
	//part_sum.reserve(bins);
	part_sum.resize(bins);

  // 1.) pile-up counts for later determination of cluster
  std::partial_sum(pop.begin(), pop.end(), part_sum.begin());
  typename std::vector<uint64_t>::const_iterator part_sum_begin = part_sum.begin();
  typename std::vector<uint64_t>::const_iterator part_sum_end = part_sum.end();
  uint64_t pop_size = part_sum.back();

  boost::dynamic_bitset<> item_inventory(pop_size, 0);
	std::uniform_int_distribution<uint64_t> distribution(0, pop_size-1);

  // 2.) start sampling
  uint64_t valid_samples = 0;
  while (valid_samples < sample_size)
  {
    uint64_t rand_item = distribution(generator);

    if (!item_inventory[rand_item])
    {
      item_inventory[rand_item] = 1;

      typename std::vector<uint64_t>::const_iterator it = std::upper_bound(part_sum_begin, part_sum_end, rand_item);
      intType ind = std::distance(part_sum_begin, it);
      ++sample[ind];

      //std::cout << "Ind: " << rand_item << " - item: " << ind << "\n";

			if (sample[ind] > min_coverage)
			{
				++valid_samples;
			}
			else
			{
			  if (sample[ind] == min_coverage)
					valid_samples += min_coverage;
			}
    }
  }
	
  return sample;
}
*/

/*
// std::set version for bookkeeping of sampled items
std::vector<intType> ran_mult_hypergeometric(const std::vector<uint64_t>& pop, intType sample_size, const gsl_rng* r, std::default_random_engine& generator, intType min_coverage = 0)
{
  std::vector<intType> sample(bins, 0);
  std::vector<uint64_t> part_sum;
	//part_sum.reserve(bins);
	part_sum.resize(bins);

  // 1.) pile-up counts for later determination of cluster
  std::partial_sum(pop.begin(), pop.end(), part_sum.begin());
  std::vector<uint64_t>::const_iterator part_sum_begin = part_sum.begin();
  std::vector<uint64_t>::const_iterator part_sum_end   = part_sum.end();
  uint64_t pop_size = part_sum.back();

	std::set<uint64_t> item_inventory;
	std::uniform_int_distribution<uint64_t> distribution(0, pop_size-1);

  // 2.) start sampling
	std::vector<uint64_t>::const_iterator it;
  uint64_t valid_samples = 0;
  while (valid_samples < sample_size)
  {
    uint64_t rand_item = distribution(generator);
		std::pair<std::set<uint64_t>::iterator, bool> _it = item_inventory.insert(rand_item);

    if (_it.second)
    {
      it = std::upper_bound(part_sum_begin, part_sum_end, rand_item);
      intType ind = std::distance(part_sum_begin, it);
      ++sample[ind];

			if (sample[ind] > min_coverage)
			{
				++valid_samples;
			}
			else
			{
			  if (sample[ind] == min_coverage)
					valid_samples += min_coverage;
			}
    }
  }
	
  return sample;
}
*/

// sorted vector version
std::vector<intType> ran_mult_hypergeometric(const std::vector<uint64_t>& pop, intType sample_size, const gsl_rng* r, std::default_random_engine& generator, intType min_coverage = 0)
{
    std::vector<intType> sample(bins, 0);
    std::vector<uint64_t> part_sum;
    part_sum.resize(bins);

    // 1.) pile-up counts for later determination of cluster
    std::partial_sum(pop.begin(), pop.end(), part_sum.begin());
    std::vector<uint64_t>::const_iterator part_sum_begin = part_sum.begin();
    std::vector<uint64_t>::const_iterator part_sum_end = part_sum.end();
    uint64_t pop_size = part_sum.back();

    std::vector<uint64_t> item_inventory;
    item_inventory.resize(2 * sample_size);
    std::uniform_int_distribution<uint64_t> distribution(0, pop_size - 1);

    for (intType i = 0; i < 2 * sample_size; ++i) {
        item_inventory[i] = distribution(generator);
    }

    // 2.) sort vector
    std::sort(item_inventory.begin(), item_inventory.end());
    item_inventory.erase(std::unique(item_inventory.begin(), item_inventory.end()), item_inventory.end());
    std::shuffle(item_inventory.begin(), item_inventory.end(), generator);

    // 3.) fill sample vector
    std::vector<uint64_t>::const_iterator it;
    intType valid_samples = 0;
    intType I = 0;
    while (valid_samples < sample_size) {
        it = std::upper_bound(part_sum_begin, part_sum_end, item_inventory[I]);
        intType ind = std::distance(part_sum_begin, it);
        ++sample[ind];

        if (sample[ind] > min_coverage) {
            ++valid_samples;
        }
        else {
            if (sample[ind] == min_coverage)
                valid_samples += min_coverage;
        }

        ++I;
    }

    return sample;
}

intType random_seed()
{
    intType random_seed;
    std::ifstream file("/dev/random", std::ios::binary);
    if (file.is_open()) {
        char* memblock;
        intType size = sizeof(intType);
        memblock = new char[size];
        file.read(memblock, size);
        file.close();
        random_seed = *reinterpret_cast<intType*>(memblock);
        delete[] memblock;
    }
    else {
        std::cerr << "Could not open /dev/random for generating a seed!\n";
        exit(EXIT_FAILURE);
    }

    return random_seed;
}

// probability vector for H_0
static std::vector<double> prob_RT;

void run_simulation(int thread_id, intVector& rank_abundance, intVector& rank_occurrence)
{
    gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, random_seed());
    std::default_random_engine generator(random_seed());

    // RT
    intVector sample_after_RT(bins);

    // PCR
    uint64_t no_PCR_molecules;
    std::vector<uint64_t> sample_after_PCR(bins);

    // Sequencing sample
    intVector sampler_after_seq(bins);

    for (int i = 1; i <= (N_trials - 1 - thread_id) / number_threads + 1; ++i) {
        if ((thread_id == 0) && (i * number_threads % 10 == 0))
            std::cout << "Trial no.: " << i* number_threads << '\n';

        // 1.) simulate RT
        gsl_ran_multinomial(r, bins, N_RT, prob_RT.data(), sample_after_RT.data());

        // 2.) simulate PCR (branching process)
        no_PCR_molecules = 0;
        for (int j = 0; j < bins; ++j) {
            if (sample_after_RT[j]) {
                sample_after_PCR[j] = sample_after_RT[j];

                // stochastic amplification
                for (int k = 0; k < cycles; ++k) {
                    sample_after_PCR[j] += gsl_ran_binomial(r, p, sample_after_PCR[j]);
                }
                no_PCR_molecules += sample_after_PCR[j];
            }
            else {
                sample_after_PCR[j] = 0;
            }
        }

        if (thread_id == 0) {
            average_PCR_molecules += no_PCR_molecules;
        }

        // 3.) sequencing sample
        if (no_PCR_molecules > 100 * N_reads) {
            // multinomial approximation is good
            std::sort(sample_after_PCR.begin(), sample_after_PCR.end(), std::greater<intType>());
            sampler_after_seq = ran_mult_multinomial(sample_after_PCR, N_reads, r, min_coverage);
            ++used_with_replacement;
        }
        else {
            // multinomial approximation is not good enough!
            sampler_after_seq = ran_mult_hypergeometric(sample_after_PCR, N_reads, r, generator, min_coverage);
            ++used_without_replacement;
        }

        std::sort(sampler_after_seq.begin(), sampler_after_seq.end(), std::greater<intType>());

        // 4.) enter rankings
        for (int j = 0; j < rank_cutoff; ++j) {
            if (sampler_after_seq[j]) {
                rank_abundance[j] += sampler_after_seq[j];
                ++rank_occurrence[j];
            }
        }
    }

    gsl_rng_free(r);
}

int main(int argc, char* argv[])
{
    boost::program_options::options_description desc("Options");
    desc.add_options()("help,h", "Print this help")(",K", boost::program_options::value<intType>(&bins)->default_value(1048576), "Number of different primer IDs")("Nrt", boost::program_options::value<intType>(&N_RT)->default_value(10000), "Number of molecules sampled in RT")(",k", boost::program_options::value<intType>(&cycles)->default_value(30), "Number of cycles in PCR")(",p", boost::program_options::value<double>(&p)->default_value(0.8), "Efficiency of PCR")("Nmc", boost::program_options::value<intType>(&N_trials)->default_value(500), "Number of Monte Carlo simulations")("Nreads", boost::program_options::value<intType>(&N_reads)->default_value(600000), "Number of simulated reads")(",c", boost::program_options::value<intType>(&rank_cutoff)->default_value(1000), "Rank cutoff beyond which samples are discarded")(",d", boost::program_options::value<intType>(&min_coverage)->default_value(0), "Minimum coverage of reads per pID")(",t", boost::program_options::value<intType>(&number_threads)->default_value(boost::thread::hardware_concurrency()), "Number of threads to use for simulation")("quiet,q", "Suppress writing result to stdout");

    // show help options
    boost::program_options::variables_map cmdline_options;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), cmdline_options);
    boost::program_options::notify(cmdline_options);

    if (cmdline_options.count("help")) {
        std::cout << desc << "\n";
        return EXIT_SUCCESS;
    }
    if (cmdline_options.count("quiet"))
        silent = true;
    else
        silent = false;

    // show parameters
    std::cout << "No primerIDs:       " << bins << '\n';
    std::cout << "No molecules in RT: " << N_RT << '\n';
    std::cout << "No cycles in PCR:   " << cycles << '\n';
    std::cout << "Efficiency of PCR:  " << p << '\n';
    std::cout << "No MC simulations:  " << N_trials << '\n';
    std::cout << "No simulated reads: " << N_reads << '\n';
    std::cout << "Rank cutoff:        " << rank_cutoff << '\n';
    std::cout << "Minimum coverage:   " << min_coverage << '\n';
    std::cout << "No threads:         " << number_threads << '\n';
    std::cout << "Ave. PCR molecules: " << N_RT* std::pow(1 + p, cycles) << '\n';

    std::vector<intVector> rank_abundances(number_threads, intVector(rank_cutoff, 0));
    std::vector<intVector> rank_occurrences(number_threads, intVector(rank_cutoff, 0));
    prob_RT.assign(bins, 1.0 / bins);

    // start simulation
    boost::thread_group all_threads;
    for (int t = 0; t < number_threads; ++t) {
        all_threads.add_thread(new boost::thread(run_simulation, t, boost::ref(rank_abundances[t]), boost::ref(rank_occurrences[t])));
    }
    all_threads.join_all();

    if (!silent)
        std::cout << "Result:\n";

    std::stringstream id;
    id << "Nrt_" << N_RT << "_Nreads_" << N_reads << "_p_" << (boost::format("%.4f") % p);
    std::cout << id.str() << '\n';

    std::ofstream output((id.str() + ".txt").c_str());

    if (!silent)
        std::cout << "Result:\n";
    output << "DATA[[\"" << id.str() << "\"]] = c(";

    double sum_abundance, sum_occurrence, expected_abundance;
    for (int i = 0; i < rank_cutoff; ++i) {
        sum_abundance = 0;
        sum_occurrence = 0;
        for (int t = 0; t < number_threads; ++t) {
            sum_abundance += rank_abundances[t][i];
            sum_occurrence += rank_occurrences[t][i];
        }

        expected_abundance = sum_abundance / sum_occurrence;

        if (i)
            output << "," << expected_abundance;
        else
            output << expected_abundance;

        if (!silent)
            std::cout << expected_abundance << ' ';
    }
    std::cout << '\n';
    output << ")\n";
    output.close();

    // show diagnostics
    if (used_with_replacement > used_without_replacement)
        std::cout << "Sampled \033[1mwith replacement\033[0m, PCR copy numbers sufficiently large\n";
    else
        std::cout << "Sampled \033[1mwithout replacement\033[0m, PCR copy numbers low\n";

    average_PCR_molecules /= ((N_trials - 1) / number_threads + 1);
    std::cout << "Average no. PCR molecules: " << average_PCR_molecules << '\n';
    return EXIT_SUCCESS;
}