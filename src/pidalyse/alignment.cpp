#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstddef>
#include <memory>

#include "alignment.hpp"
#include "statistics.hpp"

static prob_cycle PCR_prob;

#include <boost/algorithm/string.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/filesystem.hpp>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics_double.h>

/* ALIGNMENT */
// CTOR:
alignment::alignment(const std::string& fileName_)
    : m_input_fileName(fileName_), m_pid_stats(10)
{
    m_fileName = boost::filesystem::path(m_input_fileName).filename().string();
    m_fileStem = boost::filesystem::path(m_input_fileName).stem().string();

    m_reference = (m_fileName.find("3223") != std::string::npos ? ::ref3223 : ::ref3236);

    std::ifstream input;
    input.open(m_input_fileName.c_str());

    if (input.is_open()) {
        m_raw_reads.reserve(700000);

        // temporary objects for I/O
        std::string line_left, line_right;
        std::vector<std::string> split_vec_left, split_vec_right;
        int line_no = 2;

        // temporary objects for constructing raw read
        std::string final_raw_DNA;
        std::string pID;
        int num_inserts, num_dels, num_N;
        int len_overhang, len_pID;
        bool is_valid;

        while (input.good()) {
            if (line_no % 200000 == 0)
                std::cout << "Loaded " << line_no / 1000 << "k lines\n";

            // read left read mate
            getline(input, line_left);
            line_left.erase(std::remove(line_left.begin(), line_left.end(), '\n'), line_left.end());

            // skip SAM header
            if ((line_left.length() == 0) || (line_left[0] == '@'))
                continue;
            ++line_no;

            // read right read mate
            getline(input, line_right);
            line_right.erase(std::remove(line_right.begin(), line_right.end(), '\n'), line_right.end());
            ++line_no;

            // tokenize SAM lines
            boost::split(split_vec_left, line_left, boost::is_any_of("\t"), boost::token_compress_on);
            boost::split(split_vec_right, line_right, boost::is_any_of("\t"), boost::token_compress_on);

            // build merged read
            final_raw_DNA.reserve(m_reference.m_genome_length);
            final_raw_DNA.clear();
            pID.reserve(m_reference.m_pID_length);

            construct_sequence(split_vec_left[9], split_vec_left[5], stoi(split_vec_left[3]), m_reference, final_raw_DNA, pID, num_inserts, num_dels, num_N, len_overhang); // left
            construct_sequence(split_vec_right[9], split_vec_right[5], stoi(split_vec_right[3]), m_reference, final_raw_DNA, pID, num_inserts, num_dels, num_N, len_overhang); // right

            len_pID = std::max(static_cast<int>(pID.length()) + num_inserts - num_dels, 0);
            is_valid = (num_dels == 0) && (num_inserts == 0) && (num_N == 0) && (len_overhang >= m_reference.m_overhang_min_length) && (len_pID == m_reference.m_pID_length);

            m_raw_reads.emplace_back(std::move(final_raw_DNA), std::move(pID), num_inserts, num_dels, len_overhang, len_pID, is_valid);
        }
    }
    else {
        std::cout << m_input_fileName << " is not readable!\n";
        exit(EXIT_FAILURE);
    }
    input.close();

    std::sort(m_raw_reads.begin(), m_raw_reads.end(),
        [](const raw_paired_read& A, const raw_paired_read& B) {
			return (A.m_pID < B.m_pID);
        });

    std::cout << m_input_fileName << ": " << m_raw_reads.size() << " unique read identifiers\n";
}

// PREPROCESSING:
void alignment::filtering_QA()
{
    std::string curString;
    int accepted = 0;
    m_raw_pID_collection.clear();
    m_raw_pID_collection.reserve(600000);

    // 1. add to vector
    for (const auto& i : m_raw_reads) {
        m_reference.assign_counts(i.m_DNA, false);
        m_pid_stats.addLengthToHistogram(i.m_len_pID);

        if (i.m_is_valid) {
            // can add read
            if (i.m_pID != curString) {
                // new pID
                curString = i.m_pID;
                m_raw_pID_collection.push_back(std::pair<std::string, std::vector<proper_read>>(curString, std::vector<proper_read>()));
            }

            m_raw_pID_collection.back().second.emplace_back(i.m_DNA, m_reference);
        }
    }

    // 2. collect statistics
    for (const auto& i : m_raw_pID_collection) {
        m_pid_stats.addAbundanceToHistogram(i.second.size());
        accepted += i.second.size();
    }

    std::cout << m_fileName << ": " << accepted << " passed QA (" << std::fixed << std::setprecision(1) << static_cast<double>(accepted) / m_raw_reads.size() * 100 << "%)\n";
}

void alignment::remove_pID_collisions(int min_required_coverage, double min_plurality, bool report)
{
    m_min_current_coverage = min_required_coverage;
    m_min_majority_fraction = min_plurality;

    std::string consensus;
    m_number_collisions = 0;
    m_number_indecisive = 0;
    m_number_singletons = 0;

    m_collision_free_pID_collection.clear();
    m_consensus_pID_collection.clear();

    std::vector<proper_read*> tmp;

    for (auto& i : m_raw_pID_collection) {
        tmp.clear();
        consensus = call_consensus_and_remove_collisions(i.second, 50, tmp);

        switch (consensus[0]) {
        case '0':
            // collision
            ++m_number_collisions;
            break;

        case '1':
            // indecisive
            ++m_number_indecisive;
            break;

        default:
            // OK, proper consensus
            ++m_number_singletons;

            m_reference.assign_counts(consensus, true);
            m_pid_stats.addPrimer(i.first, tmp.size());

            m_collision_free_pID_collection.emplace_back(i.first, std::move(tmp));
            m_consensus_pID_collection.emplace_back(std::piecewise_construct,
                std::forward_as_tuple(i.first),
                std::forward_as_tuple(std::move(consensus), m_reference, tmp.size()));
            break;
        }
    }

    m_reference.normalise_counts();

    if (report) {
        std::cout << "Collision-free pIDs: " << m_number_singletons << '\n';
        std::cout << "    Indecisive pIDs: " << m_number_indecisive << '\n';
        std::cout << "     Collision pIDs: " << m_number_collisions << " (" << std::fixed << std::setprecision(1) << static_cast<double>(m_number_collisions) / (m_number_singletons + m_number_collisions) * 100 << "%)\n";
        std::cout << "         Total pIDs: " << m_raw_pID_collection.size() << '\n';
    }
}

// RT STUFF:
hamming_return_type alignment::count_mismatches_at_locus(int locus) const
{
    int mismatches = 0;
    int valid_trials = 0;
    int Ns = 0;

    for (const auto& i : m_consensus_pID_collection) {
        if (valid_base(i.second.m_fullConsensus[locus])) {
            ++valid_trials;
            mismatches += (i.second.m_fullConsensus[locus] != i.second.m_ref.m_all_reference_strains[0].m_DNA[locus]);
        }
        else {
            ++Ns;
        }
    }

    return hamming_return_type(mismatches, valid_trials, Ns);
}

hamming_return_type alignment::calculate_RT_mismatches() const
{
    int mismatches = 0;
    int valid_trials = 0;
    int Ns = 0;
    hamming_return_type temp;

    for (const auto& i : m_consensus_pID_collection) {
        temp = i.second.calculate_homozygous_mismatches();

        mismatches += std::get<0>(temp);
        valid_trials += std::get<1>(temp);
        Ns += std::get<2>(temp);
    }

    return hamming_return_type(mismatches, valid_trials, Ns);
}

double alignment::LogLik(double s, double r) const
{
    double total_likelihood = 0;
    double temp;

    if (m_reference.m_freq_initialised == false) {
        std::cout << "Cannot calculate Log Likelihood without estimated reference frequencies!\n";
        exit(EXIT_FAILURE);
    }

    std::map<std::string, double> string_to_prob_cache;
    std::map<std::string, double>::const_iterator it;

    for (const auto& i : m_consensus_pID_collection) {
        it = string_to_prob_cache.find(i.second.m_heterozygous_loci_string);
        if (it == string_to_prob_cache.end()) {
            // new sequence
            temp = i.second.log_prob(s, r);
            string_to_prob_cache.emplace(i.second.m_heterozygous_loci_string, temp);
            total_likelihood += temp;
        }
        else {
            total_likelihood += it->second;
        }
    }

    return total_likelihood;
}

// PCR STUFF:
d_hamming_return_type alignment::calculate_PCR_mismatches() const
{
    double mismatches = 0;
    uint64_t valid_trials = 0;
    uint64_t Ns = 0;

    int cov;
    ranked_DNA_list ranks;

    char wt, mt;
    int num_wt, num_mt;

    double p;

    for (const auto& i : m_collision_free_pID_collection) {
        cov = i.second.size();
        // PCRerror = false;

        for (int j : m_reference.m_homozygous_loci) {
            // 1.) perform pileup
            ranks.reset();

            for (const auto& k : i.second) {
                ranks.add_base(k->m_fullRead[j]);
            }

            // 2.) perform ranking
            std::tie<char, int>(wt, num_wt) = ranks.get_rank_base(0);
            std::tie<char, int>(mt, num_mt) = ranks.get_rank_base(1);

            if (ranks.get_ambig_count() >= 2) {
                ++Ns;
                continue;
            }

            ++valid_trials;

            if (ranks.get_total_mt_count() <= 2)
                continue;

            p = PCR_prob(num_mt, num_mt + num_wt, 2);
            mismatches += p;
        }
    }

    return d_hamming_return_type(mismatches, valid_trials, Ns);
}

// STATISTICS STUFF:
void alignment::add_to_merge_statistics(seq_statistics& merge_into) const
{
    merge_into.mergestatistics(this->m_pid_stats);
}

// DISPLAY I/O:
void alignment::show_recombination_patterns() const
{
    for (const auto& i : m_consensus_pID_collection) {
        if (i.second.m_best_reference == i.second.m_ref.m_K) {
            // recombinant
            std::cout << i.first << " (" << i.second.m_multiplicity << ")\n";

            for (const auto& j : i.second.m_ref.m_all_reference_strains) {
                std::cout << j.m_name << '\t';

                for (int k : i.second.m_ref.m_heterozygous_loci) {
                    if ((i.second.m_fullRead[k] == 'N') || (i.second.m_fullRead[k] == j.m_heterozygous_loci_string[k]))
                        std::cout << i.second.m_fullRead[k];
                    else
                        std::cout << ' ';
                }

                std::cout << '\n';
            }

            std::cout << '\n';
        }
    }
}

void alignment::show_clone_frequencies() const
{
    m_reference.display_strains_abridged(m_fileStem);
}

// FILE I/O:
void alignment::write_consensus_to_fasta() const
{
    std::string output_fileName = m_fileStem + "_cons.fasta";

    std::ofstream output(output_fileName.c_str());
    int k = 0;
    for (const auto& i : m_consensus_pID_collection) {
        ++k;
        output << '>' << k
               << "_pID:" << i.first
               << "_m:" << i.second.m_multiplicity
               << "_T:" << (i.second.m_best_reference == i.second.m_ref.m_K ? "RECOMB" : i.second.m_ref.m_all_reference_strains[i.second.m_best_reference].m_name)
               << "_H:" << i.second.m_hamming_distance_to_best_reference << '\n'
               << i.second.m_fullConsensus << '\n';
    }

    output.close();
}

void alignment::write_to_fasta() const
{
    std::string output_fileName = m_fileStem + ".fasta";

    std::ofstream output(output_fileName.c_str());
    int k = 0, m;
    for (const auto& i : m_raw_pID_collection) {
        ++k;
        m = 0;
        for (const auto& j : i.second) {
            ++m;
            output << '>' << k << ',' << m
                   << "_pID:" << i.first << '\n'
                   << j.m_fullRead << '\n';
        }
    }

    output.close();
}

void alignment::write_raw_to_fasta() const
{
    std::string output_fileName = m_fileStem + "_raw.fasta";

    std::ofstream output(output_fileName.c_str());
    int k = 0;
    for (const auto& i : m_raw_reads) {
        ++k;
        output << '>' << k
               << "_pID:" << (i.m_pID.length() == 0 ? "NA" : i.m_pID)
               << "_i:" << i.m_num_inserts
               << "_d:" << i.m_num_dels
               << "_l:" << i.m_len_overhang
               << "_v:" << (i.m_is_valid ? 'T' : 'F') << '\n'
               << i.m_DNA << '\n';
    }

    output.close();
}

void alignment::write_frequencies_to_MATLAB(std::ofstream& output, bool pID) const
{
    const std::vector<record>& strains = m_reference.m_all_reference_strains;

    output << std::fixed << std::setprecision(6) << (pID ? strains[0].m_pID_frequency : strains[0].m_raw_frequency);

    for (int i = 1; i < m_reference.m_K; ++i) {
        output << ',' << (pID ? strains[i].m_pID_frequency : strains[i].m_raw_frequency);
    }
    output << ";\n";
}

void alignment::write_statistics_histograms() const
{
    m_pid_stats.write_histograms(m_fileStem);
}

void alignment::write_statistics_to_csv() const
{
    m_pid_stats.write_to_csv(m_fileStem);
}

// PRIVATE FUNCTIONS:
std::string alignment::call_consensus_and_remove_collisions(std::vector<proper_read>& reads, int minDisplay, std::vector<proper_read*>& filtered_reads)
{
    static kmeans clustering;
    std::string full_consensus;

    if (reads.size() < m_min_current_coverage) {
        // indecisive - coverage too low
        full_consensus = "1";
    }
    else {
        double majority_fraction;
        int majority_size;
        std::tie(majority_fraction, majority_size, filtered_reads, full_consensus) = clustering(reads, m_min_majority_fraction);

        if (majority_fraction < m_min_majority_fraction) {
            // collision
            full_consensus = "0";
        }
        else {
            // no collision
            if (majority_size < m_min_current_coverage) {
                // indecisive - coverage too low
                full_consensus = "1";
            }
        }
    }

    return full_consensus;
}

void alignment::construct_sequence(const std::string& SEQ, const std::string& CIGAR, const int POS, const reference& ref, std::string& final_raw_sequence, std::string& pID, int& num_inserts, int& num_dels, int& num_N, int& len_overhang) const
{
    // current position in CIGAR
    int curPos = 0;

    // current position in READ
    int curOffset = 0;

    // current position on GENOME
    // SAM is 1-based
    int offsetRef = POS - 1;

    // current length of stretch in READ
    int curLen;
    num_inserts = 0;

    final_raw_sequence.append(offsetRef - final_raw_sequence.length(), '-');

    // 1.) tokenize CIGAR
    for (size_t i = 0; i < CIGAR.length(); ++i) {
        if (isalpha(CIGAR[i])) {
            curLen = stoi(CIGAR.substr(curPos, i - curPos));
            curPos = i + 1;

            switch (CIGAR[i]) {
            case 'M':
                final_raw_sequence.append(SEQ, curOffset, curLen);
                curOffset += curLen;
                offsetRef += curLen;
                break;

            case 'D':
                final_raw_sequence.append(curLen, '-');
                offsetRef += curLen;
                break;

            case 'I':
                if (offsetRef == ref.m_pID_start) {
                    num_inserts += curLen;
                }

                curOffset += curLen;
                break;

            case 'S':
                curOffset += curLen;
                break;

            case 'H':
                curOffset += curLen;
                break;
            }
        }
    }

    // 2.) calculate statistics for reads
    if (POS > 260) {
        // right read
        if (offsetRef > ref.m_pID_start) {
            // has a pID
            pID = final_raw_sequence.substr(ref.m_pID_start, ref.m_pID_length);
            num_dels = std::count(pID.begin(), pID.end(), '-');
            num_N = std::count(pID.begin(), pID.end(), 'N');
        }
        else {
            // has no pID
            pID.clear();
            num_dels = 0;
        }

        len_overhang = std::max(offsetRef - ref.m_overhang_start, 0);

        if (offsetRef > ref.m_replace_start)
            final_raw_sequence.erase(ref.m_replace_start);

        //final_raw_sequence.append(ref.genome_length - final_raw_sequence.length(), '-');
    }
}

/* ALIGNMENTS */
// CTOR:
alignments::alignments(const std::vector<std::string>& inputFiles_)
    : m_num_alignments(inputFiles_.size()), m_pid_stats(10), m_min_current_coverage(-1), m_min_majority_fraction(-1), m_is_collisions_removed(false)
{
    m_collections_alignments.reserve(m_num_alignments);

    for (const std::string& i : inputFiles_) {
        std::cout << "Loading " << i << '\n';
        m_collections_alignments.emplace_back(i);
    }
}

// PREPROCESSING:
void alignments::filtering_QA()
{
    for (alignment& i : m_collections_alignments) {
        i.filtering_QA();
    }

    m_is_collisions_removed = false;
}

void alignments::remove_pID_collisions(int min_required_coverage, double min_plurality, bool report)
{
    m_min_current_coverage = min_required_coverage;
    m_min_majority_fraction = min_plurality;

    for (alignment& i : m_collections_alignments) {
        i.remove_pID_collisions(m_min_current_coverage, m_min_majority_fraction, report);
    }

    m_is_collisions_removed = true;
}

// RT STUFF:
double alignments::estimate_RT_substitution_rate(bool report)
{
    uint64_t MLE_mismatches = 0;
    uint64_t MLE_valid_trials = 0;
    uint64_t MLE_Ns = 0;
    hamming_return_type temp;

    // 1.) MLE
    for (const alignment& i : m_collections_alignments) {
        temp = i.calculate_RT_mismatches();

        MLE_mismatches += std::get<0>(temp);
        MLE_valid_trials += std::get<1>(temp);
        MLE_Ns += std::get<2>(temp);
    }

    double RT_sub_rate_MLE = static_cast<double>(MLE_mismatches) / MLE_valid_trials;
    double RT_sub_rate_MLE_CI_lower = boost::math::binomial_distribution<>::find_lower_bound_on_p(MLE_valid_trials, MLE_mismatches, 0.025);
    double RT_sub_rate_MLE_CI_upper = boost::math::binomial_distribution<>::find_upper_bound_on_p(MLE_valid_trials, MLE_mismatches, 0.025);

    // 2.) MoM-estimator
    uint64_t mismatches;
    uint64_t valid_trials;
    uint64_t Ns;

    double temp_rate;
    std::vector<double> vector_mt_freqs, log_vector_mt_freqs;

    std::cout << std::scientific << std::setprecision(3);

    for (const auto& i : reference::m_shared_homozygous_loci) {
        mismatches = 0;
        valid_trials = 0;
        Ns = 0;

        for (const auto& j : m_collections_alignments) {
            temp = j.count_mismatches_at_locus(i);

            mismatches += std::get<0>(temp);
            valid_trials += std::get<1>(temp);
            Ns += std::get<2>(temp);
        }

        temp_rate = static_cast<double>(mismatches) / valid_trials;
        vector_mt_freqs.push_back(temp_rate);
        log_vector_mt_freqs.push_back(log(temp_rate));
    }

    double sample_size = vector_mt_freqs.size();

    // arithmetic mean:
    double MoM_mu = gsl_stats_mean(vector_mt_freqs.data(), 1, sample_size);
    double MoM_sigma = gsl_stats_variance_m(vector_mt_freqs.data(), 1, sample_size, MoM_mu);

    double RT_sub_rate_MoM = MoM_mu;
    double RT_sub_rate_MoM_CI_lower = MoM_mu - gsl_cdf_tdist_Pinv(1 - 0.05 / 2, sample_size - 1) * MoM_sigma / sqrt(sample_size);
    double RT_sub_rate_MoM_CI_upper = MoM_mu + gsl_cdf_tdist_Pinv(1 - 0.05 / 2, sample_size - 1) * MoM_sigma / sqrt(sample_size);

    // geometric mean:
    double log_MoM_mu = gsl_stats_mean(log_vector_mt_freqs.data(), 1, sample_size);
    double log_MoM_sigma = gsl_stats_variance_m(log_vector_mt_freqs.data(), 1, sample_size, log_MoM_mu);

    double RT_sub_rate_log_MoM = exp(log_MoM_mu);
    double RT_sub_rate_log_MoM_CI_lower = exp(log_MoM_mu - gsl_cdf_tdist_Pinv(1 - 0.05 / 2, sample_size - 1) * log_MoM_sigma / sqrt(sample_size));
    double RT_sub_rate_log_MoM_CI_upper = exp(log_MoM_mu + gsl_cdf_tdist_Pinv(1 - 0.05 / 2, sample_size - 1) * log_MoM_sigma / sqrt(sample_size));

    // report
    if (report) {
        std::cout << std::string(50, '=') << '\n';
        std::cout << "RT substitution rate" << '\n';
        std::cout << "-----------------------\n";
        std::cout << '\n';

        std::cout << "   1.) MLE:\n";
        std::cout << "           mt bases: " << MLE_mismatches << '\n';
        std::cout << "        Total bases: " << MLE_valid_trials << '\n';
        std::cout << "          'N' bases: " << MLE_Ns << '\n';
        std::cout << '\n';
        std::cout << "   Est. RT sub rate: " << std::scientific << std::setprecision(2) << RT_sub_rate_MLE << '\n';
        std::cout << "             95% CI: [" << RT_sub_rate_MLE_CI_lower << ", " << RT_sub_rate_MLE_CI_upper << "]\n\n";

        std::cout << "   2.) MoM:\n";
        std::cout << "   Est. RT sub rate: " << std::scientific << std::setprecision(2) << RT_sub_rate_MoM << '\n';
        std::cout << "             95% CI: [" << RT_sub_rate_MoM_CI_lower << ", " << RT_sub_rate_MoM_CI_upper << "]\n\n";
        std::cout << "   using logarithmic transform:\n";
        std::cout << "   Est. RT sub rate: " << std::scientific << std::setprecision(2) << RT_sub_rate_log_MoM << '\n';
        std::cout << "             95% CI: [" << RT_sub_rate_log_MoM_CI_lower << ", " << RT_sub_rate_log_MoM_CI_upper << "]\n\n";

        std::cout << std::string(50, '=') << '\n';
    }

    m_RT_sub_rate = RT_sub_rate_MLE;
    m_RT_sub_rate_CI_lower = RT_sub_rate_MLE_CI_lower;
    m_RT_sub_rate_CI_upper = RT_sub_rate_MLE_CI_upper;

    return m_RT_sub_rate;
}

double alignments::estimate_RT_recombination_rate(bool report)
{
    return estimate_RT_recombination_rate(m_RT_sub_rate, report);
}

double alignments::estimate_RT_recombination_rate(double s, bool report)
{
    std::pair<double, double> result;

    // 1.) calculate MLE
    result = boost::math::tools::brent_find_minima([&](double r) -> double {return this->neg_LogLik_recombination(s, r);
    }, 1E-10, 0.9, 1000);

    m_RT_recomb_rate = result.first;
    m_RT_LogLik = -result.second;

    // 2.) calculate 95% CI
    uintmax_t max_iter = 1000;
    boost::math::tools::eps_tolerance<double> tol(10000);

    // lower CI boundary
    max_iter = 1000;
    result = boost::math::tools::toms748_solve([&](double r) -> double {return this->calculate_RT_LogLik_recombination(s, r) - m_RT_LogLik + 1.920729;
    }, m_RT_recomb_rate / 50, m_RT_recomb_rate, tol, max_iter);

    m_RT_recomb_rate_CI_lower = (result.first + result.second) / 2;

    // upper CI boundary
    max_iter = 1000;
    result = boost::math::tools::toms748_solve([&](double r) -> double {return this->calculate_RT_LogLik_recombination(s, r) - m_RT_LogLik + 1.920729;
    }, m_RT_recomb_rate, m_RT_recomb_rate * 50, tol, max_iter);

    m_RT_recomb_rate_CI_upper = (result.first + result.second) / 2;

    if (report) {
        // 3.) display
        std::cout << std::string(50, '=') << '\n';
        std::cout << "RT recombination rate" << '\n';
        std::cout << "-----------------------\n";
        std::cout << '\n';

        std::cout << "Est. RT recomb rate: " << std::scientific << std::setprecision(2) << m_RT_recomb_rate << '\n';
        std::cout << "     Log-Likelihood: " << std::fixed << std::setprecision(2) << m_RT_LogLik << '\n';
        std::cout << "             95% CI: [" << std::scientific << std::setprecision(2) << m_RT_recomb_rate_CI_lower << ", " << m_RT_recomb_rate_CI_upper << "]\n";

        std::cout << std::string(50, '=') << '\n';
    }

    return m_RT_recomb_rate;
}

void alignments::plot_RT_recombination_LogLik(double upper, int n) const
{
    std::pair<double, double> result;
    double CI_LogLik = calculate_RT_LogLik_recombination(m_RT_sub_rate, upper);

    uintmax_t max_iter = 1000;
    boost::math::tools::eps_tolerance<double> tol(10000);

    // upper plot boundary
    result = boost::math::tools::toms748_solve([&](double r) -> double {return this->calculate_RT_LogLik_recombination(m_RT_sub_rate, r) - CI_LogLik;
    }, m_RT_recomb_rate / 100, m_RT_recomb_rate, tol, max_iter);

    double lower = (result.first + result.second) / 2;

    std::cout << "Drawing Plot from " << lower << " to " << upper << "\n";
    double factor = 1.0 / (n - 1) * (log(upper) - log(lower));
    factor = exp(factor);

    std::vector<double> x;
    x.reserve(n + 10);

    std::vector<double> y;
    y.reserve(n + 10);

    double temp;

    int I = 1;
    for (double i = lower; i <= upper; i *= factor, ++I) {
        x.emplace_back(i);

        temp = calculate_RT_LogLik_recombination(m_RT_sub_rate, i);
        y.emplace_back(temp);

        if (I % 10 == 0)
            std::cout << "Iteration " << std::fixed << std::setprecision(0) << I << "\tr: " << std::scientific << std::setprecision(4) << i << "\tLogLik: " << temp << '\n';
    }

    // write R file
    std::ofstream R_Data("LogLikData.R");

    R_Data << "x = c(" << std::scientific << std::setprecision(4) << x[0];
    for (size_t i = 1; i < x.size(); ++i)
        R_Data << ", " << x[i];

    R_Data << ")\n";

    R_Data << "y = c(" << std::fixed << std::setprecision(4) << y[0];
    for (size_t i = 1; i < y.size(); ++i)
        R_Data << ", " << y[i];

    R_Data << ")\n\n";

    R_Data << "LOGLIKMAX = " << std::fixed << std::setprecision(4) << m_RT_LogLik << "\n";
    R_Data << "X_CI_LOW = " << std::scientific << std::setprecision(4) << m_RT_recomb_rate_CI_lower << "\n";
    R_Data << "X_CI_HIGH = " << std::scientific << std::setprecision(4) << m_RT_recomb_rate_CI_upper << "\n";
    R_Data << "r = " << std::scientific << std::setprecision(2) << m_RT_recomb_rate << "\n";

    R_Data.close();
}

// PCR STUFF:
double alignments::estimate_PCR_substitution_rate(bool report)
{
    double mismatches = 0;
    uint64_t valid_trials = 0;
    uint64_t Ns = 0;
    d_hamming_return_type temp;

    for (const alignment& i : m_collections_alignments) {
        temp = i.calculate_PCR_mismatches();

        mismatches += std::get<0>(temp);
        valid_trials += std::get<1>(temp);
        Ns += std::get<2>(temp);
    }

    m_PCR_sub_rate = mismatches / valid_trials;
    m_PCR_sub_rate_CI_lower = boost::math::binomial_distribution<>::find_lower_bound_on_p(valid_trials, mismatches, 0.025);
    m_PCR_sub_rate_CI_upper = boost::math::binomial_distribution<>::find_upper_bound_on_p(valid_trials, mismatches, 0.025);

    if (report) {
        std::cout << std::string(50, '=') << '\n';
        std::cout << "PCR substitution rate" << '\n';
        std::cout << "-----------------------\n";
        std::cout << '\n';

        std::cout << "           mt bases: " << mismatches << '\n';
        std::cout << "        Total bases: " << valid_trials << '\n';
        std::cout << "          'N' bases: " << Ns << '\n';
        std::cout << '\n';

        std::cout << "  Est. PCR sub rate: " << std::scientific << std::setprecision(2) << m_PCR_sub_rate << '\n';
        std::cout << "             95% CI: [" << m_PCR_sub_rate_CI_lower << ", " << m_PCR_sub_rate_CI_upper << "]\n";

        std::cout << std::string(50, '=') << '\n';
    }

    return m_PCR_sub_rate;
}

// DISPLAY I/O:
void alignments::show_recombination_patterns() const
{
    for (const alignment& i : m_collections_alignments) {
        i.show_recombination_patterns();
    }
}

void alignments::show_clone_frequencies() const
{
    for (const alignment& i : m_collections_alignments) {
        i.show_clone_frequencies();
    }
}

// FILE I/O:
void alignments::write_all_consensus_to_fasta() const
{
    for (const alignment& i : m_collections_alignments) {
        i.write_consensus_to_fasta();
    }
}

void alignments::write_all_to_fasta() const
{
    for (const alignment& i : m_collections_alignments) {
        i.write_to_fasta();
    }
}

void alignments::write_all_raw_to_fasta() const
{
    for (const alignment& i : m_collections_alignments) {
        i.write_raw_to_fasta();
    }
}

void alignments::write_frequencies_to_MATLAB() const
{
    std::ofstream output("Variances.m");

    // 1.) raw frequencies
    output << "Data_raw = [...\n";
    for (const alignment& i : m_collections_alignments) {
        i.write_frequencies_to_MATLAB(output, false);
    }
    output << "];\n\n";

    // 2.) pID frequencies
    output << "Data_pID = [...\n";
    for (const alignment& i : m_collections_alignments) {
        i.write_frequencies_to_MATLAB(output, true);
    }
    output << "];\n\n";

    output.close();
}

void alignments::write_all_statistics()
{
    m_pid_stats.reset();

    for (const alignment& i : m_collections_alignments) {
        i.add_to_merge_statistics(this->m_pid_stats);
        i.write_statistics_histograms();
        i.write_statistics_to_csv();
    }

    m_pid_stats.write_to_csv("Total");
}

// PRIVATE FUNCTIONS:
double alignments::calculate_RT_LogLik_recombination(double s, double r) const
{
    double total_likelihood = 0;

    for (const alignment& i : m_collections_alignments)
        total_likelihood += i.LogLik(s, r);

    return total_likelihood;
}

double alignments::neg_LogLik_recombination(double s, double r) const
{
    return -calculate_RT_LogLik_recombination(s, r);
}
