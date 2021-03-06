#include <array>
#include <iomanip>
#include <fstream>
#include <iostream>

#include <gsl/gsl_cdf.h>

#include "statistics.hpp"

#ifndef REVERSE_ENDIAN
#define _ORDER L - 1 -
#else
#define _ORDER
#endif

std::string number_To_DNA(int number, int L)
{
    // map from numeric placeholder to ASCII char
    // 0 -> 'A'
    // 1 -> 'C'
    // 2 -> 'G'
    // 3 -> 'T'
    const static char arr[] = { 'A', 'C', 'G', 'T' };

    std::string result(L, 'A');
    int j;

    for (int i = L - 1; i >= 0; --i) {
        j = std::pow(4, i);
        result[_ORDER i] = arr[number / j];
        number %= j;
    }

    return result;
}

int DNA_to_number(const std::string& DNA)
{
    // map from ASCII char to numeric placeholder
    // 'A' -> 0
    // 'C' -> 1
    // 'G' -> 2
    // 'T' -> 3
    const static char arr[] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, -1, 1,
        -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3 };
    const int L = DNA.length();

    int result = 0;
    for (int i = 0; i < L; ++i) {
        result += std::pow(4, i) * arr[static_cast<unsigned char>(DNA[_ORDER i])];
    }

    return result;
}

/* SEQ_STATISTICS */
seq_statistics::seq_statistics(int L_)
    : m_L(L_), m_uniq_counts(0), m_repl_counts(0), m_whole_pID_counts_uniq(L_), m_whole_pID_counts_repl(L_), m_observed_max_abundance(0), m_histogram_of_abundance(0, m_max_abundance + 1), m_histogram_of_pID_length(0, m_max_len_pID + 1), m_histogram_of_min_hamming_hetero_distance(0, m_max_min_hamming_hetero_distance + 1)
{
}

void seq_statistics::reset()
{
    m_uniq_counts = 0;
    m_repl_counts = 0;

    m_whole_pID_counts_uniq = 0;
    m_whole_pID_counts_repl = 0;

    m_histogram_of_abundance = 0;
    m_observed_max_abundance = 0;

    m_histogram_of_pID_length = 0;
}

void seq_statistics::addLengthToHistogram(int lengthpID)
{
    if (lengthpID <= m_max_len_pID)
        ++m_histogram_of_pID_length[lengthpID];
}

void seq_statistics::addAbundanceToHistogram(int replicates)
{
    if (replicates <= m_max_abundance) {
        ++m_histogram_of_abundance[replicates];
        m_observed_max_abundance = std::max(m_observed_max_abundance, replicates);
    }
    else {
        throw std::range_error("PCR replicate number is too large for histogram!\n");
    }
}

void seq_statistics::addMinHeteroHammingToHistogram(int distance)
{
    if (distance <= m_max_min_hamming_hetero_distance) {
        ++m_histogram_of_min_hamming_hetero_distance[distance];
    }
    else {
        throw std::range_error("Hamming distance is too large for histogram!\n");
    }
}

void seq_statistics::addPrimer(const std::string& strPrimer, int replicates)
{
    // 1.) add to total RT count
    ++m_whole_pID_counts_uniq[strPrimer];
    m_whole_pID_counts_repl[strPrimer] += replicates;

    // 2.) keep track of sum
    ++m_uniq_counts;
    m_repl_counts += replicates;
}

void seq_statistics::mergestatistics(const seq_statistics& statisticsB)
{
    // 1.) total counts
    m_uniq_counts += statisticsB.m_uniq_counts;
    m_repl_counts += statisticsB.m_repl_counts;

    // 2.) sum up total RT count
    m_whole_pID_counts_uniq += statisticsB.m_whole_pID_counts_uniq;
    m_whole_pID_counts_repl += statisticsB.m_whole_pID_counts_repl;

    // 3.) histograms
    m_histogram_of_abundance += statisticsB.m_histogram_of_abundance;
    m_observed_max_abundance = std::max(m_observed_max_abundance, statisticsB.m_observed_max_abundance);

    m_histogram_of_pID_length += statisticsB.m_histogram_of_pID_length;

    m_histogram_of_min_hamming_hetero_distance += statisticsB.m_histogram_of_min_hamming_hetero_distance;
}

void seq_statistics::write_to_csv(const std::string& file_stem) const
{
    std::string DNAid;
    std::ofstream output((file_stem + "_pID_counts.csv").c_str());

    output << "Primer,unique,replicates\n";
    for (size_t i = 0; i < m_whole_pID_counts_uniq.size(); ++i) {
        if (m_whole_pID_counts_uniq[i]) {
            DNAid = number_To_DNA(i, m_L);

            // check for indexing scheme
            if (m_whole_pID_counts_repl[DNAid] != m_whole_pID_counts_repl[i]) {
                std::cerr << "Failure!\n";
                exit(EXIT_FAILURE);
            }

            output << DNAid << ',' << m_whole_pID_counts_uniq[i] << ',' << m_whole_pID_counts_repl[i] << '\n';
        }
    }

    output.close();
}

void seq_statistics::write_histograms(const std::string& file_stem) const
{
    // pID length histogram
    std::ofstream output_pID_length((file_stem + "_pID_length_histograms.csv").c_str());
    output_pID_length << "Length,Count\n";
    for (int i = 0; i <= m_max_len_pID; ++i)
        output_pID_length << i << ',' << m_histogram_of_pID_length[i] << '\n';
    output_pID_length.close();

    // abundance histogram
    std::ofstream output_abundance((file_stem + "_abundance_histograms.csv").c_str());
    output_abundance << "Abundance,Count\n";
    for (int i = 1; i <= m_observed_max_abundance; ++i)
        output_abundance << i << ',' << m_histogram_of_abundance[i] << '\n';
    output_abundance.close();

    // minimum heterozygous Hamming distance histogram
    std::ofstream output_hetero_hamming((file_stem + "_min_hetero_hamming_histograms.csv").c_str());
    output_hetero_hamming << "Hamming,Count\n";
    for (int i = 0; i <= m_max_min_hamming_hetero_distance; ++i)
        output_hetero_hamming << i << ',' << m_histogram_of_min_hamming_hetero_distance[i] << '\n';
    output_hetero_hamming.close();
}

/* PROB_CYCLE */
double prob_cycle::operator()(int num_mt, int num_total, int num_min, int cycle) const
{
    std::map<std::tuple<int, int, int, int>, double>::const_iterator it = m_cached_results.find(std::forward_as_tuple(num_mt, num_total, num_min, cycle));
    if (it != m_cached_results.end()) {
        return it->second;
    }
    else {
        double p = p_cycle_given_X_and_constraint(num_mt, num_total, num_min, cycle);
        m_cached_results.emplace(std::make_tuple(num_mt, num_total, num_min, cycle), p);

        return p;
    }
}

int prob_cycle::getSize() const
{
    return m_cached_results.size();
}

double prob_cycle::p_cycle(int cycle) const
{
    return exp2(static_cast<double>(cycle) - 1);
}

double prob_cycle::p_X_given_cycle_and_constraint(int num_mt, int num_total, int num_min, int cycle) const
{
    double p = exp2(-static_cast<double>(cycle));

    // return gsl_ran_binomial_pdf(no_mt, p, no_total) / (gsl_cdf_binomial_P(no_total/2, p, no_total) - gsl_cdf_binomial_P(no_min, p, no_total));
    return pdf(boost::math::binomial(num_total, p), num_mt) / (cdf(boost::math::binomial(num_total, p), num_total / 2) - cdf(boost::math::binomial(num_total, p), num_min));
}

double prob_cycle::p_cycle_given_X_and_constraint(int num_mt, int num_total, int num_min, int cycle) const
{
    double sum = 0;

    for (int j = 1; j < prob_cycle::cycle_depth; ++j) {
        sum += p_X_given_cycle_and_constraint(num_mt, num_total, num_min, j) * p_cycle(j);
    }

    return p_X_given_cycle_and_constraint(num_mt, num_total, num_min, cycle) * p_cycle(cycle) / sum;
}

/* RANKED_DNA_LIST */
// INITIALIZERS:
ranked_DNA_list::ranked_DNA_list()
{
    reset();
}

void ranked_DNA_list::reset()
{
    for (const auto i : { 'A', 'C', 'G', 'T', '-', 'N' }) {
        m_base_ranks.left[i] = 0;
    }

    m_is_sorted = false;

    m_num_total = 0;
    m_num_rest = 0;
    m_num_ambig = 0;
}

// MODIFIERS:
void ranked_DNA_list::add_base(char base)
{
    // throws error if base is unrecognized
    valid_base(base);

    m_is_sorted = false;

    ++m_base_ranks.left[base];
    ++m_num_total;
}

void ranked_DNA_list::sort()
{
    if (!m_is_sorted) {
        m_base_ranks.right.sort(std::greater<int>());
        m_is_sorted = true;

        m_num_ambig = m_base_ranks.left['N'];
        m_num_rest = m_num_total - m_num_ambig;
    }
}

// ACCESSORS:
std::pair<char, int> ranked_DNA_list::get_rank_base(unsigned int i)
{
    if (i > 5)
        throw std::string("Invalid rank, has to be <6!\n");

    this->sort();
    return std::pair<char, int>(m_base_ranks.right[i].second, m_base_ranks.right[i].first);
}

int ranked_DNA_list::get_ambig_count()
{
    this->sort();
    return m_num_ambig;
}

int ranked_DNA_list::get_total_mt_count()
{
    this->sort();
    return m_num_total - m_base_ranks.right[0].first;
}

double ranked_DNA_list::get_ambig_frac()
{
    this->sort();
    return static_cast<double>(m_num_ambig) / m_num_total;
}

double ranked_DNA_list::get_nonambig_frac()
{
    this->sort();
    return static_cast<double>(m_num_rest) / m_num_total;
}

/* KMEANS */
static std::default_random_engine generator;
static std::bernoulli_distribution distribution(0.5);

std::string call_hetero_consensus(const std::vector<proper_read>& reads, const std::vector<int>& indices)
{
    std::string new_consensus;
    ranked_DNA_list ranks;

    char wt, mt;
    int num_wt, num_mt;

    for (int i : reads[0].m_ref.m_heterozygous_loci) {
        ranks.reset();

        for (int j : indices) {
            ranks.add_base(reads[j].m_fullRead[i]);
        }

        std::tie<char, int>(wt, num_wt) = ranks.get_rank_base(0);
        std::tie<char, int>(mt, num_mt) = ranks.get_rank_base(1);

        if (num_wt > num_mt)
            new_consensus.push_back(wt);
        else {
            if (distribution(generator))
                new_consensus.push_back(wt);
            else
                new_consensus.push_back(mt);
        }
    }

    return new_consensus;
}

std::string call_hetero_consensus(const std::vector<proper_read*>& reads)
{
    std::string new_consensus;
    ranked_DNA_list ranks;

    char wt, mt;
    int num_wt, num_mt;

    for (int i : reads[0]->m_ref.m_heterozygous_loci) {
        ranks.reset();

        for (auto j : reads) {
            ranks.add_base(j->m_fullRead[i]);
        }

        std::tie<char, int>(wt, num_wt) = ranks.get_rank_base(0);
        std::tie<char, int>(mt, num_mt) = ranks.get_rank_base(1);

        if (num_wt > num_mt)
            new_consensus.push_back(wt);
        else {
            if (distribution(generator))
                new_consensus.push_back(wt);
            else
                new_consensus.push_back(mt);
        }
    }

    return new_consensus;
}

std::string call_full_consensus(const std::vector<proper_read*>& reads, double minMajorFraction)
{
    std::string new_consensus(reads[0]->m_ref.m_replace_start, 'N');
    ranked_DNA_list ranks;

    char wt;
    int num_wt;
    const int cov = reads.size();

    for (int i = 0; i < reads[0]->m_ref.m_replace_start; ++i) {
        ranks.reset();

        for (auto j : reads) {
            ranks.add_base(j->m_fullRead[i]);
        }

        std::tie<char, int>(wt, num_wt) = ranks.get_rank_base(0);

        if (num_wt >= minMajorFraction * cov) {
            new_consensus[i] = wt;
        }
    }

    return new_consensus;
}

// functor
std::tuple<double, int, std::vector<proper_read*>, std::string> kmeans::operator()(std::vector<proper_read>& reads, double minMajorFraction)
{
    m_total_size = reads.size();

    // 1.) create graph of reads
    Graph G(m_total_size);
    int mismatches, valid_trials, Ns;

    for (int i = 0; i < m_total_size - 1; ++i) {
        for (int j = i + 1; j < m_total_size; ++j) {
            std::tie(mismatches, valid_trials, Ns) = reads[i].hetero_hamming_distance(reads[j]);

            if (mismatches <= max_hamming) {
                boost::add_edge(i, j, G);
            }
        }
    }

    // 2.) determine number of connected components
    std::vector<int> component(m_total_size);
    int num = boost::connected_components(G, &component[0]);

    std::vector<std::pair<int, std::vector<int>>> comp_sizes(num);
    for (int i = 0; i < m_total_size; ++i) {
        comp_sizes[component[i]].first = component[i];
        comp_sizes[component[i]].second.emplace_back(i);
    }

    // 3.) sort connected components by size
    std::sort(comp_sizes.begin(), comp_sizes.end(),
        [](const std::pair<int, std::vector<int>>& left,
                  const std::pair<int, std::vector<int>>& right) {
  									return left.second.size() > right.second.size();
        });

    // 4.) initialize cluster centers
    m_size1 = comp_sizes[0].second.size();
    m_c_mean1 = call_hetero_consensus(reads, comp_sizes[0].second);

    if (comp_sizes.size() == 1) {
        // only one connected component, definitely no collision
        m_size1 = m_total_size;
        m_size2 = 0;

        m_frac1 = 1;
        m_frac2 = 0;

        m_inter_cluster_ham = 0;
        m_c_mean2.clear();

        // copy all reads into pointers
        m_cluster1.clear();
        for (int i = 0; i < m_total_size; ++i) {
            m_cluster1.push_back(&reads[i]);
        }

        m_full_consensus1 = call_full_consensus(m_cluster1, minMajorFraction);
        return std::tuple<double, int, std::vector<proper_read*>, std::string>(m_frac1, m_size1, m_cluster1, m_full_consensus1);
    }
    else {
        m_size2 = comp_sizes[1].second.size();
        m_c_mean2 = call_hetero_consensus(reads, comp_sizes[1].second);
    }

    // 5.) perform k-means clustering
    int dist1, dist2;
    for (int i = 0; i <= 3; ++i) {
        m_cluster1.clear();
        m_cluster2.clear();

        for (int j = 0; j < m_total_size; ++j) {
            std::tie(dist1, valid_trials, Ns) = reads[j].hetero_hamming_distance(m_c_mean1);
            std::tie(dist2, valid_trials, Ns) = reads[j].hetero_hamming_distance(m_c_mean2);
            (dist1 <= dist2 ? m_cluster1 : m_cluster2).push_back(&reads[j]);
        }

        m_size1 = m_cluster1.size();
        m_size2 = m_cluster2.size();

        m_c_mean1 = call_hetero_consensus(m_cluster1);
        m_c_mean2 = call_hetero_consensus(m_cluster2);
    }

    m_frac1 = static_cast<double>(m_size1) / (m_size1 + m_size2);
    m_frac2 = static_cast<double>(m_size2) / (m_size1 + m_size2);

    std::tie(m_inter_cluster_ham, valid_trials, Ns) = hamming_distance(m_c_mean1, m_c_mean2);
    m_full_consensus1 = call_full_consensus(m_cluster1, minMajorFraction);

    return std::tuple<double, int, std::vector<proper_read*>, std::string>(m_frac1, m_size1, m_cluster1, m_full_consensus1);
}