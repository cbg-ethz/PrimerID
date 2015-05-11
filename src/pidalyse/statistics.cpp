#include <array>
#include <iomanip>
#include <fstream>
#include <iostream>

#include <gsl/gsl_cdf.h>

#include "statistics.hpp"

#ifdef REVERSE_ENDIAN
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
    : m_L(L_), m_unique_counts(0), m_repl_counts(0), m_marginal_counts(L_, DNAvector<int>(1)), m_marginal_counts_repl(L_, DNAvector<int>(1)), m_pairwise_counts(L_ - 1, DNAvector<int>(2)), m_pairwise_counts_repl(L_ - 1, DNAvector<int>(2)), m_whole_pID_counts(L_), m_whole_pID_counts_repl(L_), m_collision_free_whole_pID_counts(L_), m_collision_free_whole_pID_counts_repl(L_), m_histogram_of_abundance(0, 10000), m_max_abundance(0), m_histogram_of_pID_length(0, m_max_len_pID + 1)
{
}

void seq_statistics::reset()
{
    m_marginal_counts.assign(m_L, DNAvector<int>(1));
    m_marginal_counts_repl.assign(m_L, DNAvector<int>(1));

    m_pairwise_counts.assign(m_L - 1, DNAvector<int>(2));
    m_pairwise_counts_repl.assign(m_L - 1, DNAvector<int>(2));

    m_whole_pID_counts = DNAvector<int>(1);
    m_whole_pID_counts_repl = DNAvector<int>(1);

    m_histogram_of_abundance = 0;
}

void seq_statistics::show_statistics() const
{
    const static std::array<char, 4> DNA{ { 'A', 'C', 'G', 'T' } };

    for (size_t i = 0; i < 4; ++i) {
        std::cout << std::fixed << std::setprecision(2) << DNA[i] << ": " << m_marginal_counts[0][i];

        for (size_t j = 1; j < m_L; ++j) {
            std::cout << '\t' << m_marginal_counts[j][i];
        }

        std::cout << '\n';
    }
}

void seq_statistics::addLengthToHistogram(int lengthpID)
{
    if (lengthpID <= m_max_len_pID)
        ++m_histogram_of_pID_length[lengthpID];
}

void seq_statistics::addPrimer(const std::string& strPrimer, int replicates)
{
    // 1.) add to main effects
    for (size_t i = 0; i < m_L; ++i) {
        ++m_marginal_counts[i][strPrimer.substr(i, 1)];
        m_marginal_counts_repl[i][strPrimer.substr(i, 1)] += replicates;
    }

    // 2.) add to pairwise effects
    for (size_t i = 0; i < m_L - 1; ++i) {
        ++m_pairwise_counts[i][strPrimer.substr(i, 2)];
        m_pairwise_counts_repl[i][strPrimer.substr(i, 2)] += replicates;
    }

    // 3.) add to total RT count
    ++m_whole_pID_counts[strPrimer];
    m_whole_pID_counts_repl[strPrimer] += replicates;

    // 4.) keep track of sum
    ++m_unique_counts;
    m_repl_counts += replicates;

    // 5.) add to histogram
    ++m_histogram_of_abundance[replicates];
    m_max_abundance = std::max(m_max_abundance, replicates);
}

void seq_statistics::addPrimer_collisionFree(const std::string& strPrimer, int replicates)
{
    // 1.) add to total RT count
    ++m_collision_free_whole_pID_counts[strPrimer];
    m_collision_free_whole_pID_counts_repl[strPrimer] += replicates;
}

void seq_statistics::mergestatistics(const seq_statistics& statisticsB)
{
    // 1.) sum up main effects
    for (size_t i = 0; i < m_L; ++i) {
        m_marginal_counts[i] += statisticsB.m_marginal_counts[i];
        m_marginal_counts_repl[i] += statisticsB.m_marginal_counts_repl[i];
    }

    // 2.) sum up pairwise effects
    for (size_t i = 0; i < m_L - 1; ++i) {
        m_pairwise_counts[i] += statisticsB.m_pairwise_counts[i];
        m_pairwise_counts_repl[i] += statisticsB.m_pairwise_counts_repl[i];
    }

    // 3.) sum up total RT count
    m_whole_pID_counts += statisticsB.m_whole_pID_counts;
    m_whole_pID_counts_repl += statisticsB.m_whole_pID_counts_repl;

    m_collision_free_whole_pID_counts += statisticsB.m_collision_free_whole_pID_counts;
    m_collision_free_whole_pID_counts_repl += statisticsB.m_collision_free_whole_pID_counts_repl;

    // 4.) histograms
    m_histogram_of_abundance += statisticsB.m_histogram_of_abundance;
    m_max_abundance = std::max(m_max_abundance, statisticsB.m_max_abundance);

    m_histogram_of_pID_length += statisticsB.m_histogram_of_pID_length;

    // 5.) total counts
    m_unique_counts += statisticsB.m_unique_counts;
    m_repl_counts += statisticsB.m_repl_counts;
}

void seq_statistics::write_to_csv(const std::string& file_stem) const
{
    //const static std::array<char, 4> DNA{{'A','C','G','T'}};

    std::ofstream output((file_stem + "_unique.csv").c_str());
    std::ofstream output_repl((file_stem + "_replicates.csv").c_str());

    for (size_t i = 0; i < 4; ++i) {
        output << std::fixed << std::setprecision(5) << m_marginal_counts[0][i] / static_cast<double>(m_unique_counts);
        output_repl << std::fixed << std::setprecision(5) << m_marginal_counts_repl[0][i] / static_cast<double>(m_repl_counts);

        for (size_t j = 1; j < m_L; ++j) {
            output << ',' << m_marginal_counts[j][i] / static_cast<double>(m_unique_counts);
            output_repl << ',' << m_marginal_counts_repl[j][i] / static_cast<double>(m_repl_counts);
        }

        output << '\n';
        output_repl << '\n';
    }

    output.close();
    output_repl.close();

    // write list of primers and their PCR abundances
    std::ofstream output_primers((file_stem + "_primers.csv").c_str());
    std::ofstream output_primers_collision_free((file_stem + "_primers_collision_free.csv").c_str());
    std::string tmp;

    output_primers << "Primer,Count\n";
    output_primers_collision_free << "Primer,Count\n";
    for (size_t i = 0; i < m_whole_pID_counts_repl.size(); ++i) {
        if (m_whole_pID_counts_repl[i]) {
            tmp = number_To_DNA(i, m_L);
            if (m_whole_pID_counts_repl[tmp] != m_whole_pID_counts_repl[i]) {
                std::cerr << "Failure!\n";
                exit(EXIT_FAILURE);
            }
            output_primers << tmp << ',' << m_whole_pID_counts_repl[i] << '\n';
        }

        if (m_collision_free_whole_pID_counts_repl[i]) {
            tmp = number_To_DNA(i, m_L);
            if (m_collision_free_whole_pID_counts_repl[tmp] != m_collision_free_whole_pID_counts_repl[i]) {
                std::cerr << "Failure!\n";
                exit(EXIT_FAILURE);
            }
            output_primers_collision_free << tmp << ',' << m_collision_free_whole_pID_counts_repl[i] << '\n';
        }
    }
    output_primers.close();
    output_primers_collision_free.close();
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
    for (int i = 1; i <= m_max_abundance; ++i)
        output_abundance << i << ',' << m_histogram_of_abundance[i] << '\n';
    output_abundance.close();
}

void seq_statistics::calculate_comprehensive_statistics(const std::string& strLabel) const
{
    const static std::array<std::string, 4> DNA{ { "A", "C", "G", "T" } };
    const static std::array<std::string, 3> trun_DNA{ { "C", "G", "T" } };

    std::vector<DNAvector<double>> _marginal_probabilities(m_L, DNAvector<double>(1));

    std::vector<DNAvector<double>> _pairwise_probabilities(m_L - 1, DNAvector<double>(2));

    // 1.) calculate (MLE) probabilities
    for (size_t i = 0; i < m_L; ++i) {
        // marginals
        for (size_t j = 0; j < m_marginal_counts[i].size(); ++j) {
            _marginal_probabilities[i][j] = m_marginal_counts[i][j] / static_cast<double>(m_unique_counts);
        }

        // pairwise
        if (i != m_L - 1) {
            for (size_t j = 0; j < m_pairwise_counts[i].size(); ++j) {
                _pairwise_probabilities[i][j] = m_pairwise_counts[i][j] / static_cast<double>(m_unique_counts);
            }
        }
    }

    // 2.) print results
    if (!strLabel.empty()) {
        std::cout << strLabel << ":";
    }

    // numbering
    for (size_t i = 0; i < m_L; ++i)
        std::cout << "\t\t" << i + 1 << "\t";
    std::cout << '\n';

    for (size_t i = 0; i < m_L; ++i)
        std::cout << "\tp=\tbeta=\t";
    std::cout << '\n';

    double beta;
    /* marginals */
    for (const auto& j : DNA) {
        std::cout << j << ":";
        for (size_t i = 0; i < m_L; ++i) {
            beta = std::log(_marginal_probabilities[i][j] / _marginal_probabilities[i]["A"]);

            std::cout << "\t" << std::fixed << std::setprecision(4) << _marginal_probabilities[i][j] << "\t" << (beta < 0 ? "" : " ") << beta << "\t";
        }
        std::cout << '\n';
    }

    std::vector<DNAvector<double>> beta_1(m_L - 1, DNAvector<double>(1));
    std::vector<DNAvector<double>> beta_2(m_L - 1, DNAvector<double>(1));

    // main effects
    // 1st locus
    for (const auto& j : trun_DNA) {
        std::cout << "b_1(" << j << "):\t";

        for (size_t i = 0; i < m_L - 1; ++i) {
            beta = std::log(_pairwise_probabilities[i][j + "A"] / _pairwise_probabilities[i]["AA"]);
            beta_1[i][j] = beta;
            std::cout << "\t\t" << (beta < 0 ? "" : " ") << beta << "\t";
        }

        std::cout << '\n';
    }

    // 2nd locus
    for (const auto& j : trun_DNA) {
        std::cout << "b_2(" << j << "):\t";

        for (size_t i = 0; i < m_L - 1; ++i) {
            beta = std::log(_pairwise_probabilities[i]["A" + j] / _pairwise_probabilities[i]["AA"]);
            beta_2[i][j] = beta;
            std::cout << "\t\t" << (beta < 0 ? "" : " ") << beta << "\t";
        }

        std::cout << '\n';
    }

    // interaction effects
    for (const auto& j2 : trun_DNA) {
        for (const auto& j1 : trun_DNA) {
            std::cout << "b_12(" << j1 << j2 << "):";

            for (size_t i = 0; i < m_L - 1; ++i) {
                beta = std::log(_pairwise_probabilities[i][j1 + j2] / _pairwise_probabilities[i]["AA"]) - beta_1[i][j1] - beta_2[i][j2];
                std::cout << "\t\t" << (beta < 0 ? "" : " ") << beta << "\t";
            }

            std::cout << '\n';
        }
    }

    // perform independence test
    std::vector<double> LRs(9, 0);
    std::cout << "G:\t";
    for (size_t i = 0; i < m_L - 1; ++i) {
        for (const auto& j2 : DNA) {
            for (const auto& j1 : DNA) {
                LRs[i] += m_pairwise_counts[i][j1 + j2] * std::log(_pairwise_probabilities[i][j1 + j2] / (_marginal_probabilities[i][j1] * _marginal_probabilities[i + 1][j2]));
            }
        }

        LRs[i] *= 2;
        std::cout << "\t" << LRs[i] << " (p=" << 1 - gsl_cdf_chisq_P(LRs[i], 9) << ")";
    }
    std::cout << "\n\n";
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

    for (int i = 0; i < reads[0]->m_ref.m_replace_start; ++i) {
        ranks.reset();

        for (auto j : reads) {
            ranks.add_base(j->m_fullRead[i]);
        }

        std::tie<char, int>(wt, num_wt) = ranks.get_rank_base(0);

        if (ranks.get_nonambig_frac() >= minMajorFraction)
            new_consensus[i] = wt;
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