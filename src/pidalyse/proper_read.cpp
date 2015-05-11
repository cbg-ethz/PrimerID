#include <cmath>
#include <tuple>
#include <functional>
#include <exception>

#include "proper_read.hpp"

hamming_return_type
hamming_distance(const std::string& sA, const std::string& sB)
{
    int mismatches = 0;
    int valid_trials = 0;
    int Ns = 0;

    if (sA.length() != sB.length())
        throw std::range_error("String lengths do not match!\n");

    for (size_t i = 0; i < sA.length(); ++i) {
        if (non_ambiguous_base(sA[i]) && non_ambiguous_base(sB[i])) {
            ++valid_trials;
            mismatches += (sA[i] != sB[i]);
        }
        else {
            ++Ns;
        }
    }

    return hamming_return_type(mismatches, valid_trials, Ns);
}

// proper_read
proper_read::proper_read(const std::string& DNA_, const reference& reference_)
    : m_ref(reference_), m_fullRead(DNA_), m_num_N(0)
{
    // 1.) construct heterozygous string
    for (int i : m_ref.m_heterozygous_loci) {
        m_num_N += (!non_ambiguous_base(m_fullRead[i]));
        m_heterozygous_loci_string.push_back(m_fullRead[i]);
    }

    m_num_of_valid_heterozygous_bases = m_ref.m_num_heterozygous_loci - m_num_N;

    // 2.) find closest reference strain
    int best_index = m_ref.m_K;
    int best_ham = 100000;

    int cur_ham;
    int valid_trials;
    int Ns;

    for (int i = 0; i < m_ref.m_K; ++i) {
        std::tie(cur_ham, valid_trials, Ns) = hamming_distance(
            m_heterozygous_loci_string,
            m_ref.m_all_reference_strains[i].m_heterozygous_loci_string);

        if (cur_ham < best_ham) {
            best_ham = cur_ham;
            best_index = i;
        }
    }

    // 3.) assign RECOMBINANT if hamming distance to known is too large
    m_hamming_distance_to_best_reference = best_ham;
    if (m_hamming_distance_to_best_reference >= 2) {
        // suspect, possibly recombinant
        m_best_reference = m_ref.m_K;
    }
    else {
        // probably no recombinant
        m_best_reference = best_index;
    }
}

hamming_return_type
proper_read::hetero_hamming_distance(const std::string& other_string) const
{
    return hamming_distance(this->m_heterozygous_loci_string, other_string);
}

hamming_return_type
proper_read::hetero_hamming_distance(const proper_read& other_read) const
{
    return hetero_hamming_distance(other_read.m_heterozygous_loci_string);
}

// consensus_read
consensus_read::consensus_read(std::string&& DNA_, const reference& reference_, int multiplicity_)
    : proper_read(DNA_, reference_), m_fullConsensus(std::move(DNA_)), m_multiplicity(multiplicity_) {}

hamming_return_type consensus_read::calculate_homozygous_mismatches() const
{
    int mismatches = 0;
    int valid_trials = 0;
    int Ns = 0;

    for (int i : m_ref.m_homozygous_loci) {
        if (valid_base(m_fullConsensus[i])) {
            ++valid_trials;
            mismatches += (m_fullConsensus[i] != m_ref.m_all_reference_strains[0].m_DNA[i]);
        }
        else {
            ++Ns;
        }
    }

    return hamming_return_type(mismatches, valid_trials, Ns);
}

inline long double emission_probability(char X_i, char Z_i, long double s)
{
    if (non_ambiguous_base(X_i))
        return (X_i == Z_i ? 1 - s : s / 3);
    else
        return 1;
}

inline long double transition_probability(char Z_i, char Z_i_min_1, int n_steps,
    long double r,
    const reference& ref)
{
    long double p_norecombination = pow(1 - r, n_steps);
    long double p_ij = ref.m_all_reference_strains[Z_i].m_pID_frequency + p_norecombination * ((Z_i == Z_i_min_1) - ref.m_all_reference_strains[Z_i].m_pID_frequency);

    return p_ij;
}

double consensus_read::log_prob(double s, double r) const
{
    // performs the full HMM calculation for one sequence
    // FORWARD algorithm: P(X) = \sum_{Z} P(X, Z)
    // complexity: O(L*K^2), i.e., approx 25*L

    // probability is a function of r, i.e. the likelihood
    // of this sequence can be used to make inference on r
    std::vector<long double> f_old(m_ref.m_K), f_new(m_ref.m_K);

    // 1.) initialize
    for (int j = 0; j < m_ref.m_K; ++j)
        f_new[j] = m_ref.m_all_reference_strains[j].m_pID_frequency * emission_probability(
                                                                          m_heterozygous_loci_string[0],
                                                                          m_ref.m_all_reference_strains[j].m_heterozygous_loci_string[0], s);

    // 2.) recursion
    long double sum;
    int n_jump;

    for (int i = 1; i < m_ref.m_num_heterozygous_loci; ++i) {
        f_new.swap(f_old);
        n_jump = m_ref.m_heterozygous_loci[i] - m_ref.m_heterozygous_loci[i - 1];

        for (int j = 0; j < m_ref.m_K; ++j) {
            sum = 0;

            for (int k = 0; k < m_ref.m_K; ++k) {
                sum += transition_probability(j, k, n_jump, r, m_ref) * f_old[k];
            }

            f_new[j] = emission_probability(
                           m_heterozygous_loci_string[i],
                           m_ref.m_all_reference_strains[j].m_heterozygous_loci_string[i], s) * sum;
        }
    }

    // 3.) termination
    sum = 0;
    for (int j = 0; j < m_ref.m_K; ++j)
        sum += f_new[j];

    return log(sum);
}