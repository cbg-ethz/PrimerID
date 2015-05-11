#ifndef PROPER_READ_HPP
#define PROPER_READ_HPP

#include <vector>
#include <string>

#include "reference.hpp"

class proper_read {
public:
    // REFERENCE STUFF:
    const reference& m_ref;
    const std::string& m_fullRead;

    int m_best_reference;
    int m_hamming_distance_to_best_reference;

    int m_num_N;
    int m_num_of_valid_heterozygous_bases;

    // HETEROZYGOUS LOCI STRING:
    std::string m_heterozygous_loci_string;

    // MEMBER FUNCTIONS:
    proper_read(const std::string& DNA_, const reference& reference_);

    hamming_return_type hetero_hamming_distance(const std::string& other_string) const;
    hamming_return_type hetero_hamming_distance(const proper_read& other_read) const;
};

class consensus_read : public proper_read {
public:
    // CONTAINS FULL CONSENSUS SEQUENCE:
    const std::string m_fullConsensus;
    int m_multiplicity;

    // MEMBER FUNCTIONS:
    using proper_read::proper_read;
    consensus_read(std::string&& input_, const reference& reference_, int multiplicity_);

    hamming_return_type calculate_homozygous_mismatches() const;

    // s:  substitution rate
    // r: recombination rate
    double log_prob(double s, double r) const;
};

#endif /* PROPER_READ_HPP */