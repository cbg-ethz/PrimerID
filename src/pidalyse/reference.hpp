#ifndef REFERENCE_HPP
#define REFERENCE_HPP

#include <vector>
#include <string>
#include <stdexcept>

enum class referenceType {
    START_3223,
    START_3236
};

typedef std::tuple<uint64_t, uint64_t, uint64_t> hamming_return_type;
typedef std::tuple<double, uint64_t, uint64_t> d_hamming_return_type;

hamming_return_type hamming_distance(const std::string& sA, const std::string& sB);

class record {
public:
    // SEQUENCE PROPERTIES:
    std::string m_DNA;
    std::string m_name;

    std::string m_heterozygous_loci_string;

    // STATISTICS:
    int m_raw_counts;
    double m_raw_frequency;

    int m_pID_counts;
    double m_pID_frequency;

    // MEMBER FUNCTIONS:
    record(const std::string& DNA_, const std::string& name_);
};

class reference {
public:
    // GLOBAL PARAMETERS:
    static const int m_min_valid = 22;
    static const int m_max_mismatches = 0;

    // MEMBER VARIABLES:
    std::string m_referenceFile;
    std::string m_fileName;
    std::string m_fileStem;

    referenceType m_reference_variant;

    int m_K;
    int m_genome_length;

    int m_raw_total;
    int m_pID_total;

    std::vector<record> m_all_reference_strains;
    bool m_freq_initialised;
	int attempted_consensus_assignments;

    // LOCI INFORMATION:
    std::vector<int> m_included_loci;
    std::vector<int> m_not_included_loci;

    // HETEROZYGOUS INFORMATION:
    int m_num_heterozygous_loci;
    std::vector<int> m_heterozygous_loci;

    // HOMOZYGOUS INFORMATION:
    int m_num_homozygous_loci;
    std::vector<int> m_homozygous_loci;
    std::string m_homozygous_string;
    static std::vector<int> m_shared_homozygous_loci;

    // INFORMATION FOR TRUNCATING:
    int m_replace_start;
    int m_replace_length;

    int m_pID_start;
    int m_pID_length;

    int m_overhang_start;
    int m_overhang_min_length;

    // MEMBER FUNCTIONS:
    reference(const std::string& fileName_, referenceType variant_);
    reference() = default;

    void display_hamming_distance() const;
    void display_strains_hetero() const;
    void display_strains_verbose() const;
    void display_strains_abridged(const std::string& fileName) const;

    void assign_counts(const std::string& read, bool consensus = false);
    void normalise_counts();
    void reset_consensus_counts();
};

void set_reference_intersection(reference& A, reference& B);

extern reference ref3223;
extern reference ref3236;

inline bool valid_base(char base)
{
    switch (base) {
    case 'A':
    case 'C':
    case 'G':
    case 'T':
        return true;

    case '-':
    case 'N':
        return false;

    default:
        std::string errmsg("Unrecognized base: ");
        errmsg.push_back(base);
        throw std::runtime_error(errmsg);
    }
}

inline bool non_ambiguous_base(char base)
{
    switch (base) {
    case 'A':
    case 'C':
    case 'G':
    case 'T':
    case '-':
        return true;

    case 'N':
        return false;

    default:
        std::string errmsg("Unrecognized base: ");
        errmsg.push_back(base);
        throw std::runtime_error(errmsg);
    }
}

#endif /* REFERENCE_HPP */