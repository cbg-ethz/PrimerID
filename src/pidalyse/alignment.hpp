#ifndef ALIGNMENT_HPP
#define ALIGNMENT_HPP

#include <string>
#include <vector>
#include <map>
#include <tuple>

#include <cstddef>
#include <cstdint>

#include "reference.hpp"
#include "proper_read.hpp"
#include "statistics.hpp"

struct raw_paired_read {
    std::string m_DNA;
    std::string m_pID;
    int m_num_inserts;
    int m_num_dels;
    int m_len_overhang;
    int m_len_pID;
    bool m_is_valid;

    raw_paired_read(const raw_paired_read&) = default; // Copy constructor
    raw_paired_read(raw_paired_read&&) = default; // Move constructor
    raw_paired_read& operator=(const raw_paired_read&)& = default; // Copy assignment operator
    raw_paired_read& operator=(raw_paired_read&&)& = default; // Move assignment operator
    ~raw_paired_read() = default;

    // practically used ctor:
    raw_paired_read(
        std::string&& DNA_,
        std::string&& pID_,
        int num_inserts_,
        int num_dels_,
        int len_overhang_,
        int len_pID_,
        bool is_valid_)
        : m_DNA(std::move(DNA_)), m_pID(std::move(pID_)), m_num_inserts(num_inserts_), m_num_dels(num_dels_), m_len_overhang(len_overhang_), m_len_pID(len_pID_), m_is_valid(is_valid_) {}
}; // store raw reads from SAM data

class alignment {
public:
    // CTOR:
    alignment(const std::string& fileName_);

    // PREPROCESSING:
    void filtering_QA();
    void remove_pID_collisions(int min_required_coverage, double min_plurality, bool report = false);

    // RT STUFF:
    hamming_return_type count_mismatches_at_locus(int locus) const;
    hamming_return_type calculate_RT_mismatches() const;
    double LogLik(double s, double r) const;

    // PCR STUFF:
    d_hamming_return_type calculate_PCR_mismatches() const;

    // STATISTICS STUFF:
    void add_to_merge_statistics(seq_statistics&) const;

    // DISPLAY I/O:
    void show_recombination_patterns() const;
    void show_clone_frequencies() const;

    // FILE I/O:
    void write_consensus_to_fasta() const;
    void write_to_fasta() const;
    void write_raw_to_fasta() const;
    void write_frequencies_to_MATLAB(std::ofstream& output, bool pID) const;

    void write_statistics_histograms() const;
    void write_statistics_to_csv() const;

private:
    // I/O PARAMETERS:
    std::string m_input_fileName;
    std::string m_fileName;
    std::string m_fileStem;

    // DATA:
    reference m_reference; // reference copy, 5VM frequencies
    std::vector<raw_paired_read> m_raw_reads; // collection of raw PAIRED reads
    std::vector<std::pair<std::string, std::vector<proper_read>>> m_raw_pID_collection; // map of pID -> PCR ensembl
    std::vector<std::pair<std::string, std::vector<proper_read*>>> m_collision_free_pID_collection; // map of pID -> collision-free PCR ensembl
    std::vector<std::pair<std::string, consensus_read>> m_consensus_pID_collection; // map of pID -> consensus sequence of collision-free PCR ensembl

    // STATISTICS:
    seq_statistics m_pid_stats;

    int m_number_singletons;
    int m_number_collisions;
    int m_number_indecisive;

    int m_mismatches;
    int m_valid_trials;
    int m_Ns;

    // PARAMETERS:
    size_t m_min_current_coverage;
    double m_min_majority_fraction;

    // PRIVATE FUNCTIONS:
    std::string call_consensus_and_remove_collisions(std::vector<proper_read>& reads, int minDisplay, std::vector<proper_read*>& filtered_reads);
    void construct_sequence(const std::string& SEQ, const std::string& CIGAR, const int POS, const reference& ref, std::string& final_raw_sequence, std::string& pID, int& num_inserts, int& num_dels, int& num_N, int& len_overhang) const;
};

class alignments {
public:
    // CTOR:
    alignments(const std::vector<std::string>& inputFiles_);

    // PREPROCESSING:
    void filtering_QA();
    void remove_pID_collisions(int min_required_coverage, double min_plurality, bool report = false);

    // RT STUFF:
    double estimate_RT_substitution_rate(bool report = false);

    double estimate_RT_recombination_rate(bool report = false);
    double estimate_RT_recombination_rate(double s, bool report = false);

    void plot_RT_recombination_LogLik(double upper = 1E-5, int n = 100) const;

    // PCR STUFF:
    double estimate_PCR_substitution_rate(bool report = false);

    // DISPLAY I/O:
    void show_recombination_patterns() const;
    void show_clone_frequencies() const;

    // FILE I/O:
    void write_all_consensus_to_fasta() const;
    void write_all_to_fasta() const;
    void write_all_raw_to_fasta() const;
    void write_frequencies_to_MATLAB() const;
    void write_all_statistics();

private:
    // DATA:
    size_t m_num_alignments; // number of input alignments
    std::vector<alignment> m_collections_alignments; // holds all alignment objects

    // STATISTICS:
    seq_statistics m_pid_stats;

    double m_RT_sub_rate;
    double m_RT_sub_rate_CI_lower;
    double m_RT_sub_rate_CI_upper;

    double m_RT_recomb_rate;
    double m_RT_recomb_rate_CI_lower;
    double m_RT_recomb_rate_CI_upper;
    double m_RT_LogLik;

    double m_PCR_sub_rate;
    double m_PCR_sub_rate_CI_lower;
    double m_PCR_sub_rate_CI_upper;

    // PARAMETERS:
    int m_min_current_coverage;
    double m_min_majority_fraction;

    bool m_is_collisions_removed;

    // PRIVATE FUNCTIONS:
    double calculate_RT_LogLik_recombination(double s, double r) const;
    double neg_LogLik_recombination(double s, double r) const;
};

#endif /* ALIGNMENT_HPP */