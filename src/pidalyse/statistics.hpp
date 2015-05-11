#ifndef STATISTICS_HPP
#define STATISTICS_HPP

#include <string>
#include <vector>
#include <valarray>
#include <tuple>
#include <map>

#include <random>
#include <utility>
#include <string>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include <boost/bimap.hpp>
#include <boost/bimap/vector_of.hpp>

#include <boost/math/distributions/binomial.hpp>

#include "proper_read.hpp"

std::string number_To_DNA(int number, int L);
int DNA_to_number(const std::string& DNA);

/* DNAVECTOR */
template <typename T = int>
class DNAvector : public std::valarray<T> {
public:
    size_t m_L;

    DNAvector(int L_)
        : std::valarray<T>(static_cast<T>(0), static_cast<int>(std::pow(4, L_))), m_L(L_) {}

    using std::valarray<T>::operator[];
    using std::valarray<T>::operator+=;
    using std::valarray<T>::operator/=;

    // indexing with DNA sequences
    T& operator[](const std::string& str)
    {
        return this->operator[](DNA_to_number(str));
    }

    const T& operator[](const std::string& str) const
    {
        return this->operator[](DNA_to_number(str));
    }
};

// fix for libc++ bug
#ifdef _LIBCPP_VERSION
template <class _Tp>
struct std::__is_val_expr<DNAvector<_Tp>> : std::true_type {
};
#endif /* _LIBCPP_VERSION */

class seq_statistics {
public:
    seq_statistics(int L_);

    void reset();
    void show_statistics() const;
    void addLengthToHistogram(int lengthpID);
    void addPrimer(const std::string& strPrimer, int replicates);
    void addPrimer_collisionFree(const std::string& strPrimer, int replicates);
    void mergestatistics(const seq_statistics& statisticsB);
    void write_to_csv(const std::string& file_stem) const;
    void write_histograms(const std::string& file_stem) const;
    void calculate_comprehensive_statistics(const std::string& strLabel) const;

private:
    size_t m_L;
    int m_unique_counts;
    int m_repl_counts;

    std::vector<DNAvector<int>> m_marginal_counts; // position-wise counts of bases in pID
    std::vector<DNAvector<int>> m_marginal_counts_repl; // position-wise counts of bases in pID times PCR multiplicity

    std::vector<DNAvector<int>> m_pairwise_counts; // position-wise counts of bases in pID
    std::vector<DNAvector<int>> m_pairwise_counts_repl; // position-wise counts of bases in pID times PCR multiplicity

    DNAvector<int> m_whole_pID_counts; // pID counts
    DNAvector<int> m_whole_pID_counts_repl; // pID counts times PCR multiplicity

    DNAvector<int> m_collision_free_whole_pID_counts; // pID counts without collisions
    DNAvector<int> m_collision_free_whole_pID_counts_repl; // pID counts without collisions times PCR multiplicity

    std::valarray<int> m_histogram_of_abundance; // histogram of abundances
    int m_max_abundance;

    std::valarray<int> m_histogram_of_pID_length; // histogram of pID lengths
    const static int m_max_len_pID = 20;
};

class prob_cycle {
public:
    double operator()(int num_mt, int num_total, int num_min, int cycle = 1) const;
    int getSize() const;

    double p_cycle(int cycle) const;
    double p_X_given_cycle_and_constraint(int num_mt, int num_total, int num_min, int cycle) const;
    double p_cycle_given_X_and_constraint(int num_mt, int num_total, int num_min, int cycle) const;

private:
    static const int cycle_depth = 20;
    mutable std::map<std::tuple<int, int, int, int>, double> m_cached_results;
};

class ranked_DNA_list {
public:
    // INITIALIZERS:
    ranked_DNA_list();
    void reset();

    // MODIFIERS:
    void add_base(char base);
    void sort();

    // ACCESSORS:
    std::pair<char, int> get_rank_base(unsigned int i = 0);
    int get_ambig_count();
    int get_total_mt_count();
    double get_ambig_frac();
    double get_nonambig_frac();

private:
    boost::bimap<boost::bimaps::set_of<char>, boost::bimaps::vector_of<int>> m_base_ranks;

    bool m_is_sorted;

    int m_num_total;
    int m_num_rest;
    int m_num_ambig;
};

std::string call_hetero_consensus(const std::vector<proper_read>& reads, const std::vector<int>& indices);
std::string call_hetero_consensus(const std::vector<proper_read*>& reads);
std::string call_full_consensus(const std::vector<proper_read*>& reads, double minMajorFraction);

struct kmeans {
public:
    // functor
    std::tuple<double, int, std::vector<proper_read*>, std::string> operator()(std::vector<proper_read>& reads, double minMajorFraction);

private:
    int m_total_size;

    int m_size1;
    int m_size2;

    double m_frac1;
    double m_frac2;

    std::string m_c_mean1;
    std::string m_c_mean2;

    std::string m_full_consensus1;

    int m_inter_cluster_ham;

    std::vector<proper_read*> m_cluster1;
    std::vector<proper_read*> m_cluster2;

    const static int max_hamming = 1;

    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
};

#endif /* STATISTICS_HPP */