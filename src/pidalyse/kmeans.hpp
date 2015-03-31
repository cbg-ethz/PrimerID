#ifndef _KMEANS_HPP_
#define _KMEANS_HPP_

#include <random>
#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>
Graph;

struct ranked_DNA_list
{
  char base_wt;
  int no_wt;

  char base_mt;
  int no_mt;

  int no_total;
  int no_rest;
  int no_N;

  double frac_mt_wt;

  void operator()(const std::map<char, int>& pileup)
  {
    ranking.clear();

    ranking.emplace_back('A', pileup.at('A'));
    ranking.emplace_back('C', pileup.at('C'));
    ranking.emplace_back('G', pileup.at('G'));
    ranking.emplace_back('T', pileup.at('T'));
    no_N = pileup.at('N');

    no_rest = ranking[0].second + ranking[1].second + ranking[2].second + ranking[3].second;
    no_total = no_rest + no_N;

    std::sort(ranking.begin(), ranking.end(),
              [](const std::pair<char, int>& left,
                 const std::pair<char, int>& right)
    {
      return left.second > right.second;
    });

    base_wt = ranking[0].first;
    no_wt = ranking[0].second;

    base_mt = ranking[1].first;
    no_mt = ranking[1].second;

    frac_mt_wt = static_cast<double>(no_mt + no_wt) / no_total;
  }

private:
  std::vector<std::pair<char, int>> ranking;
};

static std::default_random_engine generator;
static std::bernoulli_distribution distribution(0.5);

void call_consensus(const std::vector<proper_read>& reads,
                    const std::vector<int>& indices,
                    std::string& new_consensus)
{
  new_consensus.clear();
  std::map<char, int> current_loci;
  ranked_DNA_list ranks;

  for (int i : reads[0]._ref.heterozygous_loci)
  {
    current_loci = {
      { 'A', 0 }, { 'C', 0 }, { 'G', 0 }, { 'T', 0 }, { 'N', 0 }
    };

    for (int j : indices)
    {
      ++current_loci[reads[j]._fullRead[i]];
    }

    ranks(current_loci);

    if (ranks.no_wt > ranks.no_mt)
      new_consensus.push_back(ranks.base_wt);
    else
    {
      if (distribution(generator))
        new_consensus.push_back(ranks.base_wt);
      else
        new_consensus.push_back(ranks.base_mt);
    }
  }
}

void call_consensus(const std::vector<proper_read*>& reads,
                    std::string& new_consensus)
{
  new_consensus.clear();
  std::map<char, int> current_loci;
  ranked_DNA_list ranks;

	for (int i : reads[0]->_ref.heterozygous_loci)
  {
    current_loci = {
      { 'A', 0 }, { 'C', 0 }, { 'G', 0 }, { 'T', 0 }, { 'N', 0 }
    };

    for (auto j : reads)
    {
      ++current_loci[j->_fullRead[i]];
    }

    ranks(current_loci);

    if (ranks.no_wt > ranks.no_mt)
      new_consensus.push_back(ranks.base_wt);
    else
    {
      if (distribution(generator))
        new_consensus.push_back(ranks.base_wt);
      else
        new_consensus.push_back(ranks.base_mt);
    }
  }
}

void call_full_consensus(const std::vector<proper_read*>& reads,
                         std::string& new_consensus, double minMajorFraction)
{
  new_consensus.assign(reads[0]->_ref.genome_length, 'N');
  std::map<char, int> current_loci;
  ranked_DNA_list ranks;

  for (int i : reads[0]->_ref.included_loci)
  {
    current_loci = {
      { 'A', 0 }, { 'C', 0 }, { 'G', 0 }, { 'T', 0 }, { 'N', 0 }
    };

    for (auto j : reads)
    {
      ++current_loci[j->_fullRead[i]];
    }

    ranks(current_loci);

    if (static_cast<double>(ranks.no_wt) / ranks.no_rest > minMajorFraction)
      new_consensus[i] = ranks.base_wt;
    else
      new_consensus[i] = 'N';
  }
}

struct kmeans
{
public:
  int _n;

  int _size1;
  int _size2;

  double _frac1;
  double _frac2;

  std::string _c_mean1;
  std::string _c_mean2;

  std::string _full_consensus1;

  int _inter_cluster_ham;

  std::vector<proper_read*> _cluster1;
  std::vector<proper_read*> _cluster2;

  const static int max_hamming = 1;

  // functor
  void operator()(std::vector<proper_read>& reads, double minMajorFraction)
  {
    _n = reads.size();

    // 1.) create graph of reads
    Graph G(_n);
    int mismatches, valid_trials, Ns;

    for (int i = 0; i < _n - 1; ++i)
    {
      for (int j = i + 1; j < _n; ++j)
      {
        std::tie(mismatches, valid_trials, Ns) = reads[i].hetero_hamming_distance(reads[j]);

        if ((mismatches <= max_hamming)/* && (Ns < 3)*/)
        {
          boost::add_edge(i, j, G);
        }
      }
    }

    // 2.) determine number of connected components
    std::vector<int> component(_n);
    int num = boost::connected_components(G, &component[0]);

    std::vector<std::pair<int, std::vector<int>>> comp_sizes(num);
    for (int i = 0; i < _n; ++i)
    {
      comp_sizes[component[i]].first = component[i];
      comp_sizes[component[i]].second.emplace_back(i);
    }

    // 3.) sort connected components by size
    std::sort(comp_sizes.begin(), comp_sizes.end(),
              [](const std::pair<int, std::vector<int>>& left,
                 const std::pair<int, std::vector<int>>& right)
    {
      return left.second.size() > right.second.size();
    });

    // 4.) initialize cluster centers
    _size1 = comp_sizes[0].second.size();
    call_consensus(reads, comp_sizes[0].second, _c_mean1);

    if (comp_sizes.size() == 1)
    {
			// only one connected component, definitely no collision
      _size2 = 0;
      _inter_cluster_ham = 0;
      _c_mean2.clear();
			
			// copy all reads into pointers
			_cluster1.clear();
			for (int i = 0; i < _n; ++i)
			{
				_cluster1.push_back(&reads[i]);
			}
			
      return;
    }
    else
    {
      _size2 = comp_sizes[1].second.size();
      call_consensus(reads, comp_sizes[1].second, _c_mean2);
    }
		
		/*
		if (reads.size() == 93)
		{
			for(int i = 0; i < comp_sizes.size(); ++i)
			{
				std::cout << comp_sizes[i].first << '\n';
				
				for(const auto& j : comp_sizes[i].second)
				{
					std::cout << reads[j].heterozygous_loci_string << '\n';
				}
				
				std::cout << '\n' << '\n';
			}
		}
		*/
		
		//std::cout << "Number reads: " << reads.size() << '\n';
		//std::cout << "S1:\t" << _size1 << '\t' << _c_mean1 << '\n';
		//std::cout << "S2:\t" << _size2 << '\t' << _c_mean2 << '\n';

    // 5.) perform k-means clustering
    int dist1, dist2;
    for (int i = 0; i <= 3; ++i)
    {
      _cluster1.clear();
      _cluster2.clear();

      for (int j = 0; j < _n; ++j)
      {
        std::tie(dist1, valid_trials, Ns) = reads[j].hetero_hamming_distance(_c_mean1);
        std::tie(dist2, valid_trials, Ns) = reads[j].hetero_hamming_distance(_c_mean2);

        if (dist1 <= dist2)
          _cluster1.push_back(&reads[j]);
        else
          _cluster2.push_back(&reads[j]);
      }

      _size1 = _cluster1.size();
      _size2 = _cluster2.size();
			
			/*
			if (reads.size() == 93)
			{
				std::cout << "S1:\t" << _size1 << '\t' << _c_mean1 << '\n';
				std::cout << "S2:\t" << _size2 << '\t' << _c_mean2 << '\n';
				std::cout << '\n';
			}
			*/
			
      call_consensus(_cluster1, _c_mean1);
      call_consensus(_cluster2, _c_mean2);
    }

    _frac1 = static_cast<double>(_size1) / (_size1 + _size2);
    _frac2 = static_cast<double>(_size2) / (_size1 + _size2);

    std::tie(_inter_cluster_ham, valid_trials, Ns) = hamming_distance(_c_mean1, _c_mean2);
    call_full_consensus(_cluster1, _full_consensus1, minMajorFraction);
  }
};

#endif /* _KMEANS_HPP_ */