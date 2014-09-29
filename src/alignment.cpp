#define _ALIGNMENT_CPP_

#include <iostream>
#include <iomanip>
#include <fstream>

#include <alignment.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/roots.hpp>

#include <prob_cycle.hpp>

// alignment
struct SAMentry
{
	std::string QNAME;
	int         FLAG;
	std::string RNAME;
	int         POS;
	int         MAPQ;
	std::string CIGAR;
	std::string RNEXT;
	int         PNEXT;
	int         TLEN;
	std::string SEQ;
	std::string QUAL;

	SAMentry () = default;

	SAMentry(
		const std::string& _QNAME,
		const std::string& _FLAG,
		const std::string& _RNAME,
		const std::string& _POS,
		const std::string& _MAPQ,
		const std::string& _CIGAR,
		const std::string& _RNEXT,
		const std::string& _PNEXT,
		const std::string& _TLEN,
		const std::string& _SEQ,
		const std::string& _QUAL
	) :
		QNAME (_QNAME),
		FLAG  (stoi(_FLAG)),
		RNAME (_RNAME),
		POS   (stoi(_POS)),
		MAPQ  (stoi(_MAPQ)),
		CIGAR (_CIGAR),
		RNEXT (_RNEXT),
		PNEXT (stoi(_PNEXT)),
		TLEN  (stoi(_TLEN)),
		SEQ   (_SEQ),
		QUAL  (_SEQ) {}
};

typedef std::pair<std::string, std::pair<SAMentry, SAMentry> > raw_sequence_pair;
typedef std::map<std::string, std::pair<SAMentry, SAMentry> > raw_sequence_pairs;

bool construct_sequence(const SAMentry& READ, std::string& DNA, std::string& PrimerID, const reference& ref)
{
	const std::string& read_sequence = READ.SEQ;
	const std::string& CIGAR = READ.CIGAR;
	
	int curPos = 0;
	int curOffset = 0;
	int offsetRef = READ.POS;
	int curLen;
	int lengthPID = 0;
	
	bool valid = true;
	
	DNA.clear();
	PrimerID.clear();
	
	// 1.) tokenize CIGAR
	for (int i = 0; i < CIGAR.length(); ++i)
	{
		if (isalpha(CIGAR[i]))
		{
			curLen = stoi(CIGAR.substr(curPos, i-curPos));
			curPos = i+1;
			
			switch (CIGAR[i]) {
				case 'M':
					DNA.append(read_sequence, curOffset, curLen);
					curOffset += curLen;
					offsetRef += curLen;
					break;
					
				case 'D':
					DNA.append(curLen, 'N');
					offsetRef += curLen;
					break;
					
				case 'I':
					if (offsetRef == ref.PID_start)
					{
						valid = false;
						lengthPID += curLen;
					}
					
					curOffset += curLen;
					break;
					
				case 'S':
					curOffset += curLen;
					break;
			}
		}
	}
	
	// 2.) check all strings
	if (READ.POS > 260)
	{
		// right sequence
		if (offsetRef - ref.overhang_start >= ref.overhang_min_length)
		{
			// overhang is good
			PrimerID = DNA.substr(ref.PID_start - READ.POS, ref.PID_length);
			lengthPID += std::count(PrimerID.begin(), PrimerID.end(), 'N');
			
			if (lengthPID > 0)
			{
				// contains some N's, PrimerID shorter than 10 nucleotides
				valid = false;
			}
		}
		else
		{
			// overhang too short
			valid = false;
		}
		
		if (READ.POS + DNA.length() > ref.replace_start)
			DNA.erase(ref.replace_start - READ.POS);
		
		DNA.append(ref.PID_start - READ.POS - DNA.length(), 'N');
	}
	
	return valid;
}

void merge_reads(const std::string& left_Read, const std::string& right_Read, int left_start, int right_start, const reference& ref, std::string& new_final_read)
{
	int global_start, global_end;
	
	global_start = left_start;
	global_end   = right_start + right_Read.length();
	
	new_final_read.clear();
	new_final_read.append(std::max(left_start-1, 0), 'N');
	new_final_read.append(left_Read);
	new_final_read.append(std::max(static_cast<int>(right_start-1 -new_final_read.length()), 0), 'N');
	new_final_read.append(right_Read);
	new_final_read.append(std::max(static_cast<int>(ref.genome_length -new_final_read.length()), 0), 'N');
}

alignment::alignment(const std::string& fileName) : input_fileName(fileName)
{
	int slash_pos = input_fileName.find_last_of('/');
	slash_pos = (slash_pos == std::string::npos ? 0 : slash_pos + 1);
	just_fileName = input_fileName.substr(slash_pos);
	
	if (input_fileName.find("3223") != std::string::npos)
		_reference = ref3223;
	else
		_reference = ref3236;	
	
	// 1.) read raw SAM input
	raw_sequence_pairs sam_data;
	
	std::ifstream input;
	input.open(fileName.c_str());
	raw_sequence_pairs::iterator it;
	
	if (input.is_open())
	{
		std::string temp, QNAME;
		std::vector<std::string> SplitVec;
		int line_no = 0;
		
		while (input.good())
		{
			++line_no;
			if (line_no % 100000 == 0)
				std::cout << "Loaded " << line_no << " lines\n";
			
			getline(input, temp);
			temp.erase(std::remove(temp.begin(), temp.end(), '\n'), temp.end());
			
			if ((temp.length() == 0) || (temp[0] == '@'))
				continue;
			
			boost::split(SplitVec, temp, boost::is_any_of("\t"), boost::token_compress_on);
			
			QNAME = SplitVec[0].substr(0, SplitVec[0].find(' '));
			it = sam_data.find(QNAME);
			if (it == sam_data.end())
			{
				// not found
				sam_data.insert(raw_sequence_pair(QNAME, std::pair<SAMentry, SAMentry>(
					SAMentry(
						QNAME,
						SplitVec[1],
						SplitVec[2],
						SplitVec[3],
						SplitVec[4],
						SplitVec[5],
						SplitVec[6],
						SplitVec[7],
						SplitVec[8],
						SplitVec[9],
						SplitVec[10]
						),
					SAMentry()
					)
				));
			}
			else
			{
				// found
				it->second.second = SAMentry(
						QNAME,
						SplitVec[1],
						SplitVec[2],
						SplitVec[3],
						SplitVec[4],
						SplitVec[5],
						SplitVec[6],
						SplitVec[7],
						SplitVec[8],
						SplitVec[9],
						SplitVec[10]
					);
				
				// swap elements such that first element is left sequence and second element is right sequence
				if (it->second.first.POS > 260)
				{
					std::swap(it->second.first, it->second.second);
				}
			}
		}
	}
	else
	{
		std::cout << fileName << " is not readable!\n";
		exit(EXIT_FAILURE);
	}
	input.close();
		
	// 2.) convert raw SAM data to proper sequencing reads
	int accepted = 0;
	
	std::string primerID;
	std::string leftRead, rightRead;
	std::string temp_final_read;
	std::pair<std::map<std::string, std::vector<proper_read> >::iterator, bool> insert_it;
	bool valid;
	
	for (raw_sequence_pairs::const_iterator i = sam_data.begin(); i != sam_data.end(); ++i)
	{	
		// 1.) check whether reads have a matching mate
		if ((i->second.second.SEQ.length() == 0) || (i->second.first.SEQ.length() == 0))
			continue;
		
		// 2.) parse left and right read and check whether primerID is valid
		        construct_sequence(i->second.first, leftRead, primerID, _reference);
		valid = construct_sequence(i->second.second, rightRead, primerID, _reference);

		// 3.) merge left and right read
		merge_reads(leftRead, rightRead, i->second.first.POS, i->second.second.POS, _reference, temp_final_read);
		
		_reference.assign_counts(temp_final_read);
		
		if (!(valid))
			continue;
		
		// 4.) add to map
		insert_it = raw_primerID_map.emplace(primerID, std::vector<proper_read>());
		insert_it.first->second.emplace_back(temp_final_read, _reference);		
		++accepted;
	}
	
	std::cout << input_fileName << ": " << sam_data.size() << " pairs, " << accepted << " passed initial QA (" << std::fixed << std::setprecision(1) << static_cast<double>(accepted)/sam_data.size()*100 << "%)\n";
}

std::string call_consensus(const std::vector<proper_read>& strings, int minC, double ambigFrac)
{
	std::map<char, int> current_loci;
	std::string consensus;
	
	for (int i = 0; i < strings[0].fullRead.length(); ++i)
	{
		current_loci = {{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}, {'N', 0}};
		
		for (int j = 0; j < strings.size(); ++j)
		{
			++current_loci[strings[j].fullRead[i]];
		}
		
		std::map<char, int>::const_iterator it = std::max_element(current_loci.begin(), current_loci.end(),
			[](const std::pair<char, int>& p1, const std::pair<char, int>& p2)
			{
				return p1.second < p2.second;
			}
		);
		
		if (it->first == 'N')
			consensus.push_back('N');
		else
		{
			if (strings.size()-current_loci['N'] >= minC)
			{
				if (static_cast<double>(it->second)/strings.size() >= ambigFrac)
					consensus.push_back(it->first);
				else
					consensus.push_back('N');
			}
			else
				consensus.push_back('N');
		}
	}
	
	return consensus;
}

/*
void __verbose(const std::string& PrimerID, const std::vector<proper_read>& reads, int valid_reads)
{
	std::cout << PrimerID << " (" << reads.size() << ", valid: " << valid_reads << ") - " << std::fixed << std::setprecision(3) << mismatch_rate << "\n";
	for (int i = 0; i < reads.size(); ++i)
	{
		std::cout << reads[i].heterozygous_loci_string;
		std::cout << '\t';
		if (reads[i].best_reference == _reference.K)
			std::cout << "RECOMBINANT";
		else
			std::cout << _reference.all_reference_strains[reads[i].best_reference].name;
	
		std::cout << " (" << reads[i].hamming_distance_to_best_reference << " / " << reads[i].no_of_valid_heterozygous_bases << ")\n";
	}
	std::cout << "\n";
}
*/

std::string alignment::call_consensus_and_remove_collisions(std::vector<proper_read>& reads, int minDisplay, const std::string& PrimerID)
{
	// do pileup to call most represented strain
	int valid_reads = 0;
	std::vector<int> Ref_pileup(_reference.K + 1, 0);
	for (int i = 0; i < reads.size(); ++i)
	{
		if (reads[i].no_of_valid_heterozygous_bases >= reference::min_valid)
		{
			++valid_reads;
			
			if (reads[i].hamming_distance_to_best_reference >= 2)
			{
				// suspect, possibly recombinant
				reads[i].best_reference = _reference.K;
			}
			
			++Ref_pileup[reads[i].best_reference];
		}
	}
	
	// determine most represented virus
	int max = 0;
	int index_max = -1;
	
	for (int k = 0; k < _reference.K+1; ++k)
	{
		if (Ref_pileup[k] > max)
		{
			max = Ref_pileup[k];
			index_max = k;
		}
	}
	
	bool collision = false;
	double match_rate = 0;
	if (index_max != _reference.K)
	{
		// major variant is not a recombinant
		if (valid_reads >= _min_current_coverage)
			collision = (static_cast<double>(max)/valid_reads < _min_majority_fraction);
	}
	else
	{
		// major variant is likely a recombinant
		// have to check whether recombinant is not a collision, i.e., a pure RT recombinant
		int no_invalid_pairs = 0;
		
		for (int i = 0; i < reads.size()-1; ++i)
		{
			if (reads[i].no_of_valid_heterozygous_bases >= reference::min_valid)
			{
				for (int j = i+1; j < reads.size(); ++j)
				{
					if (reads[j].no_of_valid_heterozygous_bases >= reference::min_valid)
					{
						no_invalid_pairs += (hamming_distance(reads[i].heterozygous_loci_string, reads[j].heterozygous_loci_string) >= 2);
					}
				}
			}
		}
		match_rate = 1.0 - 2*static_cast<double>(no_invalid_pairs)/(valid_reads*(valid_reads-1));
		
		collision = (match_rate < _min_majority_fraction);
	}
	
	std::string new_consensus;
	if (collision == false)
	{
		if (valid_reads >= _min_current_coverage)
		{
			//__verbose(const std::string& PrimerID, const std::vector<proper_read>& reads, int valid_reads)
			
			// can call consensus safely
			return call_consensus(reads, _min_current_coverage, _min_majority_fraction);
		}
		else
		{
			// too few sequences to call consensus sequence, indecisive
			return std::string("1");
		}
	}
	else
	{
		// collision
		if ((!(PrimerID.empty())) && (reads.size() > minDisplay) && (minDisplay)/* && (index_max == _reference.K)*/)
		{
			//__verbose(const std::string& PrimerID, const std::vector<proper_read>& reads, int valid_reads)
		}
		
		return std::string("0");
	}
}

void alignment::remove_primerID_collisions(int minC, double minMajorFraction, bool report, int minDisplay)
{
	_min_current_coverage = minC;
	_min_majority_fraction = minMajorFraction;
	
	std::string consensus;
	number_collisions = 0;
	number_indecisive = 0;
	number_singletons = 0;
	
	collision_free_primerID_map.clear();
	consensus_primerID_map.clear();
	
	std::pair<std::map<std::string, consensus_read>::iterator, bool> it;
	
	for (std::map<std::string, std::vector<proper_read> >::iterator i = raw_primerID_map.begin(); i != raw_primerID_map.end(); ++i)
	{
		if (i->second.size() >= minC)
		{
			consensus = call_consensus_and_remove_collisions(i->second, minDisplay, i->first);
			
			switch(consensus[0])
			{
				case '0':
					// collision
					++number_collisions;
					break;
					
				case '1':
					// indecisive
					++number_indecisive;
					break;
					
				default:
					// OK, proper consensus
					++number_singletons;
					
					collision_free_primerID_map.emplace(i->first, &(i->second));
					it = consensus_primerID_map.emplace(std::piecewise_construct,
						std::forward_as_tuple(i->first),
						std::forward_as_tuple(consensus, _reference, i->second.size())
					);
						
					_reference.assign_counts(consensus, true);
					break;
			}
		}
	}
	
	_reference.normalise_counts();
	
	if (report)
	{
		std::cout << "Indecisive PrimerIDs: " << number_indecisive << '\n';
		std::cout << " Collision PrimerIDs: " << number_collisions << " (" << std::fixed << std::setprecision(1) << static_cast<double>(number_collisions)/(number_singletons+number_collisions)*100 << "%)\n";
		std::cout << "     Total PrimerIDs: " << number_singletons+number_collisions << '\n';
	}
}

void alignment::show_primerIDs_with_min_coverage(int minC) const
{
	for (std::map<std::string, std::vector<proper_read> >::const_iterator i = raw_primerID_map.begin(); i != raw_primerID_map.end(); ++i)
	{
		if (i->second.size() >= minC)
		{
			std::cout << i->first << '\n' << std::string(10, '-') << '\n';
			for (std::vector<proper_read>::const_iterator j = i->second.begin(); j != i->second.end(); ++j)
				std::cout << j->heterozygous_loci_string << '\n';
			std::cout << '\n';
		}
	}
}


std::tuple<uint64_t, uint64_t, uint64_t> alignment::calculate_RT_mismatches() const
{
	uint64_t mismatches = 0;
	uint64_t trials = 0;
	uint64_t numberN = 0;
	
	for (std::map<std::string, consensus_read>::const_iterator i = consensus_primerID_map.begin(); i != consensus_primerID_map.end(); ++i)
	{
		for (int j = 0; j < _reference.no_homozygous_loci; ++j)
		{
			if (i->second.homozygous_loci_string[j] != 'N')
			{
				++trials;
				mismatches += (i->second.homozygous_loci_string[j] != _reference.homozygous_string[j]);
			}
			else
			{
				++numberN;
			}
		}
	}
	
	return std::tuple<uint64_t, uint64_t, uint64_t>(mismatches, trials, numberN);
}


long double alignment::LogLik(long double s, long double r) const
{
	long double total_likelihood = 0;
	long double temp;
	
	if (_reference.freq_initialised == false)
	{
		std::cout << "Cannot calculate Log Likelihood without estimated reference frequencies!\n";
		exit(EXIT_FAILURE);
	}
	
	std::map<std::string, long double> string_to_prob_cache;
	std::map<std::string, long double>::const_iterator it;
	
	//int no_unique = 0;
	//int no_total = 0;
	
	for(std::map<std::string, consensus_read>::const_iterator i = consensus_primerID_map.begin(); i != consensus_primerID_map.end(); ++i)
	{
		//++no_total;		
		it = string_to_prob_cache.find(i->second.heterozygous_loci_string);
		if (it == string_to_prob_cache.end())
		{
			// new sequence
			temp = i->second.log_prob(s, r);
			string_to_prob_cache.emplace(i->second.heterozygous_loci_string, temp);
			total_likelihood += temp;
			//++no_unique;
			
			//std::cout << temp << '\n';
		}
		else
		{
			total_likelihood += it->second;
		}
	}
	
	//std::cout << "Unique: " << no_unique << "\tTotal: " << no_total << '\n';
	return total_likelihood;
}


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
	
	void operator() (const std::map<char, int>& pileup)
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
			[](const std::pair<char, int>& a, const std::pair<char, int>& b)
			{
				return a.second >= b.second;
			}
		);
		
		base_wt = ranking[0].first;
		  no_wt = ranking[0].second;
		
		base_mt = ranking[1].first;
		  no_mt = ranking[1].second;
		
		frac_mt_wt = static_cast<double>(no_mt+no_wt) / no_total;
		
		/*
		for(int i = 0; i < 4; ++i)
		{
			std::cout << ranking[i].first << ':' << ranking[i].second << '\n';
		}
		
		std::cout << '\n' << "N:" << no_N << '\n';
		std::cout << "Frac: " << frac_mt_wt << '\n';
		*/
	}
	
	private:
		std::vector<std::pair<char, int> > ranking;
};

std::tuple<double, uint64_t, uint64_t> alignment::calculate_PCR_mismatches() const
{
	double mismatches = 0;
	uint64_t trials = 0;
	uint64_t numberN = 0;
	
	std::map<char, int> current_loci;
	
	bool PCRerror;
	int where;
		
	int cov;
	//double fracMajor;
	ranked_DNA_list ranker;
	double p;
	
	for (std::map<std::string, std::vector<proper_read>* >::const_iterator i = collision_free_primerID_map.begin(); i != collision_free_primerID_map.end(); ++i)
	{
		cov = i->second->size();
		PCRerror = false;
		
		for (int j = 0; j < _reference.no_homozygous_loci; ++j)
		{
			// 1.) perform pileup
			current_loci = {{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}, {'N', 0}};
			
			for (int k = 0; k < cov; ++k)
			{
				++current_loci[(*(i->second))[k].homozygous_loci_string[j]];
			}
			
			// 2.) perform ranking
			ranker(current_loci);
			
			if (ranker.no_N >= 2)
				continue;
			
			++trials;
			
			if (ranker.no_mt <= 2)
				continue;

			p = PCR_prob(ranker.no_mt, ranker.no_mt + ranker.no_wt, 2);
			mismatches += p;

			if (p > 0.995)
			{
				++mismatches;
				
				PCRerror = true;
				where = j;
			}				
		}
		
		/*
		if ((PCRerror) && (cov > 10))
		{
			std::cout << i->first << " (" << where << ")\n";
			
			for (std::vector<proper_read>::const_iterator x = i->second->begin(); x != i->second->end(); ++x)
			{
				std::cout << x->homozygous_loci_string << '\n';
			}
			std::cout << '\n';
		}
		*/
		
	}
	
	return std::tuple<double, uint64_t, uint64_t>(mismatches, trials, numberN);
}

void alignment::calculate_effective_RNA_number() const
{
	int eff_N = raw_primerID_map.size() + round(0.5 * number_singletons * number_singletons / number_collisions);
	
	std::cout << "        Number RNAs: " << eff_N;
}

void alignment::display_raw_and_primerID_counts() const
{
	for (int i = 0; i < _reference.K; ++i)
	{
		std::cout << _reference.all_reference_strains[i].name << '\t';
		std::cout << _reference.all_reference_strains[i].PID_counts << '\t';
		std::cout << _reference.all_reference_strains[i].raw_counts << '\n';
	}
	std::cout << "-----------------------\n";
	std::cout << '\t' << _reference.PID_total << '\t' << _reference.raw_total;
}

void alignment::write_consensus_to_fasta() const
{
	int startPos = input_fileName.find_last_of('/')+1;
	int endPos = input_fileName.find_last_of(".sam")-3;
	
	std::string output_fileName = input_fileName.substr(startPos, endPos-startPos) + "_cons.fasta";
	
	std::ofstream output(output_fileName.c_str());
	for (std::map<std::string, consensus_read>::const_iterator i = consensus_primerID_map.begin(); i != consensus_primerID_map.end(); ++i)
	{
		output << '>' << i->first << '_' << i->second.multiplicity << '\n' << i->second.fullRead << '\n';
	}
	output.close();
}

void alignment::write_to_fasta() const
{
	int startPos = input_fileName.find_last_of('/')+1;
	int endPos = input_fileName.find_last_of(".sam")-3;
	
	std::string output_fileName = input_fileName.substr(startPos, endPos-startPos) + ".fasta";
	
	std::ofstream output(output_fileName.c_str());
	for (std::map<std::string, std::vector<proper_read>* >::const_iterator i = collision_free_primerID_map.begin(); i != collision_free_primerID_map.end(); ++i)
	{
		int k = 0;
		for (std::vector<proper_read>::const_iterator j = i->second->begin(); j != i->second->end(); ++j)
		{
			++k;
			output << '>' << i->first << '_' << k << '\n' << j->fullRead << '\n';
		}
	}
	output.close();
}


// alignments
alignments::alignments(const std::vector<std::string>& inputFiles) : 
	_min_current_coverage(-1),
	_min_majority_fraction(-1),
	is_collisions_removed(false)
{
	for (std::vector<std::string>::const_iterator i = inputFiles.begin(); i != inputFiles.end(); ++i)
	{
		std::cout << "Loading " << *i << '\n';
		collections_alignments.emplace_back(*i);
	}
}

void alignments::remove_primerID_collisions(int minC, double minMajorFraction, bool printOut)
{
	_min_current_coverage = minC;
	_min_majority_fraction = minMajorFraction;
	
	for (std::vector<alignment>::iterator i = collections_alignments.begin(); i != collections_alignments.end(); ++i)
	{
		i->remove_primerID_collisions(_min_current_coverage, _min_majority_fraction, printOut);
	}
}

/* *** RT stuff *** */
double alignments::estimate_RT_substitution_rate(bool report)
{
	double mismatches = 0;
	uint64_t trials = 0;
	uint64_t numberN = 0;
	
	std::tuple<uint64_t, uint64_t, uint64_t> result;
	
	for (std::vector<alignment>::const_iterator i = collections_alignments.begin(); i != collections_alignments.end(); ++i)
	{
		result = i->calculate_RT_mismatches();
		
		mismatches += std::get<0>(result);
		trials += std::get<1>(result);
		numberN += std::get<2>(result);
	}
	
	double sub_rate = mismatches/trials;
	
	if (report)
	{
		double CI_lower = boost::math::binomial_distribution<>::find_lower_bound_on_p(trials, mismatches, 0.025);
		double CI_upper = boost::math::binomial_distribution<>::find_upper_bound_on_p(trials, mismatches, 0.025);
	
		std::cout << std::string(50, '=') << '\n';
		std::cout << "RT substitution rate" << '\n';
		std::cout << "-----------------------\n";
		std::cout << '\n';
	
		std::cout << "           mt bases: " << static_cast<int>(mismatches) << '\n';
		std::cout << "        Total bases: " << trials << '\n';
		std::cout << "          'N' bases: " << numberN << '\n';
		std::cout << '\n';
	
		std::cout << "   Est. RT sub rate: " << std::scientific << std::setprecision(2) << sub_rate << '\n';
		std::cout << "             95% CI: [" << CI_lower << ", " << CI_upper << "]\n";
	
		std::cout << std::string(50, '=') << '\n';
	}
	
	_RT_sub_rate = sub_rate;
	return sub_rate;
}

long double alignments::calculate_RT_LogLik_recombination(double s, double r) const
{
	long double total_likelihood = 0;
	
	for (std::vector<alignment>::const_iterator i = collections_alignments.begin(); i != collections_alignments.end(); ++i)
		total_likelihood += i->LogLik(s, r);
	
	return total_likelihood;
}

long double alignments::neg_LogLik_recombination(double s, double r) const
{
	return -calculate_RT_LogLik_recombination(s, r);
}

double alignments::estimate_RT_recombination_rate(bool report)
{
	return estimate_RT_recombination_rate(_RT_sub_rate, report);
}

double alignments::estimate_RT_recombination_rate(double s, bool report)
{	
	std::pair<long double, long double> result;
	
	// 1.) calculate MLE
	result = boost::math::tools::brent_find_minima([&](double r) -> long double {return this->neg_LogLik_recombination(s, r);},
	1E-10, 0.9, 1000);
	
	long double recomb_rate = result.first;
	long double LogLik = -result.second;
	
	if (report)
	{
		// 2.) calculate 95% CI
		uintmax_t max_iter = 1000;
		boost::math::tools::eps_tolerance<long double> tol(10000);
		
		// lower CI boundary
		max_iter = 1000;
		result = boost::math::tools::toms748_solve([&](double r) -> long double {return this->calculate_RT_LogLik_recombination(s, r) - LogLik + 1.920729;},
			recomb_rate/10,
			recomb_rate,
			tol,
			max_iter);
	
		double CI_lower = (result.first + result.second) / 2;
	
		// upper CI boundary
		max_iter = 1000;
		result = boost::math::tools::toms748_solve([&](double r) -> long double {return this->calculate_RT_LogLik_recombination(s, r) - LogLik + 1.920729;},
			recomb_rate,
			recomb_rate*10,
			tol,
			max_iter);
	
		double CI_upper = (result.first + result.second) / 2;
	
		// 3.) display	
		std::cout << std::string(50, '=') << '\n';
		std::cout << "RT recombination rate" << '\n';
		std::cout << "-----------------------\n";
		std::cout << '\n';
	
		std::cout << "Est. RT recomb rate: " << std::scientific << std::setprecision(2) << recomb_rate << '\n';
		std::cout << "     Log-Likelihood: " << std::fixed << std::setprecision(2) << LogLik << '\n';
		std::cout << "             95% CI: [" << std::scientific << std::setprecision(2) << CI_lower << ", " << CI_upper << "]\n";
	
		std::cout << std::string(50, '=') << '\n';
	}
	
	_RT_recomb_rate = recomb_rate;
	return recomb_rate;
}

/* *** PCR stuff *** */
double alignments::estimate_PCR_substitution_rate(bool report) const
{
	uint64_t mismatches = 0;
	uint64_t trials = 0;
	uint64_t numberN = 0;
	
	std::tuple<uint64_t, uint64_t, uint64_t> result;
	
	for (std::vector<alignment>::const_iterator i = collections_alignments.begin(); i != collections_alignments.end(); ++i)
	{
		result = i->calculate_PCR_mismatches();
		
		mismatches += std::get<0>(result);
		trials += std::get<1>(result);
		numberN += std::get<2>(result);
	}
	
	double sub_rate = static_cast<double>(mismatches)/trials;
	if (report)
	{
		double CI_lower = boost::math::binomial_distribution<>::find_lower_bound_on_p(trials, mismatches, 0.025);
		double CI_upper = boost::math::binomial_distribution<>::find_upper_bound_on_p(trials, mismatches, 0.025);
	
		std::cout << std::string(50, '=') << '\n';
		std::cout << "PCR substitution rate" << '\n';
		std::cout << "-----------------------\n";
		std::cout << '\n';
	
		std::cout << "           mt bases: " << mismatches << '\n';
		std::cout << "        Total bases: " << trials << '\n';
		std::cout << "          'N' bases: " << numberN << '\n';
		std::cout << '\n';
	
		std::cout << "  Est. PCR sub rate: " << std::scientific << std::setprecision(2) << sub_rate << '\n';
		std::cout << "             95% CI: [" << CI_lower << ", " << CI_upper << "]\n";
	
		std::cout << std::string(50, '=') << '\n';
	}
	
	return sub_rate;
}

/* *** variance/overdispersion stuff *** */
void alignments::display_raw_and_primerID_counts() const
{
	std::cout << std::string(50, '=') << '\n';
	std::cout << "Overdispersion/Variance" << '\n';
	std::cout << "-----------------------";
	
	for (std::vector<alignment>::const_iterator i = collections_alignments.begin(); i != collections_alignments.end(); ++i)
	{
		std::cout << '\n' << i->just_fileName << '\n';
		i->display_raw_and_primerID_counts();
		std::cout << '\n';
	}
	
	std::cout << std::string(50, '=') << '\n';
}

/* *** chao estimator *** */
void alignments::estimate_effective_RNA_number(bool report) const
{
	if (report)
	{
		std::cout << std::string(50, '=') << '\n';
		std::cout << "RNA numbers" << '\n';
		std::cout << "-----------------------";
	
		for (std::vector<alignment>::const_iterator i = collections_alignments.begin(); i != collections_alignments.end(); ++i)
		{
			std::cout << '\n' << i->just_fileName << '\n';
			i->calculate_effective_RNA_number();
			std::cout << '\n';
		}
	
		std::cout << std::string(50, '=') << '\n';
	}
}

/* *** diagnostics *** */
void alignments::show_primerIDs_with_min_coverage(int minC) const
{
	for (std::vector<alignment>::const_iterator i = collections_alignments.begin(); i != collections_alignments.end(); ++i)
	{
		std::cout << i->input_fileName << '\n';
		
		i->show_primerIDs_with_min_coverage(minC);
	}
}

void alignments::write_all_consensus_to_fasta() const
{
	for (std::vector<alignment>::const_iterator i = collections_alignments.begin(); i != collections_alignments.end(); ++i)
	{
		i->write_consensus_to_fasta();
	}
}

void alignments::write_all_to_fasta() const
{
	for (std::vector<alignment>::const_iterator i = collections_alignments.begin(); i != collections_alignments.end(); ++i)
	{
		i->write_to_fasta();
	}
}