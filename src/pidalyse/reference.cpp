#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include "reference.hpp"

const std::vector<int> included_loci_3223 = { 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246,
    294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494 };

const std::vector<int> included_loci_3236 = { 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246,
    312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494 };

// record
record::record(const std::string& DNA_, const std::string& name_)
    : m_DNA(DNA_), m_name(name_), m_raw_counts(0), m_raw_frequency(-1), m_pID_counts(0), m_pID_frequency(-1){};

// reference
reference::reference(const std::string& fileName_, referenceType variant_)
    : m_referenceFile(fileName_), m_reference_variant(variant_), m_raw_total(0), m_pID_total(0), m_freq_initialised(false), attempted_consensus_assignments(0), m_num_heterozygous_loci(0), m_num_homozygous_loci(0)
{
    m_fileName = boost::filesystem::path(m_referenceFile).filename().string();
    m_fileStem = boost::filesystem::path(m_referenceFile).stem().string();

    // 1.) load FASTA file
    std::ifstream input;
    input.open(m_referenceFile.c_str());

    if (input.is_open()) {
        std::string temp, id;

        while (input.good()) {
            getline(input, temp);
            temp.erase(std::remove(temp.begin(), temp.end(), '\n'), temp.end());

            if (temp.length() == 0)
                continue;

            if (temp[0] == '>') {
                // header
                id = temp.substr(1, std::string::npos);
            }
            else {
                // sequence
                m_all_reference_strains.emplace_back(record(temp, id));
            }
        }
    }
    else {
        std::cout << m_fileName << " is not readable!\n";
        exit(EXIT_FAILURE);
    }

    input.close();

    // 2.) initialise property variable
    // this depends on the locus under consideration
    m_replace_start = 495;
    if (m_reference_variant == referenceType::START_3223)
        m_replace_length = 14;
    else
        m_replace_length = 27;

    m_pID_length = 10;
    m_pID_start = m_replace_start + m_replace_length;

    m_overhang_start = m_pID_start + m_pID_length;
    m_overhang_min_length = 23;

    // 3.) set loci to include
    // hard-coded at the moment, will be replaced at a later stage
    // 0-based indexing used here, makes usage with C++ easier
    m_included_loci = (m_reference_variant == referenceType::START_3223 ? ::included_loci_3223 : ::included_loci_3236);
    if (!std::is_sorted(m_included_loci.begin(), m_included_loci.end())) {
        std::cout << "The included loci are not sorted for " << m_referenceFile << '\n';
        exit(EXIT_FAILURE);
    }

    // determine complement of included_loci
    for (int i = 0; i < m_replace_start; ++i) {
        if (!std::binary_search(m_included_loci.begin(), m_included_loci.end(), i)) {
            m_not_included_loci.push_back(i);
        }
    }

    // 4.) determine homozygous/heterozygous loci
    m_K = m_all_reference_strains.size();
    m_genome_length = m_all_reference_strains[0].m_DNA.length();

    char allele;
    bool homozygous;

    // decide whether locus i is homo or heterozygous
    for (int i : m_included_loci) {
        allele = m_all_reference_strains[0].m_DNA[i];
        homozygous = true;

        for (int j = 1; j < m_K; ++j) {
            homozygous = (allele == m_all_reference_strains[j].m_DNA[i]);
            if (!(homozygous))
                break;
        }

        if (homozygous) {
            m_homozygous_loci.push_back(i);
            m_homozygous_string.push_back(allele);
        }
        else {
            m_heterozygous_loci.push_back(i);
        }
    }

    m_num_heterozygous_loci = m_heterozygous_loci.size();
    m_num_homozygous_loci = m_homozygous_loci.size();

    // construct heterozygous strings
    for (record& i : m_all_reference_strains) {
        for (int j : m_heterozygous_loci) {
            i.m_heterozygous_loci_string.push_back(i.m_DNA[j]);
        }
    }
}

// list of all shared homozygous loci
std::vector<int> reference::m_shared_homozygous_loci;

void set_reference_intersection(reference& A, reference& B)
{
    reference::m_shared_homozygous_loci.clear();
    std::set_intersection(A.m_homozygous_loci.begin(), A.m_homozygous_loci.end(), B.m_homozygous_loci.begin(), B.m_homozygous_loci.end(), back_inserter(reference::m_shared_homozygous_loci));
}

void reference::display_strains_hetero() const
{
    std::cout << std::string(50, '=') << '\n';
    std::cout << "Reference strains" << '\n';
    std::cout << "-----------------------\n";

    std::cout << "Reference: " << m_fileName << '\n';
    std::cout << '\n';

    std::cout << "Heterozygous loci: ";

    // show positions of heterozygous loci
    for (int j : m_heterozygous_loci)
        std::cout << j << ' ';
    std::cout << '\n';

    // show heterozygous strings
    for (const record& i : m_all_reference_strains) {
        std::cout << i.m_name << '\t' << i.m_heterozygous_loci_string << '\n';
    }
    std::cout << '\n';

    std::cout << "Number included loci:     " << m_included_loci.size() << '\n';
    std::cout << "Number homozygous loci:   " << m_num_homozygous_loci << '\n';
    std::cout << "Number heterozygous loci: " << m_num_heterozygous_loci << '\n';

    std::cout << std::string(50, '=') << '\n';
}

void reference::display_strains_verbose() const
{
    const static int width = 100;

    std::cout << std::string(50, '=') << '\n';
    std::cout << "Reference strains (verbose)" << '\n';
    std::cout << "-----------------------\n";

    std::cout << "Reference: " << m_fileName << '\n';
    std::cout << '\n';

    for (int i = 0; i < m_genome_length; i += width) {
        // print 5VM DNA
        for (const record& j : m_all_reference_strains)
            std::cout << j.m_name << '\t' << j.m_DNA.substr(i, width) << '\n';

        // print annotation
        std::cout << std::string(7 - std::to_string(i + 1).length(), ' ') << i + 1;
        std::cout << '\t';
        for (int j = i; j < std::min(i + width, m_genome_length); ++j) {
            if (j < m_replace_start) {
                // in main read
                if (std::binary_search(m_heterozygous_loci.begin(), m_heterozygous_loci.end(), j)) {
                    // heterozygous locus
                    if (std::binary_search(m_included_loci.begin(), m_included_loci.end(), j)) {
                        // included locus
                        std::cout << 'X';
                    }
                    else {
                        std::cout << '*';
                    }
                }
                else {
                    // homozygous locus
                    if (std::binary_search(m_included_loci.begin(), m_included_loci.end(), j)) {
                        // included locus
                        std::cout << '.';
                    }
                    else {
                        std::cout << ' ';
                    }
                }
            }
            else {
                // remainder, i.e. right primer
                std::cout << ' ';
            }
        }

        std::cout << ' ' << std::min(i + width, m_genome_length) << "\n\n";
    }

    std::cout << "Number included loci:     " << m_included_loci.size() << '\n';
    std::cout << "Number homozygous loci:   " << m_num_homozygous_loci << '\n';
    std::cout << "Number heterozygous loci: " << m_num_heterozygous_loci << '\n';

    std::cout << std::string(50, '=') << '\n';
}

void reference::display_strains_abridged(const std::string& input_fileName) const
{
    std::cout << std::fixed << std::setprecision(2) << std::string(50, '=') << '\n';

    if (!input_fileName.empty())
        std::cout << input_fileName << ":\n";

    std::printf("%17s%22s\n", "Raw", "pID");
    for (const record& i : m_all_reference_strains) {
        std::printf("%-7s%6.2f%% (%6d)%12.2f%% (%5d)\n", i.m_name.c_str(), i.m_raw_frequency * 100, i.m_raw_counts, i.m_pID_frequency * 100, i.m_pID_counts);
    }
    std::printf("\nTOTAL:%16d%21d\n", m_raw_total, m_pID_total);
	std::printf("\nAttempted pID Assignments: %d (%.1f%%)\n", attempted_consensus_assignments, 100.0 * m_pID_total / attempted_consensus_assignments );

    std::cout << std::string(50, '=') << '\n';
}

void reference::assign_counts(const std::string& read, bool consensus)
{
	attempted_consensus_assignments += consensus;
	
    int min_dist = 1000;
    int index = m_K;

    int hamming;
    int valid_trials;

    for (int i = 0; i < m_K; ++i) {
        hamming = 0;
        valid_trials = 0;

        for (int j : m_heterozygous_loci) {
            if (valid_base(read[j])) {
                ++valid_trials;
                hamming += (read[j] != m_all_reference_strains[i].m_DNA[j]);
            }
        }

        // minimum number of valid bases to call a reference
        if (valid_trials < reference::m_min_valid)
            break;

        if (hamming < min_dist) {
            min_dist = hamming;

            if (min_dist < 2) {
                index = i;
            }
        }
    }

    if (index != m_K) {
        if (min_dist <= reference::m_max_mismatches) {
            ++(consensus ? m_all_reference_strains[index].m_pID_counts : m_all_reference_strains[index].m_raw_counts);
        }
    }
}

void reference::normalise_counts()
{
    m_raw_total = 0;
    m_pID_total = 0;

    for (const record& i : m_all_reference_strains) {
        m_raw_total += i.m_raw_counts;
        m_pID_total += i.m_pID_counts;
    }

    for (record& i : m_all_reference_strains) {
        i.m_raw_frequency = static_cast<double>(i.m_raw_counts) / m_raw_total;
        i.m_pID_frequency = static_cast<double>(i.m_pID_counts) / m_pID_total;
    }

    m_freq_initialised = true;
}

void reference::reset_consensus_counts()
{
    m_raw_total = 0;
    m_pID_total = 0;

    for (record& i : m_all_reference_strains) {
        i.m_raw_counts = 0;
        i.m_raw_frequency = -1;

        i.m_pID_counts = 0;
        i.m_pID_frequency = -1;
    }

    m_freq_initialised = false;
	attempted_consensus_assignments = 0;
}