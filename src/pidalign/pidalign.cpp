/**
 * Copyright (c) 2014-2015 Jochen Singer, David Seifert
 *
 * This file is part of PIDalign
 *
 * PIDalign is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or any later version.
 *
 * PIDalign is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * PIDalign. If not, see <http://www.gnu.org/licenses/>.
 */

#include <boost/progress.hpp>
#include <src/pidalign/threadpool/threadpool.hpp>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <seqan/modifier.h>

#include "ednafull.hpp"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <thread>
#include <chrono>

using namespace seqan;
typedef String<Iupac> RefSeq_t;

uint32_t total;
std::atomic_uint_fast32_t progress(0);

// stats
std::atomic_uint_fast32_t discard_score(0);
std::atomic_uint_fast32_t discard_del(0);
uint32_t discard_bad_mate(0);
uint32_t discard_singles(0);
uint32_t discard_tlen(0);
uint32_t accepted(0);

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.

struct AppOptions {
    StringSet<CharString> patternFileNames;
    CharString refFileName;
    CharString outputFileName;
    uint32_t numThreads;
    uint32_t blockSize;
    int minScore;
    uint32_t numDeletion;
    uint32_t verbosity_level;
    bool Phred;
    bool quiet;
    int32_t minimum_tlen;

    AppOptions()
        : numThreads(std::thread::hardware_concurrency()), blockSize(500), minScore(-1000), numDeletion(4), verbosity_level(1), Phred(false), quiet(false), minimum_tlen(520)
    {
    }
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(AppOptions& options, int argc, char const** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("pidalign");
    // Set short description, version, and date.
    setShortDescription(parser, "Needleman-Wunsch (Gotoh variant) aligner with output in SAM");
    setVersion(parser, "0.1");
    setDate(parser, "May 2015");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-ref\\fP <reference.fasta> \\fB-r1\\fP <R1.fastq> \\fB-r2\\fP <R2.fastq> \\fB-o\\fP <out.sam>");
    addDescription(parser, "pidalign performs a full-exhaustive Needleman-Wunsch alignment with affine gap costs. For performance reasons, pidalign is multithreaded, and it defaults to using the number of logical cores as number of threads.");

    addOption(parser, ArgParseOption("ref", "referenceFileName", "Name of the reference file.", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "ref");
    setValidValues(parser, "ref", "fasta fa fna");

    addOption(parser, ArgParseOption("r1", "reads1FileName", "Name of the file containing the first mate pair reads.", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "r1");
    setValidValues(parser, "r1", "fastq fq fnq");

    addOption(parser, ArgParseOption("r2", "reads2FileName", "Name of the file containing the second mate pair reads.", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "r2");
    setValidValues(parser, "r2", "fastq fq fnq");

    addOption(parser, ArgParseOption("o", "outputFileName", "Name of the output file.", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "o");
    setValidValues(parser, "o", "sam");

    addOption(parser, ArgParseOption("t", "threads", "The number of threads to be used.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("b", "blockSize", "Block size of per thread processing.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("s", "minScore", "MinScore of alignment in order to be considered.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("d", "deletionCutOff", "Maximum # consecutive deletions/gaps that are allowed.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("T", "minTlen", "Minimum fragment size (TLEN) of paired-end reads.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("v", "verbosityLevel", "Level of detail to print, with 0 being silent, 1 being the default and 2 displaying a progress bar.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("Q", "phred", "Include Phred scores in final alignment."));
    addOption(parser, ArgParseOption("q", "quiet", "Do not print to stdout, equivalent to -v 0."));

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBpidalign\\fP \\fB-ref\\fP \\fIreference.fasta\\fP\\fP \\fB-r1\\fP \\fP\\fIR1.fastq\\fP \\fB-r2\\fP \\fP\\fIR2.fastq\\fP \\fB-o\\fP \\fIoutput.sam\\fP",
        "Use \\fIreference.fasta\\fP to align \\fP\\fIR1.fastq\\fP and \\fP\\fIR2.fastq\\fP with the SAM alignment in \\fIoutput.sam\\fP.");

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    resize(options.patternFileNames, 2);
    getOptionValue(options.patternFileNames[0], parser, "r1");
    getOptionValue(options.patternFileNames[1], parser, "r2");
    getOptionValue(options.refFileName, parser, "ref");
    getOptionValue(options.outputFileName, parser, "o");
    getOptionValue(options.numThreads, parser, "t");
    getOptionValue(options.blockSize, parser, "b");
    getOptionValue(options.minScore, parser, "s");
    getOptionValue(options.numDeletion, parser, "d");
    getOptionValue(options.minimum_tlen, parser, "T");
    getOptionValue(options.verbosity_level, parser, "v");
    getOptionValue(options.Phred, parser, "Q");
    getOptionValue(options.quiet, parser, "q");

    return ArgumentParser::PARSE_OK;
}

void storeRecord(String<BamAlignmentRecord>& bamRecords,
    uint16_t flag,
    Align<RefSeq_t>& align,
    StringSet<CharString> const& readIDs,
    StringSet<CharString> const& readQUALs,
    uint32_t refId,
    int score,
    uint32_t bamIndex,
    uint32_t numDeletion,
    bool phred)
{
    int clipBegin = row(align, 1)._array[0];
    int clipEnd = length(row(align, 1)) - row(align, 1)._array[length(row(align, 1)._array) - 1];
    setClippedBeginPosition(row(align, 0), clipBegin);
    setClippedBeginPosition(row(align, 1), clipBegin);
    setClippedEndPosition(row(align, 0), clipEnd);
    setClippedEndPosition(row(align, 1), clipEnd);

    getCigarString(bamRecords[bamIndex].cigar, row(align, 0), row(align, 1), 1000); // CIGAR (6)
    String<CigarElement<>> const& cigarRef = bamRecords[bamIndex].cigar;
    bool deletion_too_long = false;
    for (const auto& i : cigarRef) {
        if ((i.operation == 'D') && (i.count >= numDeletion)) {
            deletion_too_long = true;
            break;
        }
    }

    if (((flag & BamFlags::BAM_FLAG_UNMAPPED) == 0) && (deletion_too_long)) {
        ++discard_del;
        flag = BamFlags::BAM_FLAG_UNMAPPED;
    }

    bamRecords[bamIndex].qName = readIDs[bamIndex]; // QNAME (1)
    bamRecords[bamIndex].flag = flag; // FLAG (2)
    bamRecords[bamIndex].beginPos = row(align, 1)._array[0]; // POS (4)

    if (bamRecords[bamIndex].flag & BamFlags::BAM_FLAG_UNMAPPED)
        return;

    bamRecords[bamIndex].rID = refId; // RNAME (3)
    bamRecords[bamIndex].mapQ = 60; // MAPQ (5)
    bamRecords[bamIndex].tLen = length(row(align, 1)); // TLEN (9)
    bamRecords[bamIndex].seq = source(row(align, 1)); // SEQ (10)
    if (phred) {
        bamRecords[bamIndex].qual = readQUALs[bamIndex]; // QUAL (11)
    }
    bamRecords[bamIndex].tags = "ASi"; // Tags (12+)
    appendRawPod(bamRecords[bamIndex].tags, score);
}

bool compRecords(BamAlignmentRecord const& left, BamAlignmentRecord const& right)
{
    if (left.qName < right.qName)
        return true;

    if (left.qName > right.qName)
        return false;

    if (left.beginPos < right.beginPos)
        return true;

    return false;
}

void sortRecords(String<BamAlignmentRecord>& bamRecords)
{
    std::sort(begin(bamRecords), end(bamRecords), compRecords);
}

void trimName(String<char>& readName)
{
    for (uint32_t i = 0; i < length(readName); ++i) {
        if (readName[i] == ' ') {
            resize(readName, i);
            return;
        }
    }
}

void adjustFlags(String<BamAlignmentRecord>& bamRecords, const AppOptions& options)
{
    int32_t TLEN;

    for (uint32_t i = 0; i < length(bamRecords) - 1; ++i) {
        if (bamRecords[i].qName == bamRecords[i + 1].qName) {
            bool left_read_bad = (bamRecords[i].flag & BamFlags::BAM_FLAG_UNMAPPED);
            bool right_read_bad = (bamRecords[i + 1].flag & BamFlags::BAM_FLAG_UNMAPPED);

            if (left_read_bad || right_read_bad) {
                bamRecords[i].flag = BamFlags::BAM_FLAG_UNMAPPED;
                bamRecords[i + 1].flag = BamFlags::BAM_FLAG_UNMAPPED;

                discard_bad_mate += (left_read_bad ^ right_read_bad);
            }
            else {
                // left read
                bamRecords[i].flag |= BamFlags::BAM_FLAG_MULTIPLE | BamFlags::BAM_FLAG_ALL_PROPER | BAM_FLAG_FIRST;
                if (bamRecords[i + 1].flag & BamFlags::BAM_FLAG_RC)
                    bamRecords[i].flag |= BamFlags::BAM_FLAG_NEXT_RC;
                bamRecords[i].pNext = bamRecords[i + 1].beginPos;
                bamRecords[i].rNextId = bamRecords[i + 1].rID;

                // right read
                bamRecords[i + 1].flag |= BamFlags::BAM_FLAG_MULTIPLE | BamFlags::BAM_FLAG_ALL_PROPER | BAM_FLAG_LAST;
                if (bamRecords[i].flag & BamFlags::BAM_FLAG_RC)
                    bamRecords[i + 1].flag |= BamFlags::BAM_FLAG_NEXT_RC;
                bamRecords[i + 1].pNext = bamRecords[i].beginPos;
                bamRecords[i + 1].rNextId = bamRecords[i].rID;

                // determine length of template
                TLEN = bamRecords[i + 1].beginPos + bamRecords[i + 1].tLen - bamRecords[i].beginPos;
                bamRecords[i].tLen = TLEN;
                bamRecords[i + 1].tLen = -TLEN;

                if (TLEN < options.minimum_tlen) {
                    bamRecords[i].flag = BamFlags::BAM_FLAG_UNMAPPED;
                    bamRecords[i + 1].flag = BamFlags::BAM_FLAG_UNMAPPED;

                    discard_tlen += 2;
                }
                else {
                    accepted += 2;
                }
            }
            ++i;
        }
        else {
            if ((bamRecords[i].flag & BamFlags::BAM_FLAG_UNMAPPED) == 0) {
                ++discard_singles;
            }
            bamRecords[i].flag = BamFlags::BAM_FLAG_UNMAPPED;
        }
    }
}

void show_progress()
{
    boost::progress_display show_progress(total);
    uint32_t now = 0, previous = 0;

    while (progress != total) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        now = progress.load();
        show_progress += (now - previous);
        previous = now;
    }
}

void align(String<BamAlignmentRecord>& bamRecords,
    StringSet<RefSeq_t> const& refSeqs,
    StringSet<CharString> const& readIDs,
    StringSet<String<Dna5>> const& readSeqs,
    StringSet<CharString> const& readQUALs,
    Ednafull const& scoringScheme,
    AlignConfig<true, false, false, true> const& alignConfig,
    AppOptions const& options,
    uint32_t start,
    uint32_t end)
{
    uint16_t bamFlag = 0;
    String<Dna5> revComPatternSeq;

    // Preparing the align Objects
    String<Align<RefSeq_t>> alignObjs, alignObjsRevComp;
    resize(alignObjs, length(refSeqs));
    resize(alignObjsRevComp, length(refSeqs));

    int maxScore, maxScoreRevComp, maxId, maxIdRevComp;

    for (uint32_t i = 0; i < length(refSeqs); ++i) {
        resize(rows(alignObjs[i]), 2);
        resize(rows(alignObjsRevComp[i]), 2);
        assignSource(row(alignObjs[i], 0), refSeqs[i]);
        assignSource(row(alignObjsRevComp[i], 0), refSeqs[i]);
    }

    if (length(readSeqs) < end)
        end = length(readSeqs);

    progress += (end - start);

    for (uint32_t i = start; i < end; ++i) {
        // align forward
        maxScore = minValue<int>();
        maxId = -1;
        for (uint32_t j = 0; j < length(alignObjs); ++j) {
            assignSource(row(alignObjs[j], 1), readSeqs[i]);
            int result = globalAlignment(alignObjs[j], scoringScheme, alignConfig, maxScore);
            if (maxScore < result) {
                maxScore = result;
                maxId = j;
            }
        }

        // align backward
        revComPatternSeq = readSeqs[i];
        reverseComplement(revComPatternSeq);
        maxScoreRevComp = minValue<int>();
        maxIdRevComp = -1;
        for (uint32_t j = 0; j < length(alignObjsRevComp); ++j) {
            assignSource(row(alignObjsRevComp[j], 1), revComPatternSeq);
            int result = globalAlignment(alignObjsRevComp[j], scoringScheme, alignConfig, maxScoreRevComp);
            if (maxScoreRevComp < result) {
                maxScoreRevComp = result;
                maxIdRevComp = j;
            }
        }

        // check which score is higher and choose the associated alignment
        int bestScore = std::max(maxScore, maxScoreRevComp);
        int normScore = bestScore / (length(revComPatternSeq) * 5.0) * 255.0;
        int bestId;
        Align<RefSeq_t>& bestAlign = (maxScore > maxScoreRevComp) ? (bamFlag = 0, bestId = maxId, alignObjs[maxId]) : (bamFlag = BAM_FLAG_RC, bestId = maxIdRevComp, alignObjsRevComp[maxIdRevComp]);

        if (normScore <= options.minScore) {
            ++discard_score;
            bamFlag = BamFlags::BAM_FLAG_UNMAPPED;
        }
        storeRecord(bamRecords, bamFlag, bestAlign, readIDs, readQUALs, bestId, normScore, i, options.numDeletion, options.Phred);
    }
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const** argv)
{
    // Parse the command line.
    ArgumentParser parser;
    AppOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // Abort if there was an error parsing
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    if (options.quiet)
        options.verbosity_level = 0;

    if (options.verbosity_level) {
        std::cout << "pidalign 0.1\n";
        std::cout << "Parameters:  -t: " << options.numThreads << '\n';
        std::cout << "             -b: " << options.blockSize << '\n';
        std::cout << "             -s: " << options.minScore << '\n';
        std::cout << "             -d: " << options.numDeletion << '\n';
        std::cout << "             -T: " << options.minimum_tlen << '\n';
        std::cout << "             -Q: " << (options.Phred ? "Phred Scores written" : "Phred Scores NOT written") << '\n';
    }

    // Reference
    SeqFileIn seqInRef(toCString(options.refFileName));
    StringSet<String<char>> refIds;
    StringSet<RefSeq_t> refSeqs;
    readRecords(refIds, refSeqs, seqInRef);

    // Preparation of header of resulting sam file
    BamHeader bamHeader;
    BamHeaderRecord bamHeaderRecord;
    setTagValue("VN", "1.5", bamHeaderRecord);
    setTagValue("SO", "queryname", bamHeaderRecord);
    appendValue(bamHeader, bamHeaderRecord);
    StringSet<String<char>> nameStore;
    NameStoreCache<StringSet<String<char>>> nameStoreCache(nameStore);
    BamIOContext<> bamContext(nameStore, nameStoreCache);
    for (unsigned i = 0; i < length(refSeqs); ++i) {
        clear(bamHeaderRecord);
        bamHeaderRecord.type = BAM_HEADER_REFERENCE;
        setTagValue("SN", refIds[i], bamHeaderRecord);
        setTagValue("LN", std::to_string(length(refSeqs[i])), bamHeaderRecord);
        appendValue(bamHeader, bamHeaderRecord);
        nameToId(contigNamesCache(bamContext), refIds[i]);
    }

    // Alignment helpers
    Ednafull scoringScheme(-1, -20);
    AlignConfig<true, false, false, true> alignConfig;

    // Pattern
    if (options.verbosity_level)
        std::cout << "(1)\tReading input\n";

    StringSet<CharString> readIDs;
    StringSet<String<Dna5>> readSeqs;
    StringSet<CharString> readQUALs;

    for (const auto& i : options.patternFileNames) {
        // reading the fastq files
        SeqFileIn seqInPattern(toCString(i));
        while (!atEnd(seqInPattern)) {
            uint32_t newLength = length(readIDs) + 1;
            resize(readIDs, newLength);
            resize(readSeqs, newLength);

            if (options.Phred) {
                // with Phred scores
                resize(readQUALs, newLength);
                readRecord(back(readIDs), back(readSeqs), back(readQUALs), seqInPattern);
            }
            else {
                readRecord(back(readIDs), back(readSeqs), seqInPattern);
            }

            trimName(back(readIDs));
        }
    }

    // prepare the resulting bam record
    total = length(readIDs);
    String<BamAlignmentRecord> bamRecords;
    resize(bamRecords, total);

    // parallelized alignment
    if (options.verbosity_level)
        std::cout << "(2)\tAligning reads\n";
    std::thread progress_bar_thread((options.verbosity_level >= 2 ? show_progress : [] {}));
    boost::threadpool::pool threadPool(options.numThreads);
    uint32_t blockSize = options.blockSize;

    auto start = std::chrono::high_resolution_clock::now();
    for (uint32_t currentReadId = 0; currentReadId < total; currentReadId += blockSize) {
        threadPool.schedule([&, currentReadId, blockSize] {
				   align(bamRecords, refSeqs, readIDs, readSeqs, readQUALs, scoringScheme, alignConfig, options, currentReadId, currentReadId + blockSize);
        });
    }
    progress_bar_thread.join();
    threadPool.wait();

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    if (options.verbosity_level)
        std::cout << "\tElapsed time: " << duration.count() << " seconds\n";

    // post-alignment processing
    if (options.verbosity_level)
        std::cout << "(3)\tPost-processing\n";
    sortRecords(bamRecords);
    adjustFlags(bamRecords, options);

    // output
    if (options.verbosity_level)
        std::cout << "(4)\tWriting output\n\n";

    String<char, MMap<>> outFile;
    open(outFile, toCString(options.outputFileName), OPEN_WRONLY | OPEN_CREATE);
    open(outFile, toCString(options.outputFileName));
    write(outFile, bamHeader, bamContext, Sam());
    for (const auto& i : bamRecords) {
        if (!(i.flag & BamFlags::BAM_FLAG_UNMAPPED)) {
            write(outFile, i, bamContext, Sam());
        }
    }

    // statistics
    if (options.verbosity_level) {
        std::cout << std::fixed << std::setprecision(1);
        std::cout << "Input Reads           : " << total << "\n";
        std::cout << "Discarded (bad score) : " << discard_score << " (" << static_cast<double>(discard_score) / total * 100 << "%)\n";
        std::cout << "Discarded (deletions) : " << discard_del << " (" << static_cast<double>(discard_del) / total * 100 << "%)\n";
        std::cout << "Discarded (bad mate)  : " << discard_bad_mate << " (" << static_cast<double>(discard_bad_mate) / total * 100 << "%)\n";
        std::cout << "Discarded (singles)   : " << discard_singles << " (" << static_cast<double>(discard_singles) / total * 100 << "%)\n";
        std::cout << "Discarded (TLEN>)     : " << discard_tlen << " (" << static_cast<double>(discard_tlen) / total * 100 << "%)\n";
        std::cout << "Accepted Reads        : " << accepted << " (" << static_cast<double>(accepted) / total * 100 << "%)\n";
    }

    return 0;
}
