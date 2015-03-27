/**
 * Copyright (c) 2014-2015 Jochen Singer, David Seifert
 *
 * This file is part of PIDalign
 *
 * PIDalign is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or any later version.
 *
 * PIDalign is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * PIDalign. If not, see <http://www.gnu.org/licenses/>.
 */

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <seqan/modifier.h>

#include <src/ednafull.hpp>

#include <iostream>
#include <fstream>
#include <algorithm>

using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class SamStruct
// --------------------------------------------------------------------------

//struct SamStruct
//{
//    String<char> readName;
//    String<char> contigName;
//    String<CigarElement<> > cigar;
//    String<char> seq;
//    unsigned pos;
//    int score;
//};

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

struct AppOptions
{
    // The first (and only) argument of the program is stored here.
    StringSet<CharString> patternFileNames;
    CharString refFileName;
    CharString outputFileName;
    int numThreads;
    unsigned bufferSize;
    int minScore;

    AppOptions() :
        numThreads(1),
        bufferSize(100000),
        minScore(minValue<int>())
    {}
};


// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("pidalign");
    // Set short description, version, and date.
    setShortDescription(parser, "Needleman-Wunsch (Gotoh variant) aligner that produces a SAM output");
    setVersion(parser, "0.1");
    setDate(parser, "Februrary 2015");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB<R1.fastq>\\fP [\\fB<R2.fastq>\\fP]");
    addDescription(parser, "pidalign performs a full-exhaustive Needleman-Wunsch alignment with affine gap costs. For performance reasons, pidalign is multithreaded, and it is generally advised to use the -t option to max out the parallelization.");

    addOption(parser, ArgParseOption("r1", "reads1FileName", "Name of the file containg the first mate pair reads.", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "r1");
    setValidValues(parser, "r1", "fasta fa fna fastq fq fnq");

    addOption(parser, ArgParseOption("r2", "reads2FileName", "Name of the file containg the second mate pair reads.", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "r2");
    setValidValues(parser, "r2", "fasta fa fna fastq fq fnq");

    addOption(parser, ArgParseOption("ref", "referenceFileName", "Name of the reference file.", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "ref");
    setValidValues(parser, "ref", "fasta fa fna");

    addOption(parser, ArgParseOption("o", "outputFileName", "Name of the output file.", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "o");
    setValidValues(parser, "o", "sam");

    addOption(parser, ArgParseOption("t", "threads", "The number of threads to be used.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("b", "bufferSize", "The number of hits stored in a buffer before writing them to disk.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("s", "minScore", "MinScore of alignment in order to be considered.", ArgParseArgument::INTEGER));

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBpidalign\\fP \\fB-r\\fP \\fIreference.fasta\\fP\\fP \\fB-o\\fP \\fIoutput.sam\\fP \\fP\\fIR*.fastq\\fP",
                "Use \\fIreference.fasta\\fP to align \\fP\\fIR*.fastq\\fP with the SAM alignment in \\fIoutput.sam\\fP.");

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
    getOptionValue(options.bufferSize, parser, "b");
    getOptionValue(options.minScore, parser, "s");

    return ArgumentParser::PARSE_OK;
}

template <typename TAlign>
void storeRecord(String<BamAlignmentRecord> & bamRecords,
                unsigned short flag,
                TAlign & align,
                StringSet<String<char> > const & patternIds,
                unsigned refId,
                int score,
                unsigned bamIndex)
{
    int clipBegin = row(align, 1)._array[0];
    int clipEnd = length(row(align, 1)) - row(align, 1)._array[length(row(align, 1)._array) - 1];
    setClippedBeginPosition(row(align, 0), clipBegin);
    setClippedBeginPosition(row(align, 1), clipBegin);
    setClippedEndPosition(row(align, 0), clipEnd);
    setClippedEndPosition(row(align, 1), clipEnd);
    getCigarString(bamRecords[bamIndex].cigar,row(align, 0), row(align, 1), 1000);
    String<CigarElement<> > const & cigarRef = bamRecords[bamIndex].cigar;
    for (unsigned i = 0; i < length(cigarRef); ++i)
        if (cigarRef[i].operation == 'D' && cigarRef[i].count >= 4)
            bamRecords[bamIndex].flag = BamFlags::BAM_FLAG_UNMAPPED;

    if (bamRecords[bamIndex].flag & BamFlags::BAM_FLAG_UNMAPPED)
        return;

    bamRecords[bamIndex].rID = refId;                            // position of ref in header
    bamRecords[bamIndex].qName = patternIds[bamIndex];              //qName
    bamRecords[bamIndex].flag = flag;                               //FLAG
    bamRecords[bamIndex].beginPos = row(align, 1)._array[0] + 1;    //POS
    bamRecords[bamIndex].mapQ = 60;                              //MapQual
    //stream << 0 << "\t";                                            //PNEXT
    //stream << 0 << "\t";                                            //TLen
    bamRecords[bamIndex].seq = source(row(align, 1));               //SEQ
    //stream << "*" << "\n";                                          //QUAL
    bamRecords[bamIndex].tags = "ASi";
    appendRawPod(bamRecords[bamIndex].tags, score);
}

bool compRecords(BamAlignmentRecord const & left, BamAlignmentRecord const & right)
{
    if (left.qName < right.qName)
        return true;

    if (left.qName > right.qName)
        return false;

    if (left.beginPos < right.beginPos)
        return true;

    return false;
}

void
sortRecords(String<BamAlignmentRecord> & bamRecords)
{
    std::sort(begin(bamRecords), end(bamRecords), compRecords);
}

void
trimName(String<char> & readName)
{
    for (unsigned i = 0; i < length(readName); ++i)
        if (readName[i] == ' ')
        {
            resize(readName, i + 1);
            return;
        }
}

void adjustFlags(String<BamAlignmentRecord> & bamRecords)
{
    for (unsigned i = 0; i < length(bamRecords) - 1; ++i)
    {
        if (bamRecords[i].qName == bamRecords[i+1].qName)
        {
            if ( (bamRecords[i].flag & BamFlags::BAM_FLAG_UNMAPPED) || 
                    (bamRecords[i+1].flag & BamFlags::BAM_FLAG_UNMAPPED) )
            {
                bamRecords[i].flag = BamFlags::BAM_FLAG_UNMAPPED;
                bamRecords[i+1].flag = BamFlags::BAM_FLAG_UNMAPPED;
            }
            else
            {
                bamRecords[i].flag |= BamFlags::BAM_FLAG_MULTIPLE | BamFlags::BAM_FLAG_ALL_PROPER | BAM_FLAG_FIRST;
                bamRecords[i].pNext = bamRecords[i+1].rID;
                bamRecords[i+1].flag |= BamFlags::BAM_FLAG_MULTIPLE | BamFlags::BAM_FLAG_ALL_PROPER | BAM_FLAG_LAST;
                bamRecords[i+1].pNext = bamRecords[i].rID;
            }
            ++i;
        }
        else
        {
            bamRecords[i].flag = BamFlags::BAM_FLAG_UNMAPPED;
        }
    }
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    typedef String<Iupac> TRefSeq;

    // Parse the command line.
    ArgumentParser parser;
    AppOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    std::cout << "pidalign\n"
              << "========\n\n";

    // Reference
    SeqFileIn seqInRef(toCString(options.refFileName));
    StringSet<String<char> > refIds;
    StringSet<TRefSeq> refSeqs;
    readRecords(refIds, refSeqs, seqInRef);

    // Preparation of header of resulting sam file
    //std::ofstream stream(toCString(options.outputFileName), std::ios::out | std::ios::app);
    //stream << "@HD\tVN:1.5\tSO:coordinate\n";
    //for (unsigned i = 0; i < length(refSeqs); ++i)
    //  stream << "@SQ\tSN:" << refIds[i] << "\tLN:" << length(refSeqs[i]) << "\n";
    //stream.close();
    BamHeader bamHeader;
    BamHeaderRecord bamHeaderRecord;
    setTagValue("VN", "1.5", bamHeaderRecord);
    setTagValue("SO", "queryname", bamHeaderRecord);
    appendValue(bamHeader, bamHeaderRecord);
    StringSet<String<char> > nameStore;
    NameStoreCache<StringSet<String<char> > > nameStoreCache(nameStore);
    BamIOContext<> bamContext(nameStore, nameStoreCache);
    for (unsigned i = 0; i < length(refSeqs); ++i)
    {
        clear(bamHeaderRecord);
        bamHeaderRecord.type = BAM_HEADER_REFERENCE;
        setTagValue("SN", refIds[i], bamHeaderRecord);
        setTagValue("LN", std::to_string(length(refSeqs[i])), bamHeaderRecord);
        appendValue(bamHeader, bamHeaderRecord);
        //unsigned globalRefId = nameToId(contigNamesCache(bamContext), refIds[i]);
        nameToId(contigNamesCache(bamContext), refIds[i]);
        //contigLengths(bamContext)[globalRefId] = length(refSeqs[i]);
    }

    // Alignment helpers
    Ednafull scoringScheme(-1, -20);
    AlignConfig<true, false, false, true> alignConfig;

    // Pattern
    StringSet<CharString> patternIds;
    StringSet<String<Dna5> > patternSeqs;
    String<Dna5> revComPatternSeq;
    for (unsigned fileCounter = 0; fileCounter < length(options.patternFileNames); ++fileCounter)
    {
        // reading the fastq files
        SeqFileIn seqInPattern(toCString(options.patternFileNames[fileCounter]));
        while (!atEnd(seqInPattern))
        {
            unsigned newLength = length(patternIds) + 1;
            resize(patternIds, newLength);
            resize(patternSeqs, newLength);
            readRecord(back(patternIds), back(patternSeqs), seqInPattern);
            trimName(back(patternIds));
        }
        //readRecords(patternIds, patternSeqs, seqInPattern);
    }

    // prepare the resulting bam record
    String<BamAlignmentRecord> bamRecords;
    resize(bamRecords, length(patternIds));

    // This is the parallelization starting point
    unsigned readIndex= 0;
    unsigned short bamFlag = 0;
    SEQAN_OMP_PRAGMA(parallel num_threads(options.numThreads))
    for (; ;)
    {
        unsigned threadReadIndex;
        SEQAN_OMP_PRAGMA(critical (adjust readIndex))
        {
            if (readIndex >= length(patternSeqs))
                break;

            threadReadIndex = readIndex;
            readIndex += options.bufferSize;
        }

        // Preparing the align Objects
        String<Align<TRefSeq> > alignObjs, alignObjsRevComp;
        resize(alignObjs, length(refSeqs));
        resize(alignObjsRevComp, length(refSeqs));

        int maxScore, maxScoreRevComp, maxId, maxIdRevComp;

        for (unsigned i = 0; i < length(refSeqs); ++i)
        {
            resize(rows(alignObjs[i]), 2);
            resize(rows(alignObjsRevComp[i]), 2);
            assignSource(row(alignObjs[i], 0), refSeqs[i]);
            assignSource(row(alignObjsRevComp[i], 0), refSeqs[i]);
        }

        unsigned threadIndexEnd = (length(patternSeqs) < readIndex + options.bufferSize - 1) ? length(patternSeqs) : readIndex + options.bufferSize - 1;
        for (unsigned i = threadReadIndex; i < threadIndexEnd; ++i)
        {
            // align forward
            maxScore = minValue<int>();
            maxId = -1;
            for (unsigned j = 0; j < length(alignObjs); ++j)
            {
                assignSource(row(alignObjs[j], 1), patternSeqs[i]);
                int result = globalAlignment(alignObjs[j], scoringScheme, alignConfig, maxScore);
                if (maxScore < result)
                {
                    maxScore = result;
                    maxId = j;
                }
            }

            // align backward
            revComPatternSeq = patternSeqs[i];
            reverseComplement(revComPatternSeq);
            maxScoreRevComp = minValue<int>();
            maxIdRevComp = -1;
            for (unsigned j = 0; j < length(alignObjsRevComp); ++j)
            {
                assignSource(row(alignObjsRevComp[j], 1), revComPatternSeq);
                int result = globalAlignment(alignObjsRevComp[j], scoringScheme, alignConfig, maxScoreRevComp);
                if (maxScoreRevComp < result)
                {
                    maxScoreRevComp = result;
                    maxIdRevComp = j;
                }
            }

            // check whehter the score for the forward alignment is
            // better than the one of the backward
            if (maxScore > maxScoreRevComp)
            {
                int normScore = maxScore / (length(revComPatternSeq) * 5.0) * 255.0;
                if (normScore > options.minScore)
                    bamFlag = 0;
                else 
                    bamFlag = BamFlags::BAM_FLAG_UNMAPPED;

                storeRecord(bamRecords, bamFlag, alignObjs[maxId], patternIds, maxId, normScore, i);
            }
            else
            {
                int normScore = maxScoreRevComp / (length(revComPatternSeq) * 5.0) * 255.0;
                if (normScore > options.minScore)
                    bamFlag = BAM_FLAG_RC;
                else 
                    bamFlag = BamFlags::BAM_FLAG_UNMAPPED;

                storeRecord(bamRecords, bamFlag, alignObjsRevComp[maxIdRevComp], patternIds, maxIdRevComp, normScore, i);
            }
        }
    }
    sortRecords(bamRecords);
    adjustFlags(bamRecords);

    String<char, MMap<> > outFile;
    open(outFile, toCString(options.outputFileName));
    write(outFile, bamHeader, bamContext, Sam());
    for (unsigned i = 0; i < length(bamRecords); ++i)
    {
        if (!(bamRecords[i].flag & BamFlags::BAM_FLAG_UNMAPPED))
            write(outFile, bamRecords[i], bamContext, Sam());
    }

    return 0;
}


