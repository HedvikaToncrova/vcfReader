#pragma once

#include <map>
#include <vector>
#include <string>
#include <memory>

#include "vcfParser.hpp"

struct testGenomeData;

namespace vcf {

enum class MutationType
{
    IDEN, SVN, INS, DEL, MVN
};

/**
 *  Struct for counting mutation types within.  Has a nice << operator
 */
struct MutationTypeCounter
{
    MutationTypeCounter();
    std::map<MutationType, int> counter;

    // sum of all mutations (excl. identity)
    size_t totalNumberOfMutations();
};

std::ostream& operator<< ( std::ostream& os, const MutationTypeCounter& mutationTypeCounter);

/**
 *  Processes a position record: evaluates the mutation types and updates total and gene mutation counters
 */
class PositionRecordProcessor
{
public:
    PositionRecordProcessor(PositionRecord positionRecord) :
        m_positionRecord(positionRecord)
    {}
    
    /**
     *  Updates total and gene mutation counters for the mutations recorded in the m_positionRecord
     */
    void updateCounters(std::shared_ptr<MutationTypeCounter> totalMutationCounter,
                        std::shared_ptr< std::map<std::string, MutationTypeCounter> > geneMutationCounter);
    
    /**
     * From two strings (reference and alternate base(es)) determines the mutation type
     */
    static MutationType evaluateMutationType(std::string ref, std::string alt);
    /**
     * Is substring included in string?
     */
    static bool isIncluded( std::string str, std::string substr);
    
private:
    PositionRecord m_positionRecord;
};


/**
 *  Class that stores mutationType counts on gene and total levels for all data in the vcf file
 */
class GenomeData
{
public:
    GenomeData( std::string vcfFilePath );
    
    /**
     *  Outputs to the screen:
     *     1) Total number of mutations represented in the file
     *     2) Breakdown of this count into a mutation type (SVN, INS, DEL, MVN and IDEN)
     *        (no identity should be stored in the vcf format.  It's printed for the completenes and error checking)
     *     3) A list of genes within which at least one variant has been found, the total number of variants
     *        found in that gene, and a breakdown by variant type as in 2).  Variants with no gene name are not printed
     */
    void outputResults();

private:
    friend struct ::testGenomeData;
    
    std::string                                                     m_vcfFilePath;
    std::shared_ptr<MutationTypeCounter>                            m_totalMutationCounter;
    std::shared_ptr<std::map<std::string, MutationTypeCounter>>     m_geneMutationCounter;

};

} // namespace vcf
