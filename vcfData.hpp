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
 *  Struct for counting mution types within a higher structure (Gene or Genome)
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
 * Representation of a gene and all information about the mutations on it
 */
class GeneData
{
public:
    GeneData() {}

    /**
     *  Adds position record.  Evaluates a mutation type(s) and constructs respective variant records.
     *  Updates mutation counters
     */
    MutationTypeCounter addPositionRecord(std::shared_ptr<PositionRecord> record);
    
    /**
     *  Returns broken down count of mutations for this gene
     */
    MutationTypeCounter geneMutationTypeCounter() const { return m_geneMutationCounter; }
    
    /**
     * From two strings (reference and alternate base(es)) determines the mutation type
     */
    static MutationType evaluateMutationType(std::string ref, std::string alt);
    /**
     * Is substring included in string?
     */
    static bool isIncluded( std::string str, std::string substr);
    
private:
    /**
     * Representation of a variant record.  References position record which stores additional details
     */
    struct VariantRecord
    {
        std::string ref;
        std::string alt;
        std::shared_ptr<PositionRecord> positionRecordPtr;
    };
    
    MutationTypeCounter                                     m_geneMutationCounter;
    std::map< MutationType, std::vector<VariantRecord> >    m_variantRecords;

};

/**
 *  Class that stores all the data for a the vcf file
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
     *        found in that gene, and a breakdown by variant type as in 2).  Variants with no gene name are printed too with
     *        an empty string instead of a gene name
     */
    void outputResults();

private:
    friend struct ::testGenomeData;
    
    std::string                                     m_vcfFilePath;
    std::vector<std::shared_ptr<PositionRecord>>    m_positionRecords;
    MutationTypeCounter                             m_mutationCounter;
    std::map<std::string, GeneData>                 m_geneData;

};

} // namespace vcf
