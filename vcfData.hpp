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

    MutationTypeCounter addPositionRecord(std::shared_ptr<PositionRecord> record);
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

class GenomeData
{
public:
    GenomeData( std::string vcfFilePath );
    
    void outputResults();

private:
    friend struct ::testGenomeData;
    
    std::string                                     m_vcfFilePath;
    std::vector<std::shared_ptr<PositionRecord>>    m_positionRecords;
    MutationTypeCounter                             m_mutationCounter;
    std::map<std::string, GeneData>                 m_geneData; // maybe would be better as a map from the name

};

} // namespace vcf
