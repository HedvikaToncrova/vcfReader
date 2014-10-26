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

struct MutationTypeCounter
{
    MutationTypeCounter();
    std::map<MutationType, int> counter;

    size_t totalNumberOfMutations();
};

std::ostream& operator<< ( std::ostream& os, const MutationTypeCounter& mutationTypeCounter);
    
class GeneData
{
public:
    GeneData() {}


    struct VariantRecord
    {
        std::string ref;
        std::string alt;
        std::shared_ptr<PositionRecord> positionRecordPtr;
    };

    MutationTypeCounter addPositionRecord(std::shared_ptr<PositionRecord> record);
    MutationTypeCounter geneMutationTypeCounter() const { return m_geneMutationCounter; }
    
    static MutationType evaluateMutationType(std::string ref, std::string alt);
    static bool isIncluded( std::string str, std::string substr);
private:
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
