#pragma once

#include <map>
#include <vector>
#include <string>
#include <memory>

#include "vcfParser.hpp"

namespace vcf {

enum class MutationType
{
    IDEN, SVN, INS, DEL, MVN
};

typedef std::map<MutationType, int> mutationTypeCounter_t;
    
class GeneData
{
public:
    GeneData();


    struct VariantRecord
    {
        std::string ref;
        std::string alt;
        std::shared_ptr<PositionRecord> positionRecordPtr;
    };

    mutationTypeCounter_t addPositionRecord(std::shared_ptr<PositionRecord> record);
    std::map<MutationType, size_t> getMutationCounts;
    
    static MutationType evaluateMutationType(std::string ref, std::string alt);
    static bool isIncluded( std::string str, std::string substr);
private:
    std::map< MutationType, std::vector<VariantRecord> > m_variantRecords;

};

class GenomeData
{
public:
    GenomeData( std::string vcfFilePath );

private:
    std::vector<std::shared_ptr<PositionRecord>>    m_positionRecords;
    mutationTypeCounter_t                           m_mutationCounter;
    std::map<std::string, GeneData>                 m_geneData; // maybe would be better as a map from the name

};

} // namespace vcf
