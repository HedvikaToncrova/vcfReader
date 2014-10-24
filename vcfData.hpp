#include <map>
#include <vector>
#include <string>

namespace vcf {

struct PositionRecord
{
    std::string chrom;
    size_t pos;
    std::string id;
    std::vector<std::string> ref;
    std::vector<std::string> alt;
    bool pass;
    std::string geneName;
};

enum class MutationType
{
    SVN, INS, DEL, MVN, IDEN
};

class GeneData
{
public:
    GeneData();


    struct VariantRecord
    {
        std::string chrom;
        size_t pos;
        std::string id;
        std::string ref;
        std::string alt;
    };

    MutationType addVarianRecord(VariantRecord record);
    std::map<MutationType, size_t> getMutationCounts; 
private:
    std::map< MutationType, std::vector<VariantRecord> > m_variantRecords;

};

class GenomeData
{
public:
  //  GenomeData( std::string vcfFilePath );
  GenomeData( const char * vcfFilePath );

private:
    std::map<MutationType, int>   m_mutationCounter;
    std::map<std::string, GeneData>  m_geneData; // maybe would be better as a map from the name

};

} // namespace vcf
