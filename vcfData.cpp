#include <iostream>

#include "vcfData.hpp"
#include "vcfParser.hpp"

namespace vcf
{

GenomeData::GenomeData( std::string vcfFilePath )
{
    auto parser = VcfParser(vcfFilePath);
    while (parser.hasNextValidRecord() )
    {
        m_positionRecords.push_back(std::make_shared<PositionRecord>(parser.getNextValidRecord()));
        std::cout << "next valid record with position " << m_positionRecords.back()->pos << std::endl;
    }
}
    
mutationTypeCounter_t GeneData::addPositionRecord(std::shared_ptr<PositionRecord> record)
{
    mutationTypeCounter_t positioMutations;
    for( auto ref : record->ref )
    {
        for(auto alt : record->alt )
        {
            //auto mutationType = evaluateMutationType(ref, alt);
            VariantRecord varRec = {ref, alt, record};
            //m_variantRecords.insert( {mutationType, varRec} );
        }
    }
    return positioMutations;
}
    
MutationType GeneData::evaluateMutationType(std::string ref, std::string alt)
{
    if(ref == alt)
    {
	       return MutationType::IDEN;
    }
    size_t refLength = ref.length();
    size_t altLength = alt.length();
    
    if( refLength == altLength && refLength == 1)
    {
        return MutationType::SVN;
    }
    
    if( refLength < altLength && isIncluded( alt, ref ) )
    {
        return MutationType::INS;
    }
    
    if( refLength > altLength && isIncluded( ref, alt ) )
    {
        return MutationType::DEL;
    }
    
    //else
    return MutationType::MVN;
}
    
bool GeneData::isIncluded( std::string str, std::string substr)
{
    size_t strLength = str.length();
    size_t substrLength = substr.length();
    
    if( strLength < substrLength)
    {
        return false;
    }
    
    for(size_t i = 0; i <= (strLength - substrLength); i++)
    {
        if( str.compare(i, substrLength, substr) == 0 )
        {
            return true;
        }
    }
    return false;
}

} // namespace vcf
