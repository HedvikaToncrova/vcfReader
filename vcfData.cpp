#include <iostream>

#include "vcfData.hpp"
#include "vcfParser.hpp"

namespace vcf
{
    
MutationTypeCounter::MutationTypeCounter()
{
    counter[MutationType::IDEN] = 0;
    counter[MutationType::SVN] = 0;
    counter[MutationType::INS] = 0;
    counter[MutationType::DEL] = 0;
    counter[MutationType::MVN] = 0;
}

size_t MutationTypeCounter::totalNumberOfMutations()
{
    size_t total = 0;
    for( auto mutation : counter )
    {
        total += mutation.second;
    }
    return total - counter[MutationType::IDEN];
}
    
std::ostream& operator<< ( std::ostream& os, const MutationTypeCounter& mutationTypeCounter)
{
    os << "\t SVN: " <<  mutationTypeCounter.counter.at(MutationType::SVN) << std::endl
       << "\t INS: " <<  mutationTypeCounter.counter.at(MutationType::INS) << std::endl
       << "\t DEL: " <<  mutationTypeCounter.counter.at(MutationType::DEL) << std::endl
       << "\t MVN: " <<  mutationTypeCounter.counter.at(MutationType::MVN) << std::endl
       << "\t IDENTITY: " <<  mutationTypeCounter.counter.at(MutationType::IDEN) << std::endl;
    
    return os;
}
    
GenomeData::GenomeData( std::string vcfFilePath ) :
    m_vcfFilePath(vcfFilePath)
{
    VcfParser parser(m_vcfFilePath);
    while (parser.hasNextValidRecord() )
    {
        m_positionRecords.push_back(std::make_shared<PositionRecord>(parser.getNextValidRecord()));
        auto geneName = m_positionRecords.back()->geneName;
        auto geneIt = m_geneData.find(geneName);
        if( geneIt == m_geneData.end() )
        {
            GeneData gd;
            geneIt = m_geneData.insert({geneName, gd}).first;
        }
        auto positionMutations = geneIt->second.addPositionRecord(m_positionRecords.back());

        for(auto mutationType : positionMutations.counter)
        {
            m_mutationCounter.counter[mutationType.first] += mutationType.second;
        }
    }
}
    
void GenomeData::outputResults()
{
    std::cout << "**********  RESULTS IN TOTAL  ************" << std::endl;
    std::cout << "Total number of mutations in " << m_vcfFilePath << ": "
              << m_mutationCounter.totalNumberOfMutations() << std::endl
              << m_mutationCounter << std::endl;

    std::cout << "\n \n" << "*********  RESULTS PER GENE ************* " << std::endl;
    for( auto gene : m_geneData )
    {
        auto geneMutationCounter = gene.second.geneMutationTypeCounter();
        std::cout << gene.first << ": " << geneMutationCounter.totalNumberOfMutations() << std::endl
            << geneMutationCounter << std::endl;
        
    }
}
    
MutationTypeCounter GeneData::addPositionRecord(std::shared_ptr<PositionRecord> record)
{
    MutationTypeCounter positionMutations;
    for(auto alt : record->alt )
    {
        auto mutationType = evaluateMutationType(record->ref, alt);
        positionMutations.counter[mutationType]++;
        m_geneMutationCounter.counter[mutationType]++;
        VariantRecord varRec = {record->ref, alt, record};
        if( m_variantRecords.count(mutationType) )
        {
            m_variantRecords[mutationType].push_back(varRec);
        }
        else
        {
            m_variantRecords.insert( {mutationType, {varRec}} );
        }
    }
    return positionMutations;
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
