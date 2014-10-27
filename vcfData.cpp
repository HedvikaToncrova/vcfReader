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
    m_vcfFilePath(vcfFilePath),
    m_totalMutationCounter(std::make_shared<MutationTypeCounter>()),
    m_geneMutationCounter(std::make_shared<std::map<std::string, MutationTypeCounter>>())

{
    VcfParser parser(m_vcfFilePath);
    while (parser.hasNextValidRecord() )
    {
        PositionRecordProcessor recordProcessor( parser.getNextValidRecord() );
        recordProcessor.updateCounters(m_totalMutationCounter, m_geneMutationCounter);
    }
}
    
void GenomeData::outputResults()
{
    std::cout << "**********  RESULTS IN TOTAL  ************" << std::endl;
    std::cout << "Total number of mutations in " << m_vcfFilePath << ": "
              << m_totalMutationCounter->totalNumberOfMutations() << std::endl
              << *m_totalMutationCounter << std::endl;
    
    std::cout << "\n \n" << "*********  RESULTS PER GENE ************* " << std::endl;
    for( auto gene : *m_geneMutationCounter )
    {
        std::cout << gene.first << ": " << gene.second.totalNumberOfMutations() << std::endl
            << gene.second << std::endl;
        
    }
}
    
void PositionRecordProcessor::updateCounters(std::shared_ptr<MutationTypeCounter> totalMutationCounter,
                                             std::shared_ptr< std::map<std::string, MutationTypeCounter> > geneMutationCounter)
{
    for(auto alt : m_positionRecord.alt )
    {
        // update total counter
        auto mutationType = evaluateMutationType(m_positionRecord.ref, alt);
        totalMutationCounter->counter[mutationType]++;

        // update gene counters (one position record can be associated with multiple genes)
        for( auto geneName : m_positionRecord.geneNames)
        {
            auto geneIt = geneMutationCounter->find(geneName);
            if( geneIt == geneMutationCounter->end() )
            {
                MutationTypeCounter mtc;
                mtc.counter[mutationType]++;
                geneMutationCounter->insert({geneName, mtc});
            }
            else
            {
                geneIt->second.counter[mutationType]++;
            }
        }
    }
}
    
MutationType PositionRecordProcessor::evaluateMutationType(std::string ref, std::string alt)
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
    
bool PositionRecordProcessor::isIncluded( std::string str, std::string substr)
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
