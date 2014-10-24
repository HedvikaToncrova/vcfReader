#include <fstream>
#include <sstream>

#include "vcfParser.hpp"

namespace vcf{

VcfParser::VcfParser( const char * vcfFilePath ) :
    m_vcfFile(vcfFilePath)
{
    assignNextRecord();
    
}
    
    PositionRecord VcfParser::getNextValidRecord()
{
    PositionRecord currentRecord = m_nextValidRecord;
    assignNextRecord();
    return currentRecord;
}

void VcfParser::assignNextRecord()
{
    m_hasNextRow = false;
    while(m_vcfFile.good() && !m_hasNextRow)
    {
        std::string line;
        std::getline(m_vcfFile,line);
        
        //do not process metadata for the moment
        if(line.front() != '#')
        {
            std::stringstream lineStream(line);
            std::vector<std::string> parsedLine;
            std::string cell;
            
            while(std::getline(lineStream,cell,'\t'))
            {
                parsedLine.push_back(cell);
            }
            m_hasNextRow = constructPositionRecord(parsedLine);
        }
    }
}
    
bool VcfParser::constructPositionRecord( const std::vector<std::string>& parsedLine)
{
    return true;
}
    
} // namespace vcf