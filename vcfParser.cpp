#include <fstream>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

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
                boost::algorithm::trim(cell);
                parsedLine.push_back(cell);
            }
            m_hasNextRow = constructPositionRecord(parsedLine);
        }
    }
}
    
bool VcfParser::constructPositionRecord( const std::vector<std::string>& parsedLine)
{
    if (parsedLine.size() != numberOfVcfFields)
    {
        throw VcfParserError("Incorrect number of fields in this row");
    }
    
    std::cout << parsedLine[1] << "  " << parsedLine[6] << "    " << parsedLine[6].compare("PASS") << std::endl;
    
    if( parsedLine[6].compare("PASS") == 0)
    {
        PositionRecord result;
        try
        {
            result.chrom = parsedLine[0];
            result.pos = boost::lexical_cast<size_t>(parsedLine[1]);
            result.id = parsedLine[2];
            result.ref = splitAndCapitalise(parsedLine[3]);
            result.alt = splitAndCapitalise(parsedLine[4]);
            result.pass = true;
            result.geneName = extractGeneName(parsedLine[7]);
        }
        catch( boost::bad_lexical_cast const& )
        {
            throw VcfParserError("Cannot parse " + parsedLine[1] + " as integer");
        }
        std::cout << "Returning true for position " << result.pos << std::endl;
        return true;
    }
    

    return false;
}
    
} // namespace vcf