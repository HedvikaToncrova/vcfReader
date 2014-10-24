#include <fstream>
#include <sstream>
#include <algorithm>
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
            auto parsedLine = splitWithDelimiter(line, '\t');
            m_hasNextRow = constructNextValidRecord(parsedLine);
        }
    }
}
    
bool VcfParser::constructNextValidRecord( const std::vector<std::string>& parsedLine)
{
    if (parsedLine.size() != numberOfVcfFields)
    {
        throw VcfParserError("Incorrect number of fields in this row");
    }

    if( parsedLine[6].compare("PASS") == 0)
    {
        try
        {
            m_nextValidRecord.chrom = parsedLine[0];
            m_nextValidRecord.pos = boost::lexical_cast<size_t>(parsedLine[1]);
            m_nextValidRecord.id = parsedLine[2];
            m_nextValidRecord.ref = splitAndCapitalise(parsedLine[3]);
            m_nextValidRecord.alt = splitAndCapitalise(parsedLine[4]);
            m_nextValidRecord.pass = true;
            m_nextValidRecord.geneName = extractGeneName(parsedLine[7]);
            
        }
        catch( boost::bad_lexical_cast const& )
        {
            throw VcfParserError("Cannot parse " + parsedLine[1] + " as integer");
        }
        return true;
    }
    

    return false;
}
    
std::vector<std::string> VcfParser::splitWithDelimiter(std::string str ,char delim) const
{
    std::stringstream lineStream(str);
    std::vector<std::string> parsedLine;
    std::string cell;
    
    while(std::getline(lineStream,cell,delim))
    {
        boost::algorithm::trim(cell);
        parsedLine.push_back(cell);
    }
    return parsedLine;
}
    
std::vector<std::string> VcfParser::splitAndCapitalise(std::string str) const
{
    std::transform(str.begin(), str.end(),str.begin(), ::toupper);
    return splitWithDelimiter(str, ',');
}
    
std::string VcfParser::extractGeneName(std::string str) const
{
    auto splitString = splitWithDelimiter(str, ';');
    std::string geneName;
    for( auto s : splitString )
    {
        if(boost::starts_with(s, "Gene="))
        {
            geneName = s.substr(5);
        }
    }
    return geneName;
}
    
} // namespace vcf