#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>

namespace vcf
{
    
namespace{
    const size_t numberOfVcfFields = 10;
}
    
class VcfParserError : public std::runtime_error {
public:
    explicit VcfParserError(const std::string &msg) : std::runtime_error(msg) {}
};

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
    
class VcfParser
{
public:
    explicit VcfParser( const char * vcfFilePath );
    
    bool hasNextValidRecord() const { return m_hasNextRow; }
    PositionRecord getNextValidRecord();
private:
    void assignNextRecord();
    bool constructPositionRecord( const std::vector<std::string>& parsedLine);
    std::vector<std::string> splitAndCapitalise(std::string str);
    std:string extractGeneName(std::string str);
    
    
    std::ifstream   m_vcfFile;
    bool            m_hasNextRow;
    PositionRecord  m_nextValidRecord;
    
};

}