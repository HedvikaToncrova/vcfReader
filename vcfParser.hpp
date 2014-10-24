#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

namespace vcf
{

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
    
    std::ifstream   m_vcfFile;
    bool            m_hasNextRow;
    PositionRecord  m_nextValidRecord;
    
};

}