#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>

struct testSplitWithDelimiter;
struct testSplitAndCapitalise;
struct testExtractGeneName;

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

/**
 * Simple parser of the vcf file format.  Parses only data records, no metadata or headers
 */
class VcfParser
{
public:
    explicit VcfParser( const char * vcfFilePath );
    
    /**
     * Checks for a next valid record in the vcf file stream
     *
     *  Definition of valid: 1) line does not start with #
     *                       2) record has all 10 specified, tab delimited field
     *                       3) record has a PASS quality filtering
     */
    bool hasNextValidRecord() const { return m_hasNextRow; }
    
    /**
     * Returns the next valid record in the vcf file stream
     * For definition of valid see the declaration of hasNextValidRecord
     */
    PositionRecord getNextValidRecord();
private:
    void assignNextRecord();
    bool constructNextValidRecord( const std::vector<std::string>& parsedLine);
    std::vector<std::string> splitWithDelimiter(std::string str, char delim) const;
    std::vector<std::string> splitAndCapitalise(std::string str) const;
    /**
     * Extracts gene name from the INFO field.  If gene name not specified returns a
     * empty string
     */
    std::string extractGeneName(std::string str) const;
    
    friend struct ::testSplitWithDelimiter;
    friend struct ::testSplitAndCapitalise;
    friend struct ::testExtractGeneName;
    
    std::ifstream   m_vcfFile;
    bool            m_hasNextRow;
    PositionRecord  m_nextValidRecord;
    
};

}