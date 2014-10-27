#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>

struct testSplitWithDelimiter;
struct testSplitAndCapitalise;
struct testExtractGeneNames;

namespace vcf
{
    
namespace{
    const size_t numberOfVcfFields = 10;
}
    
class VcfParserError : public std::runtime_error {
public:
    explicit VcfParserError(const std::string &msg) : std::runtime_error(msg) {}
};

/**
 *  Data for a single entry in the vcf file.
 */
struct PositionRecord
{
    std::string chrom;
    size_t pos;
    std::string id;
    std::string ref;
    std::vector<std::string> alt;
    bool pass;
    std::vector<std::string> geneNames;
};

/**
 * Simple parser of the vcf file format.  Parses only data records, no metadata or headers
 */
class VcfParser
{
public:
    explicit VcfParser( std::string vcfFilePath );
    
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
    std::string capitalise(std::string str) const;
    /**
     * Extracts gene names from the INFO field (comma delimited list of gene names).
     * If gene name(s) not specified returns an empty list
     */
    std::vector<std::string> extractGeneNames(std::string str) const;
    
    friend struct ::testSplitWithDelimiter;
    friend struct ::testSplitAndCapitalise;
    friend struct ::testExtractGeneNames;
    
    std::ifstream   m_vcfFile;
    bool            m_hasNextRow;
    PositionRecord  m_nextValidRecord;
    
};

}
