#include <iostream>

#include "vcfData.hpp"
#include "vcfParser.hpp"

namespace vcf
{

  GenomeData::GenomeData( const char * vcfFilePath )
  {
      auto parser = VcfParser(vcfFilePath);
      while (parser.hasNextValidRecord() )
      {
          PositionRecord record = parser.getNextValidRecord();
      }
  }

} // namespace vcf
