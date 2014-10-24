#include <iostream>
#include <fstream>

#include "vcfData.hpp"

namespace vcf
{

  GenomeData::GenomeData( const char * vcfFilePath )
  {
    std::ifstream vcfFile(vcfFilePath);
    if(vcfFile.good())
    {
      std::string line;
      std::getline( vcfFile, line, '\t');
      while( std::getline( vcfFile, line, '\t') )
        {
          std::cout << line.length() << "~~~~~~~" << line[0] << std::endl;
        }
    }


    vcfFile.close();
  }

} // namespace vcf
