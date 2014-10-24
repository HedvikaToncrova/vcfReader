#include <iostream>

#include "vcfData.hpp"

int main(int argc, const char * argv[])
{
    if( argc == 2 )
    {
      //std::string vcfInputPath( argv[1], strlen(argv[1])); 
      // vcf::GenomeData genomeData( vcfInputPath );
         vcf::GenomeData genomeData( argv[1] );
    }
    else
    {
        vcf::GenomeData genomeData( "data/testData.vcf" );
        return 1;
    }
    return 0;
}

