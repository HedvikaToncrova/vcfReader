#include <iostream>

#include "vcfData.hpp"

int main(int argc, const char * argv[])
{
    if( argc == 2 )
    {
        std::string vcfInputPath( argv[1], strlen(argv[1]));
        vcf::GenomeData genomeData( vcfInputPath );
        genomeData.outputResults();
    }
    else
    {
        vcf::GenomeData genomeData( "testData/testVcfFile.vcf" );
        genomeData.outputResults();
        return 1;
    }
    return 0;
}

