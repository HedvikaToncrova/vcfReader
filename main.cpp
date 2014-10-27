#include <iostream>
#include <string.h>

#include "vcfData.hpp"

int main(int argc, const char * argv[])
{
    if( argc == 2 )
    {
        std::string vcfInputPath( argv[1], strlen(argv[1]) );
        vcf::GenomeData genomeData( vcfInputPath );
        genomeData.outputResults();
    }
    else
    {
        std::cout << "Invalid program options.  Usage: [ input vcf file ]" << std::endl;
        return 1;
    }
    return 0;
}

