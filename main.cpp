#include <fstream>
#include <iostream>
#include "Phasing.h"
#include "Haplotag.h"
#include "ModCall.h"
#include "SomaticHaplotag.h"
#include "PurityEstimation.h"


#define PROGRAM_BIN "main"
#define VERSION "1.0.0"


static std::string version = VERSION;

static const char *STRIDE_USAGE_MESSAGE =
"Version: " VERSION " \n"
"Usage: " PROGRAM_BIN " <command> [options]\n"  
"[Germline Mode]\n"
"    phase                  run phasing algorithm.\n"
"    haplotag               tag reads by haplotype.\n\n"
//"    modcall                convert bam file to modification vcf file.\n\n"
"[Somatic Mode]\n"
"    somatic_haplotag       tag reads by somatic haplotype algorithm.\n"
"    estimate_purity        estimate tumor purity.\n"
"\n";

int main(int argc, char** argv)
{
    if(argc <= 1)
    {
        std::cout << STRIDE_USAGE_MESSAGE;
        return 0;
    }
    
    std::string command(argv[1]);
    
    if(command=="phase")
    {
        PhasingMain(argc - 1, argv + 1, version);
    }
    else if(command=="haplotag")
    {
        HaplotagMain(argc - 1, argv + 1, version);
    }
    else if(command=="somatic_haplotag")
    {
        SomaticHaplotagMain(argc - 1, argv + 1, version);
    }
    else if(command=="estimate_purity")
    {
        PurityEstimMain(argc - 1, argv + 1, version);
    }
    else if(command=="modcall")
    {
         ModCallMain(argc - 1, argv + 1, version);
    }
    else{
        std::cout << STRIDE_USAGE_MESSAGE;
        return 0;
    }

    return 0;
}
