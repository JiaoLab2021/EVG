// g++ vcf_breakpoint.cpp -o vcf_breakpoint -lz -O2
#include "../include/vcf_breakpoint.hpp"

using namespace std;

int main_breakpoint(int argc, char** argv)
{
    // Input file
    string vcfFileName;
    
    // Output file prefix
    string prefix = "breakpoint";

    // Breakpoint error magnitude
    int32_t breakpointErrorSize_ = 1;

    // Input parameter
    int c;
    while (true)
    {
        static const struct option long_options[] = {
            {"vcf", required_argument, 0, 'v'},

            {"prefix", required_argument, 0, 'p'},
            {"size", required_argument, 0, 's'},

            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "v:p:s:h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
        case 'v':
            vcfFileName = optarg;
            break;
        case 'p':
            prefix = optarg;
            break;;
        case 's':
            breakpointErrorSize_ = stol(optarg);
            break;
        case 'h':
        case '?':
            help_breakpoint(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    // Check whether the parameters are correct
    if (vcfFileName.empty())
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
             << "Parameter error.\n";
        help_breakpoint(argv);
        return 1;
    }
    
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running ...\n";

    // init
    VCFBreakpoint VCFBreakpointClass(vcfFileName, prefix, breakpointErrorSize_);
    // set breakpoint
    VCFBreakpointClass.vcf_breakpoint();
    
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done ...\n";

    return 0;
}

// Help document
void help_breakpoint(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -v [options]" << endl
       << "set error for variants breakpoints." << endl
       << endl
       << "required arguments (Note: vcf files must be sorted):" << endl
       << "    -v, --vcf           FILE       vcf file to set breakpoint error" << endl
       << endl
	   << "optional arguments:" << endl
       << "    -p, --prefix        STRING     output prefix [breakpoint]" << endl
       << "    -s, --size          INT        breakpoint error size [1]" << endl
       << endl
       << "    -h, --help                     print this help document" << endl;
}


/**
 * init
 *
 * @param vcfFileName                     tinput VCF  file name
 * @param prefix                          output file prefix
 * @param breakpointErrorSize             breakpoint error
 * 
**/
VCFBreakpoint::VCFBreakpoint(
    const string & vcfFileName, 
    const string & prefix, 
    const int32_t & breakpointErrorSize
) : vcfFileName_(vcfFileName), prefix_(prefix), breakpointErrorSize_(breakpointErrorSize) {}


void VCFBreakpoint::vcf_breakpoint()
{
    // input file stream
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(vcfFileName_);

    // Record the result of the conversion
    string outTxt;

    // ouyput file stream
    SAVE outFile(prefix_ + ".vcf.gz");

    while(VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        if (INFOSTRUCTTMP.line.find("#") != string::npos) {  // Check for comment lines
            string headLine = INFOSTRUCTTMP.line + "\n";
            outFile.save(headLine);
        } else {
            // Use regular expressions to extract length information
            std::regex reg("END=(\\d+)");
            
            // Additive error
            INFOSTRUCTTMP.lineVec[1] = (static_cast<int32_t>(INFOSTRUCTTMP.POS) > breakpointErrorSize_) ? to_string(static_cast<int32_t>(INFOSTRUCTTMP.POS) - breakpointErrorSize_) : "1";

            // end position
            uint32_t refEnd = stoul(INFOSTRUCTTMP.lineVec[1]) + INFOSTRUCTTMP.LEN - 1;
            INFOSTRUCTTMP.lineVec[7] = regex_replace(string(INFOSTRUCTTMP.lineVec[7]), regex(reg), string("END=" + to_string(refEnd)));

            // Adds to the output string
            outTxt += join(INFOSTRUCTTMP.lineVec, "\t") + "\n";

            // save result
            if (outTxt.size() > 10 * 1024 * 1024) {  // It is written every 10Mb
                outFile.save(outTxt);

                // Empty string
                outTxt.clear();
            }
        }
    }

    // Last write
    outFile.save(outTxt);

    // Empty string
    string().swap(outTxt);
}