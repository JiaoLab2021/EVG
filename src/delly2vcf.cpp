#include "../include/delly2vcf.hpp"

using namespace std;

int main_delly2vcf(int argc, char** argv)
{
    string dellyFilename;
    string refFileName;
    string outputFilename;

    // Input parameter
    int c;
    while (true)
    {
        static const struct option long_options[] = 
        {
            {"vcf", required_argument, 0, 'v'},
            {"reference", required_argument, 0, 'r'},
            {"out", no_argument, 0, 'o'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "v:r:o:h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
        case 'v':
            dellyFilename = optarg;
            break;
        case 'r':
            refFileName = optarg;
            break;
        case 'o':
            outputFilename = optarg;
            break;
        case 'h':
        case '?':
            help_delly2vcf(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 2)
    {
        help_delly2vcf(argv);
        return 1;
    }

    // Set the name of the output file
    if (outputFilename.empty()) {
        outputFilename = split(dellyFilename, ".")[0] + ".convert.vcf";
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running ...\n";

    // init
    Delly2VCF Delly2VCFClass(refFileName, dellyFilename, outputFilename);

    // build index
    Delly2VCFClass.build_fasta_index();

    // convert
    Delly2VCFClass.vcf_convert();

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done ...\n";
    return 0;
}

// Help document
void help_delly2vcf(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -v -r [options]" << endl
       << "convert delly-generated VCF file to standard format." << endl
       << endl
       << "required arguments:" << endl
       << "    -v, --vcf         FILE    VCF file output by delly" << endl
       << "    -r, --reference   FILE    reference genome (fasta)" << endl
       << endl
	   << "optional arguments:" << endl
       << "    -o, --out         FILE    output file name [xxx.convert.vcf]" << endl
       << endl
       << "    -h, --help              print this help document" << endl;
}


/**
 * init
 *
 * @param refFileName
 * @param dellyFilename
 * @param outputFilename
 * 
**/
Delly2VCF::Delly2VCF(
    const string& refFileName,
    const string& dellyFilename, 
    const string& outputFilename
) : refFileName_(refFileName), dellyFilename_(dellyFilename), outputFilename_(outputFilename) {}


void Delly2VCF::build_fasta_index()
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Building refFileName_ index: " << refFileName_ << endl;

    uint64_t genomeSize = 0;

    // open fasta file
    gzFile gzfp = gzopen(refFileName_.c_str(), "rb");

    // input file stream
    if(!gzfp) {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
            << "'" << refFileName_ << "': No such file or directory." << endl;
        exit(1);
    }
    else {
        kseq_t *ks;
        ks = kseq_init(gzfp);

        // Define local variables to avoid frequent application and release of memory inside the loop
        string chromosome;
        string sequence;

        while( kseq_read(ks) >= 0 ) {
            // ks->name.s  name
            // ks->seq.s   sequence
            chromosome = ks->name.s;
            sequence = ks->seq.s;

            // record the length of sequence
            genomeSize += ks->seq.l;

            // build fasta index
            FastaMap_.emplace(chromosome, sequence);
        }

        // free memory and close file
        kseq_destroy(ks);
        gzclose(gzfp);
    }
}

void Delly2VCF::vcf_convert()
{
    SAVE SAVEClass(outputFilename_);

    stringstream outStream; // Use stringstream instead of string concatenation
    static const int32_t CACHE_SIZE = 1024 * 1024 * 10; // Cache size is 10mb
    outStream.str().reserve(CACHE_SIZE);

    // The regular expression is used to replace END again. If there is [], it cannot be used in the regex library, so delete it first
    std::regex end_reg("END=(\\d+)");
    std::regex reg1("\\[");
    std::regex reg2("\\]");

    smatch end_search; // Look for END= in the informations field
    
    srand((int)time(0));  // Generate random seeds, or replace 0 with NULL.

    // input file stream
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(dellyFilename_);

    // open file
    ifstream dellyFile;
    dellyFile.open(dellyFilename_, ios::in);

    // If not traversed, continue
    while (VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        // empty line, skip
        if (INFOSTRUCTTMP.line.empty()) {
            continue;
        }

        // commant line
        if (INFOSTRUCTTMP.line.find("#") != string::npos) {
            outStream << INFOSTRUCTTMP.line << endl;
            continue;
        }

        // low qual, skip
        if (INFOSTRUCTTMP.FILTER != "PASS") {
            continue;
        }

        // Remove square brackets
        INFOSTRUCTTMP.ALT = regex_replace(string(INFOSTRUCTTMP.ALT), regex(reg1), string(""));
        INFOSTRUCTTMP.ALT = regex_replace(string(INFOSTRUCTTMP.ALT), regex(reg2), string(""));

        // Look for the END location in the INFO field
        uint32_t refEnd = 0;
        if (regex_search(INFOSTRUCTTMP.INFO, end_search, end_reg)) {
            refEnd = stoul(end_search.str(1));
        } else {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: INFO column is missing END information, skipping this position: " << INFOSTRUCTTMP.line << endl;
            continue;
        }

        uint32_t refLen = refEnd - INFOSTRUCTTMP.POS + 1;
        if (refLen > 0.1 * 1000 * 1000) {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: The variant is too long, skipping this position: " << INFOSTRUCTTMP.line << endl;
            continue;
        }
        
        string refSeq = FastaMap_[INFOSTRUCTTMP.CHROM].substr(INFOSTRUCTTMP.POS - 1, refLen);
        string qrySeq = "";

        if (INFOSTRUCTTMP.ALT == "<INV>") {
            qrySeq = refSeq;
            reverse(qrySeq.begin(), qrySeq.end());
        } else if (INFOSTRUCTTMP.ALT == "<DEL>") {
            qrySeq = refSeq[0];
        } else if (INFOSTRUCTTMP.ALT == "<DUP>") {
            qrySeq = refSeq + refSeq;
        } else if (INFOSTRUCTTMP.ALT == "<INS>") {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: ALT == <INS>, skipping this position: " << INFOSTRUCTTMP.line << endl;
            continue;
        } else if (INFOSTRUCTTMP.ALT.find(":") != string::npos) {  // Translocation.
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: Skipping translocation: " << INFOSTRUCTTMP.line << endl;
            continue;
        } else {
            qrySeq = INFOSTRUCTTMP.ALT;
        }

        INFOSTRUCTTMP.lineVec[3] = refSeq;
        INFOSTRUCTTMP.lineVec[4] = qrySeq;

        outStream << join(INFOSTRUCTTMP.lineVec, "\t") + "\n";;

        if (outStream.tellp() >= CACHE_SIZE) {  // Cache size is 10mb
            string outTxt = outStream.str();
            SAVEClass.save(outTxt);
            // Clear stringstream
            outStream.str(string());
            outStream.clear();
        }
    }

    if (outStream.tellp() > 0) {  //Write for the last time
        string outTxt = outStream.str();
        SAVEClass.save(outTxt);
        // Clear stringstream
        outStream.str(string());
        outStream.clear();
    }
}