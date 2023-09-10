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
    if (outputFilename.empty())
    {
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
    // 打开文件
    ifstream dellyFile;
    dellyFile.open(dellyFilename_, ios::in);

    string informations;
    string outTxt;

    string chromosome;
    uint32_t start;
    string refSeq;
    uint32_t end;
    string qrySeq;

    // The regular expression is used to replace END again. If there is [], it cannot be used in the regex library, so delete it first
    std::regex end_reg("END=(\\d+)");
    std::regex reg1("\\[");
    std::regex reg2("\\]");

    smatch end_search; // Look for END= in the informations field
    
    srand((int)time(0));  // Generate random seeds, or replace 0 with NULL.

    if(!dellyFile.is_open()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "'" << dellyFilename_ << "': No such file or directory." << endl;
        exit(1);
    } else {
        while(getline(dellyFile,informations)) {
            if (informations.find("#") != string::npos) {
                std::regex reg_INFO("INFO\t.*");
                outTxt += regex_replace(string(informations), regex(reg_INFO), string("INFO")) + "\n";
            } else {
                // Remove square brackets
                informations = regex_replace(string(informations), regex(reg1), string(""));
                informations = regex_replace(string(informations), regex(reg2), string(""));
                
                chromosome = split(informations, "\t")[0];
                start = stoul(split(informations, "\t")[1]);
                refSeq = split(informations, "\t")[3];
                qrySeq = split(informations, "\t")[4];

                std::regex reg_ref_seq("\t" + refSeq + "\t");
                std::regex reg_qry_seq("\t" + qrySeq + "\t");
                std::regex reg_PASS("PASS.*");

                
                // Look for the END location in the INFO field
                if (regex_search(informations, end_search, end_reg)) {
                    // end = stoi(end_search.format("$1")) - 1;
                    end = stoul(end_search.str(1)) - 1;
                } else {
                    cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: skip ->" << informations << endl;
                    continue;
                }

                // Some of the structural variations delly finds are going to be very, very large, narrow that variation down.
                if ((end - start) > 10000) {
                    end = start + 10000;
                }        
                
                refSeq = FastaMap_[chromosome].substr(start-1, end-start+1);
                
                if (qrySeq == "<INV>") {
                    qrySeq = refSeq;
                    reverse(qrySeq.begin(), qrySeq.end());
                } else if (qrySeq == "<DEL>") {
                    qrySeq = refSeq[0];
                } else if (qrySeq == "<DUP>") {
                    qrySeq = refSeq + refSeq;
                } else if (qrySeq == "<INS>") {
                    uint32_t chromosome_len = FastaMap_[chromosome].length();
                    uint32_t qry_start = rand()%(chromosome_len-10001)+0;
                    uint32_t qryLen = rand()%10000+50;
                    qrySeq = refSeq[0] + FastaMap_[chromosome].substr(qry_start, qryLen);
                } else if (qrySeq.find(chromosome) == string::npos) {  // Translocation.
                
                    if (refSeq.length() == 1) {
                        cerr << "[" << __func__ << "::" << getTime() << "] " << "Skip this short variations: " << chromosome << " " << start << endl;
                        continue;
                    } else {
                        qrySeq = refSeq[0];
                    }
                } else {
                    cerr << "[" << __func__ << "::" << getTime() << "] " << "This variations could not be converted, please check code: " << informations << endl;
                }

                
                informations = regex_replace(string(informations), regex(reg_ref_seq), string("\t" + refSeq + "\t"));
                informations = regex_replace(string(informations), regex(reg_qry_seq), string("\t" + qrySeq + "\t"));
                informations = regex_replace(string(informations), regex(reg_PASS), string("PASS\tEND=" + to_string(end)));

                outTxt += informations + "\n";
            }
        }
    }

    // save result
    ofstream outFile;

    outFile.open(outputFilename_, ios::app);

    if(!outFile.is_open()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "'" << outputFilename_ << "': No such file or directory." << endl;
        exit(1);
    }
    outFile << outTxt;

    // free memory
    dellyFile.close();
    outFile.close();
}