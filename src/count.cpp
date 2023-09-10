#include "../include/count.hpp"

using namespace std;

int main_count(int argc, char** argv)
{
    // Input file
    string FileName1 = "";
    string FileName2 = "";
    string fastaFileName = "";

    // Output file
    string outputFileName_ = "";

    // Input parameter
    int c;
    while (true)
    {
        static const struct option long_options[] = 
        {
            {"input", required_argument, 0, 'i'},
            {"out", required_argument, 0, 'o'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "i:o:h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
        case 'i':
            if (FileName1.empty())
            {
                FileName1 = optarg;
            }
            else
            {
                FileName2 = optarg;
            }
            break;
        case 'o':
            outputFileName_ = optarg;
            break;
        case 'h':
        case '?':
            help_count(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 2)
    {
        help_count(argv);
        return 1;
    }

    // print log
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running ...\n";
    
    // init
    Count CountClass(
        FileName1, 
        FileName2, 
        outputFileName_
    );
    // count
    CountClass.fastq_a_count();
    // save
    CountClass.save_result();

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done ...\n";

    return 0;
}

// help document
void help_count(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -i FILE [options]" << endl
       << "calculate the number of bases and reads of fasta/q" << endl
       << endl
       << "required arguments:" << endl
       << "    -i, --input     FILE    input fastq/a, possibly compressed, two are allowed, one for each mate" << endl
       << endl
	   << "optional arguments:" << endl
       << "    -o, --out       FILE    output file name [stdout]" << endl
       << "    -h, --help              print this help document" << endl;
}

Count::Count(
    const string& inputFileName1, 
    const string& inputFileName2, 
    const string & outputFileName
) : inputFileName1_(inputFileName1), inputFileName2_(inputFileName2), outputFileName_(outputFileName) {}


// Open the fastq/a.g file
void Count::fastq_a_count()
{
    if (inputFileName1_.length() > 0 && inputFileName2_.length() > 0) {  // Double-ended sequencing
        // Input file stream
        gzFile gzfp1 = gzopen(inputFileName1_.c_str(), "rb");
        gzFile gzfp2 = gzopen(inputFileName2_.c_str(), "rb");

        // Open file
        if(!gzfp1 || !gzfp2) {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'" << inputFileName1_ << "' or '" << inputFileName2_
                << "': No such file or directory.\n";
            exit(1);
        }
        else {
            kseq_t *ks1;
            kseq_t *ks2;
            ks1 = kseq_init(gzfp1);
            ks2 = kseq_init(gzfp2);
        
            while(kseq_read(ks1) >= 0 && kseq_read(ks2) >= 0) {
                countOut_.readBase += ks1->seq.l;
                countOut_.readBase += ks2->seq.l;
                countOut_.readNum += 2;
            }

            // Free the memory and close the file
            kseq_destroy(ks1);
            kseq_destroy(ks2);
            gzclose(gzfp1);
            gzclose(gzfp2);
        }
    } else if (inputFileName1_.length() > 0 && inputFileName2_.length() == 0) {  // Single-ended sequencing
        // Input file stream
        gzFile gzfp = gzopen(inputFileName1_.c_str(), "rb");

        // Open file
        if(!gzfp) {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << inputFileName1_
                << "': No such file or directory.\n";
            exit(1);
        } else {
            kseq_t *ks;
            ks = kseq_init(gzfp);
        
            while(kseq_read(ks) >= 0) {
                countOut_.readBase +=  ks->seq.l;
                countOut_.readNum++;
            }

            // Free the memory and close the file
            kseq_destroy(ks);
            gzclose(gzfp);
        }
    }
}

void Count::save_result()
{
    // Output result
    string outTxt = "readBase:" + to_string(countOut_.readBase) + "\n" + 
                    "readNum:" + to_string(countOut_.readNum) + "\n" + 
                    "readLen:" + to_string(countOut_.readBase/countOut_.readNum) + "\n";

    SAVE saveClass(outputFileName_);
    saveClass.save(outTxt);
}
