#ifndef SAMPLE_HPP
#define SAMPLE_HPP

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <getopt.h>

#include "count.hpp"
#include "get_time.hpp"
#include "strip_split_join.hpp"

using namespace std;

namespace sample
{
    void save(string outputFileName, string outTxt);

    // Open the fastq.gz file
    void sample(string inputFileName1, string inputFileName2, const vector<long long int> & randVecSel, string prefix)
    {
        if (inputFileName1.size() > 0 && inputFileName2.size() > 0) // Double-ended sequencing
        {
            // Output file name
            string outputFileName1 = prefix + "." + split(inputFileName1, "/")[split(inputFileName1, "/").size() - 1];
            string outputFileName2 = prefix + "." + split(inputFileName2, "/")[split(inputFileName2, "/").size() - 1];
            if (outputFileName1.rfind(".gz") == string::npos && outputFileName1.rfind(".GZ") == string::npos)
            {
                outputFileName1 += ".gz";
            }
            if (outputFileName2.rfind(".gz") == string::npos && outputFileName2.rfind(".GZ") == string::npos)
            {
                outputFileName2 += ".gz";
            }

            // Check whether the file exists. If yes, delete the file
            // Output file stream
            gzFile gzfpO1 = gzopen(outputFileName1.c_str(), "wb");
            gzFile gzfpO2 = gzopen(outputFileName2.c_str(), "wb");
            gzclose(gzfpO1);
            gzclose(gzfpO2);
            
            // Input file stream
            gzFile gzfpI1 = gzopen(inputFileName1.c_str(), "rb");
            gzFile gzfpI2 = gzopen(inputFileName2.c_str(), "rb");

            // Record read
            long long int readNum = 0;

            // Output string
            string outTxt1;
            string outTxt2;

            // The index recorded in the vec
            long long int index = 0;

            // Open file
            if(!gzfpI1 || !gzfpI2)
            {
                cerr << "[" << __func__ << "::" << getTime() << "] '" 
                     << inputFileName1 << "' or '" << inputFileName2
                     << "': No such file or directory.\n";
                exit(1);
            }
            else
            {
                kseq_t *ks1;
                kseq_t *ks2;
                ks1 = kseq_init(gzfpI1);
                ks2 = kseq_init(gzfpI2);
            
                while(kseq_read(ks1) >= 0 && kseq_read(ks2) >= 0)
                {
                    if (randVecSel[index] == readNum) // If it is a selected read, it is added to outTxt
                    {
                        string readName1 = ks1->name.s;
                        readName1 = "@" + readName1;
                        string readSeq1 = ks1->seq.s;
                        string readQual1 = ks1->qual.s;

                        string readName2 = ks2->name.s;
                        readName2 = "@" + readName2;
                        string readSeq2 = ks2->seq.s;
                        string readQual2 = ks2->qual.s;

                        // Check whether fastq has joint sequence information
                        string readNameLast1;
                        string readNameLast2;
                        if (ks1->comment.s)
                        {
                            readNameLast1 = ' ' + string(ks1->comment.s);
                            readName1 += readNameLast1;
                        }
                        else
                        {
                            readNameLast1 = "";
                            readName1 += " 1";
                        }

                        if (ks2->comment.s)
                        {
                            readNameLast2 = ' ' + string(ks2->comment.s);
                            readName2 += readNameLast2;
                        }
                        else
                        {
                            readNameLast2 = "";
                            readName2 += " 1";
                        }

                        outTxt1 += readName1 + 
                                readNameLast1 + 
                                "\n" + 
                                readSeq1 + 
                                "\n+\n" + 
                                readQual1 +
                                "\n";

                        outTxt2 += readName2 + 
                                readNameLast2 + 
                                "\n" + 
                                readSeq2 + 
                                "\n+\n" + 
                                readQual2 +
                                "\n";

                        if (outTxt1.length() > 100000000) // Write once every 100Mb
                        {
                            // Write to disk
                            save(outputFileName1, outTxt1);
                            save(outputFileName2, outTxt2);

                            // Empty string
                            outTxt1.clear();
                            outTxt2.clear();
                            string().swap(outTxt1);
                            string().swap(outTxt2);
                        }

                        // index points to the next element
                        index++;
                    }
                    // read count + 1
                    readNum++;
                }

                // Writes to disk
                save(outputFileName1, outTxt1);
                save(outputFileName2, outTxt2);

                // Empty string
                outTxt1.clear();
                outTxt2.clear();
                string().swap(outTxt1);
                string().swap(outTxt2);

                // Free the memory and close the file
                kseq_destroy(ks1);
                kseq_destroy(ks2);
                gzclose(gzfpI1);
                gzclose(gzfpI2);
            }
        }
        else if (inputFileName1.size() > 0 && inputFileName2.size() == 0) // Single-ended sequencing
        {
            // Output file name
            string outputFileName1 = prefix + "." + split(inputFileName1, "/")[split(inputFileName1, "/").size() - 1];
            if (outputFileName1.rfind(".gz") == string::npos && outputFileName1.rfind(".GZ") == string::npos)
            {
                outputFileName1 += ".gz";
            }

            // Check whether the file exists. If yes, delete the file
            // Output file stream
            gzFile gzfpO1 = gzopen(outputFileName1.c_str(), "wb");
            gzclose(gzfpO1);
            
            // Input file stream
            gzFile gzfpI1 = gzopen(inputFileName1.c_str(), "rb");

            // Record read
            long long int readNum = 0;

            // Output string
            string outTxt1;

            // The index recorded in the vec
            long long int index = 0;

            // Open file
            if(!gzfpI1)
            {
                cerr << "[" << __func__ << "::" << getTime() << "] '" 
                     << inputFileName1
                     << "': No such file or directory.\n";
                    << endl;
                exit(1);
            }
            else
            {
                kseq_t *ks1;
                ks1 = kseq_init(gzfpI1);
            
                while(kseq_read(ks1) >= 0)
                {
                    if (randVecSel[index] == readNum) // If it is a selected read, it is added to outTxt
                    {
                        string readName1 = ks1->name.s;
                        readName1 = "@" + readName1;
                        string readSeq1 = ks1->seq.s;
                        string readQual1 = ks1->qual.s;

                        // Check whether fastq has joint sequence information
                        string readNameLast1;
                        if (ks1->comment.s)
                        {
                            readNameLast1 = ' ' + string(ks1->comment.s);
                            readName1 += readNameLast1;
                        }
                        else
                        {
                            readNameLast1 = "";
                            readName1 += " 1";
                        }

                        outTxt1 += readName1 + 
                                readNameLast1 + 
                                "\n" + 
                                readSeq1 + 
                                "\n+\n" + 
                                readQual1 +
                                "\n";

                        if (outTxt1.length() > 100000000) // Write once every 100Mb
                        {
                            // Write to disk
                            save(outputFileName1, outTxt1);

                            // Empty string
                            outTxt1.clear();
                            string().swap(outTxt1);
                        }

                        // index points to the next element
                        index++;
                    }
                    // read count +1
                    readNum++;
                }

                // Writes to disk
                save(outputFileName1, outTxt1);

                // Empty string
                outTxt1.clear();
                string().swap(outTxt1);

                // Free the memory and close the file
                kseq_destroy(ks1);
                gzclose(gzfpI1);
            }
        }
        else
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                    << " please input sequencing file." 
                    << endl;
            exit(1);
        }
    }

    // Save the fastq.gz file
    void save(string outputFileName, string outTxt)
    {
        // Output file stream
        gzFile gzfp = gzopen(outputFileName.c_str(), "ab");

        // Open file
        if(!gzfp)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] '" 
                 << outputFileName
                 << "': No such file or directory.\n";
            exit(1);
        }
        else
        {
            gzwrite(gzfp, outTxt.c_str(), outTxt.length());
        }

        // Free the memory and close the file
        gzclose(gzfp);
    }

    // A random set of data
    vector<long long int> randVector(long long int num)
    {
        vector<long long int> result;
        result.clear();
        result.reserve(num);
        srand((long long int)time(0));
        for (size_t i = 0; i < num; i++)
        {
            result.push_back(i);
        }
        long long int p1;
        long long int p2;
        long long int temp;
        long long int count = num;

        while (--num)
        {
            p1 = num;
            p2 = rand() % num;
            temp = result[p1];
            result[p1] = result[p2];
            result[p2] = temp;
        }
        return result;
    }
}

#endif