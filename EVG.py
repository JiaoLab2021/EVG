#!/usr/bin/env python3

# -*- coding: utf-8 -*-


__data__ = "2024/06/24"
__version__ = "1.2.0"
__author__ = "Zezhen Du"
__email__ = "dzz0539@gmail.com or dzz0539@163.com"

import os
import argparse
import sys
import shutil
import concurrent.futures
import threading
import time
import logging

# code path
code_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(code_dir, "src/"))
import convert, select_software, merge, bwa, GraphTyper2, \
    GraphAligner, BayesTyper, Paragraph, VG_MAP_Giraffe, PanGenie, check_software

# get parser
class MyParser:
    def __init__(self):
        self.flag = "-"*65

        self.setup_logging()
        self.initialize_parser()

        # Parse arguments
        self.args = self.parser.parse_args()
        self.setup_paths()

    def setup_logging(self):
        """Configure the logging"""
        self.logger = logging.getLogger('MyParser')
        # log
        self.logger.error(f"data: {__data__}")
        self.logger.error(f"version: {__version__}")
        self.logger.error(f"author: {__author__}")
        self.logger.error(f"If you encounter any issues related to the code, please don't hesitate to contact us via email at {__email__}.\n")
        formatter = logging.Formatter('[%(asctime)s] %(message)s')
        handler = logging.StreamHandler()  # output to the console
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)

    def initialize_parser(self):
        """Setup command line argument parser"""
        self.parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
        # Input options
        self.input_option = self.parser.add_argument_group(title="required input")
        self.input_option.add_argument("-r", dest="reference", help="input FASTA reference", type=str, required=True)
        self.input_option.add_argument("-v", dest="vcf", help="input the merged vcf file", type=str, required=True)
        self.input_option.add_argument("-s", dest="samples", help="samples file (name read1_path read2_path)", type=str, required=True)
        # Depth option
        self.depth_option = self.parser.add_argument_group(title="depth")
        self.depth_option.add_argument("--depth", dest="depth", help="read depth for genotyping [15]", type=float, default=15)
        # Algorithm option
        self.algorithm_option = self.parser.add_argument_group(title="algorithm")
        self.algorithm_option.add_argument("--software", dest="software", help="genotyping software [auto]", type=str, nargs='+', choices=['VG-MAP', 'VG-Giraffe', 'GraphAligner', 'Paragraph', 'BayesTyper', 'GraphTyper2', 'PanGenie'], default="auto")
        self.algorithm_option.add_argument("--mode", dest="mode", help="software mode [precise]", type=str, choices=['precise', 'fast'], default="precise")
        # Parallel option
        self.parallel_option = self.parser.add_argument_group(title="parallel")
        self.parallel_option.add_argument("--threads", dest="threads", help="number of threads [10]", type=int, default=10)
        self.parallel_option.add_argument("--jobs", dest="jobs", help="run n jobs in parallel [3]", type=int, default=3)
        self.parallel_option.add_argument("--number", dest="number", help="BayesTyper and Paragraph process N samples at once [5]", type=int, default=5)
        # Resume option
        self.resume_option = self.parser.add_argument_group(title="resume")
        self.resume_option.add_argument("--force", action="store_true", help="erase already existing output directory [False]")
        self.resume_option.add_argument("--restart", action="store_true", help="attempt to restart existing workflow based on the state of output files [False]")

    def setup_paths(self):
        """Setup paths based on parsed arguments"""
        # input
        self.reference_file = os.path.abspath(self.args.reference)
        self.vcf_file = os.path.abspath(self.args.vcf)
        self.samples_file = os.path.abspath(self.args.samples)

        # Check if the file exists
        if not os.path.exists(self.reference_file):
            self.logger.error(f"'{self.reference_file}' No such file or directory.")
            sys.exit(1)
        if not os.path.exists(self.vcf_file):
            self.logger.error(f"'{self.vcf_file}' No such file or directory.")
            sys.exit(1)
        if not os.path.exists(self.samples_file):
            self.logger.error(f"'{self.samples_file}' No such file or directory.")
            sys.exit(1)


# EVG
class MyEVG(MyParser):
    def __init__(self):
        super().__init__()

        ## Initialize other attributes
        # Lock
        self.lock = threading.Lock()
        # environment variable
        self.env_path = {'PATH': os.environ.get('PATH')}
        # work directory
        self.base_work_dir = os.getcwd()
        # reads
        self.samples_map = {}  # map<sample_name, list<read_path> >
        self.read_infos_map = {}  # map<sample_name, map<read1: "read1", read2: "read2", depth: "depth", length: "length", bam: "bam"> >
        # reference genome information
        self.refgenome_base = 0
        # Converted file path
        self.convert_file_map = {}  # map<sample_name, reference: "refernece", vcf: "vcf", sv: "sv", region: "region", map_index: "map_index", giraffe_index: "giraffe_index", GraphAligner_index "GraphAligner_index", PanGenie_index: "PanGenie_index">
        # genotype vcf path
        self.bayestyper_vcf_map = {}  # map<sample_name, vcf_path>
        self.paragraph_vcf_map = {}  # map<sample_name, vcf_path>
        self.genotype_vcf_map = {}  # map<sample_name, map<software, vcf_path> >
        # merge vcf path
        self.merge_vcf_map = {}

    def check_softwares(self):
        # check software
        for software in ['graphvcf', 'fastAQ', 'tabix', 'bwa', 'samtools', 'vg', 'GraphAligner', 'paragraph', 'bayesTyper', 'graphtyper', 'PanGenie']:
            stdout, stderr, log_out = check_software.check_software_existence(software, self.env_path)
            if log_out:
                raise SystemExit(log_out.strip())

    def execute_and_check(self, function, *args):
        stdout, stderr, log_out = function(*args)
        if log_out:
            raise SystemExit(log_out.strip())
        return stdout, stderr
    
    # capture subprocess exit status
    def throw_exception(self, name):
        # print log
        log = '%s' % name.__cause__
        raise SystemExit(log)

    # Create a directory
    def makedir(self, path_dir: str):
        """Create a directory with handling for existing directories"""
        if self.args.force and self.args.restart:
            raise ValueError("Both `force` and `restart` arguments were provided. Please provide only one of them.")
        
        if os.path.isdir(path_dir):
            if self.args.force:  # If forced to delete, empty the directory and create a new one
                shutil.rmtree(path_dir)
                os.makedirs(path_dir)
                self.logger.error(f"'{path_dir}' already exists, cleared and recreated.")
            elif self.args.restart:
                self.logger.error(f"'{path_dir}' already exists, using restart, skipping folder clearance.")
            else:  # Print a warning and exit code if not forced to delete
                raise SystemExit(f"'{path_dir}' already exists. Use '--force' to overwrite or '--restart' to restart workflow.")
        else:
            os.makedirs(path_dir)
            self.logger.error(f"Created directory: '{path_dir}'")

    # Parse the samples file
    def get_sample_path(self):
        try:
            with open(self.samples_file, 'r') as f:
                samples_list = f.readlines()
        except IOError as e:
            raise ValueError(f'Error occurred while reading the file {self.samples_file}: {e}')

        if not samples_list:
            raise ValueError(f'Empty file: {self.samples_file}')

        for sample in samples_list:
            if not sample.strip() or "#" in sample or "read1_path" in sample:
                continue

            sample_list = sample.split()
            if len(sample_list) >= 2:
                name = sample_list[0]
                read1_path = os.path.abspath(sample_list[1])
                read2_path = os.path.abspath(sample_list[2]) if len(sample_list) > 2 else ""
                self.samples_map[name] = [read1_path, read2_path]
            else:
                raise ValueError(f"Invalid entry in samples file. Expected 2 or 3 fields, got: {sample}")

    # Convert reads files
    def read_convert(self, work_path: str, sample_name: str, read1: str, read2: str):
        # log
        stdout = stderr = log_out = ""

        # switch paths
        os.chdir(work_path)

        # ################################### Transform the second-generation data ###################################
        # Evaluate fastq files
        if read1:
            stdout, stderr, log_out, read_base, read_num, length = convert.fastq_count(code_dir, self.env_path, read1, read2)
            # Report an error if there is a problem with the exit code
            if log_out:
                return stdout, stderr, log_out

            # downsample
            stdout, stderr, log_out, read1_out, read2_out = convert.downsample(code_dir, read1, read2, read_base, self.refgenome_base, self.args.depth, self.env_path, self.args.restart)
            # Report an error if there is a problem with the exit code
            if log_out:
                return stdout, stderr, log_out

            # actual sequencing depth
            depth = min(int(read_base) / int(self.refgenome_base), self.args.depth)
        else:
            length = 0
            read1_out, read2_out = "", ""
            depth = 0

        # Determine whether the sequencing depth is 0
        if depth == 0:
            log_out = f'Sequencing depth is zero. Please verify the input files: {read1} {read2}'
            return stdout, stderr, log_out

        # Assigned to the total hash table
        with self.lock:  # Use locks to protect the following operations
            self.read_infos_map[sample_name] = {
                "read1": read1_out,
                "read2": read2_out,
                "depth": depth,
                "length": length
            }

        return stdout, stderr, log_out
    
    # Convert reference genome
    def refgenome_convert(self, work_path):
        # switch paths
        os.chdir(work_path)

        # convert
        stdout, stderr, log_out, self.reference_file = convert.convert_reference(self.reference_file, self.env_path, self.args.restart)
        # Report an error if there is a problem with the exit code
        if log_out:
            return stdout, stderr, log_out
        
        # count
        stdout, stderr, log_out, self.refgenome_base = convert.fasta_count(code_dir, self.env_path, self.reference_file)

        return stdout, stderr, log_out
    
    # vcf convert
    def vcf_convert(self, work_path):
        # switch paths
        os.chdir(work_path)

        # convert
        stdout, stderr, log_out, vcf_file, sv_file, region_file = convert.vcf_convert(code_dir, self.reference_file, self.vcf_file, 350, self.args.mode, self.env_path, self.args.restart)
        # Report an error if there is a problem with the exit code
        if log_out:
            return stdout, stderr, log_out

        # bgzip
        if sv_file != vcf_file:  # If the two files are the same, sv_file will not be compressed, otherwise an error will be reported
            stdout, stderr, log_out, sv_file_gz = convert.bgzip_vcf(sv_file, self.env_path, self.args.threads, self.args.restart)
            # Report an error if there is a problem with the exit code
            if log_out:
                return stdout, stderr, log_out
            stdout, stderr, log_out, vcf_file_gz = convert.bgzip_vcf(vcf_file, self.env_path, self.args.threads, self.args.restart)
            # Report an error if there is a problem with the exit code
            if log_out:
                return stdout, stderr, log_out
        else:
            stdout, stderr, log_out, vcf_file_gz = convert.bgzip_vcf(vcf_file, self.env_path, self.args.threads, self.args.restart)
            sv_file_gz = vcf_file_gz
            # Report an error if there is a problem with the exit code
            if log_out:
                return stdout, stderr, log_out

        # The path to store the converted file
        with self.lock:  # Use locks to protect the following operations
            self.convert_file_map["vcf"] = vcf_file_gz
            self.convert_file_map["sv"] = sv_file_gz
            self.convert_file_map["region"] = region_file

        return stdout, stderr, log_out
    
    # build index (bwa, vg, PanGenie)
    def build_indexes(self, work_path):
        stdout = stderr = log_out = ""
        
        # Minimum sequencing depth and sequencing length for obtaining sequencing data
        depth_min = 1000
        read_len_min = 100000
        for _, value in self.read_infos_map.items():
            depth_min = min(depth_min, value["depth"])
            read_len_min = min(read_len_min, value["length"])

        # run select_software
        if isinstance(self.args.software, str):  # Program Automatic Judgment Software
            self.args.software = select_software.main(depth_min, read_len_min, self.refgenome_base)

        # The temporary list is used to determine whether to build an index
        select_software_tmp_list = list(set(self.args.software) & {"VG-MAP", "VG-Giraffe"})
        select_software_tmp_list = ["map" if x == "VG-MAP" else "giraffe" if x == "VG-Giraffe" else x for x in select_software_tmp_list]

        # multi-process process pool
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.args.jobs) as executor:
            # Back to the main working path
            os.chdir(work_path)

            futures = []
            # bwa index
            future = executor.submit(
                bwa.index, 
                work_path, 
                self.reference_file,
                self.env_path, 
                self.args.restart
            )
            # save result
            futures.append(future)
            time.sleep(0.05)
            
            # Map and giraffe
            for software in select_software_tmp_list:
                # switch paths
                os.chdir(work_path)

                index_dir = os.path.join(work_path, software)
                self.makedir(index_dir)

                # submit task
                future = executor.submit(
                    VG_MAP_Giraffe.vg_autoindex,
                    self.reference_file,
                    self.convert_file_map["vcf"],
                    software, 
                    self.env_path, 
                    self.args.threads,
                    index_dir,
                    self.args.restart
                )
                # save result
                futures.append(future)

                time.sleep(0.05)

            # Back to the main working path
            os.chdir(work_path)

            # GraphAligner build index
            if "GraphAligner" in self.args.software:
                # switch paths
                index_dir = os.path.join(work_path, "GraphAligner")
                self.makedir(index_dir)

                # submit task
                future = executor.submit(
                    GraphAligner.vg_index,
                    self.reference_file,
                    self.convert_file_map["vcf"],
                    self.args.threads,
                    index_dir, 
                    self.env_path, 
                    self.args.restart
                )
                # save result
                futures.append(future)

                time.sleep(0.05)

            # return to main path
            os.chdir(work_path)

            # PanGenie build index
            if "PanGenie" in self.args.software:
                # switch paths
                index_dir = os.path.join(work_path, "PanGenie")
                self.makedir(index_dir)

                # submit task
                future = executor.submit(
                    PanGenie.run_index, 
                    self.reference_file, 
                    self.convert_file_map["vcf"],
                    index_dir, 
                    self.env_path, 
                    self.args.threads,
                    self.args.restart
                )
                # save result
                futures.append(future)

                time.sleep(0.05)

            # return to main path
            os.chdir(work_path)

            # Obtain the return value of multithreading. If log_out exists, it indicates that the operation failed and the exit code
            for future in concurrent.futures.as_completed(futures):  # concurrent execution
                try:
                    stdout, stderr, log_out = future.result()  # check return value
                    if log_out:
                        return stdout, stderr, log_out
                except Exception as e:
                    self.logger.error(f"An error occurred: {e}")

        # The path to store the converted file
        with self.lock:  # Use locks to protect the following operations
            self.convert_file_map["map_index"] = os.path.join(work_path, "map")
            self.convert_file_map["giraffe_index"] = os.path.join(work_path, "giraffe")
            self.convert_file_map["GraphAligner_index"] = os.path.join(work_path, "GraphAligner")
            self.convert_file_map["PanGenie_index"] = os.path.join(work_path, "PanGenie")

        return stdout, stderr, log_out


    # Input file format conversion and build index
    def file_convert_index(self):
        # Return to the main working path
        os.chdir(self.base_work_dir)

        # reference and vcf, Convert storage paths and create folders
        convert_dir = os.path.join(self.base_work_dir, "convert")  # Conversion file storage path
        self.makedir(convert_dir)

        ## fastAQ count -i refgenome.fa
        self.logger.error("")
        self.logger.error(f"{self.flag} fastAQ count {self.flag}")
        stdout, stderr, log_out, self.refgenome_base = convert.fasta_count(code_dir, self.env_path, self.reference_file)
        if log_out:
            raise SystemExit(log_out.strip())

        # fastAQ sample
        self.logger.error("")
        self.logger.error(f"{self.flag} fastAQ sample {self.flag}")
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.args.threads) as executor:
            futures = []
            for sample_name, value in self.samples_map.items():
                future = executor.submit(
                    self.read_convert,
                    convert_dir,
                    sample_name,
                    value[0],
                    value[1]
                )
                futures.append(future)

                # Multi-thread interval
                time.sleep(0.05)

            for future in concurrent.futures.as_completed(futures):
                try:
                    stdout, stderr, log_out = future.result()
                    if log_out:
                        raise SystemExit(log_out.strip())
                except Exception as e:
                    self.logger.error(f"An error occurred: {e}")

        ## refgenome and VCF convert
        self.logger.error("")
        self.logger.error(f"{self.flag} Convert reference genome and VCF file, and build index. {self.flag}")
        self.execute_and_check(self.refgenome_convert, convert_dir)
        self.execute_and_check(self.vcf_convert, convert_dir)

        # build index
        self.execute_and_check(self.build_indexes, convert_dir)

    # bwa
    def run_bwa(self):
        # Check whether bwa needs to be run
        if not any(s in self.args.software for s in ["GraphTyper2", "BayesTyper", "Paragraph"]):
            return

        self.logger.error("")
        self.logger.error(f"{self.flag} BWA MEM {self.flag}")

        # directory
        os.chdir(self.base_work_dir)
        bwa_dir = os.path.join(self.base_work_dir, "bwa")
        self.makedir(bwa_dir)

        # Execute bwa mem using ThreadPoolExecutor
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.args.jobs) as executor:
            futures = []
            for sample_name, value in self.read_infos_map.items():
                # Submit to process pool
                future = executor.submit(
                    bwa.mem,
                    bwa_dir, 
                    self.reference_file,
                    self.args.threads,
                    sample_name, 
                    self.env_path, 
                    self.args.restart,
                    value["read1"],
                    value["read2"]
                )
                futures.append(future)
                time.sleep(0.05)

        for future in concurrent.futures.as_completed(futures):
            try:
                stdout, stderr, log_out, sample_name, bam_file = future.result()
                if log_out:
                    raise SystemExit(log_out.strip())
                # Record the bam file path
                self.read_infos_map[sample_name]["bam"] = bam_file
            except Exception as e:
                self.logger.error(f"An error occurred: {e}")

    # Each line runs independently
    def run_genotype(self):
        self.logger.error("")
        self.logger.error(f"{self.flag} Genotyping {self.flag}")

        # directory
        genotype_dir = os.path.join(self.base_work_dir, "genotype")
        self.makedir(genotype_dir)
        os.chdir(genotype_dir)

        # ################################### config ###################################
        # Generate configuration files for GraphTyper2, BayesTyper and Paragraph
        bam2graphtyper_file = ""
        bam2bayestyper_file_list = []
        bam2paragraph_file_list = []
        if "GraphTyper2" in self.args.software:
            bam2graphtyper_file = bwa.bam2graphtyper(self.read_infos_map)
        if "BayesTyper" in self.args.software:
            self.bayestyper_vcf_map, bam2bayestyper_file_list = bwa.bam2bayestyper(self.read_infos_map, self.args.number)
        if "Paragraph" in self.args.software:
            self.paragraph_vcf_map, bam2paragraph_file_list = bwa.bam2paragraph(self.read_infos_map, self.args.number)

        # ################################### genotype ###################################
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.args.jobs) as executor:
            future_to_software_sample = {}  # Used to track the software and sample name corresponding to each future
            for software in self.args.software:
                # working path
                software_work_path = os.path.join(genotype_dir, software)
                self.makedir(software_work_path)

                if software in ["VG-MAP", "VG-Giraffe", "GraphAligner", "PanGenie"]:
                    for sample_name, value in self.read_infos_map.items():
                        if "VG-MAP" == software:
                            future = executor.submit(
                                VG_MAP_Giraffe.main, 
                                software_work_path, 
                                sample_name,
                                value["read1"],
                                value["read2"], 
                                value["depth"],
                                software,
                                self.convert_file_map["map_index"],
                                self.args.threads,
                                self.env_path, 
                                self.args.restart
                            )
                        elif "VG-Giraffe" == software:
                            future = executor.submit(
                                VG_MAP_Giraffe.main,
                                software_work_path, 
                                sample_name,
                                value["read1"],
                                value["read2"], 
                                value["depth"],
                                software,
                                self.convert_file_map["giraffe_index"],
                                self.args.threads,
                                self.env_path, 
                                self.args.restart
                            )
                        elif "GraphAligner" == software:
                            future = executor.submit(
                                GraphAligner.main, 
                                software_work_path, 
                                sample_name,
                                value["read1"],
                                value["depth"],
                                self.convert_file_map["GraphAligner_index"],
                                self.args.threads,
                                self.env_path, 
                                self.args.restart
                            )
                        else:
                            future = executor.submit(
                                PanGenie.main,
                                software_work_path, 
                                sample_name,
                                value["read1"],
                                value["read2"],
                                self.convert_file_map["PanGenie_index"],
                                self.args.threads,
                                self.env_path, 
                                self.args.restart
                            )
                        # record
                        future_to_software_sample[future] = (software, sample_name)
                        time.sleep(0.05)
                elif software == "GraphTyper2":
                    future = executor.submit(
                        GraphTyper2.main, 
                        software_work_path, 
                        self.reference_file, 
                        self.convert_file_map["vcf"], 
                        bam2graphtyper_file, 
                        self.convert_file_map["region"], 
                        self.args.threads,
                        self.env_path, 
                        self.args.restart
                    )
                    # record
                    future_to_software_sample[future] = (software, "")
                elif software == "Paragraph":
                    num = 0
                    for bam2paragraph_file in bam2paragraph_file_list:
                        # Create folder and switch paths
                        paragraph_sample_path = os.path.join(software_work_path, str(num))
                        self.makedir(paragraph_sample_path)
                        # executor
                        future = executor.submit(
                            Paragraph.main,
                            paragraph_sample_path,  
                            self.reference_file, 
                            self.convert_file_map["sv"], 
                            bam2paragraph_file, 
                            self.args.threads,
                            self.env_path, 
                            self.args.restart
                        )
                        # record
                        future_to_software_sample[future] = (software, str(num))
                        num += 1
                elif software == "BayesTyper":
                    num = 0
                    for bam2bayestyper_file in bam2bayestyper_file_list:
                        # Create folder and switch paths
                        bayestyper_sample_path = os.path.join(software_work_path, str(num))
                        self.makedir(bayestyper_sample_path)
                        # executor
                        future = executor.submit(
                            BayesTyper.main,
                            bayestyper_sample_path, 
                            self.reference_file, 
                            self.convert_file_map["vcf"], 
                            bam2bayestyper_file, 
                            self.read_infos_map, 
                            self.args.threads,
                            self.env_path, 
                            self.args.restart
                        )
                        # record
                        future_to_software_sample[future] = (software, str(num))
                        num += 1
                time.sleep(0.05)

            # Get and process the return value
            for future in concurrent.futures.as_completed(future_to_software_sample):
                software, sample_name = future_to_software_sample[future]
                try:
                    stdout, stderr, log_out, vcf_file = future.result()
                    if log_out:
                        raise SystemExit(log_out.strip())
                    # update genotype_vcf_map
                    if software in ["VG-MAP", "VG-Giraffe", "GraphAligner", "PanGenie", "BayesTyper", "Paragraph"]:
                        if software not in self.genotype_vcf_map:
                            self.genotype_vcf_map[software] = {}
                        self.genotype_vcf_map[software][sample_name] = vcf_file
                    elif software == "GraphTyper2":
                        self.genotype_vcf_map[software] = vcf_file
                except Exception as exc:
                    self.logger.error(f'{software} execution for {sample_name} generated an exception: {exc}')

    # graphvcf merge
    def vcf_merge(self):
        self.logger.error("")
        self.logger.error(f"{self.flag} Merge the result {self.flag}")

        # directory
        merge_dir = os.path.join(self.base_work_dir, "merge")
        self.makedir(merge_dir)
        os.chdir(merge_dir)

        with concurrent.futures.ThreadPoolExecutor(max_workers=self.args.jobs) as executor:
            future_to_sample = {}  # Used to track the sample name corresponding to each future

            sample_name_tmp_list = sorted(self.read_infos_map.keys())  # Get all samples_name first
            for sample_name_tmp in sample_name_tmp_list:
                vcf_out_tmp_list = []  # Temporary list for submitting tasks
                software_tmp_list = []  # Temporary list for submitting tasks
                for software in self.args.software:  # software not in results file
                    if software not in self.genotype_vcf_map.keys():
                        self.logger.error(f"Warning: Results for '{software}' are missing for sample '{sample_name_tmp}', so it has been skipped.")
                        continue
                    software_tmp_list.append(software)
                    if isinstance(self.genotype_vcf_map[software], dict):
                        if software == "BayesTyper":
                            if sample_name_tmp not in self.bayestyper_vcf_map.keys():
                                self.logger.error(f"Warning: Results for '{software}' are missing for the sample '{sample_name_tmp}', so it has been skipped.")
                                continue
                            vcf_out_tmp_list.append(self.bayestyper_vcf_map[sample_name_tmp])
                        elif software == "Paragraph":
                            if sample_name_tmp not in self.paragraph_vcf_map.keys():
                                self.logger.error(f"Warning: Results for '{software}' are missing for the sample '{sample_name_tmp}', so it has been skipped.")
                                continue
                            vcf_out_tmp_list.append(self.paragraph_vcf_map[sample_name_tmp])
                        else:
                            if sample_name_tmp not in self.genotype_vcf_map[software].keys():  # sample_name is not in the result file
                                self.logger.error(f"Warning: Results for '{software}' are missing for the sample '{sample_name_tmp}', so it has been skipped.")
                                continue
                            vcf_out_tmp_list.append(self.genotype_vcf_map[software][sample_name_tmp])
                    else:
                        vcf_out_tmp_list.append(self.genotype_vcf_map[software])

                # Submit the merge task to the thread pool
                future = executor.submit(
                    merge.main,
                    code_dir,
                    merge_dir,
                    self.convert_file_map["vcf"],
                    sample_name_tmp,
                    vcf_out_tmp_list,
                    software_tmp_list,
                    self.env_path,
                    self.args.restart
                )
                future_to_sample[future] = sample_name_tmp
                time.sleep(0.05)

        # return
        for future in concurrent.futures.as_completed(future_to_sample):
            sample_name = future_to_sample[future]
            try:
                stdout, stderr, log_out, merge_vcf_file = future.result()
                if log_out:
                    raise SystemExit(log_out.strip())
                self.merge_vcf_map[sample_name] = merge_vcf_file
            except Exception as exc:
                self.logger.error(f'Merge execution for {sample_name} generated an exception: {exc}')

    # Print the result of merge
    def print_result(self):
        self.logger.error("")
        self.logger.error(f"{self.flag} Result {self.flag}")
        for key, value in self.merge_vcf_map.items():
            self.logger.error( f"{key}: {value}")


def main():
    # EVG
    EVGClass = MyEVG()
    # check
    EVGClass.check_softwares()
    # Get sequencing file path and output pathname
    EVGClass.get_sample_path()
    # convert
    EVGClass.file_convert_index()
    # bwa
    EVGClass.run_bwa()
    # genotype
    EVGClass.run_genotype()
    # merge
    EVGClass.vcf_merge()
    # print result
    EVGClass.print_result()


if __name__ == '__main__':
    main()
