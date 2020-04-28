#coding:utf-8

"""
#==========================================================
Contact: Jie Li
v1.0	20181207	Use this to download sra data from NCBI with sratoolkit
#==========================================================
"""

import os
import sys
import argparse as ap
import pandas as pd
import basictools as basic
import pipelinetools as pipe

def read_paras(args):
	p = ap.ArgumentParser(description=__doc__)
	ars = p.add_argument
	ars("--accession", dest="accession", required=True,
		help="accession number of data you want to download")
	ars("--list", dest="list", action="store_true",
		help="set if your --accession input is a list of acc nums.")
	ars("--rm", dest="rm", action="store_true",
		help="set if you want to delete all *.sh.* files")
	ars("--outdir", dest="outdir", default=os.getcwd(),
		help="output dir to output data")
	return vars(p.parse_args())

def download_sra(result, shell, sra, softwares, rm, rawdata):	
	rawdata.write(sra + "\t" + result + sra+"/"+sra+"_1.fastq.gz\t"+result + sra+"/"+sra+"_2.fastq.gz\n")
	shell = shell + "sra_download_" + sra + ".sh"
	result = result + sra + "/"
	basic.mkdir([result])
	sl = softwares['prefetch'] + " --min-size 1 --max-size 100000000 " + sra + " &&\n"
	sl += "mv " + softwares['cache'] + "sra/"+ sra + ".sra " + result + " &&\n"
	sl += softwares['fastq-dump']+" --helicos -I -O "+result+" --gzip "+"--split-files "+result+sra +".sra &&\n"
	sl += softwares['fastqc']+" -o "+result+" "+result+sra +"_1.fastq.gz "+result+sra+"_2.fastq.gz &&\n"
	sl += "ls " + result + "*_fastqc.zip|xargs -t -i unzip -u -d " + result + " {} &&\n"
	for i in ['_1_', '_2_']:
		sl += "cp " + result + sra + i + "fastqc/Images/per_base_quality.png " +\
			  result + sra + i + "per_base_quality.png &&\n"
		sl += "cp " + result + sra + i + "fastqc/Images/duplication_levels.png " +\
			  result + sra + i + "duplication_levels.png &&\n"
		sl += "cp " + result + sra + i + "fastqc/Images/per_base_sequence_content.png " +\
			  result + sra + i + "per_base_sequence_content.png &&\n"
	pipe.creat_shell(sl, shell, rm, "select=1:ncpus=2:mem=2GB", "02:99:00")

if __name__ == "__main__":
	pars = read_paras(sys.argv)
	outdir = basic.outdir(pars['outdir']) + "/"
	shell = outdir + "shell/SRA_download/"
	result = outdir + "result/SRA_download/"
	basic.mkdir([outdir, outdir+"shell", outdir +"result", shell, result])
	softwares = basic.parse_config(os.path.dirname(sys.argv[0]) +"/softwares.txt")
	sras = pd.read_csv(pars['accession'], header=None)[0] if pars['list'] else [pars['accession']]
	with open(outdir + "result/SRA_download/rawdata.list", "w") as rawdata:
		rawdata.write("#sample\tFQ1\tFQ2\n")
		for sra in sras:
			download_sra(result, shell, sra, softwares, pars['rm'], rawdata)
	sl = softwares['python36'] + " " + \
		 os.path.dirname(sys.argv[0]) + "/stat_fastqc_result-v1.0.py" +\
		 " --fqlist " + result + "rawdata.list" +\
		 " --FastQC_dir " + result +\
		 " --outfile " + result + "stats_of_FastQC_result.xls &&\n" 
	pipe.creat_shell(sl, shell + "stat_FastQC.sh", pars['rm'], "select=1:ncpus=4:mem=5GB", "00:99:00")




