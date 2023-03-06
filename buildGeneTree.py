#take all fastas in specified directory and check all for specific gene
from os import listdir, getcwd, makedirs, chdir, path
from os.path import isfile, join, splitext, exists
import shutil
import subprocess
import sys
import runBlast
import multiprocessing as mp
from Bio import SeqIO
import argparse

# maxLengthDifference=0.1
# refFasta=sys.argv[1]
# fastaDir=sys.argv[2]
# minIdentity=70
# proteinSearch=False

##### Args parser #######
parser = argparse.ArgumentParser(description='Take nucleotide/amino acid sequence, search files for homologues, build a phylo tree')
parser.add_argument('-r','--ref_sequence', help="Nucleotide/amino acid sequence homologues of which will be used to build a tree", required=True)
parser.add_argument('-d', '--fasta_dir', help="Location of fasta files to search for homologues", required=True)
parser.add_argument('-f', '--file_to_search', default=[], help="By default, all .fna and .fasta file in fasta_dir are search, but a subset can be specified by giving a plain text file", required=False)
parser.add_argument('-l', '--max_length_diff', default=10, help="Maximum difference in length between supplied gene/protein sequence and a BLASTn hit. Default is 10 (i.e +/-10%)", required=False)
parser.add_argument('-i', '--min_identity', default=70, help="Minimum BLASTn identity between supplied gene/protein sequence and a BLASTn hit. Default is 70 i.e. 70%", required=False)
parser.add_argument('-n', '--nucleotide_search', default=True, help="Not currently supported. Reference and targets are nucleotide (True) or proteins (False)", required=False)
parser.add_argument('-o', '--output_dir', default=True, help="Directory in which output will be stored", required=True)
parser.add_argument('-t', '--cpu_threads', default=min(mp.cpu_count(),20), help="Number of CPUs to use for BLASTn search, default is lower of machine total and 20", required=False)

args = parser.parse_args()
if isfile(join(getcwd(), args.ref_sequence)):
    ref_fasta=getcwd()+"/"+args.ref_sequence
else:
    ref_fasta=args.ref_sequence
fasta_dir=args.fasta_dir
wd=args.output_dir

cpu_threads=args.cpu_threads
if not exists(args.fasta_dir):
    print(f'The fasta directory {args.fasta_dir} does not exits')
    sys.exit()

files = [f for f in listdir(args.fasta_dir) if isfile(join(args.fasta_dir, f)) and (splitext(f)[-1]==".fasta" or splitext(f)[-1]==".fna")]

if len(files)==0:
    print(f'No files ending in .fna or .fasta found in fasta directory {args.fasta_dir}')
    sys.exit()

if "file_to_search" in args and len(args.file_to_search)>0: #remove from files those that are not present in file_to_search
    files_to_search=[f.replace(".fasta","").replace(".fna","") for f in args.file_to_search]
    files=[f for f in files if f.replace(".fasta","").replace(".fna","") in files_to_search]
    if len(files)==0:
        print(f'No files from file_to_search are present in fasta directory {args.fasta_dir}')
        sys.exit()
min_identity=args.min_identity/100
max_length_difference=args.max_length_diff/100

proteinSearch=False

##### Args parser END #######

# SamplesToExclude=[]
# with open("/mnt/storage5/anton/MGE/scripts/buildGeneTree/excludeSamples.txt") as file:
#     for line in file:
#         SamplesToExclude.append(line.strip())

#files = [f for f in listdir(fastaDir) if isfile(join(fastaDir, f)) and (splitext(f)[1]==".fasta" or splitext(f)[1]==".fna")]
# with open("/mnt/storage5/anton/KpReplicons/selectedAssemblies.tsv") as file:
#     files=[]
#     for line in file:
#         files.append(line.strip())

if not exists(wd):
    makedirs(wd)
chdir(wd)



#create blast DB
if exists("./tempBlastDB/"):
    shutil.rmtree("./tempBlastDB/")

makedirs("./tempBlastDB/")
subprocess.call("makeblastdb -in "+ref_fasta+" -title temp -out "+"./tempBlastDB/temp -dbtype nucl \
      -blastdb_version 4 1>/dev/null", shell=True)

ref_fasta_length=0
for i, record in enumerate(SeqIO.parse(ref_fasta, "fasta")): #get reference lengh
    ref_fasta_length=len(str(record.seq))
    if i>0:
        print(f'Reference file {ref_fasta} contains more then one sequence. This could be an error')

sequences={} #key=sample_counterForBlastMatch, value=sequence

i=0
while i<len(files):
    # if files[i].replace(".fna","").replace(".fasta","") in SamplesToExclude:
    #     files.remove(files[i])
    # else:
    files[i]=[files[i], fasta_dir, int(ref_fasta_length*(1-max_length_difference)), \
            int(ref_fasta_length*(1+max_length_difference)), min_identity, proteinSearch]
    i+=1

#sys.exit()
print("Collecting gene sequences")
progressCounter=0
pool=mp.Pool(cpu_threads)
results=pool.map(runBlast.main, files)
pool.close()
pool.join()
for file in results:
    for sequence in file:
        sequences[sequence]=file[sequence]
print("Completed collecting gene sequences")

#output hits
# if exists(wd+"/treetemp/"):
#     shutil.rmtree(wd+"/treetemp/")
# if not exists(wd+"/treetemp/"):
#     makedirs(wd+"/treetemp/")
# chdir(wd+"/treetemp/") 
# outputFile=wd+"unaligned.fasta"
with open("unaligned.fasta", "w") as output:
    for sequence in sequences:
        output.write(">"+sequence+"\n")
        output.write(sequences[sequence]+"\n")

print("Aligning")
head, tail = path.split(ref_fasta)
if proteinSearch:
    subprocess.call(f'probcons unaligned.fasta > {tail.replace(".fasta","")}_aligned.fasta', shell=True)
else:
    subprocess.call(f'mafft --auto unaligned.fasta > {tail.replace(".fasta","")}_aligned.fasta', shell=True)
print("Finished aligning")

print("Building tree")
#subprocess.call("iqtree -nt AUTO -pre treetemp -s aligned.fasta -bb 1000", shell=True)
#subprocess.call("iqtree -nt AUTO -m GTR+F+I+G4 -pre "+tail.replace(".fasta","")+" -s "+tail.replace(".fasta","")+"_aligned.fasta -bb 1000", shell=True)
print("Done")