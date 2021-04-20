#take all fastas in specified directory and check all for specific gene
from os import listdir, getcwd, makedirs, remove, chdir
from os.path import isfile, join, splitext, exists
import shutil
import pandas as pd
import subprocess
import sys
import runBlast
import multiprocessing as mp
from Bio import SeqIO

maxLengthDifference=0.2
wd=getcwd()
chdir(wd)
SamplesToExclude=["3042","5145","Kp_4240","Kp4151", "3766", "5170","Kp4189","Kp4173","3763","Kp4179","Kp931"]
refFasta=sys.argv[1]
fastaDir=sys.argv[2]
minIdentity=80

files = [f for f in listdir(fastaDir) if isfile(join(fastaDir, f)) and (splitext(f)[1]==".fasta" or splitext(f)[1]==".fna")]

#create blast DB
if exists(wd+"/tempBlastDB/"):
    shutil.rmtree(wd+"/tempBlastDB/")
makedirs(wd+"/tempBlastDB/")
subprocess.call("makeblastdb -in "+refFasta+" -title temp -out "+wd+"/tempBlastDB/temp -dbtype nucl \
      -blastdb_version 4 1>/dev/null", shell=True)

refFastaLength=0
for record in SeqIO.parse(refFasta, "fasta"): #get reference lengh
    refFastaLength=len(str(record.seq))

sequences={} #key=sample_counterForBlastMatch, value=sequence
for i in range(0,len(files)):
    files[i]=[files[i], fastaDir, int(refFastaLength*(1-maxLengthDifference)), \
            int(refFastaLength*(1+maxLengthDifference)), minIdentity]

print("Collecting gene sequences")
progressCounter=0
pool=mp.Pool(20)
results=pool.map(runBlast.main, files)
pool.close()
pool.join()
for file in results:
    for sequence in file:
        sequences[sequence]=file[sequence]
print("Completed collecting gene sequences")

#output hits
if exists(wd+"/treetemp/"):
    shutil.rmtree(wd+"/treetemp/")
makedirs(wd+"/treetemp/")
chdir(wd+"/treetemp/") 
outputFile=wd+"/treetemp/"+"unaligned.fasta"
with open(outputFile, "w") as output:
    for sequence in sequences:
        output.write(">"+sequence+"\n")
        output.write(sequences[sequence]+"\n")
print("Aligning")
subprocess.call("mafft --auto "+wd+"/treetemp/unaligned.fasta > "+wd+"/treetemp/aligned.fasta", shell=True)
print("Finished aligning")

print("Building tree")
subprocess.call("iqtree -nt AUTO -pre treetemp -s aligned.fasta -bb 1000", shell=True)
print("Done")