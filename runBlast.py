import subprocess
from Bio.Seq import Seq
from os.path import split

sequences={} #key=sample_counterForBlastMatch, value=sequence
def main(inputArgs):
    file, minLength, maxLength, minIdentity, wd, proteinSearch=inputArgs
    inSampleGeneCounter=1
    #quickblast: reference, query, prints query
    sequence=subprocess.run(f'blastn -query {file} -task \'megablast\' \
            -max_target_seqs 1000000000 -db {wd}/tempBlastDB/temp \
            -num_threads 1 -evalue 1.0E-5 -word_size 20 \
            -outfmt \"6 delim=  qseqid qstart qend sseqid sstart send pident evalue qseq\"'
            , shell=True, executable="/bin/bash", stdout=subprocess.PIPE)
    blastHits=sequence.stdout.decode().split("\n")
    for blasthit in blastHits:
        blasthit=blasthit.split("\t")
        if len(blasthit)>1 and \
                abs(int(blasthit[4])-int(blasthit[5]))>=minLength and \
                abs(int(blasthit[4])-int(blasthit[5]))<=maxLength and \
                float(blasthit[6])>minIdentity: #otherwise nothing found
            #check orientation of sequence and flip to match reference sequence orientation
            if inSampleGeneCounter==1:
                suffix=""
            else:
                suffix="_"+str(inSampleGeneCounter)
            header=split(file)[1].replace(".fasta","").replace(".fna","")
            if int(blasthit[4])>int(blasthit[5]): #4 and 5 are query coordinates
                #get reverse complement
                temp_sequence = Seq(blasthit[-1].strip())
                sequences[header+suffix]=str(temp_sequence.reverse_complement())
            else:
                sequences[header+suffix]=blasthit[-1].strip()
            inSampleGeneCounter+=1
    return sequences