import subprocess
from Bio.Seq import Seq
from os import getcwd

sequences={} #key=sample_counterForBlastMatch, value=sequence
def main(inputArgs):
    file, fastaDir, minLength, maxLength, minIdentity=inputArgs
    inSampleGeneCounter=1
    #quickblast: reference, query, prints query
    sequence=subprocess.run("blastn -query "+fastaDir+file+" -task 'megablast' \
            -max_target_seqs 1000000000 -db "+getcwd()+"/tempBlastDB/temp \
            -num_threads 8 -evalue 1.0E-5 -word_size 20 \
            -outfmt \"6 delim=  qseqid qstart qend sseqid sstart send pident evalue qseq\""
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
            if int(blasthit[4])>int(blasthit[5]): #4 and 5 are query coordinates
                #get reverse complement
                temp_sequence = Seq(blasthit[-1].strip())
                sequences[file.replace(".fasta","")+suffix]=str(temp_sequence.reverse_complement())
            else:
                sequences[file.replace(".fasta","")+suffix]=blasthit[-1].strip()
            inSampleGeneCounter+=1
    return sequences