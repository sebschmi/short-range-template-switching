import numpy as np
import pandas as pd
from itertools import compress
from pyfaidx import Faidx
import subprocess
from io import StringIO
     
import sys, getopt, os

######################################################

def shift(arr, num, fill_value=np.nan):
    result = np.empty_like(arr)
    if num > 0:
        result[:num] = fill_value
        result[num:] = arr[:-num]
    elif num < 0:
        result[num:] = fill_value
        result[:num] = arr[-num:]
    else:
        result[:] = arr
    return result

######################################################


def main(argv):

    create = False
    align = False
    score = False
    fafile = ''
    vcffile = ''
    locifile = ''
    bamfile = ''
    alleledir = ''
    bamdir = ''
    indid = ''
    untrimmed = False
    trimlen = 50
    oldformat = False
    oldindex = False
    iscram = False

    cmdopts = 'tsm_alleles.py -c -a -s -v <vcffile> -l <locifile> -r <readbam> -d <alleledir> -b <bamdir> -u -t <trimlen> -i <individual> -m [-o -x]'

    try:
        opts, args = getopt.getopt(argv,"hcasf:v:l:r:d:b:ut:i:oxm",["create","align","score","fafile=","vcffile=","locifile=","readbam=","alleledir=","bamdir=","untrimmed","trimlen=","indiv=","oldformat","oldindex","cram"])
        
    except getopt.GetoptError:
        print (cmdopts)
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-h':
            print (cmdopts)
            sys.exit()
        elif opt in ("-c", "--create"):
            create = True
        elif opt in ("-a", "--align"):
            align = True
        elif opt in ("-s", "--score"):
            score = True
        elif opt in ("-f", "--fafile"):
            fafile = arg
        elif opt in ("-v", "--vcffile"):
            vcffile = arg
        elif opt in ("-l", "--locifile"):
            locifile = arg
        elif opt in ("-r", "--readbam"):
            bamfile = arg
        elif opt in ("-d", "--alleledir"):
            alleledir = arg
        elif opt in ("-b", "--bamdir"):
            bamdir = arg
        elif opt in ("-u", "--untrimmed"):
            untrimmed = True
        elif opt in ("-t", "--trimlen"):
            trimlen = int(arg)
        elif opt in ("-i", "--indiv"):
            indid = arg
        elif opt in ("-o", "--oldformat"):
            oldformat = True
        elif opt in ("-x", "--oldindex"):
            oldindex = True
        elif opt in ("-m", "--cram"):
            iscram = True
    
    ######################################################

    max_dist = 50
    seq_offs = 500

    ######################################################

    # condition that alleles are created
    if create:

        if(fafile == '' or vcffile == '' or locifile == '' or alleledir == '' or indid == ''):
            print (cmdopts)
            sys.exit()

        fasta = Faidx(fafile)
        
        # loop for loci
        with open(locifile, 'r') as f:

            # loop for loci
            for line in f:
                
                chrm,start,stop,ind = line.strip().split("\t")
                ind = int(ind)

                cmd = "bcftools view -s "+indid+" -r "+chrm+":"+start+"-"+stop+" "+vcffile+" | grep -v '0|0' | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\t%INFO/VARTYPE\n'"    
                gts = subprocess.check_output(cmd, shell=True)
        
                data = pd.read_csv(StringIO(gts.decode('UTF-8')), sep="\t|\|", dtype=object, header=None, engine='python')
        
                chrm = np.array(data.iloc[:,0])
                strt = np.array(data.iloc[:,1].astype(int))
                refc = np.array(data.iloc[:,2])
                altc = np.array(data.iloc[:,3])
                stop = strt + [len(s)-1 for s in refc]
                svtp = np.array(data.iloc[:,6])

                geno = np.empty((data.shape[0],2), dtype=bool)
                for i in range(4,6):
                    geno[:,(i-4)] = np.array(data.iloc[:,i] == "1")

    

                # skip the rest if there are no cases
                if sum(geno[:,ind])==0:
                    break
                
                chm = chrm[geno[:,ind]]
                sta = strt[geno[:,ind]]
                sto = stop[geno[:,ind]]
                ref = refc[geno[:,ind]]
                alt = altc[geno[:,ind]]
                svt = svtp[geno[:,ind]]

                ic = chm[0]


                keep = np.logical_and(np.equal(chm,ic), sta-np.array(shift(sto,1,-1000)) < max_dist)

                #print(ind,ic,len(keep))
                #print(sta[keep])

                # loop over the positions to find clusters
                for p in [i for i, x in enumerate(keep) if x]:
                        
                    # condition that the position has variant
                    if keep[p]:

                        clst = []
                        if(chm[p-1] == ic):
                            clst = [p-1]

                        j=0
                        while p+j<len(keep) and keep[p+j]:
                            clst.append(p+j)
                            keep[p+j] = False
                            j+=1

                        #print("sample",ind,chm[p],sta[clst[0]],sta[clst[-1]],"[",len(clst),"]", flush= True)

                        seq_offs_extra = 0

                        for c in clst:
                            seq_offs_extra +=  len(ref[c])

                        seq = str(fasta.fetch(ic,sta[clst[0]]-seq_offs,sta[clst[-1]]+seq_offs+seq_offs_extra))

                        seq_st = sta[clst[0]]-seq_offs
                        seq_po = 0
                        ref_seq = str()
                        alt_seq = str()

                        # loop over the variants in a cluster
                        for c in clst:

                            ref_al = ref[c]                    
                            alt_al = alt[c]

                            if not untrimmed: 
                                if(len(ref_al)>2*trimlen):
                                    ref_al = ref[c][0:trimlen]+"NNNNN"+ref[c][-trimlen:] 

                                if(len(alt_al)>2*trimlen):
                                    alt_al = alt[c][0:trimlen]+"NNNNN"+ref[c][-trimlen:] 

                            ld = len(ref_al)-len(alt_al)

                            ref_seq += seq[seq_po:(sta[c]-seq_st)] #+ "|"
                            ref_seq += ref_al 

                            if ld < 0:
                                ref_seq += "-"*(-1*ld)
                            #ref_seq += "|"

                            alt_seq += seq[seq_po:(sta[c]-seq_st)] #+ "|"                    
                            alt_seq += alt_al

                            if ld > 0:
                                alt_seq += "-"*(ld)
                            #alt_seq += "|"
                                
                            seq_po = sta[c]-seq_st+len(ref[c])
                                
                        # loop over the variants in a cluster
                            
                        ref_seq += seq[seq_po:]
                        alt_seq += seq[seq_po:]


                        locname = str(chm[p])+"_"+str(sta[clst[0]])+"_"+str(sta[clst[-1]])+".fa"

                        with open(alleledir+"/"+locname, 'w') as f:
                            f.write(">REF\n"+ref_seq)
                            f.write("\n>ALT\n"+alt_seq)
                                
                            
                        cmd = "bwa index " + alleledir+"/"+locname + " 2> /dev/null"
                        os.system(cmd)
                            
                    # condition that the position has variant
                # loop over the positions to find clusters
            # loop for loci
        # loop for loci
    # condition that alleles are created
  
    # condition that reads are mapped
    if align:
        
        if(bamfile == '' or bamdir == '' or alleledir == ''):
            print (cmdopts)
            sys.exit()
        
        files = []
        files += [each for each in os.listdir(alleledir) if each.endswith('.fa')]

        # loop over allele files
        for file in files: 
            chrm,start,stop = None, None, None
            if oldformat:
                chrm,start,p2,p3,stop,suf = file.split(".")
            else:
                chrm,start,stop = file.split("_")
                stop = stop[:-3]
            file2 = file[:-3]+".fq"
            file3 = file[:-3]+".bam"

            startp = str(int(start)-150)
            stopp = str(int(stop)+150)
            
            cmd= "samtools view -h "+bamfile+" "+chrm+":"+startp+"-"+stopp+" | samtools sort -n - | bedtools bamtofastq -i - -fq "+bamdir+"/"+file2 
            if iscram:
                cmd= "samtools view -h -T "+fafile+" "+bamfile+" "+chrm+":"+startp+"-"+stopp+" | samtools sort -n - | bedtools bamtofastq -i - -fq "+bamdir+"/"+file2 
            #print(cmd)
            os.system(cmd)

            if(oldindex):
                cmd = "bwa index "+alleledir+"/"+file+" "+bamdir+"/"+file2+" 2> /dev/null"
                os.system(cmd)
                    
            cmd = "bwa mem "+alleledir+"/"+file+" "+bamdir+"/"+file2+" 2> /dev/null | " + \
                    "awk '(length($10) > 74 || $1 ~ /^@/ )' | samtools view -F 4 -b | samtools sort - > "+ bamdir+"/"+file3 + "&&" + \
                    " samtools index "+bamdir+"/"+file3
            os.system(cmd)
            
        # loop over allele files
    # condition that reads are mapped

    # condition that bams are scored
    if score:

        # loop over allele files
        if(bamdir == ''):
            print (cmdopts)
            sys.exit()

        files = []
        files += [each for each in os.listdir(bamdir) if each.endswith('.bam')]

        ref,alt ="REF","ALT"
        
        for file in files: 
            chrm,start,stop = None, None, None
            if oldformat:
                chrm,start,p2,p3,stop,suf = file.split(".")
                ref,alt ="ref","var"
            else:
                chrm,start,stop = file.split("_")
                stop = stop[:-4]
            
            cmd = "samtools view -h "+ bamdir+"/"+file + " | head -2 | tail -1 | cut -f3"
            len1 = int(subprocess.check_output(cmd, shell=True).decode('UTF-8').split(":")[1])
            mid1 = int(len1/2)
            
            cmd = "samtools view -h "+ bamdir+"/"+file + " | head -3 | tail -1 | cut -f3"
            len2 = int(subprocess.check_output(cmd, shell=True).decode('UTF-8').split(":")[1])
            mid2 = int(len2/2)
            
            cmd = "samtools depth "+ bamdir+"/"+file +" -r "+ref+":" + str(seq_offs-20) + "-"+ str(seq_offs-10) + " | awk '{s+=$3}END{if(NR>0)print s/NR;else print 0}'"
            dpufr = subprocess.check_output(cmd, shell=True).decode('UTF-8').strip()
            
            cmd = "samtools depth "+ bamdir+"/"+file +" -r "+alt+":" + str(seq_offs-20) + "-"+ str(seq_offs-10) + " | awk '{s+=$3}END{if(NR>0)print s/NR;else print 0}'"
            dpufa = subprocess.check_output(cmd, shell=True).decode('UTF-8').strip()
            
            #cmd = "samtools depth "+ bamdir+"/"+file +" -r REF:" + str(mid1-1) + "-"+ str(mid1+1) + " | awk '{s+=$3}END{if(NR>0)print s/NR;else print 0}'"
            #dpmpr = subprocess.check_output(cmd, shell=True).decode('UTF-8').strip()
            
            #cmd = "samtools depth "+ bamdir+"/"+file +" -r ALT:" + str(mid2-1) + "-"+ str(mid2+1) + " | awk '{s+=$3}END{if(NR>0)print s/NR;else print 0}'"
            #dpmpa = subprocess.check_output(cmd, shell=True).decode('UTF-8').strip()

            cmd = "samtools depth "+ bamdir+"/"+file +" -r "+ref+":" + str(seq_offs+1) + "-"+ str(seq_offs+2) + " | awk '{s+=$3}END{if(NR>0)print s/NR;else print 0}'"
            dpiur = subprocess.check_output(cmd, shell=True).decode('UTF-8').strip()

            cmd = "samtools depth "+ bamdir+"/"+file +" -r "+alt+":" + str(seq_offs+1) + "-"+ str(seq_offs+2) + " | awk '{s+=$3}END{if(NR>0)print s/NR;else print 0}'"
            dpiua = subprocess.check_output(cmd, shell=True).decode('UTF-8').strip()

            cmd = "samtools depth "+ bamdir+"/"+file +" -r "+ref+":" + str(len1-seq_offs-2) + "-"+ str(len1-seq_offs-1) + " | awk '{s+=$3}END{if(NR>0)print s/NR;else print 0}'"
            dpidr = subprocess.check_output(cmd, shell=True).decode('UTF-8').strip()

            cmd = "samtools depth "+ bamdir+"/"+file +" -r "+alt+":" + str(len2-seq_offs-2) + "-"+ str(len2-seq_offs-1) + " | awk '{s+=$3}END{if(NR>0)print s/NR;else print 0}'"
            dpida = subprocess.check_output(cmd, shell=True).decode('UTF-8').strip()

            cmd = "samtools depth "+ bamdir+"/"+file +" -r "+ref+":" + str(len1-seq_offs+10) + "-"+ str(len1-seq_offs+20) + " | awk '{s+=$3}END{if(NR>0)print s/NR;else print 0}'"
            dpdfr = subprocess.check_output(cmd, shell=True).decode('UTF-8').strip()
            
            cmd = "samtools depth "+ bamdir+"/"+file +" -r "+alt+":" + str(len2-seq_offs+10) + "-"+ str(len2-seq_offs+20) + " | awk '{s+=$3}END{if(NR>0)print s/NR;else print 0}'"
            dpdfa = subprocess.check_output(cmd, shell=True).decode('UTF-8').strip()

            print(chrm,start,stop,dpufr,str((float(dpiur)+float(dpidr))/2),dpdfr,dpufa,str((float(dpiua)+float(dpida))/2),dpdfa)

        # loop over allele files
    # condition that bams are scored

######################################################

if __name__ == "__main__":
   main(sys.argv[1:])

######################################################
