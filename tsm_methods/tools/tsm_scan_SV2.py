import numpy as np
import pandas as pd
from itertools import compress
from pyfaidx import Faidx

import sys, getopt

sys.path.insert(0,'./tools/fpa')
import fpa_ext2

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

    chrnum = ''
    fafile = ''
    gtfile = ''
    indiv = ''
    untrimmed = False
    trimlen = 50
    alltypes = True
    nbases = 4
    
    try:
        opts, args = getopt.getopt(argv,"hc:f:g:i:ut:n:S",["chrom=","fafile=","gtfile=","indiv=","untrimmed","trimlen=","nbases=","svonly"])
        
    except getopt.GetoptError:
        print ('tsm_scan.py -c chr -f <fastafile> -g <genotypefile> -i <individual> -u -t <trimlen> -n <nbases> -S')
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-h':
            print ('tsm_scan.py -c chr -f <fastafile> -g <genotypefile> -i <individual> -u -t <trimlen> -n <nbases> -S')
            sys.exit()
        elif opt in ("-c", "--chrom"):
            chrnum = arg
        elif opt in ("-f", "--fafile"):
            fafile = arg
        elif opt in ("-g", "--gtfile"):
            gtfile = arg
        elif opt in ("-i", "--indiv"):
            indiv = int(arg)
        elif opt in ("-u", "--untrimmed"):
            untrimmed = True
        elif opt in ("-t", "--trimlen"):
            trimlen = int(arg)
        elif opt in ("-n", "--nbases"):
            nbases = int(arg)
        elif opt in ("-S", "--svonly"):
            alltypes = False
            
    if(fafile == '' or gtfile == ''):
        print ('tsm_scan.py -c chr -f <fastafile> -g <genotypefile> -i <individual> -u -t <trimlen> -n <nbases> -S')
        sys.exit()
        
    fasta = Faidx(fafile)
    data = pd.read_csv(gtfile, sep="\t|\|", dtype=object, header=None, engine='python')

######################################################

    chrm = np.array(data.iloc[:,0])
    strt = np.array(data.iloc[:,1].astype(int))
    refc = np.array(data.iloc[:,2])
    altc = np.array(data.iloc[:,3])
    stop = strt + [len(s)-1 for s in refc]
    svtp = np.array(data.iloc[:,data.shape[1]-1])
    
    numind = data.shape[1]-5
    geno = np.empty((data.shape[0],numind), dtype=bool)
    for i in range(4,data.shape[1]-1):
        geno[:,(i-4)] = np.array(data.iloc[:,i] == "1")

######################################################

    max_dist = 50
    seq_offs = 150

######################################################

    f = fpa_ext2.FPA2()
    f.set_int("min_length",8)
    f.set_int("scan_flank",100)

    iset = range(0,numind)
    
    if(indiv != ''):
        iset = [indiv]

    for ind in iset:
    
        #print("haplotype ",ind)

        if sum(geno[:,ind])>1:
            chm = chrm[geno[:,ind]]
            sta = strt[geno[:,ind]]
            sto = stop[geno[:,ind]]
            ref = refc[geno[:,ind]]
            alt = altc[geno[:,ind]]
            svt = svtp[geno[:,ind]]
            
            cset = ["chr" + str(i) for i in range(1,22)]+['chrX','chrY']
            
            if chrnum != '':
                cset = [chrnum]
            
            for ic in cset:
                
                keep = np.logical_and(np.equal(chm,ic), sta-np.array(shift(sto,1,-1000)) < max_dist)

                #print(ind,ic,len(keep))
                #print(sta[keep])
                
                for p in [i for i, x in enumerate(keep) if x]:
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
                        has_sv = False

                        for c in clst:
                            seq_offs_extra +=  len(ref[c])
                            if(svt[c]=="SV"):
                                has_sv = True


                        if(has_sv or alltypes):

                            seq = str(fasta.fetch(ic,sta[clst[0]]-seq_offs,sta[clst[-1]]+seq_offs+seq_offs_extra))

                            seq_st = sta[clst[0]]-seq_offs
                            seq_po = 0
                            ref_seq = str()
                            alt_seq = str()

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

                            ref_seq += seq[seq_po:]
                            alt_seq += seq[seq_po:]

                            locus = seq[sta[clst[0]]-seq_st:sta[clst[-1]]-seq_st+1]
                            mask = "masked"
                            if locus.upper() == locus:
                                mask = "unmasked"
                                

                            hits1 = f.scan_two(ref_seq,alt_seq,False)
                            hits2 = f.scan_two(alt_seq,ref_seq,False)

                            score1 = score2 = -1

                            for h in hits1:
                                a = h.split(",")
                                for i in range(1,len(a)):
                                    a[i] = float(a[i])

                                if a[10] > 0.9 and a[11] > 0.95 and a[12] > 0.9 and a[13] > a[15] and a[17]+a[18]+a[19] > 1 and a[19] > 0 and a[20]>=nbases:
                                    score1 = a[13]

                            for h in hits2:
                                a = h.split(",")
                                for i in range(1,len(a)):
                                    a[i] = float(a[i])

                                if a[10] > 0.9 and a[11] > 0.95 and a[12] > 0.9 and a[13] > a[15] and a[17]+a[18]+a[19] > 1 and a[19] > 0 and a[20]>=nbases:
                                    score2 = a[13]


                            if score1>score2 or score1==score2:

                                scannum = 1
                                if(score1==score2):
                                    scannum = 0
                                
                                for h in hits1:
                                    a = h.split(",")
                                    for i in range(1,len(a)):
                                        a[i] = float(a[i])

                                    if a[10] > 0.9 and a[11] > 0.95 and a[12] > 0.9 and a[13] > a[15] and a[17]+a[18]+a[19] > 1 and a[19] > 0 and a[20]>=nbases:
                                        print("sample",ind,chm[p],sta[clst[0]],sta[clst[-1]],"[",len(clst),ld,"] ",scannum, mask, flush= True)
                                        print(h,ref_seq,alt_seq,"",sep="\n", flush= True)

                                        f.print_hit(ref_seq,alt_seq,h,True)

                            if score2>score1:

                                scannum = 2

                                for h in hits2:
                                    a = h.split(",")
                                    for i in range(1,len(a)):
                                        a[i] = float(a[i])

                                    if a[10] > 0.9 and a[11] > 0.95 and a[12] > 0.9 and a[13] > a[15] and a[17]+a[18]+a[19] > 1 and a[19] > 0 and a[20]>=nbases:
                                        print("sample",ind,chm[p],sta[clst[0]],sta[clst[-1]],"[",len(clst),ld,"] ",scannum, mask, flush= True)
                                        print(h,ref_seq,alt_seq,"",sep="\n", flush= True)

                                        f.print_hit(alt_seq,ref_seq,h,True)


######################################################

if __name__ == "__main__":
   main(sys.argv[1:])

######################################################

# ~ chrom #0
# ~ clus_start_chrom
# ~ clus_start_align
# ~ clust_start1
# ~ clust_end1
# ~ sp1_qry
# ~ sp1_ref
# ~ sp2_ref
# ~ sp3_ref
# ~ sp4_ref
# ~ iden_up #10
# ~ ident_rep
# ~ ident_down
# ~ ident_inv
# ~ ident_fwd
# ~ ident_org
# ~ masked
# ~ delta sum_ins
# ~ delta sum_del
# ~ delta sum_mis
# ~ sum_nuc #20
# ~ CpG
# ~ clus_ins
# ~ clus_del
# ~ clus_mis
