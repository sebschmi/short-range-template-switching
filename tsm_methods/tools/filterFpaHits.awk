BEGIN{
    cmd1="zcat "inf"| sort -k2,2n -k17,17n -k26,26n -k27,27n"
    min52=5
    minsim=0.95
    while(cmd1 | getline) {
        # keep longest contig for each locus 
        hit_site=$1"."$2"."$3"."$4"."$5
        if(keep[hit_site]=="") {
            keep[hit_site]=$0;
            keep_length[hit_site]=length($6);
            #print "y "$0
        } else {
            if(length($6)>keep_length[hit_site]) {
                keep[hit_site]=$0;
                keep_length[hit_site]=length($6);
                #print "x "$0
            }
        }
    }
    close(cmd1)

    for(hit_site in keep) {
        $0=keep[hit_site]
        if($0!="") {
            OFS="\t";
            c=$1;s=$2;e=$5;
            if(s>$4){s=$4}
            if(e<$3){e=$3}
            print c,s,e,hit_site > chr"_tmp_fpa.bed"
            print $0 >> outf"_0";
        }
    }
    system("bedtools sort -i "chr"_tmp_fpa.bed > "chr"_tmp_fpa1.bed")
    system("sort --version-sort "outf"_0 > "outf"_1")
    
    # remove those intersecting with a bed file
    cmd2="bedtools intersect -a "chr"_tmp_fpa1.bed -b "bed" -v | sort| uniq > "chr"_tmp_fpa2.bed"
    system(cmd2)
    filein=chr"_tmp_fpa2.bed"
    # ensure that flanking sequences match perfectly
    while(( getline < filein) > 0 ) {
        $0=keep[$4]
        if($0!="") {
            print $0 >> outf"_2";        
            fwd1=substr($6,$12-19,20)
            fwd3=substr($6,$15,20)
            cmd4="samtools faidx -n200 "ref" "$1":"$2-19"-"$2"|grep -v '>'|tr -d '\n'"
            cmd4 | getline frag1
            close(cmd4)
            cmd5="samtools faidx -n200 "ref" "$1":"$5"-"$5+19"|grep -v '>'"
            cmd5 | getline frag3
            close(cmd5)
            if(fwd1==frag1 && fwd3==frag3) {
                hit_site=$1"."$2"."$3"."$4"."$5
                print ">"hit_site"\n"$6 | "cat > "chr"_tmp_fpa.fa"
                print $0 >> outf"_3";
            }
        }
    }
    
    # require that complexity > 0.25
    path="./tools/SeqComplex"
    cmd3="perl -I "path" "path"/profileComplexSeq.pl "chr"_tmp_fpa.fa &> /dev/null"
    while(( cmd3 | getline)>0){ }
    close(cmd3)
    cmd3="sort --version-sort -k1,1 "chr"_tmp_fpa.complex > "chr"_tmp_fpa2.complex"
    while(( cmd3 | getline)>0){ }
    close(cmd3)
    filein=chr"_tmp_fpa2.complex"
    while((getline < filein)>0){
        if($1!~/seq/ && $18>0.25){ 
            $0=keep[$1]
            print $0 >> outf"_4";

            # those have (1)-(4) < 250 bp
            if ( $5-$2<=250 && 
            # minimum distance for (1)-(4)
                 $5-$2>=min52 &&
            # minimum total two differences
                 ($19+$20+$21)-($22+$23+$24)>=2 &&
            # and high similarity for up/downstream
                 $25>=minsim && $26>=minsim && $27>=minsim ) {  
                print $0 >> outf"_5";
                
                # minimum one mismatch
                if($21>=1) {
                    print $0 >> outf"_6";
                    
                    if(keep2[$6]=="") {
                        keep2[$6]=$0;
                        keep2_ident[$6]=$26;
                    } else {
                        if($26>keep2_ident[$6]) {
                            keep2[$6]=$0;
                            keep2_ident[$6]=$26;
                        }
                    }
                }
            }
        }
    }
    for(hit_frag in keep2) {
        $0=keep2[hit_frag]
        #hit_site=$1"."$2"."$3"."$4"."$5
        keep3[$2]=$0;
    }
    n = asort(keep3,keep4,"@ind_num_asc")
    for (i = 1; i <= n; i++) {
        $0=keep4[i]
        print $0 >> outf"_7";
    }
    
    system("rm -f "chr"_tmp_fpa.bed "chr"_tmp_fpa1.bed "chr"_tmp_fpa2.bed")
    system("rm -f "chr"_tmp_fpa.fa "chr"_tmp_fpa.complex "chr"_tmp_fpa2.complex")
}
