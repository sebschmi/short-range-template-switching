{
    len1=500
    cmd1="samtools faidx -n600 "ref" "$1":"$2-len1"-"$2"|grep -v '>'|tr -d '\n'"
    cmd1 | getline frag1
    close(cmd1)

    cmd3="samtools faidx -n600 "ref" "$1":"$5"-"$5+len1"|grep -v '>'"
    cmd3 | getline frag3
    close(cmd3)

    if($31~/S/) {
        #cmd2="samtools faidx -n300 "ref" "$1":"$2+1"-"$5-1"|grep -v '>'|tr -d '\n'"
        cmd2="samtools faidx -n600 "ref" "$1":"$4"-"$3"|grep -v '>'|tr -d '\n'";  
        frag2=""
        if($2+1<=$5-1)
            cmd2 | getline frag2
        close(cmd2)

        fwd2=substr($6,$12+1,$15-$12-1)
    } else {
        cmd2="samtools faidx -n600 "ref" "$1":"$2+1"-"$5-1"|grep -v '>'|tr -d '\n'"
        frag2=""
        if($2+1<=$5-1)
            cmd2 | getline frag2
        close(cmd2)

        fwd2=substr($6,$13,$14-$13+1)
    }
    fasfile=fas"/"$1"."$2"."$3"."$4"."$5".fa"
    OFS=""
    printf(">ref\n"frag1""frag2""frag3"\n") > fasfile
    printf(">var\n"frag1""fwd2""frag3"\n") > fasfile
    close(fasfile)
}
