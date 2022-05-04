/***************************************************************************
 *   Copyright (C) 2010-2018 by Ari Löytynoja                              *
 *   ari.loytynoja@gmail.com                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "fpa.h"

/********************    alignment stuff    ******************************/

/*
void FPA::build_indeces(std::string *s1,std::string *s2)
{

    if(s1->length() != s2->length())
    {
        std::cout<<"expecting aligned sequences. exiting.\n\n";
        exit(0);
    }

    slg = s1->length();

    index1.reserve(slg);
    index2.reserve(slg);
    rindex1.reserve(slg);
    rindex2.reserve(slg);

    fseq1.reserve(slg);
    fseq2.reserve(slg);

    int p1=0; int p2=0;

    for(int i=0;i<slg;i++)
    {
        index1.push_back(p1);
        index2.push_back(p2);

        if(s1->at(i) != '-')
        {
            rindex1.push_back(i);
            p1++;
        }
        if(s2->at(i) != '-')
        {
            rindex2.push_back(i);
            p2++;
        }
    }

    fsl1 = p1;
    fsl2 = p2;

    rindex1.resize(p1);
    rindex2.resize(p2);


    std::string alpha = "ACGTN";

    int ci;

    std::string::iterator sit = s1->begin();
    for(;sit != s1->end();sit++)
    {
        ci = alpha.find(toupper(*sit));
        if (ci>=0 && ci<5)
            fseq1.push_back(ci);
        else
            fseq1.push_back(-1);
    }

    sit = s2->begin();
    for(;sit != s2->end();sit++)
    {
        ci = alpha.find(toupper(*sit));
        if (ci>=0 && ci<5)
            fseq2.push_back(ci);
        else
            fseq2.push_back(-1);
    }

}
*/

void FPA::build_sequences(std::string *s1,std::string *s2, bool reverse)
{

    seq1.reserve(s1->length());
    rev1.reserve(s1->length());
    seq2.reserve(s2->length());
    mask1.reserve(s1->length());

    std::string alpha = "ACGTN";

    int ci;
    int p1=0; int p2=0;

    std::string::iterator sit = s1->begin();
    for(;sit != s1->end() ;sit++)
    {
        ci = alpha.find(toupper(*sit));
        if (ci>=0 && ci<5)
        {
            seq1.push_back(ci);

            int m = 0;
            if(islower(*sit))
                m += 1;

            mask1.push_back(m);

            p1++;
        }
    }


    sit = s1->begin();
    if(reverse)
    {
        for(;sit != s1->end();sit++)
        {
            ci = alpha.find(toupper(*sit));
            if (ci>=0 && ci<5)
                rev1.push_back(ci);

        }
    }
    else
    {
        for(;sit != s1->end();sit++)
        {
            ci = alpha.find(toupper(*sit));
            if (ci>=0 && ci<4)
                rev1.push_back(3-ci);
            else if (ci==4)
                rev1.push_back(ci);

        }
    }

    sit = s2->begin();
    for(;sit != s2->end();sit++)
    {
        ci = alpha.find(toupper(*sit));
        if (ci>=0 && ci<5)
        {
            seq2.push_back(ci);
            p2++;
        }
    }

    seq1.resize(p1);
    rev1.resize(p1);
    seq2.resize(p2);
    mask1.resize(p1);

    start1 = 0;
    end1 = seq1.size();
    start2 = 0;
    end2 = seq2.size();

    sl1 = end1-start1;
    sl2 = end2-start2;

    true_start2 = start2;
    true_end2 = end2;

}



void FPA::fwd_align_sequences(std::vector<seqCoordinate> *path)
{
    int gapscore = -2;

    int ei = seq1.size();
    int ej = seq2.size();

    Array2D<int> mat1(ei+1,ej+1);
    Array2D<int> ptr1(ei+1,ej+1);

    int large_neg = -10000;

    for(int j=0;j<=ej;j++)
    {
        // mat1
        for(int i=0;i<=ei;i++)
        {
            int ptr = none;
            int score = large_neg;

            if(i==0 && j==0)
//            if(i==0 || j==0)
                score = 0;

            else
            {
                if(i>0 && j>0)
                {
                    score = mat1(i-1,j-1) + substitution_score(i,j);
                    ptr = match;
                }

                if(i>0)
                {
                    int this_score = mat1(i-1,j) + gapscore;
                    if(this_score > score)
                    {
                        score = this_score;
                        ptr = xgap;
                    }
                }

                if(j>0)
                {
                    int this_score = mat1(i,j-1) + gapscore;
                    if(this_score > score)
                    {
                        score = this_score;
                        ptr = ygap;
                    }
                }

                if(i>0 && j==ej)
                {
                    int this_score = mat1(i-1,j);
                    if(this_score > score)
                    {
                        score = this_score;
                        ptr = xgap;
                    }
                }

                if(j>0 && i==ei)
                {
                    int this_score = mat1(i,j-1);
                    if(this_score > score)
                    {
                        score = this_score;
                        ptr = ygap;
                    }
                }
            }

            mat1(i,j) = score;
            ptr1(i,j) = ptr;

//            std::cout<<ptr<<" ";
        }
//        std::cout<<"\n";
    }

    int i=ei;
    int j=ej;

    bool has_match=true;

//    for(;i>=0 || j>=0;)
    for(;i>=0 && j>=0;)
    {
        seqCoordinate c;
        c.matrix = 1;
        c.pos_x = i;
        c.pos_y = j;

        if(ptr1(i,j)==match)
        {
            path->push_back(c);
            i--;j--;
            has_match=true;
        }

        else if(ptr1(i,j)==xgap)
        {
            c.pos_y = -1;
            if(has_match)
                path->push_back(c);
            i--;
        }

        else if(ptr1(i,j)==ygap)
        {
            c.pos_x = -1;
            if(has_match)
                path->push_back(c);
            j--;
        }

        else if(ptr1(i,j)==none)
        {
            break;
        }

    }
}

void FPA::align_sequences(std::vector<seqCoordinate> *path,std::vector<switchPoint> *points,bool local)
{

    int maxtermgap = 0;
    int gapscore = -2;

    Array2D<int> mat1(sl1+1,sl2+1);
    Array2D<int> mat2(sl1+1,sl2+1);
    Array2D<int> mat3(sl1+1,sl2+1);

    Array2D<int> ptr1(sl1+1,sl2+1);
    Array2D<int> ptr2(sl1+1,sl2+1);
    Array2D<int> ptr3(sl1+1,sl2+1);

    Array2D<int> sco2(sl1+1,sl2+1);
    Array2D<int> len2(sl1+1,sl2+1);
    Array2D<int> dist2(sl1+1,sl2+1);

    Array2D<int> sco3(sl1+1,sl2+1);
    Array2D<int> len3(sl1+1,sl2+1);

    int large_neg = -10000;

    if(debug_matrix)
    {
        std::cout<<start1<<" "<<end1<<"; "<<start2<<" "<<end2<<"; "<<true_start2<<" "<<true_end2<<" | "<<sl1<<" "<<sl2<<" | "<<seq1.size()<<" "<<seq2.size()<<"\n";
        for(int i=0;i<sl1;i++)
            std::cout<<" "<<std::string("-ACGTN").at(seq1.at(i)+1);
        std::cout<<std::endl;
    }


    int maxs = 0;
    int maxi = 0;
    int maxj = 0;

    for(int j=0;j<=sl2;j++)
    {
        if(debug_matrix && j<sl2) std::cout<<start2+j<<" "<<std::string("-ACGTN").at(seq2.at(j)+1);
        // mat1
        for(int i=0;i<=sl1;i++)
        {
            int ptr = none;
            int score = large_neg;

            if(i==0 && j==0)
                score = 0;
            else if( (i==0 && j<=maxtermgap) || (j==0 && i<=maxtermgap) )
                score = 0;
            else if(not local && force_overlap && j+start2<true_start2)
                ;
            else if(not local && force_overlap )
                ;
            else if(not local && force_overlap )
                ;
            else if(not local &&  i==0 && j+start2==true_start2 )
                score = 0;

            else
            {
                if(i>0 && j>0)
                {
                    score = mat1(i-1,j-1) + substitution_score(start1+i,start2+j);
                    ptr = match;
                }

                if(i>0)
                {
                    int this_score = mat1(i-1,j) + gapscore;
                    if(this_score > score)
                    {
                        score = this_score;
                        ptr = xgap;
                    }
                }

                if(j>0)
                {
                    int this_score = mat1(i,j-1) + gapscore;
                    if(this_score > score)
                    {
                        score = this_score;
                        ptr = ygap;
                    }
                }
            }

            if(score>0 || not local)
            {
                mat1(i,j) = score;
                ptr1(i,j) = ptr;
            }
            else if(local)
            {
                mat1(i,j) = 0;
                ptr1(i,j) = none;
            }

            if(score>maxs)
            {
                maxs=score;
                maxi=i;
                maxj=j;
            }

            if(debug_matrix) std::cout<<" "<<score;
        }
        if(debug_matrix) std::cout<<"\n";
    }
    if(debug_matrix) std::cout<<"\n";
    if(debug_matrix) std::cout<<maxs<<" "<<maxi<<" "<<maxj<<"\n";


    for(int i=0;i<=sl1;i++  )
        {

        int pscore = large_neg;
        int pk = -1;

        // from mat1
        for(int k=0;k<=sl2;k++)
        {
            if(mat1(i-1,k) > pscore)
            {
                 pscore = mat1(i-1,k);
                 pk = k;
            }
        }

        // mat2
        for(int j=sl2;j>=0;j--)
        {
            if(debug_matrix && j<sl2) std::cout<<j<<" "<<std::string("-ACGTN").at(seq2.at(j)+1);

            int score = large_neg;
            int subst = 0;
            int ptr = none;
            int len = 0;
            int sco = 0;
            int dist = sl2;

            if(i>0 && j>0)
            {
                ptr = match;
                subst = rev_substitution_score(start1+i,start2+j);

                if(j<sl2)
                {
                    score = subst + mat2(i-1,j+1);
                    len = len2(i-1,j+1)+1;
                    sco = sco2(i-1,j+1)+subst;
                    dist = dist2(i-1,j+1);
                }

                // from mat1
                if(subst + pscore > score)
                {
                    score = subst + pscore;
                    ptr = pk;
                    len = 1;
                    sco = subst;
                    dist = abs(j-pk);
                }
            }

            mat2(i,j) = score;
            ptr2(i,j) = ptr;

            len2(i,j) = len;
            sco2(i,j) = sco;

            dist2(i,j) = dist;

            if(debug_matrix) std::cout<<" "<<score<<"; ";
        }
        if(debug_matrix) std::cout<<"\n";
    }
    if(debug_matrix) std::cout<<"\n";

    maxs = 0;
    maxi = 0;
    maxj = 0;

    for(int i=0;i<=sl1;i++)
    {
        for(int j=0;j<=sl2;j++)
        {
            if(debug_matrix && j<sl2) std::cout<<j<<" "<<std::string("-ACGTN").at(seq2.at(j)+1);
            // mat3
            int score = large_neg;
            int ptr = none;
            int len = 0;
            int sco = 0;

            if(i>0 && j>0)
            {
                if(not local && force_overlap)
                        ;
                else if(not local && force_overlap)
                    ;
                else if(not local && force_overlap && j+start2>true_end2)
                    ;
                else
                {
                    ptr = match;
                    int subst = substitution_score(start1+i,start2+j);
                    score = subst + mat3(i-1,j-1);

                    len = len3(i-1,j-1);
                    sco = sco3(i-1,j-1);

                    int this_score = mat3(i-1,j) + gapscore;
                    if(this_score > score)
                    {
                        score = this_score;
                        ptr = xgap;

                        len = len3(i-1,j);
                        sco = sco3(i-1,j);
                    }

                    this_score = mat3(i,j-1) + gapscore;
                    if(this_score > score)
                    {
                        score = this_score;
                        ptr = ygap;

                        len = len3(i,j-1);
                        sco = sco3(i,j-1);
                    }
                    if(i==sl1 && j>=true_end2)
                    {
                        this_score = mat3(i,j-1);
                        if(this_score > score)
                        {
                            score = this_score;
                            ptr = ygap;

                            len = len3(i,j-1);
                            sco = sco3(i,j-1);
                        }
                    }

                    // from mat2
                    int dist = sl2;
                    for(int k=0;k<=sl2;k++)
                    {
                        if(subst + mat2(i-1,k) > score ||
                           ( subst + mat2(i-1,k) == score && dist2(i-1,k)<dist) )
                        {
                            score = subst + mat2(i-1,k);
                            ptr = k;

                            len = len2(i-1,k);
                            sco = sco2(i-1,k);
                            dist = dist2(i-1,k);

                            if(score>maxs)
                            {
                                maxs=score;
                                maxi=i;
                                maxj=k;
                            }

                        }
                    }
                }
            }
            if(score>0 || not local)
            {
                mat3(i,j) = score;
                ptr3(i,j) = ptr;

                len3(i,j) = len;
                sco3(i,j) = sco;
            }
            else if(local)
            {
                mat3(i,j) = 0;
                ptr3(i,j) = none;

                len3(i,j) = 0;
                sco3(i,j) = 0;
            }
            if(debug_matrix) std::cout<<" "<<score<<"; ";
        }
        if(debug_matrix) std::cout<<"\n";
    }
    if(debug_matrix) std::cout<<"\n";

    if(debug_matrix) std::cout<<maxs<<" "<<maxi<<" "<<maxj<<"\n";

    int max_i=-1;
    int max_j=-1;
    int max_end=large_neg;
    int max_len = -1;
    int max_sco = -1;

    for(int i=sl1;i>0;i--)
    {
        for(int j=sl2;j>0;j--)
        {
            if(maximise_score)
            {
                if(sco3(i,j)>max_sco)
                {
                    max_i = i;
                    max_j = j;
                    max_end = mat3(i,j);
                    max_len = len3(i,j);
                    max_sco = sco3(i,j);
                }
            }
            else if (maximise_length)
            {
                if(len3(i,j)>max_len)
                {
                    max_i = i;
                    max_j = j;
                    max_end = mat3(i,j);
                    max_len = len3(i,j);
                    max_sco = sco3(i,j);
                }
            }
            else
            {
                if(mat3(i,j)>max_end)
                {
                    max_i = i;
                    max_j = j;
                    max_end = mat3(i,j);
                    max_len = len3(i,j);
                    max_sco = sco3(i,j);
                }
            }
        }
    }


    int i=max_i;
    int j=max_j;
    bool ptr_found = false;

    for(;i>=0 || j>=0;)
    {
        seqCoordinate c;
        c.matrix = 3;
        c.pos_x = start1+i;
        c.pos_y = start2+j;

        if(ptr3(i,j)==match)
        {
            path->push_back(c);

            i--;j--;
        }

        else if(ptr3(i,j)==xgap)
        {
            c.pos_y = -1;
            path->push_back(c);

            i--;
        }

        else if(ptr3(i,j)==ygap)
        {
            c.pos_x = -1;
            path->push_back(c);

            j--;
        }

        else if(ptr3(i,j)==none)
        {
            std::cout<<"ptr3 = none. weird. exiting.\n";
            exit(0);
        }

        else
        {
            path->push_back(c);

            /*here -> correct*/
            points->at(3).i = start1+i;
            points->at(3).j = start2+j;

            j = ptr3(i,j);
            i--;

            points->at(2).i = start1+i;
            points->at(2).j = start2+j;

            ptr_found = true;
            break;
        }
    }

    if(!ptr_found)
    {
        std::cout<<"backtracking 1st path failed. exiting.\n";
        exit(0);
    }

    ptr_found = false;

    for(;i>=0 || j<=sl2;)
    {
        seqCoordinate c;
        c.matrix = 2;
        c.pos_x = start1+i;
        c.pos_y = start2+j;
        path->push_back(c);
        if(ptr2(i,j)==match)
        {
            i--;j++;
        }

        else if(ptr2(i,j)==none)
        {
            std::cout<<"ptr2 = none. weird. exiting.\n";
            exit(0);
        }

        else
        {
            points->at(1).i = start1+i;
            points->at(1).j = start2+j;

            j = ptr2(i,j);
            i--;

            points->at(0).i = start1+i;
            points->at(0).j = start2+j;

            ptr_found = true;
            break;
        }
    }

    if(!ptr_found)
    {
        std::cout<<"backtracking 2nd path failed. exiting.\n";
        exit(0);
    }

    for(;i>=0 || j>=0;)
    {
        seqCoordinate c;
        c.matrix = 1;
        c.pos_x = start1+i;
        c.pos_y = start2+j;

        if(ptr1(i,j)==match)
        {
            path->push_back(c);
            i--;j--;
        }

        else if(ptr1(i,j)==xgap)
        {
            c.pos_y = -1;
            path->push_back(c);
            i--;
        }

        else if(ptr1(i,j)==ygap)
        {
            c.pos_x = -1;
            path->push_back(c);
            j--;
        }

        else if(ptr1(i,j)==none)
        {
            break;
        }
    }
}


/********************    alignment stuff    ******************************/


/********************    alignment output   ******************************/

void FPA::print_switch_process(std::vector<seqCoordinate> *path,std::vector<switchPoint> *points)
{
    std::string alpha = "ACGTN";

    std::string out1("");
    std::string out2("");
    std::string out2_gaps;
    std::string out3(" ");


    int point2 = points->at(1).j;
    int point3 = points->at(2).j;

    int p=path->size()-1;
    int site1 = path->at(p).pos_x;
    int site2 = path->at(p).pos_y;
    int site_mat = path->at(p).matrix;

    int first_site2 = site2;
    if(first_site2<0)
    {
        for(int i=p-1;i>0 && first_site2<0;i--)
            first_site2 = path->at(i).pos_y;
    }


    std::string qry(" ");
    std::string ref(" ");
    std::string rref(" ");

    int ref_pos = site2+1;
    int ref_end = 0;
    if(point3<site2)
    {
        qry += " ";
        ref += " ";
        rref += " ";
        out1 += " ";
        out3 += " ";

        for(int i=point3;i<site2;i++)
        {
            out1 += " ";
            out3 += " ";
            ref += alpha.at(seq2.at(i-1));
            if(seq2.at(i-1)<4)
                rref += alpha.at(3-seq2.at(i-1));
            else
                rref += "N";
            ref_end = i-1;
        }
    }

    out1 += "\bL ";

    while(site_mat==1)
    {
        if(site1>=0 && site1<=(int)seq1.size())
        {
            out1+=alpha.at(seq1.at(site1-1));
            qry+=alpha.at(seq1.at(site1-1));
        }
        else
        {
            out1+="-";
            qry+="-";
        }
        if(site2>=0 && site2<=(int)seq2.size())
        {
            ref+=alpha.at(seq2.at(site2-1));
            if(seq2.at(site2-1)<4)
                rref+=alpha.at(3-seq2.at(site2-1));
            else
                rref += "N";
            ref_pos = site2+1;
        }
        else
        {
            ref+="-";
            rref+="-";
            out2_gaps+=" ";
            out3+=" ";
        }
        ref_end = site2-1;

        p--;
        site1 = path->at(p).pos_x;
        site2 = path->at(p).pos_y;
        site_mat = path->at(p).matrix;
    }
    out1+=std::string(" 1");

    int first_site3 = site2;
    if(first_site3<0)
    {
        for(int i=p-1;i>0 && first_site3<0;i--)
            first_site3 = path->at(i).pos_y;
    }

    qry+="1 3";

    int last_y=0;

    while(site_mat==2)
    {
        if(site1>=0 && site1<=(int)seq1.size())
        {
            out2=alpha.at(seq1.at(site1-1))+out2;
            qry+=alpha.at(seq1.at(site1-1));
        }
        else
        {
            std::cout<<"error!\n";
        }

        last_y = site2;

        p--;
        site1 = path->at(p).pos_x;
        site2 = path->at(p).pos_y;
        site_mat = path->at(p).matrix;
    }


    if(last_y>1)
        out2_gaps+=" ";

    if(ref_pos<=0)
        ref_pos=1;
    for(int i=ref_pos;i<site2 && i<(int)seq2.size();i++)
    {
        ref+=alpha.at(seq2.at(i-1));
        if(seq2.at(i-1)<4)
            rref += alpha.at(3-seq2.at(i-1));
        else
            rref += "N";
        ref_end = i-1;
    }

    out2=out2_gaps+std::string("\b3 ")+out2+std::string(" 2");
    for(int i=first_site2;i<last_y-1;i++)
    {
        out2=std::string(" ")+out2;
    }


    for(int i=first_site2;i<site2-1;i++)
    {
        out3+=std::string(" ");
        if(seq2.at(i)<0)
            out3+=std::string(" ");
    }

    out3+=std::string("\b4 ");

    qry+="2 4";

    while(site_mat==3 )
    {
        if(site1>=0 && site1<=(int)seq1.size())
        {
            out3+=alpha.at(seq1.at(site1-1));
            qry+=alpha.at(seq1.at(site1-1));
        }
        else
        {
            out3+="-";
            qry+="-";
        }
        if(site2>=0 && site2<=(int)seq2.size() && site2-1>ref_end)
        {
            ref+=alpha.at(seq2.at(site2-1));
            if(seq2.at(site2-1)<4)
                rref+=alpha.at(3-seq2.at(site2-1));
            else
                rref += "N";
        }
        else if( site2-1>ref_end )
        {
            ref+="-";
            rref+="-";
        }

        p--;
        if(p>=0)
        {
            site1 = path->at(p).pos_x;
            site2 = path->at(p).pos_y;
            site_mat = path->at(p).matrix;
        }
        else
            break;
    }
    site2++;
    if(point2>site2)
    {
        for(;site2<=point2;site2++)
        {
            ref+=alpha.at(seq2.at(site2-1));
            if(seq2.at(site2-1)<4)
                rref+=alpha.at(3-seq2.at(site2-1));
            else
                rref += "N";
        }
    }
    out3+=std::string(" R");

    std::string qry0;
    for(int i=start1;i<end1;i++)
        qry0+=alpha.at(seq1.at(i));

    std::string ref0;
    for(int i=start2;i<end2;i++)
        ref0+=alpha.at(seq2.at(i));

    std::string ref00;
    for(int i=true_start2;i<true_end2;i++)
        ref00+=alpha.at(seq2.at(i));

    std::cout<<"Switch process:"<<"\nF1:  "<<out1<<"\nF3:  "<<out3<<"\nRF:  "<<ref<<"\nRR:  "<<rref<<"\nF2:  "<<out2<<"\n\n";

}



void FPA::compute_fpa_score(std::vector<seqCoordinate> *path,std::vector<switchPoint> *points,std::vector<seqCoordinate> *fwd_path,int flanking)
{

    int start_i = std::max(points->at(0).i-flanking,0);
    int stop_i = points->at(0).i;

    int p=path->size()-1;

    int site1 = path->at(p).pos_x;
    int site2 = path->at(p).pos_y;

    while(site1<=start_i)
    {
        p--;
        site1 = path->at(p).pos_x;
        site2 = path->at(p).pos_y;
    }
//    int start_site1 = site1;
//    int start_site2 = site2;

    seq1_frag1.clear();
    seq2_frag1.clear();
    low_frag1.clear();

    while(site1<=stop_i)
    {
        if(site1>0 && site1<=(int)seq1.size() &&
           site2>0 && site2<=(int)seq2.size() &&
           seq1.at(site1-1) != seq2.at(site2-1) )
            low_frag1.push_back(true);
        else if(site2<0)
            low_frag1.push_back(true);
        else
            low_frag1.push_back(false);

        if(site1>0 && site1<=(int)seq1.size())
            seq1_frag1.push_back(seq1.at(site1-1));
        else
            seq1_frag1.push_back(-1);

        if(site2>0 && site2<=(int)seq2.size())
            seq2_frag1.push_back(seq2.at(site2-1));
        else
            seq2_frag1.push_back(-1);

        p--;
        site1 = path->at(p).pos_x;
        site2 = path->at(p).pos_y;
    }

    seq1_frag2.clear();
    seq2_frag2.clear();
    low_frag2.clear();

    stop_i = points->at(2).i;

    int sA=0; int sC=0; int sG=0; int sT=0;
    while(site1<=stop_i)
    {
        if(site1>0 && site1<=(int)seq1.size() &&
           site2>0 && site2<=(int)seq2.size() &&
          ( (reverse && seq1.at(site1-1) != seq2.at(site2-1) ) ||
            ( not reverse && seq1.at(site1-1) != 3-seq2.at(site2-1) ) ) )
            low_frag2.push_back(true);
        else if(site2<0)
            low_frag2.push_back(true);
        else
            low_frag2.push_back(false);

        if(site1>0 && site1<=(int)seq1.size())
        {
            seq1_frag2.push_back(seq1.at(site1-1));
            int c = seq1.at(site1-1);
            if(c==0) sA=1;
            if(c==1) sC=1;
            if(c==2) sG=1;
            if(c==3) sT=1;
        }
        else
            seq1_frag2.push_back(-1);

        if(site2>0 && site2<=(int)seq2.size())
            seq2_frag2.push_back(3-seq2.at(site2-1));
        else
            seq2_frag2.push_back(-1);

        p--;
        site1 = path->at(p).pos_x;
        site2 = path->at(p).pos_y;
    }
    sumNuc = sA+sC+sG+sT;

    seq1_frag3.clear();
    seq2_frag3.clear();
    low_frag3.clear();

    stop_i = std::min(points->at(3).i+flanking,(int)seq1.size());

    while(site1<stop_i && p>=0)
    {
        if(site1>0 && site1<=(int)seq1.size() &&
           site2>0 && site2<=(int)seq2.size() &&
           seq1.at(site1-1) != seq2.at(site2-1) )
            low_frag3.push_back(true);
        else if(site2<0)
            low_frag3.push_back(true);
        else
            low_frag3.push_back(false);

        if(site1>0 && site1<=(int)seq1.size())
            seq1_frag3.push_back(seq1.at(site1-1));
        else
            seq1_frag3.push_back(-1);

        if(site2>0 && site2<=(int)seq2.size())
            seq2_frag3.push_back(seq2.at(site2-1));
        else
            seq2_frag3.push_back(-1);

        p--;
        if(p>=0)
        {
            site1 = path->at(p).pos_x;
            site2 = path->at(p).pos_y;
        }
    }

//    int stop_site1 =site1;
//    int stop_site2 =site2;

    seq1_fwd.clear();
    seq2_fwd.clear();
    low_fwd.clear();

    p=fwd_path->size()-1;

    site1 = fwd_path->at(p).pos_x;
    site2 = fwd_path->at(p).pos_y;

//    std::cout<<site1<<" "<<site2<<" "<<start_i<<" "<<stop_i<<" : "<<fwd_path->size()<<std::endl;
    /*
    std::string alpha = "-ACGTN";

    for(int i=fwd_path->size()-1;i>=0;--i)
    {
        site1 = fwd_path->at(i).pos_x;
        site2 = fwd_path->at(i).pos_y;

        if(site1-1>=0)
            std::cout<<alpha.at(seq1.at(site1-1)+1)<<" ";
        else
            std::cout<<"- ";
        if(site2-1>=0)
            std::cout<<alpha.at(seq2.at(site2-1)+1)<<"\n";
        else
            std::cout<<"-\n";

    }
    */

    while(site1<=start_i)
    {
        p--;
        site1 = fwd_path->at(p).pos_x;
        site2 = fwd_path->at(p).pos_y;
    }
    stop_i = std::min(points->at(3).i+flanking,(int)seq1.size());

//    std::cout<<site1<<" "<<site2<<" "<<stop_i<<std::endl;


//    std::cout<<site1<<" "<<site2<<" "<<stop_i<<" : "<<start_site1<<" "<<start_site2<<" "<<stop_site1<<" "<<stop_site2<<std::endl;
    while(site1<stop_i && p>=0)
    {
        if(site1>0 && site1<=(int)seq1.size() &&
           site2>0 && site2<=(int)seq2.size() &&
           seq1.at(site1-1) != seq2.at(site2-1) )
            low_fwd.push_back(true);
        else if(site2<0)
            low_fwd.push_back(true);
        else
            low_fwd.push_back(false);

        if(site1>0 && site1<=(int)seq1.size())
            seq1_fwd.push_back(seq1.at(site1-1));
        else
            seq1_fwd.push_back(-1);

        if(site2>0 && site2<=(int)seq2.size())
            seq2_fwd.push_back(seq2.at(site2-1));
        else
            seq2_fwd.push_back(-1);

        /*
        if(site1-1>=0)
            std::cout<<seq1.at(site1-1)<<" ";
        else
            std::cout<<"- ";
        if(site2-1>=0)
            std::cout<<seq2.at(site2-1)<<"\n";
        else
            std::cout<<"-\n";
        //*/

        p--;
        if(p>=0)
        {
            site1 = fwd_path->at(p).pos_x;
            site2 = fwd_path->at(p).pos_y;
        }
    }


    up_ident = 0;
    repeat_ident = 0;
    down_ident = 0;
    inv_ident = 0;
    fwd_ident = 0;

    inv_sum_length = 0;
    inv_sum_ins = 0;
    inv_sum_del = 0;
    inv_sum_mis = 0;

    int sum1=0;
    for(int i=0;i<(int)seq1_frag1.size();i++)
    {
        if(seq1_frag1.at(i)==seq2_frag1.at(i))
            sum1++;
        else if(seq1_frag1.at(i)<0)
            inv_sum_del++;
        else if(seq2_frag1.at(i)<0)
            inv_sum_ins++;
        else
            inv_sum_mis++;

        inv_sum_length++;
     }

    if(seq1_frag1.size()>0)
        up_ident = float(sum1)/int(seq1_frag1.size());
    else
        up_ident = 0;

    int sum2=0;
    int pState=-1;
    bool hasCG=false;
    bool hasGC=false;

    for(int i=0;i<(int)seq1_frag2.size();i++)
    {
        if(not reverse && seq1_frag2.at(i)==seq2_frag2.at(i))
            sum2++;
        else if(reverse && seq1_frag2.at(i)>=0 && seq1_frag2.at(i)<4 && seq1_frag2.at(i)==3-seq2_frag2.at(i))
            sum2++;
        else if(seq1_frag2.at(i)<0)
            inv_sum_del++;
        else if(seq2_frag2.at(i)<0)
            inv_sum_ins++;
        else
            inv_sum_mis++;

        inv_sum_length++;

        if(pState==1 && seq1_frag2.at(i)==2)
            hasCG=true;
        else if(pState==2 && seq1_frag2.at(i)==1)
            hasGC=true;

        pState=seq1_frag2.at(i);
    }
    CpG=0;
    if(hasCG)
        CpG+=1;
    if(hasGC)
        CpG+=2;

    repeat_ident = float(sum2)/int(seq1_frag2.size());

    int sum3=0;
    for(int i=0;i<(int)seq1_frag3.size();i++)
    {
        if(seq1_frag3.at(i)==seq2_frag3.at(i))
            sum3++;
        else if(seq1_frag3.at(i)<0)
            inv_sum_del++;
        else if(seq2_frag3.at(i)<0)
            inv_sum_ins++;
        else
            inv_sum_mis++;

        inv_sum_length++;
     }

    if(seq1_frag3.size()>0)
        down_ident = float(sum3)/int(seq1_frag3.size());
    else
        down_ident = 0;

    inv_ident = float(sum1+sum2+sum3)/int(seq1_frag1.size()+seq1_frag2.size()+seq1_frag3.size());
//    std::cout<<sum1<<" "<<sum2<<" "<<sum3<<" : "<<seq1_frag1.size()<<" "<<seq1_frag2.size()<<" "<<seq1_frag3.size()<<" : "<<inv_sum_del<<" "<<inv_sum_ins<<" "<<inv_sum_mis<<std::endl;

    //Forward alignment
    //
    fwd_sum_mis = 0;
    fwd_sum_ins = 0;
    fwd_sum_del = 0;

    int sum4=0;
    for(int i=0;i<(int)seq1_fwd.size();i++)
    {
        if(seq1_fwd.at(i)==seq2_fwd.at(i))
            sum4++;
        else if(seq1_fwd.at(i)<0)
            fwd_sum_del++;
        else if(seq2_fwd.at(i)<0)
            fwd_sum_ins++;
        else
            fwd_sum_mis++;
     }


/*
    std::vector<FPA::seqCoordinate>::iterator it = fwd_path->begin();
    for(;it!=fwd_path->end();it++)
    {
        if(it->pos_x<0) {
            fwd_sum_del++;
        }
        else if(it->pos_y<0) {
            fwd_sum_ins++;
        }
        else
        {
            if(seq1.at(it->pos_x-1)==seq2.at(it->pos_y-1))
                sum4++;
            else
                fwd_sum_mis++;
        }
    }
    */

    fwd_ident = float(sum4)/(fwd_sum_del+fwd_sum_ins+fwd_sum_mis+sum4);
//    std::cout<<sum4<<" "<<fwd_sum_del<<" "<<fwd_sum_ins<<" "<<fwd_sum_mis<<std::endl;


}

void FPA::print_inversion_fragment(std::vector<switchPoint> *points)
{

    std::string alpha = "-ACGTN";
    std::string lowalpha = "-acgtn";

//    std::cout<<std::setprecision(3);

//    std::cout<<points->at(0).i<<","<<points->at(0).j<<","<<points->at(1).j<<","<<points->at(2).j<<","<<points->at(3).j<<","
//        <<up_ident<<","<<repeat_ident<<","<<down_ident<<","<<inv_ident<<","<<fwd_ident<<","<<sum_ins<<","<<sum_del<<","<<sum_mis<<","<<sumNuc<<","<<CpG<<"\n";


    if(long_output)
    {
        std::cout<<"\nSwitch point 1: "<<points->at(0).j<<" ("<<points->at(0).i<<")\n";
        std::cout<<  "       point 2: "<<points->at(1).j<<" ("<<points->at(1).i<<")\n";
        std::cout<<  "       point 3: "<<points->at(2).j<<" ("<<points->at(2).i<<")\n";
        std::cout<<  "       point 4: "<<points->at(3).j<<" ("<<points->at(3).i<<")\n";

        std::cout<<"\nIdentity upstream: "<<up_ident<<std::endl;
        std::cout<<"           repeat: "<<repeat_ident<<std::endl;
        std::cout<<"       downstream: "<<down_ident<<std::endl;
        std::cout<<"        inversion: "<<inv_ident<<std::endl;
        std::cout<<"          forward: "<<fwd_ident<<std::endl;
    }

    std::stringstream qry;
    std::stringstream ref;
    for(int i=0;i<(int)seq1_frag1.size();i++)
    {
        if(low_frag1.at(i))
            qry<<lowalpha.at(seq1_frag1.at(i)+1);
        else
            qry<<alpha.at(seq1_frag1.at(i)+1);
        ref<<alpha.at(seq2_frag1.at(i)+1);
    }
    ref<<"|";
    qry<<"|";
    for(int i=0;i<(int)seq1_frag2.size();i++)
    {
        if(low_frag2.at(i))
            qry<<lowalpha.at(seq1_frag2.at(i)+1);
        else
            qry<<alpha.at(seq1_frag2.at(i)+1);
        ref<<alpha.at(seq2_frag2.at(i)+1);
    }
    ref<<"|";
    qry<<"|";
    for(int i=0;i<(int)seq1_frag3.size();i++)
    {
        if(low_frag3.at(i))
            qry<<lowalpha.at(seq1_frag3.at(i)+1);
        else
            qry<<alpha.at(seq1_frag3.at(i)+1);
        ref<<alpha.at(seq2_frag3.at(i)+1);
    }
    ref<<"\n";
    qry<<"\n";

    std::stringstream fwd_qry;
    std::stringstream fwd_ref;
    for(int i=0;i<(int)seq1_fwd.size();i++)
    {
        if(low_fwd.at(i))
            fwd_qry<<lowalpha.at(seq1_fwd.at(i)+1);
        else
            fwd_qry<<alpha.at(seq1_fwd.at(i)+1);
        fwd_ref<<alpha.at(seq2_fwd.at(i)+1);
    }
    fwd_ref<<"\n";
    fwd_qry<<"\n";

    if(long_output)
        std::cout<<"\nForward alignment:\n"<<qry_name<<" "<<fwd_qry.str()<<ref_name<<" "<<fwd_ref.str()<<std::endl;

    std::cout<<"Template-switch alignment:\n"<<qry_name<<" "<<qry.str()<<ref_name<<" "<<ref.str()<<std::endl;

}

/********************    alignment output   ******************************/



