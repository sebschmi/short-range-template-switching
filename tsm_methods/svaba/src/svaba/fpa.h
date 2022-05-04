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

#ifndef FPA_H
#define FPA_H

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <sstream>

class FPA
{
public:
    struct Fasta_entry
    {
        std::string name;
        std::string sequence;
        int length;
    };

    struct switchPoint
    {
        int i;
        int j;
    };

    struct seqCoordinate
    {
        int pos_x;
        int pos_y;
        int matrix;
    };

private:

    /******************** widely used variables ******************************/

    std::vector<int> index1;
    std::vector<int> index2;
    std::vector<int> rindex1;
    std::vector<int> rindex2;

    std::vector<int> fseq1;
    std::vector<int> fseq2;

    std::vector<int> seq1;
    std::vector<int> rev1;
    std::vector<int> seq2;

    std::vector<int> mask1;

    int slg;
    int fsl1;
    int fsl2;


    int sl1;
    int sl2;
    int start1;
    int end1;
    int start2;
    int end2;
    int true_start2;
    int true_end2;


    std::string qry_name;
    std::string ref_name;


    template<typename T>
    struct Array2D
    {
        private:
            int width;
            int org_width;
            int org_height;
            T * data;
        public:
            T& operator() (int x, int y) { return data[y*width + x]; }
            Array2D(const int w, const int h) : width(w), org_width(w), org_height(h) { data = new T[w*h]; }
            void resize(int nw, int nh) { if(nw*nh <= org_width*org_height) { width=nw; } else { delete [] data; data = new T[nw*nh]; org_width = nw; org_height = nh; } }
            ~Array2D() { delete [] data; }
    };

    enum Move_ptr {match=-1, xgap=-2, ygap=-3, none=-4};

    /******************** widely used variables ******************************/

    /******************** command-line argument ******************************/

    bool maximise_score = false;
    bool maximise_length = false;

    bool reverse = false;
    bool force_overlap = false;
    bool perfect_copy = false;

    bool long_output = false;
    bool debug_matrix = false;

    /******************** command-line argument ******************************/

    float up_ident = 0;
    float repeat_ident = 0;
    float down_ident = 0;
    float inv_ident = 0;
    float fwd_ident = 0;

    int inv_sum_length = 0;
    int inv_sum_ins = 0;
    int inv_sum_del = 0;
    int inv_sum_mis = 0;

    int CpG=0;
    int fwd_sum_mis = 0;
    int fwd_sum_ins = 0;
    int fwd_sum_del = 0;
    int sumNuc = 0;

    std::vector<int> seq1_frag1;
    std::vector<int> seq2_frag1;
    std::vector<bool> low_frag1;

    std::vector<int> seq1_frag2;
    std::vector<int> seq2_frag2;
    std::vector<bool> low_frag2;

    std::vector<int> seq1_frag3;
    std::vector<int> seq2_frag3;
    std::vector<bool> low_frag3;

    std::vector<int> seq1_fwd;
    std::vector<int> seq2_fwd;
    std::vector<bool> low_fwd;

    int substitution_score(int i, int j)
    {
        if(seq1.at(i-1)<4 && seq1.at(i-1) == seq2.at(j-1))
            return 1;
        else
            return -1;
    }

    int rev_substitution_score(int i, int j)
    {
        if(perfect_copy)
        {
            if(rev1.at(i-1)<4 && rev1.at(i-1) == seq2.at(j-1))
                return 1;
            else
                return -10000;
        }
        else
        {
            if(rev1.at(i-1)<4 && rev1.at(i-1) == seq2.at(j-1))
                return 1;
            else
                return -1;
        }
    }

    void fwd_compare_sequences(int *identical,int *length,int epo_start1,int epo_stop1)
    {
        *identical = 0;
        *length = 0;

        for(int i=rindex1.at(epo_start1);i<rindex1.at(epo_stop1);i++)
        {
            if(fseq1.at(i) == fseq2.at(i))
                (*identical)++;
            (*length)++;
        }
    }

//    void build_indeces(std::string *s1,std::string *s2);

public:
    FPA() {}

    void read_data(std::string *datafile,Fasta_entry *entry);
    void read_pair_data(std::string *datafile, Fasta_entry *first, Fasta_entry *second, bool swap_pair);

    void build_sequences(std::string *s1, std::string *s2, bool reverse=false);

    void fwd_align_sequences(std::vector<seqCoordinate> *path);
    void align_sequences(std::vector<seqCoordinate> *path,std::vector<switchPoint> *points,bool local=true);

    void compute_fpa_score(std::vector<seqCoordinate> *path, std::vector<switchPoint> *points, std::vector<seqCoordinate> *fwd_path, int flanking=20);

    void print_switch_process(std::vector<seqCoordinate> *path,std::vector<switchPoint> *points);
    void print_inversion_fragment(std::vector<switchPoint> *points);
//    void print_events(std::vector<seqCoordinate> *path, std::vector<switchPoint> *points, std::vector<seqCoordinate> *fwd_path, std::string line);

    void set_fpa_long_output(bool i) { long_output= i; }
    void set_fpa_debug(bool i) { debug_matrix = i; }

    bool get_fpa_long_output() { return long_output; }
    bool get_fpa_debug() { return debug_matrix; }

    int get_qry_sp1(std::vector<switchPoint> *points) { return points->at(0).i; }
    int get_qry_sp2(std::vector<switchPoint> *points) { return points->at(1).i; }
    int get_qry_sp3(std::vector<switchPoint> *points) { return points->at(2).i; }
    int get_qry_sp4(std::vector<switchPoint> *points) { return points->at(3).i; }
    int get_ref_sp1(std::vector<switchPoint> *points) { return points->at(0).j; }
    int get_ref_sp2(std::vector<switchPoint> *points) { return points->at(1).j; }
    int get_ref_sp3(std::vector<switchPoint> *points) { return points->at(2).j; }
    int get_ref_sp4(std::vector<switchPoint> *points) { return points->at(3).j; }
    int get_inv_length(std::vector<switchPoint> *points) { return points->at(1).j-points->at(2).j+1; }
    float get_up_ident() { return up_ident; }
    float get_repeat_ident() { return repeat_ident; }
    float get_down_ident() { return down_ident; }
    float get_inv_ident() { return inv_ident; }
    float get_fwd_ident() { return fwd_ident; }
    int get_fwd_sum_ins() { return fwd_sum_ins; }
    int get_fwd_sum_del() { return fwd_sum_del; }
    int get_fwd_sum_mis() { return fwd_sum_mis; }
    int get_inv_sum_ins() { return inv_sum_ins; }
    int get_inv_sum_del() { return inv_sum_del; }
    int get_inv_sum_mis() { return inv_sum_mis; }
    int get_sum_nuc() { return sumNuc; }
    int get_cpg() { return CpG; }

};

#endif // FPA_H
