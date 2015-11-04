//  BRAT is a whole genome bisulfite sequence mapping program
//  Copyright (C) 2009  Elena Yavorska Harris
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
// 
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
/*
Same seeds as in Brat-1.1.4, but less space since it treats chroms as one continuous genome
  */

#include<iomanip>
#include<map>
#include<bitset>
#include<ctime>
#include<list>
#include<cmath>
#include<vector>
#include<string.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<assert.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

using namespace std;

//lookup8.c, by Bob Jenkins, January 4 1997, Public Domain.
#define mix64(a,b) \
{ \
  a -= b;  a ^= (b>>43); \
  b -= a; b ^= (a<<9); \
  a -= b; a ^= (b>>38); \
  b -= a; b ^= (a<<23); \
  a -= b; a ^= (b>>35); \
  b -= a; b ^= (a<<49); \
  a -= b; a ^= (b>>12); \
  b -= a; b ^= (a<<18); \
}

/*
#define mix64(a,b,c) \
{ \
  a -= b; a -= c; a ^= (c>>43); \
  b -= c; b -= a; b ^= (a<<9); \
  c -= a; c -= b; c ^= (b>>8); \
  a -= b; a -= c; a ^= (c>>38); \
  b -= c; b -= a; b ^= (a<<23); \
  c -= a; c -= b; c ^= (b>>5); \
  a -= b; a -= c; a ^= (c>>35); \
  b -= c; b -= a; b ^= (a<<49); \
  c -= a; c -= b; c ^= (b>>11); \
  a -= b; a -= c; a ^= (c>>12); \
  b -= c; b -= a; b ^= (a<<18); \
  c -= a; c -= b; c ^= (b>>22); \
}
*/
//End lookup8.c, by Bob Jenkins, January 4 1997, Public Domain.


class Stat{
public:


	long int mapped_pairs;
	long int ambiguous;
	long int unmapped;
	long int single_mapped;
	long int single_ambiguous;
	long int broken_dif_chroms;
	long int broken_wrong_insert;
	long int broken_inverse;
	long int broken_mate_unmapped;
	long int broken_mate_amb; //if pair is unmapped and mate is ambiguous, means something else is wrong
	//dif chroms, inverse or wrong insert, otherwise it would be counted as amb pair
	long int broken_both_unmapped;
	long int invalid;
	Stat();
};
Stat::Stat(){

	mapped_pairs = 0; 
	ambiguous = 0; 
	unmapped = 0; 
	single_mapped = 0;
	single_ambiguous = 0;
	 broken_dif_chroms = 0;
	 broken_wrong_insert = 0;
	 broken_inverse = 0;
	 broken_mate_unmapped = 0;
	 broken_mate_amb = 0; //if pair is unmapped and mate is ambiguous, means something else is wrong
	//dif chroms, inverse or wrong insert, otherwise it would be counted as amb pair
	 broken_both_unmapped = 0;
	invalid = 0;
}
class Mates{
public:
	Mates(unsigned int x, unsigned int y, unsigned int z, char st, short int mark) : chrom_id(x), pos1(y), pos2(z), strand(st), mark_unique(mark), mism1(1000), mism2(1000), min_mism(1000) {}; 
	Mates(): chrom_id(-1), pos1(0), pos2(0), strand('+'), mark_unique(0), mism1(100), mism2(100), min_mism(100) {};
	void operator =(Mates a){ chrom_id = a.chrom_id; pos1 = a.pos1 ; pos2 = a.pos2; strand = a.strand; mark_unique = a.mark_unique; };
	unsigned int chrom_id;
	unsigned int pos1;
	unsigned int pos2;
	char strand;
	short int mark_unique;
	unsigned short mism1: 8;
	unsigned short mism2: 8;
	unsigned short min_mism: 8;

};

long int Reverse(unsigned long int i, int bits){//by Nicholas Allen

i = ((i >> 1) & 0x5555555555555555) | ((i & 0x5555555555555555) << 1);
i = ((i >> 2) & 0x3333333333333333) | ((i & 0x3333333333333333) << 2);
i = ((i >> 4) & 0x0F0F0F0F0F0F0F0F) | ((i & 0x0F0F0F0F0F0F0F0F) << 4);
i = ((i >> 8) & 0x00FF00FF00FF00FF) | ((i & 0x00FF00FF00FF00FF) << 8);
i = ((i >> 16) & 0x0000FFFF0000FFFF) | ((i & 0x0000FFFF0000FFFF) << 16);
i = ( i >> 32 ) | ( i << 32);
i = i >> (64 - bits);

return i;
}
void convert_lmer_ta_long(string lmer, vector<unsigned long int> &ta_read, int asize){
	long int j, i;
	j = 0;
	while( j < asize){
		unsigned long int gb = 0; 
		for(i = 0; i < 64 && j < asize; i++){
			gb <<= 1;
			if(lmer[j] == 'T' || lmer[j] == 't' ||  lmer[j] == 'C' || lmer[j] == 'c')
				gb |= 1;
			j++;
		}//for i
		
		ta_read.push_back(gb);

	}//while j
}//convert_ta

void convert_lmer_c_long(string current_chrom, vector<unsigned long int> &cg_read,
						 int num_iter)
{//this function is used only once when there is a beginning of a chromosome
	//aread is an lmer in a genome
	//gb is binary representative of this genome lmer
	long int j =0, i = 0;

	while(i < num_iter){
		unsigned long int gb = 0;

		for(j = 0; j < 64 && i < num_iter; j++){
			gb <<= 1;
			if(current_chrom[i] == 'C' || current_chrom[i] == 'c' ){
				gb |= 1;
			}
			else if(current_chrom[i] == 'G' || current_chrom[i] == 'g'){
				gb |= 1;
			}
			i++;
		}//for j
		cg_read.push_back(gb);
	}//while

}//connvert_lmer_c()


void convert_lmer_ta(string lmer, unsigned long int &gb, unsigned long int &ns, int asize){
	long int j;
	gb = 0; 
	ns = 0;
	for(j = 0; j < asize; j++){
		gb <<= 1;
		ns <<= 1;

		if(lmer[j] == 'T' || lmer[j] == 't' ||  lmer[j] == 'C' || lmer[j] == 'c')
			gb |= 1;
		else if(lmer[j] == 'N' || lmer[j] == 'n')
			ns |= 1;

	}//for j
}//convert_ta

void next_lmer_ta(char ch, unsigned long int &gb, unsigned long int &ns, unsigned long int MASK_READ){

	gb <<= 1;
	ns <<= 1;

	if(ch == 'T' || ch == 't'  || ch == 'C' || ch == 'c')
		gb |= 1;
	else if(ch == 'N' || ch == 'n')
		ns |= 1;

	gb &= MASK_READ;
	ns &= MASK_READ;

}//next_lmer_ta()


long int str_to_int(string s)
{
	istringstream is(s);
	long int result = -1;
	is >> result;
	return result;
}

string int_to_str(long int x)
{
	stringstream outstr;
	outstr << x;
	return outstr.str();

}
string double_to_str(double x)
{
	stringstream outstr;
	outstr << x;
	return outstr.str();

}

long int max(long int x, long int y)
{	
	if(x > y)
		return x;
	else
		return y;

}
unsigned long int min(unsigned long int x, unsigned long int y){ 
	if(x < y)
		return x;
	else
		return y;
}



void convert_lmer_c(string current_chrom, unsigned long int &gb, int num_iter)
{//this function is used only once when there is a beginning of a chromosome
	//aread is an lmer in a genome
	//gb is binary representative of this genome lmer
	long int j ;

	gb = 0;
		for(j = 0; j < num_iter; j++){
			gb <<= 1;
			if(current_chrom[j] == 'C' || current_chrom[j] == 'c' ){
				gb |= 1;
			}
			else if(current_chrom[j] == 'G' || current_chrom[j] == 'g'){
				gb |= 1;
			}

		}//for j

}//connvert_lmer_c()

void next_lmer_cg(char ch, unsigned long int &gb,  unsigned long int MASK_READ){

	gb <<= 1;

	if(ch == 'G' || ch == 'g'  || ch == 'C' || ch == 'c')
		gb |= 1;

	gb &= MASK_READ;

}//next_lmer_ta()
void next_lmer(char ch, unsigned long int &gb, unsigned long int &ns, 
			   unsigned long int &cg, unsigned long int MASK_READ)
{
	gb <<= 1;
	ns <<= 1;
	cg <<=1;

	if(ch == 'T' || ch == 't')  
		gb |= 1;
	else if(ch == 'G' || ch == 'g'){
		cg |= 1;
	}
	else if(ch == 'C' || ch == 'c'){
		gb |= 1;
		cg |= 1;
	}
	else if(ch == 'N' || ch == 'n')
		ns |= 1;

	gb &= MASK_READ;
	cg &= MASK_READ;

	ns &= MASK_READ;

}//next_lmer()


void usage(){

	string mes("\nUSAGE: brat -r <references> -s <input-reads> -o <output> [OPTIONS]\n\nOptions:\n");
	string gen_ref_opt("  -r <references-file>       file with file names containing references (one\n                          reference per file)");
	string singles(    "  -s <input-reads-file>      file with single reads (or queries)");
	string output(     "  -o <output-file>");
	string is_pair(    "  -pe                        set Pair-End mapping (if not specified,\n                          mapping for single reads is done)");
	string is_bs(      "  -bs                        set BS-mapping (if not specified,\n                          normal mapping is done)");
	string first(      "  -1 <paired-ends-file1>     mates 1");
	string second(     "  -2 <paired-ends-file2>     mates 2");
	string min_insert( "  -i <min>                   min insert");
	string max_insert( "  -a <max>                   max insert");
	string arich(      "  -A                         singles that are pair2 file,\n                          that are mapped to an A-rich strand");
	string unmapped(   "  -u                         output unmapped reads/pairs");
	string out_amb(    "  -M                         output ambiguous reads/pairs");
	string mism(	   "  -m <non-negative integer>  the number of non-BS-mismatches allowed");
	string first_bits( "  -f <integer >= 24>         the first <int> bases, within which only one non-BS-mismatch is allowed");
	cout << mes << endl;
	cout << gen_ref_opt << endl << singles << endl
        	<< output << endl
		<< is_pair << endl << is_bs << endl
		<< first << endl << second << endl << min_insert << endl << max_insert
		<< endl << arich << endl << unmapped << endl << out_amb << endl 
		<< mism << endl << first_bits << endl << endl;

	string mes2(      "Limitations:\n");
	string ref_size(  "  a single reference size <= 2^45");
	string chrom_ids( "  number of references <= 2^18");
	string space(     "  space <= (total references sizes)*4.5*8 + 16*(2^24)");
	
	cout << mes2 << endl;
	cout << ref_size << endl << chrom_ids << endl << space << endl;

	string mes3(      "Output Format:\n");
	string ri(        "  read id <int>           the order of a read/pair in the input file (starts\n                          with 0)");
	string r1(        "  read1 <string>          mate 1 as in the input file");
	string r2(        "  read2 <string>          (pair-end only) mate 2 as in the input file");
	string chr(       "  reference name <string>");
	string str(       "  strand                  (+) if mate 1 is found in forward strand,\n                          (-) if mate 1 is in reverse strand");
	string pos1(      "  position1 <int>         mate 1 leftmost position in the reference");
	string pos2(      "  position2 <int>         (pair-end only) mate 2 leftmost pos. in the reference");
	string mism1(     "  mism1 <int>             number of mismatches in mate 1, 5' mate");
	string mism2(     "  mism2 <int>             (pair-end only) number of mismatches in mate 2, 3' mate");	
	cout << endl << mes3 << endl;
	cout << ri << endl << r1 << endl << r2 << endl << chr << endl << str << endl
		<< pos1 << endl << pos2 << endl << mism1 << endl << mism2 << endl;
}//usage()


struct gen_info
{
	gen_info() : pos(0), chrom(0) {}
	~gen_info(){}
	gen_info(unsigned int p, unsigned short int c) : pos(p), chrom(c) {}
	void operator =(gen_info a){ chrom = a.chrom; pos = a.pos; };

	unsigned int pos ;
	unsigned short int chrom ;
};

class Anchor{
public:
	Anchor(int argc, char* argv[]);
	~Anchor();

	vector<string> names;
	
	vector< unsigned int* > pos_index;
	vector<long int> pos_sizes;

	vector< vector<unsigned long int> > seeds;
	vector< vector<unsigned long int> > seeds2;
	vector<unsigned long int> gen_seeds;
	vector<unsigned long int> gen_seeds2;

	void make_seeds();

	unsigned long int len_mer;
	unsigned long int read_len;
	vector<unsigned long int> chrom_blocks;
	vector<unsigned long int> size_chrom;
	vector<unsigned long int> starts_size_chrom;
	map<int, string> chroms_names;
	vector<unsigned long int> chroms_starts;//for ta/cg/ns-genomes only
	vector<unsigned long int> starts_chroms_starts;

	long int binary_total;//size of binary representation of the genome
	
	unsigned long int* ta_genome;
	unsigned long int* cg_genome;
	unsigned long int* ns_genome; //presence of Ns

	unsigned long int hash_table_size;

	int parse_options(int argc, char* argv[]);

	void read_genome();

	void find_pos_sizes(int t);
	void make_pos_index(int t);

	void print(Stat &astat);

	vector<Mates> mates_pairs;


	void convert_gen_mer_ta(unsigned long int &ta_gen_lmer,
					 unsigned long int gen_pos,
					 unsigned long int aread_len);
	void convert_gen_mer_cg(unsigned long int &ta_gen_lmer,
					 unsigned long int gen_pos,
					 unsigned long int aread_len);
	void convert_gen_mer_ns(unsigned long int &ta_gen_lmer,
					 unsigned long int gen_pos,
					 unsigned long int aread_len);

	void convert_gen_mer_ta_long(vector<unsigned long int> &ta_gen_lmer,
					 unsigned long int gen_pos,
					 unsigned long int aread_len, int asize);
	void convert_gen_mer_cg_long(vector<unsigned long int> &ta_gen_lmer,
					 unsigned long int gen_pos,
					 unsigned long int aread_len, int asize);

	void convert_gen_mer_ns_long(vector<unsigned long int> &ta_gen_lmer,
					 unsigned long int pos,
					 unsigned long int aread_len, int asize);

	int count_mism(unsigned long int bm);
	int count_mism_long(vector<unsigned long int> &reads, vector<unsigned long int> &genmers,
		vector<unsigned long int> &ta_check, int asize);

	void pair_reads(unsigned long int read_id, string aread1, string aread2,
				long int &contain_Ns, Stat &astat, short int len1, short int len2, int t, int rs);

	void pair_reads_long(unsigned long int read_id, string aread1, string aread2,
				long int &contain_Ns, Stat &astat, short int len1, short int len2, int t, int rs);

	void pair_reads_normal(unsigned long int read_id, string aread1, string aread2,
				long int &contain_Ns, Stat &astat, short int len1, short int len2, int t, int rs);
	void pair_reads_normal_long(unsigned long int read_id, string aread1, string aread2,
				long int &contain_Ns, Stat &astat, short int len1, short int len2, int t, int rs);

	void pair_reads_singles(unsigned long int read_id, string aread1,
				long int &contain_Ns, Stat &astat, short int len1, int t, int rs);
	void pair_reads_singles_long(unsigned long int read_id, string aread1,
				long int &contain_Ns, Stat &astat, short int read_len1,  int t, int rs );

	void pair_reads_singlesA(unsigned long int read_id, string aread1,
				long int &contain_Ns, Stat &astat, short int len1, int t, int rs);
	void pair_reads_singlesA_long(unsigned long int read_id, string aread1,
				long int &contain_Ns, Stat &astat, short int len1, int t, int rs);

	void pair_reads_normal_singles(unsigned long int read_id, string aread1,
				long int &contain_Ns, Stat &astat, short int len1, int t, int rs);

	void pair_reads_normal_singles_long(unsigned long int read_id, string aread1,
				long int &contain_Ns, Stat &astat, short int len1, int t, int rs);


	string genome_names;
	long int min_insert;
	long int max_insert;
	unsigned long int MASK_BYTE;
	unsigned long int MASK_READ;
	unsigned long int MASK_LMER;
	unsigned long int MASK_HALF;
	unsigned long int MASK_SEED;
	unsigned long int byte_size;
	unsigned long int abyte ;
	unsigned long int SIZE ;//genome char per an unsigned long int
	unsigned long int two_in_Size; //log_2(SIZE)
	unsigned long int MASK_SIZE;
	unsigned long int block_size ;
	unsigned long int double_block ;//2 end of lines
	unsigned long int two_in_byte ;//log_2(byte_size)
	unsigned long int block_byte_prod_divide ;//=log_2(256*8)
	unsigned long int block_divide ;//log_2(256)
	unsigned long int bytes;//number of bytes per unsigned long int

	int mism_thresh;
	int first_bits;//within these first bits only one non-BS-mism is allowed

	unsigned long int MASK_chrom;
	ofstream out_pairs;

	string is_pe;
	unsigned long int cur_block;
	long int num_chroms;
	string output_pairs;
	unsigned long int CHROMS_BITS;
	string reads_file1;
	string reads_file2;
	string is_bs;
	bool cascade_amb;
	bool A_rich;
	bool unmapped;
	ofstream out_unmapped1;
	ofstream out_unmapped2;

	long int total_reads;

	vector<short int> reads_len1;
	vector<short int> reads_len2;

	void release_space();
	unsigned long int dif_len_mer;

	bool out_amb;

	unsigned long int hash2m(unsigned long int ta);
	unsigned long int hash2c(unsigned long int ta, unsigned long int cg);
	unsigned long int hash2gs(unsigned long int ta, unsigned long int cg);
	unsigned long int hash2cs(unsigned long int ta, unsigned long int cg);
	unsigned long int hash_aread(unsigned long int ta, unsigned long int cg, 
		bool is_tac);

};	

void Anchor::release_space()
{
	long int i;
	for(i = 0; i < hash_table_size; i++)
		delete [] pos_index[i];

}//release pos index
Anchor::~Anchor(){

	long int i;

		delete [] ta_genome;
		delete [] cg_genome;
		delete [] ns_genome;

	out_unmapped1.close();
	out_unmapped2.close();
}

int Anchor::count_mism(unsigned long int bits){

/////////////////Andrew Smith author of the following code:
  bits = ((bits & 0xAAAAAAAAAAAAAAAA) >> 1)  + (bits & 0x5555555555555555);
  bits = ((bits & 0xCCCCCCCCCCCCCCCC) >> 2)  + (bits & 0x3333333333333333);
  bits = ((bits & 0xF0F0F0F0F0F0F0F0) >> 4)  + (bits & 0x0F0F0F0F0F0F0F0F);
  bits = ((bits & 0xFF00FF00FF00FF00) >> 8)  + (bits & 0x00FF00FF00FF00FF);
  bits = ((bits & 0xFFFF0000FFFF0000) >> 16) + (bits & 0x0000FFFF0000FFFF);
  return ((bits & 0xFFFFFFFF00000000) >> 32) + (bits & 0x00000000FFFFFFFF);
/////////////////////END Andrew Smith

}//count_mism()

int Anchor::count_mism_long(vector<unsigned long int> &read, vector<unsigned long int> &gen_mer,
							vector<unsigned long int> &ta_check, int asize){

	int i, mism = 0;
	for(i = 0; i < asize; i++){
		unsigned long int res = read[i] ^ gen_mer[i];
		ta_check.push_back(res);
		mism += count_mism(res);
	}//for i
	
	return mism;
}//count_mism_long()

Anchor::Anchor(int argc, char* argv[])
{
	first_bits = 36;
	out_amb = false;
	unmapped = false;
	A_rich = false;
	mism_thresh = 0;
	min_insert = 100;
	max_insert = 300;
	is_bs = "0";
	is_pe = "0";
	int res = parse_options(argc, argv);
	if(res == -1){
		usage();
		exit(0);
	}
	if(is_bs != "1" && A_rich == true){
		cout << "\nOption A works only with BS-mapping. Ignoring option A" << endl;
		A_rich = false;
	
	}
	string unm1 =  reads_file1 + ".unm";
	string unm2 = reads_file2 + ".unm";
	if(unmapped == true){
		out_unmapped1.open(unm1.c_str(), ios::out);
		if(is_pe == "1"){
			out_unmapped2.open(unm2.c_str(), ios::out);
		}
	}//if unm
	int min_len = 10000, cur_len;

	ifstream in;
	long int st, en;
	in.open(reads_file1.c_str(), ios::in);
	if(!in){
		cout << "\nERROR: could not open " << reads_file1 << endl;
		exit(0);
	}
	string aread;
	in >> aread;

		if(aread.length() == 0){
			cout << "\nERROR: no reads in " << reads_file1 << endl;
			in.close();
			exit(0);
		}
		if(aread[0] == '@' || aread[0] == '>'){
			cout << "\nERROR: Wrong format of reads file: only one read per a line format is supported." << endl;
			in.close();
			exit(0);
		}
	total_reads = 0;
	while(!in.eof()){
		in >> st;
		in >> en;
		total_reads++;
		in >> aread;
	
	}
	in.close(); in.clear();

	mates_pairs.resize(total_reads);
	reads_len1.resize(total_reads);

	in.open(reads_file1.c_str(), ios::in);
	if(!in){
		cout << "\nERROR: could not open " << reads_file1 << endl;
		exit(0);
	}
	long int i;

	in >> aread;
	long int read_id = 0;
	while(!in.eof()){
		in >> st;
		in >> en;
		cur_len = aread.length();
		reads_len1[read_id] = cur_len;
		if(cur_len < 24 )
			mates_pairs[read_id].mark_unique = -1;
		else{
			if(min_len > cur_len)
				min_len = cur_len;
		}
		int Ns = 0;
		for(i = 0; i < cur_len; i++){
			if(aread[i] == 'N' || aread[i] == 'n')
				Ns++;
		}
		if(Ns > mism_thresh)
			mates_pairs[read_id].mark_unique = -1;

		read_id++;
		in >> aread;
	
	}
	in.close(); in.clear();


	if(is_pe == "1"){
		reads_len2.resize(total_reads);

		in.open(reads_file2.c_str(), ios::in);
		if(!in){
			cout << "\nERROR: could not open " << reads_file2 << endl;
			exit(0);
		}


	in >> aread;
		if(aread.length() == 0){
			cout << "\nERROR: no reads in " << reads_file1 << endl;
			in.close();
			exit(0);
		}
		if(aread[0] == '@' || aread[0] == '>'){
			cout << "\nERROR: Wrong format of reads file: only one read per a line format is supported." << endl;
			in.close();
			exit(0);
		}
	read_id = 0;
	while(!in.eof()){
		in >> st;
		in >> en;
		cur_len = aread.length();
		reads_len2[read_id] = cur_len;
		if(cur_len < 24 )
			mates_pairs[read_id].mark_unique = -1;
		else{
			if(min_len > cur_len)
				min_len = cur_len;
		}
		int Ns = 0;
		for(i = 0; i < cur_len; i++){
			if(aread[i] == 'N' || aread[i] == 'n')
				Ns++;
		}
		if(Ns > mism_thresh)
			mates_pairs[read_id].mark_unique = -1;
		
		read_id++;
		in >> aread;
	
	}
	in.close(); in.clear();

	
	
	
	}//is pe

	CHROMS_BITS = 18;
	out_pairs.open(output_pairs.c_str(), ios::out);
	unsigned long int one = 1;

	abyte = 8;
	SIZE = 64;//genome char per an unsigned long int
	two_in_Size = 6;
	MASK_SIZE = (one << two_in_Size) - 1;
	MASK_HALF = (one << 32) -1 ;
	block_size = 256;
	double_block = 256 + 256 + 2;//2 end of lines
	byte_size = 8;//byte size is characters(genome) per a byte
	two_in_byte = 3;//log_2(byte_size)

	block_byte_prod_divide = 11;//=log_2(256*8)
	block_divide = 8;//log_2(256)
	bytes = SIZE >> two_in_byte;//number of bytes per unsigned long int
	MASK_chrom = (one << CHROMS_BITS) - 1;

	MASK_SEED = (one << 24) -1;
	len_mer = min(min_len, SIZE);

	if(mism_thresh > 0)
		len_mer = min(len_mer, first_bits);
	len_mer = max(len_mer, 24);

	unsigned long int ta_ind = 0;

	float base = 2;
	float apow = 24;

	hash_table_size = pow(base, apow);

	MASK_BYTE = (one << (byte_size)) -1;
	if(len_mer == SIZE)
		MASK_LMER = 0xFFFFFFFFFFFFFFFF;
	else if(len_mer < 32)
		MASK_LMER = (one << len_mer) - 1;//seed mask
	else{
		int dist = len_mer - 32;
		MASK_LMER = (((one << dist) - 1) << 32) | 0x00000000FFFFFFFF;
		
	}

	make_seeds();
	dif_len_mer = SIZE - len_mer;


	unsigned long int u, j;

	in.open(genome_names.c_str(), ios::in);
	if(!in){
		cout << "ERROR: cannot open " << genome_names << endl;
		exit(0);
	}

	string aname;
	num_chroms = 0;
	in >> aname;
	while(!in.eof()){
	
		if(aname.length() == 0){
			cout << "\nERROR: No input in the file " << genome_names << endl;
			exit(0);
		}
		if(aname[0] == '>' || aname[0] == '@'){
			cout << "\nERROR: " << genome_names << " must contain the names of the files with references." << endl;
			exit(0);
		}
		names.push_back(aname);
		num_chroms++;
		in >> aname;
	}//while
	in.close(); in.clear();

	string line;

	//find out sizes of chromosomes (or references)
//	ta_genome.resize(1);

	long int cur_size;
	for(u = 0; u < num_chroms; u++){
		in.open(names[u].c_str(), ios::in);
		if(!in){
			cout << "ERROR: cannot open " << names[u] << endl;
			exit(0);
		}
		cur_size = 0;

		getline(in, line);

		if(line.length() == 0){
			cout << "\nERROR: No input in the file " << names[u] << endl;
			exit(0);
		}
		if(line[0] != '>' ){
			cout << "\nERROR: " << names[u] << " must be in FASTA format." << endl;
			exit(0);
		}

		getline(in, line);
		long int line_len = line.length();
		while(!in.eof()){
			cur_size++;

			getline(in, line);
		}//while

		cur_size = cur_size * line_len;

		unsigned long int rem = SIZE - (MASK_SIZE & cur_size);
		cur_size += rem ;//now cur_size is divisible by SIZE

		//add max_insert between chroms so that we don't map mates to different chroms
		cur_size += ((max_insert >> two_in_Size) + 1) * SIZE;

		size_chrom.push_back(cur_size);
		long int new_size = (cur_size >> two_in_Size) ;// =/(8*64) 1pb takes 1 bit, and 1 long takes 64 bits
		
		chroms_starts.push_back(new_size);


		in.close(); in.clear();
	}


	unsigned long int cur_total = 0;
	unsigned long int starts_total = 0;
	for(i = 0; i < size_chrom.size(); i++){
		starts_size_chrom.push_back(cur_total);
		cur_total += size_chrom[i];

		starts_chroms_starts.push_back(starts_total);
		starts_total += chroms_starts[i];
	}
	starts_chroms_starts.push_back(starts_total);//the end of the last chrom

	cur_total = cur_total >> two_in_Size;
	
	binary_total = cur_total;

	ta_genome = new unsigned long int[cur_total];
	cg_genome = new unsigned long int[cur_total];
	ns_genome = new unsigned long int[cur_total];

	assert(ta_genome != 0);
	assert(cg_genome != 0);
	assert(ns_genome != 0);


	read_genome();

	pos_sizes.resize(hash_table_size);
	pos_index.resize(hash_table_size);


	MASK_READ = MASK_LMER;//for initialization only

}//Anchor()

//lookup8.c, by Bob Jenkins, January 4 1997, Public Domain.


unsigned long int Anchor::hash2m(unsigned long int ta)
{
	// in ta ones are Ts and Cs, similarly in cg ones are Cs and Gs, and in tg Ts and Gs
	unsigned long int a,b;

	a = 0x9e3779b97f4a7c13 ;//LL; 
	b = 0x9e3779b97f4a7c13;//LL; 
    a += ta & MASK_HALF;//mask least significant 32 bits
    b += ta >> 32;
    mix64(a,b);

  return ( b & MASK_SEED);
}

unsigned long int Anchor::hash2c(unsigned long int ta, unsigned long int cg)
{
	// in ta ones are Ts and Cs, similarly in cg ones are Cs and Gs, and in tg Ts and Gs
	unsigned long int a,b;

	a = 0x9e3779b97f4a7c13 ;//LL; 
	b = 0x9e3779b97f4a7c13;//LL; 

	//Gs info
    a += cg;//mask least significant 32 bits
    b += ta;//Cs info only
    mix64(a,b);

  return ( b & MASK_SEED);
}

unsigned long int Anchor::hash2gs(unsigned long int ta, unsigned long int cg)
{
	// in ta ones are Ts and Cs, similarly in cg ones are Cs and Gs, and in tg Ts and Gs
	unsigned long int a,b;

	a = 0x9e3779b97f4a7c13 ;//LL; 
	b = 0x9e3779b97f4a7c13;//LL; 

	//Gs info
    a += ((~ta) & MASK_LMER) & cg;//mask least significant 32 bits
    b += ta;//Cs info only
    mix64(a,b);

  return ( b & MASK_SEED);
}

unsigned long int Anchor::hash2cs(unsigned long int ta, unsigned long int cg)
{
	
	// in ta ones are Ts and Cs, similarly in cg ones are Cs and Gs, and in tg Ts and Gs
	unsigned long int a,b;

	a = 0x9e3779b97f4a7c13 ;//LL; 
	b = 0x9e3779b97f4a7c13;//LL; 

	//Gs info
    a += ta & cg;//mask least significant 32 bits
    b += ta;//Cs info only
    mix64(a,b);

  return ( b & MASK_SEED);
}
//End lookup8.c, by Bob Jenkins, January 4 1997, Public Domain.



unsigned long int Anchor::hash_aread(unsigned long int ta, unsigned long int cg, 
						 bool is_tac)
{
//is_tac = true, then call hash2cs; else call hash2gs

	if(mism_thresh > 0)
		return(hash2m(ta));
	else if(is_tac == true)
		return hash2cs(ta, cg);
	else
		return hash2gs(ta, cg);

}//hash_aread

void Anchor::make_seeds(){
	
	unsigned long int one =1 , MASK_quarter, MQ1, MQ2, MQ3, MASK_half_lenmer, new_len_mer;
	long int quarter, i, ahalf;
	quarter = len_mer/3;
	ahalf = len_mer >> 1;//divide by 2
	MASK_quarter = (one << quarter) - 1;
	long int double_quarter = quarter + quarter;
	MQ1 = MASK_quarter << double_quarter ;
	MQ2 = MASK_quarter << quarter;

	if(is_pe == "1"){
		if(mism_thresh == 0){
			gen_seeds.push_back(MASK_LMER);
			gen_seeds2.push_back(MASK_LMER);

			vector<unsigned long int> dum;
			dum.push_back(MASK_LMER);
			seeds.push_back(dum);
			seeds2.push_back(dum);
		}
		else{//mism > 0 and pairs


				new_len_mer = double_quarter + quarter;
				int shift = len_mer -  new_len_mer;
				MASK_LMER = MASK_LMER >> shift;

				gen_seeds.resize(4);
				gen_seeds2.resize(4);

				gen_seeds[0] = MASK_LMER;
				gen_seeds[1] = MQ1 | MQ2;
				gen_seeds[2] = MQ1 | MQ2;
				gen_seeds[3] = MQ1 | MASK_quarter;

				gen_seeds2[0] = MASK_LMER;
				gen_seeds2[1] = MQ2 | MASK_quarter;
				gen_seeds2[2] = MQ1 | MASK_quarter;
				gen_seeds2[3] = MQ2 | MASK_quarter;

				seeds.resize(4);
				seeds2.resize(4);

				seeds[0].push_back(MASK_LMER);
				seeds[1].push_back(MQ1 | MQ2);
				seeds[1].push_back(MQ2 | MASK_quarter);
				seeds[1].push_back(MQ1 | MQ2);
				seeds[1].push_back(MQ2 | MASK_quarter);

				seeds[2].push_back(MQ1 | MQ2);
				seeds[2].push_back(MQ1 | MASK_quarter);
				seeds[2].push_back(MQ1 | MASK_quarter);

				seeds[3].push_back(MQ2 | MASK_quarter);
				seeds[3].push_back(MQ1 | MASK_quarter);



				seeds2[0].push_back(MASK_LMER);
				seeds2[1].push_back(MQ2 | MASK_quarter);
				seeds2[1].push_back(MQ1 | MQ2);
				seeds2[1].push_back(MQ1 | MQ2);
				seeds2[1].push_back(MQ2 | MASK_quarter);

				seeds2[2].push_back(MQ1 | MASK_quarter);
				seeds2[2].push_back(MQ1 | MQ2);
				seeds2[2].push_back(MQ1 | MASK_quarter);

				seeds2[3].push_back(MQ1 | MASK_quarter);
				seeds2[3].push_back(MQ2 | MASK_quarter);



			len_mer = new_len_mer;
		}//mism > 0
	}//pairs
	else{
		if(mism_thresh == 0){
			gen_seeds.push_back(MASK_LMER);
			vector<unsigned long int> dum;
			dum.push_back(MASK_LMER);

			seeds.push_back(dum);
		}
		else{
			new_len_mer = double_quarter + quarter;

			int shift = len_mer -  new_len_mer;

			MASK_LMER = MASK_LMER >> shift;
			len_mer = new_len_mer;

			gen_seeds.resize(4);
			gen_seeds[0] = MASK_LMER;
			gen_seeds[1] = MQ1 | MQ2;
			gen_seeds[2] = MQ2 | MASK_quarter;
			gen_seeds[3] = MQ1 | MASK_quarter;


			seeds.resize(4);

			seeds[0].push_back(MASK_LMER);
			seeds[1].push_back(MQ1 | MQ2);
			seeds[2].push_back(MQ2 | MASK_quarter);
			seeds[3].push_back(MQ1 | MASK_quarter);
		}
	
	}//singles

}//make seeds


void Anchor::find_pos_sizes(int t)
{

	long int i, u;
	string line;
	ifstream in;

	for(i = 0; i < hash_table_size; i++)
		pos_sizes[i] = 0;

	//find out sizes of chromosomes (or references)
	const long int buffer_size = 10000000;
	char *buffer = new char[buffer_size];
	assert( buffer != 0);
	long int lmer_pos = 0;
	long int real_size = 0;

	for(u = 0; u < num_chroms; u++){

		long int cur_size = 0;
		long int chrom_id = u;

		lmer_pos = starts_size_chrom[chrom_id];

		in.open(names[u].c_str(), ios::in);
		if(!in){
			cout << "\nERROR: cannot open " << names[u] << endl;
			exit(0);
		}

		getline(in, line);//first line is not counted has >
		if(line[0] != '>'){
		
			cout << "\nERROR: " << names[u] << " is not in FASTA format. " << endl;
			exit(0);
		}

		in.read(buffer, buffer_size);
		real_size = in.gcount();

		if(real_size == 0){
			cout << "\nERROR: no input in " << names[u] << endl;
			exit(0);
		}
		unsigned long int lmer_1 = len_mer - 1;
		unsigned long int ta_seed = 0, ns_seed = 0, ta_ind2, ta_ind, ta_lmer=0, cg_lmer=0;
		long int last_lmer_char;

		string temp_seed("");
		last_lmer_char = 0;
		i = 0;
		while(i < lmer_1){
			if(buffer[last_lmer_char] != '\n'){
				temp_seed = temp_seed + buffer[last_lmer_char];
				i++;
			}
			last_lmer_char++;
		}//for i

		convert_lmer_ta(temp_seed, ta_seed, ns_seed, lmer_1);
		unsigned long int cg_seed = 0;
		convert_lmer_c(temp_seed, cg_seed,  lmer_1);
		long int ends_of_lines = 0;
		long int stop ;
		while(true){

			cur_size += real_size;

			stop = real_size - lmer_1;
			for(; last_lmer_char < real_size; last_lmer_char++){

				char ch_seed = buffer[last_lmer_char];
				if(ch_seed != '\n'){

					next_lmer(ch_seed, ta_seed, ns_seed, cg_seed, MASK_LMER);

	if(ns_seed == 0){
				

						ta_lmer = ta_seed & gen_seeds[t];
						cg_lmer = cg_seed & gen_seeds[t];

					if(is_bs == "1"){

						if(mism_thresh == 0){

							ta_ind = hash2gs(ta_lmer, cg_lmer);
							ta_ind2 = hash2cs(ta_lmer, cg_lmer);

							pos_sizes[ta_ind]++;
							if(ta_ind != ta_ind2)
								pos_sizes[ta_ind2]++;
							if(is_pe == "1" && gen_seeds[t] != gen_seeds2[t]){
								ta_lmer = ta_seed & gen_seeds2[t];
								cg_lmer = cg_seed & gen_seeds2[t];
								ta_ind = hash2gs(ta_lmer, cg_lmer);
								ta_ind2 = hash2cs(ta_lmer, cg_lmer);
								pos_sizes[ta_ind]++;
								if(ta_ind != ta_ind2)
									pos_sizes[ta_ind2]++;
							}
						}//if thresh_mism = 0
						else{//thresh_mism > 0
							ta_ind = hash2m(ta_lmer);
							pos_sizes[ta_ind]++;
							if(is_pe == "1" && gen_seeds[t] != gen_seeds2[t]){
								ta_lmer = ta_seed & gen_seeds2[t];
								ta_ind2 = hash2m(ta_lmer);
								if(ta_ind != ta_ind2)
									pos_sizes[ta_ind2]++;
							}
						
						}//thresh_mism > 0

						}
						else{
							ta_ind = hash2c(ta_lmer, cg_lmer);
							pos_sizes[ta_ind]++;

							if(is_pe == "1" && gen_seeds[t] != gen_seeds2[t]){
								ta_lmer = ta_seed & gen_seeds2[t];
								cg_lmer = cg_seed & gen_seeds2[t];
								ta_ind2 = hash2c(ta_lmer, cg_lmer);
								if(ta_ind != ta_ind2)
									pos_sizes[ta_ind2]++;
							}

						}//else
	}//if no Nx
						lmer_pos++;
				}//if char is not end of line
				else{
					ends_of_lines++;
				}
//				last_lmer_char++;

			}//for i

			last_lmer_char = 0;

			if(in.eof())
				break;

			in.read(buffer, buffer_size);
			real_size = in.gcount();

			if(real_size == 0)
				break;

		}//while
		in.close(); in.clear();
		
		cur_size -= ends_of_lines;

		i = u+1;
	}//for u

	for(i = 0; i < hash_table_size; i++){
		long int cur_size = pos_sizes[i];
		pos_index[i] = new unsigned int[cur_size];
		assert( pos_index[i] != 0);

		pos_sizes[i] = 0; 
	}//for i


	delete [] buffer;
//////////////////////////
	cout << "Space for hash table is allocated" << endl;
}//find pos sizes

void Anchor::make_pos_index(int t)
{
	unsigned long int i, ch, j, u, y;


	//out2 is a file handler for Pos Index
	unsigned long int rows = hash_table_size;
	unsigned int lmer_pos =0;
	string aname;
	
	string next, line, chrom;
	unsigned short int chrom_id = 0;

	//the first chromosome has starts 0/0
	unsigned long int cur_chrom_cg_size = 0;
	unsigned long int cur_chrom_ta_size = 0;
	unsigned long int count_bytes = 0;
	unsigned long int count_bytes_ta = 0;
	unsigned long int count_blocks_ta = 0;
		unsigned long int count_blocks = 0;

	unsigned long int one = 1;
	unsigned long int double_byte = 16;
	unsigned long int quadra_byte = 32;
	unsigned long int triple_byte = 24;
	unsigned long int x ;
	long int index;
	ifstream in;
	long int last_lmer_char;
	const long int buffer_size = 10000000;
	char *buffer = new char[buffer_size];
	assert( buffer != 0);

	long int real_size = 0;

	for(u = 0; u < names.size(); u++){ 


		chrom_id = u;
		lmer_pos = starts_size_chrom[chrom_id];
		index = 0;

		in.open(names[u].c_str(), ios::in);
		if(!in){
			cout << "can't open " << names[u] << " names[u] reading genome" << endl;
			exit(0);
		}

		getline(in, line); //the first line has >


		in.read(buffer, buffer_size);
		real_size = in.gcount();
		if(real_size == 0){
			cout << "\nERROR: " << names[u] << " has no input " << endl;
			exit(0);
		}
		unsigned long int cmer = 0;
		unsigned long int gmer = 0;
		unsigned long int ta_seed = 0;
		
		unsigned long int count_bytes = 0;
		unsigned long int count_bits = 0;
		unsigned long int count_bytes_ta = 0;
		unsigned long int count_blocks_ta = 0;
		unsigned long int count_blocks = 0;

		unsigned long int lmer_1 = len_mer - 1;
		unsigned long int gta = 0;
		unsigned long int ns = 0, ns_seed = 0, cg_seed = 0;
		unsigned long int cg = 0;
		unsigned long int ta_ind = 0, ta_ind2 = 0;
		unsigned long int ta = 0;
		long int stop;
		string temp_seed("");
		last_lmer_char = 0;
		i = 0;
		while(i < lmer_1){
			if(buffer[last_lmer_char] != '\n'){
				temp_seed = temp_seed + buffer[last_lmer_char];
				i++;
			}
			last_lmer_char++;
		}//for i

			for(i = 0; i < lmer_1; i++){
	
				char next_char = temp_seed[i];

					cg <<= 1;

					ta <<= 1;

					ns <<= 1;
					if(next_char == 'T' || next_char == 't'){
							ta |= 1;
					}
					else if(next_char == 'C' || next_char == 'c' ){
						cg |= 1;
						ta |= 1;//for c
					}
					else if(next_char == 'G' || next_char == 'g')
						cg |=1;
					else if(next_char == 'N' || next_char == 'n')
						ns |= 1;

					count_bits++;

					ns_seed = ns & MASK_LMER;
					ta_seed = ta & MASK_LMER;
					cg_seed = cg & MASK_LMER;


			}//for i



		unsigned int chroms_bits = CHROMS_BITS;


		while(true){

			stop = real_size - lmer_1;
			for(; last_lmer_char < real_size; last_lmer_char++){

				char next_char = buffer[last_lmer_char];

				if(next_char != '\n'){

					cg <<= 1;
					ta <<= 1;
					ns <<= 1;
					if(next_char == 'T' || next_char == 't')
							ta |= 1;
					else if(next_char == 'C' || next_char == 'c' ){
						cg |= 1;
						ta |= 1;//for c
					}
					else if(next_char == 'G' || next_char == 'g')
						cg |=1;
					else if(next_char == 'N' || next_char == 'n' )
						ns |=1;

					ns_seed = ns & MASK_LMER;

					count_bits++;

					if(ns_seed == 0){

						ta_seed = ta & gen_seeds[t];
						cg_seed = cg & gen_seeds[t];

						unsigned int gi = lmer_pos ;

					if(is_bs == "1"){

						if(mism_thresh == 0){
							ta_ind = hash2gs(ta_seed, cg_seed);
							ta_ind2 = hash2cs(ta_seed, cg_seed);
							long int cur_pointer = pos_sizes[ta_ind];

							pos_index[ta_ind][cur_pointer] = gi;//y;
							pos_sizes[ta_ind] = cur_pointer + 1;

							if(ta_ind != ta_ind2){
								cur_pointer =  pos_sizes[ta_ind2];
								pos_index[ta_ind2][cur_pointer] = gi;//y;
								pos_sizes[ta_ind2] = cur_pointer + 1;
								
							}//if

							if(is_pe == "1" && gen_seeds[t] != gen_seeds2[t]){
								ta_seed = ta & gen_seeds2[t];
								cg_seed = cg & gen_seeds2[t];
								ta_ind = hash2gs(ta_seed, cg_seed);
								ta_ind2 = hash2cs(ta_seed, cg_seed);
								cur_pointer = pos_sizes[ta_ind];

								pos_index[ta_ind][cur_pointer] = gi;//y;
								pos_sizes[ta_ind] = cur_pointer + 1;

								if(ta_ind != ta_ind2){
									cur_pointer =  pos_sizes[ta_ind2];
									pos_index[ta_ind2][cur_pointer] = gi;//y;
									pos_sizes[ta_ind2] = cur_pointer + 1;
									
								}//if
							}//if pair-end
						}//if thresh-mism = 0
						else{
							ta_ind = hash2m(ta_seed);
							long int cur_pointer = pos_sizes[ta_ind];

							pos_index[ta_ind][cur_pointer] = gi;
							pos_sizes[ta_ind]++;

							if(is_pe == "1" && gen_seeds[t] != gen_seeds2[t]){
								ta_seed = ta & gen_seeds2[t];
								ta_ind2 = hash2m(ta_seed);
								if(ta_ind != ta_ind2){
									cur_pointer = pos_sizes[ta_ind2];
									pos_index[ta_ind2][cur_pointer] = gi;
									pos_sizes[ta_ind2]++;
								}
							}//if pe
						}//thresh-mism > 0


						}//if BS
						else{ //normal mapping 
							ta_ind = hash2c(ta_seed, cg_seed);
							long int cur_pointer = pos_sizes[ta_ind];
							pos_index[ta_ind][cur_pointer] = gi;//y;
							pos_sizes[ta_ind] = cur_pointer + 1;
							if(is_pe == "1"  && gen_seeds[t] != gen_seeds2[t]){
								ta_seed = ta & gen_seeds2[t];
								cg_seed = cg & gen_seeds2[t];
								ta_ind2 = hash2c(ta_seed, cg_seed);
								if(ta_ind != ta_ind2){
									cur_pointer = pos_sizes[ta_ind2];
									pos_index[ta_ind2][cur_pointer] = gi;//y;
									pos_sizes[ta_ind2] = cur_pointer + 1;
								}
							}//if pair-end

						}//else
					}//if No Ns
					lmer_pos++;

				}//if not end of line char
	//			last_lmer_char++;
			}//for i

			last_lmer_char = 0;

			if(in.eof())
				break;

			in.read(buffer, buffer_size);
			real_size = in.gcount();

			if(real_size == 0)
				break;

		}//while
		 in.close(); in.clear();
		int t1 = t+1;
		cout << "round " << t1 << ": " << names[u] << " is hashed " << endl;

	}//for u all chroms

	delete [] buffer;

}//make_pos_index

//read genome
void Anchor::read_genome()
{
	//out2 is a file handler for Pos Index
	unsigned long int i, ch, j, u, y;
	string aname;
	
	unsigned long int Mask_All = 0xFFFFFFFFFFFFFFFF;

	string next, line, chrom;
	unsigned int chrom_id = 0;

	//the first chromosome has starts 0/0
	unsigned long int cur_chrom_cg_size = 0;
	unsigned long int cur_chrom_ta_size = 0;
	unsigned long int count_bytes = 0;
	unsigned long int count_bytes_ta = 0;
	unsigned long int count_blocks_ta = 0;
		unsigned long int count_blocks = 0;

	unsigned long int one = 1;
	unsigned long int double_byte = 16;
	unsigned long int quadra_byte = 32;
	unsigned long int triple_byte = 24;
	unsigned long int x ;
	long int index;
	ifstream in;
	long int last_lmer_char;
	const long int buffer_size = 10000000;
	char *buffer = new char[buffer_size];
	assert( buffer != 0);

	long int real_size = 0;
	long int lmer_pos ;

	for(u = 0; u < names.size(); u++){ 

		long int cur_size = 0;

		chrom_id = u;
		index = starts_chroms_starts[chrom_id];
		lmer_pos = 0;
		in.open(names[u].c_str(), ios::in);
		if(!in){
			cout << "can't open " << names[u] << " names[u] reading genome" << endl;
			exit(0);
		}

		getline(in, line); //the first line has >

		line = line.substr(1, line.length());
		istringstream istr(line);
		istr >> aname;
		chroms_names.insert(make_pair(chrom_id, aname));

		in.read(buffer, buffer_size);
		real_size = in.gcount();

		if(real_size == 0){
			cout << "\nERROR: " << names[u] << " has no input " << endl;
			exit(0);
		}
		unsigned long int cmer = 0;
		unsigned long int gmer = 0;
		unsigned long int ta_seed = 0;
		
		unsigned long int count_bytes = 0;
		unsigned long int count_bits = 0;
		unsigned long int count_bytes_ta = 0;
		unsigned long int count_blocks_ta = 0;
		unsigned long int count_blocks = 0;

		unsigned long int lmer_1 = len_mer - 1;
		unsigned long int gta = 0;
		unsigned long int ns = 0, ns_seed = 0, cg_seed = 0;
		unsigned long int cg = 0;
		unsigned long int ta_ind = 0, ta_ind2 = 0;
		unsigned long int ta = 0;
		long int stop;
		string temp_seed("");
		last_lmer_char = 0;
		i = 0;
		while(i < lmer_1){
			if(buffer[last_lmer_char] != '\n'){
				temp_seed = temp_seed + buffer[last_lmer_char];
				i++;
			}
			last_lmer_char++;
		}//for i

			for(i = 0; i < lmer_1; i++){
	
				char next_char = temp_seed[i];

					cg <<= 1;

					ta <<= 1;
					ns <<= 1;

					if(next_char == 'T' || next_char == 't'){
							ta |= 1;
					}
					else if(next_char == 'C' || next_char == 'c' ){
						cg |= 1;
						ta |= 1;//for c
					}
					else if(next_char == 'G' || next_char == 'g')
						cg |=1;
					else if(next_char == 'N' || next_char == 'n')
						ns |= 1;

					count_bits++;

					ns_seed = ns & MASK_LMER;
					ta_seed = ta & MASK_LMER;
					cg_seed = cg & MASK_LMER;

					if(count_bits % SIZE == 0){//BUG possible: change byte_size = 4

						ta_genome[index] =ta ;
						cg_genome[index] = cg ;
						ns_genome[index] = ns;
						index++;

					}//if

			}//for i



		unsigned int chroms_bits = CHROMS_BITS;


		while(true){

			stop = real_size - lmer_1;
			for(; last_lmer_char < real_size; last_lmer_char++){

				char next_char = buffer[last_lmer_char];

				if(next_char != '\n'){

					cg <<= 1;
					ta <<= 1;
					ns <<= 1;

					if(next_char == 'T' || next_char == 't')
							ta |= 1;
					else if(next_char == 'C' || next_char == 'c' ){
						cg |= 1;
						ta |= 1;//for c
					}
					else if(next_char == 'G' || next_char == 'g')
						cg |=1;
					else if(next_char == 'N' || next_char == 'n')
						ns |= 1;

					count_bits++;

					if(count_bits % SIZE == 0){//BUG possible: change byte_size = 4

						ta_genome[index] = ta;
						cg_genome[index] = cg;
						ns_genome[index] = ns;

						index++;
						lmer_pos++;
					}//if
				}//if not end of line char
			}//for i

			last_lmer_char = 0;

			if(in.eof())
				break;

			in.read(buffer, buffer_size);
			real_size = in.gcount();

			if(real_size == 0)
				break;


		}//while
		 in.close(); in.clear();



		int remainder_cg = count_bits % SIZE;//same for ta
		
		unsigned long int Mask_remainder = 0xFFFFFFFFFFFFFFFF;
		if(remainder_cg > 0){
			Mask_remainder = (one << (SIZE - remainder_cg)) - 1;
			cg = (cg << (SIZE - remainder_cg));
			ta = (ta << (SIZE - remainder_cg)) ;
			ns = (ns << (SIZE - remainder_cg)) | Mask_remainder;
	


			ta_genome[index] =ta;
			cg_genome[index] = cg;
			ns_genome[index] = ns;
			index++;
		}
		//starts_chroms_starts contains starts for chroms
		for(; index < starts_chroms_starts[chrom_id + 1]; index++)
			ns_genome[index] = Mask_All;

		cout << names[u] << " is preprocessed " << endl;

	}//for u all chroms

	delete [] buffer;

}//
//read genome
////////////////////////// LONG READS

void Anchor::convert_gen_mer_ta_long(vector<unsigned long int> &ta_gen_lmer,
					 unsigned long int pos,
					 unsigned long int aread_len, int asize)
{

	unsigned long int cur_byte, first, second;
	//current binary char inside current block: there 8 characters fit in one binary characters
	//(pos >> two_in_byte) = pos / byte_size = number of binary characters needed to cover this pos
	cur_byte = (pos >> two_in_Size) ;

	unsigned long int char_pos = pos & MASK_SIZE;//pos of genome char inside byte-block
	long int i, j;
	//this is equivalent to bytes_in_lmer = ceil(aread_len/byte_size)

	int len_rest = aread_len & MASK_SIZE;

	if(char_pos == 0){
		j = cur_byte;
		for(i = 0; i < asize  && j < binary_total; i++){
			ta_gen_lmer.push_back(ta_genome[j]);

			j++;
		}//for i

		if(len_rest < SIZE)
			ta_gen_lmer[asize - 1] >>= (SIZE - len_rest);

	}
	else{

		//This is not done
		for(i = 0; i < asize && cur_byte < binary_total; i++){
			first = ta_genome[cur_byte];
			second = ta_genome[cur_byte + 1];
			
			first <<= char_pos;
			second >>= (SIZE - char_pos);
			ta_gen_lmer.push_back(first | second);
			cur_byte++;

		}
		if(len_rest < SIZE)
			ta_gen_lmer[asize - 1] >>= (SIZE - len_rest);
	}//else


}//convert_gen_mer_ta

void Anchor::convert_gen_mer_cg_long(vector<unsigned long int> &ta_gen_lmer,
					 unsigned long int pos,
					 unsigned long int aread_len, int asize)
{

	unsigned long int cur_byte, first, second;
	//current binary char inside current block: there 8 characters fit in one binary characters
	//(pos >> two_in_byte) = pos / byte_size = number of binary characters needed to cover this pos
	cur_byte = (pos >> two_in_Size) ;

	unsigned long int char_pos = pos & MASK_SIZE;//pos of genome char inside byte-block
	long int i, j;
	//this is equivalent to bytes_in_lmer = ceil(aread_len/byte_size)

	int len_rest = aread_len & MASK_SIZE;

	if(char_pos == 0){
		j = cur_byte;
		for(i = 0; i < asize - 1 && j < binary_total; i++){
			ta_gen_lmer.push_back(cg_genome[j]);
			j++;
		}//for i
		if(len_rest < SIZE)
			ta_gen_lmer.push_back(cg_genome[j] >> (SIZE - len_rest));
	}
	else{

		for(i = 0; i < asize && cur_byte < binary_total; i++){
			first = cg_genome[cur_byte];
			second = cg_genome[cur_byte + 1];
			
			first <<= char_pos;
			second >>= (SIZE - char_pos);
			ta_gen_lmer.push_back(first | second);
			cur_byte++;
		}
		if(len_rest < SIZE)
			ta_gen_lmer[asize - 1] >>= (SIZE - len_rest);
	}//else


}//convert_gen_mer_ta

void Anchor::convert_gen_mer_ns_long(vector<unsigned long int> &ta_gen_lmer,
					 unsigned long int pos,
					 unsigned long int aread_len, int asize)
{

	unsigned long int cur_byte, first, second;
	//current binary char inside current block: there 8 characters fit in one binary characters
	//(pos >> two_in_byte) = pos / byte_size = number of binary characters needed to cover this pos
	cur_byte = (pos >> two_in_Size) ;

	unsigned long int char_pos = pos & MASK_SIZE;//pos of genome char inside byte-block
	long int i, j;
	//this is equivalent to bytes_in_lmer = ceil(aread_len/byte_size)

	int len_rest = aread_len & MASK_SIZE;

	if(char_pos == 0){
		j = cur_byte;
		for(i = 0; i < asize  && j < binary_total; i++){
			ta_gen_lmer.push_back(ns_genome[j]);

			j++;
		}//for i
		if(len_rest < SIZE)
			ta_gen_lmer[asize - 1] >>= (SIZE - len_rest);

	}
	else{


		for(i = 0; i < asize && cur_byte < binary_total; i++){
			first = ns_genome[cur_byte];
			second = ns_genome[cur_byte + 1];
			
			first <<= char_pos;
			second >>= (SIZE - char_pos);
			ta_gen_lmer.push_back(first | second);
			cur_byte++;

		}
		if(len_rest < SIZE)
			ta_gen_lmer[asize - 1] >>= (SIZE - len_rest);
	}//else


}//convert_gen_mer_ta

///////////////////////////// END LONG

void Anchor::convert_gen_mer_ta(unsigned long int &ta_gen_lmer,
					 unsigned long int pos,
					 unsigned long int aread_len)
{

	unsigned long int cur_byte, first,  second;
	//current binary char inside current block: there 8 characters fit in one binary characters
	//(pos >> two_in_byte) = pos / byte_size = number of binary characters needed to cover this pos
	cur_byte = (pos >> two_in_Size) ;

	unsigned long int char_pos = pos & MASK_SIZE;//pos of genome char inside byte-block

	//this is equivalent to bytes_in_lmer = ceil(aread_len/byte_size)
	first = ta_genome[cur_byte];

	if(char_pos == 0)
		ta_gen_lmer = first >> (SIZE - aread_len);
	else{
		cur_byte++;
		second = ta_genome[cur_byte];
		
		first <<= char_pos;
		second >>= (SIZE - char_pos);
		ta_gen_lmer = (first | second) >> (SIZE - aread_len);
	}


}//convert_gen_mer_ta


void Anchor::convert_gen_mer_cg(unsigned long int &ta_gen_lmer,
					 unsigned long int pos,
					 unsigned long int aread_len)
{
	unsigned long int cur_byte, first,  second;
	//current binary char inside current block: there 8 characters fit in one binary characters
	//(pos >> two_in_byte) = pos / byte_size = number of binary characters needed to cover this pos
	cur_byte = (pos >> two_in_Size) ;

	unsigned long int char_pos = pos & MASK_SIZE;//pos of genome char inside byte-block

	//this is equivalent to bytes_in_lmer = ceil(aread_len/byte_size)
	first = cg_genome[cur_byte];
	if(char_pos == 0)
		ta_gen_lmer = first >> (SIZE - aread_len);
	else{
		cur_byte++;
		second = cg_genome[cur_byte];
		
		first <<= char_pos;
		second >>= (SIZE - char_pos);
		ta_gen_lmer = (first | second) >> (SIZE - aread_len);
	}


}//convert_gen_mer_cg

void Anchor::convert_gen_mer_ns(unsigned long int &ta_gen_lmer,
					 unsigned long int pos,
					 unsigned long int aread_len)
{

	unsigned long int cur_byte, first,  second;
	//current binary char inside current block: there 8 characters fit in one binary characters
	//(pos >> two_in_byte) = pos / byte_size = number of binary characters needed to cover this pos
	cur_byte = (pos >> two_in_Size) ;

	unsigned long int char_pos = pos & MASK_SIZE;//pos of genome char inside byte-block

	//this is equivalent to bytes_in_lmer = ceil(aread_len/byte_size)
	first = ns_genome[cur_byte];

	if(char_pos == 0)
		ta_gen_lmer = first >> (SIZE - aread_len);
	else{
		cur_byte++;
		second = ns_genome[cur_byte];
		
		first <<= char_pos;
		second >>= (SIZE - char_pos);
		ta_gen_lmer = (first | second) >> (SIZE - aread_len);
	}


}//convert_gen_mer_ta


void Anchor::pair_reads(unsigned long int read_id, string aread1, string aread2,
				long int &contain_Ns, Stat &astat, short int read_len1, short int read_len2, int t, int rs )
{


	unsigned long int one = 1;

	long int asize,  pos1, pos2, dist, i;
	string strand, strand2;
	unsigned short int mism1, mism2;
	unsigned long int ta_lmer1, ta_lmer2, cg_lmer1, cg_lmer2;
	string line;

	unsigned long int MASK_READ1 = ( 1 << read_len1) - 1;
	unsigned long int MASK_READ2 = (1 << read_len2) -1;


	if(read_len1 > 31){
		dist = read_len1 - 32;
		MASK_READ1 = (((one << dist) - 1) << 32) | 0x00000000FFFFFFFF;
	}
	if(read_len2 > 31){
		dist = read_len2 - 32;
		MASK_READ2 = (((one << dist) - 1) << 32) | 0x00000000FFFFFFFF;
		
	}

	unsigned long int ta_rc1=0, ta_rc2=0, ind1 = 0, ind2 = 0;//hash-indices 
	unsigned long int c_rc_read1, c_rc_read2, c_read1, c_read2;
	unsigned long int  ta_read1, ta_read2, ta_rc_read1, ta_rc_read2, ns=0, ns2=0;
	unsigned long int ta_gen_mer1 = 0, cg_gen_mer1 = 0, ta_gen_mer2 =0, cg_gen_mer2 = 0;

	bool pair_ambiguous = false;

	unsigned long int tc_mism1, ag_mism1, tc_mism2, ag_mism2;

	//ta-reads
	convert_lmer_ta(aread1, ta_read1, ns, read_len1);
	convert_lmer_ta(aread2, ta_read2, ns2, read_len2);
	
	//c-reads
	convert_lmer_c(aread1, c_read1,  read_len1);
	convert_lmer_c(aread2, c_read2,  read_len2);

  //rev-complement of c-reads

	int dif1 = read_len1 - len_mer;
	int dif2 = read_len2 - len_mer;
	ta_lmer1 = (ta_read1 >> dif1) & seeds[t][rs];
	cg_lmer1 = (c_read1 >> dif1) & seeds[t][rs];

	cg_lmer2 = (c_read2 >> dif2) & seeds2[t][rs];
	ta_lmer2 = (ta_read2 >> dif2) & seeds2[t][rs];

	ind1 = hash_aread(ta_lmer1, cg_lmer1, false);// hash2gs(ta_lmer1,  cg_lmer1, MASK_SEED, MASK_LMER);

	ind2 = hash_aread(ta_lmer2, cg_lmer2, true);




	//rev-complement of g-reads


	ta_rc_read2 = MASK_READ2 & (~Reverse(ta_read2, read_len2));
	c_rc_read2 =  Reverse(c_read2, read_len2);//shows positions of Gs in reverse strand

	ta_lmer2 = (ta_rc_read2 & MASK_LMER) & seeds2[t][rs];
	cg_lmer2 = (c_rc_read2 & MASK_LMER) & seeds2[t][rs];
	ta_rc2 =  hash_aread(ta_lmer2, cg_lmer2, false);



	//process forward read1 and rev read2
	//start and end of list of genome positions for forw1
	long int start1 = 0;
	long int end1 = pos_sizes[ind1];

	//start and end of list of genome positions for rev-compl2

	long int start2 = 0;
	long int end2 = pos_sizes[ta_rc2];
	long int cur = start2;

	while(start1 < end1){
		if(start2 >= end2)
			break;

		pos1 = pos_index[ind1][start1];


			cur = start2;
			for(; cur != end2 ; cur++){
				pos2 = pos_index[ta_rc2][cur];
				if(dif2 <= pos2){
					pos2 -= dif2;
					if(A_rich == true)
						dist = pos1 - pos2;//to map reads from PCR1- and PCR2+
					else
						dist = pos2 - pos1;
					if(dist >= min_insert ){
						break;
					}
				}//if dif2 <= pos2
			}//for cur

			start2 = cur;


			for(; cur != end2; cur++){
				pos2 = pos_index[ta_rc2][cur];
				if(dif2 <= pos2){
					pos2 -= dif2;

				if(A_rich == true)
					dist = pos1 - pos2;
				else
					dist = pos2 - pos1;
				if(dist >= min_insert && dist <= max_insert ){
					//cg check
					bool cg_mism1 = false, cg_mism2 = false;
					convert_gen_mer_ta(ta_gen_mer1, pos1,  read_len1);

					unsigned long int ta_check = ta_read1 ^ ta_gen_mer1, ta_check2;
					mism1 = count_mism(ta_check);

					if(mism1 <= mism_thresh){

						convert_gen_mer_ta(ta_gen_mer2, pos2, read_len2);
						ta_check2 = ta_rc_read2 ^ ta_gen_mer2;
						mism2 = count_mism(ta_check2);

						if(mism2 <= mism_thresh){
							unsigned long int ns_gen_mer1 = 0, ns_gen_mer2 =  0;
							convert_gen_mer_cg(cg_gen_mer1, pos1,  read_len1);
							convert_gen_mer_cg(cg_gen_mer2, pos2,  read_len2);
							convert_gen_mer_ns(ns_gen_mer1, pos1,  read_len1);
							convert_gen_mer_ns(ns_gen_mer2, pos2,  read_len2);

							ag_mism1 = (~(ta_read1 | c_read1)) & cg_gen_mer1;
							ta_check |= (ag_mism1 | ns_gen_mer1);
							ag_mism2 = (~(ta_rc_read2 | c_rc_read2)) & cg_gen_mer2;
							ta_check2 |= (ag_mism2 | ns_gen_mer2);

							unsigned long int cg_check = c_read1 ^ (c_read1 & cg_gen_mer1);

							ta_check |=  cg_check;
							unsigned long int cg_check2 = c_rc_read2 ^ (c_rc_read2 & cg_gen_mer2);
							ta_check2 |= cg_check2 ;

								mism1 = count_mism(ta_check);
								if(mism1 > mism_thresh)
									cg_mism1 = true;

								mism2 = count_mism(ta_check2);
								if(mism2 > mism_thresh)
									cg_mism2 = true;



						}//if passed second ta-check
						else
							cg_mism2 = true;//not passed 2 ta-check
					}//if passed ta=check
					else{
						cg_mism1 = true;
					}//else



					if(cg_mism1 == false && cg_mism2 == false){
						//T->C and A->G not in the same read check
						//read1 check
						char astrand = '+';
						//check tc/ag mism don't happen in mates
							int min_mism = mism1;
							int second_mism = mism2;
							if(mism2 < mism1){
								min_mism = mism2;
								second_mism = mism1;
							}
							int cur_second = mates_pairs[read_id].min_mism;
							if(cur_second <  mates_pairs[read_id].mism1)
								cur_second = mates_pairs[read_id].mism1;
							else if(cur_second < mates_pairs[read_id].mism2)
								cur_second = mates_pairs[read_id].mism2;

							if(mates_pairs[read_id].mark_unique == 0){
								mates_pairs[read_id].mark_unique = 1;
								mates_pairs[read_id].pos1 = pos1;
								mates_pairs[read_id].mism1 = mism1;
								mates_pairs[read_id].pos2 = pos2;
								mates_pairs[read_id].mism2 = mism2;
								mates_pairs[read_id].min_mism = min_mism;
								mates_pairs[read_id].strand = '+';//by the first read
							}//if
							else if(mates_pairs[read_id].pos1 != pos1 || 
								mates_pairs[read_id].pos2 != pos2 || mates_pairs[read_id].strand != astrand){


								if((min_mism < mates_pairs[read_id].min_mism) || 
									(min_mism == mates_pairs[read_id].min_mism && 
									second_mism < cur_second)){
									//unique
		
									mates_pairs[read_id].mark_unique = 1;
									mates_pairs[read_id].pos1 = pos1;
									mates_pairs[read_id].strand = '+';//by the first read
									mates_pairs[read_id].mism1 = mism1;
									mates_pairs[read_id].pos2 = pos2;
									mates_pairs[read_id].mism2 = mism2;
									mates_pairs[read_id].min_mism = min_mism;
										
								}
								else if(min_mism == mates_pairs[read_id].min_mism &&
									second_mism == cur_second){

									mates_pairs[read_id].mark_unique++;
									if(min_mism == 0 && second_mism == 0)
										pair_ambiguous = true;
								}
							}//else if

					}//both reads passed cg-check
				}//if correct insert range

				else if(dist > max_insert)
					break;
				}//if dif2 <= pos2
				if(pair_ambiguous == true){					
					break;
				}
			}//while cur < end2 && cur <= max_insert
			
			start1++;

		if(pair_ambiguous == true){
			break;
		}//if
	}//for forward read's list

	if(pair_ambiguous == false){
		//take rev-complement of the first read
		ta_rc_read1 = MASK_READ1 & (~Reverse(ta_read1, read_len1));


		c_rc_read1 =  Reverse(c_read1, read_len1);//shows positions of Gs in reverse strand
		ta_lmer1 = (ta_rc_read1 & MASK_LMER) & seeds[t][rs];
		cg_lmer1 = (c_rc_read1  & MASK_LMER) & seeds[t][rs];
		ta_rc1 = hash_aread(ta_lmer1, cg_lmer1, true);

		//process rev-complemnet of read1 and forward of read2

		start1 = 0;
		end1 = pos_sizes[ind2];

		start2 = 0;
		end2 = pos_sizes[ta_rc1];

		while(start1 < end1){

			if(start2 >= end2)
				break;

			pos1 =  pos_index[ind2][start1];

				cur = start2;
				for(; cur != end2 ; cur++){
					pos2 =  pos_index[ta_rc1][cur] ;
					if(	dif1 <= pos2){		
						pos2 -= dif1;
						if(A_rich == true)
							dist = pos1 - pos2;
						else
							dist = pos2 - pos1;
						
						if(dist >= min_insert ){
							break;
						}//if
					}//if
				}//for cur

				start2 = cur;

					for(; cur != end2; cur++){
						pos2 =  pos_index[ta_rc1][cur] ;

						if(dif1 <= pos2){
						pos2 -= dif1;
						if(A_rich == true)
							dist = pos1 - pos2;
						else
							dist = pos2 - pos1;
						if(dist >= min_insert && dist <= max_insert ){
					//cg check
					
							bool cg_mism1 = false, cg_mism2 = false;
							convert_gen_mer_ta(ta_gen_mer2, pos1,  read_len2);


					unsigned long int ta_check2 = ta_read2 ^ ta_gen_mer2, ta_check;
					mism2 = count_mism(ta_check2);

					if(mism2 <= mism_thresh){
						convert_gen_mer_ta(ta_gen_mer1, pos2, read_len1);
						ta_check = ta_rc_read1 ^ ta_gen_mer1;
						mism1 = count_mism(ta_check);

						if(mism1 <= mism_thresh){
							unsigned long int ns_gen_mer1 = 0, ns_gen_mer2 = 0;
							convert_gen_mer_cg( cg_gen_mer2, pos1,  read_len2);
							convert_gen_mer_cg( cg_gen_mer1, pos2,  read_len1);
							convert_gen_mer_ns(ns_gen_mer1, pos2, read_len1);
							convert_gen_mer_ns(ns_gen_mer2, pos1, read_len2);

							tc_mism2 = (ta_read2 & (~c_read2)) & cg_gen_mer2;
							tc_mism1 = (ta_rc_read1 &(~c_rc_read1)) & cg_gen_mer1;

							ta_check2 |= (tc_mism2 | ns_gen_mer2);
							ta_check |= (tc_mism1 | ns_gen_mer1);

								unsigned long int cg_check2 = c_read2 ^ (c_read2 & cg_gen_mer2);
								ta_check2 |= cg_check2 ;
								unsigned long int cg_check =c_rc_read1 ^ (c_rc_read1 & cg_gen_mer1);
								ta_check |=  cg_check;

								
							mism2 = count_mism(ta_check2);	
							if(mism2 > mism_thresh)
								cg_mism2 = true;

							mism1 = count_mism(ta_check);
							if(mism1 > mism_thresh)
								cg_mism1 = true;



						}//if second ta-check passed
						else
							cg_mism2 = true;
						}//if first ta_check passed
						else{
							cg_mism1 = true;
						}//else


							if(cg_mism1 == false && cg_mism2 == false){

							char astrand = '-';
							int min_mism = mism1;
							int second_mism = mism2;
							if(mism2 < mism1){
								min_mism = mism2;
								second_mism = mism1;
							}
							int cur_second = mates_pairs[read_id].min_mism;
							if(cur_second <  mates_pairs[read_id].mism1)
								cur_second = mates_pairs[read_id].mism1;
							else if(cur_second < mates_pairs[read_id].mism2)
								cur_second = mates_pairs[read_id].mism2;

							if(mates_pairs[read_id].mark_unique == 0){
								mates_pairs[read_id].mark_unique = 1;
								mates_pairs[read_id].pos1 = pos2;
								mates_pairs[read_id].mism1 = mism1;
								mates_pairs[read_id].pos2 = pos1;
								mates_pairs[read_id].mism2 = mism2;
								mates_pairs[read_id].min_mism = min_mism;
								mates_pairs[read_id].strand = '-';//by the first read
							}//if
							else if(mates_pairs[read_id].pos1 != pos2 || 
								mates_pairs[read_id].pos2 != pos1 || mates_pairs[read_id].strand != astrand){

								if((min_mism < mates_pairs[read_id].min_mism) || 
									(min_mism == mates_pairs[read_id].min_mism && 
									second_mism < cur_second)){
									//unique
									mates_pairs[read_id].mark_unique = 1;
									mates_pairs[read_id].pos1 = pos2;
									mates_pairs[read_id].strand = '-';//by the first read
									mates_pairs[read_id].mism1 = mism1;
									mates_pairs[read_id].pos2 = pos1;
									mates_pairs[read_id].mism2 = mism2;
									mates_pairs[read_id].min_mism = min_mism;
										
								}
								else if(min_mism == mates_pairs[read_id].min_mism &&
									second_mism == cur_second){
									mates_pairs[read_id].mark_unique++;
									if(min_mism == 0 && second_mism == 0)
										pair_ambiguous = true;
								}
							}//else if


					
						}//both reads passed cg-check
					}//if correct insert range
					else if(dist > max_insert)
						break;
					}//if dif1 <= pos2
					if(pair_ambiguous == true)
						break;
				}//while cur < end2 && cur <= max_insert
			start1++;

		if(pair_ambiguous == true){
			break;
		}//if
	}//while
	}//if pair is not ambiguous



}//pair_reads new

////////////////////////////////LONG READs////////////////////////////////////////////

void Anchor::pair_reads_long(unsigned long int read_id, string aread1, string aread2,
				long int &contain_Ns, Stat &astat, short int read_len1, short int read_len2, int t, int rs )
{



	unsigned long int one = 1;

	long int asize,   pos1, pos2, dist, i, j;
	string strand, strand2;
	unsigned short int mism1, mism2;
	unsigned long int ta_lmer1, ta_lmer2, cg_lmer1, cg_lmer2;
	string line;

	short int len1 = read_len1 & MASK_SIZE;//= (read_len1 mod SIZE) 
	short int len2 = read_len2 & MASK_SIZE;

	if(len1 == 0)
		len1 = SIZE;
	if(len2 == 0)
		len2 = SIZE;

	unsigned long int MASK_READ1, MASK_READ2;
	if(len1 < SIZE){
		MASK_READ1 = ( 1 << len1) - 1;
		if(len1 > 31){
			dist = len1 - 32;
			MASK_READ1 = (((one << dist) - 1) << 32) | 0x00000000FFFFFFFF;
		}

	}
	else
		MASK_READ1 = 0xFFFFFFFFFFFFFFFF;

	if(len2 < SIZE){
		MASK_READ2 = (1 << len2) -1;
		if(len2 > 31){
			dist = len2 - 32;
			MASK_READ2 = (((one << dist) - 1) << 32) | 0x00000000FFFFFFFF;
			
		}//if
	}//if
	else
		MASK_READ2 = 0xFFFFFFFFFFFFFFFF;


	unsigned long int ta_rc1=0, ta_rc2=0, ind1 = 0, ind2 = 0;//hash-indices 
	vector<unsigned long int> c_rc_read1, c_rc_read2, c_read1, c_read2;//cg-representation
	vector<unsigned long int>  ta_read1, ta_read2, ta_rc_read1, ta_rc_read2;//ta-representation


	bool pair_ambiguous = false;


	//ta-reads
	convert_lmer_ta_long(aread1, ta_read1, read_len1);
	convert_lmer_ta_long(aread2, ta_read2, read_len2);

	int ta_size1 = ta_read1.size();
	int ta_size2 = ta_read2.size();

	
	//c-reads
	convert_lmer_c_long(aread1, c_read1,  read_len1);
	convert_lmer_c_long(aread2, c_read2,  read_len2);

	


  //rev-complement of c-reads

	int dif1 = read_len1 - len_mer;
	int dif2 = read_len2 - len_mer;//to find genome kmers


	//we need the first len_mer bases for each read
	if(ta_size1 > 1){


	if(dif_len_mer > 0){
		ta_lmer1 = (ta_read1[0] >> dif_len_mer) & seeds[t][rs];
		cg_lmer1 = (c_read1[0] >> dif_len_mer) & seeds[t][rs];
	}
	else{
		ta_lmer1 = ta_read1[0]  & seeds[t][rs];
		cg_lmer1 = c_read1[0]  & seeds[t][rs];
	
	}//else

	}//if
	else{
		ta_lmer1 = (ta_read1[0] >> dif1) & seeds[t][rs];//dif1 < 64
		cg_lmer1 = (c_read1[0] >> dif1) & seeds[t][rs];
	
	}
	if(ta_size2 > 1){
		if(dif_len_mer > 0){
			cg_lmer2 = (c_read2[0] >> dif_len_mer) & seeds2[t][rs];
			ta_lmer2 = (ta_read2[0] >> dif_len_mer) & seeds2[t][rs];
		}
		else{
			cg_lmer2 = c_read2[0]  & seeds2[t][rs];
			ta_lmer2 = ta_read2[0]  & seeds2[t][rs];

		}//else
	}//if
	else{
		cg_lmer2 = (c_read2[0] >> dif2) & seeds2[t][rs];
		ta_lmer2 = (ta_read2[0] >> dif2) & seeds2[t][rs];
	
	}
	ind1 = hash_aread(ta_lmer1,  cg_lmer1, false);
	ind2 = hash_aread(ta_lmer2, cg_lmer2, true);




	//rev-complement of g-reads

	//rev-complement hashing on the first len-mer bases
	ta_rc_read2.resize(ta_size2);
	c_rc_read2.resize(ta_size2);




	int stop = ta_size2 - 1;

	j= stop;
	for(i = 0; i < stop; i++){
		ta_rc_read2[j] = ~Reverse(ta_read2[i], SIZE);
		c_rc_read2[j] =  Reverse(c_read2[i], SIZE);//shows positions of Gs in reverse strand
		j--;
	}
	ta_rc_read2[0] = MASK_READ2 & (~Reverse(ta_read2[stop], len2));
	c_rc_read2[0] =  MASK_READ2 & Reverse(c_read2[stop], len2);//shows positions of Gs in reverse strand

	ta_lmer2 = (ta_rc_read2[stop] & MASK_LMER) & seeds2[t][rs];
	cg_lmer2 = (c_rc_read2[stop] & MASK_LMER) & seeds2[t][rs];

	if(ta_size2 > 1){
		if(len2 < SIZE){
			int residue = SIZE - len2;
			for(i = 0; i < stop; i++){
				ta_rc_read2[i] = (ta_rc_read2[i] << residue) | (ta_rc_read2[i+1] >> len2);
				c_rc_read2[i] = (c_rc_read2[i] << residue) | (c_rc_read2[i+1] >> len2);
			}//for i
			unsigned long int MASK2 = (one << len2) - 1 ;
			ta_rc_read2[stop] = ta_rc_read2[stop] & MASK2;
			c_rc_read2[stop] = c_rc_read2[stop] & MASK2;
		}//if len2 < SIZe
	
	}//ta_size2 > 1


	ta_rc2 =  hash_aread(ta_lmer2, cg_lmer2, false);



	//process forward read1 and rev read2
	//start and end of list of genome positions for forw1
	long int start1 = 0;
	long int end1 = pos_sizes[ind1];

	//start and end of list of genome positions for rev-compl2

	long int start2 = 0;
	long int end2 = pos_sizes[ta_rc2];
	long int cur = start2;

	while(start1 < end1){
		if(start2 >= end2)
			break;

		pos1 = pos_index[ind1][start1];

			cur = start2;
			for(; cur != end2 ; cur++){
				pos2 = pos_index[ta_rc2][cur];
				if(dif2 <= pos2){
					pos2 -= dif2;
					if(A_rich == true)
						dist = pos1 - pos2;
					else
						dist = pos2 - pos1;
					
					if(dist >= min_insert ){
						break;
					}
				}//if dif2 <= pos2
			}//for cur

			start2 = cur;


			for(; cur != end2; cur++){
				pos2 = pos_index[ta_rc2][cur] ;
				if(dif2 <= pos2){

					pos2 -= dif2;

					if(A_rich == true)
						dist = pos1 - pos2;
					else
						dist = pos2 - pos1;

				if(dist >= min_insert && dist <= max_insert ){
					//cg check
					bool cg_mism1 = false, cg_mism2 = false;
					vector<unsigned long int> ta_gen_mer1;

					convert_gen_mer_ta_long(ta_gen_mer1, pos1, read_len1, ta_size1);


					vector<unsigned long int> ta_check , ta_check2;
					mism1 = count_mism_long(ta_read1, ta_gen_mer1, ta_check, ta_size1);

					if(mism1 <= mism_thresh){
						vector<unsigned long int> ta_gen_mer2;
						
						convert_gen_mer_ta_long(ta_gen_mer2, pos2, read_len2, ta_size2);

						mism2 = count_mism_long(ta_rc_read2, ta_gen_mer2, ta_check2, ta_size2);

						if(mism2 <= mism_thresh){
							vector<unsigned long int> cg_gen_mer1, cg_gen_mer2, ns_gen_mer1, ns_gen_mer2;

							convert_gen_mer_cg_long(cg_gen_mer1, pos1,  read_len1, ta_size1);
							convert_gen_mer_cg_long(cg_gen_mer2, pos2, read_len2, ta_size2);

							convert_gen_mer_ns_long(ns_gen_mer1, pos1, read_len1, ta_size1);
							convert_gen_mer_ns_long(ns_gen_mer2, pos2, read_len2, ta_size2);

							for(j = 0; j < ta_size1; j++){
								ta_check[j] |= ((~(ta_read1[j] | c_read1[j])) & cg_gen_mer1[j]) | ns_gen_mer1[j];
							}
							for(j = 0; j < ta_size2; j++){
								ta_check2[j] |=  ((~(ta_rc_read2[j] | c_rc_read2[j])) & cg_gen_mer2[j]) | ns_gen_mer2[j];

							}

								for(j = 0; j < ta_size1; j++){
									ta_check[j] |= c_read1[j] ^ (c_read1[j] & cg_gen_mer1[j]);

								}//for j

								for(j = 0; j < ta_size2; j++)
										ta_check2[j] |= c_rc_read2[j] ^ (c_rc_read2[j] & cg_gen_mer2[j]);

							mism1 = 0;
							for(j = 0; j < ta_size1; j++){
								mism1 += count_mism(ta_check[j]);
							}

							if(mism1 > mism_thresh)
								cg_mism1 = true;

							mism2 = 0;
							for(j = 0; j < ta_size2; j++){
								mism2 += count_mism(ta_check2[j]);
							}//for j
							if(mism2 > mism_thresh)
								cg_mism2 = true;

						}//if passed second ta-check
						else
							cg_mism2 = true;//not passed 2 ta-check
					}//if passed ta=check
					else{
						cg_mism1 = true;
					}//else



					if(cg_mism1 == false && cg_mism2 == false){
						//T->C and A->G not in the same read check
						//read1 check
						char astrand = '+';
						//check tc/ag mism don't happen in mates
							int min_mism = mism1;
							int second_mism = mism2;
							if(mism2 < mism1){
								min_mism = mism2;
								second_mism = mism1;
							}
							int cur_second = mates_pairs[read_id].min_mism;
							if(cur_second <  mates_pairs[read_id].mism1)
								cur_second = mates_pairs[read_id].mism1;
							else if(cur_second < mates_pairs[read_id].mism2)
								cur_second = mates_pairs[read_id].mism2;

							if(mates_pairs[read_id].mark_unique == 0){
								mates_pairs[read_id].mark_unique = 1;
								mates_pairs[read_id].pos1 = pos1;
								mates_pairs[read_id].mism1 = mism1;
								mates_pairs[read_id].pos2 = pos2;
								mates_pairs[read_id].mism2 = mism2;
								mates_pairs[read_id].min_mism = min_mism;
								mates_pairs[read_id].strand = '+';//by the first read
							}//if
							else if(mates_pairs[read_id].pos1 != pos1 || 
								mates_pairs[read_id].pos2 != pos2 || mates_pairs[read_id].strand != astrand){


								if((min_mism < mates_pairs[read_id].min_mism) || 
									(min_mism == mates_pairs[read_id].min_mism && 
									second_mism < cur_second)){
									//unique
		
									mates_pairs[read_id].mark_unique = 1;
									mates_pairs[read_id].pos1 = pos1;
									mates_pairs[read_id].strand = '+';//by the first read
									mates_pairs[read_id].mism1 = mism1;
									mates_pairs[read_id].pos2 = pos2;
									mates_pairs[read_id].mism2 = mism2;
									mates_pairs[read_id].min_mism = min_mism;
										
								}
								else if(min_mism == mates_pairs[read_id].min_mism &&
									second_mism == cur_second){

									mates_pairs[read_id].mark_unique++;
									if(min_mism == 0 && second_mism == 0)
										pair_ambiguous = true;
								}
							}//else if

					}//both reads passed cg-check
				}//if correct insert range

				else if(dist > max_insert)
					break;

				}//if dif2 <= pos2

				if(pair_ambiguous == true){					
					break;
				}
			}//while cur < end2 && cur <= max_insert
			
			start1++;

		if(pair_ambiguous == true){
			break;
		}//if


	}//for forward read's list




	if(pair_ambiguous == false){
		//take rev-complement of the first read

		//rev-complement hashing on the first len-mer bases
		ta_rc_read1.resize(ta_size1);
		c_rc_read1.resize(ta_size1);

		int stop = ta_size1 - 1;
		j= stop;
		for(i = 0; i < stop; i++){
			ta_rc_read1[j] = ~Reverse(ta_read1[i], SIZE);
			c_rc_read1[j] =  Reverse(c_read1[i], SIZE);//shows positions of Gs in reverse strand
			j--;
		}
		ta_rc_read1[0] = MASK_READ1 & (~Reverse(ta_read1[stop], len1));
		c_rc_read1[0] =  MASK_READ1 & Reverse(c_read1[stop], len1);//shows positions of Gs in reverse strand

		ta_lmer1 = (ta_rc_read1[stop] & MASK_LMER) & seeds[t][rs];
		cg_lmer1 = (c_rc_read1[stop] & MASK_LMER) & seeds[t][rs];

		if(ta_size1 > 1){
			int residue = SIZE - len1;
			if(len1 < SIZE){
				for(i = 0; i < ta_size1 -1; i++){
					ta_rc_read1[i] = (ta_rc_read1[i] << residue) | (ta_rc_read1[i+1] >> len1);
					c_rc_read1[i] = (c_rc_read1[i] << residue) | (c_rc_read1[i+1] >> len1);
				}//for i
				unsigned long int MASK1 = (one << len1) - 1 ;
				ta_rc_read1[stop] = ta_rc_read1[stop] & MASK1;
				c_rc_read1[stop] = c_rc_read1[stop] & MASK1;
			}//if len1 < SIZE
		
		}//ta size1


		ta_rc1 = hash_aread(ta_lmer1, cg_lmer1, true);

		//process rev-complemnet of read1 and forward of read2

		start1 = 0;
		end1 = pos_sizes[ind2];

		start2 = 0;
		end2 = pos_sizes[ta_rc1];

		while(start1 < end1){

			if(start2 >= end2)
				break;

			pos1 =  pos_index[ind2][start1];

				cur = start2;
				for(; cur != end2 ; cur++){
					pos2 =  pos_index[ta_rc1][cur] ;
					if(	dif1 <= pos2){		
						pos2 -= dif1;
						if(A_rich == true)
							dist = pos1 - pos2;
						else
							dist = pos2 - pos1;
		
						
						if(dist >= min_insert ){
							break;
						}//if
					}//if
				}//for cur

				start2 = cur;

					for(; cur != end2; cur++){
						pos2 =  pos_index[ta_rc1][cur] ;

						if(dif1 <= pos2){
						pos2 -= dif1;
						if(A_rich == true)
							dist = pos1 - pos2;
						else
							dist = pos2 - pos1;
					
						if(dist >= min_insert && dist <= max_insert ){
					//cg check
						vector<unsigned long int> ta_gen_mer2, ta_gen_mer1;

						bool cg_mism1 = false, cg_mism2 = false;
						convert_gen_mer_ta_long(ta_gen_mer2, pos1, read_len2, ta_size2);


					vector<unsigned long int> ta_check2 , ta_check;
					mism2 = count_mism_long(ta_read2, ta_gen_mer2, ta_check2, ta_size2);

					if(mism2 <= mism_thresh){
						convert_gen_mer_ta_long(ta_gen_mer1, pos2, read_len1, ta_size1);

						mism1 = count_mism_long(ta_rc_read1, ta_gen_mer1, ta_check, ta_size1);

						if(mism1 <= mism_thresh){
							vector<unsigned long int> cg_gen_mer1, cg_gen_mer2, ns_gen_mer1, ns_gen_mer2;

							convert_gen_mer_cg_long( cg_gen_mer2, pos1,  read_len2, ta_size2);
							convert_gen_mer_cg_long( cg_gen_mer1, pos2,  read_len1, ta_size1);
						
							convert_gen_mer_ns_long(ns_gen_mer2, pos1,  read_len2, ta_size2);
							convert_gen_mer_ns_long(ns_gen_mer1, pos2, read_len1, ta_size1);

							for(j = 0; j < ta_size2; j++){
				ta_check2[j] |= ((ta_read2[j] & (~c_read2[j])) & cg_gen_mer2[j]) | ns_gen_mer2[j];
							}
								

							for(j = 0; j < ta_size1; j++){
			ta_check[j] |= ((ta_rc_read1[j] &(~c_rc_read1[j])) & cg_gen_mer1[j]) | ns_gen_mer1[j];
							
							}


								for(j = 0; j < ta_size2; j++)
									ta_check2[j] |= c_read2[j] ^ (c_read2[j] & cg_gen_mer2[j]);

								for(j = 0; j < ta_size1; j++)
									ta_check[j] |= c_rc_read1[j] ^ (c_rc_read1[j] & cg_gen_mer1[j]);


								
							mism2 = 0;
							for(j = 0; j < ta_size2; j++)
								mism2 += count_mism(ta_check2[j]);
							
							if(mism2 > mism_thresh)
								cg_mism2 = true;

							mism1 = 0;
							for(j = 0; j < ta_size1; j++)
								mism1 += count_mism(ta_check[j]);

							if(mism1 > mism_thresh)
								cg_mism1 = true;



						}//if second ta-check passed
						else
							cg_mism2 = true;
						}//if first ta_check passed
						else{
							cg_mism1 = true;
						}//else


							if(cg_mism1 == false && cg_mism2 == false){

							char astrand = '-';
							int min_mism = mism1;
							int second_mism = mism2;
							if(mism2 < mism1){
								min_mism = mism2;
								second_mism = mism1;
							}
							int cur_second = mates_pairs[read_id].min_mism;
							if(cur_second <  mates_pairs[read_id].mism1)
								cur_second = mates_pairs[read_id].mism1;
							else if(cur_second < mates_pairs[read_id].mism2)
								cur_second = mates_pairs[read_id].mism2;

							if(mates_pairs[read_id].mark_unique == 0){
								mates_pairs[read_id].mark_unique = 1;
								mates_pairs[read_id].pos1 = pos2;
								mates_pairs[read_id].mism1 = mism1;
								mates_pairs[read_id].pos2 = pos1;
								mates_pairs[read_id].mism2 = mism2;
								mates_pairs[read_id].min_mism = min_mism;
								mates_pairs[read_id].strand = '-';//by the first read
							}//if
							else if(mates_pairs[read_id].pos1 != pos2 || 
								mates_pairs[read_id].pos2 != pos1 || mates_pairs[read_id].strand != astrand){

								if((min_mism < mates_pairs[read_id].min_mism) || 
									(min_mism == mates_pairs[read_id].min_mism && 
									second_mism < cur_second)){
									//unique
									mates_pairs[read_id].mark_unique = 1;
									mates_pairs[read_id].pos1 = pos2;
									mates_pairs[read_id].strand = '-';//by the first read
									mates_pairs[read_id].mism1 = mism1;
									mates_pairs[read_id].pos2 = pos1;
									mates_pairs[read_id].mism2 = mism2;
									mates_pairs[read_id].min_mism = min_mism;
										
								}
								else if(min_mism == mates_pairs[read_id].min_mism &&
									second_mism == cur_second){
									mates_pairs[read_id].mark_unique++;
									if(min_mism == 0 && second_mism == 0)
										pair_ambiguous = true;
								}
							}//else if


					
						}//both reads passed cg-check
					}//if correct insert range
					else if(dist > max_insert)
						break;
					}//if dif1 <= pos2
					if(pair_ambiguous == true)
						break;
				}//while cur < end2 && cur <= max_insert
			start1++;


		if(pair_ambiguous == true){
			break;
		}//if


	}//while
	}//if pair is not ambiguous



}//pair_reads_long new



void Anchor::pair_reads_normal_long(unsigned long int read_id, string aread1, string aread2,
				long int &contain_Ns, Stat &astat, short int read_len1, short int read_len2, int t, int rs )
{



	unsigned long int one = 1;

	long int asize,   pos1, pos2, dist, i, j;
	string strand, strand2;
	unsigned short int mism1, mism2;
	unsigned long int ta_lmer1, ta_lmer2, cg_lmer1, cg_lmer2;
	string line;

	short int len1 = read_len1 & MASK_SIZE;//= (read_len1 mod SIZE) 
	short int len2 = read_len2 & MASK_SIZE;

	if(len1 == 0)
		len1 = SIZE;
	if(len2 == 0)
		len2 = SIZE;

	unsigned long int MASK_READ1, MASK_READ2;
	if(len1 < SIZE){
		MASK_READ1 = ( 1 << len1) - 1;
		if(len1 > 31){
			dist = len1 - 32;
			MASK_READ1 = (((one << dist) - 1) << 32) | 0x00000000FFFFFFFF;
		}

	}
	else
		MASK_READ1 = 0xFFFFFFFFFFFFFFFF;

	if(len2 < SIZE){
		MASK_READ2 = (1 << len2) -1;
		if(len2 > 31){
			dist = len2 - 32;
			MASK_READ2 = (((one << dist) - 1) << 32) | 0x00000000FFFFFFFF;
			
		}//if
	}//if
	else
		MASK_READ2 = 0xFFFFFFFFFFFFFFFF;


	unsigned long int ta_rc1=0, ta_rc2=0, ind1 = 0, ind2 = 0;//hash-indices 
	vector<unsigned long int> c_rc_read1, c_rc_read2, c_read1, c_read2;//cg-representation
	vector<unsigned long int>  ta_read1, ta_read2, ta_rc_read1, ta_rc_read2;//ta-representation


	bool pair_ambiguous = false;


	//ta-reads
	convert_lmer_ta_long(aread1, ta_read1, read_len1);
	convert_lmer_ta_long(aread2, ta_read2, read_len2);

	int ta_size1 = ta_read1.size();
	int ta_size2 = ta_read2.size();

	
	//c-reads
	convert_lmer_c_long(aread1, c_read1,  read_len1);
	convert_lmer_c_long(aread2, c_read2,  read_len2);

	


  //rev-complement of c-reads

	int dif1 = read_len1 - len_mer;
	int dif2 = read_len2 - len_mer;//to find genome kmers


	//we need the first len_mer bases for each read
	if(ta_size1 > 1){


	if(dif_len_mer > 0){
		ta_lmer1 = (ta_read1[0] >> dif_len_mer) & seeds[t][rs];
		cg_lmer1 = (c_read1[0] >> dif_len_mer) & seeds[t][rs];
	}
	else{
		ta_lmer1 = ta_read1[0]  & seeds[t][rs];
		cg_lmer1 = c_read1[0]  & seeds[t][rs];
	
	}//else

	}//if
	else{
		ta_lmer1 = (ta_read1[0] >> dif1) & seeds[t][rs];//dif1 < 64
		cg_lmer1 = (c_read1[0] >> dif1) & seeds[t][rs];
	
	}
	if(ta_size2 > 1){
		if(dif_len_mer > 0){
			cg_lmer2 = (c_read2[0] >> dif_len_mer) & seeds2[t][rs];
			ta_lmer2 = (ta_read2[0] >> dif_len_mer) & seeds2[t][rs];
		}
		else{
			cg_lmer2 = c_read2[0]  & seeds2[t][rs];
			ta_lmer2 = ta_read2[0]  & seeds2[t][rs];

		}//else
	}//if
	else{
		cg_lmer2 = (c_read2[0] >> dif2) & seeds2[t][rs];
		ta_lmer2 = (ta_read2[0] >> dif2) & seeds2[t][rs];
	
	}

	ind1 = hash2c(ta_lmer1, cg_lmer1 );
	ind2 = hash2c(ta_lmer2, cg_lmer2);


	//rev-complement of g-reads

	//rev-complement hashing on the first len-mer bases
	ta_rc_read2.resize(ta_size2);
	c_rc_read2.resize(ta_size2);




	int stop = ta_size2 - 1;

	j= stop;
	for(i = 0; i < stop; i++){
		ta_rc_read2[j] = ~Reverse(ta_read2[i], SIZE);
		c_rc_read2[j] =  Reverse(c_read2[i], SIZE);//shows positions of Gs in reverse strand
		j--;
	}
	ta_rc_read2[0] = MASK_READ2 & (~Reverse(ta_read2[stop], len2));
	c_rc_read2[0] =  MASK_READ2 & Reverse(c_read2[stop], len2);//shows positions of Gs in reverse strand

	ta_lmer2 = (ta_rc_read2[stop] & MASK_LMER) & seeds2[t][rs];
	cg_lmer2 = (c_rc_read2[stop] & MASK_LMER) & seeds2[t][rs];

	if(ta_size2 > 1){
		if(len2 < SIZE){
			int residue = SIZE - len2;
			for(i = 0; i < stop; i++){
				ta_rc_read2[i] = (ta_rc_read2[i] << residue) | (ta_rc_read2[i+1] >> len2);
				c_rc_read2[i] = (c_rc_read2[i] << residue) | (c_rc_read2[i+1] >> len2);
			}//for i
			unsigned long int MASK2 = (one << len2) - 1 ;
			ta_rc_read2[stop] = ta_rc_read2[stop] & MASK2;
			c_rc_read2[stop] = c_rc_read2[stop] & MASK2;
		}//if len2 < SIZe
	
	}//ta_size2 > 1


	ta_rc2 = hash2c(ta_lmer2, cg_lmer2);



	//process forward read1 and rev read2
	//start and end of list of genome positions for forw1
	long int start1 = 0;
	long int end1 = pos_sizes[ind1];

	//start and end of list of genome positions for rev-compl2

	long int start2 = 0;
	long int end2 = pos_sizes[ta_rc2];
	long int cur = start2;

	while(start1 < end1){
		if(start2 >= end2)
			break;

		pos1 = pos_index[ind1][start1];
			cur = start2;
			for(; cur != end2 ; cur++){
				pos2 = pos_index[ta_rc2][cur];
				if(dif2 <= pos2){
					pos2 -= dif2;

					dist = pos2 - pos1;
					if(dist >= min_insert ){
						break;
					}
				}//if dif2 <= pos2
			}//for cur

			start2 = cur;


			for(; cur != end2; cur++){
				pos2 = pos_index[ta_rc2][cur] ;

			

				if(dif2 <= pos2){

					pos2 -= dif2;


				dist = pos2 - pos1;
				if(dist >= min_insert && dist <= max_insert ){
					//cg check
					bool cg_mism1 = false, cg_mism2 = false;
					vector<unsigned long int> ta_gen_mer1;

					convert_gen_mer_ta_long(ta_gen_mer1, pos1, read_len1, ta_size1);


					vector<unsigned long int> ta_check , ta_check2;
					mism1 = count_mism_long(ta_read1, ta_gen_mer1, ta_check, ta_size1);

					if(mism1 <= mism_thresh){
						vector<unsigned long int> ta_gen_mer2;
						
						convert_gen_mer_ta_long(ta_gen_mer2, pos2,  read_len2, ta_size2);

						mism2 = count_mism_long(ta_rc_read2, ta_gen_mer2, ta_check2, ta_size2);

						if(mism2 <= mism_thresh){
							vector<unsigned long int> cg_gen_mer1, cg_gen_mer2, ns_gen_mer1, ns_gen_mer2;

							convert_gen_mer_cg_long(cg_gen_mer1, pos1, read_len1, ta_size1);
							convert_gen_mer_cg_long(cg_gen_mer2, pos2,  read_len2, ta_size2);

							convert_gen_mer_ns_long(ns_gen_mer1, pos1,  read_len1, ta_size1);
							convert_gen_mer_ns_long(ns_gen_mer2, pos2,  read_len2, ta_size2);

							for(j = 0; j < ta_size1; j++){
								ta_check[j] |= ns_gen_mer1[j];
							}
							for(j = 0; j < ta_size2; j++){
								ta_check2[j] |= ns_gen_mer2[j];

							}

								for(j = 0; j < ta_size1; j++){
									ta_check[j] |= c_read1[j] ^  cg_gen_mer1[j];

								}//for j
								for(j = 0; j < ta_size2; j++)
										ta_check2[j] |= c_rc_read2[j] ^  cg_gen_mer2[j];

							mism1 = 0;
							for(j = 0; j < ta_size1; j++){
								mism1 += count_mism(ta_check[j]);
							}

							if(mism1 > mism_thresh)
								cg_mism1 = true;

							mism2 = 0;
							for(j = 0; j < ta_size2; j++){
								mism2 += count_mism(ta_check2[j]);
							}//for j
							if(mism2 > mism_thresh)
								cg_mism2 = true;

						}//if passed second ta-check
						else
							cg_mism2 = true;//not passed 2 ta-check
					}//if passed ta=check
					else{
						cg_mism1 = true;
					}//else



					if(cg_mism1 == false && cg_mism2 == false){
						//T->C and A->G not in the same read check
						//read1 check
						char astrand = '+';
						//check tc/ag mism don't happen in mates
							int min_mism = mism1;
							int second_mism = mism2;
							if(mism2 < mism1){
								min_mism = mism2;
								second_mism = mism1;
							}
							int cur_second = mates_pairs[read_id].min_mism;
							if(cur_second <  mates_pairs[read_id].mism1)
								cur_second = mates_pairs[read_id].mism1;
							else if(cur_second < mates_pairs[read_id].mism2)
								cur_second = mates_pairs[read_id].mism2;

							if(mates_pairs[read_id].mark_unique == 0){
								mates_pairs[read_id].mark_unique = 1;
								mates_pairs[read_id].pos1 = pos1;
								mates_pairs[read_id].mism1 = mism1;
								mates_pairs[read_id].pos2 = pos2;
								mates_pairs[read_id].mism2 = mism2;
								mates_pairs[read_id].min_mism = min_mism;
								mates_pairs[read_id].strand = '+';//by the first read
							}//if
							else if(mates_pairs[read_id].pos1 != pos1 || 
								mates_pairs[read_id].pos2 != pos2 || mates_pairs[read_id].strand != astrand){


								if((min_mism < mates_pairs[read_id].min_mism) || 
									(min_mism == mates_pairs[read_id].min_mism && 
									second_mism < cur_second)){
									//unique
		
									mates_pairs[read_id].mark_unique = 1;
									mates_pairs[read_id].pos1 = pos1;
									mates_pairs[read_id].strand = '+';//by the first read
									mates_pairs[read_id].mism1 = mism1;
									mates_pairs[read_id].pos2 = pos2;
									mates_pairs[read_id].mism2 = mism2;
									mates_pairs[read_id].min_mism = min_mism;
										
								}
								else if(min_mism == mates_pairs[read_id].min_mism &&
									second_mism == cur_second){

									mates_pairs[read_id].mark_unique++;
									if(min_mism == 0 && second_mism == 0)
										pair_ambiguous = true;
								}
							}//else if

					}//both reads passed cg-check
				}//if correct insert range

				else if(dist > max_insert)
					break;

				}//if dif2 <= pos2

				if(pair_ambiguous == true){					
					break;
				}
			}//while cur < end2 && cur <= max_insert
			
			start1++;



		if(pair_ambiguous == true){
			break;
		}//if


	}//for forward read's list




	if(pair_ambiguous == false){
		//take rev-complement of the first read

		//rev-complement hashing on the first len-mer bases
		ta_rc_read1.resize(ta_size1);
		c_rc_read1.resize(ta_size1);

		int stop = ta_size1 - 1;
		j= stop;
		for(i = 0; i < stop; i++){
			ta_rc_read1[j] = ~Reverse(ta_read1[i], SIZE);
			c_rc_read1[j] =  Reverse(c_read1[i], SIZE);//shows positions of Gs in reverse strand
			j--;
		}
		ta_rc_read1[0] = MASK_READ1 & (~Reverse(ta_read1[stop], len1));
		c_rc_read1[0] =  MASK_READ1 & Reverse(c_read1[stop], len1);//shows positions of Gs in reverse strand

		ta_lmer1 = (ta_rc_read1[stop] & MASK_LMER) & seeds[t][rs];
		cg_lmer1 = (c_rc_read1[stop] & MASK_LMER) & seeds[t][rs];

		if(ta_size1 > 1){
			int residue = SIZE - len1;
			if(len1 < SIZE){
				for(i = 0; i < ta_size1 -1; i++){
					ta_rc_read1[i] = (ta_rc_read1[i] << residue) | (ta_rc_read1[i+1] >> len1);
					c_rc_read1[i] = (c_rc_read1[i] << residue) | (c_rc_read1[i+1] >> len1);
				}//for i
				unsigned long int MASK1 = (one << len1) - 1 ;
				ta_rc_read1[stop] = ta_rc_read1[stop] & MASK1;
				c_rc_read1[stop] = c_rc_read1[stop] & MASK1;
			}//if len1 < SIZE
		
		}//ta size1


		ta_rc1 = hash2c(ta_lmer1, cg_lmer1);

		//process rev-complemnet of read1 and forward of read2

		start1 = 0;
		end1 = pos_sizes[ind2];

		start2 = 0;
		end2 = pos_sizes[ta_rc1];

		while(start1 < end1){

			if(start2 >= end2)
				break;

			pos1 =  pos_index[ind2][start1];


				cur = start2;
				for(; cur != end2 ; cur++){
					pos2 =  pos_index[ta_rc1][cur];
					if(	dif1 <= pos2){		
						pos2 -= dif1;
						dist = pos2 - pos1;
						
						if(dist >= min_insert ){
							break;
						}//if
					}//if
				}//for cur

				start2 = cur;

					for(; cur != end2; cur++){
						pos2 =  pos_index[ta_rc1][cur] ;

						if(dif1 <= pos2){
						pos2 -= dif1;
						dist = pos2 - pos1;
						if(dist >= min_insert && dist <= max_insert ){
					//cg check
						vector<unsigned long int> ta_gen_mer2, ta_gen_mer1;

						bool cg_mism1 = false, cg_mism2 = false;
						convert_gen_mer_ta_long(ta_gen_mer2, pos1, read_len2, ta_size2);


					vector<unsigned long int> ta_check2 , ta_check;
					mism2 = count_mism_long(ta_read2, ta_gen_mer2, ta_check2, ta_size2);

					if(mism2 <= mism_thresh){
						convert_gen_mer_ta_long(ta_gen_mer1, pos2,  read_len1, ta_size1);

						mism1 = count_mism_long(ta_rc_read1, ta_gen_mer1, ta_check, ta_size1);

						if(mism1 <= mism_thresh){
							vector<unsigned long int> cg_gen_mer1, cg_gen_mer2, ns_gen_mer1, ns_gen_mer2;

							convert_gen_mer_cg_long( cg_gen_mer2, pos1,  read_len2, ta_size2);
							convert_gen_mer_cg_long( cg_gen_mer1, pos2,  read_len1, ta_size1);
						
							convert_gen_mer_ns_long(ns_gen_mer2, pos1, read_len2, ta_size2);
							convert_gen_mer_ns_long(ns_gen_mer1, pos2,  read_len1, ta_size1);

							for(j = 0; j < ta_size2; j++){
								ta_check2[j] |=  ns_gen_mer2[j];
							}
								

							for(j = 0; j < ta_size1; j++){
								ta_check[j] |=  ns_gen_mer1[j];
							
							}


								for(j = 0; j < ta_size2; j++)
									ta_check2[j] |= c_read2[j] ^ cg_gen_mer2[j];

								for(j = 0; j < ta_size1; j++)
									ta_check[j] |= c_rc_read1[j] ^ cg_gen_mer1[j];


								
							mism2 = 0;
							for(j = 0; j < ta_size2; j++)
								mism2 += count_mism(ta_check2[j]);
							
							if(mism2 > mism_thresh)
								cg_mism2 = true;

							mism1 = 0;
							for(j = 0; j < ta_size1; j++)
								mism1 += count_mism(ta_check[j]);

							if(mism1 > mism_thresh)
								cg_mism1 = true;



						}//if second ta-check passed
						else
							cg_mism2 = true;
						}//if first ta_check passed
						else{
							cg_mism1 = true;
						}//else


							if(cg_mism1 == false && cg_mism2 == false){

							char astrand = '-';
							int min_mism = mism1;
							int second_mism = mism2;
							if(mism2 < mism1){
								min_mism = mism2;
								second_mism = mism1;
							}
							int cur_second = mates_pairs[read_id].min_mism;
							if(cur_second <  mates_pairs[read_id].mism1)
								cur_second = mates_pairs[read_id].mism1;
							else if(cur_second < mates_pairs[read_id].mism2)
								cur_second = mates_pairs[read_id].mism2;

							if(mates_pairs[read_id].mark_unique == 0){
								mates_pairs[read_id].mark_unique = 1;
								mates_pairs[read_id].pos1 = pos2;
								mates_pairs[read_id].mism1 = mism1;
								mates_pairs[read_id].pos2 = pos1;
								mates_pairs[read_id].mism2 = mism2;
								mates_pairs[read_id].min_mism = min_mism;
								mates_pairs[read_id].strand = '-';//by the first read
							}//if
							else if(mates_pairs[read_id].pos1 != pos2 || 
								mates_pairs[read_id].pos2 != pos1 || mates_pairs[read_id].strand != astrand){

								if((min_mism < mates_pairs[read_id].min_mism) || 
									(min_mism == mates_pairs[read_id].min_mism && 
									second_mism < cur_second)){
									//unique
									mates_pairs[read_id].mark_unique = 1;
									mates_pairs[read_id].pos1 = pos2;
									mates_pairs[read_id].strand = '-';//by the first read
									mates_pairs[read_id].mism1 = mism1;
									mates_pairs[read_id].pos2 = pos1;
									mates_pairs[read_id].mism2 = mism2;
									mates_pairs[read_id].min_mism = min_mism;
										
								}
								else if(min_mism == mates_pairs[read_id].min_mism &&
									second_mism == cur_second){
									mates_pairs[read_id].mark_unique++;
									if(min_mism == 0 && second_mism == 0)
										pair_ambiguous = true;
								}
							}//else if


					
						}//both reads passed cg-check
					}//if correct insert range
					else if(dist > max_insert)
						break;
					}//if dif1 <= pos2
					if(pair_ambiguous == true)
						break;
				}//while cur < end2 && cur <= max_insert
			start1++;

		if(pair_ambiguous == true){
			break;
		}//if



	}//while
	}//if pair is not ambiguous



}//pair_reads_normal_long new



void Anchor::pair_reads_singles_long(unsigned long int read_id, string aread1,
				long int &contain_Ns, Stat &astat, short int read_len1,  int t, int rs )
{



	unsigned long int one = 1;

	long int asize,  pos1, dist, i, j;
	string strand;
	unsigned short int mism1;
	unsigned long int ta_lmer1, cg_lmer1;
	string line;

	short int len1 = read_len1 & MASK_SIZE;//= (read_len1 mod SIZE) 

	if(len1 == 0)
		len1 = SIZE;

	unsigned long int MASK_READ1, MASK_READ2;
	if(len1 < SIZE){
		MASK_READ1 = ( 1 << len1) - 1;
		if(len1 > 31){
			dist = len1 - 32;
			MASK_READ1 = (((one << dist) - 1) << 32) | 0x00000000FFFFFFFF;
		}

	}
	else
		MASK_READ1 = 0xFFFFFFFFFFFFFFFF;


	unsigned long int ta_rc1=0, ind1 = 0;//hash-indices 
	vector<unsigned long int> c_rc_read1, c_read1;//cg-representation
	vector<unsigned long int>  ta_read1, ta_rc_read1;//ta-representation

	bool pair_ambiguous = false;
	//ta-reads
	convert_lmer_ta_long(aread1, ta_read1, read_len1);

	int ta_size1 = ta_read1.size();


	int chr_id;
	
	//c-reads
	convert_lmer_c_long(aread1, c_read1,  read_len1);

  //rev-complement of c-reads

	int dif1 = read_len1 - len_mer;

	//we need the first len_mer bases for each read
	if(ta_size1 > 1){

		if(dif_len_mer > 0){
			ta_lmer1 = (ta_read1[0] >> dif_len_mer) & seeds[t][rs];
			cg_lmer1 = (c_read1[0] >> dif_len_mer) & seeds[t][rs];
		}
		else{
			ta_lmer1 = ta_read1[0]  & seeds[t][rs];
			cg_lmer1 = c_read1[0]  & seeds[t][rs];
		
		}//else

	}//if
	else{
		ta_lmer1 = (ta_read1[0] >> dif1) & seeds[t][rs];//dif1 < 64
		cg_lmer1 = (c_read1[0] >> dif1) & seeds[t][rs];
	
	}

	ind1 = hash_aread(ta_lmer1,  cg_lmer1, false);


	//rev-complement of g-reads


	//process forward read1 and rev read2
	//start and end of list of genome positions for forw1
	long int start1 = 0;
	long int end1 = pos_sizes[ind1];

	//start and end of list of genome positions for rev-compl2

	while(start1 < end1){

		pos1 = pos_index[ind1][start1];

			//cg check
			bool cg_mism1 = false;
			vector<unsigned long int> ta_gen_mer1;

			convert_gen_mer_ta_long(ta_gen_mer1, pos1, read_len1, ta_size1);


			vector<unsigned long int> ta_check;
			mism1 = count_mism_long(ta_read1, ta_gen_mer1, ta_check, ta_size1);


			if(mism1 <= mism_thresh){

				vector<unsigned long int> cg_gen_mer1, ns_gen_mer1;

				convert_gen_mer_cg_long(cg_gen_mer1, pos1,  read_len1, ta_size1);

				convert_gen_mer_ns_long(ns_gen_mer1, pos1,  read_len1, ta_size1);

				for(j = 0; j < ta_size1; j++)
					ta_check[j] |= ((~(ta_read1[j] | c_read1[j])) & cg_gen_mer1[j]) | ns_gen_mer1[j];
				

					for(j = 0; j < ta_size1; j++)
						ta_check[j] |= c_read1[j] ^ (c_read1[j] & cg_gen_mer1[j]);
					
				mism1 = 0;
				for(j = 0; j < ta_size1; j++)
					mism1 += count_mism(ta_check[j]);
			
				if(mism1 > mism_thresh)
					cg_mism1 = true;
			}//if passed ta=check
			else{
				cg_mism1 = true;
			}//else

			if(cg_mism1 == false ){
				//T->C and A->G not in the same read check
				//read1 check
				char astrand = '+';
				if(mates_pairs[read_id].mark_unique == 0){
					mates_pairs[read_id].mark_unique = 1;
					mates_pairs[read_id].pos1 = pos1;
					mates_pairs[read_id].mism1 = mism1;
					mates_pairs[read_id].strand = '+';//by the first read
				}//if
				else if(mates_pairs[read_id].pos1 != pos1 || 
					mates_pairs[read_id].strand != astrand){
						if(mism1 < mates_pairs[read_id].mism1){
						//unique
							mates_pairs[read_id].mark_unique = 1;
							mates_pairs[read_id].pos1 = pos1;
							mates_pairs[read_id].strand = '+';//by the first read
							mates_pairs[read_id].mism1 = mism1;
													
						}
						else if(mism1 == mates_pairs[read_id].mism1){
							mates_pairs[read_id].mark_unique++;
							if(mism1 == 0 )
								pair_ambiguous = true;
						}
				}//else if
			}//both reads passed cg-check


			if(pair_ambiguous == true)					
				break;

		start1++;

	}//for forward read's list




	if(pair_ambiguous == false){
		//take rev-complement of the first read

		//rev-complement hashing on the first len-mer bases
		ta_rc_read1.resize(ta_size1);
		c_rc_read1.resize(ta_size1);

		int stop = ta_size1 - 1;
		j= stop;
		for(i = 0; i < stop; i++){
			ta_rc_read1[j] = ~Reverse(ta_read1[i], SIZE);
			c_rc_read1[j] =  Reverse(c_read1[i], SIZE);//shows positions of Gs in reverse strand
			j--;

		}
		ta_rc_read1[0] = MASK_READ1 & (~Reverse(ta_read1[stop], len1));
		c_rc_read1[0] =  MASK_READ1 & Reverse(c_read1[stop], len1);//shows positions of Gs in reverse strand

		ta_lmer1 = (ta_rc_read1[stop] & MASK_LMER) & seeds[t][rs];
		cg_lmer1 = (c_rc_read1[stop] & MASK_LMER) & seeds[t][rs];


		if(ta_size1 > 1){
			int residue = SIZE - len1;
			if(len1 < SIZE){
				for(i = 0; i < ta_size1 -1; i++){
					ta_rc_read1[i] = (ta_rc_read1[i] << residue) | (ta_rc_read1[i+1] >> len1);
					c_rc_read1[i] = (c_rc_read1[i] << residue) | (c_rc_read1[i+1] >> len1);

				}//for i
				unsigned long int MASK1 = (one << len1) - 1 ;
				ta_rc_read1[stop] = ta_rc_read1[stop] & MASK1;
				c_rc_read1[stop] = c_rc_read1[stop] & MASK1;
			}//if len1 < SIZE
		
		}//ta size1


		ta_rc1 = hash_aread(ta_lmer1, cg_lmer1, true);

		//process rev-complemnet of read1 and forward of read2

		start1 = 0;
		end1 = pos_sizes[ta_rc1];


		while(start1 < end1){


			pos1 =  pos_index[ta_rc1][start1];

			if(dif1 <= pos1){
				pos1 -= dif1;

				vector<unsigned long int>  ta_gen_mer1;

				bool cg_mism1 = false;
				vector<unsigned long int>  ta_check;
				convert_gen_mer_ta_long(ta_gen_mer1, pos1,  read_len1, ta_size1);
				mism1 = count_mism_long(ta_rc_read1, ta_gen_mer1, ta_check, ta_size1);

				if(mism1 <= mism_thresh){
					vector<unsigned long int> cg_gen_mer1, ns_gen_mer1;
					convert_gen_mer_cg_long( cg_gen_mer1, pos1,  read_len1, ta_size1);
					convert_gen_mer_ns_long(ns_gen_mer1, pos1,  read_len1, ta_size1);

					for(j = 0; j < ta_size1; j++)
ta_check[j] |= ((ta_rc_read1[j] &(~c_rc_read1[j])) & cg_gen_mer1[j]) | ns_gen_mer1[j];
						
						for(j = 0; j < ta_size1; j++)
							ta_check[j] |= c_rc_read1[j] ^ (c_rc_read1[j] & cg_gen_mer1[j]);
					mism1 = 0;
					for(j = 0; j < ta_size1; j++)
						mism1 += count_mism(ta_check[j]);

					if(mism1 > mism_thresh)
						cg_mism1 = true;
				}//if second ta-check passed
				else
					cg_mism1 = true;

				if(cg_mism1 == false ){
					char astrand = '-';
					if(mates_pairs[read_id].mark_unique == 0){
						mates_pairs[read_id].mark_unique = 1;
						mates_pairs[read_id].pos1 = pos1;
						mates_pairs[read_id].mism1 = mism1;
						mates_pairs[read_id].strand = '-';//by the first read
					}//if
					else if(mates_pairs[read_id].pos1 != pos1 || 
						mates_pairs[read_id].strand != astrand){
						if(mism1 < mates_pairs[read_id].mism1){
						//unique
							mates_pairs[read_id].mark_unique = 1;
							mates_pairs[read_id].pos1 = pos1;
							mates_pairs[read_id].strand = '-';//by the first read
							mates_pairs[read_id].mism1 = mism1;
						}
						else if(mism1 == mates_pairs[read_id].mism1){
							mates_pairs[read_id].mark_unique++;
							if(mism1 == 0)
								pair_ambiguous = true;
							}
						}//else if
					}//both reads passed cg-check

					if(pair_ambiguous == true)
						break;
			}
		start1++;

	}//while
	}//if pair is not ambiguous



}//pair_reads_singles_long new




void Anchor::pair_reads_singlesA_long(unsigned long int read_id, string aread1,
				long int &contain_Ns, Stat &astat, short int read_len1, int t, int rs )
{


	unsigned long int one = 1;

	long int asize,  pos1, dist, i, j;
	string strand;
	unsigned short int mism1;
	unsigned long int ta_lmer1, cg_lmer1;
	string line;

	short int len1 = read_len1 & MASK_SIZE;//= (read_len1 mod SIZE) 

	if(len1 == 0)
		len1 = SIZE;

	unsigned long int MASK_READ1, MASK_READ2;
	if(len1 < SIZE){
		MASK_READ1 = ( 1 << len1) - 1;
		if(len1 > 31){
			dist = len1 - 32;
			MASK_READ1 = (((one << dist) - 1) << 32) | 0x00000000FFFFFFFF;
		}

	}
	else
		MASK_READ1 = 0xFFFFFFFFFFFFFFFF;


	unsigned long int ta_rc1=0, ind1 = 0;//hash-indices 
	vector<unsigned long int> c_rc_read1, c_read1;//cg-representation
	vector<unsigned long int>  ta_read1, ta_rc_read1;//ta-representation

	bool pair_ambiguous = false;
	//ta-reads
	convert_lmer_ta_long(aread1, ta_read1, read_len1);

	int ta_size1 = ta_read1.size();


	int chr_id;
	
	//c-reads
	convert_lmer_c_long(aread1, c_read1,  read_len1);

  //rev-complement of c-reads

	int dif1 = read_len1 - len_mer;

	//we need the first len_mer bases for each read
	if(ta_size1 > 1){

		if(dif_len_mer > 0){
			ta_lmer1 = (ta_read1[0] >> dif_len_mer) & seeds[t][rs];
			cg_lmer1 = (c_read1[0] >> dif_len_mer) & seeds[t][rs];
		}
		else{
			ta_lmer1 = ta_read1[0]  & seeds[t][rs];
			cg_lmer1 = c_read1[0]  & seeds[t][rs];
		
		}//else

	}//if
	else{
		ta_lmer1 = (ta_read1[0] >> dif1) & seeds[t][rs];//dif1 < 64
		cg_lmer1 = (c_read1[0] >> dif1) & seeds[t][rs];
	
	}

	ind1 = hash_aread(ta_lmer1, cg_lmer1, true);

	if(ind1 >= hash_table_size){
cout << "ind1 is out of boundary" << endl;
exit(0);
	}
	//rev-complement of g-reads


	//process forward read1 and rev read2
	//start and end of list of genome positions for forw1
	long int start1 = 0;
	long int end1 = pos_sizes[ind1];

	//start and end of list of genome positions for rev-compl2

	while(start1 < end1){

		pos1 = pos_index[ind1][start1];

			//cg check
			bool cg_mism1 = false;
			vector<unsigned long int> ta_gen_mer1;

			convert_gen_mer_ta_long(ta_gen_mer1, pos1, read_len1, ta_size1);


			vector<unsigned long int> ta_check;
			mism1 = count_mism_long(ta_read1, ta_gen_mer1, ta_check, ta_size1);


			if(mism1 <= mism_thresh){

				vector<unsigned long int> cg_gen_mer1, ns_gen_mer1;

				convert_gen_mer_cg_long(cg_gen_mer1, pos1,  read_len1, ta_size1);

				convert_gen_mer_ns_long(ns_gen_mer1, pos1, read_len1, ta_size1);

				for(j = 0; j < ta_size1; j++)
ta_check[j] |= ((ta_read1[j] & (~c_read1[j])) & cg_gen_mer1[j]) | ns_gen_mer1[j];

				

					for(j = 0; j < ta_size1; j++)
						ta_check[j] |= c_read1[j] ^ (c_read1[j] & cg_gen_mer1[j]);
					
				mism1 = 0;
				for(j = 0; j < ta_size1; j++)
					mism1 += count_mism(ta_check[j]);
			
				if(mism1 > mism_thresh)
					cg_mism1 = true;
			}//if passed ta=check
			else{
				cg_mism1 = true;
			}//else

			if(cg_mism1 == false ){
				//T->C and A->G not in the same read check
				//read1 check
				char astrand = '+';
				if(mates_pairs[read_id].mark_unique == 0){
					mates_pairs[read_id].mark_unique = 1;
					mates_pairs[read_id].pos1 = pos1;
					mates_pairs[read_id].mism1 = mism1;
					mates_pairs[read_id].strand = '+';//by the first read
				}//if
				else if(mates_pairs[read_id].pos1 != pos1 || 
					mates_pairs[read_id].strand != astrand){
						if(mism1 < mates_pairs[read_id].mism1){
						//unique
							mates_pairs[read_id].mark_unique = 1;
							mates_pairs[read_id].pos1 = pos1;
							mates_pairs[read_id].strand = '+';//by the first read
							mates_pairs[read_id].mism1 = mism1;
													
						}
						else if(mism1 == mates_pairs[read_id].mism1){
							mates_pairs[read_id].mark_unique++;
							if(mism1 == 0 )
								pair_ambiguous = true;
						}
				}//else if
			}//both reads passed cg-check


			if(pair_ambiguous == true)					
				break;

		start1++;

	}//for forward read's list




	if(pair_ambiguous == false){
		//take rev-complement of the first read



		//rev-complement hashing on the first len-mer bases
		ta_rc_read1.resize(ta_size1);
		c_rc_read1.resize(ta_size1);

		int stop = ta_size1 - 1;
		j= stop;
		for(i = 0; i < stop; i++){
			ta_rc_read1[j] = ~Reverse(ta_read1[i], SIZE);
			c_rc_read1[j] =  Reverse(c_read1[i], SIZE);//shows positions of Gs in reverse strand
			j--;
		}
		ta_rc_read1[0] = MASK_READ1 & (~Reverse(ta_read1[stop], len1));
		c_rc_read1[0] =  MASK_READ1 & Reverse(c_read1[stop], len1);//shows positions of Gs in reverse strand

		ta_lmer1 = (ta_rc_read1[stop] & MASK_LMER) & seeds[t][rs];
		cg_lmer1 = (c_rc_read1[stop] & MASK_LMER) & seeds[t][rs];

		if(ta_size1 > 1){
			int residue = SIZE - len1;
			if(len1 < SIZE){
				for(i = 0; i < ta_size1 -1; i++){
					ta_rc_read1[i] = (ta_rc_read1[i] << residue) | (ta_rc_read1[i+1] >> len1);
					c_rc_read1[i] = (c_rc_read1[i] << residue) | (c_rc_read1[i+1] >> len1);
				}//for i
				unsigned long int MASK1 = (one << len1) - 1 ;
				ta_rc_read1[stop] = ta_rc_read1[stop] & MASK1;
				c_rc_read1[stop] = c_rc_read1[stop] & MASK1;
			}//if len1 < SIZE
		
		}//ta size1


		ta_rc1 = hash_aread(ta_lmer1, cg_lmer1, false);

		//process rev-complemnet of read1 and forward of read2

		start1 = 0;
		end1 = pos_sizes[ta_rc1];


		while(start1 < end1){


			pos1 =  pos_index[ta_rc1][start1];
			if(dif1 <= pos1 ){
				pos1 -= dif1;

				vector<unsigned long int>  ta_gen_mer1;

				bool cg_mism1 = false;
				vector<unsigned long int>  ta_check;
				convert_gen_mer_ta_long(ta_gen_mer1, pos1, read_len1, ta_size1);
				mism1 = count_mism_long(ta_rc_read1, ta_gen_mer1, ta_check, ta_size1);

				if(mism1 <= mism_thresh){
					vector<unsigned long int> cg_gen_mer1,ns_gen_mer1;
					convert_gen_mer_cg_long( cg_gen_mer1, pos1,  read_len1, ta_size1);
					convert_gen_mer_ns_long(ns_gen_mer1, pos1,  read_len1, ta_size1);

					for(j = 0; j < ta_size1; j++)
	ta_check[j] |= ((~(ta_rc_read1[j] | c_rc_read1[j])) & cg_gen_mer1[j]) | ns_gen_mer1[j];


						
						for(j = 0; j < ta_size1; j++)
							ta_check[j] |= c_rc_read1[j] ^ (c_rc_read1[j] & cg_gen_mer1[j]);
					mism1 = 0;
					for(j = 0; j < ta_size1; j++)
						mism1 += count_mism(ta_check[j]);

					if(mism1 > mism_thresh)
						cg_mism1 = true;
				}//if second ta-check passed
				else
					cg_mism1 = true;

				if(cg_mism1 == false ){
					char astrand = '-';
					if(mates_pairs[read_id].mark_unique == 0){
						mates_pairs[read_id].mark_unique = 1;
						mates_pairs[read_id].pos1 = pos1;
						mates_pairs[read_id].mism1 = mism1;
						mates_pairs[read_id].strand = '-';//by the first read
					}//if
					else if(mates_pairs[read_id].pos1 != pos1 || 
						mates_pairs[read_id].strand != astrand){
						if(mism1 < mates_pairs[read_id].mism1){
						//unique
							mates_pairs[read_id].mark_unique = 1;
							mates_pairs[read_id].pos1 = pos1;
							mates_pairs[read_id].strand = '-';//by the first read
							mates_pairs[read_id].mism1 = mism1;
						}
						else if(mism1 == mates_pairs[read_id].mism1){
							mates_pairs[read_id].mark_unique++;
							if(mism1 == 0)
								pair_ambiguous = true;
							}
						}//else if
					}//both reads passed cg-check

					if(pair_ambiguous == true)
						break;

		}//if read2 fits into a ref
		start1++;

	}//while
	}//if pair is not ambiguous



}//pair_reads_singles_long new


void Anchor::pair_reads_normal_singles_long(unsigned long int read_id, string aread1,
				long int &contain_Ns, Stat &astat, short int read_len1, int t, int rs )
{



	unsigned long int one = 1;

	long int asize,   pos1, dist, i, j;
	string strand;
	unsigned short int mism1;
	unsigned long int ta_lmer1, cg_lmer1;
	string line;

	short int len1 = read_len1 & MASK_SIZE;//= (read_len1 mod SIZE) 

	if(len1 == 0)
		len1 = SIZE;

	unsigned long int MASK_READ1, MASK_READ2;
	if(len1 < SIZE){
		MASK_READ1 = ( 1 << len1) - 1;
		if(len1 > 31){
			dist = len1 - 32;
			MASK_READ1 = (((one << dist) - 1) << 32) | 0x00000000FFFFFFFF;
		}

	}
	else
		MASK_READ1 = 0xFFFFFFFFFFFFFFFF;


	unsigned long int ta_rc1=0, ind1 = 0;//hash-indices 
	vector<unsigned long int> c_rc_read1, c_read1;//cg-representation
	vector<unsigned long int>  ta_read1, ta_rc_read1;//ta-representation

	bool pair_ambiguous = false;
	//ta-reads
	convert_lmer_ta_long(aread1, ta_read1, read_len1);

	int ta_size1 = ta_read1.size();


	int chr_id;
	
	//c-reads
	convert_lmer_c_long(aread1, c_read1,  read_len1);

  //rev-complement of c-reads

	int dif1 = read_len1 - len_mer;

	//we need the first len_mer bases for each read
	if(ta_size1 > 1){

		if(dif_len_mer > 0){
			ta_lmer1 = (ta_read1[0] >> dif_len_mer) & seeds[t][rs];
			cg_lmer1 = (c_read1[0] >> dif_len_mer) & seeds[t][rs];
		}
		else{
			ta_lmer1 = ta_read1[0]  & seeds[t][rs];
			cg_lmer1 = c_read1[0]  & seeds[t][rs];
		
		}//else

	}//if
	else{
		ta_lmer1 = (ta_read1[0] >> dif1) & seeds[t][rs];//dif1 < 64
		cg_lmer1 = (c_read1[0] >> dif1) & seeds[t][rs];
	
	}

	ind1 = hash2c(ta_lmer1,  cg_lmer1);

	//rev-complement of g-reads


	//process forward read1 and rev read2
	//start and end of list of genome positions for forw1
	long int start1 = 0;
	long int end1 = pos_sizes[ind1];

	//start and end of list of genome positions for rev-compl2

	while(start1 < end1){

		pos1 = pos_index[ind1][start1];

			//cg check
			bool cg_mism1 = false;
			vector<unsigned long int> ta_gen_mer1;

			convert_gen_mer_ta_long(ta_gen_mer1, pos1, read_len1, ta_size1);


			vector<unsigned long int> ta_check;
			mism1 = count_mism_long(ta_read1, ta_gen_mer1, ta_check, ta_size1);


			if(mism1 <= mism_thresh){

				vector<unsigned long int> cg_gen_mer1, ns_gen_mer1;

				convert_gen_mer_cg_long(cg_gen_mer1, pos1,  read_len1, ta_size1);

				convert_gen_mer_ns_long(ns_gen_mer1, pos1, read_len1, ta_size1);

				for(j = 0; j < ta_size1; j++)
					ta_check[j] |=  ns_gen_mer1[j];
				

					for(j = 0; j < ta_size1; j++)
						ta_check[j] |= c_read1[j] ^ cg_gen_mer1[j];
					
				mism1 = 0;
				for(j = 0; j < ta_size1; j++)
					mism1 += count_mism(ta_check[j]);
			
				if(mism1 > mism_thresh)
					cg_mism1 = true;
			}//if passed ta=check
			else{
				cg_mism1 = true;
			}//else

			if(cg_mism1 == false ){
				//T->C and A->G not in the same read check
				//read1 check
				char astrand = '+';
				if(mates_pairs[read_id].mark_unique == 0){
					mates_pairs[read_id].mark_unique = 1;
					mates_pairs[read_id].pos1 = pos1;
					mates_pairs[read_id].mism1 = mism1;
					mates_pairs[read_id].strand = '+';//by the first read
				}//if
				else if(mates_pairs[read_id].pos1 != pos1 || 
					mates_pairs[read_id].strand != astrand){
						if(mism1 < mates_pairs[read_id].mism1){
						//unique
							mates_pairs[read_id].mark_unique = 1;
							mates_pairs[read_id].pos1 = pos1;
							mates_pairs[read_id].strand = '+';//by the first read
							mates_pairs[read_id].mism1 = mism1;
													
						}
						else if(mism1 == mates_pairs[read_id].mism1){
							mates_pairs[read_id].mark_unique++;
							if(mism1 == 0 )
								pair_ambiguous = true;
						}
				}//else if
			}//both reads passed cg-check


			if(pair_ambiguous == true)					
				break;

		start1++;

	}//for forward read's list




	if(pair_ambiguous == false){
		//take rev-complement of the first read

		//rev-complement hashing on the first len-mer bases
		ta_rc_read1.resize(ta_size1);
		c_rc_read1.resize(ta_size1);

		int stop = ta_size1 - 1;
		j= stop;
		for(i = 0; i < stop; i++){
			ta_rc_read1[j] = ~Reverse(ta_read1[i], SIZE);
			c_rc_read1[j] =  Reverse(c_read1[i], SIZE);//shows positions of Gs in reverse strand
			j--;

		}
		ta_rc_read1[0] = MASK_READ1 & (~Reverse(ta_read1[stop], len1));
		c_rc_read1[0] =  MASK_READ1 & Reverse(c_read1[stop], len1);//shows positions of Gs in reverse strand

		ta_lmer1 = (ta_rc_read1[stop] & MASK_LMER) & seeds[t][rs];
		cg_lmer1 = (c_rc_read1[stop] & MASK_LMER) & seeds[t][rs];


		if(ta_size1 > 1){
			int residue = SIZE - len1;
			if(len1 < SIZE){
				for(i = 0; i < ta_size1 -1; i++){
					ta_rc_read1[i] = (ta_rc_read1[i] << residue) | (ta_rc_read1[i+1] >> len1);
					c_rc_read1[i] = (c_rc_read1[i] << residue) | (c_rc_read1[i+1] >> len1);

				}//for i
				unsigned long int MASK1 = (one << len1) - 1 ;
				ta_rc_read1[stop] = ta_rc_read1[stop] & MASK1;
				c_rc_read1[stop] = c_rc_read1[stop] & MASK1;
			}//if len1 < SIZE
		
		}//ta size1

		ta_rc1 = hash2c(ta_lmer1, cg_lmer1);


		//process rev-complemnet of read1 and forward of read2

		start1 = 0;
		end1 = pos_sizes[ta_rc1];


		while(start1 < end1){


			pos1 =  pos_index[ta_rc1][start1];

			if(dif1 <= pos1){
				pos1 -= dif1;

				vector<unsigned long int>  ta_gen_mer1;

				bool cg_mism1 = false;
				vector<unsigned long int>  ta_check;
				convert_gen_mer_ta_long(ta_gen_mer1, pos1, read_len1, ta_size1);
				mism1 = count_mism_long(ta_rc_read1, ta_gen_mer1, ta_check, ta_size1);

				if(mism1 <= mism_thresh){
					vector<unsigned long int> cg_gen_mer1, ns_gen_mer1;
					convert_gen_mer_cg_long( cg_gen_mer1, pos1, read_len1, ta_size1);
					convert_gen_mer_ns_long(ns_gen_mer1, pos1,  read_len1, ta_size1);

					for(j = 0; j < ta_size1; j++)
						ta_check[j] |= ns_gen_mer1[j];
						
						for(j = 0; j < ta_size1; j++)
							ta_check[j] |= c_rc_read1[j] ^  cg_gen_mer1[j];
					mism1 = 0;
					for(j = 0; j < ta_size1; j++)
						mism1 += count_mism(ta_check[j]);

					if(mism1 > mism_thresh)
						cg_mism1 = true;
				}//if second ta-check passed
				else
					cg_mism1 = true;

				if(cg_mism1 == false ){
					char astrand = '-';
					if(mates_pairs[read_id].mark_unique == 0){
						mates_pairs[read_id].mark_unique = 1;
						mates_pairs[read_id].pos1 = pos1;
						mates_pairs[read_id].mism1 = mism1;
						mates_pairs[read_id].strand = '-';//by the first read
					}//if
					else if(mates_pairs[read_id].pos1 != pos1 || 
						mates_pairs[read_id].strand != astrand){
						if(mism1 < mates_pairs[read_id].mism1){
						//unique
							mates_pairs[read_id].mark_unique = 1;
							mates_pairs[read_id].pos1 = pos1;
							mates_pairs[read_id].strand = '-';//by the first read
							mates_pairs[read_id].mism1 = mism1;
						}
						else if(mism1 == mates_pairs[read_id].mism1){
							mates_pairs[read_id].mark_unique++;
							if(mism1 == 0)
								pair_ambiguous = true;
							}
						}//else if
					}//both reads passed cg-check

					if(pair_ambiguous == true)
						break;

		}//if read2 fits into a ref
		start1++;

	}//while
	}//if pair is not ambiguous



}//pair_reads_singles_long new


/////////////////////////////// END LONG READS///////////////////////////////////////


void Anchor::pair_reads_normal(unsigned long int read_id, string aread1, string aread2,
				long int &contain_Ns, Stat &astat, short int read_len1, short int read_len2, int t, int rs )
{

	unsigned long int one = 1;

	long int asize, pos1, pos2, dist, i;
	string strand, strand2;

	string line;


	unsigned long int MASK_READ1 = ( 1 << read_len1) - 1;
	unsigned long int MASK_READ2 = (1 << read_len2) -1;


	if(read_len1 > 31){
		dist = read_len1 - 32;
		MASK_READ1 = (((one << dist) - 1) << 32) | 0x00000000FFFFFFFF;
	}
	if(read_len2 > 31){
		dist = read_len2 - 32;
		MASK_READ2 = (((one << dist) - 1) << 32) | 0x00000000FFFFFFFF;
		
	}
	unsigned short int mism1, mism2;
	unsigned long int ta_rc1=0, ta_rc2=0, ind1 = 0, ind2 = 0;//hash-indices 
	unsigned long int c_rc_read1, c_rc_read2, c_read1, c_read2;
	unsigned long int  ta_read1, ta_read2, ta_rc_read1, ta_rc_read2, ns=0, ns2=0;
	unsigned long int ta_gen_mer1 = 0, cg_gen_mer1 = 0, ta_gen_mer2 =0, cg_gen_mer2 = 0;
	unsigned long int ta_lmer1, ta_lmer2, cg_lmer1, cg_lmer2, ta_check, ta_check2;
	
	bool pair_ambiguous = false;
	//ta-reads
	convert_lmer_ta(aread1, ta_read1, ns, read_len1);
	convert_lmer_ta(aread2, ta_read2, ns2, read_len2);

	
	//c-reads
	convert_lmer_c(aread1, c_read1, read_len1);
	convert_lmer_c(aread2, c_read2, read_len2);

	int dif1 = read_len1 - len_mer;
	int dif2 = read_len2 - len_mer;

	ta_lmer1 = (ta_read1 >> dif1) & seeds[t][rs];
	cg_lmer1 = (c_read1 >> dif1) & seeds[t][rs];

	ta_lmer2 = (ta_read2 >> dif2) & seeds2[t][rs];
	cg_lmer2 = (c_read2 >> dif2) & seeds2[t][rs];

	ind1 = hash2c(ta_lmer1, cg_lmer1);
	ind2 = hash2c(ta_lmer2, cg_lmer2);


  //rev-complement of c-reads


	//rev-complement of g-reads


	ta_rc_read2 = MASK_READ2 & (~Reverse(ta_read2, read_len2));
	c_rc_read2 =  Reverse(c_read2, read_len2);//shows positions of Gs in reverse strand
	ta_lmer2 = (ta_rc_read2 & MASK_LMER) & seeds2[t][rs];
	cg_lmer2 = (c_rc_read2  & MASK_LMER) & seeds2[t][rs];

	ta_rc2 = hash2c(ta_lmer2, cg_lmer2);

	//process forward read1 and rev read2
	//start and end of list of genome positions for forw1
	long int start1 = 0;
	long int end1 = pos_sizes[ind1];

	//start and end of list of genome positions for rev-compl2

	long int start2 = 0;
	long int end2 = pos_sizes[ta_rc2];
	long int cur = start2;

	while(start1 < end1){
		if(start2 >= end2)
			break;

		pos1 = pos_index[ind1][start1];


			cur = start2;
			for(; cur != end2 ; cur++){
				pos2 = pos_index[ta_rc2][cur];
				if(dif2 <= pos2){

				pos2 -= dif2;
				dist = pos2 - pos1;
				if(dist >= min_insert ){
					break;
				}
				}//if dif2 <= pos2
			}//for cur

			start2 = cur;


			for(; cur != end2; cur++){
				pos2 = pos_index[ta_rc2][cur] ;

				if(dif2 <= pos2){

				pos2 -= dif2;
				dist = pos2 - pos1;
				if(dist >= min_insert && dist <= max_insert ){
					//cg check
					bool cg_mism1 = false, cg_mism2 = false;
					convert_gen_mer_ta(ta_gen_mer1, pos1, read_len1);

					unsigned long int ta_check = ta_read1 ^ ta_gen_mer1, ta_check2;
					mism1 = count_mism(ta_check);

					if(mism1 <= mism_thresh){

						convert_gen_mer_ta(ta_gen_mer2, pos2,  read_len2);
						ta_check2 = ta_rc_read2 ^ ta_gen_mer2;
						mism2 = count_mism(ta_check2);

						if(mism2 <= mism_thresh){
							convert_gen_mer_cg(cg_gen_mer1, pos1,  read_len1);
							convert_gen_mer_cg(cg_gen_mer2, pos2,  read_len2);

							unsigned long int ns_gen_mer1 = 0, ns_gen_mer2 = 0;
							convert_gen_mer_ns(ns_gen_mer1, pos1, read_len1);
							convert_gen_mer_ns(ns_gen_mer2, pos2,  read_len2);
							ta_check |= ns_gen_mer1;
							ta_check2 |= ns_gen_mer2;
								unsigned long int cg_check = c_read1 ^ cg_gen_mer1 ;
								ta_check |= cg_check;
								mism1 = count_mism(ta_check);
								if(mism1 > mism_thresh)
									cg_mism1 = true;
									unsigned long int cg_check2 = c_rc_read2 ^  cg_gen_mer2;
									ta_check2 |= cg_check2;
									mism2 = count_mism(ta_check2);
									if(mism2 > mism_thresh)
										cg_mism2 = true;
						}//if passed second ta-check
						else
							cg_mism2 = true;//not passed 2 ta-check
					}//if passed ta=check
					else{
						cg_mism1 = true;
					}//else



					if(cg_mism1 == false && cg_mism2 == false){
							char astrand = '+';
							int min_mism = mism1;
							int second_mism = mism2;
							if(mism2 < mism1){
								min_mism = mism2;
								second_mism = mism1;
							}
							int cur_second = mates_pairs[read_id].min_mism;
							if(cur_second < mates_pairs[read_id].mism1)
								cur_second = mates_pairs[read_id].mism1;
							else if(cur_second < mates_pairs[read_id].mism2)
								cur_second = mates_pairs[read_id].mism2;

							if(mates_pairs[read_id].mark_unique == 0){
								mates_pairs[read_id].mark_unique = 1;
								mates_pairs[read_id].pos1 = pos1;
								mates_pairs[read_id].mism1 = mism1;
								mates_pairs[read_id].pos2 = pos2;
								mates_pairs[read_id].mism2 = mism2;
								mates_pairs[read_id].min_mism = min_mism;
								mates_pairs[read_id].strand = '+';//by the first read
							}//if
							else if(mates_pairs[read_id].pos1 != pos1 || 
								mates_pairs[read_id].pos2 != pos2 || mates_pairs[read_id].strand != astrand){

								if((min_mism < mates_pairs[read_id].min_mism) || 
									(min_mism == mates_pairs[read_id].min_mism && 
									second_mism < cur_second)){
									//unique
									mates_pairs[read_id].mark_unique = 1;
									mates_pairs[read_id].pos1 = pos1;
									mates_pairs[read_id].strand = '+';//by the first read
									mates_pairs[read_id].mism1 = mism1;
									mates_pairs[read_id].pos2 = pos2;
									mates_pairs[read_id].mism2 = mism2;
									mates_pairs[read_id].min_mism = min_mism;
										
								}
								else if(min_mism == mates_pairs[read_id].min_mism &&
									second_mism == cur_second){
									mates_pairs[read_id].mark_unique++;
									if(min_mism == 0 && second_mism == 0)
										pair_ambiguous = true;
								}
							}//else if

						

					}//both reads passed cg-check
				}//if correct insert range
				else if(dist > max_insert)
					break;
				}//if dif2 <= pos2
				if(pair_ambiguous == true){					
					break;
				}
			}//while cur < end2 && cur <= max_insert
			
			start1++;


		if(pair_ambiguous == true){
			break;
		}//if
	}//for forward read's list

	if(pair_ambiguous == false){
		//take rev-complement of the first read
		ta_rc_read1 = MASK_READ1 & (~Reverse(ta_read1, read_len1));
		c_rc_read1 =  Reverse(c_read1, read_len1);//shows positions of Gs in reverse strand

		ta_lmer1 = (ta_rc_read1  & MASK_LMER) & seeds[t][rs];
		cg_lmer1 = (c_rc_read1  & MASK_LMER) & seeds[t][rs];

		ta_rc1 = hash2c(ta_lmer1, cg_lmer1);


		//process rev-complemnet of read1 and forward of read2

		start1 = 0;
		end1 = pos_sizes[ind2];

		start2 = 0;
		end2 = pos_sizes[ta_rc1];

		while(start1 < end1){

			if(start2 >= end2)
				break;

			pos1 =  pos_index[ind2][start1];

				cur = start2;
				for(; cur != end2 ; cur++){
					pos2 =  pos_index[ta_rc1][cur] ;
					if(dif1 <= pos2){
							pos2 -= dif1;
						dist = pos2 - pos1;
						if(dist >= min_insert ){
							break;
						}//if
					}//if dif1 <= pos2
				}//for cur

				start2 = cur;

					for(; cur != end2; cur++){
						pos2 =  pos_index[ta_rc1][cur];
						if(dif1 <= pos2){

							pos2 -= dif1;
						dist = pos2 - pos1;
						if(dist >= min_insert && dist <= max_insert ){
					//cg check
					
							bool cg_mism1 = false, cg_mism2 = false;
							convert_gen_mer_ta(ta_gen_mer2, pos1, read_len2);

					unsigned long int ta_check2 = ta_read2 ^ ta_gen_mer2, ta_check;
					mism2 = count_mism(ta_check2);

					if(mism2 <= mism_thresh){
						convert_gen_mer_ta(ta_gen_mer1, pos2,  read_len1);
						ta_check = ta_rc_read1 ^ ta_gen_mer1;
						mism1 = count_mism(ta_check);

						if(mism1 <= mism_thresh){
							convert_gen_mer_cg( cg_gen_mer2, pos1,  read_len2);
							convert_gen_mer_cg( cg_gen_mer1, pos2,  read_len1);

							unsigned long int ns_gen_mer1 = 0, ns_gen_mer2 = 0;
							convert_gen_mer_ns(ns_gen_mer1, pos2,  read_len1);
							convert_gen_mer_ns(ns_gen_mer2, pos1,  read_len2);
							ta_check |= ns_gen_mer1;
							ta_check2 |= ns_gen_mer2;

								unsigned long int cg_check2 = c_read2 ^ cg_gen_mer2;
								ta_check2 |= cg_check2;
								mism2 = count_mism(ta_check2);
								if(mism2 > mism_thresh)
									cg_mism2 = true;
								unsigned long int cg_check =c_rc_read1 ^ cg_gen_mer1;
								ta_check |= cg_check;
								mism1 = count_mism(ta_check);
								if(mism1 > mism_thresh)
									cg_mism1 = true;

						}//if second ta-check passed
						else
							cg_mism2 = true;
						}//if first ta_check passed
						else{
							cg_mism1 = true;
						}//else


							if(cg_mism1 == false && cg_mism2 == false){

							char astrand = '-';
							int min_mism = mism1;
							int second_mism = mism2;
							if(mism2 < mism1){
								min_mism = mism2;
								second_mism = mism1;
							}
							int cur_second = mates_pairs[read_id].min_mism;
							if(cur_second <  mates_pairs[read_id].mism1)
								cur_second = mates_pairs[read_id].mism1;
							else if(cur_second < mates_pairs[read_id].mism2)
								cur_second = mates_pairs[read_id].mism2;

							if(mates_pairs[read_id].mark_unique == 0){
								mates_pairs[read_id].mark_unique = 1;
								mates_pairs[read_id].pos1 = pos2;
								mates_pairs[read_id].mism1 = mism1;
								mates_pairs[read_id].pos2 = pos1;
								mates_pairs[read_id].mism2 = mism2;
								mates_pairs[read_id].min_mism = min_mism;
								mates_pairs[read_id].strand = '-';//by the first read
							}//if
							else if(mates_pairs[read_id].pos1 != pos2 || 
								mates_pairs[read_id].pos2 != pos1 || mates_pairs[read_id].strand != astrand){

								if((min_mism < mates_pairs[read_id].min_mism) || 
									(min_mism == mates_pairs[read_id].min_mism && 
									second_mism < cur_second)){
									//unique
									mates_pairs[read_id].mark_unique = 1;
									mates_pairs[read_id].pos1 = pos2;
									mates_pairs[read_id].strand = '-';//by the first read
									mates_pairs[read_id].mism1 = mism1;
									mates_pairs[read_id].pos2 = pos1;
									mates_pairs[read_id].mism2 = mism2;
									mates_pairs[read_id].min_mism = min_mism;
										
								}
								else if(min_mism == mates_pairs[read_id].min_mism &&
									second_mism == cur_second){
									mates_pairs[read_id].mark_unique++;
									if(min_mism == 0 && second_mism == 0)
										pair_ambiguous = true;
								}
							}//else if
		

						}//both reads passed cg-check
					}//if correct insert range
					else if(dist > max_insert)
						break;
					}//if dif1 <= pos2
					if(pair_ambiguous == true)
						break;
				}//while cur < end2 && cur <= max_insert
			start1++;

		if(pair_ambiguous == true){
			break;
		}//if
	}//while
	}//if pair is not ambiguous




}//pair_reads_normal

//START SINGLES
void Anchor::pair_reads_singles(unsigned long int read_id, string aread1,
				long int &contain_Ns, Stat &astat, short int read_len1,  int t, int rs )
{


	unsigned long int one = 1, mism1, ta_lmer1, cg_lmer1;

	long int asize,   pos1, dist, i;
	string strand;

	string line;
	unsigned long int MASK_READ1 = ( 1 << read_len1) - 1;


	if(read_len1 > 31){
		dist = read_len1 - 32;
		MASK_READ1 = (((one << dist) - 1) << 32) | 0x00000000FFFFFFFF;
	}

	unsigned long int ta_rc1=0,  ind1 = 0;//hash-indices 
	unsigned long int c_rc_read1,  c_read1;
	unsigned long int  ta_read1,  ta_rc_read1,  ns=0;
	unsigned long int ta_gen_mer1 = 0, cg_gen_mer1 = 0;

	bool pair_ambiguous = false;
	//ta-reads
	convert_lmer_ta(aread1, ta_read1, ns, read_len1);


	int chr_id, chr_id2;
	
	//c-reads
	convert_lmer_c(aread1, c_read1, read_len1);
	
	int dif1 = read_len1 - len_mer;
	
	ta_lmer1 = (ta_read1 >> dif1) & seeds[t][rs];
	cg_lmer1 = (c_read1 >> dif1) & seeds[t][rs];

	ind1 = hash_aread(ta_lmer1, cg_lmer1, false);//hash2m(ta_read1,  MASK_SEED, MASK_HALF);//ta_read1 >> dif_read_seed1;

  //rev-complement of c-reads


	unsigned long int tc_mism1, ag_mism1;

	//process forward read1 and rev read2
	//start and end of list of genome positions for forw1
	long int start1 =0;
	long int end1 = pos_sizes[ind1];


	while(start1 < end1){

		pos1 = pos_index[ind1][start1];

		bool cg_mism1 = false;
		convert_gen_mer_ta(ta_gen_mer1, pos1, read_len1);

		unsigned long int ta_check = ta_read1 ^ ta_gen_mer1;
		mism1 = count_mism(ta_check);

		if(mism1 <= mism_thresh){
			unsigned long int ns_gen_mer1 = 0;
			convert_gen_mer_ns(ns_gen_mer1, pos1, read_len1);
			ta_check |= ns_gen_mer1;

			convert_gen_mer_cg(cg_gen_mer1, pos1, read_len1);
			ag_mism1 = (~(ta_read1 | c_read1)) & cg_gen_mer1;
			ta_check |= ag_mism1;

				unsigned long int cg_check = c_read1 ^ (c_read1 & cg_gen_mer1);
				//count total mism
				ta_check |= cg_check;
				
			mism1 = count_mism(ta_check);
			if(mism1 > mism_thresh)
				cg_mism1 = true;


		}//if passed ta=check
		else{
			cg_mism1 = true;
		}//else
		if(cg_mism1 == false){
			char astrand = '+';
			//T->C and A->G not in the same read check
				if(mates_pairs[read_id].mark_unique == 0){
					mates_pairs[read_id].mark_unique++;
					mates_pairs[read_id].pos1 = pos1;
					mates_pairs[read_id].strand = '+';//by the first read
					mates_pairs[read_id].mism1 = mism1;
				}//if
				else if(mates_pairs[read_id].pos1 != pos1 || mates_pairs[read_id].strand != astrand){
					if(mism1 < mates_pairs[read_id].mism1){
						//unique
						mates_pairs[read_id].mark_unique = 1;
						mates_pairs[read_id].pos1 = pos1;
						mates_pairs[read_id].strand = '+';//by the first read
						mates_pairs[read_id].mism1 = mism1;
						
					}
					else if(mism1 == mates_pairs[read_id].mism1){
						mates_pairs[read_id].mark_unique++;
						if(mism1 == 0)
							pair_ambiguous = true;
					}
				}//else ambiguous pair							

		}//both reads passed cg-check
		if(pair_ambiguous == true)			
			break;
		start1++;
	}//for forward read's list

	if(pair_ambiguous == false){

		//take rev-complement of the first read
		ta_rc_read1 = MASK_READ1 & (~Reverse(ta_read1, read_len1));

		c_rc_read1 =  Reverse(c_read1, read_len1);//shows positions of Gs in reverse strand

		ta_lmer1 = (ta_rc_read1  & MASK_LMER) & seeds[t][rs];
		cg_lmer1 = (c_rc_read1  & MASK_LMER) & seeds[t][rs];

		ta_rc1  = hash_aread(ta_lmer1, cg_lmer1, true);// hash2m(ta_rc_read1,  MASK_SEED, MASK_HALF);

		//process rev-complemnet of read1 
		start1 = 0;
		end1 = pos_sizes[ta_rc1];

		while(start1 < end1){

			pos1 = pos_index[ta_rc1][start1] ;

			if(dif1 <= pos1){
				pos1 -= dif1;

			bool cg_mism1 = false;

			convert_gen_mer_ta(ta_gen_mer1, pos1, read_len1);
			unsigned long int ta_check = ta_rc_read1 ^ ta_gen_mer1;
			mism1 = count_mism(ta_check);
			if(mism1 <= mism_thresh){
				unsigned long int ns_gen_mer1 = 0;
				convert_gen_mer_ns(ns_gen_mer1, pos1,  read_len1);
				ta_check |= ns_gen_mer1;

				convert_gen_mer_cg(cg_gen_mer1, pos1,  read_len1);
				tc_mism1 = (ta_rc_read1 &(~c_rc_read1)) & cg_gen_mer1;
				ta_check |= tc_mism1;
					unsigned long int cg_check = c_rc_read1 ^ (c_rc_read1 & cg_gen_mer1);
					ta_check |= cg_check;
				mism1 = count_mism(ta_check);
				if(mism1 > mism_thresh)
					cg_mism1 = true;
	
			}//if first ta_check passed
			else{
				cg_mism1 = true;
			}
			if(cg_mism1 == false){
				//T->C and A->G not in the same read check
				char astrand = '-';
				if(mates_pairs[read_id].mark_unique == 0){
					mates_pairs[read_id].mark_unique++;
					mates_pairs[read_id].pos1 = pos1;
					mates_pairs[read_id].strand = '-';//by the first read
					mates_pairs[read_id].mism1 = mism1;
				}//if
				else if(mates_pairs[read_id].pos1 != pos1 || mates_pairs[read_id].strand != astrand){
					if(mism1 < mates_pairs[read_id].mism1){
						//unique
						mates_pairs[read_id].mark_unique = 1;
						mates_pairs[read_id].pos1 = pos1;
						mates_pairs[read_id].strand = '-';//by the first read
						mates_pairs[read_id].mism1 = mism1;
						
					}
					else if(mism1 == mates_pairs[read_id].mism1){
						mates_pairs[read_id].mark_unique++;
						if(mism1 == 0)
							pair_ambiguous = true;
					}
				}//else					

			}//passed cg-check
			}//if dif1 <= pos1
			if(pair_ambiguous == true)
				break;
			start1++;
		}//while

	}//if pair is not ambiguous


}//pair_reads_singles new



///////////////////////////////////////////////////////////
//A-rich strands
void Anchor::pair_reads_singlesA(unsigned long int read_id, string aread1,
				long int &contain_Ns, Stat &astat, short int read_len1, int t, int rs )
{


//cout << read_id << "\t" << aread1 << "\t" << aread2 << endl;
	unsigned long int one = 1;
	short int mism1;
	unsigned long int ta_lmer1, cg_lmer1;
	long int asize,  pos1, dist, i;
	string strand;

	string line;
	unsigned long int MASK_READ1 = ( 1 << read_len1) - 1;

	int dif1 = read_len1 - len_mer;
	if(read_len1 > 31){
		dist = read_len1 - 32;
		MASK_READ1 = (((one << dist) - 1) << 32) | 0x00000000FFFFFFFF;
	}

	unsigned long int ta_rc1=0,  ind1 = 0;//hash-indices 
	unsigned long int c_rc_read1,  c_read1;
	unsigned long int  ta_read1,  ta_rc_read1,  ns=0;
	unsigned long int ta_gen_mer1 = 0, cg_gen_mer1 = 0;

	bool pair_ambiguous = false;
	//ta-reads
	convert_lmer_ta(aread1, ta_read1, ns, read_len1);
	int chr_id, chr_id2;
	
	//c-reads
	convert_lmer_c(aread1, c_read1,  read_len1);

  //rev-complement of c-reads
	ta_lmer1 = (ta_read1 >> dif1) & seeds[t][rs];
	cg_lmer1 = (c_read1 >> dif1) & seeds[t][rs];
	ind1 = hash_aread(ta_lmer1, cg_lmer1, true);

	unsigned long int tc_mism1, ag_mism1;

	//process forward read1 and rev read2
	//start and end of list of genome positions for forw1
	long int start1 =0;
	long int end1 = pos_sizes[ind1];


	while(start1 < end1){

		pos1 = pos_index[ind1][start1];

		bool cg_mism1 = false;
		convert_gen_mer_ta(ta_gen_mer1, pos1,  read_len1);

		unsigned long int ta_check = ta_read1 ^ ta_gen_mer1;
		mism1 = count_mism(ta_check);

		if(mism1 <= mism_thresh){
			convert_gen_mer_cg(cg_gen_mer1, pos1, read_len1);
			tc_mism1 = (ta_read1 & (~c_read1)) & cg_gen_mer1 ;
			ta_check |= tc_mism1;
			unsigned long int ns_gen_mer1 = 0;
			convert_gen_mer_ns(ns_gen_mer1, pos1,  read_len1);
			ta_check |= ns_gen_mer1;

				unsigned long int cg_check = c_read1 ^ (c_read1 & cg_gen_mer1);
				//count total mism
				ta_check |= cg_check;
			mism1 = count_mism(ta_check);
			if(mism1 > mism_thresh)
				cg_mism1 = true;

		}//if passed ta=check
		else{
			cg_mism1 = true;
		}//else
		if(cg_mism1 == false){
				char astrand = '+';
				if(mates_pairs[read_id].mark_unique == 0){
					mates_pairs[read_id].mark_unique++;
					mates_pairs[read_id].pos1 = pos1;
					mates_pairs[read_id].strand = '+';//by the first read
					mates_pairs[read_id].mism1 = mism1;
				}//if
				else if(mates_pairs[read_id].pos1 != pos1 || mates_pairs[read_id].strand != astrand){
					if(mism1 < mates_pairs[read_id].mism1){
						//unique
						mates_pairs[read_id].mark_unique = 1;
						mates_pairs[read_id].pos1 = pos1;
						mates_pairs[read_id].strand = '+';//by the first read
						mates_pairs[read_id].mism1 = mism1;
						
					}
					else if(mism1 == mates_pairs[read_id].mism1){
						mates_pairs[read_id].mark_unique++;
						if(mism1 == 0)
							pair_ambiguous = true;
					}
				}//else ambiguous pair							

		}//both reads passed cg-check
		if(pair_ambiguous == true)			
			break;
		start1++;
	}//for forward read's list

	if(pair_ambiguous == false){

		//take rev-complement of the first read
		ta_rc_read1 = MASK_READ1 & (~Reverse(ta_read1, read_len1));
		c_rc_read1 =  Reverse(c_read1, read_len1);//shows positions of Gs in reverse strand

		ta_lmer1 = (ta_rc_read1  & MASK_LMER) & seeds[t][rs];
		cg_lmer1 = (c_rc_read1  & MASK_LMER) & seeds[t][rs];

		ta_rc1 = hash_aread(ta_lmer1, cg_lmer1, false);


		//process rev-complemnet of read1 
		start1 = 0;
		end1 = pos_sizes[ta_rc1];

		while(start1 < end1){

			pos1 = pos_index[ta_rc1][start1] ;
			bool cg_mism1 = false;

			if(dif1 <= pos1){
				pos1 -= dif1;
			convert_gen_mer_ta(ta_gen_mer1, pos1, read_len1);

			unsigned long int ta_check = ta_rc_read1 ^ ta_gen_mer1;
			mism1 = count_mism(ta_check);
			if(mism1 <= mism_thresh){
				unsigned long int ns_gen_mer1 = 0;
				convert_gen_mer_ns(ns_gen_mer1, pos1,  read_len1);
				ta_check |= ns_gen_mer1;

				convert_gen_mer_cg(cg_gen_mer1, pos1, read_len1);

				ag_mism1 = (~(ta_rc_read1 | c_rc_read1)) & cg_gen_mer1;
				ta_check |= ag_mism1;

					unsigned long int cg_check = c_rc_read1 ^ (c_rc_read1 & cg_gen_mer1);
					ta_check |= cg_check;

				mism1 = count_mism(ta_check);
				if(mism1 > mism_thresh)
					cg_mism1 = true;


				}//if first ta_check passed
			else{
				cg_mism1 = true;
			}
			if(cg_mism1 == false){
				char astrand = '-';

				if(mates_pairs[read_id].mark_unique == 0){
					mates_pairs[read_id].mark_unique++;
					mates_pairs[read_id].pos1 = pos1;
					mates_pairs[read_id].strand = '-';//by the first read
					mates_pairs[read_id].mism1 = mism1;
				}//if
				else if(mates_pairs[read_id].pos1 != pos1 || mates_pairs[read_id].strand != astrand){
					if(mism1 < mates_pairs[read_id].mism1){
						//unique
						mates_pairs[read_id].mark_unique = 1;
						mates_pairs[read_id].pos1 = pos1;
						mates_pairs[read_id].strand = '-';//by the first read
						mates_pairs[read_id].mism1 = mism1;
						
					}
					else if(mism1 == mates_pairs[read_id].mism1){
						mates_pairs[read_id].mark_unique++;
						if(mism1 == 0)
							pair_ambiguous = true;
					}
				}//else							
	
			}//passed cg-check
			}//if dif1 <= pos1
			if(pair_ambiguous == true)
				break;
			start1++;
		}//while

	}//if pair is not ambiguous




}//pair_reads_singles new

////////////////////////////////////////////////////////////
//Normal
void Anchor::pair_reads_normal_singles(unsigned long int read_id, string aread1,
				long int &contain_Ns, Stat &astat, short int read_len1, int t, int rs )
{

//cout << read_id << "\t" << aread1 << "\t" << aread2 << endl;
	unsigned long int one = 1, ta_lmer1, cg_lmer1;

	long int asize,    pos1, dist, i;
	string strand;

	string line;

	unsigned long int MASK_READ1 = ( 1 << read_len1) - 1;


	if(read_len1 > 31){
		dist = read_len1 - 32;
		MASK_READ1 = (((one << dist) - 1) << 32) | 0x00000000FFFFFFFF;
	}

	unsigned long int ta_rc1=0,  ind1 = 0;//hash-indices 
	unsigned long int c_rc_read1,  c_read1;
	unsigned long int  ta_read1,  ta_rc_read1,  ns=0;
	unsigned long int ta_gen_mer1 = 0, cg_gen_mer1 = 0, ta_check;

	bool pair_ambiguous = false;
	//ta-reads
	convert_lmer_ta(aread1, ta_read1, ns, read_len1);
	int chr_id, chr_id2;
	unsigned short int mism1;
	//c-reads
	convert_lmer_c(aread1, c_read1,  read_len1);

	int dif1 = read_len1 - len_mer;
	ta_lmer1 = (ta_read1 >> dif1) & seeds[t][rs];
	cg_lmer1 = (c_read1 >> dif1) & seeds[t][rs];


	ind1 = hash2c(ta_lmer1,  cg_lmer1);



  //rev-complement of c-reads



	//process forward read1 and rev read2
	//start and end of list of genome positions for forw1
	long int start1 =0;
	long int end1 = pos_sizes[ind1];


	while(start1 < end1){

		pos1 = pos_index[ind1][start1];

		bool cg_mism1 = false;
		convert_gen_mer_ta(ta_gen_mer1, pos1,  read_len1);

		unsigned long int ta_check = ta_read1 ^ ta_gen_mer1;
		mism1 = count_mism(ta_check);

		if(mism1 <= mism_thresh){
			unsigned long int ns_gen_mer1 = 0;
			convert_gen_mer_ns(ns_gen_mer1, pos1,  read_len1);
			ta_check |= ns_gen_mer1;


			convert_gen_mer_cg(cg_gen_mer1, pos1,  read_len1);
				unsigned long int cg_check = c_read1 ^ cg_gen_mer1;
				//count total mism
				ta_check |= cg_check;
				mism1 = count_mism(ta_check);

				if(mism1 > mism_thresh)
					cg_mism1 = true;
		}//if passed ta=check
		else{
			cg_mism1 = true;
		}//else
		if(cg_mism1 == false){
			char astrand = '+';
			//T->C and A->G not in the same read check
				if(mates_pairs[read_id].mark_unique == 0){
					mates_pairs[read_id].mark_unique++;
					mates_pairs[read_id].pos1 = pos1;
					mates_pairs[read_id].strand = '+';//by the first read
					mates_pairs[read_id].mism1 = mism1;
				}//if
				else if(mates_pairs[read_id].pos1 != pos1 || mates_pairs[read_id].strand != astrand){
					if(mism1 < mates_pairs[read_id].mism1){
						//unique
						mates_pairs[read_id].mark_unique = 1;
						mates_pairs[read_id].pos1 = pos1;
						mates_pairs[read_id].strand = '+';//by the first read
						mates_pairs[read_id].mism1 = mism1;
						
					}
					else if(mism1 == mates_pairs[read_id].mism1){
						mates_pairs[read_id].mark_unique++;
						if(mism1 == 0)
							pair_ambiguous = true;
					}
				}//else ambiguous pair											

		}//both reads passed cg-check
		if(pair_ambiguous == true)			
			break;
		start1++;
	}//for forward read's list

	if(pair_ambiguous == false){

		//take rev-complement of the first read
		ta_rc_read1 = MASK_READ1 & (~Reverse(ta_read1, read_len1));
		c_rc_read1 =  Reverse(c_read1, read_len1);//shows positions of Gs in reverse strand
		ta_lmer1 = (ta_rc_read1  & MASK_LMER) & seeds[t][rs];
		cg_lmer1 = (c_rc_read1  & MASK_LMER) & seeds[t][rs];

		ta_rc1 = hash2c(ta_lmer1, cg_lmer1);


		//process rev-complemnet of read1 
		start1 = 0;
		end1 = pos_sizes[ta_rc1];

		while(start1 < end1){

			pos1 = pos_index[ta_rc1][start1];
			bool cg_mism1 = false;

			if(dif1 <= pos1){
				pos1 -= dif1;

			convert_gen_mer_ta(ta_gen_mer1, pos1,  read_len1);

			unsigned long int ta_check = ta_rc_read1 ^ ta_gen_mer1;
			mism1 = count_mism(ta_check);
			if(mism1 <= mism_thresh){
				unsigned long int ns_gen_mer1 = 0;
				convert_gen_mer_ns(ns_gen_mer1, pos1, read_len1);
				ta_check |= ns_gen_mer1;

				convert_gen_mer_cg(cg_gen_mer1, pos1,  read_len1);
					unsigned long int cg_check = c_rc_read1 ^ cg_gen_mer1;
					ta_check |= cg_check;
					mism1 = count_mism(ta_check);
					if(mism1 > mism_thresh)
						cg_mism1 = true;
				}//if first ta_check passed
			else{
				cg_mism1 = true;
			}
			if(cg_mism1 == false){
				char astrand = '-';
				if(mates_pairs[read_id].mark_unique == 0){
					mates_pairs[read_id].mark_unique++;
					mates_pairs[read_id].pos1 = pos1;
					mates_pairs[read_id].strand = '-';//by the first read
					mates_pairs[read_id].mism1 = mism1;
				}//if
				else if(mates_pairs[read_id].pos1 != pos1 || mates_pairs[read_id].strand != astrand){
					if(mism1 < mates_pairs[read_id].mism1){
						//unique
						mates_pairs[read_id].mark_unique = 1;
						mates_pairs[read_id].pos1 = pos1;
						mates_pairs[read_id].strand = '-';//by the first read
						mates_pairs[read_id].mism1 = mism1;
						
					}
					else if(mism1 == mates_pairs[read_id].mism1){
						mates_pairs[read_id].mark_unique++;
						if(mism1 == 0)
							pair_ambiguous = true;
					}
				}//else		
			}//passed cg-check
			}//if dif1 <= pos1
			if(pair_ambiguous == true)
				break;
			start1++;
		}//while

	}//if pair is not ambiguous


}//pair_reads_singles new


int Anchor::parse_options(int argc, char* argv[])
{
	int res = 0;

	long int i, j;
	for(i = 1; i < argc; i++){
		if(argv[i][0] != '-')
			return -1;
		switch(argv[i][1]){
		case 'r': genome_names = argv[++i]; break;
		case 's': reads_file1 = argv[++i]; break;
		case '1': reads_file1 = argv[++i]; break;
		case '2': reads_file2 = argv[++i]; break;
		case 'o': output_pairs = argv[++i]; break;
		case 'p': is_pe = "1"; break;
		case 'b': is_bs = "1"; break;
		case 'i': min_insert = atoi(argv[++i]); break;
		case 'a': max_insert = atoi(argv[++i]); break;
		case 'A': A_rich = true; break;
		case 'u': unmapped = true; break;
		case 'm': mism_thresh = atoi(argv[++i]); break;
		case 'M': out_amb = true; break;
		case 'f': first_bits = atoi(argv[++i]); break;
		case 'h': usage(); break;
		default: cerr << "ERROR: Invalid char=" << argv[i][1] << "=this" << endl; usage(); exit(0); break;
		}
	}//for i
	return res;
}//parse_opt()

void Anchor::print(Stat &astat){
	long int i, j;
	astat.mapped_pairs = 0;
	astat.ambiguous = 0;
	astat.invalid = 0;
	long int st1, en1, st2, en2;
	long int orig_start1, orig_start2;

		ofstream out_amb1, out_amb2;
		string output_amb1 = reads_file1 + ".amb";
		string output_amb2 =  reads_file2 + ".amb";


		if(out_amb == true){
			out_amb1.open(output_amb1.c_str(), ios::out);
		
			if(is_pe == "1")
				out_amb2.open(output_amb2.c_str(), ios::out);
	
		}	//if to output ambiguous option

	ifstream in;
	in.open(reads_file1.c_str(), ios::in);
	if(!in){
		cout << "\nERROR: could not open " << reads_file1 << " to print output " << endl;
		exit(0);
	}
	if(is_pe == "1"){

		ifstream in2;
		in2.open(reads_file2.c_str(), ios::in);
		if(!in2){
			cout << "\nERROR: could not open " << reads_file2 << " to print output " << endl;
			exit(0);
		}



		long int read_id = 0;
		string aread1, aread2;
		in >> aread1;
		in2 >> aread2;
		while(!in.eof()){
			in >> st1 >> en1;
			in2 >> st2 >> en2;

			if(mates_pairs[read_id].mark_unique == 1){

				//determine chrom id for this read
				int chrom_id = 0;
				long int asize = size_chrom.size();
				//find the first start of chrom that is greater than pos
				long int pos1 = mates_pairs[read_id].pos1; 
				for(i = 0; i < asize; i++){
					if(pos1 >= starts_size_chrom[i])
						chrom_id++;
					else
						break;

				}//for i
				chrom_id--;

				pos1 = pos1 - starts_size_chrom[chrom_id];
				long int pos2 = mates_pairs[read_id].pos2 - starts_size_chrom[chrom_id];
				astat.mapped_pairs++;
				map<int, string>::iterator afind = chroms_names.find(chrom_id);
				string chr_name = afind->second;
				
				if(mates_pairs[read_id].strand == '+'){
					orig_start1 = pos1 - st1;
					orig_start2 = pos2 - en2;
				}
				else{
					orig_start1 = pos1 - en1;
					orig_start2 = pos2 - st2;
				}
				if(orig_start1 < 0)
					orig_start1 = pos1;
				if(orig_start2 < 0)
					orig_start2 = pos2;

				out_pairs << read_id << "\t" << aread1 << "\t" << aread2 << "\t"
						<< chr_name << "\t" << mates_pairs[read_id].strand <<"\t"
						<< pos1 << "\t" 
						<< pos2 << "\t"
						<< mates_pairs[read_id].mism1 << "\t"
						<< mates_pairs[read_id].mism2 << "\t"
						<< orig_start1 << "\t" << orig_start2 << endl;
			}//if
			else{
				if(mates_pairs[read_id].mark_unique > 1){
					astat.ambiguous++;
					if(out_amb == true){
						out_amb1 << aread1 << endl;
						out_amb2 << aread2 << endl;
					}
				}//if amb
		
				if(unmapped == true && mates_pairs[read_id].mark_unique < 1){
					out_unmapped1 << aread1 << "\t" << st1 << "\t" << en1 << endl;
					out_unmapped2 << aread2 << "\t" << st2 << "\t" << en2 << endl;
				}//if unmap

				if(mates_pairs[read_id].mark_unique < 0)
					astat.invalid++;
			}
			read_id++;
			in >> aread1;
			in2 >> aread2;
		}//while		

		in2.close(); in2.clear();
	}//print pair-ends
	else{

		long int read_id = 0;
		string aread1;
		in >> aread1;
		while(!in.eof()){
			in >> st1 >> en1;
			if(mates_pairs[read_id].mark_unique == 1){
				int chrom_id = 0;
				long int asize = size_chrom.size();
				//find the first start of chrom that is greater than pos
				long int pos1 = mates_pairs[read_id].pos1; 
				for(i = 0; i < asize; i++){
					if(pos1 >= starts_size_chrom[i])
						chrom_id++;
					else
						break;

				}//for i
				chrom_id--;

				pos1 = pos1 - starts_size_chrom[chrom_id];

				astat.mapped_pairs++;
				map<int, string>::iterator afind = chroms_names.find(chrom_id);
				string chr_name = afind->second;
				if(mates_pairs[read_id].strand == '+')
					orig_start1 = pos1 - st1;
				else
					orig_start1 = pos1 - en1;
				if(orig_start1 < 0)
					orig_start1 = pos1;
				out_pairs << read_id << "\t" << aread1 << "\t" 
						<< chr_name << "\t" << mates_pairs[read_id].strand <<"\t"
						<< pos1 << "\t" << mates_pairs[read_id].mism1 << "\t" << orig_start1 << endl;
			}//if
			else if(mates_pairs[read_id].mark_unique > 1){

				astat.ambiguous++;

				if(out_amb == true)
					out_amb1 << aread1 << endl;
			}
			else
			{
				if(unmapped == true)
					out_unmapped1 << aread1 << "\t" << st1 << "\t" << en1 << endl;

				if(mates_pairs[read_id].mark_unique < 0)
					astat.invalid++;
			}
			read_id++;
			in >> aread1;
		}//while
			
	}//print singles
	in.close(); in.clear();

	if(out_amb == true){
		out_amb1.close();
		if(is_pe == "1")
			out_amb2.close();
	}
}//print



int main(int argc, char* argv[])
{
	/**************************************
	This program calculates distribution of different fingerprints in the genome given
	*******************************************/
		/**************************************
		bitset [11000] is translated into a string as 00011 and into long int as 3
		so, least significant bits in the bitset are the leftmost
  *******************************************/
	
	const int SIZE = 64;
	const float alphabet_size = 4;
	long int i, j;
	long int asize , u, y;

	if(argc < 6){
		usage();
		return 0;
	}//if
	

	unsigned long int chrom_id = 0;

	long int total = 0;

	unsigned long int len_mer = 0;
	long int   pos, some;
	string strand, aread, ascore;

	long int id,  pos1, pos2;
	string  aread2, strand2, line;

	ifstream in ;

	Stat astat;

	ifstream in2;
	string read, read2;
	unsigned long int gta, ns, gta2, ns2, one = 1, byte_size = 8;
	int len1 = 0, len2 = 0;

	long int read_id = 0;
	Anchor anc(argc, argv);//reads genome builds ta-, cg- binary representations once

	int t = 0, rs = 0;
	long int st1, en1, st2, en2;
	int seeds_size = anc.gen_seeds.size();

for(t = 0; t < seeds_size; t++){

		anc.find_pos_sizes(t);
		anc.make_pos_index(t);

	int reads_seeds_size = anc.seeds[t].size();
	for(rs = 0; rs < reads_seeds_size; rs++){

	//find the total number of input pairs or single reads
	if(anc.is_pe == "1"){

		one = 1;
		long int thresh_mism = 0;


		int r;

		int chr_id;
		long int contain_Ns = 0;
		in.open(anc.reads_file1.c_str(), ios::in);
		if(!in){
			cout << "can't open " << anc.reads_file1 << " reads_file " << endl;
			exit(0);
		}
		in2.open(anc.reads_file2.c_str(), ios::in);
		if(!in2){
			cout << "can't open " << anc.reads_file2 << " reads_file2" << endl;
			exit(0);
		}
		

	if(anc.is_bs == "1"){
		in >> read;
		in2 >> read2;

		read_id = 0;

		while(!in.eof()){
			in >> st1 >> en1;
			in2 >> st2 >> en2;
			if(anc.mates_pairs[read_id].mark_unique >= 0){
				if(anc.mates_pairs[read_id].mark_unique == 0 ){//unmapped to any locations: map always
					if(anc.reads_len1[read_id] > 64 || anc.reads_len2[read_id] > 64)
						anc.pair_reads_long(read_id, read, read2, contain_Ns, astat,  anc.reads_len1[read_id], anc.reads_len2[read_id], t, rs);
					else
						anc.pair_reads(read_id, read, read2, contain_Ns, astat,  anc.reads_len1[read_id], anc.reads_len2[read_id], t, rs);

				}
				else if(anc.mates_pairs[read_id].mark_unique >= 1 &&
					(anc.mates_pairs[read_id].mism1 > 0 || anc.mates_pairs[read_id].mism2 > 0)){//unmapped to any locations: map always
					if(anc.reads_len1[read_id] > 64 || anc.reads_len2[read_id] > 64)
						anc.pair_reads_long(read_id, read, read2, contain_Ns, astat,  anc.reads_len1[read_id], anc.reads_len2[read_id], t, rs);
					else
						anc.pair_reads(read_id, read, read2, contain_Ns, astat,  anc.reads_len1[read_id], anc.reads_len2[read_id], t, rs);

				}

			}//if valid
			read_id++;
if(read_id % 100000 == 0)
cout << read_id << endl;
			in >> read;
			in2 >> read2;
		
		}//while
	}//if bisulfited
	else{
		in >> read;
		in2 >> read2;

		read_id = 0;

		while(!in.eof()){
			in >> st1 >> en1;
			in2 >> st2 >> en2;
			if(anc.mates_pairs[read_id].mark_unique >= 0){
				if(anc.mates_pairs[read_id].mark_unique == 0 ){//unmapped to any locations: map always
					if(anc.reads_len1[read_id] > 64 || anc.reads_len2[read_id] > 64)
						anc.pair_reads_normal_long(read_id, read, read2, contain_Ns, astat,  anc.reads_len1[read_id], anc.reads_len2[read_id], t, rs);
					else
						anc.pair_reads_normal(read_id, read, read2, contain_Ns, astat,  anc.reads_len1[read_id], anc.reads_len2[read_id], t, rs);

				}
				else if(anc.mates_pairs[read_id].mark_unique >= 1 &&
					(anc.mates_pairs[read_id].mism1 > 0 || anc.mates_pairs[read_id].mism2 > 0)){//unmapped to any locations: map always
					if(anc.reads_len1[read_id] > 64 || anc.reads_len2[read_id] > 64)
						anc.pair_reads_normal_long(read_id, read, read2, contain_Ns, astat,  anc.reads_len1[read_id], anc.reads_len2[read_id], t, rs);
					else
						anc.pair_reads_normal(read_id, read, read2, contain_Ns, astat,  anc.reads_len1[read_id], anc.reads_len2[read_id], t, rs);


				}
		

			}//if valid

			read_id++;
if(read_id % 100000 == 0)
cout << read_id << endl;

			in >> read;
			in2 >> read2;
		
		}//while
	
	
	}//normal reads
	in.close(); in.clear();
	in2.close(); in2.clear();


	}//if pair-end
	else{//singles

		long int contain_Ns = 0;
		long int thresh_mism = 0;

		int r;

		int chr_id;
		in.open(anc.reads_file1.c_str(), ios::in);
		if(!in){
			cout << "can't open " << anc.reads_file1 << " reads_file " << endl;
			exit(0);
		}
long int check = 0;
		if(anc.is_bs == "1"){
			in >> read;
			read_id = 0;
			while(!in.eof()){
				in >> st1 >> en1;

				if(anc.mates_pairs[read_id].mark_unique >= 0){
					if(anc.mates_pairs[read_id].mark_unique == 0 )
					{//if have not been mapped, try to map
						if(anc.A_rich == false){
							if(anc.reads_len1[read_id] > 64)
								anc.pair_reads_singles_long(read_id, read, contain_Ns, astat, anc.reads_len1[read_id], t, rs);
							else	
								anc.pair_reads_singles(read_id, read, contain_Ns, astat, anc.reads_len1[read_id], t, rs);
						}
						else{
							if(anc.reads_len1[read_id] > 64)
								anc.pair_reads_singlesA_long(read_id, read, contain_Ns, astat, anc.reads_len1[read_id], t, rs);
							else	
								anc.pair_reads_singlesA(read_id, read, contain_Ns, astat, anc.reads_len1[read_id], t, rs);

					
						}//A-rich is true
					
					}//if not mapped yet
					else if(anc.mates_pairs[read_id].mark_unique >= 1 && anc.mates_pairs[read_id].mism1 > 0){
						if(anc.A_rich == false){
							if(anc.reads_len1[read_id] > 64)
								anc.pair_reads_singles_long(read_id, read, contain_Ns, astat, anc.reads_len1[read_id], t, rs);
							else	
								anc.pair_reads_singles(read_id, read, contain_Ns, astat, anc.reads_len1[read_id], t, rs);
						}
						else{
							if(anc.reads_len1[read_id] > 64)
								anc.pair_reads_singlesA_long(read_id, read, contain_Ns, astat, anc.reads_len1[read_id], t, rs);
							else	
								anc.pair_reads_singlesA(read_id, read, contain_Ns, astat, anc.reads_len1[read_id], t, rs);

					
						}//A-rich is true
					
					}
				}//if valid read
				read_id++;
if(read_id % 100000 == 0)
cout << read_id << endl;
				
				in >> read;
			
			}//while
			in.close(); in.clear();

		}//if bs
		else{
			in >> read;
			read_id = 0;
			while(!in.eof()){
				in >> st1 >> en1;
				if(anc.mates_pairs[read_id].mark_unique >= 0){
					if(anc.mates_pairs[read_id].mark_unique == 0 ){
							if(anc.reads_len1[read_id] > 64)
								anc.pair_reads_normal_singles_long(read_id, read, contain_Ns, astat, anc.reads_len1[read_id], t, rs);
							else	
								anc.pair_reads_normal_singles(read_id, read, contain_Ns, astat, anc.reads_len1[read_id], t, rs);
					}
					else if(anc.mates_pairs[read_id].mark_unique >= 1 && anc.mates_pairs[read_id].mism1 > 0){
							if(anc.reads_len1[read_id] > 64)
								anc.pair_reads_normal_singles_long(read_id, read, contain_Ns, astat, anc.reads_len1[read_id], t, rs);
							else	
								anc.pair_reads_normal_singles(read_id, read, contain_Ns, astat, anc.reads_len1[read_id], t, rs);
					
					}
				}//if valid read
				read_id++;	
if(read_id % 100000 == 0)
cout << read_id << endl;
				
				in >> read;
			
			}//while
			in.close(); in.clear();
		
		}//normal

	}//else singles

	}//for rs, all reads seeds for current gen-seeds

		anc.release_space();
	}//for t, all seeds

		anc.print(astat);
		anc.out_pairs.close(); anc.out_pairs.clear();
		astat.unmapped = (read_id) - (astat.mapped_pairs) - (astat.ambiguous) ;
		cout << "uniquely mapped, ambiguous singles, unmapped:" << endl;

		cout << astat.mapped_pairs << endl;
		cout << astat.ambiguous << endl;
		cout << astat.unmapped << endl;

		cout << "The number of reads not processed because of Ns or wrong length is " << astat.invalid << endl;

	

	return 0;
}//main


