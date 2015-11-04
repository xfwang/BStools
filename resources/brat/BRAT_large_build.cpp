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



void convert_lmer_c(string current_chrom, unsigned long int &gb,  int num_iter)
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

	string mes("\nUSAGE: brat-large-build -r <references> -P <directory name> [OPTIONS]\n\nOptions:\n");
	string gen_ref_opt("  -r <references-file>    file with file names containing references (one\n                          reference per file)");
	string is_bs(      "  -bs                     set BS-mapping (if not specified,\n                          normal mapping is done)");
	string first_bits( "  -f <integer >= 24>      the first <int> bases, within which non-BS-mismatches are NOT allowed");
	string pref(       "  -P <directory name>     the output files will be saved in this directory");
	string speed_mode( "  -S                      set Speed-Mode: space will be doubled, but running time will significantly decrease");

	cout << mes << endl;
	cout << gen_ref_opt << endl  
			<< is_bs << endl
			<< first_bits << endl 
			<< pref << endl
			<< speed_mode << endl << endl;

	string mes2(      "Limitations:\n");
	string ref_size(  "  a single reference size <= 2^32");
	string chrom_ids( "  number of references <= 2^32");
	string space(     "  space <= (the largest reference size)*4.5 + 16*(2^24) + 269 + 24Reads");
	
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



class Anchor{
public:
	Anchor(int argc, char* argv[]);
	~Anchor();

	vector<string> names;
	map<unsigned int, string> chroms_names;
	vector<unsigned long int> chroms_starts;
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

	
	vector< unsigned long int* > ta_genome;
	vector< unsigned long int* > cg_genome;
	vector< unsigned long int* > ns_genome; //presence of Ns

	unsigned long int hash_table_size;

	int parse_options(int argc, char* argv[]);

	void read_genome(int chr);

	void find_pos_sizes(int t);
	void make_pos_index(int t);

	void print(int u);
	void print_info();

	vector<Mates> mates_pairs;


	void convert_gen_mer_ta(unsigned long int &ta_gen_lmer,
					 unsigned long int gen_pos,
					 unsigned long int chrom_id,
					 unsigned long int aread_len);
	void convert_gen_mer_cg(unsigned long int &ta_gen_lmer,
					 unsigned long int gen_pos,
					 unsigned long int chrom_id,
					 unsigned long int aread_len);
	void convert_gen_mer_ns(unsigned long int &ta_gen_lmer,
					 unsigned long int gen_pos,
					 unsigned long int chrom_id,
					 unsigned long int aread_len);

	void convert_gen_mer_ta_long(vector<unsigned long int> &ta_gen_lmer,
					 unsigned long int gen_pos,
					 unsigned long int chrom_id,
					 unsigned long int aread_len, int asize);
	void convert_gen_mer_cg_long(vector<unsigned long int> &ta_gen_lmer,
					 unsigned long int gen_pos,
					 unsigned long int chrom_id,
					 unsigned long int aread_len, int asize);

	void convert_gen_mer_ns_long(vector<unsigned long int> &ta_gen_lmer,
					 unsigned long int pos,
					 unsigned long int chrom_id, 
					 unsigned long int aread_len, int asize);

	int count_mism(unsigned long int bm);
	int count_mism_long(vector<unsigned long int> &reads, vector<unsigned long int> &genmers,
		vector<unsigned long int> &ta_check, int asize);

	void pair_reads(unsigned long int read_id, string aread1, string aread2,
				long int &contain_Ns, Stat &astat, short int len1, short int len2, int t, int rs, int chr);

	void pair_reads_long(unsigned long int read_id, string aread1, string aread2,
				long int &contain_Ns, Stat &astat, short int len1, short int len2, int t, int rs, int chr);

	void pair_reads_normal(unsigned long int read_id, string aread1, string aread2,
				long int &contain_Ns, Stat &astat, short int len1, short int len2, int t, int rs, int chr);
	void pair_reads_normal_long(unsigned long int read_id, string aread1, string aread2,
				long int &contain_Ns, Stat &astat, short int len1, short int len2, int t, int rs, int chr);

	void pair_reads_singles(unsigned long int read_id, string aread1,
				long int &contain_Ns, Stat &astat, short int len1, int t, int rs, int chr);
	void pair_reads_singles_long(unsigned long int read_id, string aread1,
				long int &contain_Ns, Stat &astat, short int read_len1,  int t, int rs, int chr );

	void pair_reads_singlesA(unsigned long int read_id, string aread1,
				long int &contain_Ns, Stat &astat, short int len1, int t, int rs, int chr);
	void pair_reads_singlesA_long(unsigned long int read_id, string aread1,
				long int &contain_Ns, Stat &astat, short int len1, int t, int rs, int chr);

	void pair_reads_normal_singles(unsigned long int read_id, string aread1,
				long int &contain_Ns, Stat &astat, short int len1, int t, int rs, int chr);

	void pair_reads_normal_singles_long(unsigned long int read_id, string aread1,
				long int &contain_Ns, Stat &astat, short int len1, int t, int rs, int chr);

	unsigned long int hash2m(unsigned long int ta);
	unsigned long int hash2c(unsigned long int ta, unsigned long int cg);
	unsigned long int hash2gs(unsigned long int ta, unsigned long int cg);
	unsigned long int hash2cs(unsigned long int ta, unsigned long int cg);
	unsigned long int hash_aread(unsigned long int ta, unsigned long int cg, 
		bool is_tac);

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

	string prefix;
	string speed_mode;

	int mism_thresh;

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

	void release_space(int chr);
	unsigned long int dif_len_mer;

	bool out_amb;
	void free_pos_index();

	int first_bits;
	long int cur_chrom;

	bool fast;
};	

void Anchor::release_space(int chr)
{
	long int i;
		delete [] ta_genome[chr];
		delete [] cg_genome[chr];
		delete [] ns_genome[chr];

	for(i = 0; i < hash_table_size; i++)
		delete [] pos_index[i];

}//release pos index
void Anchor::free_pos_index(){
	long int i;
	if(unmapped == true){
		out_unmapped1.close();
		out_unmapped2.close();
	}
	for(i = 0; i < hash_table_size; i++)
		delete [] pos_index[i];

}

Anchor::~Anchor(){

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
	fast = false;
	speed_mode = "0";
	out_amb = false;
	unmapped = false;
	A_rich = false;
	mism_thresh = 0;
	min_insert = 100;
	max_insert = 300;
	is_bs = "0";
	is_pe = "0";
	first_bits = 24;
	
	
	int res = parse_options(argc, argv);

        if(first_bits < 24 || first_bits > 64){                
		cout << "\nERROR: option f must NOT be less than 24 or greater than 64" << endl;
                        exit(0);
        }

	if(speed_mode == "1")
		fast = true;
	if(res == -1){
		usage();
		exit(0);
	}
	if(is_bs != "1" && A_rich == true){
		cout << "\nOption A works only with BS-mapping. Ignoring option A" << endl;
		A_rich = false;
	
	}
	string unm1 = "unm_" + reads_file1;
	string unm2 = "unm_" + reads_file2;
	if(unmapped == true){
		out_unmapped1.open(unm1.c_str(), ios::out);
		if(is_pe == "1"){
			out_unmapped2.open(unm2.c_str(), ios::out);
		}
	}//if unm
	int min_len = 10000, cur_len;

	ifstream in;
	cur_chrom = 0;
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
	len_mer = max(24, first_bits);

	len_mer = min(len_mer, SIZE);


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
	ta_genome.resize(num_chroms);
	cg_genome.resize(num_chroms);
	ns_genome.resize(num_chroms);

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

		size_chrom.push_back(cur_size);
		long int new_size = (cur_size >> (two_in_Size)) + 1;//(SIZE*3);// =/(8*64) 1pb takes 1 bit, and 1 long takes 64 bits
		
		chroms_starts.push_back(new_size);

		in.close(); in.clear();
	}


	pos_sizes.resize(hash_table_size);
	pos_index.resize(hash_table_size);


	MASK_READ = MASK_LMER;//for initialization only

}//Anchor()

void Anchor::make_seeds(){
	


	if(is_pe == "1"){
			gen_seeds.push_back(MASK_LMER);
			gen_seeds2.push_back(MASK_LMER);

			vector<unsigned long int> dum;
			dum.push_back(MASK_LMER);
			seeds.push_back(dum);
			seeds2.push_back(dum);

	}//pairs
	else{
			gen_seeds.push_back(MASK_LMER);
			vector<unsigned long int> dum;
			dum.push_back(MASK_LMER);

			seeds.push_back(dum);

	
	}//singles

}//make seeds



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

	if(fast == true){
		if(is_tac == true)
			return hash2cs(ta, cg);
		else
			return hash2gs(ta, cg);
	}
	else
		return hash2m(ta);
}//hash_aread


void Anchor::find_pos_sizes(int u)
{
	int t = 0;

	long int i;
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

		long int cur_size = 0;
		long int chrom_id = u;
		lmer_pos = 0;
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


			stop = real_size - lmer_1;
			for(; last_lmer_char < real_size; last_lmer_char++){

				char ch_seed = buffer[last_lmer_char];
				if(ch_seed != '\n'){

					next_lmer(ch_seed, ta_seed, ns_seed, cg_seed, MASK_LMER);

					if(ns_seed == 0){


						ta_lmer = ta_seed & gen_seeds[t];
						cg_lmer = cg_seed & gen_seeds[t];
						if(is_bs == "1"){
							if(fast == true){
								ta_ind = hash_aread(ta_lmer, cg_lmer, true);
								ta_ind2 = hash_aread(ta_lmer, cg_lmer, false);
								if(ta_ind2 != ta_ind)
									pos_sizes[ta_ind2]++;
							}
							else
								ta_ind = hash2m(ta_lmer);
if(ta_ind 	>= hash_table_size)
{
cout << "ta_ind is out of boundary " << endl;
exit(0);
}		
							pos_sizes[ta_ind]++;
						}
						else{
							ta_ind = hash2c(ta_lmer, cg_lmer);
							pos_sizes[ta_ind]++;

						
						}//else
					}//if No Ns
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

void Anchor::make_pos_index(int u)
{
	int t = 0;

	unsigned long int i, ch, j,  y;


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



		chrom_id = u;
		lmer_pos = 0;
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
		
		unsigned long int count_bits = 0;

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

						if(is_bs == "1"){
							if(fast == true){
								ta_ind = hash_aread(ta_seed, cg_seed, true);
								ta_ind2 = hash_aread(ta_seed, cg_seed, false);

								long int cur_pointer = pos_sizes[ta_ind];

								pos_index[ta_ind][cur_pointer] = lmer_pos;
								pos_sizes[ta_ind] = cur_pointer + 1;

								if(ta_ind2 != ta_ind){
									cur_pointer = pos_sizes[ta_ind2];
									pos_index[ta_ind2][cur_pointer] = lmer_pos;
									pos_sizes[ta_ind2] = cur_pointer + 1;
								}
							}
							else{
								ta_ind  = hash2m(ta_seed);

								long int cur_pointer = pos_sizes[ta_ind];

								pos_index[ta_ind][cur_pointer] = lmer_pos;
								pos_sizes[ta_ind] = cur_pointer + 1;

							}

if(ta_ind 	>= hash_table_size)
{
cout << "make index: ta_ind is out of boundary " << endl;
exit(0);
}		

						}//if
						else{ //normal mapping 
							ta_ind = hash2c(ta_seed, cg_seed);
							long int cur_pointer = pos_sizes[ta_ind];
							pos_index[ta_ind][cur_pointer] = lmer_pos;							pos_sizes[ta_ind] = cur_pointer + 1;
						}//else
					}//if no Ns
					lmer_pos++;

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
		int t1 = t+1;
		cout << "round " << t1 << ": " << names[u] << " is hashed " << endl;


	delete [] buffer;

}//make_pos_index

//read genome
void Anchor::read_genome(int u)
{
	long int new_size = chroms_starts[u];
		ta_genome[u] = new unsigned long int[new_size];
		cg_genome[u] = new unsigned long int[new_size];
		ns_genome[u] = new unsigned long int[new_size];

		assert(ta_genome[u] != 0);
		assert(cg_genome[u] != 0);
		assert(ns_genome[u] != 0);


	//out2 is a file handler for Pos Index
	unsigned long int i, ch, j, y;
	string aname;
	
	string next, line, chrom;
	unsigned int chrom_id = 0;

	//the first chromosome has starts 0/0
	unsigned long int cur_chrom_cg_size = 0;
	unsigned long int cur_chrom_ta_size = 0;

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



		chrom_id = u;
		index = 0;
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

						ta_genome[chrom_id][index] =ta ;
						cg_genome[chrom_id][index] = cg ;
						ns_genome[chrom_id][index] = ns;
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
	if(chroms_starts[chrom_id] <= index){

		cout << "index is out of boundary at lmer pos " << lmer_pos << endl;
		exit(0);
	}
						ta_genome[chrom_id][index] = ta;
						cg_genome[chrom_id][index] = cg;
						ns_genome[chrom_id][index] = ns;

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

		cg = (cg << (SIZE - remainder_cg));
		ta = (ta << (SIZE - remainder_cg)) ;
		ns = (ns << (SIZE - remainder_cg));

if(chroms_starts[chrom_id] <= index){

	cout << "index is out of boundary at lmer pos " << lmer_pos << endl;
	exit(0);
}

		ta_genome[chrom_id][index] =ta;
		cg_genome[chrom_id][index] = cg ;
		ns_genome[chrom_id][index] = ns;

		index++;

		cout << names[u] << " is preprocessed " << endl;


	delete [] buffer;

}//
//read genome
////////////////////////// LONG READS

int Anchor::parse_options(int argc, char* argv[])
{
	int res = 0;

	long int i, j;
	for(i = 1; i < argc; i++){
		if(argv[i][0] != '-')
			return -1;
		switch(argv[i][1]){
		case 'r': genome_names = argv[++i]; break;
		case 'b': is_bs = "1"; break;
		case 'f': first_bits = atoi(argv[++i]); break;
		case 'P': prefix = argv[++i]; break;
		case 'S': speed_mode = "1"; break;
		case 'h': usage(); break;
		default: cerr << "ERROR: Invalid char=" << argv[i][1] << "=this" << endl; usage(); exit(0); break;
		}
	}//for i
	return res;
}//parse_opt()
void Anchor::print_info(){
	long int i, j, y, u;

	string output4 = prefix + "/" + "INFO.txt";
	ofstream out2;
	
	out2.open(output4.c_str(), ios::out);
	out2 << num_chroms << endl;
	out2 << hash_table_size << endl;
	out2 << first_bits << endl;
	out2 << is_bs << endl;
	out2 << speed_mode << endl;

	map<unsigned int, string>::iterator beg = chroms_names.begin();
	map<unsigned int, string>::iterator fin = chroms_names.end();

	for(; beg != fin; beg++){
		out2 << beg->first << "\t" << beg->second << endl;
	}
	long int asize = num_chroms;
	for(i = 0; i < asize; i++)
		out2 << names[i] << endl;

	for(i = 0; i < asize; i++){
		out2 << size_chrom[i] << endl;
	}

	for(i = 0; i < asize; i++)
		out2 << chroms_starts[i] << endl;

	out2.close(); out2.clear();

}

void Anchor::print(int u){
	long int i, j, y;
	long int double_byte = byte_size << 1;
	long int triple_byte = double_byte + byte_size;
	long int quadra = double_byte + double_byte;
	long int five = quadra + byte_size;
	long int six = five + byte_size;
	long int seven = six + byte_size;

	string cur_chrom;
	map<unsigned int, string>::iterator fin = chroms_names.end();
	map<unsigned int, string>::iterator afind = chroms_names.find(u);
	cur_chrom = afind->second;

	ofstream out1, out2, out3;

	string output1 = prefix + "/" + cur_chrom + ".ht";//hash table
	string output2 = prefix + "/" + cur_chrom + ".hs";//hash sizes
	string output3 = prefix + "/" + cur_chrom + ".bg";//binary genome
	out1.open(output1.c_str(), ios::out);

if(!out1){
cout << "Can't open " << output1 << " for output." << endl;

}
else
cout << "Opening " << output1 << " to output a hash table." << endl;

	long int asize = 4;
	char abuffer[4];
	for(i = 0; i < hash_table_size; i++){
		
		for(j = 0; j < pos_sizes[i]; j++){
			unsigned int cur = pos_index[i][j];

			unsigned int x = (cur >> triple_byte) & MASK_BYTE;
			abuffer[0] = char(x);
			x = (cur >> double_byte) & MASK_BYTE;
			abuffer[1] = char(x);
			x = (cur >> byte_size) & MASK_BYTE;
			abuffer[2] = char(x);
			x = cur & MASK_BYTE;
			abuffer[3] = char(x);

			out1.write(&abuffer[0], asize);

		}
	}
	out1.close(); out1.clear();

	out1.open(output3.c_str(), ios::out);
if(!out1){
cout << "Can't open " << output3 << " for output." << endl;

}
else
cout << "Opening " << output3 << " to output binary genome." << endl;

	char abuffer2[8];
	for(i = 0; i < chroms_starts[u]; i++){
		unsigned long int cur = ta_genome[u][i];
		unsigned int x = (cur >> seven) & MASK_BYTE;
		abuffer2[0] = char(x);
		x = (cur >> six) & MASK_BYTE;
		abuffer2[1] = char(x);
		x = (cur >> five) & MASK_BYTE;
		abuffer2[2] = char(x);
		x = (cur >> quadra) & MASK_BYTE;
		abuffer2[3] = char(x);
		x = (cur >> triple_byte) & MASK_BYTE;
		abuffer2[4] = char(x);
		x = (cur >> double_byte) & MASK_BYTE;
		abuffer2[5] = char(x);
		x = (cur >> byte_size) & MASK_BYTE;
		abuffer2[6] = char(x);
		x = cur & MASK_BYTE;
		abuffer2[7] = char(x);

		out1.write(&abuffer2[0], 8);
	}//for i
	for(i = 0; i < chroms_starts[u]; i++){
		unsigned long int cur = cg_genome[u][i];

		unsigned int x = (cur >> seven) & MASK_BYTE;
		abuffer2[0] = char(x);
		x = (cur >> six) & MASK_BYTE;
		abuffer2[1] = char(x);
		x = (cur >> five) & MASK_BYTE;
		abuffer2[2] = char(x);
		x = (cur >> quadra) & MASK_BYTE;
		abuffer2[3] = char(x);
		x = (cur >> triple_byte) & MASK_BYTE;
		abuffer2[4] = char(x);
		x = (cur >> double_byte) & MASK_BYTE;
		abuffer2[5] = char(x);
		x = (cur >> byte_size) & MASK_BYTE;
		abuffer2[6] = char(x);
		x = cur & MASK_BYTE;
		abuffer2[7] = char(x);

		out1.write(&abuffer2[0], 8);
	}//for i
	for(i = 0; i < chroms_starts[u]; i++){
		unsigned long int cur = ns_genome[u][i];

		unsigned int x = (cur >> seven) & MASK_BYTE;
		abuffer2[0] = char(x);
		x = (cur >> six) & MASK_BYTE;
		abuffer2[1] = char(x);
		x = (cur >> five) & MASK_BYTE;
		abuffer2[2] = char(x);
		x = (cur >> quadra) & MASK_BYTE;
		abuffer2[3] = char(x);
		x = (cur >> triple_byte) & MASK_BYTE;
		abuffer2[4] = char(x);
		x = (cur >> double_byte) & MASK_BYTE;
		abuffer2[5] = char(x);
		x = (cur >> byte_size) & MASK_BYTE;
		abuffer2[6] = char(x);
		x = cur & MASK_BYTE;
		abuffer2[7] = char(x);

		out1.write(&abuffer2[0], 8);
	}//for i	

	out1.close(); out1.clear();


	out2.open(output2.c_str(), ios::out);
if(!out2){
cout << "Can't open " << output2 << " for output." << endl;

}
else
cout << "Opening " << output2 << " to output sizes of hash table entries." << endl;


	for(i = 0; i < hash_table_size; i++){
		unsigned long int cur = pos_sizes[i];

		unsigned int x = (cur >> seven) & MASK_BYTE;
		abuffer2[0] = char(x);
		x = (cur >> six) & MASK_BYTE;
		abuffer2[1] = char(x);
		x = (cur >> five) & MASK_BYTE;
		abuffer2[2] = char(x);
		x = (cur >> quadra) & MASK_BYTE;
		abuffer2[3] = char(x);
		x = (cur >> triple_byte) & MASK_BYTE;
		abuffer2[4] = char(x);
		x = (cur >> double_byte) & MASK_BYTE;
		abuffer2[5] = char(x);
		x = (cur >> byte_size) & MASK_BYTE;
		abuffer2[6] = char(x);
		x = cur & MASK_BYTE;
		abuffer2[7] = char(x);

		out2.write(&abuffer2[0], 8);
	}
	out2.close(); out2.clear();

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

	if(argc < 3){
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

	int seeds_size = anc.gen_seeds.size();

	long int num_chroms = anc.num_chroms;

	for(u = 0; u < num_chroms; u++){

		anc.cur_chrom = u;
		anc.read_genome(u);
		anc.find_pos_sizes(u);
		anc.make_pos_index(u);
		anc.print(u);

		anc.release_space(u);
	}//for u, all chroms
	anc.print_info();


	return 0;
}//main


