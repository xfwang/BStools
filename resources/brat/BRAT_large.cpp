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

	string mes("\nUSAGE: brat-large -r <references> -s <input single reads> -o <output> [OPTIONS]\n\nOptions:\n");
	string gen_ref_opt("  -r <references-file>    file with file names containing references (one\n reference per file)");
	string singles(    "  -s <input-reads-file>   file with single reads (or queries)");
	string output(     "  -o <output-file>");
	string is_pair(    "  -pe                     set Pair-End mapping (if not specified,\n mapping for single reads is done)");
	string is_bs(      "  -bs                     set BS-mapping (if not specified,\n normal mapping is done)");
	string first(      "  -1 <paired-ends-file1>  5`mates");
	string second(     "  -2 <paired-ends-file2>  3`mates");
	string min_insert( "  -i <min>                min insert");
	string max_insert( "  -a <max>                max insert");
	string arich(      "  -A                      singles from paired-ends-file2,\n that are mapped to an A-rich strand");
	string unmapped(   "  -u                      output unmapped reads/pairs");
	string out_amb(    "  -M                      output ambiguous reads/pairs");
	string mism(	   "  -m <non-negative integer>  the number of non-BS-mismatches allowed");
	string first_bits( "  -f <integer >= 24>         the first f bases, within which non-BS-mismatches are NOT allowed");
	string pref(       "  -P <directory name for index files> with this option need to use brat-large-build first");
	string speed_mode( "  -S                      set Speed-Mode that uses double space, but speeds running time more than double");

	cout << mes << endl;
	cout << gen_ref_opt << endl << singles << endl
        	<< output << endl
		<< is_pair << endl << is_bs << endl
		<< first << endl << second << endl << min_insert << endl << max_insert
		<< endl << arich << endl << unmapped << endl << out_amb << endl 
		<< pref << endl << speed_mode << endl;

	string mes2(      "Limitations:\n");
	string ref_size(  "  a single reference size <= 2^32");
	string chrom_ids( "  number of references <= 2^32");
	string space(     "  space <= (largest reference size)*4 + 16*(2^24) + 269 + 24*Reads");
	
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
	void find_pos_sizes(int chr);
	void make_pos_index(int chrt);

	void read_genome2(int chr);
	void find_pos_sizes2(int chr);
	void make_pos_index2(int chrt);

	void print(Stat &astat);

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


	string genome_names;

	string prefix;

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

	unsigned long int MASK_chrom;
	ofstream out_pairs;

	string is_pe;
	string speed_mode;

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

	unsigned long int hash2m(unsigned long int ta);
	unsigned long int hash2c(unsigned long int ta, unsigned long int cg);
	unsigned long int hash2gs(unsigned long int ta, unsigned long int cg);
	unsigned long int hash2cs(unsigned long int ta, unsigned long int cg);
	unsigned long int hash_aread(unsigned long int ta, unsigned long int cg, 
		bool is_tac);

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

	long int i;

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
	speed_mode = "0";
	fast = false;
	prefix = "";
	out_amb = false;
	unmapped = false;
	A_rich = false;
	mism_thresh = 0;
	min_insert = 100;
	max_insert = 300;
	is_bs = "0";
	is_pe = "0";
	first_bits = 0;

	int res = parse_options(argc, argv);

	if(res == -1){
		usage();
		exit(0);
	}
	if(is_bs != "1" && A_rich == true){
		cout << "\nOption A works only with BS-mapping. Ignoring option A" << endl;
		A_rich = false;
	
	}
	string unm1 = reads_file1 + ".unm";
	string unm2 = reads_file2 + ".unm";
	if(unmapped == true){
		out_unmapped1.open(unm1.c_str(), ios::out);
		if(is_pe == "1"){
			out_unmapped2.open(unm2.c_str(), ios::out);
		}
	}//if unm
	int min_len = 10000, cur_len;
	
	long int st1, en1, st2, en2;

	ifstream in;
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
		in >> st1 >> en1;
		total_reads++;
		in >> aread;
	
	}
	in.close(); in.clear();

	mates_pairs.resize(total_reads);
	reads_len1.resize(total_reads);
/////////////////////////////////////////////
	if(speed_mode == "1")
		fast = true;

if(prefix.length() > 0){

	//Since index is built to map reads of various length, index is built for the min-length, f=24,
	//unless user provided a different value for option f


	int temp_first_bits;
	string input = prefix + "/" + "INFO.txt";
	in.open(input.c_str(), ios::in);

	in >> num_chroms ;
	in >> hash_table_size;
	in >> temp_first_bits;
	string temp_bs;
	in >> temp_bs;
	in >> speed_mode;

	if(speed_mode == "1")
		fast = true;
	else
		fast = false;

	if(temp_bs != is_bs){
		cout << "\nERROR: While building index the option bs was " << temp_bs << endl;
		cout << "Now bs option was = " << is_bs << endl;
		cout << "The option bs must be same for building index and for mapping." << endl;
		exit(0);
	}
	if(first_bits > 0){//user provided f option
		if(temp_first_bits != first_bits){
			cout << "\nWARNING: Option f is different from option f provided with brat-large-build." << endl;
			cout << "Continue with option f provided with brat-large-build" << endl;
		}
	}

	first_bits = temp_first_bits;

	unsigned int chr_id;
	string aname;
	long int i;
	for(i = 0; i < num_chroms; i++){
		in >> chr_id;
		in >> aname;
		chroms_names.insert(make_pair(chr_id, aname));
	}

	for(i = 0; i < num_chroms; i++){
		in >> aname;
		names.push_back(aname);
	}
	long int asize;
	for(i = 0; i < num_chroms; i++){
		in >> asize;
		size_chrom.push_back(asize);
	}

	for(i = 0; i < num_chroms; i++){
		in >> asize;
		chroms_starts.push_back(asize);
	}

	in.close(); in.clear();

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
	len_mer =  first_bits;// rage first_bits is [24...64] is checked when index is built 


	unsigned long int ta_ind = 0;

	float base = 2;
	float apow = 24;


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

	string line;

	//find out sizes of chromosomes (or references)
	ta_genome.resize(num_chroms);
	cg_genome.resize(num_chroms);
	ns_genome.resize(num_chroms);

	long int cur_size;

	pos_sizes.resize(hash_table_size);
	pos_index.resize(hash_table_size);


	in.open(reads_file1.c_str(), ios::in);
	if(!in){
		cout << "\nERROR: could not open " << reads_file1 << endl;
		exit(0);
	}


	in >> aread;
	long int read_id = 0;
	while(!in.eof()){
		in >> st1 >> en1;
		cur_len = aread.length();
		reads_len1[read_id] = cur_len;
		if(cur_len < first_bits )
			mates_pairs[read_id].mark_unique = -1;
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
		in >> st2 >> en2;
		cur_len = aread.length();
		reads_len2[read_id] = cur_len;
		if(cur_len < first_bits )
			mates_pairs[read_id].mark_unique = -1;
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

	cout << "\nPlease note that with option P, reads of length less than f will not be mapped (where f is option set during brat-large-build)" << endl;


}//if prefix is provided
else{

	if(first_bits == 0)
		first_bits = 24;

	in.open(reads_file1.c_str(), ios::in);
	if(!in){
		cout << "\nERROR: could not open " << reads_file1 << endl;
		exit(0);
	}
	long int i;

	in >> aread;
	long int read_id = 0;
	while(!in.eof()){
		in >> st1 >> en1;
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
		in >> st1 >> en1;
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

		size_chrom.push_back(cur_size);
		long int new_size = (cur_size >> (two_in_Size)) + 1;//(SIZE*3);// =/(8*64) 1pb takes 1 bit, and 1 long takes 64 bits
		
		chroms_starts.push_back(new_size);

		in.close(); in.clear();
	}


	pos_sizes.resize(hash_table_size);
	pos_index.resize(hash_table_size);
}
/////////////////////////////////////////////

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

void Anchor::find_pos_sizes2(int u)
{
	ifstream in;
	map<unsigned int, string>::iterator fin2 = chroms_names.end();
	map<unsigned int, string>::iterator afind2 = chroms_names.find(u);
	string cur_chrom = afind2->second;

	string input = prefix + "/" + cur_chrom + ".hs";

	in.open(input.c_str(), ios::in);
	if(!in){
		cout << "ERROR: cannot open " << input << endl;
		exit(0);
	}

	long int cur_size = hash_table_size << 3;//*8
	char *abuf = new char[cur_size];
	in.read(abuf, cur_size);
	long int read_char = in.gcount();
	if(read_char < cur_size){
		cout << "\nERROR: wrong input in the file " << input << endl;
		cout << "\nThe number of characters in this file should be not less than " << cur_size << endl;
		exit(0);
	}
	long int j = 0;
	long int y, i;
	for(i = 0; i < hash_table_size; i++){
		unsigned long int res = 0;

		for(y = 0; y < 8; y++){
			res <<= byte_size;

			char x = abuf[j];
			j++;
			res |= (int(x) & MASK_BYTE);

		}
		pos_sizes[i] = res;
	}//for i
	in.close(); in.clear();
	delete [] abuf;

}//find pos sizes2

void Anchor::make_pos_index2(int chr)
{
	unsigned long int i, ch, j, u, y;

	long int cur_size;

	long int buffer_size = 16777216;//divisible by 4
	char *buffer = new char[buffer_size];
	assert( buffer != 0);

	for(i = 0; i < hash_table_size; i++){
		cur_size = pos_sizes[i];
		pos_index[i] = new unsigned int[cur_size];
		assert( pos_index[i] != 0);

	}//for i
	int cur_chr = chr + 1;
	cout << "Space is allocated for ref " << cur_chr << endl;

	ifstream in;
	map<unsigned int, string>::iterator fin2 = chroms_names.end();
	map<unsigned int, string>::iterator afind2 = chroms_names.find(chr);
	string cur_chrom = afind2->second;

	string input = prefix + "/" + cur_chrom + ".ht";


	in.open(input.c_str(), ios::in);
	in.read(buffer, buffer_size);
	long int read_char = in.gcount();
	buffer_size = read_char;
/*
	if(read_char < buffer_size){
		cout << "\nERROR: Unexpected HASH TABLE size" << endl;
		exit(0);
	}
*/
  u = 0; i = 0;
	long int cur_count = 0;
	while(true){


		long int cur_ind = 0;
			
		cur_size = pos_sizes[i];
		cur_ind = 0;
		j = 0; 
		while(j < cur_size){
			unsigned int res = 0;
			for(y = 0; y < 4; y++){
				char x = buffer[u];
				u++;
				res <<= byte_size;
				res |= (int(x) & MASK_BYTE);
			}
			pos_index[i][cur_ind] = res;
			cur_ind++;
			if(u >= buffer_size){
				in.read(buffer, buffer_size);
				u = 0;
				if(in.gcount() == 0 ){
					break;
				}
			}


			j++;
		}//while
		if(in.gcount() == 0){
			break;
		}
		i++;
		if(i >= hash_table_size)
			break;
	}//while
	in.close(); in.clear();

	delete [] buffer;
	chr++;
	cout << "HASH TABLE is read for ref " << chr << endl;
}//make_pos_index2

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
						}//if
						else{ //normal mapping 
							ta_ind = hash2c(ta_seed, cg_seed);
							long int cur_pointer = pos_sizes[ta_ind];
							pos_index[ta_ind][cur_pointer] = lmer_pos;	
							pos_sizes[ta_ind] = cur_pointer + 1;
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

void Anchor::read_genome2(int u)
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
	

	ifstream in;
	long int last_lmer_char;
	long int binary_total = new_size;

	const long int buffer_size = binary_total << 3;//*8
	char *buffer = new char[buffer_size];
	assert(buffer != 0);

	map<unsigned int, string>::iterator fin2 = chroms_names.end();
	map<unsigned int, string>::iterator afind2 = chroms_names.find(u);
	string cur_chrom = afind2->second;

	string input = prefix + "/" + cur_chrom + ".bg";

	in.open(input.c_str(), ios::in);
	if(!in){
		cout << "\nERROR: could not open " << input << endl;
		exit(0);
	}
	in.read(buffer, buffer_size);
	long int read_char = in.gcount();
	if(read_char < buffer_size){
		cout << "\nERROR: Unexpected size of a binary genome" << endl;
		exit(0);
	}

	long int cur_size = binary_total;//*8
	j = 0;
	for(i = 0; i < cur_size; i++){
		
		unsigned long int res = 0;
		for(y = 0; y < 8; y++){
			res <<= byte_size;

			char x = buffer[j];
			j++;
			res |= (int(x) & MASK_BYTE);
		}
		ta_genome[u][i] = res;
	} //for i

	in.read(buffer, buffer_size);
	read_char = in.gcount();
	if(read_char < buffer_size){
		cout << "\nERROR: Unexpected size of a binary genome" << endl;
		exit(0);
	}

	j = 0;
	for(i = 0; i < cur_size; i++){
		
		unsigned long int res = 0;
		for(y = 0; y < 8; y++){
			res <<= byte_size;

			char x = buffer[j];
			j++;
			res |= (int(x) & MASK_BYTE);
		}
		cg_genome[u][i] = res;

	} //for i

	in.read(buffer, buffer_size);
	read_char = in.gcount();
	if(read_char < buffer_size){
		cout << "\nERROR: Unexpected size of a binary genome" << endl;
		exit(0);
	}

	j = 0;
	for(i = 0; i < cur_size; i++){
		
		unsigned long int res = 0;
		for(y = 0; y < 8; y++){
			res <<= byte_size;

			char x = buffer[j];
			j++;
			res |= (int(x) & MASK_BYTE);
		}
		ns_genome[u][i] = res;

	} //for i
	in.close(); in.clear();
	delete [] buffer;
	u++;
	cout << "\nRef " << u << " is read" << endl;
}//
//read genome
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

////////////////////////// LONG READS

void Anchor::convert_gen_mer_ta_long(vector<unsigned long int> &ta_gen_lmer,
					 unsigned long int pos,
					 unsigned long int chrom_id, 
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
		for(i = 0; i < asize  && j < chroms_starts[chrom_id]; i++){
			ta_gen_lmer.push_back(ta_genome[chrom_id][j]);

			j++;
		}//for i

		if(len_rest < SIZE)
			ta_gen_lmer[asize - 1] >>= (SIZE - len_rest);

	}
	else{

		//This is not done
		for(i = 0; i < asize && cur_byte < chroms_starts[chrom_id]; i++){
			first = ta_genome[chrom_id][cur_byte];
			second = ta_genome[chrom_id][cur_byte + 1];
			
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
					 unsigned long int chrom_id, 
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
		for(i = 0; i < asize - 1 && j < chroms_starts[chrom_id]; i++){
			ta_gen_lmer.push_back(cg_genome[chrom_id][j]);
			j++;
		}//for i
		if(len_rest < SIZE)
			ta_gen_lmer.push_back(cg_genome[chrom_id][j] >> (SIZE - len_rest));
	}
	else{

		for(i = 0; i < asize && cur_byte < chroms_starts[chrom_id]; i++){
			first = cg_genome[chrom_id][cur_byte];
			second = cg_genome[chrom_id][cur_byte + 1];
			
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
					 unsigned long int chrom_id, 
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
		for(i = 0; i < asize  && j < chroms_starts[chrom_id]; i++){
			ta_gen_lmer.push_back(ns_genome[chrom_id][j]);

			j++;
		}//for i
		if(len_rest < SIZE)
			ta_gen_lmer[asize - 1] >>= (SIZE - len_rest);

	}
	else{


		for(i = 0; i < asize && cur_byte < chroms_starts[chrom_id]; i++){
			first = ns_genome[chrom_id][cur_byte];
			second = ns_genome[chrom_id][cur_byte + 1];
			
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
					 unsigned long int chrom_id, 
					 unsigned long int aread_len)
{

	unsigned long int cur_byte, first,  second;
	//current binary char inside current block: there 8 characters fit in one binary characters
	//(pos >> two_in_byte) = pos / byte_size = number of binary characters needed to cover this pos
	cur_byte = (pos >> two_in_Size) ;

	unsigned long int char_pos = pos & MASK_SIZE;//pos of genome char inside byte-block

	//this is equivalent to bytes_in_lmer = ceil(aread_len/byte_size)
	first = ta_genome[chrom_id][cur_byte];

	if(char_pos == 0)
		ta_gen_lmer = first >> (SIZE - aread_len);
	else{
		cur_byte++;
		second = ta_genome[chrom_id][cur_byte];
		
		first <<= char_pos;
		second >>= (SIZE - char_pos);
		ta_gen_lmer = (first | second) >> (SIZE - aread_len);
	}


}//convert_gen_mer_ta


void Anchor::convert_gen_mer_cg(unsigned long int &ta_gen_lmer,
					 unsigned long int pos,
					 unsigned long int chrom_id, 
					 unsigned long int aread_len)
{
	unsigned long int cur_byte, first,  second;
	//current binary char inside current block: there 8 characters fit in one binary characters
	//(pos >> two_in_byte) = pos / byte_size = number of binary characters needed to cover this pos
	cur_byte = (pos >> two_in_Size) ;

	unsigned long int char_pos = pos & MASK_SIZE;//pos of genome char inside byte-block

	//this is equivalent to bytes_in_lmer = ceil(aread_len/byte_size)
	first = cg_genome[chrom_id][cur_byte];
	if(char_pos == 0)
		ta_gen_lmer = first >> (SIZE - aread_len);
	else{
		cur_byte++;
		second = cg_genome[chrom_id][cur_byte];
		
		first <<= char_pos;
		second >>= (SIZE - char_pos);
		ta_gen_lmer = (first | second) >> (SIZE - aread_len);
	}


}//convert_gen_mer_cg

void Anchor::convert_gen_mer_ns(unsigned long int &ta_gen_lmer,
					 unsigned long int pos,
					 unsigned long int chrom_id, 
					 unsigned long int aread_len)
{

	unsigned long int cur_byte, first,  second;
	//current binary char inside current block: there 8 characters fit in one binary characters
	//(pos >> two_in_byte) = pos / byte_size = number of binary characters needed to cover this pos
	cur_byte = (pos >> two_in_Size) ;

	unsigned long int char_pos = pos & MASK_SIZE;//pos of genome char inside byte-block

	//this is equivalent to bytes_in_lmer = ceil(aread_len/byte_size)
	first = ns_genome[chrom_id][cur_byte];

	if(char_pos == 0)
		ta_gen_lmer = first >> (SIZE - aread_len);
	else{
		cur_byte++;
		second = ns_genome[chrom_id][cur_byte];
		
		first <<= char_pos;
		second >>= (SIZE - char_pos);
		ta_gen_lmer = (first | second) >> (SIZE - aread_len);
	}


}//convert_gen_mer_ta


void Anchor::pair_reads(unsigned long int read_id, string aread1, string aread2,
				long int &contain_Ns, Stat &astat, short int read_len1, short int read_len2, int t, int rs, int chr )
{


	unsigned long int one = 1;

	long int chrom_id1 = chr;
	long int chrom_id2 = chr;
	long int chr_id2 = chr;
	long int asize,   pos1, pos2, dist, i;
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
	convert_lmer_c(aread1, c_read1, read_len1);
	convert_lmer_c(aread2, c_read2, read_len2);

  //rev-complement of c-reads

	int dif1 = read_len1 - len_mer;
	int dif2 = read_len2 - len_mer;

	ta_lmer1 = (ta_read1 >> dif1) & seeds[t][rs];
	ta_lmer2 = (ta_read2 >> dif2) & seeds2[t][rs];
	cg_lmer1 = (c_read1 >> dif1) & seeds[t][rs];

	cg_lmer2 = (c_read2 >> dif2) & seeds2[t][rs];

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
		if((pos1 + read_len1) < size_chrom[chrom_id1]){

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
					convert_gen_mer_ta(ta_gen_mer1, pos1, chrom_id1, read_len1);

					unsigned long int ta_check = ta_read1 ^ ta_gen_mer1, ta_check2;
					mism1 = count_mism(ta_check);

					if(mism1 <= mism_thresh){

						convert_gen_mer_ta(ta_gen_mer2, pos2, chr_id2, read_len2);
						ta_check2 = ta_rc_read2 ^ ta_gen_mer2;
						mism2 = count_mism(ta_check2);

						if(mism2 <= mism_thresh){
							unsigned long int ns_gen_mer1 = 0, ns_gen_mer2 =  0;
							convert_gen_mer_cg(cg_gen_mer1, pos1, chrom_id1, read_len1);
							convert_gen_mer_cg(cg_gen_mer2, pos2, chr_id2, read_len2);
							convert_gen_mer_ns(ns_gen_mer1, pos1, chrom_id1, read_len1);
							convert_gen_mer_ns(ns_gen_mer2, pos2, chr_id2, read_len2);

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
								mates_pairs[read_id].chrom_id = chrom_id1 ;
								mates_pairs[read_id].pos1 = pos1;
								mates_pairs[read_id].mism1 = mism1;
								mates_pairs[read_id].pos2 = pos2;
								mates_pairs[read_id].mism2 = mism2;
								mates_pairs[read_id].min_mism = min_mism;
								mates_pairs[read_id].strand = '+';//by the first read
							}//if
							else if(mates_pairs[read_id].chrom_id != chrom_id1 ||
								mates_pairs[read_id].pos1 != pos1 || 
								mates_pairs[read_id].pos2 != pos2 || mates_pairs[read_id].strand != astrand){


								if((min_mism < mates_pairs[read_id].min_mism) || 
									(min_mism == mates_pairs[read_id].min_mism && 
									second_mism < cur_second)){
									//unique
		
									mates_pairs[read_id].mark_unique = 1;
									mates_pairs[read_id].chrom_id = chrom_id1 ;
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
			
			}//if first read fits into genome
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
			if(pos1 + read_len2 < size_chrom[chrom_id1]){

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
							convert_gen_mer_ta(ta_gen_mer2, pos1, chrom_id1, read_len2);


					unsigned long int ta_check2 = ta_read2 ^ ta_gen_mer2, ta_check;
					mism2 = count_mism(ta_check2);

					if(mism2 <= mism_thresh){
						convert_gen_mer_ta(ta_gen_mer1, pos2, chr_id2, read_len1);
						ta_check = ta_rc_read1 ^ ta_gen_mer1;
						mism1 = count_mism(ta_check);

						if(mism1 <= mism_thresh){
							unsigned long int ns_gen_mer1 = 0, ns_gen_mer2 = 0;
							convert_gen_mer_cg( cg_gen_mer2, pos1, chrom_id1, read_len2);
							convert_gen_mer_cg( cg_gen_mer1, pos2, chr_id2, read_len1);
							convert_gen_mer_ns(ns_gen_mer1, pos2, chr_id2, read_len1);
							convert_gen_mer_ns(ns_gen_mer2, pos1, chrom_id1, read_len2);

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
								mates_pairs[read_id].chrom_id = chrom_id1 ;
								mates_pairs[read_id].pos1 = pos2;
								mates_pairs[read_id].mism1 = mism1;
								mates_pairs[read_id].pos2 = pos1;
								mates_pairs[read_id].mism2 = mism2;
								mates_pairs[read_id].min_mism = min_mism;
								mates_pairs[read_id].strand = '-';//by the first read
							}//if
							else if(mates_pairs[read_id].chrom_id != chrom_id1 ||
								mates_pairs[read_id].pos1 != pos2 || 
								mates_pairs[read_id].pos2 != pos1 || mates_pairs[read_id].strand != astrand){

								if((min_mism < mates_pairs[read_id].min_mism) || 
									(min_mism == mates_pairs[read_id].min_mism && 
									second_mism < cur_second)){
									//unique
									mates_pairs[read_id].mark_unique = 1;
									mates_pairs[read_id].chrom_id = chrom_id1 ;
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
				}//if read fits into genome
			start1++;


		if(pair_ambiguous == true){
			break;
		}//if
	}//while
	}//if pair is not ambiguous



}//pair_reads new

////////////////////////////////LONG READs////////////////////////////////////////////
void Anchor::pair_reads_long(unsigned long int read_id, string aread1, string aread2,
				long int &contain_Ns, Stat &astat, short int read_len1, short int read_len2, int t, int rs, int chr )
{



	unsigned long int one = 1;
	long int chrom_id1 = chr;
	long int chrom_id2 = chr;
	long int chr_id2 = chr;

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
		if((pos1 + read_len1) < size_chrom[chrom_id1]){




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

					convert_gen_mer_ta_long(ta_gen_mer1, pos1, chrom_id1, read_len1, ta_size1);


					vector<unsigned long int> ta_check , ta_check2;
					mism1 = count_mism_long(ta_read1, ta_gen_mer1, ta_check, ta_size1);

					if(mism1 <= mism_thresh){
						vector<unsigned long int> ta_gen_mer2;
						
						convert_gen_mer_ta_long(ta_gen_mer2, pos2, chr_id2, read_len2, ta_size2);

						mism2 = count_mism_long(ta_rc_read2, ta_gen_mer2, ta_check2, ta_size2);

						if(mism2 <= mism_thresh){
							vector<unsigned long int> cg_gen_mer1, cg_gen_mer2, ns_gen_mer1, ns_gen_mer2;

							convert_gen_mer_cg_long(cg_gen_mer1, pos1, chrom_id1, read_len1, ta_size1);
							convert_gen_mer_cg_long(cg_gen_mer2, pos2, chr_id2, read_len2, ta_size2);

							convert_gen_mer_ns_long(ns_gen_mer1, pos1, chrom_id1, read_len1, ta_size1);
							convert_gen_mer_ns_long(ns_gen_mer2, pos2, chr_id2, read_len2, ta_size2);

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
								mates_pairs[read_id].chrom_id = chrom_id1 ;
								mates_pairs[read_id].pos1 = pos1;
								mates_pairs[read_id].mism1 = mism1;
								mates_pairs[read_id].pos2 = pos2;
								mates_pairs[read_id].mism2 = mism2;
								mates_pairs[read_id].min_mism = min_mism;
								mates_pairs[read_id].strand = '+';//by the first read
							}//if
							else if(mates_pairs[read_id].chrom_id != chrom_id1 ||
								mates_pairs[read_id].pos1 != pos1 || 
								mates_pairs[read_id].pos2 != pos2 || mates_pairs[read_id].strand != astrand){


								if((min_mism < mates_pairs[read_id].min_mism) || 
									(min_mism == mates_pairs[read_id].min_mism && 
									second_mism < cur_second)){
									//unique
		
									mates_pairs[read_id].mark_unique = 1;
									mates_pairs[read_id].chrom_id = chrom_id1 ;
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

		}//if read1 fits into a ref
		else
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
		end1 = pos_sizes[ind2];

		start2 = 0;
		end2 = pos_sizes[ta_rc1];

		while(start1 < end1){

			if(start2 >= end2)
				break;

			pos1 =  pos_index[ind2][start1];

			if(pos1 + read_len2 < size_chrom[chrom_id1]){


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
						pos2 =  pos_index[ta_rc1][cur];

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
						convert_gen_mer_ta_long(ta_gen_mer2, pos1, chrom_id1, read_len2, ta_size2);


					vector<unsigned long int> ta_check2 , ta_check;
					mism2 = count_mism_long(ta_read2, ta_gen_mer2, ta_check2, ta_size2);

					if(mism2 <= mism_thresh){
						convert_gen_mer_ta_long(ta_gen_mer1, pos2, chr_id2, read_len1, ta_size1);

						mism1 = count_mism_long(ta_rc_read1, ta_gen_mer1, ta_check, ta_size1);

						if(mism1 <= mism_thresh){
							vector<unsigned long int> cg_gen_mer1, cg_gen_mer2, ns_gen_mer1, ns_gen_mer2;

							convert_gen_mer_cg_long( cg_gen_mer2, pos1, chrom_id1, read_len2, ta_size2);
							convert_gen_mer_cg_long( cg_gen_mer1, pos2, chr_id2, read_len1, ta_size1);
						
							convert_gen_mer_ns_long(ns_gen_mer2, pos1, chrom_id1, read_len2, ta_size2);
							convert_gen_mer_ns_long(ns_gen_mer1, pos2, chr_id2, read_len1, ta_size1);

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
								mates_pairs[read_id].chrom_id = chrom_id1 ;
								mates_pairs[read_id].pos1 = pos2;
								mates_pairs[read_id].mism1 = mism1;
								mates_pairs[read_id].pos2 = pos1;
								mates_pairs[read_id].mism2 = mism2;
								mates_pairs[read_id].min_mism = min_mism;
								mates_pairs[read_id].strand = '-';//by the first read
							}//if
							else if(mates_pairs[read_id].chrom_id != chrom_id1 ||
								mates_pairs[read_id].pos1 != pos2 || 
								mates_pairs[read_id].pos2 != pos1 || mates_pairs[read_id].strand != astrand){

								if((min_mism < mates_pairs[read_id].min_mism) || 
									(min_mism == mates_pairs[read_id].min_mism && 
									second_mism < cur_second)){
									//unique
									mates_pairs[read_id].mark_unique = 1;
									mates_pairs[read_id].chrom_id = chrom_id1 ;
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

		}//if read2 fits into a ref
		else
			start1++;

	}//while
	}//if pair is not ambiguous



}//pair_reads_long new



void Anchor::pair_reads_normal_long(unsigned long int read_id, string aread1, string aread2,
				long int &contain_Ns, Stat &astat, short int read_len1, short int read_len2, int t, int rs, int chr )
{


	long int chrom_id1 = chr;
	long int chrom_id2 = chr;
	long int chr_id2 = chr;

	unsigned long int one = 1;

	long int asize,  pos1, pos2, dist, i, j;
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
		if((pos1 + read_len1) < size_chrom[chrom_id1]){




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
				pos2 = pos_index[ta_rc2][cur];
				

				if(dif2 <= pos2){

					pos2 -= dif2;


				dist = pos2 - pos1;
				if(dist >= min_insert && dist <= max_insert ){
					//cg check
					bool cg_mism1 = false, cg_mism2 = false;
					vector<unsigned long int> ta_gen_mer1;

					convert_gen_mer_ta_long(ta_gen_mer1, pos1, chrom_id1, read_len1, ta_size1);


					vector<unsigned long int> ta_check , ta_check2;
					mism1 = count_mism_long(ta_read1, ta_gen_mer1, ta_check, ta_size1);

					if(mism1 <= mism_thresh){
						vector<unsigned long int> ta_gen_mer2;
						
						convert_gen_mer_ta_long(ta_gen_mer2, pos2, chr_id2, read_len2, ta_size2);

						mism2 = count_mism_long(ta_rc_read2, ta_gen_mer2, ta_check2, ta_size2);

						if(mism2 <= mism_thresh){
							vector<unsigned long int> cg_gen_mer1, cg_gen_mer2, ns_gen_mer1, ns_gen_mer2;

							convert_gen_mer_cg_long(cg_gen_mer1, pos1, chrom_id1, read_len1, ta_size1);
							convert_gen_mer_cg_long(cg_gen_mer2, pos2, chr_id2, read_len2, ta_size2);

							convert_gen_mer_ns_long(ns_gen_mer1, pos1, chrom_id1, read_len1, ta_size1);
							convert_gen_mer_ns_long(ns_gen_mer2, pos2, chr_id2, read_len2, ta_size2);

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
								mates_pairs[read_id].chrom_id = chrom_id1 ;
								mates_pairs[read_id].pos1 = pos1;
								mates_pairs[read_id].mism1 = mism1;
								mates_pairs[read_id].pos2 = pos2;
								mates_pairs[read_id].mism2 = mism2;
								mates_pairs[read_id].min_mism = min_mism;
								mates_pairs[read_id].strand = '+';//by the first read
							}//if
							else if(mates_pairs[read_id].chrom_id != chrom_id1 ||
								mates_pairs[read_id].pos1 != pos1 || 
								mates_pairs[read_id].pos2 != pos2 || mates_pairs[read_id].strand != astrand){


								if((min_mism < mates_pairs[read_id].min_mism) || 
									(min_mism == mates_pairs[read_id].min_mism && 
									second_mism < cur_second)){
									//unique
		
									mates_pairs[read_id].mark_unique = 1;
									mates_pairs[read_id].chrom_id = chrom_id1 ;
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

		}//if read1 fits into a ref
		else
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
		end1 = pos_sizes[ind2];

		start2 = 0;
		end2 = pos_sizes[ta_rc1];

		while(start1 < end1){

			if(start2 >= end2)
				break;

			pos1 =  pos_index[ind2][start1];
			if(pos1 + read_len2 < size_chrom[chrom_id1]){


				cur = start2;
				for(; cur != end2 ; cur++){
					pos2 =  pos_index[ta_rc1][cur] ;
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
						convert_gen_mer_ta_long(ta_gen_mer2, pos1, chrom_id1, read_len2, ta_size2);


					vector<unsigned long int> ta_check2 , ta_check;
					mism2 = count_mism_long(ta_read2, ta_gen_mer2, ta_check2, ta_size2);

					if(mism2 <= mism_thresh){
						convert_gen_mer_ta_long(ta_gen_mer1, pos2, chr_id2, read_len1, ta_size1);

						mism1 = count_mism_long(ta_rc_read1, ta_gen_mer1, ta_check, ta_size1);

						if(mism1 <= mism_thresh){
							vector<unsigned long int> cg_gen_mer1, cg_gen_mer2, ns_gen_mer1, ns_gen_mer2;

							convert_gen_mer_cg_long( cg_gen_mer2, pos1, chrom_id1, read_len2, ta_size2);
							convert_gen_mer_cg_long( cg_gen_mer1, pos2, chr_id2, read_len1, ta_size1);
						
							convert_gen_mer_ns_long(ns_gen_mer2, pos1, chrom_id1, read_len2, ta_size2);
							convert_gen_mer_ns_long(ns_gen_mer1, pos2, chr_id2, read_len1, ta_size1);

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
								mates_pairs[read_id].chrom_id = chrom_id1 ;
								mates_pairs[read_id].pos1 = pos2;
								mates_pairs[read_id].mism1 = mism1;
								mates_pairs[read_id].pos2 = pos1;
								mates_pairs[read_id].mism2 = mism2;
								mates_pairs[read_id].min_mism = min_mism;
								mates_pairs[read_id].strand = '-';//by the first read
							}//if
							else if(mates_pairs[read_id].chrom_id != chrom_id1 ||
								mates_pairs[read_id].pos1 != pos2 || 
								mates_pairs[read_id].pos2 != pos1 || mates_pairs[read_id].strand != astrand){

								if((min_mism < mates_pairs[read_id].min_mism) || 
									(min_mism == mates_pairs[read_id].min_mism && 
									second_mism < cur_second)){
									//unique
									mates_pairs[read_id].mark_unique = 1;
									mates_pairs[read_id].chrom_id = chrom_id1 ;
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

		}//if read2 fits into a ref
		else
			start1++;

	}//while
	}//if pair is not ambiguous



}//pair_reads_normal_long new




void Anchor::pair_reads_singles_long(unsigned long int read_id, string aread1,
				long int &contain_Ns, Stat &astat, short int read_len1,  int t, int rs, int chr )
{


	long int chrom_id1 = chr;

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
	convert_lmer_c_long(aread1, c_read1, read_len1);

  //rev-complement of c-reads

	int dif1 = read_len1 - len_mer;

	//we need the first len_mer bases for each read
	if(ta_size1 > 1){

		if(dif_len_mer > 0){
			cg_lmer1 = (c_read1[0] >> dif_len_mer) & seeds[t][rs];
			ta_lmer1 = (ta_read1[0] >> dif_len_mer) & seeds[t][rs];
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
		if((pos1 + read_len1) < size_chrom[chrom_id1]){

			//cg check
			bool cg_mism1 = false;
			vector<unsigned long int> ta_gen_mer1;

			convert_gen_mer_ta_long(ta_gen_mer1, pos1, chrom_id1, read_len1, ta_size1);


			vector<unsigned long int> ta_check;
			mism1 = count_mism_long(ta_read1, ta_gen_mer1, ta_check, ta_size1);


			if(mism1 <= mism_thresh){

				vector<unsigned long int> cg_gen_mer1, ns_gen_mer1;

				convert_gen_mer_cg_long(cg_gen_mer1, pos1, chrom_id1, read_len1, ta_size1);

				convert_gen_mer_ns_long(ns_gen_mer1, pos1, chrom_id1, read_len1, ta_size1);

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
					mates_pairs[read_id].chrom_id = chrom_id1 ;
					mates_pairs[read_id].pos1 = pos1;
					mates_pairs[read_id].mism1 = mism1;
					mates_pairs[read_id].strand = '+';//by the first read
				}//if
				else if(mates_pairs[read_id].chrom_id != chrom_id1 ||
					mates_pairs[read_id].pos1 != pos1 || 
					mates_pairs[read_id].strand != astrand){
						if(mism1 < mates_pairs[read_id].mism1){
						//unique
							mates_pairs[read_id].mark_unique = 1;
							mates_pairs[read_id].chrom_id = chrom_id1 ;
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

		}//if read1 fits into a ref
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
			if(pos1 + read_len1 < size_chrom[chrom_id1]){

			if(dif1 <= pos1){
				pos1 -= dif1;

				vector<unsigned long int>  ta_gen_mer1;

				bool cg_mism1 = false;
				vector<unsigned long int>  ta_check;
				convert_gen_mer_ta_long(ta_gen_mer1, pos1, chrom_id1, read_len1, ta_size1);
				mism1 = count_mism_long(ta_rc_read1, ta_gen_mer1, ta_check, ta_size1);

				if(mism1 <= mism_thresh){
					vector<unsigned long int> cg_gen_mer1, ns_gen_mer1;
					convert_gen_mer_cg_long( cg_gen_mer1, pos1, chrom_id1, read_len1, ta_size1);
					convert_gen_mer_ns_long(ns_gen_mer1, pos1, chrom_id1, read_len1, ta_size1);

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
						mates_pairs[read_id].chrom_id = chrom_id1 ;
						mates_pairs[read_id].pos1 = pos1;
						mates_pairs[read_id].mism1 = mism1;
						mates_pairs[read_id].strand = '-';//by the first read
					}//if
					else if(mates_pairs[read_id].chrom_id != chrom_id1 ||
						mates_pairs[read_id].pos1 != pos1 || 
						mates_pairs[read_id].strand != astrand){
						if(mism1 < mates_pairs[read_id].mism1){
						//unique
							mates_pairs[read_id].mark_unique = 1;
							mates_pairs[read_id].chrom_id = chrom_id1 ;
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
		}//if fits into ref
		start1++;

	}//while
	}//if pair is not ambiguous



}//pair_reads_singles_long new



void Anchor::pair_reads_singlesA_long(unsigned long int read_id, string aread1,
				long int &contain_Ns, Stat &astat, short int read_len1, int t, int rs , int chr)
{

	long int chrom_id1 = chr;

	unsigned long int one = 1;

	long int asize,    pos1, dist, i, j;
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
	convert_lmer_c_long(aread1, c_read1, read_len1);

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


	//rev-complement of g-reads


	//process forward read1 and rev read2
	//start and end of list of genome positions for forw1
	long int start1 = 0;
	long int end1 = pos_sizes[ind1];

	//start and end of list of genome positions for rev-compl2

	while(start1 < end1){

		pos1 = pos_index[ind1][start1];
		if((pos1 + read_len1) < size_chrom[chrom_id1]){

			//cg check
			bool cg_mism1 = false;
			vector<unsigned long int> ta_gen_mer1;

			convert_gen_mer_ta_long(ta_gen_mer1, pos1, chrom_id1, read_len1, ta_size1);


			vector<unsigned long int> ta_check;
			mism1 = count_mism_long(ta_read1, ta_gen_mer1, ta_check, ta_size1);


			if(mism1 <= mism_thresh){

				vector<unsigned long int> cg_gen_mer1, ns_gen_mer1;

				convert_gen_mer_cg_long(cg_gen_mer1, pos1, chrom_id1, read_len1, ta_size1);

				convert_gen_mer_ns_long(ns_gen_mer1, pos1, chrom_id1, read_len1, ta_size1);

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
					mates_pairs[read_id].chrom_id = chrom_id1 ;
					mates_pairs[read_id].pos1 = pos1;
					mates_pairs[read_id].mism1 = mism1;
					mates_pairs[read_id].strand = '+';//by the first read
				}//if
				else if(mates_pairs[read_id].chrom_id != chrom_id1 ||
					mates_pairs[read_id].pos1 != pos1 || 
					mates_pairs[read_id].strand != astrand){
						if(mism1 < mates_pairs[read_id].mism1){
						//unique
							mates_pairs[read_id].mark_unique = 1;
							mates_pairs[read_id].chrom_id = chrom_id1 ;
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

		}//if read1 fits into a ref
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
			if(pos1 + read_len1 < size_chrom[chrom_id1]){

			if(dif1 <= pos1 ){
				pos1 -= dif1;

				vector<unsigned long int>  ta_gen_mer1;

				bool cg_mism1 = false;
				vector<unsigned long int>  ta_check;
				convert_gen_mer_ta_long(ta_gen_mer1, pos1, chrom_id1, read_len1, ta_size1);
				mism1 = count_mism_long(ta_rc_read1, ta_gen_mer1, ta_check, ta_size1);

				if(mism1 <= mism_thresh){
					vector<unsigned long int> cg_gen_mer1,ns_gen_mer1;
					convert_gen_mer_cg_long( cg_gen_mer1, pos1, chrom_id1, read_len1, ta_size1);
					convert_gen_mer_ns_long(ns_gen_mer1, pos1, chrom_id1, read_len1, ta_size1);

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
						mates_pairs[read_id].chrom_id = chrom_id1 ;
						mates_pairs[read_id].pos1 = pos1;
						mates_pairs[read_id].mism1 = mism1;
						mates_pairs[read_id].strand = '-';//by the first read
					}//if
					else if(mates_pairs[read_id].chrom_id != chrom_id1 ||
						mates_pairs[read_id].pos1 != pos1 || 
						mates_pairs[read_id].strand != astrand){
						if(mism1 < mates_pairs[read_id].mism1){
						//unique
							mates_pairs[read_id].mark_unique = 1;
							mates_pairs[read_id].chrom_id = chrom_id1 ;
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

		}//if read2 fits into a ref: this is dif
		}//if fits into ref
		start1++;

	}//while
	}//if pair is not ambiguous



}//pair_reads_singles_long new


void Anchor::pair_reads_normal_singles_long(unsigned long int read_id, string aread1,
				long int &contain_Ns, Stat &astat, short int read_len1, int t, int rs, int chr )
{


	long int chrom_id1 = chr;

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
	convert_lmer_c_long(aread1, c_read1, read_len1);

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
		if((pos1 + read_len1) < size_chrom[chrom_id1]){

			//cg check
			bool cg_mism1 = false;
			vector<unsigned long int> ta_gen_mer1;

			convert_gen_mer_ta_long(ta_gen_mer1, pos1, chrom_id1, read_len1, ta_size1);


			vector<unsigned long int> ta_check;
			mism1 = count_mism_long(ta_read1, ta_gen_mer1, ta_check, ta_size1);


			if(mism1 <= mism_thresh){

				vector<unsigned long int> cg_gen_mer1, ns_gen_mer1;

				convert_gen_mer_cg_long(cg_gen_mer1, pos1, chrom_id1, read_len1, ta_size1);

				convert_gen_mer_ns_long(ns_gen_mer1, pos1, chrom_id1, read_len1, ta_size1);

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
					mates_pairs[read_id].chrom_id = chrom_id1 ;
					mates_pairs[read_id].pos1 = pos1;
					mates_pairs[read_id].mism1 = mism1;
					mates_pairs[read_id].strand = '+';//by the first read
				}//if
				else if(mates_pairs[read_id].chrom_id != chrom_id1 ||
					mates_pairs[read_id].pos1 != pos1 || 
					mates_pairs[read_id].strand != astrand){
						if(mism1 < mates_pairs[read_id].mism1){
						//unique
							mates_pairs[read_id].mark_unique = 1;
							mates_pairs[read_id].chrom_id = chrom_id1 ;
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

		}//if read1 fits into a ref
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
			if(pos1 + read_len1 < size_chrom[chrom_id1]){

			if(dif1 <= pos1){
				pos1 -= dif1;

				vector<unsigned long int>  ta_gen_mer1;

				bool cg_mism1 = false;
				vector<unsigned long int>  ta_check;
				convert_gen_mer_ta_long(ta_gen_mer1, pos1, chrom_id1, read_len1, ta_size1);
				mism1 = count_mism_long(ta_rc_read1, ta_gen_mer1, ta_check, ta_size1);

				if(mism1 <= mism_thresh){
					vector<unsigned long int> cg_gen_mer1, ns_gen_mer1;
					convert_gen_mer_cg_long( cg_gen_mer1, pos1, chrom_id1, read_len1, ta_size1);
					convert_gen_mer_ns_long(ns_gen_mer1, pos1, chrom_id1, read_len1, ta_size1);

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
						mates_pairs[read_id].chrom_id = chrom_id1 ;
						mates_pairs[read_id].pos1 = pos1;
						mates_pairs[read_id].mism1 = mism1;
						mates_pairs[read_id].strand = '-';//by the first read
					}//if
					else if(mates_pairs[read_id].chrom_id != chrom_id1 ||
						mates_pairs[read_id].pos1 != pos1 || 
						mates_pairs[read_id].strand != astrand){
						if(mism1 < mates_pairs[read_id].mism1){
						//unique
							mates_pairs[read_id].mark_unique = 1;
							mates_pairs[read_id].chrom_id = chrom_id1 ;
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
		}//if fits
		start1++;

	}//while
	}//if pair is not ambiguous



}//pair_reads_singles_long new


/////////////////////////////// END LONG READS///////////////////////////////////////


void Anchor::pair_reads_normal(unsigned long int read_id, string aread1, string aread2,
				long int &contain_Ns, Stat &astat, short int read_len1, short int read_len2, int t, int rs , int chr)
{
	long int chrom_id1 = chr;
	long int chrom_id2 = chr;
	long int chr_id2 = chr;

	unsigned long int one = 1;

	long int asize,  pos1, pos2, dist, i;
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
	convert_lmer_c(aread1, c_read1,  read_len1);
	convert_lmer_c(aread2, c_read2,  read_len2);

	int dif1 = read_len1 - len_mer;
	int dif2 = read_len2 - len_mer;

	ta_lmer1 = (ta_read1 >> dif1) & seeds[t][rs];
	cg_lmer1 = (c_read1 >> dif1) & seeds[t][rs];

	ta_lmer2 = (ta_read2 >> dif2) & seeds2[t][rs];
	cg_lmer2 = (c_read2 >> dif2) & seeds2[t][rs];

	ind1 = hash2c(ta_lmer1, cg_lmer1 );
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
		if(pos1 + read_len1 < size_chrom[chrom_id1]){

			cur = start2;
			for(; cur != end2 ; cur++){
				pos2 = pos_index[ta_rc2][cur] ;
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
				pos2 = pos_index[ta_rc2][cur];
				if(dif2 <= pos2){

				pos2 -= dif2;
				dist = pos2 - pos1;
				if(dist >= min_insert && dist <= max_insert ){
					//cg check
					bool cg_mism1 = false, cg_mism2 = false;
					convert_gen_mer_ta(ta_gen_mer1, pos1, chrom_id1, read_len1);

					unsigned long int ta_check = ta_read1 ^ ta_gen_mer1, ta_check2;
					mism1 = count_mism(ta_check);

					if(mism1 <= mism_thresh){

						convert_gen_mer_ta(ta_gen_mer2, pos2, chr_id2, read_len2);
						ta_check2 = ta_rc_read2 ^ ta_gen_mer2;
						mism2 = count_mism(ta_check2);

						if(mism2 <= mism_thresh){
							convert_gen_mer_cg(cg_gen_mer1, pos1, chrom_id1, read_len1);
							convert_gen_mer_cg(cg_gen_mer2, pos2, chr_id2, read_len2);

							unsigned long int ns_gen_mer1 = 0, ns_gen_mer2 = 0;
							convert_gen_mer_ns(ns_gen_mer1, pos1, chrom_id1, read_len1);
							convert_gen_mer_ns(ns_gen_mer2, pos2, chr_id2, read_len2);
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
								mates_pairs[read_id].chrom_id = chrom_id1 ;
								mates_pairs[read_id].pos1 = pos1;
								mates_pairs[read_id].mism1 = mism1;
								mates_pairs[read_id].pos2 = pos2;
								mates_pairs[read_id].mism2 = mism2;
								mates_pairs[read_id].min_mism = min_mism;
								mates_pairs[read_id].strand = '+';//by the first read
							}//if
							else if(mates_pairs[read_id].chrom_id != chrom_id1 ||
								mates_pairs[read_id].pos1 != pos1 || 
								mates_pairs[read_id].pos2 != pos2 || mates_pairs[read_id].strand != astrand){

								if((min_mism < mates_pairs[read_id].min_mism) || 
									(min_mism == mates_pairs[read_id].min_mism && 
									second_mism < cur_second)){
									//unique
									mates_pairs[read_id].mark_unique = 1;
									mates_pairs[read_id].chrom_id = chrom_id1 ;
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
			}//if fits into genome
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
			if(pos1 + read_len2 < size_chrom[chrom_id1]){

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
						pos2 =  pos_index[ta_rc1][cur] ;
						if(dif1 <= pos2){

							pos2 -= dif1;
						dist = pos2 - pos1;
						if(dist >= min_insert && dist <= max_insert ){
					//cg check
					
							bool cg_mism1 = false, cg_mism2 = false;
							convert_gen_mer_ta(ta_gen_mer2, pos1, chrom_id1, read_len2);

					unsigned long int ta_check2 = ta_read2 ^ ta_gen_mer2, ta_check;
					mism2 = count_mism(ta_check2);

					if(mism2 <= mism_thresh){
						convert_gen_mer_ta(ta_gen_mer1, pos2, chr_id2, read_len1);
						ta_check = ta_rc_read1 ^ ta_gen_mer1;
						mism1 = count_mism(ta_check);

						if(mism1 <= mism_thresh){
							convert_gen_mer_cg( cg_gen_mer2, pos1, chrom_id1, read_len2);
							convert_gen_mer_cg( cg_gen_mer1, pos2, chr_id2, read_len1);

							unsigned long int ns_gen_mer1 = 0, ns_gen_mer2 = 0;
							convert_gen_mer_ns(ns_gen_mer1, pos2, chr_id2, read_len1);
							convert_gen_mer_ns(ns_gen_mer2, pos1, chrom_id1 , read_len2);
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
								mates_pairs[read_id].chrom_id = chrom_id1 ;
								mates_pairs[read_id].pos1 = pos2;
								mates_pairs[read_id].mism1 = mism1;
								mates_pairs[read_id].pos2 = pos1;
								mates_pairs[read_id].mism2 = mism2;
								mates_pairs[read_id].min_mism = min_mism;
								mates_pairs[read_id].strand = '-';//by the first read
							}//if
							else if(mates_pairs[read_id].chrom_id != chrom_id1 ||
								mates_pairs[read_id].pos1 != pos2 || 
								mates_pairs[read_id].pos2 != pos1 || mates_pairs[read_id].strand != astrand){

								if((min_mism < mates_pairs[read_id].min_mism) || 
									(min_mism == mates_pairs[read_id].min_mism && 
									second_mism < cur_second)){
									//unique
									mates_pairs[read_id].mark_unique = 1;
									mates_pairs[read_id].chrom_id = chrom_id1 ;
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
				}//if fits into genome
			start1++;

		if(pair_ambiguous == true){
			break;
		}//if
	}//while
	}//if pair is not ambiguous




}//pair_reads_normal

//START SINGLES
void Anchor::pair_reads_singles(unsigned long int read_id, string aread1,
				long int &contain_Ns, Stat &astat, short int read_len1,  int t, int rs, int chr )
{

	long int chrom_id1 = chr;

	unsigned long int one = 1, mism1, ta_lmer1, cg_lmer1;

	long int asize,  pos1, dist, i;
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
		if(pos1 + read_len1 < size_chrom[chrom_id1]){

		bool cg_mism1 = false;
		convert_gen_mer_ta(ta_gen_mer1, pos1, chrom_id1, read_len1);

		unsigned long int ta_check = ta_read1 ^ ta_gen_mer1;
		mism1 = count_mism(ta_check);

		if(mism1 <= mism_thresh){
			unsigned long int ns_gen_mer1 = 0;
			convert_gen_mer_ns(ns_gen_mer1, pos1, chrom_id1, read_len1);
			ta_check |= ns_gen_mer1;

			convert_gen_mer_cg(cg_gen_mer1, pos1, chrom_id1, read_len1);
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
					mates_pairs[read_id].chrom_id = chrom_id1 ;
					mates_pairs[read_id].pos1 = pos1;
					mates_pairs[read_id].strand = '+';//by the first read
					mates_pairs[read_id].mism1 = mism1;
				}//if
				else if(mates_pairs[read_id].chrom_id != chrom_id1 ||
					mates_pairs[read_id].pos1 != pos1 || mates_pairs[read_id].strand != astrand){
					if(mism1 < mates_pairs[read_id].mism1){
						//unique
						mates_pairs[read_id].mark_unique = 1;
						mates_pairs[read_id].chrom_id = chrom_id1 ;
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
		}//if fits into genome
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

			pos1 = pos_index[ta_rc1][start1];
			if(pos1 + read_len1 < size_chrom[chrom_id1]){

			if(dif1 <= pos1){
				pos1 -= dif1;

			bool cg_mism1 = false;

			convert_gen_mer_ta(ta_gen_mer1, pos1, chrom_id1, read_len1);
			unsigned long int ta_check = ta_rc_read1 ^ ta_gen_mer1;
			mism1 = count_mism(ta_check);
			if(mism1 <= mism_thresh){
				unsigned long int ns_gen_mer1 = 0;
				convert_gen_mer_ns(ns_gen_mer1, pos1, chrom_id1, read_len1);
				ta_check |= ns_gen_mer1;

				convert_gen_mer_cg(cg_gen_mer1, pos1, chrom_id1, read_len1);
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
					mates_pairs[read_id].chrom_id = chrom_id1 ;
					mates_pairs[read_id].pos1 = pos1;
					mates_pairs[read_id].strand = '-';//by the first read
					mates_pairs[read_id].mism1 = mism1;

				}//if
				else if(mates_pairs[read_id].chrom_id != chrom_id1 ||
					mates_pairs[read_id].pos1 != pos1 || mates_pairs[read_id].strand != astrand){
					if(mism1 < mates_pairs[read_id].mism1){
						//unique
						mates_pairs[read_id].mark_unique = 1;
						mates_pairs[read_id].chrom_id = chrom_id1 ;
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
			}//if fits into genome
			if(pair_ambiguous == true)
				break;
			start1++;
		}//while

	}//if pair is not ambiguous

}//pair_reads_singles new

//cascade ambiguous bs singles


///////////////////////////////////////////////////////////
//A-rich strands
void Anchor::pair_reads_singlesA(unsigned long int read_id, string aread1,
				long int &contain_Ns, Stat &astat, short int read_len1, int t, int rs, int chr )
{

	long int chrom_id1 = chr;

	unsigned long int one = 1;
	short int mism1;
	unsigned long int ta_lmer1, cg_lmer1;
	long int asize,   pos1, dist, i;
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
		if(pos1 + read_len1 < size_chrom[chrom_id1]){

		bool cg_mism1 = false;
		convert_gen_mer_ta(ta_gen_mer1, pos1, chrom_id1, read_len1);

		unsigned long int ta_check = ta_read1 ^ ta_gen_mer1;
		mism1 = count_mism(ta_check);

		if(mism1 <= mism_thresh){
			convert_gen_mer_cg(cg_gen_mer1, pos1, chrom_id1, read_len1);
			tc_mism1 = (ta_read1 & (~c_read1)) & cg_gen_mer1 ;
			ta_check |= tc_mism1;
			unsigned long int ns_gen_mer1 = 0;
			convert_gen_mer_ns(ns_gen_mer1, pos1, chrom_id1, read_len1);
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
					mates_pairs[read_id].chrom_id = chrom_id1 ;
					mates_pairs[read_id].pos1 = pos1;
					mates_pairs[read_id].strand = '+';//by the first read
					mates_pairs[read_id].mism1 = mism1;
				}//if
				else if(mates_pairs[read_id].chrom_id != chrom_id1 ||
					mates_pairs[read_id].pos1 != pos1 || mates_pairs[read_id].strand != astrand){
					if(mism1 < mates_pairs[read_id].mism1){
						//unique
						mates_pairs[read_id].mark_unique = 1;
						mates_pairs[read_id].chrom_id = chrom_id1 ;
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
		}//if fits into genome
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
			if(pos1 + read_len1 < size_chrom[chrom_id1]){
			bool cg_mism1 = false;

			if(dif1 <= pos1){
				pos1 -= dif1;
			convert_gen_mer_ta(ta_gen_mer1, pos1, chrom_id1, read_len1);

			unsigned long int ta_check = ta_rc_read1 ^ ta_gen_mer1;
			mism1 = count_mism(ta_check);
			if(mism1 <= mism_thresh){
				unsigned long int ns_gen_mer1 = 0;
				convert_gen_mer_ns(ns_gen_mer1, pos1, chrom_id1, read_len1);
				ta_check |= ns_gen_mer1;

				convert_gen_mer_cg(cg_gen_mer1, pos1, chrom_id1, read_len1);

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
					mates_pairs[read_id].chrom_id = chrom_id1 ;
					mates_pairs[read_id].pos1 = pos1;
					mates_pairs[read_id].strand = '-';//by the first read
					mates_pairs[read_id].mism1 = mism1;
				}//if
				else if(mates_pairs[read_id].chrom_id != chrom_id1 ||
					mates_pairs[read_id].pos1 != pos1 || mates_pairs[read_id].strand != astrand){
					if(mism1 < mates_pairs[read_id].mism1){
						//unique
						mates_pairs[read_id].mark_unique = 1;
						mates_pairs[read_id].chrom_id = chrom_id1 ;
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
			}//if fits into genome
			if(pair_ambiguous == true)
				break;
			start1++;
		}//while

	}//if pair is not ambiguous




}//pair_reads_singles new

////////////////////////////////////////////////////////////
//Normal
void Anchor::pair_reads_normal_singles(unsigned long int read_id, string aread1,
				long int &contain_Ns, Stat &astat, short int read_len1, int t, int rs , int chr)
{

	long int chrom_id1 = chr;

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
	convert_lmer_c(aread1, c_read1, read_len1);

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
		if(pos1 + read_len1 < size_chrom[chrom_id1]){
		bool cg_mism1 = false;
		convert_gen_mer_ta(ta_gen_mer1, pos1, chrom_id1, read_len1);

		unsigned long int ta_check = ta_read1 ^ ta_gen_mer1;
		mism1 = count_mism(ta_check);

		if(mism1 <= mism_thresh){
			unsigned long int ns_gen_mer1 = 0;
			convert_gen_mer_ns(ns_gen_mer1, pos1, chrom_id1, read_len1);
			ta_check |= ns_gen_mer1;


			convert_gen_mer_cg(cg_gen_mer1, pos1, chrom_id1, read_len1);
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
					mates_pairs[read_id].chrom_id = chrom_id1 ;
					mates_pairs[read_id].pos1 = pos1;
					mates_pairs[read_id].strand = '+';//by the first read
					mates_pairs[read_id].mism1 = mism1;
				}//if
				else if(mates_pairs[read_id].chrom_id != chrom_id1 ||
					mates_pairs[read_id].pos1 != pos1 || mates_pairs[read_id].strand != astrand){
					if(mism1 < mates_pairs[read_id].mism1){
						//unique
						mates_pairs[read_id].mark_unique = 1;
						mates_pairs[read_id].chrom_id = chrom_id1 ;
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
		}//if fits into genome
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

			pos1 = pos_index[ta_rc1][start1] ;
			if(pos1 + read_len1 < size_chrom[chrom_id1]){

			bool cg_mism1 = false;

			if(dif1 <= pos1){
				pos1 -= dif1;

			convert_gen_mer_ta(ta_gen_mer1, pos1, chrom_id1, read_len1);

			unsigned long int ta_check = ta_rc_read1 ^ ta_gen_mer1;
			mism1 = count_mism(ta_check);
			if(mism1 <= mism_thresh){
				unsigned long int ns_gen_mer1 = 0;
				convert_gen_mer_ns(ns_gen_mer1, pos1, chrom_id1, read_len1);
				ta_check |= ns_gen_mer1;

				convert_gen_mer_cg(cg_gen_mer1, pos1, chrom_id1, read_len1);
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
					mates_pairs[read_id].chrom_id = chrom_id1 ;
					mates_pairs[read_id].pos1 = pos1;
					mates_pairs[read_id].strand = '-';//by the first read
					mates_pairs[read_id].mism1 = mism1;
				}//if
				else if(mates_pairs[read_id].chrom_id != chrom_id1 ||
					mates_pairs[read_id].pos1 != pos1 || mates_pairs[read_id].strand != astrand){
					if(mism1 < mates_pairs[read_id].mism1){
						//unique
						mates_pairs[read_id].mark_unique = 1;
						mates_pairs[read_id].chrom_id = chrom_id1 ;
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
			}//if fits into genome
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
		case 'P': prefix = argv[++i]; break;
		case 'i': min_insert = atoi(argv[++i]); break;
		case 'a': max_insert = atoi(argv[++i]); break;
		case 'A': A_rich = true; break;
		case 'u': unmapped = true; break;
		case 'm': mism_thresh = atoi(argv[++i]); break;
		case 'M': out_amb = true; break;
		case 'f': first_bits = atoi(argv[++i]); break;
		case 'S': speed_mode = "1"; break;
		case 'h': usage(); break;
		default: cerr << "ERROR: Invalid char=" << argv[i][1] << "=this" << endl; usage(); exit(0); break;
		}
	}//for i
	return res;
}//parse_opt()

void Anchor::print(Stat &astat){

	astat.mapped_pairs = 0;
	astat.ambiguous = 0;
	astat.invalid = 0;

	long int orig_start1, orig_start2, st1, en1, st2, en2;
		ofstream out_amb1, out_amb2;
		string output_amb1 = reads_file1 + ".amb";
		string output_amb2 = reads_file2 + ".amb";


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
				astat.mapped_pairs++;
				map<unsigned int, string>::iterator afind = chroms_names.find(mates_pairs[read_id].chrom_id);
				string chr_name = afind->second;
				if(mates_pairs[read_id].strand == '+')
				{
					orig_start1 = mates_pairs[read_id].pos1 - st1;
					orig_start2 = mates_pairs[read_id].pos2 - en2;
				}
				else{
					orig_start1 = mates_pairs[read_id].pos1 - en1;
					orig_start2 = mates_pairs[read_id].pos2 - st2;
				}
				if(orig_start1 < 0)
					orig_start1 = mates_pairs[read_id].pos1 ;
				if(orig_start2 < 0)
					orig_start2 = mates_pairs[read_id].pos2 ;

				out_pairs << read_id << "\t" << aread1 << "\t" << aread2 << "\t"
						<< chr_name << "\t" << mates_pairs[read_id].strand <<"\t"
						<< mates_pairs[read_id].pos1 << "\t" 
						<< mates_pairs[read_id].pos2 << "\t"
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
				astat.mapped_pairs++;
				map<unsigned int, string>::iterator afind = chroms_names.find(mates_pairs[read_id].chrom_id);
				string chr_name = afind->second;
				if(mates_pairs[read_id].strand == '+')
					orig_start1 = mates_pairs[read_id].pos1 - st1;
				else
					orig_start1 = mates_pairs[read_id].pos1 - en1;
				if(orig_start1 < 0)
					orig_start1 = mates_pairs[read_id].pos1 ;
				out_pairs << read_id << "\t" << aread1 << "\t" 
						<< chr_name << "\t" << mates_pairs[read_id].strand <<"\t"
						<< mates_pairs[read_id].pos1 << "\t" << mates_pairs[read_id].mism1 << "\t" << orig_start1 << endl;
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
}
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
	long int st1, en1, st2, en2;
	int t = 0, rs = 0;

	int seeds_size = anc.gen_seeds.size();

	long int num_chroms = anc.num_chroms;

	for(u = 0; u < num_chroms; u++){

		anc.cur_chrom = u;
		if(anc.prefix.length() == 0){
			anc.read_genome(u);
			anc.find_pos_sizes(u);
			anc.make_pos_index(u);
		}
		else{
			anc.read_genome2(u);
			anc.find_pos_sizes2(u);
			anc.make_pos_index2(u);
		
		}//else

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
				if(anc.mates_pairs[read_id].mark_unique == 0 || anc.mates_pairs[read_id].mark_unique == 1){//unmapped to any locations: map always
					if(anc.reads_len1[read_id] > 64 || anc.reads_len2[read_id] > 64)
						anc.pair_reads_long(read_id, read, read2, contain_Ns, astat,  anc.reads_len1[read_id], anc.reads_len2[read_id], t, rs, u);
					else
						anc.pair_reads(read_id, read, read2, contain_Ns, astat,  anc.reads_len1[read_id], anc.reads_len2[read_id], t, rs, u);

				}
				else if((anc.mates_pairs[read_id].mark_unique > 1) &&
					(anc.mates_pairs[read_id].mism1 > 0 || anc.mates_pairs[read_id].mism2 > 0)){
					//with all other combinations of mismatches for mates,
					//if it is ambiguous once, it is always ambiguous
					if(anc.reads_len1[read_id] > 64 || anc.reads_len2[read_id] > 64)
						anc.pair_reads_long(read_id, read, read2, contain_Ns, astat,  anc.reads_len1[read_id], anc.reads_len2[read_id], t, rs, u);
					else
						anc.pair_reads(read_id, read, read2, contain_Ns, astat,  anc.reads_len1[read_id], anc.reads_len2[read_id], t, rs, u);

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
				if(anc.mates_pairs[read_id].mark_unique == 0 || anc.mates_pairs[read_id].mark_unique == 1){//unmapped to any locations: map always
					if(anc.reads_len1[read_id] > 64 || anc.reads_len2[read_id] > 64)
						anc.pair_reads_normal_long(read_id, read, read2, contain_Ns, astat,  anc.reads_len1[read_id], anc.reads_len2[read_id], t, rs, u);
					else
						anc.pair_reads_normal(read_id, read, read2, contain_Ns, astat,  anc.reads_len1[read_id], anc.reads_len2[read_id], t, rs, u);

				}
				else if((anc.mates_pairs[read_id].mark_unique > 1) &&
					(anc.mates_pairs[read_id].mism1 > 0 || anc.mates_pairs[read_id].mism2 > 0)){
					if(anc.reads_len1[read_id] > 64 || anc.reads_len2[read_id] > 64)
						anc.pair_reads_normal_long(read_id, read, read2, contain_Ns, astat,  anc.reads_len1[read_id], anc.reads_len2[read_id], t, rs, u);
					else
						anc.pair_reads_normal(read_id, read, read2, contain_Ns, astat,  anc.reads_len1[read_id], anc.reads_len2[read_id], t, rs, u);
				}//if not perfect match


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
					if(anc.mates_pairs[read_id].mark_unique == 0 || anc.mates_pairs[read_id].mark_unique == 1)
					{//if have not been mapped, try to map
						if(anc.A_rich == false){
							if(anc.reads_len1[read_id] > 64)
								anc.pair_reads_singles_long(read_id, read, contain_Ns, astat, anc.reads_len1[read_id], t, rs, u);
							else	
								anc.pair_reads_singles(read_id, read, contain_Ns, astat, anc.reads_len1[read_id], t, rs, u);
						}
						else{
							if(anc.reads_len1[read_id] > 64)
								anc.pair_reads_singlesA_long(read_id, read, contain_Ns, astat, anc.reads_len1[read_id], t, rs, u);
							else	
								anc.pair_reads_singlesA(read_id, read, contain_Ns, astat, anc.reads_len1[read_id], t, rs, u);

					
						}//A-rich is true
					
					}//if not mapped yet
					else if(anc.mates_pairs[read_id].mism1 > 0){
						if(anc.A_rich == false){
							if(anc.reads_len1[read_id] > 64)
								anc.pair_reads_singles_long(read_id, read, contain_Ns, astat, anc.reads_len1[read_id], t, rs, u);
							else	
								anc.pair_reads_singles(read_id, read, contain_Ns, astat, anc.reads_len1[read_id], t, rs, u);
						}
						else{
							if(anc.reads_len1[read_id] > 64)
								anc.pair_reads_singlesA_long(read_id, read, contain_Ns, astat, anc.reads_len1[read_id], t, rs, u);
							else	
								anc.pair_reads_singlesA(read_id, read, contain_Ns, astat, anc.reads_len1[read_id], t, rs, u);

					
						}//A-rich is true
					
					}//if mism > 0
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
					if(anc.mates_pairs[read_id].mark_unique == 0 || anc.mates_pairs[read_id].mark_unique == 1){
							if(anc.reads_len1[read_id] > 64)
								anc.pair_reads_normal_singles_long(read_id, read, contain_Ns, astat, anc.reads_len1[read_id], t, rs, u);
							else	
								anc.pair_reads_normal_singles(read_id, read, contain_Ns, astat, anc.reads_len1[read_id], t, rs, u);

				
					}//if not mapped with 0-mism
					else if(anc.mates_pairs[read_id].mism1 > 0){
							if(anc.reads_len1[read_id] > 64)
								anc.pair_reads_singlesA_long(read_id, read, contain_Ns, astat, anc.reads_len1[read_id], t, rs, u);
							else	
								anc.pair_reads_singlesA(read_id, read, contain_Ns, astat, anc.reads_len1[read_id], t, rs, u);
					
					}//if mism > 0
					
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

		anc.release_space(u);
	}//for u, all chroms

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


