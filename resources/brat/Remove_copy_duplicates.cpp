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
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
using namespace std;

class Point
{

public:
	//constructors of this class
	Point(long int xx, long int yy, long int zz) : x(xx), y(yy), z(zz) {};
	Point() : x(0), y(0), z(0) {};
	
	
	void operator =(Point a) { x = a.x; y = a.y; z = a.z; };
	//x is index, y is tails
	bool operator <(const Point& a) const
	{ 
		if(y < a.y) return true; 
		else return false;
	}


	long int x;
	long int y;
	long int z;

};
class Cov{
public:
	Cov(unsigned short int x, unsigned short int y, unsigned short int z, unsigned short int w) : a(x), c(y), g(z), t(w) {};
	Cov() : a(0), c(0), g(0), t(0) {};

	void operator =(Cov s){ a = s.a; c = s.c; g = s.g; t = s.t; };
	unsigned short int a;
	unsigned short int c;
	unsigned short int g;
	unsigned short int t;
};
long int str_to_int(string s)
{
	istringstream is(s);
	long int result = -1;
	is >> result;
	return result;
}
double str_to_double(string s)
{
	istringstream is(s);
	double result = -1;
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
long int min(long int x, long int y){ 
	if(x < y)
		return x;
	else
		return y;
}
long int rand_res(long int total_len, double i)
{
	srand( time(0)*i );
	long int res = 0 + rand() % (total_len - 1);
	return res;
}

char rev(char x)
{
	if(x == 'a')
		return 'c';
	else if(x == 'c')
		return 'g';
	else
		return 'a';

}//rev

double choose(double top, double bottom){
	//given top and bottom it returns Log_10(choose(top, bottom)), where choose is binomial coef
	if(bottom == 0 || bottom == top)
		return 0; //choose = 1, log(choose)=0
	else if(bottom == 1)
		return log10(top);
	double sum = 0.0, current;
	long int i, j;

	long int stop = top - bottom;
	for(i = top; i > stop; i--){
		current = log10(top) - log10(bottom);
		sum += current;
		top -= 1.0;
		bottom -= 1.0;
		
	}//for
	
	 return sum; 
	
}//choose
double p_value(double total_genes, double cluster, double colors, double total_all){
	/*
	Given a set of total_genes' sequences, 
	given a set of chosen genes (cluster) sequences,
	given # occurrences of motif in cluster,
	given # occurrences of motif in all_genes
	
	Find hypergeometric log_10(p_value)
	
	*/
	double stop = total_all;
	if(cluster < total_all)
		stop = cluster;
	double N = total_genes;
	double N_cluster = N - cluster;
	
	long int i;
	double j;//total_all - i
	double sum = 0.0, cluster_choose_i, N_cluster_choose_j , N_choose_total_all, current;
	for(i = colors; i <= stop; i++){
		j = total_all - i;
		cluster_choose_i = choose(cluster, i);
		N_cluster_choose_j = choose(N_cluster, j);
		N_choose_total_all = choose(N, total_all);
		current = (cluster_choose_i + N_cluster_choose_j - N_choose_total_all);// 
		sum += pow(10.0, current);
	}
	
	
	return (-log10(sum));
//return sum;	
}

string rev_complement(string s)
{
	string res("");
	int alen = s.length() - 1; 
	for(int j = alen; j >= 0; j--)
{
		if(s[j] == 'A')
			res = res + 'T';
		else if(s[j] == 'T')
			res = res + 'A';
		else if(s[j] == 't')
			res = res + 'a';
		else if(s[j] == 'a')
			res = res + 't';
		else if(s[j] == 'C')
			res = res + 'G';
		else if(s[j] == 'G')
			res = res + 'C';
		else if(s[j] == 'c')
			res = res + 'g';
		else if(s[j] == 'g')
			res = res + 'c';
		else
			res = res + 'N';
		
	}
	return res;

}//rev_comple


void count_copies(bool is_pe, bool is_singles, bool is_singles2, vector< vector<long int> > &genome,
	vector< vector<long int> > &rev_genome, vector<string> names, vector<string> names_singles,
	vector<string> names_singles2, vector< map<string, Point> > &ref_seq, long int t, long int chrom_num)
{
	string read1, read2, strand;
	string chrom_name ;
	long int id, i, j, asize, u, y, chrom_id, pos1, pos2, mism, mism1, mism2, st1, st2, en1, en2;
	long int orig_start1, orig_start2;
	ifstream in;
		map<string, Point>::iterator end = ref_seq[t].end();
		if(is_pe == true){

cout << "Processing Pair-End results: counting copy-duplicates at each location" << endl;
			for( y = 0; y < names.size(); y++){


				in.open(names[y].c_str(), ios::in);
				if(!in){
					cerr << "ERROR: can't open " << names[y] << endl;
					exit(0);
				}
				else
					cout << "opening " << names[y] << endl;

				long int count_lines = 0;
				in >> id;


				while(!in.eof()){
					count_lines++;
					in >> read1 ;
					in >> read2;
					in >> chrom_name; //starts with 1

					in >> strand;
					in >> pos1;
					in >> pos2;
					in >> mism1;
					in >> mism2;
					in >> orig_start1;
					in >> orig_start2;


					map<string, Point>::iterator afind = ref_seq[t].find(chrom_name);

					if(afind != end){

						chrom_id = (afind->second).x ;
						asize = (afind->second).y;
						if(chrom_id < 0 || chrom_id >= chrom_num)
						{
							cout << "chrom_id is out of boundary at line " << count_lines << endl;
							exit(0);
						}

						long int forw_pos = orig_start1;

						if(orig_start2 < orig_start1){//by leftmost read of a pair
							forw_pos = orig_start2;
						}
						genome[chrom_id][forw_pos]++;
					}//if there is chromosome in a current chunk
					in >> id;
				}//while
				in.close(); in.clear();

			}//for y all files with results
		
		}//if pair end
	

		if(is_singles == true){
			cout << "Processing Single-End results: counting copy duplicates" << endl;
			for( y = 0; y < names_singles.size(); y++){
				in.open(names_singles[y].c_str(), ios::in);
				if(!in){
					cerr << "ERROR: can't open " << names_singles[y] << endl;
					exit(0);
				}
				else
					cout << "opening " << names_singles[y] << endl;
				long int count_lines = 0;
				in >> id;
				while(!in.eof()){
					in >> read1 ;
					in >> chrom_name; //starts with 1
					in >> strand;
					in >> pos1;
					in >> mism;
					in >> orig_start1;

					map<string, Point>::iterator afind = ref_seq[t].find(chrom_name);
					if(afind != end){
						chrom_id = (afind->second).x ;
						asize = (afind->second).y;
						if(chrom_id < 0 || chrom_id >= chrom_num)
						{
							cout << "chrom_id is out of boundary at line " << count_lines << endl;
							exit(0);
						}			
						rev_genome[chrom_id][orig_start1]++;
					}//if there is chrom in this chunk
					in >> id;
				}//while
				in.close(); in.clear();
			}//for y all files with results

		}//if singles

                if(is_singles2 == true){
					cout << "Processing Single-End results, Mates 2: counting copy-duplicates" << endl;
                    for( y = 0; y < names_singles2.size(); y++){                                
						in.open(names_singles2[y].c_str(), ios::in);                                
						if(!in){                                      
							cerr << "ERROR: can't open " << names_singles2[y] << endl;                                        
							exit(0);                               
						}
						else                                        
							cout << "opening " << names_singles2[y] << endl;                               
						long int count_lines = 0;                              
						in >> id;
                                
						while(!in.eof()){                                  
							in >> read1 ;                                  
							in >> chrom_name; //starts with 1                                   
							in >> strand;                                        
							in >> pos1;                                        
							in >> mism;
							in >> orig_start2;
						
							map<string, Point>::iterator afind = ref_seq[t].find(chrom_name);                                        
							if(afind != end){                                              
								chrom_id = (afind->second).x ;                                                
								asize = (afind->second).y;                                                
								if(chrom_id < 0 || chrom_id >= chrom_num)                                                
								{                                               
									cout << "chrom_id is out of boundary at line " << count_lines << endl;
                                    exit(0);
								}
								rev_genome[chrom_id][orig_start2]++;
							}//if there is chrom in this chunk

                            in >> id;
                                
						}//while                                
						in.close(); in.clear();                        
					}//for y all files with results

                }//if singles2
}//count_copies()

void choose_random(vector< vector<long int> > &genome, vector< vector<long int> > &rev_genome,
	vector<long int> &sizes_chroms_cur_map, long int amap_size)
{
	long int i, j;
		for(i = 0; i < amap_size; i++){
			long int cur_size = sizes_chroms_cur_map[i];
			for(j = 0; j < cur_size; j++){
				if(genome[i][j] > 0){
					rev_genome[i][j] = 0;//if there are pairs mapped to this location, we'll take random among them: coverage from pairs is greater than from singles			}
					if(genome[i][j] > 1){
						long int res = rand() % genome[i][j] + 1;//random pair
						genome[i][j] = res;
					}//if copy duplicates
				}//if pairs
				else{
					if(rev_genome[i][j] > 1)
					{
						long int res = rand() % rev_genome[i][j] + 1;//random single
						rev_genome[i][j] = res;
					}//if copy duplicates
				}//only singles
			}//for j
		}//for i
}//choose_random

void print(bool is_pe, bool is_singles, bool is_singles2, vector< vector<long int> > &genome,
	vector< vector<long int> > &rev_genome, vector<string> names, vector<string> names_singles,
	vector<string> names_singles2, vector< map<string, Point> > &ref_seq, long int t, long int chrom_num)
{
	string read1, read2, strand;
	string chrom_name ;
	long int id, i, j, asize, u, y, chrom_id, pos1, pos2, mism, mism1, mism2, st1, st2, en1, en2;
	long int orig_start1, orig_start2;
	ifstream in;
	ofstream out;
		map<string, Point>::iterator end = ref_seq[t].end();
		if(is_pe == true){

cout << "Processing Pair-End results: removing duplicates" << endl;
			for( y = 0; y < names.size(); y++){


				in.open(names[y].c_str(), ios::in);
				if(!in){
					cerr << "ERROR: can't open " << names[y] << endl;
					exit(0);
				}
				else
					cout << "opening " << names[y] << endl;

				string output = names[y] + ".nodupl";
				out.open(output.c_str(), ios::app);

				long int count_lines = 0;
				in >> id;


				while(!in.eof()){
					count_lines++;
					in >> read1 ;
					in >> read2;
					in >> chrom_name; //starts with 1

					in >> strand;
					in >> pos1;
					in >> pos2;
					in >> mism1;
					in >> mism2;
					in >> orig_start1;
					in >> orig_start2;


					map<string, Point>::iterator afind = ref_seq[t].find(chrom_name);

					if(afind != end){

						chrom_id = (afind->second).x ;
						asize = (afind->second).y;
						if(chrom_id < 0 || chrom_id >= chrom_num)
						{
							cout << "chrom_id is out of boundary at line " << count_lines << endl;
							exit(0);
						}

						long int forw_pos = orig_start1;

						if(orig_start2 < orig_start1){//by leftmost read of a pair
							forw_pos = orig_start2;
						}

						genome[chrom_id][forw_pos]--;
						if(genome[chrom_id][forw_pos] == 0){
							out << id << "\t" << read1 << "\t" << read2 << "\t"
							<< chrom_name << "\t" << strand << "\t"	<< pos1 << "\t"
							<< pos2 << "\t"	<< mism1 << "\t" << mism2 << "\t"
							<< orig_start1 << "\t"	<< orig_start2 << endl;						
						}//if
					}//if there is chromosome in a current chunk
					in >> id;
				}//while
				in.close(); in.clear();
				out.close(); out.clear();
			}//for y all files with results
		
		}//if pair end
	

		if(is_singles == true){
			cout << "Processing Single-End results: removing duplicates" << endl;
			for( y = 0; y < names_singles.size(); y++){
				in.open(names_singles[y].c_str(), ios::in);
				if(!in){
					cerr << "ERROR: can't open " << names_singles[y] << endl;
					exit(0);
				}
				else
					cout << "opening " << names_singles[y] << endl;
				string output = names_singles[y] + ".nodupl";
				out.open(output.c_str(), ios::app);

				long int count_lines = 0;
				in >> id;
				while(!in.eof()){
					in >> read1 ;
					in >> chrom_name; //starts with 1
					in >> strand;
					in >> pos1;
					in >> mism;
					in >> orig_start1;

					map<string, Point>::iterator afind = ref_seq[t].find(chrom_name);
					if(afind != end){
						chrom_id = (afind->second).x ;
						asize = (afind->second).y;
						if(chrom_id < 0 || chrom_id >= chrom_num)
						{
							cout << "chrom_id is out of boundary at line " << count_lines << endl;
							exit(0);
						}			
						rev_genome[chrom_id][orig_start1]--;
						if(rev_genome[chrom_id][orig_start1] == 0){
							out << id << "\t" << read1 << "\t" 
							<< chrom_name << "\t" << strand << "\t"	<< pos1 << "\t"
							<< mism << "\t" << orig_start1 << endl;						
						}//if
					}//if there is chrom in this chunk
					in >> id;
				}//while
				in.close(); in.clear();
				out.close(); out.clear();
			}//for y all files with results

		}//if singles

                if(is_singles2 == true){
					cout << "Processing Single-End results, Mates 2: removing duplicates" << endl;
                    for( y = 0; y < names_singles2.size(); y++){                                
						in.open(names_singles2[y].c_str(), ios::in);                                
						if(!in){                                      
							cerr << "ERROR: can't open " << names_singles2[y] << endl;                                        
							exit(0);                               
						}
						else                                        
							cout << "opening " << names_singles2[y] << endl;

						string output = names_singles2[y] + ".nodupl";
						out.open(output.c_str(), ios::app);		

						long int count_lines = 0;                              
						in >> id;
                                
						while(!in.eof()){                                  
							in >> read1 ;                                  
							in >> chrom_name; //starts with 1                                   
							in >> strand;                                        
							in >> pos1;                                        
							in >> mism;
							in >> orig_start2;
						
							map<string, Point>::iterator afind = ref_seq[t].find(chrom_name);                                        
							if(afind != end){                                              
								chrom_id = (afind->second).x ;                                                
								asize = (afind->second).y;                                                
								if(chrom_id < 0 || chrom_id >= chrom_num)                                                
								{                                               
									cout << "chrom_id is out of boundary at line " << count_lines << endl;
                                    exit(0);
								}
								rev_genome[chrom_id][orig_start2]--;
								if(rev_genome[chrom_id][orig_start2] == 0){
									out << id << "\t" << read1 << "\t" 
									<< chrom_name << "\t" << strand << "\t"	<< pos1 << "\t"
									<< mism << "\t" << orig_start2 << endl;						
								}//if
							}//if there is chrom in this chunk

                            in >> id;
                                
						}//while                                
						in.close(); in.clear();   
						out.close(); out.clear();
					}//for y all files with results
		
                }//if singles2
}//count_copies()

void parse_options(int argc, char* argv[], string &pairs, 
				   string &singles, string &singles2, string &genome)
{
	long int i, j;
	for(i = 1; i < argc; i++){
		if(argv[i][0] != '-'){
			cerr << "\nERROR: Invalid char=" << argv[i][0] << "=this" << endl ;		
			return;
		}
		switch(argv[i][1]){
		case 'p': pairs = argv[++i]; break;
		case 's': singles = argv[++i]; break;
		case '1': singles = argv[++i]; break;
		case '2': singles2 = argv[++i]; break;
		case 'r': genome = argv[++i]; break;
		default: cerr << "ERROR: Invalid char=" << argv[i][1] << "=this" << endl ;  exit(0); break;
		}
	}//for i

}//parse_opt()


int main(int argc, char* argv[]){

	srand( time(0) );

	string read1, read2, strand;
	long int id, i, j, asize, u, y, chrom_id, pos1, pos2, mism, mism1, mism2, st1, st2, en1, en2;
	long int orig_start1, orig_start2;
	long int count = 0;
	bool no_dupl = false;
	bool is_pe = false;
	string pairs(""), singles(""), singles2(""), pref, references, out_type("");//out_type default is matrix of ACGT-counts
	parse_options(argc, argv, pairs, singles, singles2, references);


	if(pairs.length() > 0)
		is_pe = true;

	ifstream in;


	
	in.open(references.c_str(), ios::in);
	if(!in){
		cerr << "\nERROR: cannot open " << references << endl;
		exit(0);
	}
	vector<string> gen_names;
	in >> read1;
	while(!in.eof()){
		gen_names.push_back(read1);
		in >> read1;
	}//while
	in.close(); in.clear();

	vector<long int> size_chrom;
	
	asize = gen_names.size();
	string line;
	map<string, long int> chrom_names;
	vector<string> chroms_names;
		string chrom_name ;

	for(i = 0; i < asize; i++){
	
		in.open(gen_names[i].c_str(), ios::in);

		if(!in){
			cerr << "\nERROR: cannot open " << gen_names[i] << endl;
			exit(0);
		}//if
cout << "opening " << gen_names[i] << endl;
		getline(in, line);//first line of fasta file
		istringstream istr(line);
		char ch;
		istr >> ch;
		istr >> chrom_name;
		chrom_names.insert(make_pair(chrom_name, i));
		chroms_names.push_back(chrom_name);

		getline(in, line);
		long int cur_size = 0;
		while(!in.eof()){
			cur_size += line.length();
			getline(in, line);
		}//while
		in.close(); in.clear();

		size_chrom.push_back(cur_size);
cout << "size of " << chrom_name << " is " << cur_size << endl;
	}//for

	long int num_chroms = size_chrom.size();

	//to save space, we process references sequentially
	long int thresh_size = 15000000;
	vector< map<string, Point> > ref_seq;
	map<string, Point> first_map;
	ref_seq.push_back(first_map);
	long int last_map_ind = 0;
	long int cur_index = 0;
	long int cur_size = 0;

cout << "current map: " ;
	for(i = 0; i < num_chroms; i++){
		cur_size += size_chrom[i];
		ref_seq[last_map_ind].insert(make_pair(chroms_names[i], Point(cur_index, size_chrom[i], i)));
		cur_index++;
		cout << chroms_names[i] << ", " ;
		if(cur_size >= thresh_size && i < num_chroms-1){
			map<string, Point> amap;
			ref_seq.push_back(amap);
			last_map_ind++;
			cur_size = 0;
			cur_index = 0;
cout << endl;
cout << "current map: " ;
		}//if
	
	}//for i
cout << endl;


	vector<string> names;
	vector<string> names_singles2;
	vector<string> names_singles;

	if(is_pe == true){
		string aname;
		in.open(pairs.c_str(), ios::in);
		if(!in){
			cerr << "ERROR: can't open " << pairs << endl;
			exit(0);
		}
		in >> aname;
		while(!in.eof()){
			names.push_back(aname);
			in >> aname;
		}//while
		in.close(); in.clear();
cout << "the number of files with pair-end results is " << names.size() << endl;
	}//if

	long int chrom_num = size_chrom.size();
	string aname;
	bool is_singles = false;
	bool is_singles2 = false;
	if(singles2.length() > 0){
		is_singles2 = true;
	}
	if(singles.length() > 0){
		is_singles = true;
	}
	if(is_singles == true){
		in.open(singles.c_str(), ios::in);
		if(!in){
			cerr << "ERROR: can't open " << singles << endl;
			exit(0);
		}
		in >> aname;
		while(!in.eof()){
			names_singles.push_back(aname);
			in >> aname;
		}//while
		in.close(); in.clear();

cout << "the number of files with single-end results is " << names_singles.size() << endl;

	}//if singles


        if(is_singles2 == true){
                in.open(singles2.c_str(), ios::in);
                if(!in){
                        cerr << "ERROR: can't open " << singles2 << endl;
                        exit(0);
                }
                in >> aname;
                while(!in.eof()){
                        names_singles2.push_back(aname);
                        in >> aname;
                }//while
                in.close(); in.clear();

cout << "the number of files with single-end results is " << names_singles2.size() << endl;

        }//if singles2




	long int t;
	cur_size = ref_seq.size();
	long int amap_size ;
	for(t = 0; t < cur_size; t++){

		amap_size = ref_seq[t].size();

		vector<long int> sizes_chroms_cur_map(amap_size);
		map<string, Point>::iterator start = ref_seq[t].begin();
		map<string, Point>::iterator end = ref_seq[t].end();
		for(; start != end; start++){
			long int ind = (start->second).x;
			long int asize_cur = (start->second).y;
			sizes_chroms_cur_map[ind] = asize_cur;
		}//for

		vector< vector<long int> > genome;//for paired-end reads to count the number of pairs mapped to each genome position
		vector< vector<long int> > rev_genome;//for single-end reads to count the number of reads mapped to each genome pos
		for(i = 0; i < amap_size; i++){
			vector<long int> dum(sizes_chroms_cur_map[i]);
			for(j = 0; j < sizes_chroms_cur_map[i]; j++)
				dum[j] = 0;
			genome.push_back(dum);
			rev_genome.push_back(dum);
		}//for i

		//count copy-duplicates for each location
		count_copies(is_pe, is_singles, is_singles2, genome, rev_genome, names, names_singles,
			names_singles2, ref_seq, t, chrom_num);

		//if count > 1 at some location, choose random representative
		choose_random(genome, rev_genome, sizes_chroms_cur_map, amap_size);

		//print out the mapped results without copy-duplicates
		print(is_pe, is_singles, is_singles2, genome, rev_genome, names, names_singles,
			names_singles2, ref_seq, t, chrom_num);

	}//for ref_seq size

	return 0;

}

