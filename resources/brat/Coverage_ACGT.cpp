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

void fill_line(vector<unsigned int> &mark, string &prev, string &cur, string &next, long int alen_prev, 
	long int alen, long int alen_next, long int &cur_base)
{
	long int j = 0;
			for(j = 0; j < alen; j++){					
				if(cur[j] == 'C' || cur[j] == 'c'){
					//determine CHH, CHG or CG
					if(j < alen - 2){//two following bases exist
						if(cur[j+1] == 'G' || cur[j+1] == 'g')
							mark[cur_base]=1;
						else{
							if(cur[j+2] == 'G' || cur[j + 2] == 'g')
								mark[cur_base] = 2;
							else
								mark[cur_base] = 3;
						}//else first consecutive is not G
					}//if fits into a line
					else if(j == alen - 2){
						if(cur[j+1] == 'G' || cur[j+1] == 'g')
							mark[cur_base]=1;
						else{
							if(!next.empty()){
								if(next[0] == 'G' || next[0] == 'g')
									mark[cur_base] = 2;
								else
									mark[cur_base] = 3;
							}
							else
								mark[cur_base] = 2;//CH? preference is to CHG
						}//else
					}//else if 
					else if(j == alen - 1){
						if(!next.empty()){
							if(next[0] == 'G' || next[0] == 'g')
								mark[cur_base] = 1;
							else{
								if(alen_next > 1){
									if(next[1] == 'G' || next[1] == 'g')
										mark[cur_base] = 2;
									else
										mark[cur_base] = 3;
								}
								else
									mark[cur_base] = 2;//CH? preference is to CHG
							}//else
						}//next
						else
							mark[cur_base] = 1;//we don't know, but preference is to CG
					}//alen - 1
				}//if C							
				else if(cur[j] == 'G' || cur[j] == 'g'){
					if(j > 1){
						if(cur[j-1] == 'C' || cur[j-1] == 'c')
							mark[cur_base] = 4;
						else{
							if(cur[j-2] == 'C' || cur[j-2] == 'c')
								mark[cur_base] = 5;
							else
								mark[cur_base] = 6;
						}//not C

					}//if fits into line
					else if(j == 1){
						if(cur[0] == 'C' || cur[0] == 'c')
							mark[cur_base] = 4;
						else{
							if(alen_prev > 0){
								if(prev[alen_prev - 1] == 'C' || prev[alen_prev - 1] == 'c')
									mark[cur_base] = 5;
								else
									mark[cur_base] = 6;
							}//prev not empty
							else
								mark[cur_base] = 5;
						}//else CH
					}//j = 1
					else if(j == 0){
						if(alen_prev > 0){
							if(prev[alen_prev - 1] == 'C' || prev[alen_prev - 1] == 'c')
								mark[cur_base] = 4;
							else{//CH
								if(alen_prev > 1){
									if(prev[alen_prev - 2] == 'C' || prev[alen_prev - 2] == 'c')
										mark[cur_base] = 5;
									else
										mark[cur_base] = 6;
								}
								else
									mark[cur_base] = 5;
							}//else CH
						}//prev
						else
							mark[cur_base] = 4;
					}//if j = 0
				}//if base is G						
												
				cur_base++;						
			}//for j				
}
void fill_content(vector<string> &gen_names, vector<unsigned int> &mark, long int gen_ind)
{
	ifstream in;
	string line, prev, next, cur;
	long int j;
	in.open(gen_names[gen_ind].c_str(), ios::in);					
	getline(in, line);//first FASTA line					
	getline(in, prev);//at least one line in the reference will exist
	getline(in, cur);//current line
	//NEED to process previous	
	long int alen_prev = 0, alen = 0;
	alen_prev = prev.length();
	alen = cur.length();
	long int cur_base = 0;		
	line="";
	if(!prev.empty()){
		//mark previous line's characters
		if(in.eof()){
			cur ="";
			fill_line(mark, line, prev, cur, 0, alen_prev, 0, cur_base);
			return;
		}//if second line is eof
		else
			fill_line(mark, line, prev, cur, 0, alen_prev, alen, cur_base);//current exists
	}//process prev
	else
		return;//no lines in reference

	getline(in, line);//next
	if(in.eof()){//next is eof, but current is not
		next="";
		fill_line(mark, prev, cur, next, alen_prev, alen, 0, cur_base);//process current
		return;
	}//if

		long int alen_next;
		while(!in.eof()){//if reference is more than three lines		
			next = line;
			alen_next = next.length();		
			fill_line(mark, prev, cur, next, alen_prev, alen, alen_next, cur_base);

			prev = cur;
			cur = next;
			alen_prev = alen;
			alen = alen_next;

			getline(in, line);	
		}//while

	
	in.close(); in.clear();	
	next = "";
	//process the last line
	fill_line(mark, prev, cur, next, alen_prev, alen, 0, cur_base);	
}//fill_content()

void print(vector< map<string, Point> > &ref_seq, vector< vector<Cov> > &genome,
	vector< vector<Cov> > &rev_genome, vector<string> &gen_names, string pref, string out_type,
	ofstream &outf, ofstream &outr, long int t){

		string line;
		long int i, j;

		ofstream out;
		ifstream in;

		map<string, Point>::iterator start = ref_seq[t].begin();
		map<string, Point>::iterator end = ref_seq[t].end();



			if(out_type == ""){

			start = ref_seq[t].begin();

			for(; start != end; start++){
				string output =  pref + "_forw_" + start->first + ".txt";
				long int ind = (start->second).x;
				long int asize_of_chrom = (start->second).y;

				out.open(output.c_str(), ios::out);
cout << "Outputting " << output << endl;				
				for(j = 0; j < asize_of_chrom; j++){
					out << genome[ind][j].a << "\t" << genome[ind][j].c << "\t" << genome[ind][j].g << "\t" << genome[ind][j].t << endl;
				}//for j
			out.close(); out.clear();
			}//for i

			start = ref_seq[t].begin();

			for(; start != end; start++){
				string output =  pref + "_rev_" + start->first + ".txt";
				long int ind = (start->second).x;
				long int asize_of_chrom = (start->second).y;

				out.open(output.c_str(), ios::out);
				
				for(j = 0; j < asize_of_chrom; j++){
					out << rev_genome[ind][j].a << "\t" << rev_genome[ind][j].c << "\t" << rev_genome[ind][j].g << "\t" << rev_genome[ind][j].t << endl;
				}//for j
			out.close(); out.clear();
			}//for i

			}//if out_type is default matrix of ACGT-counts
			else{//NEW format: chr, start, end, coverage=C+T (forw) or G+A(rev), meth level C/(C+T) or G/(G+A)

				start = ref_seq[t].begin();

				for(; start != end; start++){
					string achrom =  start->first;

					long int ind = (start->second).x;
					long int asize_of_chrom = (start->second).y;
					long int gen_ind = (start->second).z;

					vector<unsigned int> mark(asize_of_chrom);

					for(j = 0; j < asize_of_chrom; j++){
						mark[j] = 0;
					}
					fill_content(gen_names, mark, gen_ind);

					for(j = 0; j < asize_of_chrom; j++){
						
						if(mark[j] >= 1 && mark[j] < 4){
							double total = 0.0 + genome[ind][j].c + genome[ind][j].t;
							double meth_level = 0;
							if(total > 0)
								meth_level = (genome[ind][j].c + 0.0)/total;
							string acontent("CG:");
							if(mark[j] == 2)
								acontent = "CHG:";
							else if(mark[j] == 3)
								acontent = "CHH:";

							outf << achrom << "\t" << j << "\t" << j << "\t" << acontent << total << "\t" << meth_level << "\t" 
								<< "+" <<  endl;

						}//if C
						else if(mark[j] >= 4){
							double total = 0.0 + rev_genome[ind][j].g + rev_genome[ind][j].a;
							double meth_level = 0;
							string acontent("CG:");
							if(mark[j] == 5)
								acontent = "CHG:";
							else if(mark[j] == 6)
								acontent = "CHH:";

							if(total > 0)
								meth_level = (rev_genome[ind][j].g + 0.0)/total;
							outr << achrom << "\t" << j << "\t" << j << "\t" << acontent << total << "\t" << meth_level << "\t" 
								<< "-" <<  endl;
						}//else if G
					}//for j

				}//for i


			}//else if out_type is BED format

}//print()

void increment_forw(vector< vector<Cov> > &genome, long int chrom_id, long int asize, 
	long int forw_pos, long int stop, long int read_pos, string forw_read)
{
	long int i, j = read_pos;

	for(i = forw_pos; i < stop; i++){

		if(i >= 0 && i < asize){
				if(forw_read[j] == 'A' || forw_read[j] == 'a')
					genome[chrom_id][i].a++;
				else if(forw_read[j] == 'C' || forw_read[j] == 'c')
					genome[chrom_id][i].c++;
				else if(forw_read[j] == 'G' || forw_read[j] == 'g')
					genome[chrom_id][i].g++;
				else if(forw_read[j] == 'T' || forw_read[j] == 't')
					genome[chrom_id][i].t++;
		}//if valid position
		j++;
	}//for i
}//increment_forw()

void increment_rev(vector< vector<Cov> > &genome, long int chrom_id, long int asize, 
	long int rev_pos, long int stop, long int read_pos, string rev_read)
{
	long int i, j = read_pos;		
	
	for(i = rev_pos; i < stop; i++){

								
		if(i >= 0 && i < asize){									
			if(rev_read[j] == 'A' || rev_read[j] == 'a')							
				genome[chrom_id][i].t++;									
			else if(rev_read[j] == 'C' || rev_read[j] == 'c')									
				genome[chrom_id][i].g++;									
			else if(rev_read[j] == 'G' || rev_read[j] == 'g')									
				genome[chrom_id][i].c++;									
			else if(rev_read[j] == 'T' || rev_read[j] == 't')									
				genome[chrom_id][i].a++;							
		}
		j--;
	}//for i
}//increment_rev from reverse-complement of read
void parse_options(int argc, char* argv[], string &pairs, 
				   string &singles, string &singles2, string &pref, string &genome, string &out_type, bool &no_dupl)
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
		case 'P': pref = argv[++i]; break;
		case 'r': genome = argv[++i]; break;
		case 'B': out_type = "bed"; break;
		case 'N': no_dupl = true; break;
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
	parse_options(argc, argv, pairs, singles, singles2, pref, references, out_type, no_dupl);

	ofstream outf, outr;
	string out_forw, out_rev;

	if(out_type == "bed"){
		out_forw =  pref + "_forw" + ".txt";
		out_rev = pref + "_rev" + ".txt";
		outf.open(out_forw.c_str(), ios::out);
		outr.open(out_rev.c_str(), ios::out);
	}
	if(pairs.length() > 0)
		is_pe = true;

	ifstream in;


	if(no_dupl == false){

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

		vector< vector<Cov> > genome;
		vector< vector<Cov> > rev_genome;
		for(i = 0; i < amap_size; i++){
			vector<Cov> dum(sizes_chroms_cur_map[i]);
			genome.push_back(dum);
			rev_genome.push_back(dum);
		}//for i

		if(is_pe == true){

cout << "Processing Pair-End results" << endl;
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
						long int len_mer1 = read1.length();
						long int len_mer2 = read2.length();
						long int forw_len = len_mer1;
						long int rev_len = len_mer2;

						string forw_read = read1;
						string rev_read = read2;
						long int forw_pos = pos1;
						long int rev_pos = pos2;
						if(strand == "-"){
							forw_read = read2;
							rev_read = read1;
							forw_pos = pos2;
							rev_pos = pos1;
							forw_len = len_mer2;
							rev_len = len_mer1;
						}
						j = 0;
						long int stop = forw_pos + forw_len;
						if(strand == "+"){
							increment_forw(genome, chrom_id, asize, forw_pos, stop, j, forw_read);
						}//if first mate is mapped to positive strand 
						else{
							increment_forw(rev_genome, chrom_id, asize, forw_pos, stop, j, forw_read);
						}//else rev-complement of first map is mapped to positive strand
						//process reverse read
						//check if there is an overlap between two mates
						bool overlap = false;



						//Case when pair comes from two-strands sequencing: first mates are from original genomic strands
						if(forw_pos < rev_pos){//this condition only works if min-insert is a positive integer

						if(rev_pos < forw_pos + forw_len)
							overlap = true;
						if(overlap == false){//no overlap

							j=rev_len -1;
							if(strand == "-"){
								long int stop = rev_pos + rev_len;
								increment_rev(rev_genome, chrom_id, asize, rev_pos, stop, j, rev_read);
							}//if strand -
							else{
								long int stop = rev_pos + rev_len;
								increment_rev(genome, chrom_id, asize, rev_pos, stop, j, rev_read);
							}//else strand +
						}//no overlap
						else{//overlap
							long int pos_overlapped = forw_pos + forw_len - rev_pos; 
							long int rev_stop = rev_pos + rev_len;
                            rev_pos = forw_pos + forw_len;							
                            j=rev_len -1 - pos_overlapped;
							if(strand == "-"){
								increment_rev(rev_genome, chrom_id, asize, rev_pos, rev_stop, j, rev_read);
							}//if strand -
							else{
								increment_rev(genome, chrom_id, asize, rev_pos, rev_stop, j, rev_read);
							}//if strand +
						}//overlap
						}//if first mates are from original genomic strands
						else{//first mates are from rev-complements PCR-products of original genomic strans
							if(forw_pos < rev_pos + rev_len)
								overlap = true;
							if(overlap == false){//no overlap

								j=rev_len -1;
								if(strand == "-"){
									long int stop = rev_pos + rev_len;
									increment_rev(rev_genome, chrom_id, asize, rev_pos, stop, j, rev_read);
								}//if strand -
								else{
									long int stop = rev_pos + rev_len;
									increment_rev(genome, chrom_id, asize, rev_pos, stop, j, rev_read);
								}//else strand +
							}//no overlap
							else{//overlap
								j = rev_len - 1;
								long int stop = forw_pos;
								if(strand == "-"){
									increment_rev(rev_genome, chrom_id, asize, rev_pos, stop, j, rev_read);
								}//neg
								else{
									increment_rev(genome, chrom_id, asize, rev_pos, stop, j, rev_read);								
								}//positive strand
							}//else overlap
						}//else first mates are from rev-complement PCR-products of original genomic strands
					}//if there is chromosome in a current chunk
					in >> id;
				}//while
				in.close(); in.clear();

			}//for y all files with results
			cout << "Finished calculating ACGT-count from paired-end reads." << endl;
		}//if pair end
	

		if(is_singles == true){
			cout << "Processing Single-End results" << endl;
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
						long int len_mer = read1.length();
						j = 0;
						long int stop = pos1 + len_mer;
						if(strand == "+"){
							increment_forw(genome, chrom_id, asize, pos1, stop, j, read1);
						}else{
							j=len_mer -1;
							increment_rev(rev_genome, chrom_id, asize, pos1, stop, j, read1);
						}//else reverse
					}//if there is chrom in this chunk
					in >> id;
				}//while
				in.close(); in.clear();
			}//for y all files with results
			cout << "Finished calculating ACGT-count from single reads." << endl;
		}//if singles

                if(is_singles2 == true){
					cout << "Processing Single-End results, Mates 2" << endl;
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
                                long int len_mer = read1.length();                    
								j = 0;                                                
								long int stop = pos1 + len_mer;	
								if(strand == "+"){																									
									increment_forw(rev_genome, chrom_id, asize, pos1, stop, j, read1);                                                
								}else{                                    
									j=len_mer -1;																	
									increment_rev(genome, chrom_id, asize, pos1, stop, j, read1);
								}//else reverse
							}//if there is chrom in this chunk

                            in >> id;
                                
						}//while                                
						in.close(); in.clear();                        
					}//for y all files with results
					cout << "Finished calculating ACGT-count from single reads, Mates 2." << endl;
                }//if singles2

		if(is_pe == true || is_singles == true || is_singles2 == true){
			print(ref_seq, genome, rev_genome, gen_names, pref, out_type, outf, outr, t);				
		}//if any input
	}//for ref_seq size

	outf.close(); 
	outr.close();

	}//if duplicates are allowed
	else{//NO_DUPL--------------------------------------------------------------------------------------------------------------

	}//no duplicates NO DUPL
	return 0;

}
