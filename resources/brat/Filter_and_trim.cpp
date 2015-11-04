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

class Triple{
public:
	Triple(long int xx, long int yy, long int zz, long int ww) : x(xx), y(yy), z(zz), w(ww) {};
	Triple() : x(0), y(0), z(0), w(0) {};
	void operator =(Triple a){x = a.x; y = a.y; z = a.z; w = a.w;};
	long int x;//id of a read
	long int y;//count number of positions to which read is mapped
	long int z;//count of copies of this read
	long int w;//to hold genome position where read maps to
};


class Point
{

public:
	//constructors of this class
	Point(long int xx, long int yy) : x(xx), y(yy) {};
	Point() : x(0), y(0) {};
	
	
	void operator =(Point a) { x = a.x; y = a.y;};
	//x is index, y is tails
	bool operator <(const Point& a) const
	{ 
		if(y < a.y) return true; 
		else return false;
	}


	long int x;
	long int y;

};

string rev_complement(string s)
{
	string res("");
	for(int j = s.length() - 1; j >= 0; j--)
{
		if(s[j] == 't')
			res = res + 'a';
		else if(s[j] == 'a')
			res = res + 't';
		else if(s[j] == 'c')
			res = res + 'g';
		else if(s[j] == 'g')
			res = res + 'c';
		else if(s[j] == 'A')
			res = res + 'T';
		else if(s[j] == 'T')
			res = res + 'A';
		else if(s[j] == 'C')
			res = res + 'G';
		else if(s[j] == 'G')
			res = res + 'C';
		else
			res = res + 'N';
		
	}
	return res;

}//rev_comple

int hash(char c){
	if(c == 'a')
		return 0;
	else if(c == 'c')
		return 1;
	else if(c == 'g')
		return 2;
	else
		return 3;
}
char char_hash(int x)
{
	if(x == 0)
		return 'a';
	else if(x == 1)
		return 'c';
	else if(x == 2)
		return 'g';
	else
		return 't';

}
long int str_to_int(string s)
{
	istringstream is(s);
	long int result = -1;
	is >> result;
	return result;
}

string int_to_str(int x)
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

int max(int x, int y)
{	
	if(x > y)
		return x;
	else
		return y;

}
int min(int x, int y){ 
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
void parse_options(int argc, char* argv[], string &reads_file1, 
				   string &reads_file2, string &pref, int &thresh, 
				   int &mism_thresh, int &lower_bound )
{
	int res = 0;

	long int i;
	for(i = 1; i < argc; i++){
		if(argv[i][0] != '-'){
			cout << "\nInvalid char=" << argv[i][0] << "=this" << endl ;		
			return;
		}
		switch(argv[i][1]){
		case 's': reads_file1 = argv[++i]; break;
		case '1': reads_file1 = argv[++i]; break;
		case 'm': mism_thresh = atoi(argv[++i]); break;
		case '2': reads_file2 = argv[++i]; break;
		case 'P': pref = argv[++i]; break;
		case 'q': thresh = atoi(argv[++i]); break;
		case 'L': lower_bound = atoi(argv[++i]); break;
		default: cout << "Invalid char=" << argv[i][1] << "=this" << endl ;  exit(0); break;
		}
	}//for i

}//parse_opt()

int main(int argc, char* argv[]){

	string reads_file1(""), reads_file2(""), is_pe, pref;
	int score_thresh = 0;

	int mism_thresh = 0;

	int lower_bound = 33;

	parse_options(argc, argv, reads_file1, reads_file2, pref, score_thresh, mism_thresh, lower_bound);

	if(reads_file2.length() > 0)
		is_pe = "1";

	int min_len = 23;
	
	long int  j;
	string aread, aread2, aname, read1, read2;

	ofstream out;
	ifstream in, in2;//to read input seq and prb
	

	int len ;
	int mism = 0;
	string qual, qual2;
	//open files one at a time

	int alph = 4;
	//to count the total number of sequenced reads
	long int total_reads = 0;

	//to hold the sum of qual score at each position

		//construct file name
		

		in.open(reads_file1.c_str(), ios::in);
		if(!in){
			cout << "\nERROR: Cannot open " <<reads_file1<< endl;
			exit(0);
		}
		else
			cout << "Opening " << reads_file1 << endl;
		if(is_pe == "1"){
			in2.open(reads_file2.c_str(), ios::in);
			if(!in2){
				cout << "\nERROR: Cannot open " << reads_file2 << endl;
				exit(0);
			}
			else 
				cout << "Opening " << reads_file2 << endl;
				//check for internal Ns
		}//if pair-end

	//to output mates of dif length, to output singles
	ofstream out1, out2, out_singles1, out_singles2;
	string output1 = pref + "_reads1.txt";
	string output2 = pref + "_reads2.txt";
	string output5 = pref + "_mates1.txt" ;
	string output6 = pref + "_mates2.txt";

	ofstream outf1, outf2, outf3, outf4;

	string output_fastq1 = pref + "_pair1.fastq";
	string output_fastq2 = pref + "_pair2.fastq";
	string output_s1 = pref + "_mates1.fastq";
	string output_s2 = pref + "_mates2.fastq";
	outf3.open(output_s1.c_str(), ios::out);
	outf1.open(output_fastq1.c_str(), ios::out);

	out1.open(output1.c_str(), ios::out);
	out_singles1.open(output5.c_str(), ios::out);
	if(is_pe == "1"){
		out_singles2.open(output6.c_str(), ios::out);
		out2.open(output2.c_str(), ios::out);
		outf2.open(output_fastq2.c_str(), ios::out);

		outf4.open(output_s2.c_str(), ios::out);

	}//if pair-end

	string strand, plus, id1, id2, plus2;
	string line;

	if(is_pe == "1"){

		getline(in, id1);
		while(!in.eof()){
			getline(in, read1);
			getline(in, plus);
			getline(in, qual);

			getline(in2, id2);
			getline(in2, read2);
			getline(in2, plus2);
			getline(in2, qual2);

			int st1 =0, st2 = 0;
			len = read1.length();
			int length2 = read2.length();
			int en1 = len-1;
			int en2 = length2-1;
			for(j = 0; j < len; j++){
				int aq1 = (int)qual[j] - lower_bound;
				if(aq1 < score_thresh)
					st1++;
				else
					break;
			}//for j
			for(j = 0; j < length2; j++){
				int aq2 = (int)qual2[j] - lower_bound;
				if(aq2 < score_thresh)
					st2++;
				else
					break;
			}
			//find end points
			for(j = len -1 ; j >= st1; j--){
				int eq = (int)qual[j] - lower_bound;
				if(eq < score_thresh)
					en1--;
				else
					break;
			}
			for(j = length2 -1 ; j >= st2; j--){
				int eq = (int)qual2[j] - lower_bound;
				if(eq < score_thresh)
					en2--;
				else
					break;
			}

			//Trim the ends with Ns
			for(j = st1; j <= en1; j++){
				if(read1[j] == 'N' || read1[j] == 'n')
					st1++;
				else
					break;

			}//for j
			for(j = en1; j >= st1; j--)
			{
				if(read1[j] == 'N' || read1[j] == 'n')
					en1--;
				else
					break;
			}//for j
			for(j = st2; j <= en2; j++){
				if(read2[j] == 'N' || read2[j] == 'n')
					st2++;
				else
					break;

			}//for j
			for(j = en2; j >= st2; j--)
			{
				if(read2[j] == 'N' || read2[j] == 'n')
					en2--;
				else
					break;
			}//for j


			int len1 = en1 - st1 + 1;
			int len2 = en2 - st2 + 1;

			bool found_n1 = false, found_n2 = false;
			//check internal Ns
			long int countN1 = 0, countN2 = 0;
			if(len1 > 0){
				for(j = st1; j <= en1; j++){
					if(read1[j] == 'N' || read1[j] == 'n')
						countN1++;
					
				}
			}
			if(len2 > 0){
				for(j = st2; j <= en2; j++){
					if(read2[j] == 'N' || read2[j] == 'n')
						countN2++;
					
				}
			}
			if(countN1 > mism_thresh)
				found_n1 = true;
			if(countN2 > mism_thresh)
				found_n2 = true;
 
		if(found_n1 == false && found_n2 == false){

			if(len2 > min_len && len1 > min_len){
				read1 = read1.substr(st1, len1);
				read2 = read2.substr(st2, len2);

				qual = qual.substr(st1, len1);
				qual2 = qual2.substr(st2, len2);

				en1= len - en1 - 1;
				en2 = length2 - en2 -1;
					out1 << read1 << "\t" << st1 << "\t" << en1 << endl;
					out2 << read2 << "\t" << st2 << "\t" << en2 << endl;

					outf1 << id1 << endl;
					outf1 << read1 << endl;
					outf1 << plus << endl;
					outf1 << qual << endl;
					
					outf2 << id2 << endl;
					outf2 << read2 << endl;
					outf2 << plus2 << endl;
					outf2 << qual2 << endl;
					

			}//if dif len
			else if(len1 > min_len){
				read1 = read1.substr(st1, len1);
				en1 = len - en1 - 1 ;
				out_singles1 << read1 << "\t" << st1 << "\t" << en1 << endl;
				qual = qual.substr(st1, len1);

					outf3 << id1 << endl;
					outf3 << read1 << endl;
					outf3 << plus << endl;
					outf3 << qual << endl;


			}//if first read has legitimate read len
			else if(len2 > min_len){
				read2 = read2.substr(st2, len2);
				en2 = length2 - en2 -1;

				out_singles2 << read2 << "\t" << st2 << "\t" << en2 << endl;
				qual2 = qual2.substr(st2, len2);


					outf4 << id2 << endl;
					outf4 << read2 << endl;
					outf4 << plus2 << endl;
					outf4 << qual2 << endl;


			}//second read ok
			}//if no Ns in both reads
			else if(found_n1 == false && len1 > min_len){
				read1 = read1.substr(st1, len1);
				en1 = len - en1 - 1;
				out_singles1 << read1 << "\t" << st1 << "\t" << en1 << endl;
				qual = qual.substr(st1, len1);

					outf3 << id1 << endl;
					outf3 << read1 << endl;
					outf3 << plus << endl;
					outf3 << qual << endl;


				
			}//first read is ok
			else if(found_n2 == false && len2 > min_len){
				read2 = read2.substr(st2, len2);
				en2 = length2 - en2 - 1;
				out_singles2 << read2 << "\t" << st2 << "\t" << en2 << endl;
				qual2 = qual2.substr(st2, len2);


					outf4 << id2 << endl;
					outf4 << read2 << endl;
					outf4 << plus2 << endl;
					outf4 << qual2 << endl;
	
			}//second read has no Ns
			getline(in, id1);
		}//while
	}//if pair-end
	else{
	
		getline(in, id1);
		while(!in.eof()){
			getline(in, read1);
			getline(in, plus);
			getline(in, qual);

			len = read1.length();

			int st1 =0;
			int en1 = len-1;
			for(j = 0; j < len; j++){
				int aq1 = (int)qual[j] - lower_bound;
				if(aq1 < score_thresh)
					st1++;
				else
					break;
			}//for j
			//find end points
			for(j = len -1 ; j >= st1; j--){
				int eq = (int)qual[j] - lower_bound;
				if(eq < score_thresh)
					en1--;
				else
					break;
			}

			//Trim the ends with Ns
			for(j = st1; j <= en1; j++){
				if(read1[j] == 'N' || read1[j] == 'n')
					st1++;
				else
					break;

			}//for j
			for(j = en1; j >= st1; j--)
			{
				if(read1[j] == 'N' || read1[j] == 'n')
					en1--;
				else
					break;
			}//for j

			int len1 = en1 - st1 + 1;

			bool found_n1 = false;
			//check internal Ns
			long int countN = 0;
			if(len1 > 0){
				for(j = st1; j <= en1; j++){
					if(read1[j] == 'N' || read1[j] == 'n')
						countN++;
					
				}
			}
			if(countN  > mism_thresh)
				found_n1 = true;
			if(found_n1 == false ){

				if(len1 > min_len){
					read1 = read1.substr(st1, len1);
					qual = qual.substr(st1, len1);
					en1 = len - en1 - 1;
					out1 << read1 << "\t" << st1 << "\t" << en1 << endl;

					outf1 << id1 << endl;
					outf1 << read1 << endl;
					outf1 << plus << endl;
					outf1 << qual << endl;
					

				}//if dif len

			}//if no Ns in both reads
			getline(in, id1);
		}//while	
	}//else
	out1.close(); out2.close(); out_singles1.close(); 
	out_singles2.close();
	outf1.close(); outf2.close();
		in.close(); in.clear();
		in2.close(); in2.clear();
	//calculate the average of the scores
	outf3.close(); outf4.close();

	return 0;
}