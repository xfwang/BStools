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
#include<ctime>
#include<cmath>
#include<vector>
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
using namespace std;

void read_file_names(string from_file, vector<string> &file_names)
{		ifstream in;
		string aname;
		in.open(from_file.c_str(), ios::in);
		if(!in){
			cerr << "ERROR: can't open " << from_file << endl;
			exit(0);
		}
		in >> aname;
		while(!in.eof()){
			file_names.push_back(aname);
			in >> aname;
		}//while
		in.close(); in.clear();
}//read_file_names

void convert_pairs(ofstream &outf, ofstream &outr, vector<string> &names)
{
	ifstream in;
	long int asize = names.size();
	long int i = 0;
	for(i = 0; i < asize; i++){
		in.open(names[i].c_str(), ios::in);
		if(!in){
			cout << "Could not open " << names[i] << endl;
			outf.close(); outr.close();
			exit(0);
		}
		string read_id, read1, read2, ref;
		char strand;
		long int pos1 = 0, pos2 = 0, mism1 = 0, mism2 = 0, orig_pos1 = 0, orig_pos2 = 0;
		in >> read_id;
		while(!in.eof()){
			in >> read1 >> read2 >> ref >> strand >> pos1 >> pos2 >> mism1 >> mism2 >> orig_pos1 >> orig_pos2;
			unsigned long int bit_flag_mate1 = 1 + (1 << 1) + (strand == '-' ? (1 << 4) : 0);
			bit_flag_mate1 += (strand == '+' ? (1 << 5) : 0) + (1 << 6) ;
			unsigned long int bit_flag_mate2 = 1 + (1 << 1) + (strand == '+' ? (1 << 4) : 0);
			bit_flag_mate2 += (strand == '-' ? (1 << 5) : 0) + (1 << 7);
			long int read_len1 = read1.length();
			long int read_len2 = read2.length();
			if(strand == '+'){
				outf << read_id << "\t" << bit_flag_mate1 << "\t" << ref << "\t" << (pos1 + 1) << "\t" ;
				outf << 255 << "\t" << '*' << "\t" << '*' << "\t" << 0 << "\t" ;
				outf << (strand == '+' ? read_len1 : (- read_len1)) << "\t" << read1 << "\t" << '*' << endl;

				outf << read_id << "\t" << bit_flag_mate2 << "\t" << ref << "\t" << (pos2 + 1) << "\t" ;
				outf << 255 << "\t" << '*' << "\t" << '*' << "\t" << 0 << "\t" ;
				outf << (strand == '+' ? ( - read_len2) : read_len2)	<< "\t" << read2 << "\t" << '*' << endl;
			}//if +
			else{
				outr << read_id << "\t" << bit_flag_mate1 << "\t" << ref << "\t" << (pos1 + 1) << "\t" ;
				outr << 255 << "\t" << '*' << "\t" << '*' << "\t" << 0 << "\t" << (strand == '+' ? read_len1 : (- read_len1));
				outr << "\t" << read1 << "\t" << '*' << endl;

				outr << read_id << "\t" << bit_flag_mate2 << "\t" << ref << "\t" << (pos2 + 1) << "\t" ;
				outr << 255 << "\t" << '*' << "\t" << '*' << "\t" << 0 << "\t" << (strand == '+' ? ( - read_len2) : read_len2);
				outr << "\t" << read2 << "\t" << '*' << endl;			
			}//else -
			in >> read_id;
		}//while
		in.close(); in.clear();
	}//for i

}//convert_pairs

void convert_singles1(ofstream &outf, ofstream &outr, vector<string> &names)
{
	ifstream in;
	long int asize = names.size();
	long int i = 0;
	for(i = 0; i < asize; i++){
		in.open(names[i].c_str(), ios::in);
		if(!in){
			cout << "Could not open " << names[i] << endl;
			outf.close(); outr.close();
			exit(0);
		}
		string read_id, read1, ref;
		char strand;
		long int pos1 = 0, mism1 = 0, orig_pos1 = 0;
		in >> read_id;
		while(!in.eof()){
			in >> read1 >> ref >> strand >> pos1 >> mism1 >> orig_pos1;
			unsigned long int bit_flag_mate1 = 1 + (1 << 3) + (strand == '-' ? (1 << 4) : 0);
			bit_flag_mate1 +=  (1 << 6) ;
			if(strand == '+'){
				outf << read_id << "\t" << bit_flag_mate1 << "\t" << ref << "\t" << (pos1 + 1) << "\t" ;
				outf << 255 << "\t" << '*' << "\t" << '*' << "\t" << 0 << "\t" << 0;
				outf << "\t" << read1 << "\t" << '*' << endl;
			}//if +
			else{
				outr << read_id << "\t" << bit_flag_mate1 << "\t" << ref << "\t" << (pos1 + 1) << "\t" ;
				outr << 255 << "\t" << '*' << "\t" << '*' << "\t" << 0 << "\t" << 0;
				outr << "\t" << read1 << "\t" << '*' << endl;	
			}//else -
			in >> read_id;
		}//while
		in.close(); in.clear();
	}//for i

}//convert_singles1

void convert_singles2(ofstream &outf, ofstream &outr, vector<string> &names)
{
	ifstream in;
	long int asize = names.size();
	long int i = 0;
	for(i = 0; i < asize; i++){
		in.open(names[i].c_str(), ios::in);
		if(!in){
			cout << "Could not open " << names[i] << endl;
			outf.close(); outr.close();
			exit(0);
		}
		string read_id, read1, ref;
		char strand;
		long int pos1 = 0, mism1 = 0, orig_pos1 = 0;
		in >> read_id;
		while(!in.eof()){
			in >> read1 >> ref >> strand >> pos1 >> mism1 >> orig_pos1;
			unsigned long int bit_flag_mate1 = 1 + (1 << 3) + (strand == '-' ? (1 << 4) : 0);
			bit_flag_mate1 +=  (1 << 7) ;
			if(strand == '+'){
				outr << read_id << "\t" << bit_flag_mate1 << "\t" << ref << "\t" << (pos1 + 1) << "\t" ;
				outr << 255 << "\t" << '*' << "\t" << '*' << "\t" << 0 << "\t" << 0;
				outr << "\t" << read1 << "\t" << '*' << endl;
			}//if +
			else{
				outf << read_id << "\t" << bit_flag_mate1 << "\t" << ref << "\t" << (pos1 + 1) << "\t" ;
				outf << 255 << "\t" << '*' << "\t" << '*' << "\t" << 0 << "\t" << 0;
				outf << "\t" << read1 << "\t" << '*' << endl;	
			}//else -
			in >> read_id;
		}//while
		in.close(); in.clear();
	}//for i

}//convert_singles2

void convert_singles(ofstream &outf, ofstream &outr, vector<string> &names)
{
	ifstream in;
	long int asize = names.size();
	long int i = 0;
	for(i = 0; i < asize; i++){
		in.open(names[i].c_str(), ios::in);
		if(!in){
			cout << "Could not open " << names[i] << endl;
			outf.close(); outr.close();
			exit(0);
		}
		string read_id, read1, ref;
		char strand;
		long int pos1 = 0, mism1 = 0, orig_pos1 = 0;
		in >> read_id;
		while(!in.eof()){
			in >> read1 >> ref >> strand >> pos1 >> mism1 >> orig_pos1;
			unsigned long int bit_flag_mate1 = (strand == '-' ? (1 << 4) : 0);
		
			if(strand == '+'){
				outf << read_id << "\t" << bit_flag_mate1 << "\t" << ref << "\t" << (pos1 + 1) << "\t" ;
			outf << 255 << "\t" << '*' << "\t" << '*' << "\t" << 0 << "\t" << 0;
				outf << "\t" << read1 << "\t" << '*' << endl;
			}//if +
			else{
				outr << read_id << "\t" << bit_flag_mate1 << "\t" << ref << "\t" << (pos1 + 1) << "\t" ;
				outr << 255 << "\t" << '*' << "\t" << '*' << "\t" << 0 << "\t" << 0;
				outr << "\t" << read1 << "\t" << '*' << endl;	
			}//else -
			in >> read_id;
		}//while
		in.close(); in.clear();
	}//for i

}//convert_singles1

void parse_options(int argc, char* argv[], string &pairs, 
				   string &singles1, string &singles2, string &singles, string &pref)
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
		case '1': singles1 = argv[++i]; break;
		case '2': singles2 = argv[++i]; break;
		case 'P': pref = argv[++i]; break;
		default: cerr << "ERROR: Invalid char=" << argv[i][1] << "=this" << endl ;  exit(0); break;
		}
	}//for i

}//parse_opt()


int main(int argc, char* argv[]){

	srand( time(0) );


	string pairs(""), singles1(""), singles2(""), singles(""), pref("");//out_type default is matrix of ACGT-counts
	parse_options(argc, argv, pairs, singles1, singles2, singles, pref);

	vector<string> names;
	vector<string> names_singles2;
	vector<string> names_singles1;
	vector<string> names_singles;

	bool is_pe = false;
	if(pairs.length() > 0)
		is_pe = true;
	if(is_pe == true){
		read_file_names(pairs, names);
		cout << "the number of files with pair-end results is " << names.size() << endl;
	}//if

	bool is_singles = false;
	if(singles.length() > 0){
		is_singles = true;
	}
	if(is_singles == true){
		read_file_names(singles, names_singles);
		cout << "the number of files with single-end results is " << names_singles.size() << endl;
	}//if singles

	bool is_singles1 = false;
	if(singles.length() > 0){
		is_singles1 = true;
	}
	if(is_singles1 == true){
		read_file_names(singles1, names_singles1);
		cout << "the number of files with single-end results is " << names_singles1.size() << endl;
	}//if singles

	bool is_singles2 = false;
	if(singles2.length() > 0){
		is_singles2 = true;
	}
    if(is_singles2 == true){
		read_file_names(singles2, names_singles2);
		cout << "the number of files with single-end results is " << names_singles2.size() << endl;
    }//if singles2

	string out_forw = (pref.length() > 0 ? (pref + "_forw.sam") : ("mapped_to_forw.sam")) ;
	string out_rev = (pref.length() > 0 ? (pref + "_rev.sam") : ("mapped_to_rev.sam"));
	ofstream outf, outr;
	outf.open(out_forw.c_str(), ios::app);
	outr.open(out_rev.c_str(), ios::app);
	
	if(is_pe == true)
		convert_pairs(outf, outr, names);
	if(is_singles1 == true)
		convert_singles1(outf, outr, names_singles1);
	if(is_singles2 == true)
		convert_singles2(outf, outr, names_singles2);
	if(is_singles == true)
		convert_singles(outf, outr, names_singles);
	outf.close();
	outr.close();
	return 0;

}

