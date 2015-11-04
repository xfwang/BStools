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

void usage(){

	string mes("\nUSAGE: check_strands -r <references> -s <input-reads> -o <output> [OPTIONS]\n\nOptions:\n");
	string gen_ref_opt("  -r <references-file>       file with file names containing references (one\n                          reference per file)");
	string singles(    "  -s <input-reads-file>      file with single reads (or queries) in FASTQ format");
	string output(     "  -o <output-file>");
	string first(      "  -1 <paired-ends-file1>     mates 1 in FASTQ format");
	string second(     "  -2 <paired-ends-file2>     mates 2 in FASTQ format");

	cout << mes << endl;
	cout << gen_ref_opt << endl << singles << endl
        	<< output << endl
			<< first << endl << second << endl;

	string mes2(      "Limitations:\n");
	string ref_size(  "  a single reference size <= 2^45");
	string chrom_ids( "  number of references <= 2^18");
	string space(     "  space <= (total references sizes)*4.5*8 + 16*(2^24)");
	
	cout << mes2 << endl;
	cout << ref_size << endl << chrom_ids << endl << space << endl;


}//usage()

int parse_options(int argc, char* argv[],string &genome_names, string &reads_file1, string &reads_file2, string &output)
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
		case 'o': output = argv[++i]; break;
		case 'h': usage(); break;
		default: cerr << "ERROR: Invalid char=" << argv[i][1] << "=this" << endl; usage(); exit(0); break;
		}
	}//for i
	return res;
}//parse_opt()

long int hash_char(char c)
{
	if(c == 'A' || c == 'a')
		return 0;
	else if(c == 'C' || c == 'c')
		return 1;
	else if(c == 'G' || c == 'g')
		return 2;
	else if(c == 'T' || c == 't')
		return 3;
	else
		return -1;
}//hash_char



int main(int argc, char* argv[])
{
	/**************************************
	This program calculates distribution of different fingerprints in the genome given
	*******************************************/
		/**************************************
		bitset [11000] is translated into a string as 00011 and into long int as 3
		so, least significant bits in the bitset are the leftmost
  *******************************************/
	

	long int i, j;
	long int asize , u, y;

	if(argc < 5){
		usage();
		return 0;
	}//if
	string genome_names(""), reads_file1(""), reads_file2(""), output;
	int res = parse_options(argc, argv, genome_names,  reads_file1,  reads_file2,  output);
	if(res == -1)
	{
		usage();
		exit(0);
	}

	ifstream in;
	in.open(genome_names.c_str(), ios::in);
	if(!in){
		cout << "ERROR: could not open " << genome_names << endl;
		exit(0);
	}
	vector<string> refs;
	string aname, first, aread, third, aqual;
	in >> aname;
	while(!in.eof()){
		refs.push_back(aname);
		in >> aname;
	}//while

	in.close(); in.clear();

	long int as = 0, cs = 0, gs = 0, ts = 0;
	vector<double> gen_acgt(4);
	vector<double> first_acgt(4);
	vector<double> second_acgt(4);
	for(i = 0; i < 4; i++)
	{
		gen_acgt[i] = 0.0;
		first_acgt[i] = 0.0;
		second_acgt[i] = 0.0;
	}//for i

	asize = refs.size();
	string line;
	for(i = 0; i < asize; i++){
		in.open(refs[i].c_str(), ios::in);
		if(!in){
			cout << "ERROR: could not open " << refs[i] << endl;
			exit(0);
		}
		getline(in, line);//first line
		getline(in, line);
		while(!in.eof()){
			long int len = line.length();
			for(j = 0; j < len; j++){
				long int ind = hash_char(line[j]);
				if(ind != -1 && ind < 4){
					gen_acgt[ind] += 1.0;
				}
			}//for j
			getline(in, line);
		}//while
		in.close(); in.clear();
	}//for i
	ofstream out;
	out.open(output.c_str(), ios::out);
	out << "ACGT-distribution in the references:" << endl;
	double total = gen_acgt[0] + gen_acgt[1] + gen_acgt[2] + gen_acgt[3];
	for(j = 0; j < 4; j++){
		double ratio = gen_acgt[j]/total;
		out << ratio << "\t" ;

	}
	out << endl;

	in.open(reads_file1.c_str(), ios::in);
	if(!in){
		cout << "ERROR: could not open " << reads_file1 << endl;
		exit(0);
	}//
	getline(in, line);
	if(line[0] != '@'){
		cout << "ERROR: the files with reads must be in FASTQ format" << endl;
		exit(0);
	}
	while(!in.eof()){
		getline(in, aread);
		getline(in, third);
		getline(in, aqual);
		long int len = aread.length();
		for(j = 0; j < len; j++){
			long int ind = hash_char(aread[j]);
			if(ind != -1 && ind < 4){
				first_acgt[ind] += 1.0;;
			}
		}//for j
		getline(in, line);
	}//while
	in.close(); in.clear();

	out << "ACGT-distribution in the reads (for paired-end, from file with option -1:" << endl;
	total = first_acgt[0] + first_acgt[1] + first_acgt[2] + first_acgt[3];
	for(j = 0; j < 4; j++){
		double ratio = first_acgt[j]/total;
		out << ratio << "\t" ;
	}
	out << endl;

	if(reads_file2 != ""){

	in.open(reads_file2.c_str(), ios::in);
	if(!in){
		cout << "ERROR: could not open " << reads_file2 << endl;
		exit(0);
	}//
	getline(in, line);
	if(line[0] != '@'){
		cout << "ERROR: the files with reads must be in FASTQ format" << endl;
		exit(0);
	}
	while(!in.eof()){
		getline(in, aread);
		getline(in, third);
		getline(in, aqual);
		long int len = aread.length();
		for(j = 0; j < len; j++){
			long int ind = hash_char(aread[j]);
			if(ind != -1 && ind < 4){
				second_acgt[ind] += 1.0;;
			}
		}//for j
		getline(in, line);
	}//while
	in.close(); in.clear();
	out << "ACGT-distribution for paired-end, from file with option -2:" << endl;
	total = second_acgt[0] + second_acgt[1] + second_acgt[2] + second_acgt[3];
	for(j = 0; j < 4; j++){
		double ratio = second_acgt[j]/total;
		out << ratio << "\t" ;
	}
	out << endl;

	}//if second file
	out.close(); out.clear();
	return 0;
}//main


