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

	string mes("\nUSAGE: rev_compl -r <references> -s <input-reads> -o <output> [OPTIONS]\n\nOptions:\n");
	string singles(    "  -s <input-reads-file>      file with single reads (or queries) in FASTQ format");
	string output(     "  -o <output-file>");


	cout << mes << endl;
	cout <<  singles << endl
        	<< output << endl;



}//usage()

int parse_options(int argc, char* argv[], string &reads_file1, string &output)
{
	int res = 0;

	long int i, j;
	for(i = 1; i < argc; i++){
		if(argv[i][0] != '-')
			return -1;
		switch(argv[i][1]){
		case 's': reads_file1 = argv[++i]; break;
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

	if(argc < 3){
		usage();
		return 0;
	}//if
	string reads_file1(""),  output;
	int res = parse_options( argc, argv, reads_file1, output);
	if(res == -1){
		usage();
		exit(0);
	}

	ofstream out;
	out.open(output.c_str(), ios::out);
	ifstream in;
	in.open(reads_file1.c_str(), ios::in);
	if(!in){
		cout << "ERROR: could not open " <<reads_file1 << endl;
		exit(0);
	}
	string aread;
	long int st, en;
	in >> aread;
	while(!in.eof()){
		in >> st;
		in >> en;
		long int len = aread.length();
		j = len - 1;
		string rc("");
		for(; j >= 0; j--){
			if(aread[j] == 'A' || aread[j] == 'a')
				rc = rc + 'T';
			else if(aread[j] == 'C' || aread[j] == 'c')
				rc = rc + 'G';
			else if(aread[j] == 'G' || aread[j] == 'g')
				rc = rc + 'C';
			else if(aread[j] == 'T' || aread[j] == 't')
				rc = rc + 'A';
			else
				rc = rc + aread[j];
		}//for j
		out << rc << "\t" << en << "\t" << st << endl;
		in >> aread;
	}//while
	in.close(); in.clear();
	
	
	out.close();
	return 0;
}//main


