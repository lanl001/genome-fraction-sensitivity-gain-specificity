//============================================================================
// Name        : Evaluation.cpp
// Author      : llx
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <omp.h>
#include "Evaluation.h"
#include "parsingargs.h"

using namespace std;
using namespace __gnu_cxx;

namespace __gnu_cxx
{
size_t str_hash::operator()(const string& str) const
{
	return __stl_hash_string(str.c_str());
}
bool str_equal::operator()(const string& s1, const string& s2) const
{
	return s1 == s2;
}
}
hash_map<string, Cgenome, str_hash, str_equal> gnhm;
hash_map<string, Clongread, str_hash, str_equal> lrhm;
int TP;
int FP;
int TN;
int FN;
//const char * genomeFIFO = "/tmp/genometmp";
//const char * longreadFIFO = "/tmp/longreadtmp";
int numofthread;
int part;
pid_t pid;
string corrector;
clock_t total1;
clock_t total2;
clock_t total3;
clock_t total4;
int GetBaseNum(ifstream& file, fstream::pos_type head, fstream::pos_type tail)
{
	fstream::pos_type tempindex = file.tellg();
	char* s;
	string ss;
	char* p;
	int res = 0;
	s = (char*) malloc((tail - head + 2) * sizeof(char));
	s[tail - head + 1] = '\0';
	file.seekg(head, ios::beg);
	file.read(s, tail - head + 1);
	p = s;
	while ((*p) != '\0')
	{
		if ((*p) != '\n' && (*p) != ' ')
			res++;
		p++;
	}
	free(s);
	file.clear();
	file.seekg(tempindex, ios::beg);
	fstream::pos_type temp = file.tellg();
	return res;
}

void Hashgenome(ifstream& genomefile)
{
	char s[1001] = { 0 };
//	s[1000] = '\0';
	string ss;
	fstream::pos_type position = 0;
	unsigned long temp1 = 0, temp2 = 0;
	bool flag = false;  //if the reading of headline hasn't been completed
	bool flag2 = false;  //if the reading of conitgname hasn't been completed
	bool flag3 = false;
	string genomename;
	fstream::pos_type headindex;
	fstream::pos_type tailindex;
	fstream::pos_type tempindex;
	unsigned long gcount;
	while (genomefile.read(s, 1000), genomefile.gcount() != 0)
	{
		gcount = genomefile.gcount();
		s[gcount] = '\0';
		if (gcount != 1000)
		{
			genomefile.clear();
			genomefile.seekg(0, ios::end);
		}
		ss = s;
		if (flag2)
		{
			unsigned long p;
			if ((p = ss.find(' ', 0)) != string::npos)
			{
				genomename += ss.substr(temp1, p - temp1);
				flag2 = false;
			}
			else
			{
				if ((p = ss.find('\n', temp1)) != string::npos)
				{
					genomename += ss.substr(temp1, p - temp1);
					flag2 = false;
				}
				else
				{
					genomename += ss.substr(temp1, p - temp1);
					flag2 = true;
				}
			}
		}
		if (flag)
		{
			if ((temp2 = ss.find('\n', 0)) != string::npos)
			{
				position = (genomefile.tellg() - gcount + temp2);
				headindex = position;
				gnhm[genomename].index = position;
				flag = false;
			}
		}
		while ((temp1 = ss.find('>', temp2)) != string::npos)
		{
			if (flag3)
			{
				tailindex = (genomefile.tellg() - gcount - 1 + temp1);
				gnhm[genomename].length = GetBaseNum(genomefile, headindex, tailindex);
			}
			flag3 = true;
			unsigned long p;
			if ((p = ss.find(" ", temp1)) != string::npos)
			{
				genomename = ss.substr(temp1 + 1, p - temp1 - 1);
				flag2 = false;
			}
			else
			{
				if ((p = ss.find('\n', temp1)) != string::npos)
				{
					genomename = ss.substr(temp1 + 1, p - temp1 - 1);
					flag2 = false;
				}
				else
				{
					genomename = ss.substr(temp1 + 1, p - temp1 - 1);
					flag2 = true;
				}
			}
			if ((temp2 = ss.find('\n', temp1)) != string::npos)
			{
				position = (genomefile.tellg() - gcount + 1 + temp2);
				headindex = position;
				gnhm[genomename].index = position;
			}
			else
				flag = true;
		}
		temp1 = 0;
		temp2 = 0;
	}
	genomefile.clear();
	genomefile.seekg(0, ios::end);
	tailindex = genomefile.tellg() - 1;
	gnhm[genomename].length = GetBaseNum(genomefile, headindex, tailindex);
}

void Hashlongread(ifstream& longreadfile)
{
	char s[1001] = { 0 };
//	s[1000] = '\0';
	string ss;
	fstream::pos_type position = 0;
	unsigned long temp1 = 0, temp2 = 0;
	bool flag = false;  //if the reading of headline hasn't been completed
	bool flag2 = false;  //if the reading of conitgname hasn't been completed
	bool flag3 = false;
	string longreadname;
	fstream::pos_type headindex;
	fstream::pos_type tailindex;
	fstream::pos_type tempindex;
	unsigned long gcount;
	while (longreadfile.read(s, 1000), longreadfile.gcount() != 0)
	{
		gcount = longreadfile.gcount();
		s[gcount] = '\0';
		if (gcount != 1000)
		{
			longreadfile.clear();
			longreadfile.seekg(0, ios::end);
		}
		ss = s;
		if (flag2)
		{
			unsigned long p;
			if ((p = ss.find(' ', 0)) != string::npos)
			{
				longreadname += ss.substr(temp1, p - temp1);
				flag2 = false;
			}
			else
			{
				if ((p = ss.find('\n', temp1)) != string::npos)
				{
					longreadname += ss.substr(temp1, p - temp1);
					flag2 = false;
				}
				else
				{
					longreadname += ss.substr(temp1, p - temp1);
					flag2 = true;
				}
			}
		}
		if (flag)
		{
			if ((temp2 = ss.find('\n', 0)) != string::npos)
			{
				position = (longreadfile.tellg() - gcount + temp2);
				headindex = position;
				lrhm[longreadname].index = position;
				flag = false;
			}
		}
		while ((temp1 = ss.find('>', temp2)) != string::npos)
		{
			if (flag3)
			{
				tailindex = (longreadfile.tellg() - gcount - 1 + temp1);
				lrhm[longreadname].length = GetBaseNum(longreadfile, headindex, tailindex);
			}
			flag3 = true;
			unsigned long p;
			if ((p = ss.find(" ", temp1)) != string::npos)
			{
				longreadname = ss.substr(temp1 + 1, p - temp1 - 1);
				flag2 = false;
			}
			else
			{
				if ((p = ss.find('\n', temp1)) != string::npos)
				{
					longreadname = ss.substr(temp1 + 1, p - temp1 - 1);
					flag2 = false;
				}
				else
				{
					longreadname = ss.substr(temp1 + 1, p - temp1 - 1);
					flag2 = true;
				}
			}
			if ((temp2 = ss.find('\n', temp1)) != string::npos)
			{
				position = (longreadfile.tellg() - gcount + 1 + temp2);
				headindex = position;
				lrhm[longreadname].index = position;
			}
			else
				flag = true;
		}
		temp1 = 0;
		temp2 = 0;
	}
	longreadfile.clear();
	longreadfile.seekg(0, ios::end);
	tailindex = longreadfile.tellg() - 1;
	lrhm[longreadname].length = GetBaseNum(longreadfile, headindex, tailindex);
}

string getword(char* str, long &pos)
{
	string temp;
	while (str[pos] == '\t' || str[pos] == '\0' || str[pos] == '\n' || str[pos] == ' ')
	{
		pos++;
	}
	while (str[pos] != '\t' && str[pos] != '\0' && str[pos] != '\n' && str[pos] != ' ')
	{
		temp += str[pos];
		pos++;
	}
	return temp;
}

string Cfilebuffer::Getstring(char *argv, fstream::pos_type begin, fstream::pos_type end)
{

	string temp;
	if (strcmp(this->myfilename.c_str(), argv) == 0 && begin >= this->startoffset && begin <= this->startoffset + this->mysize - 1)
	{
		this->hited++;
		if (end <= this->startoffset + this->mysize - 1)
		{
//#pragma omp critical
//			cout << omp_get_thread_num() << ": end <= size " << "begin=" << begin << " end =" << end << " size = " << this->mysize << endl;
			temp.insert(0, buffer + begin - this->startoffset, (buffer + end - this->startoffset) - (buffer + begin - this->startoffset) + 1);
		}
		else if (end > this->startoffset + this->mysize - 1)
		{
//#pragma omp critical
//			cout << omp_get_thread_num() << ": end > size " << "begin=" << begin << " end =" << end << " size = " << this->mysize << endl;
			temp.insert(0, buffer + begin - this->startoffset, buffer + this->mysize - (buffer + begin - this->startoffset));
			temp += this->Getstring(argv, this->startoffset + this->mysize, end);
		}
	}
	else if (strcmp(this->myfilename.c_str(), argv) == 0 && (begin < this->startoffset || begin > this->startoffset + this->mysize - 1))
	{
		this->nothited++;
		ifstream file(argv);
		file.seekg(begin, ios::beg);
		file.read(this->buffer, this->expectsize);
		this->mysize = file.gcount();
		if (mysize == 0)
		{
#pragma omp critical
			{
				cout << omp_get_thread_num() << ": Err:begin pos is:'" << begin << "'which is larger than the file size ! end pos is : " << end << endl;
				file.clear();
				file.seekg(0, ios::end);
				streampos p = file.tellg();
				file.close();
				cout << "file size =" << p << endl;
				exit(-1);
			}
			string nullstring;
			return nullstring;
		}
		this->buffer[mysize] = '\0';
		file.close();
		this->startoffset = begin;
//#pragma omp critical
//		cout << omp_get_thread_num() << ": startoffset= " << begin << "size= " << mysize << endl;
		temp = this->Getstring(argv, begin, end);
	}
	else if (strcmp(this->myfilename.c_str(), argv) != 0)
	{
		string filename(argv);
		this->refreshbuffer(filename, this->mysize);
		temp = this->Getstring(argv, begin, end);
	}
	return temp;
}

string GetACut2(char *argv, fstream::pos_type position, int begin, int end, Cfilebuffer &filebuffer)
{
	string res;
	string temp;
	int remain = begin;
	fstream::pos_type p = position;
	int p2;
//	string::iterator it;
	while (remain != 0)
	{
		p2 = 0;
		temp = filebuffer.Getstring(argv, p, p + remain - 1);
		p = p + remain;
		remain = 0;
		while ((p2 = temp.find('\n', p2)) != string::npos)
		{
			remain++;
			p2++;
		}
		p2 = 0;
		while ((p2 = temp.find(' ', p2)) != string::npos)
		{
			remain++;
			p2++;
		}
	}
	remain = end - begin + 1;
	while (remain != 0)
	{
		p2 = 0;
		res += filebuffer.Getstring(argv, p, p + remain - 1);
		p = p + remain;
		remain = 0;
		while ((p2 = res.find('\n', p2)) != string::npos)
		{
			res.erase(p2, 1);
			remain++;
		}
		p2 = 0;
		while ((p2 = res.find(' ', p2)) != string::npos)
		{
			res.erase(p2, 1);
			remain++;
		}
	}
	return res;
}

string GetACut(char *argv, fstream::pos_type position, int begin, int end, Cfilebuffer &filebuffer)
{
	string res;
	string temp;
	int remain = begin;
	fstream::pos_type p = position;
	int p2 = 0;
	end += begin / 70;
	begin += begin / 70;
	remain = end - begin + 1;
	remain = remain + (remain + begin % 71) / 70;
	res = filebuffer.Getstring(argv, position + begin, position + begin + remain - 1);
	while ((p2 = res.find('\n', p2)) != string::npos)
	{
		res.erase(p2, 1);
	}
	p2 = 0;
	while ((p2 = res.find(' ', p2)) != string::npos)
	{
		res.erase(p2, 1);
	}
	return res;
}

void getErrsBeforeCorrection(string longreadname, stringstream& blasrresult, mapgroup& beforegrp, mapgroup& aftergrp)
{
	beforegrp.filename = "blasrresult.sam";
	/*ifstream beforecorrectionsamfile("blasrresult.sam");
	 char * beforebuffer;
	 filebuf *pbuf;
	 pbuf = beforecorrectionsamfile.rdbuf();
	 long size = pbuf->pubseekoff(0, ios::end, ios::in);
	 pbuf->pubseekpos(0, ios::in);
	 beforebuffer = new char[size];
	 pbuf->sgetn(beforebuffer, size);
	 long pos = 0;*/
	char c;
	long long counter = 0;
	string cigar, num;
	stringstream ss;
	long long temp;
	unsigned long long p = 0;
	bool validline = true;
	bool firstword = true;
	string gname;
	string readname;
	string tempstr;
	unsigned long long genomeindex;
	mapline line;
	while (c = blasrresult.get(), blasrresult)
	{
		if (firstword == true)
		{
			if (c == '@')
			{
				validline = false;
				firstword = false;
			}
			else
			{
				firstword = false;
				blasrresult >> readname;
				readname.insert(0, 1, c);
				size_t strpos = readname.find(longreadname);
				if (strpos == string::npos)
					validline = false;
				else
					line.longreadname = longreadname;
			}
		}
		else if (c == '	' && validline == true)
		{
			firstword = false;
			++counter;
			if (counter == 2)
			{
				blasrresult >> gname;
				if (gname == "*")
				{
					validline = false;
				}
				else
				{
					size_t strpos1 = gname.rfind('_');
					size_t strpos2;
					if (strpos1 != string::npos)
					{
						strpos2 = gname.rfind('_', strpos1 - 1);
						if (strpos2 != string::npos)
						{
							tempstr = gname.substr(strpos2 + 1, strpos1 - strpos2 - 1);
							ss.clear();
							ss << tempstr;
							ss >> genomeindex;
						}
					}
					else
						validline = false;
					line.genomename = gname.substr(0, strpos2);
				}
			}
			else if (counter == 3)
			{
				unsigned long long genomeindextemp;
				ss.clear();
				blasrresult >> tempstr;
				ss << tempstr;
				ss >> genomeindextemp;
				genomeindex += genomeindextemp;
				line.startoffset = genomeindex;
			}
			else if (counter == 5)
			{
				blasrresult >> cigar;
				errnode myerr;
				while (p < cigar.size())
				{
					while (cigar[p] >= '0' && cigar[p] <= '9')
					{
						num += cigar[p];
						p++;
					}
					if (cigar[p] == 'M')
					{
						cout << "'M' found in cigar. You need to run MDtoCigar first" << endl;
						exit(-1);
					}
					else if (cigar[p] == '=')
					{
						ss.clear();
						ss << num;
						ss >> temp;
						genomeindex += temp;
						p++;
						num.clear();
					}
					else if (cigar[p] == 'X')
					{
						ss.clear();
						ss << num;
						ss >> temp;

						for (int i = 0; i < temp; i++)
						{
							myerr.location = genomeindex++;
							myerr.err = substitution;
							line.errnodes.push_back(myerr);
						}
						p++;
						num.clear();
					}
					else if (cigar[p] == 'I')
					{
						ss.clear();
						ss << num;
						ss >> temp;
						myerr.location = genomeindex;
						for (int i = 0; i < temp; i++)
						{
							myerr.err = insertion;
							line.errnodes.push_back(myerr);
						}
						p++;
						num.clear();
						//						cout << "temp =" << temp << endl;
						//						cout << "totalnum+temp =" << totalnum << endl;
					}
					else if (cigar[p] == 'D')
					{
						ss.clear();
						ss << num;
						ss >> temp;
						for (int i = 0; i < temp; i++)
						{
							myerr.location = genomeindex++;
							myerr.err = deletion;
							line.errnodes.push_back(myerr);
						}
						p++;
						num.clear();
					}
					else if (cigar[p] == 'S' || cigar[p] == 'H' || cigar[p] == '*')
					{
						num.clear();
						p++;
					}
					else
					{
						cout << cigar[p] << endl;
					}
				}
				line.endoffset = genomeindex;
				beforegrp.maplines.push_back(line);
			}
		}
		else if (c == '\n')
		{
			counter = 0;
			p = 0;
			validline = true;
			firstword = true;
			if (c == '\0')
				break;
		}
		else
		{
			firstword = false;
		}
	}
	if (line.longreadname.compare("*") == 0)
	{
		line.longreadname = longreadname;
		line.genomename = "*";
		line.startoffset = 0;
		line.endoffset = 0;
		beforegrp.maplines.push_back(line);
		line.errnodes.clear();
	}
	else
	{
		vector<errnode>::iterator it1 = aftergrp.maplines.back().errnodes.begin();
		vector<errnode>::iterator it2 = line.errnodes.begin();
		int insertioncounter = 0;
		unsigned long long startoffset = line.startoffset < aftergrp.maplines.back().startoffset ? aftergrp.maplines.back().startoffset : line.startoffset;
		unsigned long long endoffset = line.endoffset > aftergrp.maplines.back().endoffset ? aftergrp.maplines.back().endoffset : line.endoffset;
		int tempTP = 0, tempFN = 0, tempFP = 0;
		while (it1 != aftergrp.maplines.back().errnodes.end() && it1->location < startoffset)
		{
			it1++;
		}
		while (it2 != line.errnodes.end() && it2->location < startoffset)
		{
			it2++;
		}
		while (1)
		{
			while (it2 != line.errnodes.end() && it2->location <= endoffset && (it1 == aftergrp.maplines.back().errnodes.end() || it1->location > endoffset || it2->location < it1->location))
			{
				if (it2->err == insertion)
					insertioncounter++;
				tempTP++;
				it2++;
			}
			if (it1 == aftergrp.maplines.back().errnodes.end() || (it2 != line.errnodes.end() && it2->location > endoffset) || it1->location > endoffset)
				break;
			while (it1 != aftergrp.maplines.back().errnodes.end() && it2 != line.errnodes.end() && it2->location == it1->location && it1->location <= endoffset && it2->location <= endoffset)
			{
				if (it1->err == insertion || it2->err == insertion)
					insertioncounter++;
				tempFN++;
				it1++;
				it2++;
			}
			while (it1 != aftergrp.maplines.back().errnodes.end() && it1->location <= endoffset && (it2 == line.errnodes.end() || it2->location > endoffset || it2->location > it1->location))
			{
				if (it1->err == insertion)
					insertioncounter++;
				tempFP++;
				it1++;
			}
			if (it2 == line.errnodes.end() || it2->location > endoffset || (it1 != aftergrp.maplines.back().errnodes.end() && it1->location > endoffset))
				break;
		}
		TN += endoffset - startoffset + 1 + insertioncounter - tempFP - tempTP - tempFN;
		FP += tempFP;
		TP += tempTP;
		FN += tempFN;
	}
	line.errnodes.clear();
}

void getErrsAfterCorrection(ifstream& aftercorrectionsamfile, const char* samfilename, char* gfilename, char* lrfilename)
{
	mapgroup beforegrp;
	mapgroup aftergrp;
	FILE* longreadtemp;
	FILE* genometemp;
	aftergrp.filename = samfilename;

	char * afterbuffer;
	char c;
	filebuf *pbuf;
	pbuf = aftercorrectionsamfile.rdbuf();
	long size = pbuf->pubseekoff(0, ios::end, ios::in);
	pbuf->pubseekpos(0, ios::in);
	afterbuffer = new char[size];
	pbuf->sgetn(afterbuffer, size);
	long pos = 0;
	long long counter = 0;
	string cigar, num;
	stringstream ss;
	long long temp;
	unsigned long long p = 0;
	bool validline = true;
	bool firstword = true;
	string gname;
	string readname;
	string tempstr;
	unsigned long long genomeindex;
	mapline line;
	Cfilebuffer filebuffer(gfilename, 1000);
	Cfilebuffer filebuffer2(lrfilename, 1000);
	long long filecount = omp_get_thread_num();
	string strfilecount;
	int headoffset = 0;
	int tailoffset = 0;
	while (c = afterbuffer[pos++], c != '\0')
	{
		if (firstword == true)
		{
			if (c == '@')
			{
				validline = false;
				firstword = false;
			}
			else
			{
				firstword = false;
				pos--;
				readname = getword(afterbuffer, pos);
				if (corrector.compare("LoRDEC") == 0 || corrector.compare("HALC") == 0 || corrector.compare("colormap") == 0)
				{
					size_t strpos = readname.rfind('_');
					readname = readname.substr(0, strpos);
				}
				else if (corrector.compare("ectools") == 0)
				{
					size_t strpos = readname.rfind('_');
					strpos = readname.rfind('_', --strpos);
					readname = readname.substr(0, strpos);
				}
				else if (corrector.compare("LSC") == 0)
				{
					size_t strpos = readname.rfind('|');
					readname = readname.substr(0, strpos);
				}
				else if (corrector.compare("proovread") == 0)
				{
					size_t strpos = readname.rfind('/');
					if ((strpos = readname.find('.', strpos)) != string::npos)
						readname = readname.substr(0, strpos);
				}
				else if (corrector.compare("proovread2") == 0)
				{
					size_t strpos = readname.find('.');
					strpos = readname.find('.', ++strpos);
					readname = readname.substr(0, strpos);
				}
				else if (corrector.compare("Jabba") == 0)
				{
					size_t strpos = readname.rfind('/');
					if((strpos = readname.find('_', strpos)) != string::npos)
						readname = readname.substr(0, strpos);
				}
				else if (corrector.compare("Jabba2") == 0)
				{
					size_t strpos = readname.find('.');
					strpos = readname.find('_', strpos);
					readname = readname.substr(0,strpos);
				}
				else
				{
					cout << "Invalid corrector" << endl;
					cout << "corrector can be LoRDEC/HALC/colormap/LSC/ectools/proovread/proovread2/Jabba/Jabba2" << endl;
					exit(-1);
				}
				line.longreadname = readname;
			}
		}
		else if (c == '	' && validline == true)
		{
			firstword = false;
			++counter;
			if (counter == 2)
			{
				gname = getword(afterbuffer, pos);
				if (gname == "*")
				{
					validline = false;
				}
				else
				{
					line.genomename = gname;
				}
			}
			else if (counter == 3)
			{
				ss.clear();
				tempstr = getword(afterbuffer, pos);
				ss << tempstr;
				ss >> genomeindex;
				line.startoffset = genomeindex;
			}
			else if (counter == 5)
			{
				cigar = getword(afterbuffer, pos);
				errnode myerr;
				bool head = true;
				headoffset = 0;
				tailoffset = 0;
				while (p < cigar.size())
				{
					while (cigar[p] >= '0' && cigar[p] <= '9')
					{
						num += cigar[p];
						p++;
					}
					if (cigar[p] == 'M')
					{
						cout << "'M' found in cigar. You need to run MDtoCigar first" << endl;
						exit(-1);
					}
					else if (cigar[p] == '=')
					{
						ss.clear();
						ss << num;
						ss >> temp;
						genomeindex += temp;
						p++;
						num.clear();
					}
					else if (cigar[p] == 'X')
					{
						ss.clear();
						ss << num;
						ss >> temp;

						for (int i = 0; i < temp; i++)
						{
							myerr.location = genomeindex++;
							myerr.err = substitution;
							line.errnodes.push_back(myerr);
						}
						p++;
						num.clear();
					}
					else if (cigar[p] == 'I')
					{
						ss.clear();
						ss << num;
						ss >> temp;
						myerr.location = genomeindex;
						for (int i = 0; i < temp; i++)
						{
							myerr.err = insertion;
							line.errnodes.push_back(myerr);
						}
						p++;
						num.clear();
						//						cout << "temp =" << temp << endl;
						//						cout << "totalnum+temp =" << totalnum << endl;
					}
					else if (cigar[p] == 'D')
					{
						ss.clear();
						ss << num;
						ss >> temp;
						for (int i = 0; i < temp; i++)
						{
							myerr.location = genomeindex++;
							myerr.err = deletion;
							line.errnodes.push_back(myerr);
						}
						p++;
						num.clear();
					}
					else if (cigar[p] == 'S' || cigar[p] == 'H'||cigar[p] == '*')
					{
						num.clear();
						p++;
					}
					else
					{
						num.clear();
						p++;
						cout << "\033[33m" << "Warining:Unknown cigar:" << cigar[p] << "\033[0m" << endl;
					}
				}
				line.endoffset = genomeindex;
				aftergrp.maplines.push_back(line);

				ss.clear();
				ss << filecount;
				ss >> strfilecount;
				ss.clear();
				string strpid;
				ss << pid;
				ss >> strpid;
				string tempgenomename = "/dev/shm/" + strpid + "_trashme_genometemp_" + strfilecount + ".fa";

				genometemp = fopen(tempgenomename.c_str(), "w");
				if (!genometemp)
				{
					cout << " Cannot create genometemp.fa " << endl;
					exit(-1);
				}
				if (gnhm.find(gname) == gnhm.end())
				{
					cout << "ERROR:cannot find genome in hashmap" << endl;
					exit(-1);
				}
				fprintf(genometemp, ">%s_%d_%d\n", gname.c_str(), ((signed) line.startoffset - 201 > 0 ? (signed) line.startoffset - 201 : 0),
						((signed) line.endoffset + 198 < (gnhm[gname].length - 1) ? (signed) line.endoffset + 198 : (gnhm[gname].length - 1)));
				clock_t start2 = clock();

				fprintf(genometemp, "%s\n",
						GetACut(gfilename, gnhm[gname].index, ((signed) line.startoffset - 201 > 0 ? (signed) line.startoffset - 201 : 0),
								((signed) line.endoffset + 198 < (gnhm[gname].length - 1) ? (signed) line.endoffset + 198 : (gnhm[gname].length - 1)), filebuffer).c_str());
				total3 += clock() - start2;
				fclose(genometemp);

				string templongreadname = "/dev/shm/" + strpid + "_trashme_longreadtemp_" + strfilecount + ".fa";

				longreadtemp = fopen(templongreadname.c_str(), "w");
				if (!longreadtemp)
				{
					cout << " Cannot create longreadtemp.fa " << endl;
					exit(-1);
				}
				if (lrhm.find(readname) == lrhm.end())
				{
					cout << "ERROR:cannot find longread in hashmap :" << readname << endl;
					exit(-1);
				}
				fprintf(longreadtemp, ">%s\n", readname.c_str());
				clock_t start3 = clock();
				fprintf(longreadtemp, "%s\n", GetACut(lrfilename, lrhm[readname].index, 0, lrhm[readname].length - 1, filebuffer2).c_str());
				total4 += clock() - start3;
				fclose(longreadtemp);

				stringstream blasrresult;
				string blasrcommand = "blasr " + templongreadname + " " + tempgenomename + " -sam -cigarUseSeqMatch -minMatch 8 -nCandidates 20 -bestn 1 -maxScore 20000 2>>errs"+strfilecount.c_str()+".log";
				clock_t start = clock();
				string errmsg;
				if (!myexec(blasrcommand.c_str(), blasrresult))
				{
					cout << "\033[31m" << "fail to execute blasr" << endl;
					cout << errmsg << "\033[0m" << endl;
					exit(0);
				}
				getErrsBeforeCorrection(readname, blasrresult, beforegrp, aftergrp);
				/*				string haha;
				 while (blasrresult >> haha)
				 {
				 cout << haha;
				 }
				 cout << endl;*/
				total1 += (clock() - start);
			}
		}
		else if (c == '\n')
		{
			counter = 0;
			p = 0;
			validline = true;
			firstword = true;
			line.errnodes.clear();
			if (c == '\0')
				break;
		}
		else
		{
			firstword = false;
		}
	}
	aftercorrectionsamfile.close();
	delete[] afterbuffer;
}

void parameteranalysing(int argc, char* argv[])
{
	/*
	 * parameter analyzing
	 */
	string tmpPara = "";
	for (int i = 4; i < argc; i++)
	{
		if (strlen(argv[i]) == 0)
		{
			cout << "find NULL" << endl;
			tmpPara += char(31);
		}
		else
		{
			tmpPara += argv[i];
		}
		tmpPara += " ";
	}
	std::map<std::string, std::vector<std::string> > result;
	ParsingArgs pa;
	pa.AddArgType('t', "threads", ParsingArgs::MUST_VALUE);
	pa.AddArgType('p', "part", ParsingArgs::MAYBE_VALUE);
	pa.AddArgType('c', "corrector", ParsingArgs::MUST_VALUE);
	std::string errPos;
	int iRet = pa.Parse(tmpPara, result, errPos);
	if (0 > iRet)
	{
		cout << "Invalid parameters!" << endl << iRet << errPos << endl;
		system("cat readme.txt");
		exit(-1);
	}
	else
	{
		if (result.find("p") == result.end() && result.find("part") == result.end())
		{
			cout << "parameter -p is necessary" << endl;
			exit(-1);
		}
		if (result.find("c") == result.end() && result.find("corrector") == result.end())
		{
			cout << "parameter -c is necessary" << endl;
			exit(-1);
		}
		map<std::string, std::vector<std::string> >::iterator it = result.begin();
		int argflag = 0;
		for (; it != result.end(); ++it)
		{

			if (it->first.compare("t") == 0 || it->first.compare("threads") == 0)
			{
				if (it->second.size() != 1)
				{
					cout << "Invalid parameters!" << iRet << errPos << endl;
					system("cat readme.txt");
					exit(-1);
				}
				else
				{

					std::stringstream ss;
					ss << it->second[0];
					ss >> numofthread;
					omp_set_num_threads(numofthread);
					cout << "threads = " << numofthread << endl;
				}
			}

			if (it->first.compare("p") == 0 || it->first.compare("part") == 0)
			{
				if (it->second.size() != 1)
				{
					cout << "Invalid parameters!" << iRet << errPos << endl;
					system("cat readme.txt");
					exit(-1);
				}
				else
				{

					std::stringstream ss;
					ss << it->second[0];
					ss >> part;
					cout << "part = " << part << endl;
				}
			}
			if (it->first.compare("c") == 0 || it->first.compare("corrector") == 0)
			{
				if (it->second.size() != 1)
				{
					cout << "Invalid parameters!" << iRet << errPos << endl;
					system("cat readme.txt");
					exit(-1);
				}
				else
				{
					corrector = it->second[0];
					cout << "corrector =" << corrector << endl;
				}
			}
		}
	}
}

int main(int argc, char *argv[])
{
	pid = getpid();
	numofthread = 0;
	parameteranalysing(argc, argv);
	total2 = clock();
	ios::sync_with_stdio(false);
	ifstream genomefile;
	ifstream longreadfile;
	stringstream ss;
	string afterfilename;
	long long filecounter = 0;
	long long filenum;
	TP = 0;
	FP = 0;
	TN = 0;
	FN = 0;
	cout << "\033[33m" << "warning: all fasta files should be 70 bases per line. If not, run wordsperline first" << "\033[0m" << endl;
	cout << "\033[31m" << "have you loaded blasr??????????????????" << endl;
	genomefile.open(argv[2]);
	longreadfile.open(argv[3]);
	Hashgenome(genomefile);
	Hashlongread(longreadfile);
#pragma omp parallel
	{
		ifstream aftercorrectionsamfile;
#pragma omp master
		{
			cout << omp_get_num_threads() << " threads created" << endl;
		}
#pragma omp for  private(afterfilename,aftercorrectionsamfile,ss)
		for (long long i = 0; i < part; i++)
		{
			ss.clear();
			ss << argv[1] << "_" << i << ".sam";
			ss >> afterfilename;
			aftercorrectionsamfile.open(afterfilename.c_str());
			if (!aftercorrectionsamfile)
			{
				cout << "cannot open samfile:" << afterfilename.c_str() << endl;
				exit(-1);
			}
			getErrsAfterCorrection(aftercorrectionsamfile, afterfilename.c_str(), argv[2], argv[3]);
			aftercorrectionsamfile.close();
		}
	}
	total2 = clock() - total2;
	cout << "TP= " << TP << endl; // prints !!!Hello World!!!
	cout << "FP= " << FP << endl;
	cout << "FN= " << FN << endl;
	cout << "TN= " << TN << endl;
	cout << total3 << endl;
	cout << total4 << endl;
	cout << total1 << endl;
	cout << total2 << endl;
	return 0;
}

