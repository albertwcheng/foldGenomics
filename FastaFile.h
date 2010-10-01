/***************************************************************************
 Copyright 2010 Wu Albert Cheng <albertwcheng@gmail.com>
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 *******************************************************************************/ 


#ifndef FASTAFILE_H_
#define FASTAFILE_H_

#include <fstream>
using namespace std;


class FastaFile
{
public:
	ifstream fin;
	string seqName;
	string seq;
	std::streampos last;
	bool upperCase;
	
	inline FastaFile(string filename,bool _upperCase=true): upperCase(_upperCase)
	{
		
		fin.open(filename.c_str());
		if(!fin.good())
		{
			fin.close();
			cerr<<"file "<<filename<<" cannot be opened"<<endl;
		}
	}
	
	inline bool readEntry()
	{
		seq="";
		seqName="";
		if(fin.eof())
			return false;
		
		fin>>seqName;

		if(seqName[0]!='>')
		{
			cerr<<"bad format!"<<endl;
			return false;
		}
		seqName=seqName.substr(1);
		last=fin.tellg();
		string line=" ";
		
		
		while(!fin.eof())
		{
			line="";
			fin>>line;
			//if(fin.eof())
			//	break;
			
			if(line[0]=='>')
			{
				fin.seekg(last);
				break;
				
			}
			
			if(upperCase)
				seq+=toUpperAndNoSpecialDNA(line);
			else
				seq+=line;
			last=fin.tellg();
		}
		
		return seq!="";
	}
	
	inline ~FastaFile()
	{
		fin.close();
	}
	
};

#endif /*FASTAFILE_H_*/
