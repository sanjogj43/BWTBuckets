#include<stdio.h>
#include<fstream>
#include<sstream>
#include<iostream>
#include<string>
#include<vector>
#include<map>
#include<math.h>
#include<conio.h>

using namespace std;


class BWT
{
public:
	vector<char> BWTString;
	vector<short unsigned int> LCPVal;
	vector<unsigned int> componentIds;

	int numOfBuckets;
	int compSize;
	int kmerLength;
	int numEltsInEachBucket; // Bucket size.

	string origString;
	vector<string> LCPArray;
	void QuickSort(int start, int n);
	int Partition(int start, int n);
	void swap(unsigned int &s1,unsigned int &s2);
	void findSuperMaximalRepeats();
	void findBWT();
	void findLCPArray();
	void fillUpComponentIds(int bucketId);
	unsigned int convertCharacter(char);
	unsigned int getKmerMask();
};
void BWT::findBWT()
{
	for (int i = 0; i < componentIds.size(); i++)
	{
		if (componentIds[i]>0)
			BWTString.push_back(origString[componentIds[i] - 1]);
		else
			BWTString.push_back(origString[origString.size() - 1]);
	}
}

void BWT::findLCPArray()
{
	LCPArray.push_back("\0");
	LCPVal.push_back(-1);
	for (int i = 1; i < componentIds.size(); i++)
	{
		string s1 = origString.substr(componentIds[i-1]);
		string s2 = origString.substr(componentIds[i]);
		int j = 0;
		string LCPStr="";
		while (s1[j] == s2[j])
		{
			LCPStr+=s1[j];
			j++;
		}
		LCPStr[j] += '\0';
		LCPArray.push_back(LCPStr);
		LCPVal.push_back(j);

	}
}

void BWT::QuickSort(int start, int n)
{
	if(n > 1){
		int pivotIndex = Partition(start, n);
		int n1 = pivotIndex - start;
		int n2 = n - n1 - 1;
		QuickSort(start, n1);
		QuickSort(pivotIndex + 1, n2);
	}
}

int BWT::Partition(int start, int n)
{
	string pivot = origString.substr(componentIds[start],compSize);
	int i = start + 1;
	int j = start + n - 1;
	while (i <= j){
		string iComp = origString.substr(componentIds[i], compSize);
		string jComp = origString.substr(componentIds[j], compSize);
		string iPivot = pivot, jPivot = pivot;
		int pivotPoint = componentIds[start];
		int iInc = componentIds[i], jInc=componentIds[j];
		while (iComp == iPivot)
		{
			iInc += compSize;
			pivotPoint += compSize;

			iComp = origString.substr(iInc, compSize);
			iPivot = origString.substr(pivotPoint, compSize);
		}
		//iPivot = pivot;
		pivotPoint = componentIds[start];
		while (jComp == jPivot)
		{
			jInc += compSize;
			pivotPoint += compSize;
			jPivot = origString.substr(pivotPoint, compSize);
			jComp = origString.substr(jInc, compSize);
		}
		if (iComp < iPivot)
			i++;
		else if (jComp > jPivot)
			j--;
		else{
			swap(componentIds[i], componentIds[j]);
			i++;
  			j--;
		}//else
	}//while
	swap(componentIds[start], componentIds[j]);
	return j;
}//end Partition

unsigned int BWT::convertCharacter(char c)
{
	switch (c)
	{
		case 'A':
			return 0x00;
			break;
		case 'C':
			return 0x01;
			break;
		case 'G':
			return 0x02;
			break;
		case 'T':
			return 0x03;
			break;
		default:
			return 0xffffffff;
			break;
	}
}

void BWT::swap(unsigned int &s1,unsigned int &s2)
{
	unsigned int temp = s1;
	s1 = s2;
	s2 = temp;
}

void BWT::fillUpComponentIds(int bucketId)
{
	componentIds.clear();
	LCPVal.clear();
	LCPArray.clear();
	unsigned int kmer = 0;
	unsigned int mask = getKmerMask();
	int minIndex = bucketId*numEltsInEachBucket;
	int maxIndex = (bucketId + 1)*numEltsInEachBucket-1;
	for (int i = 0; i < origString.length()-kmerLength; i++)
	{
		for (int j = 0; j < kmerLength;j++)
		{
			kmer <<= 2;
			kmer &= mask;
			unsigned int convChar = convertCharacter(origString[i+j]);
			if (convChar == 0xffffffff)
			{
				cout << "i=" << i << endl;
				cout << "j=" << j << endl;
				cout << "origString[" << i + j << "]=" << origString[i + j]<<endl;
				throw "Invalid character encountered!!";
			}
			kmer |= convChar;
		}
		if (kmer >= minIndex && kmer<maxIndex)
		{
			componentIds.push_back(i);
		}

	}
}

unsigned int BWT::getKmerMask()
{
	unsigned int x = 0xffffffff;
	x <<= kmerLength;
	x = ~x;
	return x;
}

void BWT::findSuperMaximalRepeats()
{
	fstream fout;
	fout.open("out.txt", ios::out);
	bool currentUp = false, currDown = false;
	int startInt = -1, endInt = -1;

	for (int i = 0; i < LCPVal.size() - 1; i++)
	{
		//if (!currentUp && bwt.LCPVal[i+1] < 3)
		//break;
		if (i + 1 != LCPVal.size() && LCPVal[i]<LCPVal[i + 1])
		{
			currentUp = true;
			startInt = i;
			endInt = i + 1;
		}
		if (currentUp)
		{
			if (LCPVal[i] == LCPVal[i + 1])
			{
				endInt = i + 1;
			}
			else if (LCPVal[i] > LCPVal[i + 1])
			{
				currentUp = false;
				currDown = true;
			}
		}
		if (!currentUp && currDown)
		{
			//put stint and endint in file.
			if (endInt - startInt + 1 <= 4 && LCPVal[endInt]>15) // possiblity for supermaximal repeat 4 and 15
			{
				// check for pairwise distinct
				bool pairWiseDistinct = true;
				for (int j = startInt; j < endInt; j++)
				{
					for (int k = j + 1; k <= endInt; k++)
					{
						if (BWTString[j] == BWTString[k])
						{
							pairWiseDistinct = false;
							break;
						}
					}
				}

				if (fout.is_open() && pairWiseDistinct)
				{
					fout << componentIds[startInt] << "\t" << LCPVal[endInt] << "\t" << origString.substr(componentIds[startInt], LCPVal[endInt]) << endl;
					cout << endl << "Supermaximal repeats:" << endl;
					cout << componentIds[startInt] << "\t" << LCPVal[endInt] << "\t" << origString.substr(componentIds[startInt], LCPVal[endInt]) << endl;
					pairWiseDistinct = false;
					currentUp = false;
					currDown = false;
				}
			}
		}
	}
	fout.close();
}

int main()
{ 
	string s1 = "ATATATTAG$";//"mississippi$";
	string s =s1;
	fstream myFile;
	stringstream ss;
	myFile.open("genome.txt");
	if (myFile.is_open())
	{
		string line;
		int i = 0;
		while (getline(myFile, line))
		{
			ss<<line;
			i++;
		}

		//cout << endl << i << endl;
		myFile.close();
	}
	
	s = ss.str()+"$";
	cout << s;
	BWT bwt;
	bwt.origString = s;
	int compSize = 0;
	cout << "Enter the component size:";
	cin >> bwt.compSize; 
	cout << "Enter Kmer length (<=15):";
	cin >> bwt.kmerLength;
	bwt.numOfBuckets = (bwt.compSize*bwt.compSize);
	bwt.numEltsInEachBucket = s.size() / (bwt.compSize*bwt.compSize);
	
	try
	{
		for (int i = 0; i < bwt.numOfBuckets; i++)
		{
			// fill up the components
			bwt.fillUpComponentIds(i);

			// Sort the components
			bwt.QuickSort(0, bwt.componentIds.size());
			// find BWT from sorted component ids
			bwt.findBWT();
			// find LCPs from sorted component ids
			bwt.findLCPArray();
			// Start : compute Super maximal repeats 
			bwt.findSuperMaximalRepeats();
			// End : compute Super maximal repeats 
		}
	}
	catch (const char* msg)
	{
		cout << "Exception : " << msg << endl;
	}
	cout << "over";
	_getch();
	return 0;

}