#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;
int main()
{
vector <long double> xcord;
vector <long double> ycord;
vector <long double> zcord;
vector <double> cord(3);
int n=101;  
string type1;
int check = 0;
int mol=17113;
ifstream infile("8pah0.001M_nvt22-12.xyz");
if(!infile)
{
cout << "File not exist!" << endl;
}
else
{cout << "Start importing data!" << endl;
for (int i=1; i<=n; i++)
{
string line1;
getline(infile, line1);
cout << line1 << endl;
string line2;
getline(infile, line2);
cout << line2 << endl;
for (int j=1; j<=mol; j++)
{ infile >> type1 >> cord[0] >> cord[1] >> cord[2]; 
  if (type1=="Na") 
  {
  	check++;
  xcord.push_back(cord[0]);
  ycord.push_back(cord[1]);
  zcord.push_back(cord[2]);
}
  /*cout << type1 << " "  << cord[0] << " " << cord[1] << " " << cord[2] << endl;*/
}
string line3;
getline(infile, line3);
}
cout << "End importing data!" << endl;
}
infile.close();

int num=check;
ofstream outfile1("Na.txt");
if (outfile1.is_open())
{for (int k=0; k < num; k++)
{
    outfile1 << " " << xcord[k] << " " << ycord[k] << " " << zcord[k] << endl;
}
outfile1.close();
}
}