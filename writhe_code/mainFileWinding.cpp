#include "localWrithe.h"
#include <time.h>


int main( int argc, const char* argv[] )
{
	if(argc ==3){
	std::ifstream myfile1;
	std::ifstream myfile2;
 	myfile1.open(argv[1]);
	myfile2.open(argv[2]);
	std::string output1;
	std::vector<point> points1;
	if (myfile1.is_open()) { 
	  while(std::getline(myfile1,output1)){
		points1.push_back(point(output1));
  		}
	}else{
	std::cout<<"Curve data file 1 failed to open";
	}
	myfile1.close();
       	std::string output2;
	std::vector<point> points2;
	if (myfile2.is_open()) { 
	  while(std::getline(myfile2,output2)){
		points2.push_back(point(output2));
  		}
	}else{
	std::cout<<"Curve data file 2 failed to open";
	}
	myfile2.close();
	if(points1.size()>= 3 && points2.size()>=3){
	localWrithe lw;
	double winding;
	winding  = lw.getWinding(points1,points2);
	std::cout<<winding<<"\n";
	}
	}else{
	  std::cout<<"you must supply two curve files, i.e. ./winding file1.dat file2.dat ";
	}
	return 0;
}
