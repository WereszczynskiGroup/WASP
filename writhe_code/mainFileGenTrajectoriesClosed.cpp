#include "localWrithe.h"
#include <time.h>
#include <string.h>


int main( int argc, const char* argv[] )
{
       if(argc >=2){
	 bool check=false;
	 std::ifstream myfile;
	 myfile.open(argv[1]);
	 std::string output;
	 std::vector<std::vector<point> > points;
	 int noFrames = std::atoi(argv[2]);
	 int noAtoms = std::atoi(argv[3]);
	 int startIndex = std::atoi(argv[4]);
	 // first read in and dump all the frames up to startIndex
	 if (myfile.is_open()){ 
	   for(int i=1;i<startIndex;i++){
	     for(int j=1;j<=noAtoms;j++){
	       std::getline(myfile,output);
	     }
	   }
	   // now read in the required frames
	   for(int i=1;i<=noFrames;i++){
	     std::vector<point> frame;
	     for(int k=0;k<noAtoms;k++){
	       std::getline(myfile,output);
	       frame.push_back(point(output));
	     }
	     points.push_back(frame);
	   }
	   // now calculate the set of writhe values
	   myfile.close();
	   std::ofstream ofile;
	   ofile.open(argv[5]); 
	   for(int i=0;i<noFrames;i++){
	     if(points[i].size()>3){
	       check =true;
	     }
	     //std::cout<<"frame "<<i<<"\n";
	     int index =startIndex+i;
	     if(check){
	       /*************************

                   star writhes

	       **********************/
	       localWrithe lw;
	       localWrithe lw2;
	       std::vector<double> writheOutput;
	       std::vector<point> copyPt = points[i];
	     	 writheOutput = lw.getWritheGenClosed(points[i]);
	       ofile<<writheOutput[0]<<" "<<writheOutput[1]<<" "<<writheOutput[2]<<" {";
	       if(writheOutput.size()>3){
	     	 for(int j=3;j<writheOutput.size()-1;j++){
	     	   ofile<<writheOutput[j]<<" ";
	     	 }
	     	 ofile<<writheOutput[writheOutput.size()-1]<<"} "<<index<<"\n";
	       }else{
	     	 ofile<<"} "<<index<<"\n";
	       }
	     }else{
	       ofile<<index<<" curve too short: < 3 points\n";
	     }
	   }
	   ofile.close();   
	   
	 }
	 else{
	     std::cout<<"Curve data file failed to open";
	 }
       }else{
	 std::cout<<"must supply a file containing the curve data";
       }
       return 0;
}
