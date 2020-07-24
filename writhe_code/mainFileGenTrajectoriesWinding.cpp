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
	 std::vector<std::vector<point> > points1;
   std::vector<std::vector<point> > points2;
	 int noFrames = std::atoi(argv[2]);
	 int noAtoms = std::atoi(argv[3]);
   // note that this is the number of atoms in one frame
	 int startIndex = std::atoi(argv[4]);
	 // first read in and dump all the frames up to startIndex
	 if (myfile.is_open()){ 
	   for(int i=1;i<startIndex;i++){
	     for(int j=0;j<2*noAtoms;j++){
	       std::getline(myfile,output);
	     }
	   }
	   // now read in the required frames
	   for(int i=1;i<=noFrames;i++){
	     std::vector<point> frame1;
	     std::vector<point> frame2;
	     for(int k=0;k<2*noAtoms;k++){
	       std::getline(myfile,output);
	       if(k<noAtoms){
		 frame1.push_back(point(output));
	       }else{
		 frame2.push_back(point(output));        
	       }
	     }
	     std::reverse(frame2.begin(),frame2.end());
	     points1.push_back(frame1);  
	     points2.push_back(frame2);
	     //std::cout<<frame1.size()<<" "<<frame2.size()<<"\n";
	   }
	   // now calculate the set of writhe values
	   myfile.close();
	   std::ofstream ofile;
	   ofile.open(argv[5]); 
	   for(int i=0;i<noFrames;i++){
	     if(points1[i].size()>3){
	       check =true;
	     }
	     // std::cout<<"frame "<<i<<"\n";
	     int index =startIndex+i;
	     if(check){
	       localWrithe lw2;
	       double windOutput;
	        /*************************

                   winding yo 

	       **********************/
	       if(argc >6){
	     	 // here the routine will be smoothed
	     	 windOutput = lw2.getWinding(points1[i],points2[i]);
	       }else{
	     	 windOutput = lw2.getWinding(points1[i],points2[i]);
	       }
         std::cout<<windOutput<<"\n";
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
