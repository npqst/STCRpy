#include <cstdio>
#include <stdlib.h>
#include <cstring>
#include <map>
#include <vector>
#include <cmath>
using namespace std;

struct Coor3f {
  double x;
  double y;
  double z;
};
      
vector<Coor3f> atom;
vector<char*> lines;
      
typedef float (Matrix)[4][3];
     
int main(int argc, char **argv) {
  char buf[2000];
  char *filename = argv[2];
  char *matfilename = argv[1];
  int n;
   
  FILE *matfile = fopen(matfilename, "r");
  if (matfile == NULL) {
    fprintf(stderr,"Matrix file does not exist\n");
     return 1;
  }
 
  float avgx, avgy, avgz; 
  fgets(buf, 2000, matfile);
  sscanf(buf, "%f %f %f", &avgx, &avgy, &avgz);
  
  Matrix m;
  for (int i = 0; i < 4; i++) {
    fgets(buf, 2000, matfile);
    if (sscanf(buf, "%f %f %f", m[i], m[i]+1, m[i]+2) != 3) {
      printf("%s", buf);
      fprintf(stderr, "Reading error in matrix file %s\n", matfilename);
      return 1;
    }
  }
  fclose(matfile);
                                    
  FILE *fil = fopen(filename, "r");
  char code[10];
  while (!feof(fil)) {
  if (!fgets(buf, 2000, fil)) break;
  sscanf(buf, "%s", code);
  if (!strncmp(code,"ATOM", 4)) {
    Coor3f ccoor;
    ccoor.x = atof(buf+30);
    ccoor.y = atof(buf+38);
    ccoor.z = atof(buf+46);
    atom.push_back(ccoor);
    char *x = new char[strlen(buf) + 1];
    strcpy(x, buf);
    lines.push_back(x);
  }
  }
   for (n = 0; n < atom.size(); n++) {
     atom[n].x -= avgx;
     atom[n].y -= avgy;
     atom[n].z -= avgz;
   }
                         
   //Applying matrix
   for (n = 0; n < atom.size(); n++) {
     Coor3f &c = atom[n];
     Coor3f cc;
     cc.x = c.x * m[1][0] + c.y * m[1][1] + c.z * m[1][2];
     cc.y = c.x * m[2][0] + c.y * m[2][1] + c.z * m[2][2];
     cc.z = c.x * m[3][0] + c.y * m[3][1] + c.z * m[3][2];
     c = cc;
                                                
     c.x += m[0][0];
     c.y += m[0][1];
     c.z += m[0][2];
   }
   //Printing PDB file
   for (n = 0; n < atom.size(); n++) {
     printf("%30.30s%8.3f%8.3f%8.3f%s", lines[n], atom[n].x, atom[n].y, atom[n].
     z, lines[n]+ 54);
   }
}                                                                        
/*
'FILENAME==ARGV[1]{v[3*NR-3] = $1; v[3*NR-2]=$2;v[3*NR-1]=$3}  FILENAME==A
RGV[2] && $1 == "ATOM" {x0 = substr($0,31,8)-v[0]; y0 = substr($0,39,8)-v[1]; z0
 = substr($0,47,8)-v[2]; x = v[6]*x0+v[7]*y0+v[8]*z0+v[3]; y = v[9]*x0+v[10]*y0+
 v[11]*z0+v[4]; z = v[12]*x0+v[13]*y0+v[14]*z0+v[5]; printf("%s", substr($0,1,30)
 ); printf("%8.3f%8.3f%8.3f",x,y,z); printf("%s\n",substr($0,55,7) "10.00" substr
 ($0,67))}'*/ 


 
