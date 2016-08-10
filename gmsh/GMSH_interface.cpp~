/******************************************************************************
 *  fucntion gmsh2d_  (fortran gmsh2d)                                        *
 *  2D triangular mesh generation. To minizie bandwidth of sparse matrix, the *
 *  elements needs to be reordered (eg., by Cuthill & McKee's algorithm).     *
 *  Input:                                                                    *
 *      order:          order of the element (eg. 3 nodes for a first order   *
 *                      triangle and 6 nodes for 2nd order)                   *
 *      bgm:            0, no background mesh                                 *
 *                      1, background mesh for element size is used, need to  *
 *                      specify bgmfile                                       *
 *      geofile(char*): geometry file name *.geo                              *
 *      bgmfile(char*): background mesh file name *.pos                       *
 *      ngmx:           maximum number of nodes in each element               *
 *      mcmx:           maximum number of edges in each element               *
 *      maxVertices:    maximum number of vertices                            *
 *      maxNodes:       maximum number of nodes                               *
 *      maxnic:         maximum number of boundary sections                   *
 *      maxnbd:         maximum number of boundary nodes                      *
 *  Output (mesh):                                                            *
 *      nVertices:      number of vertices                                    * 
 *      nNodes:         number of nodes                                       *
 *      nElements:      number of elements                                    *
 *      nic:            number of boundary sections                           *
 *      nbd:            number of boundary nodes                              *
 *      nbound:         number of connected boundaries                        *
 *      x,y(nNodes):    coordinates of mesh nodes                             *
 *      inod(ngmx,nElements):   element description table                     *
 *      nec(ncmx,nElements):    neighboring elements table                    *
 *      ic(nic+1):      position (in ibdnod) of the starting node in each     *
 *                      section                                               *
 *      ibdnod(nbd):    global nodal indices of boundary nodes                *
 *      nside(nbound):  number of boundary sections in each connected boundary*
 *                                                                            *
 *  Original code, 12/2014, P. Yue                                            *
 ******************************************************************************/

#include "Gmsh.h"		//GmshBatch
#include "GModel.h"
#include "MElement.h"
#include "MVertex.h"
#include "GEntity.h"
#include "CommandLine.h"	//GmshInitialize
#include "Options.h"		//GetOptions
#include "Context.h"		//CTX class
#include <iostream>
#include <string>
#include <cstring>
#include <map>

using namespace std;
int writeTecplotMESH2D(GModel* m, const std::string& name); 
int extractGMSH(GModel *m, const int ngmx, const int ncmx, 
    const int maxVertices, const int maxNodes, const int maxElements, 
    const int maxnic, const int maxnbd,
    int &nVertices, int &nNodes, int &nElements, int &nic, int &nbd, int &nbound,
    double *x, double *y, int *inod, int *nec, int *ic, int *ibdnod, int *nside);
extern "C"
{
// 1)By default, the function names in c++ file must be in lower case with 
//   "_" appended. Since fortran is case insensitive (before Fortran 95), 
//   all function names are complied into lower case.
// 2)By default, fortran only supports pass by reference, thus all dummy 
//   are defined as pointers.
void gmsh2d_(int *order, int *bgm, char* geofile, char* bgmfile, 
    int *ngmx, int *ncmx,
    int *maxVertices, int *maxNodes, int *maxElements,int *maxnic,int *maxnbd,
    int *nVertices, int *nNodes, int *nElements, int *nic, int *nbd, int *nbound,
    double *x, double *y, int *inod, int *nec, int *ic, int *ibdnod, int *nside);
}  

void gmsh2d_(int *order, int *bgm, char *geofile, char *bgmfile, 
    int *ngmx, int *ncmx,
    int *maxVertices, int *maxNodes, int *maxElements,int *maxnic,int *maxnbd,
    int *nVertices, int *nNodes, int *nElements, int *nic, int *nbd, int *nbound,
    double *x, double *y, int *inod, int *nec, int *ic, int *ibdnod, int *nside)
{
  // assemble options for gmsh
  int argc;
  char args[10][100],*argv[10];
  strcpy(args[0],"gmsh");
  strcpy(args[1],geofile);
  strcpy(args[2],"-2");
  strcpy(args[3],"-order");
//  strcpy(args[4],(char[2]){char('0'+*order),'\0'});
  args[4][0]=char('0'+*order);
  args[4][1]='\0';

  if(*bgm)
  {
    strcpy(args[5],"-bgm");
    strcpy(args[6],bgmfile);
    argc=7;
  }
  else
  {
    argc=5;
  }
  for (int i=0;i<argc;i++) argv[i]=args[i];
  
  // Note declaration of this pointer won't invoke the constructor of GModel
  GModel *m;
  //regular 2D mesh generation
  //At least one GModel is required for to hold the data
  // If GModel *m=new GModel() is used, then m=GModel::current() automatically.
  if(GModel::list.empty())      
  {
    //regular 2D mesh generation, first time calling gmsh
    //At least one GModel is required for to hold the data
    // If GModel *m=new GModel() is used, then m=GModel::current() automatically.      
    new GModel();
    //Initialization 
    GmshInitialize(argc, argv);
  }
  else
  {
    //Gmsh has been called before and there is something left in the cache
    //Destroy the current model, otherwise the size of GModel::list will 
    //increase by 1 at each call of GmshBatch. This is more economical than
    //deleting the model.
    GModel::current()->destroy();
    //Meshing has been performed before. Some parameters in static class
    //CTX::instance() need to be cleaned up.
    CTX::instance()->files.clear();
    CTX::instance()->bgmFileName.clear();
    //Initialize() doesn't work here; GetOptions needs to be called to 
    //update option parameters.
    GetOptions(argc,argv);  
  }
  //Mesh Generation
  GmshBatch();
  m=GModel::current();
  //ouput mesh (debugging purpose)
//  writeTecplotMESH2D(m,"tecmesh.dat");
  //Extract mesh information
  extractGMSH(m, *ngmx, *ncmx, 
              *maxVertices, *maxNodes, *maxElements, *maxnic, *maxnbd,         
              *nVertices, *nNodes, *nElements, *nic, *nbd, *nbound,
              x, y, inod, nec, ic, ibdnod, nside);

  GmshFinalize();  
}    

