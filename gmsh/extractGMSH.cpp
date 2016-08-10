/******************************************************************************
 *  extractGMSH2D                                                             *
 *  This program extracts the 2D mesh information for numerical simulations.  *
 *  It's developed for triangular mesh, however, it should work for all types *
 *  of 2D meshes except hybrid ones.                                          *
 *  input:                                                                    *
 *      m:      pointer to the gmesh model                                    *
 *      ngmx:   max number of nodes per element (=6 for 2nd order triangles)  *
 *      ncmx:   max number of edges per element (=3 for triangles)            *
 *      max***: max number, used to check the array bounds
 *  output:                                                                   *
 *      nVertices:  Number of vertices in the mesh                            *
 *      nNodes:     Number of nodes in the mesh                               *
 *                  If the element order is higher than 1, these nodes include*
 *                  both vertices and edge nodes.                             *
 *      nElements:  Number of elements in the mesh                            *
 *      x,y:        coordinates of the mesh nodes. All the midnodes are placed*
 *                  after the vertices.                                       *
 *      inod:       element-node table                                        *
 *      nec:        element neighbor table                                    *
 *      ic(nic+1):  starting position of boundary sections                    *
 *      ibdnod(nbd):global nodal index (not the index in gmsh) of each        *
 *                  boundary nodes                                            *
 *      nside(nbound): number of boundary sections in each connected boundary *
 *                                                                            *
 *  Original code,  7/8/2014  P. Yue                                          *
 *  nside added,    1/4/2015  P. Yue                                          *
 ******************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include "Gmsh.h"
#include "GModel.h"
#include "GEntity.h"
#include "MElement.h"
#include "MVertex.h"
#include "GmshMessage.h"
using namespace std;

class MyEdge{
  public:
    int v1,v2;
    MyEdge(int a=0,int b=0){v1=min(a,b);v2=max(a,b);}
};

class MyEdgeLessThan{
  public:
    bool operator()(const MyEdge & edge1,const MyEdge& edge2) const
    {
      if (edge1.v1<edge2.v1)
      {
        return true;
      }
      else if(edge1.v1>edge2.v1)
      {
        return false;
      }
      else
      {
        return edge1.v2<edge2.v2;
      }
    }
};

class EdgeElem{
  public:
    // This edge is the sid1-th edge of element elm1 and the sid2-th edge of 
    // element elm2. If one of elm* is zero, this edge is a boundary edge.
    // These integers start from 0 for efficiency in C++.
    int elm1,sid1,elm2,sid2;
    EdgeElem(int a, int b, int c=-1, int d=-1){elm1=a; sid1=b; elm2=c; sid2=d;}
    void setElm2(int c, int d){elm2=c;sid2=d;}
};


int extractGMSH(GModel *m, const int ngmx, const int ncmx, 
    const int maxVertices, const int maxNodes, const int maxElements, 
    const int maxnic, const int maxnbd,         
    int &nVertices, int &nNodes, int &nElements, int &nic, int &nbd, int &nbound,
    double *x, double *y, int *inod, int *nec, int *ic, int *ibdnod, int *nside)
{
  // put mesh data into entities. These entities correpsonds to the points, 
  // lines, and faces used to define the geometry for mesh generation. 
  std::vector<GEntity*> entities;
  m->getEntities(entities);

  // get the number of all vertices (including those not in the mesh),
  // and index the vertices in a continuous sequence
  const int nGNodes=m->indexMeshVertices(true);  
  // temporary variables
  int *index2node=new int[nGNodes+1];	 //working array, destroy at return
  int nEntities=entities.size();
  int i,j,k,l,num,numv,numn;
  MElement *me;
  MVertex *mv;
  int nElem,idx;
  const int nElmVert=ncmx;  //# of vertices in one element  
  // These constants will be ouput to calling programs
  nElements=0; 
  nVertices=0; 
  nNodes=0;   
  /****************************************************************************
   * Extract mesh vertices, nodes, and elements.                              *
   * x,y[nNodes]: nodal coordinates. x,y[0:nVertices-1] for  vertices and     *
   *            x,y[nVertices:nNodes] for midnodes                            * 
   * inod[nElements,ngmx]: element-node table                                 *
   ****************************************************************************/  
  for (i=1;i<=nGNodes;i++)
  {
    index2node[i]=0;
  }
  // Count the number of mesh vertices, mesh nodes, and mesh elements.
  // (The entity vertices that do not appear in mesh are not counted.)
  for (i=0;i<nEntities;i++)  
  {
    // only counts the elements of dimension 2. 
    if(entities[i]->dim()==2)
    {
      nElem=entities[i]->getNumMeshElements();
      nElements +=nElem;
      for (j=0; j<nElem;j++)
      { 
        me=entities[i]->getMeshElement(j);
        for (k=0;k<me->getNumVertices();k++)
        { 
          if(me->getNumVertices()>ngmx)
          {
            cout<<"Error in meshing. ngmx should be >= "<<me->getNumVertices()
                <<endl;
            exit(1);
          }
          idx=me->getVertex(k)->getIndex();
          if(!index2node[idx])
          {
            if(k<nElmVert)
            {
              nVertices++;
              index2node[idx]=-1;  //vertices
            }
            else
            {
              index2node[idx]=-2; //midnodes
            }
            nNodes++;
          }
        }
      }    
    }
  }

  //Check the size of nElements, nNodes,  and nElements
  if(nVertices>maxVertices||nNodes>maxNodes||nElements>maxElements)
  {
    cout<<"Error in function extraGMSH: the following requirements should be met"<<endl;
    cout<<"maxVertices>="<<nVertices<<endl;
    cout<<"maxNodes>="<<nNodes<<endl;
    cout<<"maxElements>="<<nElements<<endl;
    exit(1);
  }

  //obtain vertices coordinates
  numv=0;
  numn=nVertices;
  for (i=0;i<nEntities;i++)  
  {
    for (j=0; j<entities[i]->mesh_vertices.size();j++)
    {
      mv=entities[i]->mesh_vertices[j];
      idx=mv->getIndex();
      if(index2node[idx]==-1)
      {
        numv++;
        index2node[idx]=numv;
        x[numv-1]=mv->x();
        y[numv-1]=mv->y();
      }
      else if(index2node[idx]==-2)
      {
        numn++;
        index2node[idx]=numn;
        x[numn-1]=mv->x();
        y[numn-1]=mv->y();
      }
    }
  }  
  //obtain the element-node table
  //For higher order elements, getNumVertices()>3 and the midnodes are also
  //included in this table.
  num=0;
  for (i=0;i<nEntities;i++)  
  {
    if(entities[i]->dim()==2)
    {
      nElem=entities[i]->getNumMeshElements();
      for (j=0; j<nElem;j++)
      {
        me=entities[i]->getMeshElement(j);
        for (k=0;k<me->getNumVertices();k++)
        {
          //outFile<<" "<<index2node[me->getVertex(k)->getIndex()];
          inod[num*ngmx+k]=index2node[me->getVertex(k)->getIndex()];
        }
        num++;
      }
    }
  }
  /****************************************************************************
   * Extract element connection information                                   *
   * nec[nElements,ncmx]:  nec[i-1,j-1] stores the element index of the       *
   *        j-th neighbor (sharing the j-th edge) of the i-th element.        *
   *        The range of element index is [1,nElements].                      *
   ****************************************************************************/   
  //obtain the edge list. A map is used to store the edges. The two vertices 
  //of each edge are used as the key and the value contains two elements 
  //sharing this edge.
  map<MyEdge,EdgeElem,MyEdgeLessThan> edges;
  int v1,v2;
  map<MyEdge,EdgeElem>::iterator it;
  for (i=0;i<nElements;i++)
  {
    for (j=0;j<nElmVert;j++)
    {
      v1=inod[i*ngmx+j];
      v2=inod[i*ngmx+(j+1)%nElmVert];
      MyEdge edg(v1,v2);
      // Search whether this edge already exists in the map.
      it=edges.find(edg);
      if(it==edges.end())
        // The newly inserted edge is the (j+1)th edge of the (i+1)th element.
        edges.insert(pair<MyEdge,EdgeElem>(edg,EdgeElem(i,j)));
      else
        // This existing edge is also the (j+1)th edge of the (i+1)th element.
        it->second.setElm2(i,j);
    }
  }
  for (it=edges.begin();it!=edges.end();it++)
  {
    EdgeElem &edg=it->second;
    // nec=0 for boundary segments (since elm* are set to -1 by default)
    // 1<=nec(i,j)<=nElements, to be consistenet with Fortran
    if(edg.elm1>=0)
      nec[(edg.elm1)*ncmx+edg.sid1]=edg.elm2+1;
    if(edg.elm2>=0)
      nec[(edg.elm2)*ncmx+edg.sid2]=edg.elm1+1;
  }    
  /****************************************************************************
   * Extract boundary                                                         *
   * nbd: 	number of boundary nodes                                      *
   * nic:       number of boundary sections                                   *
   * ibdnod(nbd):  ibdnod(i-1) stores the mesh nodes index (in[1,nNodes]) of  *
   *               the i-th boundary node.                                    *
   * ic(nic+1):    the i-th boundary section starts at the ic(i-1)-th boundary*
   *               node. ic(0)=1, ic(nic)=nbd+1.                              *
   ****************************************************************************/    
  //get the physical groups. groups[i] is the physical group of dimesion i.
  map<int, vector<GEntity*> > groups[4];
  m->getPhysicalGroups(groups);
  //We only need the physical groups of 1 dimension for the boundaries
  map<int, vector<GEntity*> > &group(groups[1]);
  const int nPhysical=group.size();
  //
  nic=nPhysical;
  if(nic>maxnic)
  {
    cout<<"Error in function extraGMSH: maxic should be at least"<<endl;
    cout<<"maxnic>="<<nic<<endl;
    exit(1);
  }  
  
  int nseg=0;
  int nseg0=0;
  nbound=0;
  num=0;
  //Note that the members of group are referenced by their tags which start
  //from 1 and goes continuously to nPhysical.
  for (i=1;i<=nPhysical;i++)
  { 
    ic[nseg]=num+1;
    nseg++;
    for (j=0;j<group[i].size();j++)
    {
      GEntity *line=group[i][j];
      for (k=0;k<line->getNumMeshElements();k++)
      {
        MElement *me=line->getMeshElement(k);
        for (l=0;l<me->getNumVertices();l++)
        {
          // for higher order edges the edge vertices are organized as
          //   0--2--3--...-1
          // 1 is shared with next edge and is thus taken care of by next edge.
          if(l==1)continue;
          if(num>=maxnbd)
          {
            cout<<"Error in function extraGMSH: maxnbd should be at least"<<endl;
            cout<<"maxnbd>="<<num+1<<endl;
            exit(1);
          }
          ibdnod[num]=index2node[me->getVertex(l)->getIndex()];
          num++;
        }
      }
    }
    GEntity *lined=group[i].back();
//    GEntity *lined=group[i][group[i].size()-1];
    MElement *med=lined->getMeshElement(lined->getNumMeshElements()-1);
    if(index2node[med->getVertex(1)->getIndex()]==ibdnod[ic[nseg0]-1])
    {
      nside[nbound]=nseg-nseg0;
      nbound++;
      nseg0=nseg;
    }
  }
  nbd=num;
  ic[nseg]=num+1;     

  
  delete[] index2node;
  return 1;
}
  
  
    
