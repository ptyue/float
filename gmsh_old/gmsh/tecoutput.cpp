/* output mesh data in tecplot format
 * This file can be modified to generate mesh data for ALE code directly.
 * Currently only implemented in 2D.
 * 6/23/2014, P. Yue
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Gmsh.h"
#include "GModel.h"
#include "GEntity.h"
#include "MElement.h"
#include "MVertex.h"
#include "GmshMessage.h"
using namespace std;

// output 2d mesh data 
int writeTecplotMESH2D(GModel *m, const std::string& name)
{
  // put mesh data into entities. These entities correpsonds to the points, 
  // lines, and faces used to define the geometry for mesh generation. 
  std::vector<GEntity*> entities;
  m->getEntities(entities);

  // get the number of all vertices (including those not in the mesh),
  // and index the vertices in a continuous sequence
  const int nGNodes=m->indexMeshVertices(true);  
  //const int nGNodes=m->getNumMeshVertices();
  int *index2node=new int[nGNodes+1];	//working array, destroy at return
  int nElements=0,nVertices=0, nNodes=0;
  int nEntities=entities.size();
  int i,j,k;
  MElement *me;
  MVertex *mv;
  int nElem,nElemNodes,idx;
  const int nElmVert=3;		//Three vertices for a triangular element.
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
        nElemNodes=me->getNumVertices();
        for (k=0;k<nElemNodes;k++)
        {
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
  // At here, index2nodes[idx]=0 for dangling nodes, =-1 for vertices,
  // and =-2 for midnodes
  //cout<<"****"<<nVertices<<" "<<nNodes<<endl;
  // creat data file for output
  std::ofstream outFile(name.c_str());
  if(!outFile.good()){
    Msg::Error("Unable to open file '%s'", name.c_str());
    return 0;
  }
  
  // tecplot data file header
  outFile<<"Title = \"Gmsh grid data file\""<<std::endl;
  outFile<<"Variables = \"X\", \"Y\""<<std::endl;
  outFile<<"Zone T=\"T1\","<<"DATAPACKING=POINT,";
  outFile<<"NODES="<<nVertices<<",ELEMENTS="<<nElements<<",ZONETYPE=FETRIANGLE"<<std::endl;
  
  //output vertices coordinates
  int numv=0;
  int numn=nVertices;
  double *x=new double[nNodes],*y=new double[nNodes];
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
/*******************************************************************************
  // The following is another way to order the mesh nodes. The ordering is more
  // sporadic and thus abandoned
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
          mv=me->getVertex(k);
          idx=mv->getIndex();
          if(index2node[idx]>0)continue;
          if(k<nElmVert)
          {  
            numv++;
            index2node[idx]=numv;
            x[numv-1]=mv->x();
            y[numv-1]=mv->y();
          }
          else
          {
            numn++;
            index2node[idx]=numn;
            x[numn-1]=mv->x();
            y[numn-1]=mv->y();

          }
                       
        }
      }
    }
  }  
*******************************************************************************/
  for (i=1;i<=nVertices;i++)
  {
    outFile<<x[i-1]<<"  "<<y[i-1]<<endl;
  }

  //output elements
  for (i=0;i<nEntities;i++)  
  {
    if(entities[i]->dim()==2)
    {
      nElem=entities[i]->getNumMeshElements();
      for (j=0; j<nElem;j++)
      {
        me=entities[i]->getMeshElement(j);
        //for (k=0;k<me->getNumVertices();k++)
        for (k=0;k<nElmVert;k++)        
        {
          outFile<<" "<<index2node[me->getVertex(k)->getIndex()];
        }
        outFile<<std::endl;
      }
    }
  }
  
  
/*
  for (i=0;i<nEntities;i++)
  {
    if(entities[i]->dim()==1)
    {
      std::cout<<"Entity "<<i<<":";
      for (j=0;j<entities[i]->getNumMeshElements();j++)
      {
        MElement *me=entities[i]->getMeshElement(j);
        for ( k=0;k<me->getNumVertices();k++)
        {
          std::cout<<" "<<index2node[me->getVertex(k)->getIndex()];
        }
        std::cout<<";";
      }
      std::cout<<"\n";
    }
    std::cout<<"Entity "<<i<<": vertices"<<entities[i]->mesh_vertices.size()<<entities[i]->getTypeString();
    std::cout<<",vertices "<<entities[i]->getNumMeshVertices()<<",elements "<<entities[i]->getNumMeshElements()<<"\n";
  }
*/
  
  
  outFile.close();
  delete[] index2node	;
  return 1;
}
  
  
    
