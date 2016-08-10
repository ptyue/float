#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "GR_config.h"
#include "GR_AdaptPred.h"
#include "GR_Bdry2D.h"
#include "GR_InsertionQueue.h"
#include "GR_Mesh2D.h"
#include "GR_events.h"
#include "GR_misc.h"
#ifndef MAX
#define MAX(a,b) (a > b ? a : b)
#endif
#ifndef MIN
#define MIN(a,b) (a < b ? a : b)
#endif

extern "C" {
void  grummp_init_(int *mxnbd, int *mxvert,int *mxelem,int *ngmx,int *ncmx,
			int *nbd,int *nvert,int *nelem,int *refseg,int *bdseg,
			int *inod,int *nec,int *reft,double *coorx,double *coory,
			double *h,double *dG);
void  grummp_adp_(int *mxnbd, int *mxvert,int *mxelem,int *ngmx,int* ncmx,
		int *nvertd, double *phid, int *nbd,int *nvert,int *nelem, 
		int *bdseg,int *refseg,	int *inod,int *nec,int *reft,
		double *coorx,double *coory,double *h1,double *h2,double *h3,
		int *lrmsh,double *dG,int *lwrmesh);
	}
			
	



static bool GRUMMP_Inited=false;

bool MeshADP_2D(Mesh2D *pM2D,double *phid,double *h1,double *h2,double *h3,double dGrading,int iRUNS,int IR);
double h_length(double x,double y,double *,double *);
double phi(double x,double y,int iRUNS);
double h_length_smallest();
bool thesameedge(int na1,int na2,int nb1,int b2);
void check_clockwise(double *x1,double *x2,double *x3,int iCell);
void ReturnToSolver(Mesh2D *pM2D,int *nbd,int *nvert,int *nelem,int *ngmx,int *ncxm,
		int *refseg,int *bdseg,int *inod,int *nec,int *reft,double *coorx,double *coory);
void MarkVertLevel(Mesh2D *pM2D,double *phid);		
         

/* This subroutine to generate the initial 2D-mesh.
*/
void  grummp_init_(int *mxnbd, int *mxvert,int *mxelem,int *ngmx,int *ncmx,
			int *nbd,int *nvert,int *nelem,int *refseg,int *bdseg,
			int *inod,int *nec,int *reft,double *coorx,double *coory,
			double *h,double *dG)
{	
	double dGrading = *dG;
	double dScale = 1.0;
	char strBaseFileName[80];
	char strCanonicalName[80];
	enum eEncroachType eET = eLens;
	int iNumPasses;
	Mesh2D *pM2D;
	int debug=false;
	
	if(debug) printf("\n======== iMessageStdoutLevel = %d before setting.\n\n",iMessageStdoutLevel);
	
	sprintf(strBaseFileName,"amphi");
	sprintf(strCanonicalName,"tri");
	
	
	iMessageStdoutLevel=0;
	vGRUMMPInit(strCanonicalName);
	GRUMMP_Inited=true;
	if(debug) printf("\n======== iMessageStdoutLevel = %d after setting.\n\n",iMessageStdoutLevel);
	
	//vOpenMessageFile(strBaseFileName);

/* */
      Bdry2D B2D(strBaseFileName, false); 
      

//	dScale=2.0/(*h);
//	if(debug)printf("the dScale = %f \n", dScale);
//	printf("========before creating or adapting the mesh\n");
	
	  pM2D=new Mesh2D(B2D,2);

	  pM2D->vCheckForSmallAngles();
        pM2D->vSetEncroachmentType(eET);
        pM2D->vCheckForLargeAngles();
	  pM2D->vInitLengthScale(1.0/dScale, 1.0/dGrading);
        pM2D->vSetLSChecking(true);
        
      
        InsertionQueue IQ(pM2D);

        IQ.vQualityRefine();

        pM2D->vEvaluateQuality();
        pM2D->vSetSmoothingThreshold(40.);
        iNumPasses=0;
        vWriteFile_Mesh2D(*pM2D,"amphi");
	  delete pM2D;


// refining based on the given mesh size.
	  pM2D= new Mesh2D(strBaseFileName,2);	
        pM2D->vCheckForSmallAngles();
 
  	 if(debug){
 	 	FILE *fp=fopen("lengthscale_beforeRG.txt","w");
       	for(int iV=0;iV<pM2D->iNumVerts();iV++) {
	  		fprintf(fp,"The initial length scale of node %4d is %f\n",
	  				 iV+1, pM2D->pVVert(iV)->dLS());
	  	}
	  	fclose(fp);
	}
 
        pM2D->vInitLengthScale(1.0/dScale, 1.0/dGrading);
 
  	 if(debug){
 	 	FILE *fp=fopen("lengthscale_afterRG.txt","w");
       	for(int iV=0;iV<pM2D->iNumVerts();iV++) {
	  		fprintf(fp,"The initial length scale of node %4d is %f\n",
	  				 iV+1, pM2D->pVVert(iV)->dLS());
	  	}
	  	fclose(fp);
	}
 
 
        pM2D->vSetLSChecking(true);
	  if(debug)printf("The uniform length scale for the init mesh = %f\n", *h);
	  for(int iV=0;iV<pM2D->iNumVerts();iV++) {
	  	pM2D->pVVert(iV)->LSset=false;
	  	if(debug)printf("The initial length scale of node %4d is %f\n",
	  				 iV+1, pM2D->pVVert(iV)->dLS());
	  }
/*	  for(int iBE=0;iBE<pM2D->iNumBdryFaces();iBE++){
	  	pM2D->pVVert(pM2D->iVertIndex(pM2D->pBFBFace(iBE)->pVVert(0)))->LSset=true;
	  	pM2D->pVVert(pM2D->iVertIndex(pM2D->pBFBFace(iBE)->pVVert(1)))->LSset=true;
	  	}
*/	  for(int iV=0;iV<pM2D->iNumVerts();iV++){	 	
		Vert *pV=pM2D->pVVert(iV);
			pV->vSetLS(*h);
			pM2D->vMarkForLengthScaleCheck(pV);
	 }
	 
  	 if(debug){
 	 	FILE *fp=fopen("lengthscale_beforeUpdate.txt","w");
       	for(int iV=0;iV<pM2D->iNumVerts();iV++) {
	  		fprintf(fp,"The initial length scale of node %4d is %f\n",
	  				 iV+1, pM2D->pVVert(iV)->dLS());
	  	}
	  	fclose(fp);
	}

//	 dGrading = *dG;
//	 dGrading = 3;
//	 pM2D->vInitLengthScale(1.0/dScale, 1.0/dGrading);

  	 pM2D->vUpdateLengthScale();

  	 if(debug){
 	 	FILE *fp=fopen("lengthscale_afterUpdate.txt","w");
       	for(int iV=0;iV<pM2D->iNumVerts();iV++) {
	  		fprintf(fp,"The initial length scale of node %4d is %f\n",
	  				 iV+1, pM2D->pVVert(iV)->dLS());
	  	}
	  	fclose(fp);
	}


 	 pM2D->vAdaptToLengthScale(eET);
 	 if(debug){
 	 	FILE *fp=fopen("lengthscale_afterAdapt.txt","w");
       	for(int iV=0;iV<pM2D->iNumVerts();iV++) {
	  		fprintf(fp,"The initial length scale of node %4d is %f\n",
	  				 iV+1, pM2D->pVVert(iV)->dLS());
	  	}
	  	fclose(fp);
	}
   
      
      // Smoothing parameters for later
       pM2D->vSetSmoothingThreshold(40.);
 
        
 // return all the Mesh variables. 
 
	  *nvert=pM2D->iNumVerts();
	  *nelem=pM2D->iNumCells();
	  *nbd=pM2D->iNumBdryFaces();
	  
	  if(*nbd>*mxnbd||*nvert>*mxvert||*nelem>*mxelem){
	  	printf("Vital Error.....\n ========== Mesh is too large ===========\n");
	  	printf(" ========== More memory needed ===========\n");
	  	exit(1);
	}
	
	  ReturnToSolver(pM2D,nbd,nvert,nelem,ngmx,ncmx,refseg,bdseg,inod,nec,reft,coorx,coory);
	  
//    write out the mesh  
        vWriteFile_Mesh2D(*pM2D,strBaseFileName);
        iNumPasses = 0;
        delete pM2D;
	
 	return ;
}



/* this subroutine adapting the mesh from the given mesh file amphi.mesh

       subroutine rfnmsh(mxnbd,mxvert,mxelem,nvertd,phid,
     &                  nbd,nvert,nelem,bdseg,refseg,
     &                  inod,nec,coor,h1,h2,lrmsh)
     input:
        amphi.bdry: file specifying the domain
        old mesh info: presaved in data file or memory by GRUMMP
        maximum numbers: (for array bound checking)
           mxnbd,mxvert,mxelem                 
        old mesh related variables: 
           nvertd:  number of vertices
           phid:    phi field

        mesh size:
           h1:   mesh size along the interface(phi=0 contour)
           h2:   mesh size at infinity phi=1
           h3:   mesh size at infinity phi=-1           
    output:
        new mesh variables:
           nbd:     number of boundary nodes
           nvert:   number of vertices
           nelem:   number of elements(cells)
           bdseg:   end nodes of boundary segments
           refseg:  reference number of boundary segments
           inod:    element description table
           nec:     neighbouring elements info
           coor:    coordinates of vertices.
        lrmsh:   whether remesh has been performed
                 T, remeshed; F, the old mesh is good enough, remesh
                 not performed.

*/
void  grummp_adp_(int *mxnbd, int *mxvert,int *mxelem,int *ngmx,int* ncmx,
		int *nvertd, double *phid, int *nbd,int *nvert,int *nelem, 
		int *bdseg,int *refseg,	int *inod,int *nec,int *reft,
		double *coorx,double *coory,double *h1,double *h2,double *h3,
		int *lrmsh,double *dG,int *lwrmesh)
{	

	double dGrading;
	double dScale = 1.0;
	char strBaseFileName[80];
	char strCanonicalName[80];
	bool qLengthScaleGiven=false;
	enum eEncroachType eET = eLens;
	int iNumPasses;
	Mesh2D *pM2D;
	bool debug=false;
	
	
	sprintf(strBaseFileName,"amphi");
	sprintf(strCanonicalName,"meshopt2d");
	

	
	if(!GRUMMP_Inited){
	   vGRUMMPInit(strCanonicalName);
	   GRUMMP_Inited=true;
	   }

	bool qAllowBdryChanges = true;
	int iQualMeasure = 2, iPreviousPasses = 0;	
	
	iMessageStdoutLevel=0;
//	vOpenMessageFile(strBaseFileName);
	
      static Mesh2D M2D(strBaseFileName,2);
      static int iRUNS = 0;
      static int iREMESH = 0;
	pM2D=&M2D;
//	pM2D =new Mesh2D(strBaseFileName,2);

//	only write the mesh and return 
	if(*lwrmesh) {
	    sprintf(strBaseFileName,"amphi");
	    vWriteFile_Mesh2D(*pM2D,strBaseFileName);
	    return;			
		}
	
	dGrading = *dG;
	iRUNS=iRUNS+1;
	*lrmsh=MeshADP_2D(pM2D,phid,h1,h2,h3,dGrading,iRUNS,iREMESH);
	if(*lrmsh) iREMESH = iRUNS;
// return the value to the solver 	  
	if(*lrmsh){
 
	    *nvert=pM2D->iNumVerts();
	    *nelem=pM2D->iNumCells();
	    *nbd=pM2D->iNumBdryFaces();
	    
	    if(*nbd>*mxnbd||*nvert>*mxvert||*nelem>*mxelem){
	    	printf("Vital Error.....\n ========== Mesh is too large ===========\n");
	    	printf(" ========== More memory needed ===========\n");
	    	exit(1);
	      }

  	    ReturnToSolver(pM2D,nbd,nvert,nelem,ngmx,ncmx,refseg,bdseg,inod,nec,reft,coorx,coory);
//	    sprintf(strBaseFileName,"amphi",iRUNS);
//	    vWriteFile_Mesh2D(*pM2D,strBaseFileName);
  	
  	}
		
//	delete pM2D;
	return ;
}

// end of calling function of adapting from Fortran.
//=========================================================================

/* this function to Adapt the mesh pM2D from the parameter:
return true if remeshing performed.
phid: phi on the old mesh
h1: the lengthscale near the interface
h2: the lengthscale for the coarse mesh size. 
iRUNS: the number that MeshADP_2D is called
iREMESH: the number that Adaption was performed last time.
*/
bool MeshADP_2D(Mesh2D *pM2D,double *phid,double *h1,double *h2,double *h3,
			double dGrading,int iRUNS,int iREMESH)
{
//	 iRUNS = iRUNS + 1; // iRUNS to trace the times that the function is called. 

	 bool debug = false;
	 bool lrmsh;  // Flag to indicate if the mesh shall be remeshed.
	 char strBaseFileName[80];
	 enum eEncroachType eET = eLens;
//	 enum eEncroachType eET = eNewLens;
	 FILE *fp;
	 double dScale=1.;
	 
// Flag to indicate if the Level for the vertices shall be relabeled.	 
	 static bool lChLevel = true; 

	 assert(pM2D->qValid());
//	 pM2D->vEvaluateQuality();
	 
	 pM2D->vAllowBdryChanges();	 
	 pM2D->vCheckForSmallAngles();
//	 pM2D->vCheckForLargeAngles();
	 pM2D->vInitLengthScale(1./dScale, 1./dGrading);
	 pM2D->vSetLSChecking(true);
	 
/*  	 if(debug&&iRUNS==1){
 	 	fp=fopen("lengthscale_afterInit.txt","w");
       	for(int iV=0;iV<pM2D->iNumVerts();iV++) {
	  		fprintf(fp,"The initial length scale of node %4d is %f\n",
	  				 iV+1, pM2D->pVVert(iV)->dLS());
	  	}
	  	fclose(fp);
	}
*/	 
	 
	 Vert *pV;
	 Face *pE;
	 int iV;
	 int iE;
	 int iCell;
	 int VL,VR;
	 double VLx,VRx,VLy,VRy;
	 double lenL,lenR,lenEdge;



// if remeshing performed last time, then reMarking is needed.	 
	 if(lChLevel) {
	 	MarkVertLevel(pM2D,phid);
	 	lChLevel=false;
	 	}

//Here the remeshing criteria is applied	 
	 lrmsh=false;
	 for(iE=0;iE<pM2D->iNumFaces();iE++){
	 	pE=pM2D->pFFace(iE);
	 	VL=pM2D->iVertIndex(pE->pVVert(0));
	 	VR=pM2D->iVertIndex(pE->pVVert(1));
//	 	if(*(phid+VL)*(*(phid+VR))<= pow(*h1,4)/4){ // interface edges
	 	if(*(phid+VL)*(*(phid+VR))<= 0.0){ // interface edges
	  		if(pM2D->pVVert(VL)->dLS()>*h1*1.3||pM2D->pVVert(VR)->dLS()>*h1*1.3){
	  			if(debug)printf("==== remesh because of length============\n");
	  			lrmsh=true;
	  			break;
	  		}	
	  		 if(pM2D->pVVert(VL)->LSlevel>2||pM2D->pVVert(VR)->LSlevel>2) {
	  		 	if(debug)printf("==== remesh because of region************\n");
	  			lrmsh=true;
	  			break;
	  		}
	  	}
	  }
	 
// if remesh is needed then, level need to be remeshed. 
	 if(lrmsh) {
	 	MarkVertLevel(pM2D,phid);
	 	lChLevel=false;
	 	}

	 

	 int iLayer1_all=1;
	 int iLayer1_not=0;
	 double dLength;
	 if(debug&&lrmsh){
	 	sprintf(strBaseFileName,"LS_ADP_%4.4d.txt",iREMESH);
	 	fp=fopen(strBaseFileName,"w");
	 }
	 if(lrmsh)for(iV=0;iV<pM2D->iNumVerts();iV++){
	 	pV=pM2D->pVVert(iV);
//	 	double hlength=h_length(pV->dX(),pV->dY(),h1,h2);
	 	switch(pV->LSlevel){
	 	case 1:  // interface nodes, first layer
	 		iLayer1_all+=1;
	 		if(pV->dLS()>*h1*2) {
	 			dLength=pV->dLS()/2.0;
	 			pV->vSetLS(dLength);
	 			iLayer1_not+=1;
	 		}
	 		else	{
	 			dLength=*h1*1.0;
	 			pV->vSetLS(dLength);
	 		}
	 			
	 		pM2D->vMarkForLengthScaleCheck(pV);
	 		if(debug)fprintf(fp,"%d %10.8lf\n",iV,dLength);
	 		break;
	 	case -1: // boundary nodes
	 		dLength=MAX(*h2,*h3);
	 		pV->vSetLS(dLength);
	 		pM2D->vMarkForLengthScaleCheck(pV);
	 		if(debug)fprintf(fp,"%d %10.8lf\n",iV,dLength);
	 		break;
/*	 	case 2:  // second layer of interface
	 		if(pV->dLS()>*h1*2) {
	 			dLength= pV->dLS()/2.0;
	 			pV->vSetLS(dLength);
	 		}
	 		else	{
	 			dLength = *h1*1.0;	
	 			pV->vSetLS(dLength);
	 		}
	 		pM2D->vMarkForLengthScaleCheck(pV);
	 		if(debug)fprintf(fp,"%d %10.8lf\n",iV,dLength);
	 		break;
	 	case 3:  // third layer of interface
	 		if(pV->dLS()>*h1*2) {
	 			dLength= pV->dLS()/2.0;
	 			pV->vSetLS(dLength);
	 		}
	 		else	{
	 			dLength = *h1*1.0;	
	 			pV->vSetLS(dLength);
	 		}
	 		pM2D->vMarkForLengthScaleCheck(pV);
	 		if(debug)fprintf(fp,"%d %10.8lf\n",iV,dLength);
	 		break;
*/	 	case 10: // interior nodes. phi = +1
//	 		pV->vSetLS(sqrt(*h2*pV->dLS()));
			dLength = *h2*pow(pV->dLS()/(*h2),0.5);
			pV->vSetLS(dLength);
	 		pM2D->vMarkForLengthScaleCheck(pV);
	 		if(debug)fprintf(fp,"%d %10.8lf\n",iV,dLength);
	 		break;
	 	case 11: // interior nodes. phi = -1
//	 		pV->vSetLS(sqrt(*h2*pV->dLS()));
			dLength = *h3*pow(pV->dLS()/(*h3),0.5);
			pV->vSetLS(dLength);
	 		pM2D->vMarkForLengthScaleCheck(pV);
	 		if(debug)fprintf(fp,"%d %10.8lf\n",iV,dLength);
	 		break;
	 	default:
			dLength = *h2*pow(pV->dLS()/(*h2),0.5);
			pV->vSetLS(dLength);
	 		pM2D->vMarkForLengthScaleCheck(pV);
	 		if(debug)fprintf(fp,"%d %10.8lf\n",iV,dLength);
	 		break;

	 	}
	 	
	 }
	 
	 if(debug&&lrmsh)fclose(fp);
	 

	 if(debug)printf("The percentage of nodes in layer 1 without a good LS is %f\n",
	 				float(iLayer1_not)/iLayer1_all);			

// after set the interface

// try to get the boundary length scale to be fixed.
/*	  for(int iBE=0;iBE<pM2D->iNumBdryFaces();iBE++){
	  	pV=pM2D->pVVert(pM2D->iVertIndex(pM2D->pBFBFace(iBE)->pVVert(0)));
	  	pV->vSetLS(*h2);
	  	pM2D->vMarkForLengthScaleCheck(pV);
	  	pV=pM2D->pVVert(pM2D->iVertIndex(pM2D->pBFBFace(iBE)->pVVert(1)));
	  	pV->vSetLS(*h2);
	  	pM2D->vMarkForLengthScaleCheck(pV);
	  	}
	 

	 	sprintf(strBaseFileName,"LS_beforeUpdate_%3.3d.plt",iRUNS);
	 	FILE *fp=fopen(strBaseFileName,"w");
	 	fprintf(fp,"zone f=fepoint,N=%d,E=%d,ET=TRIANGLE\n",pM2D->iNumVerts(),pM2D->iNumCells());
	 	for(iV=0;iV<pM2D->iNumVerts();iV++)
	 		fprintf(fp,"%f %f %f\n",pM2D->pVVert(iV)->dX(),pM2D->pVVert(iV)->dY(),pM2D->pVVert(iV)->dLS());
	 	for(iCell=0;iCell<pM2D->iNumCells();iCell++){
	 		Cell *pC=pM2D->pCCell(iCell);
	 		fprintf(fp,"%d %d %d\n",pM2D->iVertIndex(pC->pVVert(0))+1,
	 			pM2D->iVertIndex(pC->pVVert(1))+1,pM2D->iVertIndex(pC->pVVert(2))+1);
	 	}
	 	fclose(fp);	
	
*/	

// Adapting the mesh based on the length scale. 
	if(debug)printf("the value of lrmsh is %d with coarse = %f, fine = %f\n",lrmsh, *h2, *h1);
	if(lrmsh){
		if(debug)printf("before remeshing applied after the length scale setting\n");	 
	 	pM2D->vUpdateLengthScale();
	 	
		if(debug){
	 	    sprintf(strBaseFileName,"LS_afterUpdate%4.4d.txt",iREMESH);
	 	    fp=fopen(strBaseFileName,"w");
	 	    for(iV=0;iV<pM2D->iNumVerts();iV++)
	 	    	fprintf(fp,"%d %f\n",iV,pM2D->pVVert(iV)->dLS());
	 	    fclose(fp);
	 	}	
	 	if(debug)printf("update the length scale is ok\n");	
	 	
//	 	if(iREMESH==24)iMessageStdoutLevel=1; 
	 	pM2D->vAdaptToLengthScale(eET);
	 	lChLevel = true;
	 	pM2D->vSetSmoothingThreshold(30.);
	 	pM2D->iSmooth(1);
//	 	iMessageStdoutLevel=0;
	 	if(debug)printf("remeshing applied after the length scale setting\n");	 
	} // end of remeshing
	
	return lrmsh;

}
//=========================================================================
void GetNodes(int *nnode, double *coord)
{
	return;
}
double h_length(double x,double y,double *h1,double *h2)
{
	double h0,r0,eps,c,r,h;
	double hh,aa;
	r0=1.0;
	h0=*h2;
	c=1.0;
	eps=*h1;
	hh=(c*eps*h0)/(c*eps+h0)/10.0;
//	aa=(1.0+hh*double(iRUNS-1))*(1.0+hh*double(iRUNS-1));
	r=sqrt(x*x+y*y);
	h=pow(cosh((r-r0)/sqrt(2)/eps),2);

	h=1.0/(1/(h*c*eps)+1/h0);

	return h;
}

double phi(double x,double y,int iRUNS)
{
	double h0,r0,eps,c,r,h;
	double hh,aa;
	r0=1.0;
	h0=0.5;
	c=1.0;
	eps=0.02;
	hh=(c*eps*h0)/(c*eps+h0)/10.0;
	aa=(1.0+hh*double(iRUNS-1))*(1.0+hh*double(iRUNS-1));
	r=sqrt(x*x/aa+y*y);
	return r-r0;
}

double h_length_smallest()
{
	double h0,r0,eps,c,r,h;
	double hh,aa;
	r0=1.0;
	h0=0.5;
	c=1.0;
	eps=0.02;
	hh=(c*eps*h0)/(c*eps+h0);
	return hh;
}
bool thesameedge(int na1,int na2,int nb1,int nb2)
{
	if(na1==nb1&&na2==nb2) {
		return true;
	}
	else if(na1==nb2&&na2==nb1) {
		return true;
	}
	else {
		return false;
	}
}

void check_clockwise(double *x1,double *x2,double *x3,int iCell)
{
	double r1x,r1y,r2x,r2y;
	r1x=*x2-*x1;
	r1y=*(x2+1)-*(x1+1);
	r2x=*x3-*x1;
	r2y=*(x3+1)-*(x1+1);
	if(r1x*r2y-r1y*r2x<0.0)
		printf("**********Not counter clockwise for Cell %d\n",iCell);
//	else
//		printf("counter clockwise for Cell %d\n",iCell);	
		
	return;	
}

// this subroutine to Return the mesh information to the solver. 
void ReturnToSolver(Mesh2D *pM2D,int *nbd,int *nvert,int *nelem,int *ngmx,int *ncmx,
		int *refseg,int *bdseg,int *inod,int *nec,int *reft,double *coorx,double *coory)
{
	  bool debug = false;
	  double sum = 0;
	  for(int iV=0;iV<*nvert;iV++){
		*(coorx+iV)=(pM2D->pVVert(iV))->dX();	
        	*(coory+iV)=(pM2D->pVVert(iV))->dY();	
        	sum = sum + pM2D->pVVert(iV)->dLS();
        	//printf("======== length scale = %g \n", pM2D->pVVert(iV)->dLS());
        	//printf("======== Y= %g \n\n", *(coord+2*iV+1));
        }
        if(debug)printf("the average length scale is %f \n\n", sum/double(*nvert));
        if(debug)printf("\n======== obtain the verts ok! ==========\n\n");
        
        Face *pF,*pF1,*pF2,*pF3;
        Cell *pC;
        int iCA,iCB;
        for(int iCell=0;iCell<*nelem;iCell++) {
        	pC=pM2D->pCCell(iCell);
        	*(reft+iCell)=1;
      	*(inod+(*ngmx)*iCell)  =pM2D->iVertIndex(pC->pVVert(0))+1;
      	*(inod+(*ngmx)*iCell+1)=pM2D->iVertIndex(pC->pVVert(1))+1;
      	*(inod+(*ngmx)*iCell+2)=pM2D->iVertIndex(pC->pVVert(2))+1;
//      	check_clockwise(coor+(*(inod+3*iCell)-1)*2,
//      	   coor+(*(inod+3*iCell+1)-1)*2,coor+(*(inod+3*iCell+2)-1)*2,iCell+1);
  
 // begin to getting the neighboring cells. 
 
 // first edge..     	
      	pF1=pM2D->pCCell(iCell)->pFFace(0);
      	pF2=pM2D->pCCell(iCell)->pFFace(1);
      	pF3=pM2D->pCCell(iCell)->pFFace(2);
      	if(thesameedge(pM2D->iVertIndex(pF1->pVVert(0))+1,
      	   pM2D->iVertIndex(pF1->pVVert(1))+1,*(inod+(*ngmx)*iCell),*(inod+(*ngmx)*iCell+1))){
      	   	pF=pF1;
      	}
      	else if(thesameedge(pM2D->iVertIndex(pF2->pVVert(0))+1,
      	   pM2D->iVertIndex(pF2->pVVert(1))+1,*(inod+(*ngmx)*iCell),*(inod+(*ngmx)*iCell+1))){
      	   	pF=pF2;
//      	   	if(debug)printf("not counter-clockwise order: cell = %d\n",iCell+1);
      	}
      	else {
      		pF=pF3;
      	}
      	
      	
      	pC=pF->pCCellLeft();
     	      switch (pC->eType()) {
             case Cell::eBdryEdge:
             case Cell::eTriBFace:
             case Cell::eQuadBFace:
             case Cell::eIntBdryEdge:
             case Cell::eIntTriBFace:
             case Cell::eIntQuadBFace:
               iCA = - ((BFace*)pC)->iBdryCond()-1;
              break;
 
             default:
               iCA = pM2D->iCellIndex(pC);
             break;
            }

      	pC=pF->pCCellRight();
     	      switch (pC->eType()) {
             case Cell::eBdryEdge:
             case Cell::eTriBFace:
             case Cell::eQuadBFace:
             case Cell::eIntBdryEdge:
             case Cell::eIntTriBFace:
             case Cell::eIntQuadBFace:
               iCB = - ((BFace*)pC)->iBdryCond()-1;
              break;
 
             default:
               iCB = pM2D->iCellIndex(pC);
             break;
            }

      	
      	*(nec+(*ncmx)*iCell)=iCA+1;
      	if(iCA==iCell) *(nec+(*ncmx)*iCell)=iCB+1;
      	if(iCell==0&&debug) printf("%d %d\n",iCA,iCB);
      	
      	
 
 //second edge..     	

      	if(thesameedge(pM2D->iVertIndex(pF2->pVVert(0))+1,
      	   pM2D->iVertIndex(pF2->pVVert(1))+1,*(inod+(*ngmx)*iCell+1),*(inod+(*ngmx)*iCell+2))){
      	   	pF=pF2;
      	}
      	else if(thesameedge(pM2D->iVertIndex(pF1->pVVert(0))+1,
      	   pM2D->iVertIndex(pF1->pVVert(1))+1,*(inod+(*ngmx)*iCell+1),*(inod+(*ngmx)*iCell+2))){
      	   	pF=pF1;
//      	   	if(debug)printf("not counter-clockwise order: cell = %d\n",iCell+1);
      	}
      	else {
      		pF=pF3;
      	}

      	pC=pF->pCCellLeft();
     	      switch (pC->eType()) {
             case Cell::eBdryEdge:
             case Cell::eTriBFace:
             case Cell::eQuadBFace:
             case Cell::eIntBdryEdge:
             case Cell::eIntTriBFace:
             case Cell::eIntQuadBFace:
               iCA = - ((BFace*)pC)->iBdryCond()-1;
              break;
 
             default:
               iCA = pM2D->iCellIndex(pC);
             break;
            }

      	pC=pF->pCCellRight();
     	      switch (pC->eType()) {
             case Cell::eBdryEdge:
             case Cell::eTriBFace:
             case Cell::eQuadBFace:
             case Cell::eIntBdryEdge:
             case Cell::eIntTriBFace:
             case Cell::eIntQuadBFace:
               iCB = - ((BFace*)pC)->iBdryCond()-1;
              break;
 
             default:
               iCB = pM2D->iCellIndex(pC);
             break;
            }

      	
      	*(nec+(*ncmx)*iCell+1)=iCA+1;
      	if(iCA==iCell) *(nec+(*ncmx)*iCell+1)=iCB+1;
      	if(iCell==0&&debug) printf("%d %d\n",iCA,iCB);
        	
 
        	
 //third edge..     	
      	if(thesameedge(pM2D->iVertIndex(pF3->pVVert(0))+1,
      	   pM2D->iVertIndex(pF3->pVVert(1))+1,*(inod+(*ngmx)*iCell+2),*(inod+(*ngmx)*iCell))){
      	   	pF=pF3;
      	}
      	else if(thesameedge(pM2D->iVertIndex(pF2->pVVert(0))+1,
      	   pM2D->iVertIndex(pF2->pVVert(1))+1,*(inod+(*ngmx)*iCell+2),*(inod+(*ngmx)*iCell))){
      	   	pF=pF2;
//      	   	if(debug)printf("not counter-clockwise order: cell = %d\n",iCell+1);
      	}
      	else {
      		pF=pF1;
      	}

      	pC=pF->pCCellLeft();
     	      switch (pC->eType()) {
             case Cell::eBdryEdge:
             case Cell::eTriBFace:
             case Cell::eQuadBFace:
             case Cell::eIntBdryEdge:
             case Cell::eIntTriBFace:
             case Cell::eIntQuadBFace:
               iCA = - ((BFace*)pC)->iBdryCond()-1;
              break;
 
             default:
               iCA = pM2D->iCellIndex(pC);
             break;
            }

      	pC=pF->pCCellRight();
     	      switch (pC->eType()) {
             case Cell::eBdryEdge:
             case Cell::eTriBFace:
             case Cell::eQuadBFace:
             case Cell::eIntBdryEdge:
             case Cell::eIntTriBFace:
             case Cell::eIntQuadBFace:
               iCB = - ((BFace*)pC)->iBdryCond()-1;
              break;
 
             default:
               iCB = pM2D->iCellIndex(pC);
             break;
            }

      	
      	*(nec+(*ncmx)*iCell+2)=iCA+1;
      	if(iCA==iCell) *(nec+(*ncmx)*iCell+2)=iCB+1;
      	if(iCell==0&&debug) printf("%d %d\n",iCA,iCB);
        
        }
        
        if(debug)printf("\n======== obtain the elements ok!========\n\n\n");
        
                
        for(int iBFace=0;iBFace<*nbd;iBFace++){
        	*(bdseg+2*iBFace) = pM2D->iVertIndex(pM2D->pBFBFace(iBFace)->pVVert(0))+1;
        	*(bdseg+2*iBFace+1) = pM2D->iVertIndex(pM2D->pBFBFace(iBFace)->pVVert(1))+1;
        	*(refseg+iBFace) = pM2D->pBFBFace(iBFace)->iBdryCond();
        }
	
}


void MarkVertLevel(Mesh2D *pM2D,double *phid)
{
	int iV,iE,iBE,VL,VR;
	Face *pE;
	
	
	 
	   for(iV=0;iV<pM2D->iNumVerts();iV++){
	   	if(*(phid+iV)>0.0)
	   		pM2D->pVVert(iV)->LSlevel=10;
	   	else 
	   		pM2D->pVVert(iV)->LSlevel=11;
	   	}	


	   for(int iBE=0;iBE<pM2D->iNumBdryFaces();iBE++){
	  	pM2D->pVVert(pM2D->iVertIndex(pM2D->pBFBFace(iBE)->pVVert(0)))->LSlevel=-1;
	  	pM2D->pVVert(pM2D->iVertIndex(pM2D->pBFBFace(iBE)->pVVert(1)))->LSlevel=-1;
	   } 
	 
	   for(iE=0;iE<pM2D->iNumFaces();iE++){
	 	pE=pM2D->pFFace(iE);
	 	VL=pM2D->iVertIndex(pE->pVVert(0));
	 	VR=pM2D->iVertIndex(pE->pVVert(1));
//	 	if(*(phid+VL)*(*(phid+VR))<= pow(*h1,4)/4){ // interface edges
	 	if(*(phid+VL)*(*(phid+VR))<= 0.0){ // interface edges
	  		pM2D->pVVert(VL)->LSlevel=1;
	  		pM2D->pVVert(VR)->LSlevel=1;	
	  	}
	  }
	  
	 
	  
	 
	   for(iE=0;iE<pM2D->iNumFaces();iE++){
	 	pE=pM2D->pFFace(iE);
	 	VL=pM2D->iVertIndex(pE->pVVert(0));
	 	VR=pM2D->iVertIndex(pE->pVVert(1)); // near interface edges
	 	if(pM2D->pVVert(VL)->LSlevel==1 && (pM2D->pVVert(VR)->LSlevel==10
	 	  ||pM2D->pVVert(VR)->LSlevel==11||pM2D->pVVert(VR)->LSlevel==-1))
	  					pM2D->pVVert(VR)->LSlevel=2;
	  	if(pM2D->pVVert(VR)->LSlevel==1 && (pM2D->pVVert(VL)->LSlevel==10
	  	  ||pM2D->pVVert(VL)->LSlevel==11||pM2D->pVVert(VL)->LSlevel==-1))
	  					pM2D->pVVert(VL)->LSlevel=2;
	   }
  	 
  	   for(iE=0;iE<pM2D->iNumFaces();iE++){
	 	pE=pM2D->pFFace(iE);
	 	VL=pM2D->iVertIndex(pE->pVVert(0));
	 	VR=pM2D->iVertIndex(pE->pVVert(1)); // near interface edges
	 	if(pM2D->pVVert(VL)->LSlevel==2 && (pM2D->pVVert(VR)->LSlevel==10
	 	 ||pM2D->pVVert(VR)->LSlevel==11||pM2D->pVVert(VR)->LSlevel==-1))
	  					pM2D->pVVert(VR)->LSlevel=3;
	  	if(pM2D->pVVert(VR)->LSlevel==2 && (pM2D->pVVert(VL)->LSlevel==10 
	  	 ||pM2D->pVVert(VL)->LSlevel==11||pM2D->pVVert(VL)->LSlevel==-1))
	  					pM2D->pVVert(VL)->LSlevel=3;
	  }
	  
	return;
}

