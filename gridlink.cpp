#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "cellarray.h"


//barebones gridlink 

// #################################
// Get Bin Size, i.e size of each bin
double get_binsize(
          const double xmin, const double xmax, const double rmax, 
          const int refine_factor, 
          const int max_ncells, 
          int *nlattice)
{


  double xdiff = xmax-xmin;

  //nmesh is number of bins. rmax is max size of each bin 
  int nmesh=(int)(refine_factor*xdiff/rmax) ;
  
  if (nmesh>max_ncells)  nmesh=max_ncells;

  double xbinsize = xdiff/nmesh;
  *nlattice = nmesh;
  return xbinsize;
}
// #################################




// #################################
// create bins given range of x, y,z and the max bin size (rmax)
//why triple pointer????
cellarray *** gridlink(unsigned int np,
		       DOUBLE *x,DOUBLE *y,DOUBLE *z,
		       double xmin, double xmax,
		       double ymin, double ymax,
		       double zmin, double zmax,
		       double rmax,
		       int *nlattice_x,
		       int *nlattice_y,
		       int *nlattice_z)
{
  cellarray ***lattice= NULL;

  int ix,iy,iz;
  int nmesh_x,nmesh_y,nmesh_z; //number of bins in x, y, and z dimensions
  int ***nallocated=NULL;
  int index;
  double xdiff,ydiff,zdiff;
  double cell_volume,box_volume;
  double xbinsize,ybinsize,zbinsize;
  int expected_n=0;
  struct timeval t0,t1;//??
  gettimeofday(&t0,NULL); //??
  
  //get bin size, defined above
  xbinsize = get_binsize(xmin,xmax,rmax,BIN_REFINE_FACTOR, NLATMAX, &nmesh_x);
  ybinsize = get_binsize(ymin,ymax,rmax,BIN_REFINE_FACTOR, NLATMAX, &nmesh_y);
  zbinsize = get_binsize(zmin,zmax,rmax,BIN_REFINE_FACTOR, NLATMAX, &nmesh_z);
  
  //get range of x, y, and z dimensions
  xdiff = xmax-xmin;
  ydiff = ymax-ymin;
  zdiff = zmax-zmin;
  
  cell_volume = xbinsize*ybinsize*zbinsize;   //volume of each bin
  box_volume=xdiff*ydiff*zdiff;    //volume of entire map
  expected_n=(int)(np*cell_volume/box_volume*MEMORY_INCREASE_FAC); //??

  fprintf(stderr,"In gridlink> Running with (nmesh = number of bins) [nmesh_x, nmesh_y, nmesh_z]  = %d,%d,%d. ",nmesh_x,nmesh_y,nmesh_z);

  //allocate memory for lattice (lattice of bins); memory of each bin = sizeof(cellarray)
  lattice = (cellarray ***) volume_malloc(sizeof(cellarray),nmesh_x,nmesh_y,nmesh_z); 

  //keep track of bins to which memory has been allocated
  nallocated = (int ***) volume_malloc(sizeof(unsigned int),nmesh_x,nmesh_y,nmesh_z);

  //allocate memory to bins *important*
  for (int i=0;i<nmesh_x;i++) {
    for (int j=0;j<nmesh_y;j++) {
      for (int k=0;k<nmesh_z;k++) {
      	lattice[i][j][k].x = my_malloc(sizeof(DOUBLE),expected_n);
      	lattice[i][j][k].y = my_malloc(sizeof(DOUBLE),expected_n);
      	lattice[i][j][k].z = my_malloc(sizeof(DOUBLE),expected_n);
      	nallocated[i][j][k] = expected_n;
      	lattice[i][j][k].nelements=0;
      }
    }
  }

  //replace all mallocing above with c plus plus


  double xinv=1.0/xbinsize;
  double yinv=1.0/ybinsize;
  double zinv=1.0/zbinsize;
  

// ####################  what is purpose of this part?
  //what is np, x, y, z? 
  for (unsigned int i=0;i<np;i++)  {
    ix=(int)((x[i]-xmin)*xinv) ; //number of xbins; ix = (x coordinate - xmin)/xbinsize
    iy=(int)((y[i]-ymin)*yinv) ;
    iz=(int)((z[i]-zmin)*zinv) ;
    if (ix>nmesh_x-1)  ix=nmesh_x-1 ;  //  make sure ix is in range of nmesh_x. what is purpose of ix?
    if (iy>nmesh_y-1)  iy=nmesh_y-1 ;
    if (iz>nmesh_z-1)  iz=nmesh_z-1 ;
    assert(ix >= 0 && ix < nmesh_x && "ix is in range");
    assert(iy >= 0 && iy < nmesh_y && "iy is in range");
    assert(iz >= 0 && iz < nmesh_z && "iz is in range");
    
    //what is expected_n?
    if(lattice[ix][iy][iz].nelements == nallocated[ix][iy][iz]) {
      expected_n = nallocated[ix][iy][iz]*MEMORY_INCREASE_FAC;

      if(expected_n == nallocated[ix][iy][iz])
	       expected_n += 3;

        lattice[ix][iy][iz].x = my_realloc(lattice[ix][iy][iz].x ,sizeof(DOUBLE),expected_n,"lattice.x");
        lattice[ix][iy][iz].y = my_realloc(lattice[ix][iy][iz].y ,sizeof(DOUBLE),expected_n,"lattice.y");
        lattice[ix][iy][iz].z = my_realloc(lattice[ix][iy][iz].z ,sizeof(DOUBLE),expected_n,"lattice.z");
        
        nallocated[ix][iy][iz] = expected_n;
    }

    assert(lattice[ix][iy][iz].nelements < nallocated[ix][iy][iz] && "Ensuring that number of particles in a cell doesn't corrupt memory");
    index=lattice[ix][iy][iz].nelements;
    lattice[ix][iy][iz].x[index] = x[i];
    lattice[ix][iy][iz].y[index] = y[i];
    lattice[ix][iy][iz].z[index] = z[i];
    lattice[ix][iy][iz].nelements++;
  }
// #################################

  //free nallocated
  volume_free((void ***) nallocated,nmesh_x,nmesh_y);
  
  //put number of bins in x,y,z dimension into memory
  *nlattice_x=nmesh_x;
  *nlattice_y=nmesh_y;
  *nlattice_z=nmesh_z;

  gettimeofday(&t1,NULL);
  fprintf(stderr," Time taken = %6.2lf sec\n",ADD_DIFF_TIME(t0,t1));
  
  return lattice;
}

