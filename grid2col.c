/* this program takes an arcinfo grid and rewrites it into
a column that can be read and plotted by gmt Jenny 7/2001   */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
  
#define MAX 20

int main (int argc, char *argv[])
{
  int i,j;
  char str;
  float lat,lon;
  float data;
  int nodata;
  int nrows,ncols;
  float cellsize;
  float xllcorner,yllcorner;
 
  

  FILE *fin,*fout;
  
  /* usage */
  if(argc!=3) {
    printf("Usage: %s  in_grid  out_file \n ",argv[0]);
    exit(0);
  }

  
  /* open the input file */
  if((fin=fopen(argv[1],"r"))==NULL)   {
    printf("ERROR: Unable to open %s\n",argv[1]);
    exit(0);
  }


  /* open the output file */
  if((fout=fopen(argv[2],"w"))==NULL)   {
    printf("ERROR: Unable to open %s\n",argv[2]);
    exit(0);
  }
  
  /* read in the header */
  fscanf(fin,"%s",&str);
  fscanf(fin,"%d",&ncols);
  fprintf(stderr,"ncols is %d\n",ncols);
  fscanf(fin,"%s",&str);
  fscanf(fin,"%d",&nrows);
  fprintf(stderr,"nrows is %d\n",nrows);
  fscanf(fin,"%s",&str);
  fscanf(fin,"%f",&xllcorner);
  fprintf(stderr,"xllcorner is %f\n",xllcorner);
  fscanf(fin,"%s",&str);
  fscanf(fin,"%f",&yllcorner);
  fprintf(stderr,"yllcorner is %f\n",yllcorner);
  fscanf(fin,"%s",&str);
  fscanf(fin,"%f",&cellsize);
  fprintf(stderr,"cellsize is %f\n",cellsize);
  fscanf(fin,"%s",&str);
  fscanf(fin,"%d",&nodata);
  fprintf(stderr,"nodata is %d\n",nodata);

  /* initialize lat and long */
  lon = xllcorner + cellsize/2;
  lat = yllcorner + cellsize*nrows - cellsize/2;

  /* read in the data */
  
  for(i=0;i<nrows;i++) {
    for(j=0;j<ncols;j++) {
      fscanf(fin,"%f",&data);
      fprintf(fout,"%.4f ",lat);
      fprintf(fout,"%.4f ",lon);
      fprintf(fout,"%.2f",data);
      fprintf(fout,"\n");
      lon = lon + cellsize;
    }
    lat = lat - cellsize;
    lon = xllcorner + cellsize/2;
  }
  

  fclose(fout);  
  fclose(fin); 
  return 0;
  
}

