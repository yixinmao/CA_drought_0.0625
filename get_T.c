#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159
#define RADIUS 6371.228 //radius of earth, km
#define MAXY 200
#define MAXC 500

double get_dist(double, double, double, double);
double calc_area(float lat, float lon, float cellsize);
int get_length(FILE *file);

int main (int argc, char *argv[])   // average T (Nov 1 - Mar 19) for each year [deg C]
{
  int i,j,nfiles,ndays;
  int syear,year,month,day,y;
  double tot_area;
  double lat,lon,cellsize,area;
  double prec,evap,runoff,baseflow;
  double airT,sm1,sm2,sm3,swe;
  double T_year[MAXY];
  int nday_year[MAXY];   // nday_year: count the days in Nov1-Mar19 in each water year
  char filename[MAXC];
  FILE *fpin,*fplist,*fpout;

  if(argc!=4){
    printf("Usage: %s latlonlist cellsize outfile\n",argv[0]);
    exit(0);
  }
  if((fplist = fopen(argv[1],"r"))==NULL){
    printf("ERROR: can't open %s\n", argv[1]);
    exit(0);
  }
  nfiles=get_length(fplist);
  cellsize=atof(argv[2]);
  for(j=0;j<MAXY;j++){
    T_year[j]=0;
  }
  tot_area=0;
  for(i=0;i<nfiles;i++){
    fscanf(fplist,"%lf %lf",&lat,&lon);
    area=calc_area(lat,lon,cellsize);
    tot_area+=area;
    sprintf(filename,"./VIC_forecast_for_CA/whole_period/fluxes_%.5f_%.5f",lat,lon);
    if((fpin = fopen(filename,"r"))==NULL){
      printf("ERROR: can't open %s\n", filename);
      exit(0);
    }
    ndays=get_length(fpin);
	  for(j=0;j<MAXY;j++){
	    nday_year[j]=0;
	  }

    for(j=0;j<ndays;j++){
      fscanf(fpin,"%d %d %d",&year,&month,&day);
      if(j==0)syear=year;
      fscanf(fpin,"%lf %lf %lf %lf",&prec,&evap,&runoff,&baseflow);
      fscanf(fpin,"%lf %lf %lf %lf %lf",&airT,&sm1,&sm2,&sm3,&swe);
      if(month>=11)   // if Nov or Dec, add it to the next water year
	{
	T_year[year-syear+1]+=airT*area;
	nday_year[year-syear+1]++;
      }
	else if(month<=3 || (month==4 && day<=1))  // if Jan - Apr1, add it to this water year
	{
		T_year[year-syear]+=airT*area;
		nday_year[year-syear]++;
	}
    }
    fclose(fpin);
  }
  fclose(fplist);

	for(y=0;y<=year-syear;y++)
	{
		T_year[y] = T_year[y] / nday_year[y] / tot_area;
	}

  if((fpout = fopen(argv[3],"w"))==NULL){
    printf("ERROR: can't open %s\n", argv[3]);
    exit(0);
  }
  for(y=syear+1;y<=year;y++){    // only print water years 1917-2014
    fprintf(fpout,"%d %.4f\n",y,T_year[y-syear]);
  }
  fclose(fpout);
  printf("%.1f\n",tot_area);

  return(0);
}
/*******************************************************************
   Function:
   double distance(double lat1, double long1, double lat2, double long2)
   Returns : distance between two locations
********************************************************************/

double get_dist(double lat1, double long1, double lat2, double long2)
{
  double theta1;
  double phi1;
  double theta2;
  double phi2;
  double dtor;
  double term1;
  double term2;
  double term3;
  double temp;
  double distance;

  /*Convert from degrees to radians*/
  dtor = 2.0*PI/360.0;
  theta1 = dtor*long1;
  phi1 = dtor*lat1;
  theta2 = dtor*long2;
  phi2 = dtor*lat2;

  term1 = cos(phi1)*cos(theta1)*cos(phi2)*cos(theta2);
  term2 = cos(phi1)*sin(theta1)*cos(phi2)*sin(theta2);
  term3 = sin(phi1)*sin(phi2);
  temp = term1+term2+term3;
  temp = (double) (1.0 < temp) ? 1.0 : temp;
  distance = RADIUS*acos(temp);

  return distance;
}

double calc_area(float lat, float lon, float cellsize)
{
  double cell_area;
  float start_lat,right_lng,left_lng;
  double delta;
  int d;
  double dist;

  lat = fabs(lat);
  lon = fabs(lon);

  start_lat = lat - cellsize / 2.;
  right_lng = lon + cellsize / 2;
  left_lng = lon - cellsize / 2;

  delta = get_dist(lat,lon,lat+cellsize/10.,lon);
  dist = 0.;
  for ( d = 0; d < 10; d++ ) {
    dist += get_dist(start_lat,left_lng,start_lat,right_lng) * delta;
    start_lat += cellsize/10;
  }

  cell_area = dist;
  return(cell_area);
}
int get_length(FILE *file)
{
  /* Calculates number of lines in an ascii file from Dave Sinkula on
     http://www.daniweb.com/code/snippet325.htm*/
  int ch, prev = '\n' /* so empty files have no lines */, lines = 0;

  while ( (ch = fgetc(file)) != EOF ) /* Read all chars in the file. */
    {
      if ( ch == '\n' )
        {
          ++lines; /* Bump the counter for every newline. */
        }
      prev = ch; /* Keep a copy to later test whether... */
    }
  if ( prev != '\n' ) /* ...the last line did not end in a newline. */
    {
      ++lines; /* If so, add one more to the total. */
    }
  rewind(file);
  return(lines);
}
