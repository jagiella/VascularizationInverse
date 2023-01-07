/*
 * extractROI.c
 *
 *  Created on: Oct 15, 2012
 *      Author: jagiella
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <unistd.h>

// OUTPUT OPTIONS
#define ASCII	1
#define BINARY	2

// ROI TYPE
//#define allPossibleRadii

#define MIN(a,b) (a<b?a:b)

//enum ROI_Type { Circlular, Horizontal, Vertical };

bool isInROI_Circular( int &x, int &y, int &mx, int &my, int r){
	if( r<0) return false;
	return ( pow(x-mx,2) + pow(y-my,2) <= pow(r,2) );
}

bool isInROI_Horizontal( int &x, int &y, int &mx, int &my, int r){
	return ( abs(y-my) <= r );
}

bool isInROI_Vertical( int &x, int &y, int &mx, int &my, int r){
	return ( abs(x-mx) <= r );
}

int shortestPathRadius(int x, int y, int z)
{ }

int main( int argc, char **argv)
{
	// ROI
	int mx=0, my=0, mz=0, r=0, d=1;
	int minR = 0;

	char output = ASCII;
	bool allPossibleRadii = false;
	bool outerRim = false;

	//int ROItype = Circlular;
	bool (*isInROI) ( int &x, int &y, int &mx, int &my, int r) = isInROI_Circular;


	char c;
	while ((c = getopt(argc, argv, "x:y:r:d:RSBV")) != -1)
		switch (c) {
		case 'x':
			mx = atoi(optarg);
			break;
		case 'y':
			my = atoi(optarg);
			break;
		case 'r':
			if(r!=0)
				minR=r;
			r = atoi(optarg);
			break;
		case 'd':
			d = atoi(optarg);
			break;
		case 'R':
			allPossibleRadii = true;
			break;
		case 'S':
			outerRim = true;
			break;
		case 'B':
			output = BINARY;
			break;
		case 'V':
			//ROItype = Vertical;
			isInROI = isInROI_Vertical;
			break;
		//case 'C':

		}


	// EXPERIMENTAL
	int maxR = MIN(mx, my);
	if( r!=0)
		maxR = r;

	for ( int index = optind; index < argc; index++){
		if( allPossibleRadii){
			float	avgValue[maxR];
			int		avgCount[maxR];

			FILE *fp = fopen( argv[index], "rb");
			if( fp){
				int X,Y,Z;
				fread( &X, sizeof(int), 1, fp);
				fread( &Y, sizeof(int), 1, fp);
				fread( &Z, sizeof(int), 1, fp);

				fprintf( stderr, "read file %s (data dimensions: %ix%ix%i)\n",
						argv[index], X, Y, Z);

				if( index == optind ){
					maxR = MIN(maxR, X-mx);
					maxR = MIN(maxR, Y-my);
				}
				fprintf( stderr, "Analyze ROI's around (%i,%i) with radii from 0 to %i\n", mx, my, maxR-1);

				for( r=minR; r<maxR; r++ ){
					avgValue[r]=0;
					avgCount[r]=0;
				}
				for( int x=0; x<X; x++)
					for( int y=0; y<Y; y++)
						for( int z=0; z<Z; z++){
							float value;
							if( 1 != fread( &value, sizeof(float), 1, fp)){
								perror("fread");
							}
							for( r=minR; r<maxR; r++ )
//							if( pow(x-mx,2) + pow(y-my,2) + pow(z-mz,2) <= pow(r,2)){
//								if(!outerRim || pow(x-mx,2) + pow(y-my,2) + pow(z-mz,2) >= pow(r-1,2)){
							if( (*isInROI)(x,y,mx,my,r)){
								if(!outerRim || !(*isInROI)(x,y,mx,my,r-d)){
									avgValue[r]+=value;
									avgCount[r]++;
								}
							}
						}

				fclose(fp);
				fprintf( stderr, "Average over %i voxels\n", avgCount[minR]);

				char filename[1024];
				if(outerRim)
					sprintf( filename, "%s.ROIRIM%ix%ir%i-%i", argv[index], mx,my, minR, maxR-1);
				else
					sprintf( filename, "%s.ROI%ix%ir%i-%i", argv[index], mx,my, minR, maxR-1);

				fp = fopen( filename, "wb+");
				if( fp==0){
					fprintf(stderr, "Can not open file %s\n", filename);
					perror("fopen()");
				}
				X=maxR-minR;
				Y=Z=1;
				fwrite( &X,  sizeof(int), 1, fp);
				fwrite( &Y,  sizeof(int), 1, fp);
				fwrite( &Z,  sizeof(int), 1, fp);

				for( r=minR; r<maxR; r++ ){
					avgValue[r]/=avgCount[r];
					fwrite( &avgValue[r],  sizeof(float), 1, fp);
				}
				fclose(fp);

				if( index == optind ){
					FILE *fp_ROI = fopen( (outerRim?"ROIRIMsize.dat":"ROIsize.dat"), "w+");
					for( r=minR; r<maxR; r++ )
						fprintf( fp_ROI, "%i %i\n", r, avgCount[r]);
					fclose(fp_ROI);
				}
			}

		}else{
			float	avgValue = 0;
			int		avgCount = 0;

			FILE *fp = fopen( argv[index], "rb");
			if( fp){
				int X,Y,Z;
				fread( &X, sizeof(int), 1, fp);
				fread( &Y, sizeof(int), 1, fp);
				fread( &Z, sizeof(int), 1, fp);

				fprintf( stderr, "read file %s (data dimensions: %ix%ix%i)\n",
						argv[index], X, Y, Z);


				for( int x=0; x<X; x++)
					for( int y=0; y<Y; y++)
						for( int z=0; z<Z; z++){
							float value;
							if( 1 != fread( &value, sizeof(float), 1, fp)){
								perror("fread");
							}
							if( pow(x-mx,2) + pow(y-my,2) + pow(z-mz,2) <= pow(r,2)){
								avgValue+=value;
								avgCount++;
							}
						}

				fclose(fp);
				fprintf( stderr, "avg: %f (of %i voxels)\n", avgValue/avgCount, avgCount);

				char filename[1024];
				sprintf( filename, "%s.ROI%ix%ir%i", argv[index], mx,my,r);
				fp = fopen( filename, "wb+");
				if( fp==0){
					fprintf(stderr, "Can not open file %s\n", filename);
					perror("fopen()");
				}
				X=Y=Z=1;
				fwrite( &X,  sizeof(int), 1, fp);
				fwrite( &Y,  sizeof(int), 1, fp);
				fwrite( &Z,  sizeof(int), 1, fp);
				avgValue/=avgCount;
				fwrite( &avgValue,  sizeof(float), 1, fp);
				fclose(fp);
			}
		}
	}

	return 0;
}
