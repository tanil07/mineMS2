#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;

//These function were originally present in xcms, and are recoded in C++.
//They are recoded in C++ for commodity purpose.

struct idxStruct
{
	int  from;
	int  to;
};


int lowerBound(double val,NumericVector mzval,int first, int length){
	int half,mid;
	while (length > 0) {
		half = length >> 1;
		mid = first;
		mid += half;
		if ( mzval[mid] < val){
			first = mid;
			first ++;
			length = length - half -1;
		}
		else length = half;
	}
	return(first);
}

int upperBound(double val,NumericVector mzval,int first, int length){
	int half,mid;
	while (length > 0) {
		half = length >> 1;
		mid = first;
		mid += half;
		if (val < mzval[mid]){
			length = half;
		}
		else {
			first = mid;
			first ++;
			length = length - half -1;
		}
	}
	return(first);
}


//Find the closest match.
// [[Rcpp::export]]
IntegerVector closeMatch(NumericVector x, NumericVector y, IntegerVector xidx,
                            IntegerVector yidx, double ppm, double dmz) {
	int lastlb=0;
	int xi,yi,lb,ub,txi,from,to;
	double dtol;
	//ppm = REAL(ppm)[0];
	//dmz = REAL(dmz)[0];
	int nx = x.size();
	int ny = y.size();

	struct idxStruct * pidxS =  (struct idxStruct *) calloc(nx,  sizeof(struct idxStruct));
	if (pidxS == NULL)
		stop("closeMatch/calloc: memory could not be allocated ! (%d bytes)\n", nx  * sizeof(struct idxStruct) );
	for (xi=0;xi < nx;xi++)
		pidxS[xi].from = ny+1;

	for (yi=0;yi < ny;yi++) {

		dtol=y[yi]*ppm*0.000001;
		if(dtol<dmz){
			dtol=dmz;
		}

		lb = lowerBound(y[yi] - dtol, x, lastlb, nx-lastlb);
		if (lb < nx-1)
			lastlb=lb;

		if (lb >= nx-1){
			lb=nx-1;
			ub=nx-1;
		} else
			ub = upperBound(y[yi] + dtol, x, lb, nx-lb);

		if (ub > nx-1)
			ub = nx -1;

		//    Rprintf("yi %d lb %d  ub %d \n",yi, lb,ub);

		for (xi=lb;xi <= ub;xi++) {
			if (fabs(y[yi] - x[xi]) <= dtol) {
				//             Rprintf("  -> Match xi %d \n",xi);
				if (yi < pidxS[xi].from)
					pidxS[xi].from = yi;
				if (yi > pidxS[xi].to)
					pidxS[xi].to = yi;
				//             Rprintf("xi %d from %d  to %d \n",xi, pidxS[xi].from, pidxS[xi].to);
			}
		}
	}

	IntegerVector res(nx,NA_INTEGER);

	for (xi=0;xi < nx;xi++) {

		// no match
		if (pidxS[xi].from == ny +1 && pidxS[xi].to == 0)
			continue;

		txi = xidx[xi] -1;
		from = pidxS[xi].from;
		to = pidxS[xi].to;

		// single match
		if (pidxS[xi].from == ny +1)
			from=pidxS[xi].to;
		if (pidxS[xi].to == 0)
			to=pidxS[xi].from;

		//Checking which point is the closest.
		double mindist=10;
		int minindex=-1;
		for (yi=from;yi <= to;yi++) {
			if(fabs(y[yi] - x[xi])<mindist){
				minindex=yi;
			}
		}
		res[txi] = yidx[minindex];
	}
	free(pidxS);
	return(res);
}

