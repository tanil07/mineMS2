#include <Rcpp.h>
#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <iostream>

using namespace Rcpp;

//Return the firstg possible smaller value for each value (xmcs inspired)
// [[Rcpp::export]]
IntegerVector FindEqualGreaterM(NumericVector inv, NumericVector values){

	IntegerVector index(values.size());

	int i, idx = 0;

	for (i = 0; i < values.size(); i++)
	{
		while (idx < inv.size() && inv[idx] < values[i])
			idx++;
		index[i] = idx;
	}
	return index;
}

int sign_d(double x)
{
	return((x>0)-(x<0));
}

//Find the limits of a peak by following the slopes on both sides.
// [[Rcpp::export]]
IntegerVector findLimDensity(NumericVector seq,int istart,int state)
{
	int size = seq.size();
	int i,difference;
	int linflex=istart;
	int rinflex=istart;

	for(i=istart; i<size-1; i++)
	{
		difference=sign_d(seq[i+1]-seq[i])-sign_d(seq[i]-seq[i-1]);
		if(difference>=1)  //inflex point
		{
			if(state==2)
			{
				state=1;
				rinflex=i;
				istart=i;
				break;
			}
			else//First encounter with an inflex
			{
				state=1;
				linflex=i;
			}
		}
		if(difference==-2&&(state)==1)
		{
			state=2;
		}
	}
	IntegerVector index = IntegerVector::create(linflex+1,rinflex+1,state);
	return index;
}


//void formulaFromString(std::string rstrformula,std::vector<std::string>& names_atoms,std::vector<int>& count_atoms)
//Parse a chemical string.
// [[Rcpp::export]]
List formulaFromString(std::string formula,std::vector<std::string> names_atoms){

	std::vector<int> count_atoms(names_atoms.size(),0);
	int bpos = -1;
	int ilen = 0;
	int bposnum = -1;
	int ilennum = 0;
	int incount = 0;
	int cpos = -1;
	size_t i;
	std::vector<std::string>::iterator atmpos;
	for (i = 0 ; i < formula.length(); i++)
	{
		if(std::isupper(formula[i]))
		{
			incount = 0;


			if(bpos != -1)
			{
				//We get the index of the atom name if existing
				std::string atmname = formula.substr(bpos,ilen);
				atmpos = std::find(names_atoms.begin(), names_atoms.end(),atmname);

				//If it does not exists we add it
				if(atmpos==names_atoms.end()){
					names_atoms.push_back(atmname);
					count_atoms.push_back(1);
					//Increamenting the counter at the same time
					atmpos = names_atoms.end() - 1;
				}
				cpos = atmpos - names_atoms.begin();
				bpos = (int)i;
				ilen = 1;
			}
			else
			{
				bpos = 0;
				ilen = 1;
			}
			if(i!=0)
			{

				//The position is furnished by the vector
				if(ilennum != 0)
				{
					//If superior to 1 we change it
					count_atoms[cpos] = atoi(formula.substr(bposnum,ilennum).c_str());
				}
				else
				{
					count_atoms[cpos] = 1;
				}
				ilennum = 0;
			}

		}
		if(std::islower(formula[i]))
		{
			ilen++;
		}
		if(std::isdigit(formula[i]))
		{
			if(incount)
			{
				ilennum++;
			}
			else
			{
				incount = 1;
				bposnum = (int)i;
				ilennum = 1;
			}
		}
	}

	std::string atmname = formula.substr(bpos,ilen);
	atmpos = std::find(names_atoms.begin(), names_atoms.end(),atmname);
	cpos = atmpos - names_atoms.begin();
	//If it does not exists we add it
	if(atmpos==names_atoms.end()){
		names_atoms.push_back(atmname);
		if(ilennum != 0)
		{
			count_atoms.push_back(atoi(formula.substr(bposnum,ilennum).c_str()));
		}else{
			count_atoms.push_back(1);
		}
		//Increamenting the counter at the same time
		atmpos = names_atoms.end() - 1;
	}else{
		if(ilennum != 0)
		{
			count_atoms[cpos] = atoi(formula.substr(bposnum,ilennum).c_str());
		}else{
			count_atoms[cpos] = 1;
		}
	}
	return List::create(Named("Atom") = wrap(names_atoms),
                     Named("Num") = wrap(count_atoms));
}


// [[Rcpp::export]]
NumericVector disjointBins(NumericVector points, NumericVector lower_lim, NumericVector upper_lim, NumericVector mean_bin) {
	int N = points.size();

	NumericVector bins(N+1);
	for(int j=0;j<(N+1);j++) bins[j]=0;

	int pi = 0;
	int pp = 0;
	int num_match = 0;
	int in_inter = 0;
	while((pi < lower_lim.size())&(pp<N)){
		if((points[pp]>=lower_lim[pi])&(points[pp]<upper_lim[pi])){
			//We check if the value is already in a bin or not.
			in_inter = 1;
			num_match++;
			if(bins[pp]!=0){
				if(abs(mean_bin[bins[pp]-1]-points[pp])>=abs(mean_bin[pi]-points[pp])){
					bins[pp]=pi+1;
				}
				pi++;
			}else{
				bins[pp]=pi+1;
				pi++;
			}
		}else{
			if(points[pp]<lower_lim[pi]){
				if(bins[pp] ==0) bins[pp]=NA_REAL;
				if(in_inter==1){
					in_inter=0;
					pi = pi - num_match;
					num_match=0;
				}
				pp++;
			}else if(points[pp]>=upper_lim[pi])
				pi++;
		}
	}
	//Case where the last element have been binned
	if(bins[pp]!=0){
		pp++;
	}

	for(int p=pp;p<(N+1);p++){
		bins[p]=NA_REAL;
	}
	//Rprintf("Out");
	return bins;
}


// [[Rcpp::export]]
List checkInter(NumericVector a_min, NumericVector a_max, NumericVector b_min, NumericVector b_max) {
	int Na = a_min.size();
	int Nb = b_min.size();
	NumericVector bin_min(Nb);
	NumericVector bin_max(Nb);
	NumericVector hbin_min(Nb);
	NumericVector hbin_max(Nb);
	NumericVector mean_a(Na);
	for(int j=0;j<Nb;j++){
		bin_min[j]=NA_REAL;
		bin_max[j]=NA_REAL;
		hbin_min[j]=NA_REAL;
		hbin_max[j]=NA_REAL;
	}
	for(int j=0;j<Na;j++){
		mean_a[j]= (a_min[j]+a_max[j])/2;
	}
	int pa = 0;
	int pb = 0;
	int in_inter = 0;
	while((pa < Na)&(pb<Nb)){
		//printf("%d %d %d %d\n", pa, pb, Na, Nb);
		//AAA*****
		//****BBB*
		if(a_max[pa]<=b_min[pb]){
			pa++;
			//*****AAA
			//*BBB****
		}else if(a_min[pa]>=b_max[pb]){ // before, it was strictly superior, do not know why
			//Last point was the last of an interval.
			if(in_inter==1){
				in_inter=0;
				bin_max[pb]=pa-1;
				pa=bin_min[pb];
			}
			pb++;
			//**AAA***
			//***BBB**
			//Case where there is an overlap to the left.
		}else if((a_min[pa]<b_max[pb])&(a_min[pa]>=b_min[pb])){
			if(in_inter==0){
				bin_min[pb]=pa;
				in_inter = 1;
			}
			pa++;
			//***AAA**
			//*BBB****
		}else if((a_max[pa]<b_max[pb])&(a_max[pa]>b_min[pb])){
			if(in_inter==0){
				bin_min[pb]=pa;
				in_inter = 1;
			}
			pa++;
			//***AA***
			//**BBBB**
		}else if((a_min[pa]>b_min[pb])&(a_max[pa]<b_max[pb])){
			if(in_inter==0){
				bin_min[pb]=pa;
				in_inter = 1;
			}
			pa++;
			//**AAAA**
			//***BB***
		}else if((a_min[pa]<=b_min[pb])&(a_max[pa]>=b_max[pb])){
			if(in_inter==0){
				bin_min[pb]=pa;
				in_inter = 1;
			}
			pa++;
		}
		//else{
			//printf("aaaaaaaaaaaahhhhh\n");
			//printf("%f %f %f %f\n", a_min[pa], b_min[pb], a_max[pa], b_max[pb]);
			//sleep(15);
			//exit(0);
		//}
	}
	if((pa==Na)&(in_inter==1)){
		bin_max[pb]=pa-1;
	}
	//We check if the mean is in the interval or not
	for(int j=0;j<Nb;j++){
		if(NumericVector::is_na(bin_min[j])){
			continue;
		}
		int beginning = 0;
		for(int i=bin_min[j];i<=bin_max[j];i++){
			if((mean_a[i]<=b_max[j])&(mean_a[i]>=b_min[j])){
				if(beginning==0){
					hbin_min[j]=i;
					beginning=1;
				}
				hbin_max[j]=i;
			}
		}
	}



	//Adjusting to R index
	for(int j=0;j<Nb;j++){
		bin_min[j]++;
		bin_max[j]++;
		hbin_min[j]++;
		hbin_max[j]++;
	}

	return(List::create(Rcpp::Named("bmin") = bin_min,
                     Rcpp::Named("bmax") = bin_max,
                     Rcpp::Named("hbmin") = hbin_min,
                     Rcpp::Named("hbmax") = hbin_max));
}



//Correct formula mistake.
bool intersect(double b1min,double b1max,double b2min,double b2max){
  double mmin = std::max(b1min,b2min);
  double mmax = std::min(b1max,b2max);
  if(mmax>=mmin){
    return(true);
  }
  return(false);
}

IntegerVector extend_match(double valmin, double valmax, int pos, NumericVector bmin, NumericVector bmax){
  int pleft = pos;
  int pright = pos;
  while(intersect(valmin,valmax,bmin[pleft],bmax[pleft]) & (pleft>=0)){
    pleft--; 
  }
  
  while(intersect(valmin,valmax,bmin[pright],bmax[pright]) & (pright<bmin.size())){
    pright++; 
  }
  IntegerVector res;
  res.push_back(pleft);
  res.push_back(pright);
  
  return res;
}

//finding the limiting range of the value.
IntegerVector bisect_search_borns(double valmin, double valmax, NumericVector bmin, NumericVector bmax){
  //Rcout <<"in" << std::endl;
  int bleft = 0;
  int bright = bmin.size()-1;
  int mid = (bleft+bright)/2;
  
  IntegerVector resn;
  resn.push_back(mid);
  resn.push_back(mid);
  
  //We try to find the first intersecting value.
  while(bleft<bright){
    if(bmax[mid]<valmin){
      bleft = mid+1;
    }else if(intersect(valmin,valmax,bmin[mid],bmax[mid])){
      return extend_match(valmin,valmax,mid,bmin,bmax);
    }else if(bmin[mid]>valmax){
      bright = mid-1;
    }
    mid = (bleft+bright)/2;
  }
  return resn;
}

// [[Rcpp::export]]
List find_combinations_ranges(NumericVector bmin, NumericVector bmax,
                              NumericVector cmzmax){
  
  int N = bmin.size();
  
  int i = 0;
  int j = 0;
  
  List to_return(N);
  for(int i=0;i<to_return.size();i++){
    
    IntegerVector v1(0);
    IntegerVector v2(0);
    List temp =List::create(Rcpp::Named("f1") = v1,Rcpp::Named("f2") = v2);
    to_return[i]=temp;
  }
  
  double mzmax = cmzmax[0];
  double hmax,hmin;
  IntegerVector res;
  while( (bmax[i] < mzmax) & (i< N)){
    j = i;
    
    while(((hmin=(bmin[j]+bmin[i])) < mzmax)& (j <N)){
      
      //We calculate the borns
      hmax = bmax[j]+bmax[i];
      
      //We find the interval limit
      res = bisect_search_borns(hmin,hmax,bmin,bmax);
      
      if(res[0]<res[1]){
        for(int ii = res[0]+1; ii < res[1]; ii++){
          List temp = to_return[ii];
          IntegerVector f1 = temp["f1"];
          IntegerVector f2 = temp["f2"];
          f1.push_back(i);
          f2.push_back(j);
          temp["f1"] = f1;
          temp["f2"] = f2;
          to_return[ii] = temp;
        }
      }
      j++;
    }
    i++;
  }
  
  return to_return;
}

