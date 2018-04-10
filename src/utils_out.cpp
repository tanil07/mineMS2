#include <Rcpp.h>
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



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
###Test ifnd equalGreaterM
FindEqualGreaterM(1:100,c(10.2,56.1,89))

###Test findLimDensity.
xsq <- seq(-10,10,by=0.01)
ds <- dnorm(xsq,mean = 2,sd = 2)+dnorm(xsq,mean = -4,1)
res <- findLimDensity(ds,400,2)
print(res)
res2 <- findLimDensity(ds,res[2],res[3])
print(res2)

plot(ds)
abline(v=res,col="red")
abline(v=res2,col="darkgreen")
*/
