#include <Rcpp.h>
using namespace Rcpp;

//C++ function used by the select function.

//
// IntegerVector select_patterns_from_spectra(List pat_list){
//
//
// }

// [[Rcpp::export]]
List patterns_from_spectra(List pat_list,int num_spectra){
	std::vector< std::set < int > > res(num_spectra);


	for(int i =0;i<num_spectra;i++){
		std::set<int> temp;
		res[i] = temp;
	}


	for(int pi = 0;pi <pat_list.size();pi++){
		S4 pat = pat_list[pi];
		IntegerMatrix occs = pat.slot("occurences");
		IntegerVector gid = occs(_,0);
		for(IntegerVector::iterator it=gid.begin();it!=gid.end();it++){
			(res[*it]).insert(pi+1);
		}
	}


	List res_list = List(num_spectra);
	for(size_t i = 0;i < res.size();i++){
		IntegerVector temp(res[i].begin(),res[i].end());
		res_list[i]=temp;
	}
	return res_list;
}

