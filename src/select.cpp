#include <Rcpp.h>
using namespace Rcpp;

//C++ function used by the select function.

// [[Rcpp::export]]
IntegerVector select_patterns_from_spectra(List pat_list, int sid){
	//IntegerVector res;
	std::set<int> res;

	for(int pi = 0;pi <pat_list.size();pi++){
		S4 pat = pat_list[pi];
		IntegerMatrix occs = pat.slot("occurrences");
		IntegerVector gid = occs(_,0);
		for(IntegerVector::iterator it=gid.begin();it!=gid.end();it++){
			if((*it)==sid) (res).insert(pi+1);
		}
	}


	IntegerVector res_int(res.begin(),res.end());
	return res_int;
}

// [[Rcpp::export]]
List patterns_from_spectra(List pat_list,int num_spectra){
	std::vector< std::set < int > > res(num_spectra);
	//Rcout << "num_spectra:"<<num_spectra<<std::endl;


	for(int i =0;i<num_spectra;i++){
		std::set<int> temp;
		res[i] = temp;
	}


	for(int pi = 0;pi <pat_list.size();pi++){
		//Rcout << pi <<"_";
		S4 pat = pat_list[pi];
		//Rcout<<"pat_sel";
		IntegerMatrix occs = pat.slot("occurrences");
		IntegerVector gid = occs(_,0);
		//Rcout << gid.size() <<"_";
		for(IntegerVector::iterator it=gid.begin();it!=gid.end();it++){
			//Rcout <<" ii"<< (*it) <<"_";
			(res[(*it)-1]).insert(pi+1);
			//Rcout <<" ins"<< (*it) <<"_";
		}
		//Rcout << std::endl;
	}

	//Rcout <<std::endl<<"res : ";

	List res_list = List(num_spectra);
	for(size_t i = 0;i < res.size();i++){
		//Rcout <<i<<"_";
		IntegerVector temp(res[i].begin(),res[i].end());
		res_list[i]=temp;
	}
	return res_list;
}

