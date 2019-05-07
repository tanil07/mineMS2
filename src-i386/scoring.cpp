#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;


//Find the MST directly from the edge data.frame uniquely for a fragmentation patterns.
// [[Rcpp::export]]
DataFrame MST(DataFrame edges,NumericVector scores, int num_vertices){
	IntegerVector efrom = edges["from"];
	IntegerVector eto = edges["to"];
	IntegerVector elab = edges["lab"];

	//IntegerVector from(num_vertices);
	IntegerVector from(num_vertices,NA_INTEGER);
	IntegerVector lab(num_vertices,NA_INTEGER);
	NumericVector bestScore(num_vertices,NA_REAL);


	for(int i=0;i<efrom.size();i++){
		if(NumericVector::is_na(bestScore[eto[i]-1])||
     (scores[elab[i]-1]<bestScore[eto[i]-1])){
			bestScore[eto[i]-1] = scores[elab[i]-1];
			from[eto[i]-1] = efrom[i];
			lab[eto[i]-1] = elab[i];
		}
	}
	return(DataFrame::create(Named("from")=from,
                   Named("to")=seq_len(num_vertices),
                   Named("lab") = lab,
                   Named("score")=bestScore));
}

//TO MODIFY TO SCORE THE PATTERN
// [[Rcpp::export]]
double scorePattern(DataFrame edges, NumericVector scores, int num_vertices){
	DataFrame mst_edges = MST(edges,scores,num_vertices);
	IntegerVector labs = mst_edges["lab"];
	double totscore = 0;

	for(IntegerVector::iterator it=labs.begin();it != labs.end(); it++){
		if(!IntegerVector::is_na(*it)){
			totscore += scores[(*it)-1];
		}
	}
	//TODO MODIFY THIS.
	return totscore;
}
