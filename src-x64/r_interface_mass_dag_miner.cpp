#include <Rcpp.h>

#include "mass_dag_miner.h"
#include "mass_graph.h"

using namespace Rcpp;
// Function which actually perform the mining.
// [[Rcpp::export]]
Rcpp::List mineClosedDags(List& vertices_list,List& edges_list,LogicalVector& processing,IntegerVector num,
                    IntegerVector k,IntegerVector size_min,LogicalVector prec_only){

	//Converting the mass graphs into an appropriate forms.
	std::vector < mass_graph > mgs;
	int num_graphs = vertices_list.size();

	for(int i = 0; i < num_graphs; i++){
		Rcpp::DataFrame df_vertices(vertices_list[i]);
		Rcpp::DataFrame df_edges(edges_list[i]);
		mass_graph temp_mg(df_vertices,df_edges);
		mgs.push_back(temp_mg);
	}
	Rcerr << "Converted the "<<num_graphs << " graphs."<<std::endl;

	int csize_min = as<int>(size_min);
	int cnum = as<int>(num);
	bool cprec_only = as<bool>(prec_only);
	int ck = as<int>(k);

	//Graph parsing.
	std::ofstream output_stream;
	mass_dag_miner mgf(mgs,ck,cprec_only,std::cout);
	Rcerr << "k-Path tree contructed."<<std::endl;
	Rcerr << "Mining frequent subgraphs. This step may be quite long"<<std::endl;
	//Frequent subgraph mining
	mgf.setSizeMin(csize_min);
	mgf.mineFrequentCompleteDag(cnum,std::cout);
	Rcerr << "Graph mining finished ."<<std::endl;

	//We return a treillis
	Rcpp::List res_list = (mgf.get_container().exportMinedPatternsRcpp());

	output_stream.close();
	return res_list;
}

