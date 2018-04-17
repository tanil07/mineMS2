#include <Rcpp.h>

#include "mass_dag_miner.h"
#include "mass_graph.h"

using namespace Rcpp;
// Function which actually perform the mining.
// [[Rcpp::export]]
void mineClosedDags(List mass_dags_df_list,LogicalVector processing,IntegerVector num,
                    IntegerVector k,IntegerVector size_min,LogicalVector prec_only){

	//Converting the mass graphs into an appropriate forms.
	std::vector < mass_graph > mgs;
	int num_graphs = mgs.size();
	for(int i = 0; i < num_graphs; i++){
		Rcpp::List temp_list(mass_dags_df_list[i]);
		Rcpp::DataFrame df_vertices(temp_list[0]);
		Rcpp::DataFrame df_edges(temp_list[1]);
		mass_graph temp_mg(df_vertices,df_edges);
		mgs.push_back(temp_mg);
	}

	int csize_min = as<int>(size_min);
	int cnum = as<int>(num);
	bool cprec_only = as<bool>(prec_only);
	int ck = as<int>(k);

//Give the minimum number of occurence of a pattern
	std::ofstream output_stream;
	output_stream.open("C:/Users/AD244905/Documents/ms-ms/res_mining.txt");
	mass_dag_miner mgf(mgs,ck,cprec_only,std::cout);//output_stream);
	mgf.setSizeMin(csize_min);
	mgf.mineFrequentCompleteDag(cnum,output_stream);
//
// auto scont = mgf.get_container();
//
// scont.exportResultingLattice(path_lattice,path_dags,output_stream);
//
    Rcpp::Rcout << "PROCESSING FINISHED !!!";
	output_stream.close();
	//return 0;
}

