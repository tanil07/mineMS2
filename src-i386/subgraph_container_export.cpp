#include <iostream>
#include <tuple>
#include <string>
#include <iomanip>

#include <Rcpp.h>


#include "subgraph_container.h"


//THis function add a node to exlat for each node of the lattice and return a dictionnary of nodes.
std::map<Vertexl,int> subgraph_container::latticeMapping(){

	std::map<Vertexl,int>  res;
    //We pass on the lattice
    graphTraitl::vertex_iterator bv,ev;
    boost::tie(bv,ev) = boost::vertices(lat);

    //Ensuring that items are lawyas the first nodes.
    int countI=0;
    int countM = num_items;


    //A node is added for each node of the old lattice.
    for(;bv != ev; bv++){
    	if(lat[*bv].item){
        	res.insert(std::make_pair(*bv,countI));
    		countI++;
    	}else{
    		res.insert(std::make_pair(*bv,countM));
    		countM++;
    	}
    }
    return res;
}


//This function define a mapping between the node of the initial lattice and the node of the new lattices.

std::string int_to_string(int val){
    std::ostringstream convert;
    convert << val;
    return convert.str();
}

Rcpp::List subgraph_container::exportMinedPatternsRcpp(){

	//latticeexport elat;
	resetIdxMap();

	//Initialisation
	std::map<Vertexl, int>  vmap = latticeMapping();

	//We generate the ids, on for the datasets, one for the nodes.
	int count_M = 0;
	int count_I = 0;
	int tot_count = 0;

	//List storing the pattern as a data.frame
	Rcpp::List list_patterns(num_patterns);

	//List sotring the id of the items.
	Rcpp::IntegerVector vec_items(num_items);

	//Supplementary nodes infos which will be returned.
	Rcpp::IntegerVector num_occs(vmap.size());
	Rcpp::IntegerVector num_occs_unique(vmap.size());
	Rcpp::IntegerVector num_losses(vmap.size());
	Rcpp::LogicalVector is_item(vmap.size());


	//The id of the pattern or the position of the pattern.
	Rcpp::IntegerVector sid(vmap.size());
	Rcpp::IntegerVector gid(vmap.size());


	//Score can be directly initalised
	Rcpp::NumericVector score(vmap.size(),0.0);



	//Vertrox used for the edge list.
	Rcpp::IntegerVector from(boost::num_edges(lat));
	Rcpp::IntegerVector to(boost::num_edges(lat));
	int count_edge = 0;


	std::map< Vertexl, patternIdx> map_l_pat = getMapping();
	//A mapping which link node is on the original

	for(auto it = vmap.begin();it != vmap.end();it++){
		int pos = it-> second;
		//The new value is generated
		if(lat[it->first].item){

			//Putting it into the vectors.
			num_occs[pos] = 1;
			num_occs_unique[pos] = 1;
			num_losses[pos] = 0;
			is_item[pos] = true;

			//The motif ID
			sid[pos] = count_I;

			//If it was extracted into the lattice.

			//We append the index of the item to the data.
			vec_items[count_I] = pos+1;
			count_I++;

		}else{
			// elat[pinf.first].item = false;
			//
			// //The ID is always added to the dataset.
			// elat[pinf.first].id = pinf.second;

			patternIdx& patId = map_l_pat[it->first];

			//Pattern is extracted.
			frag_pattern& cp  = get(patId);

			//Putting it into the vectors.
			num_occs[pos] = cp.numOccs();
			num_occs_unique[pos] = cp.numUniqueOccs().size();
			num_losses[pos] = cp.sizeGraph()-1;
			is_item[pos] = false;
			sid[pos] = count_M+1;

			//Pattern is transformed into a list and pushed ot the list of motifs.
			list_patterns[count_M] = cp.as_igraph_data_frame();
			count_M++;
		}

		// //Currently the socre is initialized to 0:
		// elat[pinf.first].score = 0;
		//In this part we add the edges to the dataset.
		graphTraitl::out_edge_iterator bo,eo;
		boost::tie(bo,eo) = boost::out_edges(it->first,lat);

		for(;bo!=eo;bo++){
			from[count_edge]=vmap[boost::source(*bo,lat)]+1;
			to[count_edge]=vmap[boost::target(*bo,lat)]+1;
			//The edge is added
			//boost::add_edge(nsource,ntarget,elat);
			count_edge++;
		}

		gid[pos] = pos+1;
		tot_count++;
	}


	//The node_infos data.frame is created
	Rcpp::DataFrame nodes_infos = Rcpp::DataFrame::create(
		Rcpp::Named("gid") = gid,
		Rcpp::Named("sid") = sid,
		Rcpp::Named("item") = is_item,
		Rcpp::Named("num_occs") = num_occs,
		Rcpp::Named("num_occs_unique") = num_occs_unique,
		Rcpp::Named("num_losses") = num_losses,
		Rcpp::Named("score") = score
	);


	//The edge infos data.frame is created
	Rcpp::DataFrame edges_infos = Rcpp::DataFrame::create(
		Rcpp::Named("from") = from,
		Rcpp::Named("to") = to
	);

	//The final resl list is created.

	Rcpp::List res_list = Rcpp::List::create(
		Rcpp::Named("nodes")=nodes_infos,
		Rcpp::Named("edges")=edges_infos,
		Rcpp::Named("patterns")=list_patterns,
		Rcpp::Named("items")=vec_items
	);
	return res_list;
}
