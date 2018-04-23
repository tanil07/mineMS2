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

//FUnciton ont used anymore.
// std::map<Vertexl,Vertexel> subgraph_container::createExportLattice(latticeexport& elat){
//
//     resetIdxMap();
//
//     //Initialisation
//     std::map<Vertexl,Vertexel> vmap = latticeMapping(elat);
//     //std::cout << "mapping done" << std::endl; //DEBUG
//
//     //We generate the ids, on for the datasets, one for the nodes.
//     int count_M = 0;
//     int count_I = 0;
//
//     std::map< Vertexl, patternIdx> map_l_pat = getMapping();
//     //std::cout << "2nd mapping done" << std::endl; //DEBUG
//     //A mapping which link node is on the original
//
//     for(auto it = vmap.begin();it != vmap.end();it++){
//         //The new value is generated
//         std::string cid = "";
//         //std::cout << "npat: "<<lat[it->first].key << " addr "<<it->first; //DEBUG
//         if(lat[it->first].item){
//             //std::cout << "item";
//             count_I++;
//             cid = "G"+int_to_string(count_I);
//             elat[it->second].item = 1;
//         }else{
//             //std::cout << "pattern";
//             count_M++;
//             cid = "M"+int_to_string(count_M);
//             elat[it->second].item = 0;
//         }
//         elat[it->second].id = cid;
//
//         //Determination of the other node information.
//         if(lat[it->first].item){
//             elat[it->second].noccs  = 1;
//             elat[it->second].num_losses  = 0;
//             elat[it->second].noccs_unique = 1;
//         }else{
//             patternIdx patId = map_l_pat[it->first];
//             elat[it->second].noccs  = get(patId).numOccs();
//             //Correcting for the loss.
//             elat[it->second].num_losses  = get(patId).sizeGraph()-1;
//             elat[it->second].noccs_unique  = get(patId).numUniqueOccs().size();
//         }
//         //std::cout <<"attr" << " "; //DEBUG
//
//         //Currently the socre is initialized to 0:
//         elat[it->second].score = 0;
//
//         //In this part we add the corresponding edge.
//         graphTraitl::out_edge_iterator bo,eo;
//         boost::tie(bo,eo) = boost::out_edges(it->first,lat);
//
//         for(;bo!=eo;bo++){
//             auto nsource = vmap[boost::source(*bo,lat)];
//             auto ntarget = vmap[boost::target(*bo,lat)];
//
//             //The edge is added
//             boost::add_edge(nsource,ntarget,elat);
//         }
//     }
//
//     return vmap;
// }


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


//
// //path_lattice is a file name indicating all the values of the dataset.
// //path_path_patterns is a directory of output ofthe data
// void subgraph_container::exportResultingLattice(std::string path_lattice,std::string path_patterns,
//                                                 std::ostream& of){
//     //The lattice is created.
//     latticeexport elat;
//     std::map<Vertexl,Vertexel> vmap = createExportLattice(elat);
//
//     std::map< Vertexl, patternIdx> map_l_pat = getMapping();
//
//     //Node properties are defined.
//     boost::dynamic_properties dp(boost::ignore_other_properties);
//     dp.property("id", boost::get(&export_node_info::id, elat));
//     dp.property("noccs", boost::get(&export_node_info::noccs, elat));
//     dp.property("num_losses", boost::get(&export_node_info::num_losses, elat));
//     dp.property("score", boost::get(&export_node_info::score, elat));
//     dp.property("item", boost::get(&export_node_info::item, elat));
//     dp.property("i", boost::get(&export_node_info::item, elat));
//     //It's not exported, the number of unique occurences is calculated in the graph.
//
//
//     std::ofstream ofile(path_lattice);
//     boost::write_graphml(ofile, elat, dp, true);
//     ofile.close();
//
//
//     //The pattern are exported using this method.
//     graphTraitl::vertex_iterator bv,ev;
//     //std::cout << "Writing patterns." << std::endl;
//     boost::tie(bv,ev) = boost::vertices(lat);
//     for( ;bv != ev ; bv++){
//         if(lat[*bv].item) continue;
//         std::string idname = elat[vmap[*bv]].id;
//
//         std::string path_file = (path_patterns+"/"+idname+".graphml");
//         //std::cout << "WWriting: "<<path_file <<std::endl;
//         frag_pattern& fp = get(map_l_pat[*bv]);
//         //The node is exported.
//         fp.write_graphml(path_file);
//         //std::cout << "Written: "<<path_file <<std::endl;
//     }
//     //of << "Patterns written." << std::endl;
// }



