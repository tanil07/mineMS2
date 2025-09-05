#include <iostream>
#include <tuple>
#include <string>
#include <iomanip>

#include <Rcpp.h>


#include "subgraph_container.h"



//This function define a mapping between the node of the initial lattice and the node of the new lattices.

std::string int_to_string(int val){
    std::ostringstream convert;
    convert << val;
    return convert.str();
}

//Second verison of the R export.
Rcpp::List subgraph_container::exportMinedPatternsRcpp(){

  //We generate the ids, on for the datasets, one for the nodes.

  //Rcpp::Rcerr << "IN";
  //We iterate on the full vector.
  int num_patterns = numPatterns();
  
  //Rcpp::Rcerr << "npat:" << num_patterns << " nitems: "<< num_items;
  //List storing the pattern as a data.frame
  Rcpp::List list_patterns(num_patterns);
  
  //List sotring the id of the items.
  Rcpp::IntegerVector vec_items(num_items);
  
  //Supplementary nodes infos which will be returned.
  Rcpp::IntegerVector num_occs(num_patterns);
  Rcpp::IntegerVector num_occs_unique(num_patterns);
  Rcpp::IntegerVector num_losses(num_patterns);
  Rcpp::LogicalVector is_item(num_patterns);
  
  //We store the iDs of the pattern evnetually becuse it is simple.
  Rcpp::IntegerVector gid(num_patterns);
  
  //Score can be directly initalised
  
  //Vertrox used for the edge list.
  int pos = 0;
  for(auto it = pmap.begin();it != pmap.end();it++){
    std::vector<frag_pattern> val = it->second;
    for(auto itt=val.begin();itt != val.end();itt++){
      frag_pattern& cp  = *itt;
      //Putting it into the vectors.
      num_occs[pos] = cp.numOccs();
      num_occs_unique[pos] = cp.numUniqueOccs().size();
      num_losses[pos] = cp.sizeGraph()-1;
      is_item[pos] = false;
      std::string id = "P" + std::to_string(pos);
      //Pattern is transformed into a list and pushed ot the list of motifs.
      list_patterns[pos] = cp.as_igraph_data_frame();
      pos++;
    }
  }
  
  
  //The node_infos data.frame of the lattice is created. is created
  Rcpp::DataFrame nodes_infos = Rcpp::DataFrame::create(
    Rcpp::Named("id") = gid,
    Rcpp::Named("item") = is_item,
    Rcpp::Named("num_occs") = num_occs,
    Rcpp::Named("num_occs_unique") = num_occs_unique,
    Rcpp::Named("num_losses") = num_losses
  );

  
  //The final results lids list is created.
  Rcpp::List res_list = Rcpp::List::create(
    Rcpp::Named("nodes")=nodes_infos,
    Rcpp::Named("patterns")=list_patterns
  );
  return res_list;
}

