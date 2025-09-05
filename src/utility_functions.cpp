#include <iostream>

#include "utility_functions.h"
#include "common.h"

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

//Function which take a vertex a visitorMap and return
//the first values in the graph if necessary.

//0: unknown 1: closed 2:visited
Vertex nextNode(VisitMap& vm,graph& G,Vertex v){
    graphTraits::out_edge_iterator bo,eo;

    for(boost::tie(bo,eo)=boost::out_edges(v,G);bo != eo; bo++){
        if(vm[*bo]!=1){
            vm[*bo]=1;
            return boost::target(*bo,G);
        }
    }

    //No successors is found
    for(boost::tie(bo,eo)=boost::out_edges(v,G);bo != eo; bo++){
            vm[*bo]=2;
    }
    return graphTraits::null_vertex();

}

void addOccs(MapOccurrences& moccs,Vertext vt,occ oc){
    bool inserted;
    MapOccurrences::iterator it;
    std::set<occ> temp_set;
    temp_set.insert(oc);
    boost::tie(it,inserted) = moccs.insert(std::pair<Vertext,std::set<occ> >(vt,temp_set));
    if(!inserted){
        moccs[vt].insert(oc);
    }
}




