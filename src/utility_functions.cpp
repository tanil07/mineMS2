#include <iostream>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "utility_functions.h"
#include "common.h"

//Function which take a vertex a visitorMap and return
//the first values in the graph if necessary.
Vertex nextNode(VisitMap& vm,Vertex v){
    std::vector<std::pair<Vertex,short> >::iterator vb,ve;
    vb = vm[v].begin();
    ve = vm[v].end();
    //We get a list of the possible sucessors
    for(std::vector<std::pair<Vertex,short> >::iterator it=vb;it!=ve;it++){
        if(!(*it).second){
            (*it).second = 1;
            return (*it).first ;
        }
    }
    vb = vm[v].begin();
    ve = vm[v].end();
    //If no sucessors is found we get a list.
    for(auto it=vb;it!=ve;it++){
        it->second = 2;
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


//DEBUG PRINTING VISITOR MAP
void printVmap(VisitMap& vm,graph& g){
    graphTraits::vertex_iterator vb,ve;
    for(boost::tie(vb,ve) = boost::vertices(g);vb!=ve;vb++){
        std::cout << "vertex " << g[*vb].mz << " : ";
        for(auto itt=vm[*vb].begin();itt!=vm[*vb].end();itt++){
            std::cout << g[(*itt).first].mz << " ";
        }
        std::cout << std::endl;
    }
}


//<template T>


