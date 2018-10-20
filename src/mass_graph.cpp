#include <fstream>
#include <exception>
#include <string>
#include <iostream>

#include <boost/array.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <string>

#include "to_string.h"
#include "mass_graph.h"
#include "common.h"


//mass_graph::mass_graph()
//{
//    //ctor
//}


mass_graph::mass_graph(Rcpp::DataFrame df_nodes,Rcpp::DataFrame df_edges){
	graph gr;
	Rcpp::NumericVector mz = df_nodes["mz"];
	Rcpp::NumericVector rel_int = df_nodes["rel_int"];
	Rcpp::LogicalVector prec = df_nodes["prec"];

	//Adding the nodes.
	std::vector<Vertex> mapped_vec(mz.size());
	for(int i=0;i<mz.size();i++){
		Vertex nvert = boost::add_vertex(gr);
		gr[nvert].intensity = rel_int[i];
		gr[nvert].mz = mz[i];
		gr[nvert].lab = i+1;
		mapped_vec[i] = nvert;
		//Precursor is added when found.
		if(prec[i]){
			precursor = nvert;
		}
	}

	//Adding the edge
	if(df_edges.nrows()!=0){
		Rcpp::IntegerVector from = df_edges["from"];
		Rcpp::IntegerVector to = df_edges["to"];
		Rcpp::IntegerVector elab = df_edges["lab"];
		for(int i=0;i<from.size();i++){
			graphTraits::edge_descriptor nedge = boost::add_edge(mapped_vec[from[i]-1],
	                                                       mapped_vec[to[i]-1],gr).first;
			gr[nedge].lab = elab[i];
		}
	}
	g = gr;
}



mass_graph::~mass_graph()
{

}

std::vector<Vertex> mass_graph::roots(){
    std::vector<Vertex> roots;
    graphTraits::vertex_iterator vbeg,vend;
    for(boost::tie(vbeg,vend)=boost::vertices(g); vbeg!=vend; vbeg++)
    {
        if(boost::in_degree(*vbeg,g)==0)
        {
            roots.push_back(*vbeg);
        }
    }
    return roots;
}


VisitMap mass_graph::buildVisitMap()
{
    //We define the map
    VisitMap vmap;

    boost::graph_traits<graph>::vertex_iterator bv,ev;
    boost::graph_traits<graph>::edge_iterator be,ee;
    //We iterate on all the vertices.


    for(boost::tie(be,ee)=edges(g);be!=ee;be++){
        vmap[*be] = 0;
    }

    return vmap;
}

//    for(boost::tie(bv,ev)=boost::vertices(g); bv!=ev; bv++)
//    {
        //We push the values
//        std::vector<std::pair<Vertex,short> > AdjacentNodes;
//        boost::graph_traits<graph>::adjacency_iterator ai,ai_end;
//        for (boost::tie(ai, ai_end) = adjacent_vertices(*bv, g);
//                ai != ai_end; ++ai)
//        {
//            std::pair<Vertex,short> temp_pair(*ai,0);
//            AdjacentNodes.push_back(temp_pair);
//        }

//        vmap[*bv]=AdjacentNodes;
//    }
//    return vmap;



//TODO ::Load the graph from graphML file
//mass_graph::mass_graph(std::string filename)
//{
//    std::ifstream datafile(filename);
//    if (!datafile)
//    {
//        std::cerr << "No ./ " <<filename <<" file" << std::endl;
//        return 0;
//    }
//
//    //We read the graph ML file
//    graph g;
//
//}

Vertexp mass_graph::get_vertex_from_gid(short gid){
    IndexMap imap = boost::get(boost::vertex_index,g);
    graphTraits::vertex_iterator bv,ev;
    for(boost::tie(bv,ev)=boost::vertices(g);bv!=ev;bv++){
        if(gid==imap[*bv])     return *bv;
    }
    return graphTraits::null_vertex();

}


//We define the getter and setter function
graph& mass_graph::get_g()
{
    return g;
}
Vertex& mass_graph::get_precursor()
{
    return precursor;
}

