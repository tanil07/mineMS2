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

//return a set of graph given and exemple number.
mass_graph::mass_graph(int i)
{
    int max_ex = 1;
    if(i>max_ex)
    {
        throw std::range_error("i should be lower than"+patch::to_string(max_ex));
        //std::cerr << "No example " <<i <<"." << std::endl;
    }
    graph Gr;
    if(i == max_ex)
    {
        enum vertices {A,B,C,D,E,F,G,H,I,N};

        int num_vertices = N;
        graph Gr(num_vertices);

        //Making a mapping between the enum and the vertex.
        typedef std::map<int,Vertex> VertexIdMap;
        VertexIdMap node_map;


        boost::graph_traits<graph>::vertex_iterator bv,ev;
        int idxn = 0;
        for(boost::tie(bv,ev) = boost::vertices(Gr); bv!=ev; bv++)
        {
            node_map.insert(std::make_pair(idxn,*bv));
            idxn++;
        }

        //We encompass the array using boost
        boost::array<node_info,9> ninfos=
        {
            {
                {8,140.034,87.93},
                {7,123.031,3.60},
                {6,122.024,3.50},
                {5,105.021,1.63},
                {4,98.984,0.77},
                {3,96.064,0.72},
                {2,96.044,0.66},
                {1,95.037,0.68},
                {0,66.034,0.5}
            }
        };



        //We add all the nodes
        for(int ie=0; ie < N; ie++)
        {
            //For each node
            Gr[ node_map[ie] ].mz=ninfos[ie].mz;
            Gr[ node_map[ie] ].intensity=ninfos[ie].intensity;
            Gr[ node_map[ie] ].lab=ninfos[ie].lab;
        }

        boost::array<edge_info,6> einfos=
        {
            {
                {3},
                {3},
                {32},
                {136},
                {62},
                {10}
            }
        };

        typedef std::pair<int,int> edge;
        boost::array<edge,6> etab= {edge(A,C),edge(B,D),
                                    edge(A,G),edge(A,I),edge(C,I),
                                    edge(G,I)
                                   };
        boost::graph_traits<graph>::edge_descriptor e;

        for(int i=0; i<etab.size(); i++)
        {

            add_edge(node_map[etab[i].first],
                     node_map[etab[i].second],einfos[i],Gr);
        }
        g = Gr;
    }               //ctor
}


mass_graph::mass_graph(Rcpp::DataFrame df_nodes,Rcpp::DataFrame df_edges){
	graph g;
	Rcpp::NumericVector mz = df_nodes["mz"];
	Rcpp::NumericVector rel_int = df_nodes["rel_int"];
	Rcpp::LogicalVector prec = df_nodes["prec"];

	//Adding the nodes.
	std::vector<Vertex> mapped_vec(mz.size());
	for(int i=0;i<mz.size();i++){
		Vertex nvert = boost::add_vertex(g);
		g[nvert].intensity = rel_int[i];
		g[nvert].mz = mz[i];
		g[nvert].lab = i+1;
		mapped_vec[i] = nvert;
		if(prec[i]){
			precursor = nvert;
		}
	}

	//We start by adding the files as nodes.
	Rcpp::IntegerVector from = df_edges["from"];
	Rcpp::IntegerVector to = df_edges["to"];
	Rcpp::IntegerVector elab = df_edges["lab"];
	for(int i=0;i<mz.size();i++){
		graphTraits::edge_descriptor nedge = boost::add_edge(mapped_vec[from[i]-1],
                                                       mapped_vec[to[i]-1],g).first;
		g[nedge].lab = elab[i];
	}
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
    //We iterate on all the vertices.

    for(boost::tie(bv,ev)=boost::vertices(g); bv!=ev; bv++)
    {
        //We push the values
        std::vector<std::pair<Vertex,bool> > AdjacentNodes;
        boost::graph_traits<graph>::adjacency_iterator ai,ai_end;
        for (boost::tie(ai, ai_end) = adjacent_vertices(*bv, g);
                ai != ai_end; ++ai)
        {
            std::pair<Vertex,bool> temp_pair(*ai,false);
            AdjacentNodes.push_back(temp_pair);
        }
        //At this step the map can be pushed in the value.
        vmap[*bv]=AdjacentNodes;
    }
    return vmap;
}




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


//We define the getter and setter function
graph& mass_graph::get_g()
{
    return g;
}
Vertex& mass_graph::get_precursor()
{
    return precursor;
}


//String repsentation of the mass graph
std::string mass_graph::to_string()
{
    //We first get the edge in bfs order.
    boost::graph_traits<graph>::vertex_iterator vgbeg,vgend;
    std::vector<std::string> res(100,"NA");
    int size_str = 0;
    res[0] = "Nodes :\n";

    for(boost::tie(vgbeg,vgend)=boost::vertices(g); vgbeg!=vgend; vgbeg++)
    {
        //We get all the vertices attributes.
        if((size_str+6)>=res.size())
        {
            res.resize(2*res.size());
        }
        res[++size_str] = patch::to_string(g[*vgbeg].lab);
        res[++size_str] = patch::to_string(g[*vgbeg].mz);
        res[++size_str] = patch::to_string(g[*vgbeg].intensity);
        res[++size_str] = patch::to_string(boost::in_degree(*vgbeg,g));
        res[++size_str] = "\n";
    }

    res[++size_str]="Edges :\n";

    boost::graph_traits<graph>::edge_iterator egbeg,egend;
    for(boost::tie(egbeg,egend)=boost::edges(g); egbeg!=egend; egbeg++)
    {
        //We get all the vertices attributes.
        if((size_str+3)>=res.size())
        {
            res.resize(2*res.size());
        }
        res[++size_str] = patch::to_string(g[*egbeg].lab);
        res[++size_str] = "\n";
    }
    res.resize(size_str);
    std::string rstr = boost::algorithm::join(res," ");

    return rstr;
}
