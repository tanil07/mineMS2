#ifndef MASS_GRAPH_H
#define MASS_GRAPH_H

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <Rcpp.h>

#include "common.h"


class mass_graph
{
    public:
    //Constructor
    	mass_graph(Rcpp::DataFrame df_nodes,Rcpp::DataFrame df_edges);


        //TODO constructor form a R object.

        virtual ~mass_graph();

        //Getter
        graph& get_g();
        Vertex& get_precursor();

        Vertexp get_vertex_from_gid(short gid);

        //Used to build the path tree.
        VisitMap buildVisitMap();
        std::vector<Vertex> roots();

    protected:
    private:
        //The boost graph object.
        graph g;
        Vertex precursor;
};

#endif // MASS_GRAPH_H
