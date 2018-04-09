#ifndef MASS_GRAPH_H
#define MASS_GRAPH_H

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "common.h"


class mass_graph
{
    public:
    //Constructors
        //Test graphs
        mass_graph(int);


        //TODO constructor form a R object.

        virtual ~mass_graph();

        //Getter
        graph& get_g();
        Vertex& get_precursor();

        //Used to build the path tree.
        VisitMap buildVisitMap();
        std::vector<Vertex> roots();

        //visu function
        std::string to_string();

        //
    protected:
    private:
        //The boost graph object.
        graph g;
        Vertex precursor;
};

#endif // MASS_GRAPH_H
