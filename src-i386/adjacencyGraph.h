#ifndef ADJACENCYGRAPH_H
#define ADJACENCYGRAPH_H

#include "common.h"
//#include "k_path_tree.h"


class k_path_tree;
//This graph as for sole purpose to handle the direct neighbouring grah.
class adjacencyGraph
{
    public:
        adjacencyGraph();
        adjacencyGraph(int n);

        //Adding some information in the graph
        void add_adj(Vertex v,graph& mg);
        void add_graph(graph& mg);
        void addKTreeVertices(k_path_tree& kt);
        void remove_nodes(std::vector<short>& labs);
        void keep_nodes(std::vector<short>& labs);
        graphadj& get_g();

        //Getting the labels ot add
        std::vector<std::pair <Vertext,short> > neighbours(short);

        virtual ~adjacencyGraph();
    protected:
    private:
        graphadj g;
        //Mapping between the label and the nodes.
        std::map<short,Vertexadj> lab_node;

        //Commodity method.
        Vertexadj getNode(short lab);
};

#endif // ADJACENCYGRAPH_H
