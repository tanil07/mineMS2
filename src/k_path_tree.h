#ifndef K_PATH_TREE_H
#define K_PATH_TREE_H

#include "common.h"

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/array.hpp>

#include "mass_graph.h"
#include "adjacencyGraph.h"


#include <vector>




class lattice_node;
class k_path_tree
{
    public:
        k_path_tree(int k);
        virtual ~k_path_tree();
        void add_graph(mass_graph& G, int gid, bool=false);
        std::vector<Vertext> find_predecessors(Vertext c);
        int get_k();

        ktree& get_t();
        short get_dist_vertex(Vertext c);
        MapOccurrences& get_occs();
        Vertext get_root();
        adjacencyGraph& get_adj();

        //Removing nodes which are non frequent
        void filter_frequent_nodes(int noccs);

        //Creating 1-edge.
        void to_string(std::ostream& of);
        std::vector<Extension> getExtensions(Vertexp v,Vertext vt);

        //All the actions which take place after reading the graph
        void post_processing();

        //Initialisation of patterns.
        std::vector<lattice_node> constructOneEdgeGraphs(std::ostream& of,bool=true);

        //Create the values

        //Getting the nodes of a kpath tree
        short getLab(Vertext v) const;
        std::vector<short> getLabSuccs(Vertext v) const;
    protected:

        //Here the number of value to store.
        void update_pos_adv(int lab, std::vector<Vertex>& pfr,std::vector<int>& plabs,
                                 int lpath, IndexMap& idx_vertex,VisitMap& vm, int gid, graph& G);

        //Here the number of value to roll_back
        void update_pos_back();


        //Given a list of labels return the correct nodes.
        Vertext find_pos(std::vector<short>);




    private:
        int k;
        Vertext root;
        //TO DO move this to the graph miner when possible.
        std::vector<mass_graph> dags;

        //pos is a list which store all the elements from the path to the root.
        std::vector<std::vector<Vertext> > pos;
        ktree t;

        //It stores the occurrences under the form of a set of struct
        MapOccurrences moccs;

        //List of traingles dependent of the graph too

        //A list giving the adjacency property of the graph
        adjacencyGraph adj;

        //Internal functions
        Vertext get_node(Vertext,int,int);

};


#endif // K_PATH_TREE_H
