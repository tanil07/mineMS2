#ifndef MASS_DAG_MINER_H
#define MASS_DAG_MINER_H

#include <stack>          // std::stack


#include "k_path_tree.h"
#include "subgraph_container.h"
#include "triangles_list.h"
#include "mass_graph.h"
#include "lattice_node.h"

class mass_dag_miner
{
    public:
        mass_dag_miner () ;
        mass_dag_miner (std::vector<mass_graph>&, int,bool, std::ostream&) ;
        virtual ~mass_dag_miner();
        void mineFrequentDag(int freq,std::ostream& os);
        void mineFrequentCompleteDag(int freq,std::ostream& of);
        void mineFrequentCompleteDag(int freq,std::set<short> vals,std::ostream& of);
        subgraph_container& get_container();

    protected:
    private:
        subgraph_container container;
        k_path_tree kt;
        triangles_list tl;
        int num_graph;
        std::stack<lattice_node> lattice_stack;
};

#endif // MASS_DAG_MINER_H
