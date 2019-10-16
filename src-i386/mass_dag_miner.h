#ifndef MASS_DAG_MINER_H
#define MASS_DAG_MINER_H

#include <stack>          // std::stack

#include "Rcpp.h"
#include "k_path_tree.h"
#include "subgraph_container.h"
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
        subgraph_container& get_container();
        void setSizeMin(int smin);
        int getSizeMin();

    protected:
    private:
        std::vector<mass_graph> mgs;
        k_path_tree kt;
        subgraph_container container;
        int num_graph;
        std::stack<lattice_node> lattice_stack;
        int sizeMin;
};

#endif // MASS_DAG_MINER_H
