#include<iostream>
#include <fstream>

//DEBUG ONLY
#include <Rcpp.h>

#include "mass_dag_miner.h"
#include "lattice_node.h"

mass_dag_miner::mass_dag_miner(): kt(1)
{
    //ctor
}

mass_dag_miner::~mass_dag_miner()
{
    //dtor
}


void mass_dag_miner::setSizeMin(int smin){
	sizeMin = smin;
}

int mass_dag_miner::getSizeMin(){
	return sizeMin;
}

mass_dag_miner::mass_dag_miner(std::vector<mass_graph>& vmasses, int k, bool prec_only, std::ostream& os) : kt(k) {

	sizeMin = 2;
    num_graph = vmasses.size();
    mgs = vmasses;
    int counter = 0;
    for(auto it=vmasses.begin();it!=vmasses.end();it++){
        //Each graph is added to the k path tree which also store the occurences.
        kt.add_graph(*it,counter,prec_only);
        counter++;
    }

    kt.post_processing();
}


void mass_dag_miner::mineFrequentCompleteDag(int freq,std::ostream& of){
    ktree& t = kt.get_t();
    int init_nodes = boost::num_vertices(t);

    kt.filter_frequent_nodes(freq);

    //Rcpp::Rcerr <<"Initial k-path tree with: "<<boost::num_vertices(t) <<" nodes down from "<< init_nodes<< std::endl;

    std::vector<lattice_node> initialPatterns = kt.constructOneEdgeGraphs(of);
    Rcpp::Rcerr <<"Mining initialized with: "<<initialPatterns.size() <<" patterns."<< std::endl;
    
    bool frequent_children = true;
    int counv = 0;
    std::stack<short> max_child_occs;

    //Used for the closeness predicates.
    Rcpp::Rcerr << "Mined patterns: "<<std::endl;
    int insert_error=0;
    int counter = 0;
    for(auto it=initialPatterns.begin();it!=initialPatterns.end();it++){

        //We empty the current stack
        while(!lattice_stack.empty()){
            lattice_stack.pop();
            max_child_occs.pop();
        }

        //We initialize the value.
        lattice_stack.push(*it);

        //We store the number of occurences of the child for each node.
        max_child_occs.push(0);

        int LIM_ITER = 1000000;
        int current_1000 = 0;
        lattice_node next_node;
        int indicator = -1;
        while(counter<LIM_ITER){
            //We always consider the last elemnt of the stack.
            lattice_node& current_node = lattice_stack.top();

            //We first check if the current pattern is root.
            bool is_root = current_node.is_root();

            //We print the label of the current node
            if(((counv%1000)==0)&(((counv/1000)!=indicator))){
                indicator = counv/1000;
                of <<"Patterns explored: "<< counv <<" closed motifs found: "<<container.numSubgraphs()<<std::endl;
            }

            //We try to get the next node.
            boost::tie(frequent_children,next_node) = current_node.get_next(kt,freq,of);
            bool inserted = false;
            
            //No child are possible
            counter++;
            if(!frequent_children){
                if(is_root){
                	if(sizeMin<=1){
                        counter++;
                    	container.insert_closed_pattern(current_node,inserted,of,insert_error);
                	}
                    break;
                }else{

                    //In this case we check that the pattern is complete.
                    //No incoherence found at this step. This ensure that no value
                    //Upper have been found directly in the tree.
                    bool inserted;
                    if(current_node.numOccs() > max_child_occs.top()){

                        //The pattern is inserted in the data structure
                        if(current_node.isCompleteD(mgs)){
                            current_node.reconstructGraphD(mgs);
                            graphp &gg = current_node.get_g();
                            graphTraitsp::adjacency_iterator ba,ea;
                            graphTraitsp::out_edge_iterator bo,eo;
                            boost::tie(ba,ea)=boost::adjacent_vertices(current_node.get_root(),gg);
                            boost::tie(bo,eo)=boost::out_edges(current_node.get_root(),gg);
                            counter++;
                            if(counter/1000!=current_1000){
                                current_1000 = counter/1000;
                                Rcpp::Rcerr << current_1000*1000 <<" ";
                            }
//                            Rcpp::Rcerr << " out " << std::endl;
                            container.insert_closed_pattern(current_node,inserted,of,insert_error);
                        }
                    }
                    max_child_occs.pop();
                    lattice_stack.pop();
                }
            }else{ //We continue
                //In this case we need to check that there is not a subgraph or a supergraph with the same values.
                counv ++;


                //We update the maxChild value before adding an element
                int num_occs_child = next_node.numOccs();
                short& new_occs = max_child_occs.top();
                if(num_occs_child<new_occs){
                    new_occs = num_occs_child;
                }

                //The node is added to both stack.
                lattice_stack.push(next_node);
                max_child_occs.push(0);
            }
        }
    }
    Rcpp::Rcerr<<"Explored "<< counv<<" graphs. "<<container.numSubgraphs()<<" closed frequent subgraphs found."<<std::endl;

    //We add the information on the found patterns
    container.postProcessing(num_graph,of);
}


//Getting the subgraph container.
subgraph_container& mass_dag_miner::get_container(){
    //returning the subgraph container.
    return container;
}
