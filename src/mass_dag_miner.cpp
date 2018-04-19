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
    //os << "Reading the "<< num_graph <<" graphs." << std::endl;
    int counter = 0;
    for(auto it=vmasses.begin();it!=vmasses.end();it++){
		//Rcpp::Rcout << "processing graph " << counter << std::endl;

        //Each graph is added to the k path tree which also store the occurences.
        kt.add_graph(*it,counter,prec_only);
        counter++;
    }

    kt.post_processing();
    //os <<"Preprocessing finished."<<std::endl;
    //kt.to_string();
}

void mass_dag_miner::mineFrequentDag(int freq,std::ostream& os){
    ktree& t = kt.get_t();
    os <<"Filtering k-tree with an initial number of "<<boost::num_vertices(t) <<" nodes."<< std::endl;
    kt.filter_frequent_nodes(freq);
    os <<"Remainings nodes: "<<boost::num_vertices(t) <<" nodes."<< std::endl;

    //We start by generating all the frequent dag.
    os <<"Constructing initial patterns...";
    std::vector<lattice_node> initialPatterns = kt.constructOneEdgeGraphs(os);
    os <<"done !"<<std::endl;
    bool frequent_children = true;
    int counv = 0;
    for(auto it=initialPatterns.begin();it!=initialPatterns.end();it++){

        //We empty the current stack
        while(!lattice_stack.empty()){
            lattice_stack.pop();
        }

        //We initialize the value.
        lattice_stack.push(*it);

        int LIM_ITER = 100000;
        int counter = 0;
        lattice_node next_node;
        //std::cout << "nroot " <<lattice_stack.size() << std::endl;
        while(counter<LIM_ITER){
            //Sleep(1000);
            //We always consider the last elemnt of the stack.
            lattice_node& current_node = lattice_stack.top();
            //std::cout << std::endl<<"old_node addr : "<< &current_node << std::endl;
            //We first check if the current pattern is root.
            bool is_root = current_node.is_root();

            //current_node.to_reduced_string();
            //We print the label of the current node

            //We try to get the next node.
            boost::tie(frequent_children,next_node) = current_node.get_next(kt,freq,os);

            //No child are possible
            counter++;
            if(!frequent_children){
                if(is_root){
                    //std::cout << "root"<<std::endl;
                    break;
                }else{
                    //Removing the last element and going back.
                    //std::cout << "popping "; //DEBUG
                    //std::cout<<"old_top "; //DEBUG
                    //lattice_stack.top().to_reduced_string(); //DEBUG
                    //std::cout<<std::endl; //DEBUG
                    lattice_stack.pop();
                    //std::cout<<"n_top "; //DEBUG
                    //lattice_stack.top().to_reduced_string(); //DEBUG
                   // std::cout<<std::endl; //DEBUG
                }
            }else{ //We continue
                //std::cout << "forward pushing "; //DEBUG
                //next_node.to_reduced_string(); //DEBUG
                //std::cout<<"old_top "; //DEBUG
                //lattice_stack.top().to_reduced_string(); //DEBUG
                counv ++;
                lattice_stack.push(next_node);
                //std::cout<<"n_top "; //DEBUG
                //lattice_stack.top().to_reduced_string(); //DEBUG
                //std::cout<<std::endl; //DEBUG
            }
        }
    }
    os<<"Explored "<< counv<<" graphs."<<std::endl;

}


void mass_dag_miner::mineFrequentCompleteDag(int freq,std::ostream& of){
    ktree& t = kt.get_t();
    //kt.to_string(of);
    of <<"Filtering k-tree with an initial number of "<<boost::num_vertices(t) <<" nodes."<< std::endl;
    //kt.to_string(of);
    kt.filter_frequent_nodes(freq);
    of <<"Remainings nodes: "<<boost::num_vertices(t) <<" nodes."<< std::endl;

        //of<<"streetree  ";
    //kt.to_string(of);
    //of<<"end  ";

    //We start by generating all the frequent dag.
    of <<"Constructing initial patterns...";
    std::vector<lattice_node> initialPatterns = kt.constructOneEdgeGraphs(of);
    of <<"done !"<<std::endl;
    bool frequent_children = true;
    int counv = 0;
    std::stack<short> max_child_occs;

    //Used for the closeness predicates.
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
        int counter = 0;
        lattice_node next_node;
        int indicator = -1;
        while(counter<LIM_ITER){
            //We always consider the last elemnt of the stack.
            lattice_node& current_node = lattice_stack.top();

            //We first check if the current pattern is root.
            bool is_root = current_node.is_root();

            //We print the label of the current node
            of << std::endl;
            current_node.to_reduced_string(of);
            if(((counv%1000)==0)&(((counv/1000)!=indicator))){
                indicator = counv/1000;
                of <<"Patterns explored: "<< counv <<" closed motifs found: "<<container.numSubgraphs()<<" in structure: " << container.numSubgraphsLattice()<<std::endl;
            }

            //We try to get the next node.
            boost::tie(frequent_children,next_node) = current_node.get_next(kt,freq,of);
            bool inserted = false;
            //No child are possible
            counter++;
            if(!frequent_children){
                if(is_root){
                	if(sizeMin<=1){
                    	container.insert_closed_pattern(current_node,inserted,of);
                	}
                    break;
                }else{

                    //In this case we check that the pattern is complete.
                    //No incoherence found at this step. This ensure that no value
                    //Upper have been found directly in the tree.
                    bool inserted;
                    if(current_node.numOccs() > max_child_occs.top()){

                        //The pattern is inserted in the data structure
                        current_node.constructFullGraph(kt.get_tl());
                        container.insert_closed_pattern(current_node,inserted,of);
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
            //of << std::endl;
        }
    }
    of<<"Explored "<< counv<<" graphs. "<<container.numSubgraphs()<<" closed frequent subgraphs found."<<
    " In lattice: "<<boost::num_vertices(container.get_lat())<<std::endl;

    //We add the information on the found patterns
    container.postProcessing(num_graph,of);


}


//All the messaging take place in this function.
void mass_dag_miner::mineFrequentCompleteDag(int freq,std::set<short> vals,std::ostream& of){
    ktree& t = kt.get_t();
    //kt.to_string(of);
    of <<"Filtering k-tree with an initial number of "<<boost::num_vertices(t) <<" nodes."<< std::endl;
    kt.filter_frequent_nodes(freq);
    of <<"Remainings nodes: "<<boost::num_vertices(t) <<" nodes."<< std::endl;
    //kt.to_string(of);
    //We start by generating all the frequent dag.
    //of <<"Constructing initial patterns...";
    std::vector<lattice_node> initialPatterns = kt.constructOneEdgeGraphs(of);
    of <<"done !"<<std::endl;


    //We filter out the possible patterns.
    initialPatterns.erase(std::remove_if(initialPatterns.begin(),initialPatterns.end(),
                   [vals](lattice_node & p)-> bool{
                        return (vals.find(p.get_key())==vals.end());
                   }),initialPatterns.end());
    //of << "Remaining: "<<initialPatterns.size()<<std::endl;


    bool frequent_children = true;
    int counv = 0;
    std::stack<short> max_child_occs;

    //Used for the closeness predicates.
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

        int LIM_ITER = 100000;
        int counter = 0;
        lattice_node next_node;
        int indicator = -1;
        while(counter<LIM_ITER){
            //We always consider the last elemnt of the stack.
            lattice_node& current_node = lattice_stack.top();

            //We first check if the current pattern is root.
            bool is_root = current_node.is_root();

            //We print the label of the current node
            //of << std::endl;
            current_node.to_reduced_string(of);
            if(((counv%1000)==0)&(((counv/1000)!=indicator))){
                indicator = counv/1000;
                //of <<"Patterns explored: "<< counv <<" closed motifs found: "<<container.numSubgraphs()<<" in structure: " << container.numSubgraphsLattice()<<std::endl;
            }

            //We try to get the next node.
            boost::tie(frequent_children,next_node) = current_node.get_next(kt,freq,of);
            bool inserted = false;
            //No child are possible
            counter++;
            if(!frequent_children){
                if(is_root){
                    //We try to add the on edge pattern.
                    container.insert_closed_pattern(current_node,inserted,of);

                    break;
                }else{
                    //In this case we check that the pattern is complete.
                    //No incoherence found at this step. This ensure that no value
                    //Upper have been found directly in the tree.
                    bool inserted;
                    if(current_node.numOccs() > max_child_occs.top()){

                        //The pattern is inserted in the data structure
                        current_node.constructFullGraph(kt.get_tl());
                        container.insert_closed_pattern(current_node,inserted,of);
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
            //of << std::endl;
        }
    }
    of<<"Explored "<< counv<<" graphs. "<<container.numSubgraphs()<<" closed frequent subgraphs found."<<
    " In lattice: "<<boost::num_vertices(container.get_lat())<<std::endl;

    //We add the information on the found patterns
    container.postProcessing(num_graph,of);


}

//Getting the subgraph container.
subgraph_container& mass_dag_miner::get_container(){
    //returning the subgraph container.
    return container;
}
