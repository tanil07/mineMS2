#include <iostream>

#include "lattice_node.h"
#include "common.h"
#include "frag_pattern.h"

lattice_node::lattice_node()
{
    //ctor
}

lattice_node::~lattice_node()
{
    //dtor
}


//The 2 standard constructors.

lattice_node::lattice_node(Vertext v,k_path_tree& kt,std::ostream& os): frag_pattern(v,kt,os),
root(true), current_ext(0){

};

//Pattern extension.
lattice_node::lattice_node(lattice_node& fp,k_path_tree& kt,int iext,
            int fthreshold,bool& created,std::ostream& os): frag_pattern(fp,kt,iext,fthreshold,created,os) {
    if(created){
        root = false;
    	current_ext = 0;
    }
}

//This function get the next pattern, except if it is a value.
std::pair<bool,lattice_node> lattice_node::get_next(k_path_tree& kt,int fthreshold,std::ostream& os){
    bool created = false;
    //Now we split the value.
    while((current_ext<numExts()) & (!created)){
        lattice_node pchild(*this,kt,current_ext,fthreshold,created,os);
        current_ext++;
        if(created){
            return std::make_pair(true,pchild);
        }
    }
    clean_up();
    //TODO see if an optimization is possible.
    return std::make_pair(false,lattice_node());
}

//Cleaning up the data.
void lattice_node::clean_up(){

}
