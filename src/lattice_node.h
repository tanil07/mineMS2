#ifndef LATTICE_NODE_H
#define LATTICE_NODE_H


#include "common.h"
#include "frag_pattern.h"

class lattice_node : public frag_pattern{
    public:
        lattice_node();

        //Pattern creation
        lattice_node(Vertext v,k_path_tree& kt,std::ostream& os);
        //Pattern extension.
        lattice_node(lattice_node& fp,k_path_tree& kt,int iext,
                      int fthreshold,bool& created,std::ostream& os);
        //this function return the next pattern to explore and
        //a boolean indicating if we better not go backwards.
        std::pair<bool,lattice_node> get_next(k_path_tree& kt,int fthreshold,std::ostream& os);
        //reinitialize the iterator of position.
        void reinitialize();

        //this function add a child to this pattern.
        //Utility function
        bool is_root(){return root;};
        virtual ~lattice_node();

    protected:
    private:
        //The root anf the current extension
        bool root;
        short current_ext;
        //The number of occurrences of the parent, and of the child.
        //Cleaning
        void clean_up();
};

#endif // LATTICE_NODE_H
