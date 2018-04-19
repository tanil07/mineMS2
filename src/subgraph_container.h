#ifndef SUBGRAPH_CONTAINER_H
#define SUBGRAPH_CONTAINER_H

#include <fstream>

#include <boost/functional/hash.hpp>
#include <Rcpp.h>

#include "common.h"
#include "frag_pattern.h"



//This class store the occurences of all the graph.
class subgraph_container
{
    public:
        //PATTERN STORE STRUCTURE.
        subgraph_container();
        virtual ~subgraph_container();

        //Method to get a pattern.
        frag_pattern& get(patternIdx idx);
        frag_pattern& get(patternKey,short);

        //Inserting a pattern in the data structure.
        void insert_pattern(frag_pattern& pat,std::ostream& of);

        void insert_closed_pattern(frag_pattern& fp, bool& inserted, std::ostream& of);

        //Finding the idx of subgraph of a graph
        std::vector<patternIdx> findSubgraph(frag_pattern& pat);

        //Finding the idx of supergraph of a graph
        std::vector<patternIdx> findSupergraph(frag_pattern& pat,std::ostream& of);

        patternIdx find_pattern_idx(frag_pattern& pat);
        std::pair<std::vector<frag_pattern>::iterator,bool> find_pattern_it(frag_pattern& pat);

        //Remove a set of patterns from the contain
        void removeIdx(std::vector<patternIdx>&);

        //Return the number of patterns.
        int numSubgraphs();
        int numSubgraphsLattice();

        std::tuple< std::vector<frag_pattern>::iterator ,std::vector<frag_pattern>::iterator ,bool > find_key_it(patternKey key);
        std::vector<int> apply_fun(std::vector<patternIdx>& idxs,int (&scorefun)(frag_pattern&));

        //Debug only
        void idx_to_string(std::ostream& of);


        //LATTICE STRUCTURE
        void removePattern_latt(std::string& key);
        void removePattern_latt(std::vector< std::string>& to_rm);

        //Cleaning the lattice
        void cleanLattice();

        //Take as input the key to the precursor,
        //The key to the succesors and perform the operation if necessary.

        void addPattern_latt(std::string& key,std::vector<std::string>& parents,std::vector<std::string>&sons);

        //This post processing step, add the network we a new value.
        void postProcessing(int num_graphs,std::ostream& of);

        void addRoot();
        int numPatterns();
        latticesub& get_lat();


        void printPatternsLL(std::ostream& of);
        void printPatternsLA(std::ostream& of);
        void printPatternsB(std::ostream& of);
        void printPatternsB();

        //Functions used for exporting purpose only
        Rcpp::List exportMinedPatternsRcpp();

        //Functions used only otuside of R.
        // void exportResultingLattice(std::string path_lattice,std::string path_patterns,
        //                                         std::ostream& of);
        // std::map<Vertexl,Vertexel> createExportLattice(latticeexport& elat);
        std::map<Vertexl, int> latticeMapping();

    protected:
    private:
        //PATTERN STORE STRUCTURE.
        //Patterns are indexed by their left-most extensions.
        std::map<patternKey,std::vector<frag_pattern> > pmap;

        int numbers;
        int num_patterns;
        int num_items;

        //Other possible implmentation using vector.
        std::map<short,std::vector<patternIdx> > lab_idx;


        void insertPatternIndex(frag_pattern& pat,patternIdx& idx,std::ostream& of);
        void removePatternIndex(frag_pattern& pat,patternIdx& idx);
        std::vector<patternIdx> getInitialIdxSupermotifs(frag_pattern& pat,std::ostream& of,int nrand=3);


        //Function used ot map the pattern in list ot other patterns.
        std::map< Vertexl, patternIdx> getMapping();

        std::vector<Vertexl> getLeafs();

        //LATTICE STRUCTURE.
        latticesub lat;

        //Map a lattice node to a dataset.
        void resetIdxMap();
        subgraph_idx idxmap;

};

#endif // SUBGRAPH_CONTAINER_H
