#include <iostream>
#include <tuple>
#include <R_ext/Random.h>

#include "subgraph_container.h"


subgraph_container::subgraph_container()
{
	num_patterns=0;
    //ctor
}

subgraph_container::~subgraph_container()
{
    //dtor
}



//Method to get an element.
frag_pattern& subgraph_container::get(patternIdx idx){
    return (pmap[idx.first])[idx.second];
}

frag_pattern& subgraph_container::get(patternKey k,short i){
    return (pmap[k])[i];
}

std::tuple< std::vector<frag_pattern>::iterator ,std::vector<frag_pattern>::iterator ,bool > subgraph_container::find_key_it(patternKey key){

    //We first try to find if the key exist.
    auto it = pmap.find(key);

    //If the pattern is found
    if(it!=pmap.end()){
        return std::make_tuple((*it).second.begin(),(*it).second.end(),true);
    }else{
        return(std::make_tuple(((*(pmap.begin())).second).end(),
                              ((*(pmap.begin())).second).end(),false));
    }
}


//Inserting a pattern in the data structure.
void subgraph_container::insert_pattern(frag_pattern& pat,std::ostream& of){
    //We get the key
    patternKey key = pat.get_key();
    auto it = pmap.find(key);

    //If the pattern is found
    if(it!=pmap.end()){
        //We insert the pattern to the righ position
        (it->second).push_back(pat);
        ((it->second).back()).clearPattern();
    }else{ //If it is not found we add the vector
        std::vector<frag_pattern> temp;
        temp.push_back(pat);
        (temp[0]).clearPattern();
        pmap.insert(std::make_pair(key,temp));
    }

}

//
void subgraph_container::insert_closed_pattern(frag_pattern& fp, bool& inserted, std::ostream& of, int& num_error){
    //The normal form is reconstructed.
    fp.calcNormalForm();

    //We check if the pattenr does not already exist, if it does, it is a labelling mistake and the
    //pattern is not inserted
    std::string & norm = fp.get_norm();
    if(norm_set.find(norm)!=norm_set.end()){
        num_error++;
        inserted = false;
        return;
    }

    std::vector<patternIdx> to_rm;
    std::vector<std::string> to_rm_norm;
    std::string current_key = fp.get_norm();
    insert_pattern(fp,of);
    inserted = true;
}



//Finding the idx of subgraph of a graph this list of score is given as.
std::vector<patternIdx> subgraph_container::findSubgraph(frag_pattern& pat){
    //We start by getting all possible subgraph from the pattern
    std::vector<Edgep> econt = pat.enumerate_possible_subpattern();
    auto & g = pat.get_g();

    //Results
    std::vector<patternIdx> res;

    //For each found possible subpatterns evaluate if there is a common subpattern.
    for(auto it = econt.begin();it!=econt.end();it++){

        //The lab is the key.
        short lab = g[*it].lab;

        //We get all the 0 exts val.
        auto eset = pat.get_out_edge_labs(*it);

        //Pattern iterator.
        std::vector<frag_pattern>::iterator bp,ep;
        bool found = false;
        auto val =  find_key_it(lab);
        bp = std::get<0>(val);
        ep = std::get<1>(val);
        found = std::get<2>(val);

        //We check if the pattern exists.
        if(!found){
            continue;
        }else{
            int idx = 0;
            for(;bp!=ep;bp++){
                if((*bp).null) continue;
                if(is_isomorphic(*bp,eset)){
                    res.push_back(std::make_pair(lab,idx));
                }
                idx++;
            }
        }
    }
    return res;
}

//We get the index to test for isomorphism.
std::vector<patternIdx> subgraph_container::getInitialIdxSupermotifs(frag_pattern& pat,std::ostream& of,int nrand){
    std::vector<patternIdx> res;
    graphp& g = pat.get_g();
    int nedges = boost::num_edges(g);
    if(nedges<1){
        return res;
    }
    int num = 0;
    //The two first vertices are always useless.
    graphTraitsp::vertex_iterator bv,ev;

    //Modificiation picking some random edges.
    graphTraitsp::edge_iterator be,ee;
    boost::tie(be,ee) = boost::edges(g);

    GetRNGstate();
    while(num<nrand){
        //We get a random summit
        auto rit = std::next(be, lround(R_unif_index(1e9))%(nedges));
            //We pick the first edge
            short elab = g[*rit].lab;

            //We get the current idx
            auto idx = lab_idx.find(elab);

            if(idx == lab_idx.end()) return res;

            //Initial value
            if(res.size()==0){
                res = idx->second;
            }else{
                std::vector<patternIdx> temp;
                std::set_intersection(res.begin(),res.end(),idx->second.begin(),idx->second.end(),
                                      std::back_inserter(temp));
                res = temp;
            }
            if(res.size()==0){
                return res;
            }
        num++;
    }
    PutRNGstate();
    return res;
}

//Function which find the supergraph of a graph.
std::vector<patternIdx> subgraph_container::findSupergraph(frag_pattern& pat,std::ostream& of){
    //TODO implement an index vector to reduce the number of ismoprhism.
    std::vector<patternIdx> to_return;
    //Necessity to reduce the possible subset of this graph.
    auto pid = getInitialIdxSupermotifs(pat,of);

    for(auto it = pid.begin();it!=pid.end();it++){
        if(is_subgraph(pat,get(*it))){
            to_return.push_back(*it);
        }
    }
    return to_return;
}



//Return the score associated to each pattern in the idx fiedl.
std::vector<int> subgraph_container::apply_fun(std::vector<patternIdx>& idxs,int (&scorefun)(frag_pattern&)){
    std::vector<int> res(idxs.size());
    size_t i = 0;
    for(auto it=idxs.begin();it<idxs.end();it++){
        res[i] = scorefun(this->get(*it));
        i++;
    }
    return res;
}
//Return the number of subgraphs.
int subgraph_container::numSubgraphs(){
    int accu=0;
    for(auto it = pmap.begin();it != pmap.end();it++){
        accu = accu + it->second.size();
    }
    return accu;
}


int subgraph_container::numSubgraphsLattice(){
    return numPatterns();
}

void subgraph_container::idx_to_string(std::ostream& of){
    for(auto it=lab_idx.begin();it != lab_idx.end(); it++){
        of << it->first <<" : ";
        for(auto itt = it->second.begin(); itt != it->second.end() ; itt++){
            of <<" "<<(*itt).first << "_"<<(*itt).second;
        }
        of << std::endl;
    }
}

int subgraph_container::numPatterns(){
  int num_patt = 0;
    for(auto it = pmap.begin(); it != pmap.end(); it++){
      num_patt += it->second.size();
    }
    return(num_patt);
}

