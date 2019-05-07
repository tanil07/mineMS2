#include <iostream>

#include<boost/graph/adjacency_list.hpp>

#include "frag_pattern.h"
#include "common.h"

#include "Rcpp.h"

//TODO eventually replace it by an iterator.
//Enumerate the possible subpattern of a pattern.
std::vector<Edgep> frag_pattern::enumerate_possible_subpattern(){
    //We generate the seed pattern.
    std::vector<Edgep> candidates;

    //We enumerate the possible subpatterns.
    graphTraitsp::vertex_iterator bv,ev;

    //We loop an all the vertices.
    for(boost::tie(bv,ev)= boost::vertices(g);bv!=ev;bv++){
        if((*bv)==root) continue;

        //For each vertices we get the adjacent list
        graphTraitsp::out_edge_iterator be,ee;
        for(boost::tie(be,ee) = boost::out_edges(*bv,g);be!=ee;be++){

            //The label on the edge is lower than the key pattern.
            if(g[*be].lab<key){
                candidates.push_back(*be);
                break;
            }
        }
    }
    return candidates;
}

//It's not possible to enumerate the possible superpatterns of the graph.


//////ISOMORPHISM FUNCTION

//Return the labels of all the out_edge from the source of e with
//A label over e.
std::vector<short> frag_pattern::get_out_edge_labs(Edgep e){
    std::vector<short> labs;

    short limlab = g[e].lab;
    graphTraitsp::out_edge_iterator be,ee;
    boost::tie(be,ee) = boost::out_edges(source(e,g),g);

    //We check that there is not extension with lower lab which works.
    for(;be!=ee;be++){
        short le = g[*be].lab;
        if(le>=limlab){
            labs.push_back(le);
        }
    }
    std::sort(labs.begin(),labs.end());

    return labs;
}

//Map a pattern to a mass graph
std::map<Vertex,Vertexp> mapPattern(graph &g,Vertex orvertex,graphp &p,Vertexp rootp){

    //First finding the mappin
    //Intermdiate map used to build the mapping
    std::map<short,Vertex> mg;
    graphTraits::out_edge_iterator bo,eo;
    //IndexMap imap = boost::get(boost::vertex_index,g);
    graphTraits::vertex_iterator bv,ev;

    boost::tie(bo,eo)=boost::out_edges(orvertex,g);
    std::transform(bo,eo,std::inserter(mg,mg.begin()),
                   [&g](graphTraits::edge_descriptor u)->std::pair<short,Vertex>{
        return std::make_pair(g[u].lab,boost::target(u,g));
    });

    //Finding the vertex corresponding otthe origin in grpah d
    //TODO stores it eventually
    std::map<Vertex,Vertexp> mapping;
    mapping.insert(std::make_pair(orvertex,rootp));
    graphTraitsp::out_edge_iterator bpo,epo;
    boost::tie(bpo,epo)=boost::out_edges(rootp,p);
//
    //Rcpp::Rcerr << "InMap :"<<std::endl;
    for(;bpo!=epo;bpo++){
        auto pos = mg.find(p[*bpo].lab);
        if(pos!=mg.end()){
            bool ins1;
            Edge e1;
            //Edgep e2;
            boost::tie(e1,ins1)=boost::edge(orvertex,(*pos).second,g);
            //boost::tie(e2,ins1)=boost::edge(rootp,(*pos).second,g);
            //Rcpp::Rcerr << g[e1].lab <<"_" << p[*bpo].lab<<std::endl;
            mapping.insert(std::make_pair((*pos).second,boost::target(*bpo,p)));
        }
    }


//    std::transform(bpo,epo,std::inserter(mapping,mapping.begin()),
//                   [&mg,&p](graphTraitsp::edge_descriptor u)->std::pair<Vertex,Vertexp>{
//        return std::make_pair((*(mg.find(p[u].lab))).second,boost::target(u,p));
//    });
    return mapping;
}



//utility function called by the ismorphism.
std::vector<short> get_labs(frag_pattern& fp,Vertexp v,bool sorted){
    std::vector<short> lda;
    graphp& g = fp.get_g();
    graphTraitsp::out_edge_iterator be,ee;
    boost::tie(be,ee)=boost::out_edges(v,g);
    std::transform(be,ee,std::back_inserter(lda),
                   [g](graphTraitsp::edge_descriptor u)->short{
        return g[u].lab;
    });
    std::sort(lda.begin(),lda.end());
    return lda;
}

//Return the vertexp of b on whic an isomorphism is found.
Vertexp find_subgraph(frag_pattern & a,frag_pattern & b){
    //We get the 0 vertex of a;
    auto a0 = get_labs(a,a.get_root(),true);
    graphTraitsp::vertex_iterator bv,ev;
    boost::tie(bv,ev) = boost::vertices(b.get_g());
    for(;bv!=ev;bv++){
        auto setv = get_labs(a,*bv,true);
        if(std::includes(setv.begin(),setv.end(),a0.begin(),a0.end())) return (*bv);
    }
    return false;
}


//Return true if the graph a is isomorphic to b, 0 otherwise.
bool is_isomorphic(frag_pattern & a,frag_pattern & b){
    auto a0 = get_labs(a,a.get_root(),true);
    auto b0 = get_labs(b,b.get_root(),true);
    return (a0==b0);
}

//Return true if the two graphs are isomorphic.
bool is_isomorphic(frag_pattern & a,std::vector<short> & b){
    auto a0 = get_labs(a,a.get_root(),true);
//    std::cout << "a:"<<std::endl;
//    for(auto it=a0.begin();it!=a0.end();it++) std::cout << *it <<"_";
//    std::cout << std::endl;
//    std::cout << "b:"<<std::endl;
//    for(auto it=b.begin();it!=b.end();it++) std::cout << *it <<"_";
//    std::cout << std::endl;
    return (a0==b);
}

//a is a subgraph of b rooted on the node v. We just have to compare the labels sets.
bool is_subgraph(frag_pattern & a,frag_pattern & b, Vertexp v){
    //We get the 0 vertex of a;
    auto a0 = get_labs(a,a.get_root(),true);
    auto bv = get_labs(a,v,true);
    return std::includes(bv.begin(),bv.end(),a0.begin(),a0.end());
}

//Return true if the set of edge derivate
bool is_subgraph_edge_set(std::vector<short> elabs,frag_pattern & p){
    auto p0 = get_labs(p,p.get_root(),true);

    return p0==elabs;
}


//DEBUG ONLY
template<class T>
void print_vec(T& vec,std::ofstream& of){
    auto e = vec.end();
    of << "vec : ";
    for(auto b=vec.begin();b!=e;b++){
        of << *b << "_";
    }
    of << std::endl;
}


//Return true if a is a vertex of b.
bool is_subgraph(frag_pattern & a,frag_pattern & b){
    //We get the 0 vertex of a;
    auto a0 = get_labs(a,a.get_root(),true);
    graphTraitsp::vertex_iterator bv,ev;
    boost::tie(bv,ev) = boost::vertices(b.get_g());
    for(;bv!=ev;bv++){
        auto setv = get_labs(b,*bv,true);
//        of << "setv: "; //DEBUG
//        print_vec(setv,of);//DEBUG
//        of << "a0: ";//DEBUG
//        print_vec(a0,of);//DEBUG

        bool iso = std::includes(setv.begin(),setv.end(),a0.begin(),a0.end());
        //std::cout << a0.size() <<" "<< setv.size() <<" "<< iso <<std::endl;
        if(iso) return iso;
    }
    return false;
}
