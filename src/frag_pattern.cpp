#include <tuple>
#include <iostream>

#include "frag_pattern.h"
#include "k_path_tree.h"
#include "common.h"
#include "triangles_list.h"

#include <boost/graph/copy.hpp>
#include <boost/graph/graphml.hpp>
#include<Rcpp.h>



//Reminder of the typedef
//typedef boost::adjacency_list<
//    boost::vecS, boost::vecS, boost::bidirectionalS,
//    node_info, edge_info> pgraph;
//typedef boost::graph_traits<pgraph> graphTraitp;
//typedef boost::graph_traits<pgraph>::vertex_descriptor Vertexp;


frag_pattern::frag_pattern()
{
    //ctor
};



std::vector<Extension>::iterator pos_extension(std::vector<Extension>& exts, Extension ext){
//We fidn the position ot remove

    auto it = exts.begin();
    auto ie = exts.end();

    for(;it!=ie;it++){
        if((std::get<0> (*it)==std::get<0>(ext)) &
           (std::get<2> (*it)==std::get<2>(ext))){
                break;
        }
    }
    //We go to the next value
    return it;
}


//short calcDist(Extension ext,graphp& g, triangles_list& tl){
//    short odist = g[(std::get<0>(ext))].lab;
//    short nlab = std::get<2>(ext);
//    return tl.get_ab(odist,nlab);
//}


void frag_pattern::filterExtensions(triangles_list& lt){
    std::vector<Extension> nexts;
    std::copy_if(extensions.begin(),extensions.end(),nexts.begin(),
    [this,lt](Extension e)-> bool{
            short odist = g[(std::get<0>(e))].lab;
        short nlab = std::get<2>(e);
        short dist = lt.get_ab(odist,nlab);

        return(dist_prec.find(dist)== dist_prec.end());
    });
}


//We chekc if the extnesion is possible without instantiating an object.
bool frag_pattern::validExtension(Extension& ext,triangles_list& tl){
    graphp parent_graph = get_g();
    Vertexp source_vertex= std::get<0>(ext);
    short new_lab = std::get<2>(ext);
    short dprec_source = parent_graph[source_vertex].lab;

    //We always check that the extnesions is not a forbidden dist.
    //Calculation of the new new dist
    short new_dist = tl.get_ab(dprec_source,new_lab);
    auto itdist = dist_prec.find(new_dist);
    if(itdist!=dist_prec.end()){
        return false;
    }
    return true;
}

//Ad the distance from the exts.begin() to the iterator
void add_dist(graphp& g,std::vector<Extension>& exts,
             std::vector<Extension>::iterator& endpos,std::unordered_set<short>& dist, triangles_list& tl){
    //Case where it is directly the root
    std::vector<short> forbidden_dist;
    auto ie = exts.begin();
    for(;ie!=endpos;ie++){
        short odist = g[std::get<0>(*ie)].lab;
        if(odist==0){
            dist.insert(std::get<2>(*ie));
        }else{
            short fdist = tl.get_ab(odist,std::get<2>(*ie));
            //std::cout << "fdi : "<<odist<<" "<<std::get<2>(*ie)<<" "<<fdist<<std::endl;
            dist.insert(fdist);
        }
    }
}

//Filter the extensions which are already in dist_prec.
void rfilterExtensions(graphp gp,std::vector<Extension>& exts,std::unordered_set<short>& dist,triangles_list& lt){
    exts.erase(std::remove_if(exts.begin(),
              exts.end(),
              [lt,dist,gp](Extension e)->bool{
                short ola = gp[std::get<0>(e)].lab;
                short la = std::get<2>(e);
                short pc = lt.get_ab(ola,la);
                if((pc==NULL_LABEL)|(dist.find(pc)!=dist.end())){
                    return true;
                }
                return false;
            }),
   exts.end());
}


void removeExtensions(std::vector<Extension>& exts, std::vector<short> labs, Vertexp origin, std::ostream& of){
//We loop through the list until we find the correct precursor.
    auto idx = exts.begin();

    //We find the begining in the extension list.
    while( idx != exts.end() && std::get<0>(*idx) != origin) idx++;
    if(idx==exts.end()) return;

    auto idx_end = idx;
    while( idx_end != exts.end() && std::get<0>(*idx_end) == origin) idx_end++;

    //We erase the value if they are present in the vector.
    auto past_the_end_it =  std::remove_if(idx,
                              idx_end,
                              [labs,origin](Extension e)-> bool{
                                    bool temp = std::find(labs.begin(),labs.end(),std::get<2>(e)) != labs.end();
//                                    if(temp) std::cout << "rr "<<origin <<"_"<<*std::find(labs.begin(),labs.end(),std::get<2>(e));
                                  return(temp);});
    //Countting the number of removeed extensions/
    int counter= 0; //DEBUG
    for(auto it = past_the_end_it; it != exts.end(); it++){ //DEBUG
        counter++; //DEBUG
    } //DEBUG

    exts.erase(past_the_end_it,idx_end);
}


//This construtor only produce pattenr with 1 edge or seed pattern.
frag_pattern::frag_pattern(Vertext v,k_path_tree& kt, std::ostream& of){
    //Setting occurences, occurences may be added :
    ktree& t = kt.get_t();
    MapOccurrences& moccs = kt.get_occs();
    triangles_list& tl = kt.get_tl();
    adjacencyGraph& adj = kt.get_adj();
    std::vector<occ> toccs(moccs[v].begin(),moccs[v].end());
    occurences = toccs;
    short lab = t[v].lab;
    short k = kt.get_k();

    //For indexing purpose at the moment.
    key = lab;

    //The graph g is initalised with a single edge.
    root = boost::add_vertex(g);
    g[root].depth=0;
    g[root].lab=0;
    Vertexp n1 = boost::add_vertex(g);
    g[n1].depth=1;
    g[n1].lab=lab;

    graphTraitsp::edge_descriptor e =boost::add_edge(root,n1,g).first;
    g[e].lab = lab;

    //The added edge is added to dist_prec
    dist_prec.insert(lab);

    //The extension list is initialized with the possible edges with single values
    auto pexts = adj.neighbours(lab);

    //The forbidden values are updated given the removed extension.
    //We create the list of extension
    std::transform(pexts.begin(),pexts.end(),std::back_inserter(extensions),
                   [this](std::pair <Vertext,short> pv)-> Extension {
                        return std::make_tuple(this->root,pv.first,pv.second);
                   });
    //We sort the list of extensions.
    std::sort(extensions.begin(),extensions.end(),
    [](const Extension & a, const Extension & b) -> bool
    {
    return ((std::get<2>(a))<(std::get<2>(b)));
    });
    //We remove the extension with a lower label.

    //an iterator is defined on the selected position.
    std::vector<Extension>::iterator limExts = extensions.begin();


    while((limExts != extensions.end())&&(lab > std::get<2>(*limExts))){
        limExts++;
    }
    //Now the removed distance are added to the precursor dist extensions.
    add_dist(g,extensions,limExts,dist_prec,tl);

    auto vfirst = limExts;
    auto vlast = extensions.end();
    std::vector<Extension> vext_origin(vfirst,vlast);
    extensions = vext_origin;

    //Now a vector of new extensions is added to the value.
    graphTraitst::adjacency_iterator bo,be;
    boost::tie(bo,be)=boost::adjacent_vertices(v,t);

    //Now we add all the extension if possible.
    std::vector<Extension> temp_exts;
    std::transform(bo,be,std::back_inserter(temp_exts),
               [n1,t](Vertext v)->Extension{
                   Extension t_ext = std::make_tuple(n1,v,t[v].lab);
                   return t_ext;
                   });

    //If k is big enough we can add the new value.
    std::vector<Extension> exts_n1 = kt.getExtensions(n1,v);

    if(REMOVE_EXTENSION){
        rfilterExtensions(g,exts_n1,dist_prec,tl);
    }

    //Now if we are nto at the last level of the arborescnece
    //we remove the extensions. We remove the extension.
    if(k>1){
        //We get all the C value for the added distance
        //We get all the extension which could be attained from the new vertex..
        std::vector<short> pc = tl.get_a_c(lab);

        //Now remove all the extensions which include a c in the predecessor node.
        removeExtensions(extensions, pc, root, of);
    }

    //We create the definitive set of extension
    extensions.insert( extensions.begin(), exts_n1.begin(), exts_n1.end() );
    calcNormalForm();
}

void printOccs(std::set<occ>& occurences){
    for(auto it=occurences.begin();it!=occurences.end();it++){
        std::cout << (*it).gid << " " << (*it).idx << std::endl;
    }
}

int frag_pattern::numExts(){return extensions.size();}
int frag_pattern::numOccs(){return occurences.size();}
int frag_pattern::sizeGraph(){return boost::num_vertices(g);}

//TODO see if it is more valuable to store the distance form the root in the patter.
//Extend a pattern given an extension (the extension is considered correct)
frag_pattern::frag_pattern(frag_pattern& fp,k_path_tree& kt,
                           int iext, int fthreshold,
                           bool& created,std::ostream& of){


    //Initialisation of the used variable
    Extension ext = fp.extensions[iext];
    Vertexp source_vertex= std::get<0>(ext);
    Vertext tnode = std::get<1>(ext);
    short lab = std::get<2>(ext);
    root = fp.root;

    //Forbbiden first
    graphp& parent_graph = fp.get_g();
    triangles_list& tl = kt.get_tl();
    auto k = kt.get_k();

    //We then check if the extension is syntaxically correct.
    short new_lab = std::get<2>(ext);
    short dprec_source = parent_graph[source_vertex].lab;
    short depth_source = parent_graph[source_vertex].depth;
    //Specific care to 0 extension
    short new_dist;
    if(depth_source==0){
        new_dist = lab;
    }else{
        new_dist = tl.get_ab(dprec_source,new_lab);
    }
    //We check if the extension is not in the forbidden dist
    if((fp.dist_prec.find(new_dist) != fp.dist_prec.end())|(new_dist==NULL_LABEL)){
        created = false;
        return;
    }

    //We first check the number of occurences
    MapOccurrences& moccs = kt.get_occs();
    std::vector<occ>noccs;

    set_intersection(moccs[tnode].begin(),moccs[tnode].end(),fp.occurences.begin(),fp.occurences.end(),
                  std::inserter(noccs,noccs.begin()));

    //Creation only if the number of occurences is sufficient.
    if(int(noccs.size())< fthreshold){
        created = false;
        return;
    }else{
        occurences = noccs;
    }


    key = fp.key;

    //Copying of the other elements.
    extensions = fp.extensions;
    dist_prec = fp.dist_prec;

    //Copy of the graph.
    boost::copy_graph(parent_graph,g);

    //We add the vertiex and the corresponding edge.
    Vertexp current_vertex =  boost::add_vertex(g);
    g[current_vertex].depth = (g[source_vertex].depth)+1;
    g[current_vertex].lab = new_dist;
    pedge_info info_new_edge = {true,new_lab};
    boost::add_edge(source_vertex,current_vertex,info_new_edge,g);
    dist_prec.insert(new_dist);


    //We handle the new extension list
    //We start by removing the extension lower than the current ones.
    auto lim_ext = extensions.begin() + iext+1;

    //Now the removed exntensions aree incorporated in dist_prec.
    add_dist(g,extensions,lim_ext,dist_prec,tl);
    //std::cout << "add_dist ";
    auto vfirst = lim_ext;
    auto vlast = extensions.end();
    std::vector<Extension> old_exts(vfirst,vlast);
    //extnesions is copied, see if there is not any better ways.
    extensions = old_exts;

    //We add the extensions originating from this new vertex
    //extracted form the k path tree.
    ktree& t = kt.get_t();
    graphTraitst::adjacency_iterator bo,be;
    boost::tie(bo,be)=boost::adjacent_vertices(tnode,t);
    //Now we add all the extension if possible.
    std::vector<Extension> temp_exts;
    std::transform(bo,be,std::back_inserter(temp_exts),
               [current_vertex,t](Vertext v)->Extension{
                   Extension t_ext = std::make_tuple(current_vertex,v,t[v].lab);
                   return t_ext;
                   });
    std::sort(temp_exts.begin(),temp_exts.end(),less_than());


    //Now we add the correct subvector
    if(k>g[current_vertex].depth){
        //We get all the extension which could be attained from the new vertex.
        std::vector<short> pc = tl.get_a_c(lab);

        //Now remove all the extensions which include a c in the predecessor node.
        removeExtensions(extensions, pc, source_vertex, of);
    }
    extensions.insert(extensions.begin(),temp_exts.begin(),temp_exts.end());
    created = true;
    calcNormalForm();
}

//Return a pair of iterator for the extension
std::pair<std::vector<Extension>::iterator,
            std::vector<Extension>::iterator> frag_pattern::get_it_exts(){
            return std::make_pair(extensions.begin(),extensions.end());
}


//Update the distance preducrosr given an extension.
std::vector<Extension>::iterator frag_pattern::update_dist_prec(std::vector<Extension>& pexts,
                                                                Extension ext,triangles_list& lt){
    //We find the position of the extension vector.
    if(pexts.size()==0){
        return pexts.begin();
    }
    //We get the label of the current node
    short vlab =  std::get<2>(ext);

    for(auto bext = pexts.begin();bext!=pexts.end();bext++){
        short ext_lab = g[std::get<0>(*bext)].lab;
        auto cval = lt.get_ab(vlab,ext_lab);

        //The C is added to the value.
        if(cval!=NULL_LABEL){
            dist_prec.insert(cval);
        }
    }

    return pexts.end();
}




frag_pattern::~frag_pattern()
{
    //dtor
}

graphp& frag_pattern::get_g(){
    return g;
}

patternKey frag_pattern::get_key() const{
    return key;
}

std::vector<occ>& frag_pattern::get_occs(){
    return occurences;
}

void frag_pattern::clearPattern(){
    extensions.clear();
    dist_prec.clear();

}


void frag_pattern::calcNormalForm(){
    //The key is calculated
    std::stringstream ss;
    graphTraitsp::out_edge_iterator bo,eo;
    boost::tie(bo,eo) = boost::out_edges(root,g);

    bool sep = false;
    for(; bo != eo ; bo++)
    {
      if(sep){
        ss << ".";
      }else{
        sep = true;
      }
      ss << g[*bo].lab;
    }
    norm =  ss.str();
}

std::string& frag_pattern::get_norm(){
    return norm;
}



Vertexp frag_pattern::get_root(){
    return root;
}

void frag_pattern::extendMotif(){

}


//This function reconstruct the full graph, given the triangle list.
//If las is true only the last element is che
void frag_pattern::constructFullGraph(triangles_list& tl){
    //For every pair of node we try to reconstruct a value.
    graphTraitsp::vertex_iterator bv,ev;
    boost::tie(bv,ev)=boost::vertices(g);
    std::vector<Vertexp> vexp(bv,ev);

    //At first an edge is added between every vertex and the precursor
    for(;bv!=ev;bv++){
        if(*bv==root) continue;
        auto eve = boost::add_edge(root,*bv,g);
        g[eve.first].lab = g[*bv].lab;
    }

    //The vector is sorted in increasing order.
    std::sort(vexp.begin(),vexp.end(),
    [this](const Vertexp & a, const Vertexp & b) -> bool
    {
        return ((this->g[a].lab)<(this->g[b].lab));
    });

    //Now we reconstruct the graph
    std::vector<Vertexp>::iterator b1,b2,e2;
    for(b1=vexp.begin();b1!=vexp.end();b1++){
        //This a
        short a = g[*b1].lab;
        for(b2=b1,b2++;b2!=vexp.end();b2++){
            short c = g[*b2].lab;
            bool inserted=false;
            graphTraitsp::edge_descriptor e;
            short b = tl.get_ac(a,c);
            if(b==NULL_LABEL) continue;

                boost::tie(e,inserted) = boost::add_edge(*b1,*b2,g);
                if(inserted){
                    g[e].lab = b;
                }
        }
    }
}

//This function reconstruct the full graph, given the triangle list.
//it only consider the edge which may be rattached to Vertexp v.
void frag_pattern::constructFullGraph(triangles_list& tl,Vertexp v){
    //For every pair of node we try to reconstruct a value.
    graphTraitsp::vertex_iterator bv,ev;
    boost::tie(bv,ev)=boost::vertices(g);
    short labv = g[v].lab;

    for(;bv!=ev;bv++){
        if((*bv)==(v)){
            continue;
        }
        short labo = g[*bv].lab;
        graphTraitsp::edge_descriptor e;
        bool inserted;

        if(labo>labv){
            //We insert the edge if needed.
            //We try to insert if it does not exist.
            //auto ex_edge = boost::edge(v,*bv,g);
            //if(!ex_edge.second){
                short b = tl.get_ac(labv,labo);
                boost::tie(e,inserted) = boost::add_edge(v,*bv,g);
                if(inserted){
                    g[e].lab = b;
                }
            //}
        }else{
            //auto ex_edge = boost::edge(v,*bv,g);
            //if(!ex_edge.second){
                short b = tl.get_ac(labo,labv);
                boost::tie(e,inserted) = boost::add_edge(*bv,v,g);
                if(inserted){
                    g[e].lab = b;
                }
            //}
        }
    }
}



/////DEBUGGUING FUNCTION
//Return true if the pattern is coherent.
bool frag_pattern::isCoherent(){
    //For each vertices we check that there is not an extension with a lower label
    graphTraitsp::vertex_iterator vb,ve;

    boost::tie(vb,ve) = boost::vertices(g);

    //We check if for very vertex his extnesion are not good.
    for(;vb!=ve;vb++){
        //We get the lab of the out extensions if there is one.
        graphTraitsp::out_edge_iterator vba,vea;
        boost::tie(vba,vea) = boost::out_edges((*vb),g);
        if(vba==vea) continue;
        short min_lab = 10000;
        for(;vba!=vea;vba++){
            short clab = g[*vba].lab;
            if(clab<min_lab){
                min_lab = clab;
            }
        }

        if(min_lab!=10000){
            for(auto it=extensions.begin();it!=extensions.end();it++){
                if((std::get<0>(*it)==(*vb))&&(std::get<2>(*it)<min_lab)){
                    return false;
                }
            }
        }
    }
    return true;
}


void frag_pattern::to_string(){
    std::cout<<"g : vertices "<< boost::num_vertices(g)<<" edges " <<
    boost::num_edges(g) <<" coherent "<<this->isCoherent()<< std::endl;
    graphTraitsp::vertex_iterator bv,ev;
    for(boost::tie(bv,ev)=boost::vertices(g);bv!=ev;bv++){
        std::cout<<g[*bv].lab<<"_"<<g[*bv].depth<<" ";
    }
    std::cout<<std::endl;

    graphTraitsp::edge_iterator bee,eee;
    for(boost::tie(bee,eee)=boost::edges(g);bee!=eee;bee++){
        std::cout<<boost::source(*bee,g)<<"|"<<g[*bee].lab<<"|"<<boost::target(*bee,g)<<" ";
    }
    std::cout<<std::endl;


    //Now we print the occurences
    auto bo=occurences.begin();
    auto eo=occurences.end();
    std::cout<<"occs : ";
    for(;bo!=eo;bo++){
        std::cout<<(*bo).gid<<"_"<<(*bo).idx<<" ";
    }
    std::cout<<std::endl;

    //Extensions
    std::vector<Extension>::iterator be,ee;
    boost::tie(be,ee) = get_it_exts();
    std::cout <<"its : "<<&(*be)<<std::endl;
    std::cout<<"exts : ";
    for(;be!=ee;be++){
        std::cout<<std::get<0>(*be)<<"_"<<std::get<1>(*be)<<"_"<<std::get<2>(*be)<<" ";
    }
    std::cout<<std::endl;

    //Prec _dist plot
    auto bd=dist_prec.begin();
    auto ed=dist_prec.end();
    std::cout<<"dist_prec : ";
    for(;bd!=ed;bd++){
        std::cout<<(*bd)<<" ";
    }
    std::cout<<std::endl;
}


void frag_pattern::to_reduced_string(){
    std::cout << "## ";
    graphTraitsp::vertex_iterator bv,ev;
    for(boost::tie(bv,ev)=boost::vertices(g);bv!=ev;bv++){
        std::cout<<g[*bv].lab<<"_"<<g[*bv].depth<<" ";
    }
    std::cout<<std::endl;
}

void frag_pattern::to_reduced_string(std::ostream& of){
    of << "## ";
    graphTraitsp::vertex_iterator bv,ev;
//    for(boost::tie(bv,ev)=boost::vertices(g);bv!=ev;bv++){
//        of<<g[*bv].lab<<"_"<<g[*bv].depth<<"_"<<g[*bv].lab<<" ";
//    }
    boost::tie(bv,ev)=boost::vertices(g);
    //Debug only
    std::vector< short> temp;
    std::transform(bv,ev,std::back_inserter(temp),
                   [this](Vertex v)->short{
                    return this->g[v].lab;
                   });
    std::sort(temp.begin(),temp.end());
    for(auto it = temp.begin(); it != temp.end();it++){
        of << *it <<" " ;
    }

//
//    //FOR DEBUGGING ONLY.
//    of << "##E ";
//    graphTraitsp::edge_iterator be,ee;
//    for(boost::tie(be,ee)=boost::edges(g);be!=ee;be++){
//        of<<g[*be].lab<<"_"<<g[boost::source(*be,g)].lab<<"_"<<g[boost::target(*be,g)].lab<<" ";
//    }
//    of<<std::endl;
//
//    auto bd=dist_prec.begin();
//    auto ed=dist_prec.end();
//    of<<"dist_prec : ";
//    for(;bd!=ed;bd++){
//        of<<(*bd)<<" ";
//    }
//    of<<std::endl;
//
//        of<<"exts : ";
//    std::vector<Extension>::iterator bex,eex;
//    boost::tie(bex,eex) = get_it_exts();
//    for(;bex!=eex;bex++){
//        of<<std::get<0>(*bex)<<"_"<<std::get<1>(*bex)<<"_"<<std::get<2>(*bex)<<" ";
//    }
//    of<<std::endl;

}

void frag_pattern::clearExts(){
    extensions.clear();
}

void frag_pattern::clearPatternFull(){
    graphp temp;
    g = temp;
    extensions.clear();
    occurences.clear();
    dist_prec.clear();
}


void frag_pattern::to_reduced_string_ext(std::ostream& of){
//    std::cout<<"g : vertices "<< boost::num_vertices(g)<<" edges " <<
//    boost::num_edges(g) <<" coherent "<<std::endl;//<<this->isCoherent()<< std::endl;
    of << "## ";
    graphTraitsp::vertex_iterator bv,ev;
    for(boost::tie(bv,ev)=boost::vertices(g);bv!=ev;bv++){
        of<<g[*bv].lab<<"_"<<g[*bv].depth<<" ";
    }
    of<<std::endl<<"edges";
    //Edges print
    graphTraitsp::edge_iterator be,ee;
    for(boost::tie(be,ee)=boost::edges(g);be!=ee;be++){
        of<<" "<<g[boost::source(*be,g)].lab<<"_" <<g[*be].lab<<"_"<<g[boost::target(*be,g)].lab<<" ";
    }


    of<<std::endl;
}


//Return the number of graph id in whic a pattern
//is found.
std::set<short> frag_pattern::numUniqueOccs(){

    std::set<short> unique_gid;
    for(auto it = occurences.begin();it!=occurences.end();it++){
        unique_gid.insert(it->gid);
    }
    return unique_gid;

}

//Fill the information of the graph.
void frag_pattern::fill_graph_info(){

    //Writing the occurences
    std::ostringstream tst;


    for(auto it = occurences.begin();it!=occurences.end();it++){
        tst <<it->gid <<"_" << it->idx << "|";
    }


    g[boost::graph_bundle].occs = tst.str();
    g[boost::graph_bundle].unique_occs = numUniqueOccs().size();

    //Writing the normal form.
    g[boost::graph_bundle].norm = norm;

}


//DIRTY FUNCTION TO EXPORT THE LATTICE AND THE PATTERNS.
void frag_pattern::write_graphml(std::string path_graphml)
{
    boost::dynamic_properties dp(boost::ignore_other_properties);
    fill_graph_info();
    //Node properties are defined.
    dp.property("dist_prec", boost::get(&pnode_info::lab, g));

    //Edge properties.
    dp.property("lab", boost::get(&pedge_info::lab, g));
    dp.property("span", boost::get(&pedge_info::span, g));
    //dp.property("occs", boost::get(&pgraph_info::occs, g));
    //dp.property("norm", boost::get(&pgraph_info::norm, g));
    //
    boost::ref_property_map<graphp*,std::string> gnorm(
    get_property(g,&pgraph_info::norm));
    dp.property("norm",gnorm);

    boost::ref_property_map<graphp*,std::string> goccs(
    get_property(g,&pgraph_info::occs));
    dp.property("occs",goccs);

    //std::cout << "Writing..." << path_graphml;
    //std::cout << "nv:"<<boost::num_vertices(g)<< "ne:"<<boost::num_edges(g)<<std::endl;
    //boost::write_graphml(std::cout, g, dp, true);
    std::ofstream ofile(path_graphml);
    boost::write_graphml(ofile, g, dp, true);
    ofile.close();
    //std::cout << "written...";
}


//Return the pattern as an edge matrix and an occurences matrix.
Rcpp::List frag_pattern::as_igraph_data_frame(){

	int nedges = boost::num_edges(g);

	//Edges matrix
	Rcpp::IntegerMatrix mat_edges(nedges,3);
	Rcpp::colnames(mat_edges) = Rcpp::CharacterVector({"from","to","lab"});

	//Filling the edge data.frame.
	graphTraitsp::edge_iterator be,ee;
	boost::tie(be,ee)=boost::edges(g);
	int count_edge = 0;
	for(;be != ee; be++){
		//Vertex descriptors are integer because of vector choise.
		mat_edges(count_edge,0) = int(boost::source(*be,g));
		mat_edges(count_edge,1) = int(boost::target(*be,g));
		mat_edges(count_edge,2) = int(g[*be].lab);
		count_edge++;
	}

	//Filling the occurences matrix
	Rcpp::IntegerMatrix mat_occs(occurences.size(),2);
	Rcpp::colnames(mat_occs) = Rcpp::CharacterVector({"gid","idx"});
	int posmat = 0;
	for(auto it = occurences.begin(); it != occurences.end(); it++){
		mat_occs(posmat,0) = (*it).gid+1;
		mat_occs(posmat,1) = (*it).idx;
		posmat++;
	}

	Rcpp::List res = Rcpp::List::create(Rcpp::Named("edges")=mat_edges,
                               Rcpp::Named("occurences")=mat_occs);
	return(res);
}


