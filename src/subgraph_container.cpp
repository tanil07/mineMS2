#include <iostream>
#include <tuple>


#include "subgraph_container.h"


subgraph_container::subgraph_container()
{
    //ctor
}

subgraph_container::~subgraph_container()
{
    //dtor
}



//Method to get an element.
frag_pattern& subgraph_container::get(patternIdx idx){
   // std::cout << "Getting " << idx.first <<" "<<idx.second;
   //if((pmap[idx.first])[idx.second].null)
    //std::cout <<"pm"<<(pmap[idx.first])[idx.second].get_norm()<<"_"<<(pmap[idx.first])[idx.second].null;
    return (pmap[idx.first])[idx.second];
}

frag_pattern& subgraph_container::get(patternKey k,short i){
    return (pmap[k])[i];
}

//Method to find an element.
patternIdx subgraph_container::find_pattern_idx(frag_pattern& pat){
    patternKey key = pat.get_key();

    //We first try to find if the key exist.
    auto it = pmap.find(key);

    //If the pattern is found
    if(it!=pmap.end()){
        for(int i = 0; i < (*it).second.size();i++){
            if((*it).second[i].null) continue;
            if(is_isomorphic(pat,(*it).second[i])){
                return(std::make_pair(key,i));
            }
        }
        return(std::make_pair(NULL_LABEL,0));
    }else{
        return(std::make_pair(NULL_LABEL,0));
    }
}

std::pair<std::vector<frag_pattern>::iterator,bool> subgraph_container::find_pattern_it(frag_pattern& pat){
    patternKey key = pat.get_key();

    //We first try to find if the key exist.
    auto it = pmap.find(key);

    //If the pattern is found
    if(it!=pmap.end()){
        for(auto i = (*it).second.begin(); i != (*it).second.end();i++){
            if((*i).null) continue;
            if(is_isomorphic(pat,(*i))){
                return(std::make_pair(i,true));
            }
        }
        return(std::make_pair(((*(pmap.begin())).second).end(),false));
    }else{
        return(std::make_pair(((*(pmap.begin())).second).end(),false));
    }
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

latticesub& subgraph_container::get_lat(){
    return lat;
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
    }else{ //If it is not found
        std::vector<frag_pattern> temp;
        temp.push_back(pat);
        (temp[0]).clearPattern();
        auto it = pmap.insert(std::make_pair(key,temp));
    }

    patternIdx pv = std::make_pair(key,(pmap[key]).size()-1);

    //In very case the pattern in inserted into the new map.
    insertPatternIndex(pat,pv,of);

}

//
void subgraph_container::insert_closed_pattern(frag_pattern& fp, bool& inserted, std::ostream& of){
    //We check the supergraph.
    //of << " Possible_Insert..";
    std::vector<patternIdx> sup_idx = findSupergraph(fp,of);
    //of<<"num_sump_"<<sup_idx.size();

    std::vector<std::string> supergraphs;
    //Number of occurences is compared for closeness
    short current_score = fp.numOccs();
    for(auto it=sup_idx.begin();it<sup_idx.end();it++){
        if(this->get(*it).numOccs()==current_score){
            //of <<" "<< this->get(*it).numOccs()<<"_vs_"<<current_score<<" ";
            inserted = false;
            return;
        }
        supergraphs.push_back(this->get(*it).get_norm());
    }

    //In this case the node is added and we can calculate the normal form associated.

    //The normal form is reconstructed.
    fp.calcNormalForm();

    //We check if the pattenr does not already exist, if it does, it is a labelling mistake and the
    //pattern is not inserted
    std::string & norm = fp.get_norm();
    if(idxmap.find(norm)!=idxmap.end()){
        inserted = false;
        return;
    }


    //of <<"patt:  "<<fp.get_norm()<<"inserted !! " << std::endl;//DEBUG
    std::vector<patternIdx> to_rm;
    std::vector<std::string> to_rm_norm;
    std::string current_key = fp.get_norm();

    //We check the subgraph. They are removed if they have the same number of occurences.
    //At this step the vector is necessarely inserted.
    int nremoved = 0;
    std::vector<std::string> subgraphs;
    //std::cout <<"fidning_sub ";
    std::vector<patternIdx> sub_idx = findSubgraph(fp);
    //std::cout <<"done! ";
    for(auto it=sub_idx.begin();it<sub_idx.end();it++){
        if(this->get(*it).null) continue;
        if(this->get(*it).numOccs()==current_score){
            to_rm.push_back(*it);
            to_rm_norm.push_back(this->get(*it).get_norm());
            nremoved++;
        }else{
            subgraphs.push_back(this->get(*it).get_norm());
        }
    }
//    of << " removed "<<nremoved; //DEBUG
    //We remove the selected patterns.


    insert_pattern(fp,of);
    addPattern_latt(current_key,subgraphs,supergraphs);
//    if(sub_idx.size()>0){
//    //std::cout << "sup: "<<supergraphs.size() <<"sub: "<< subgraphs.size() <<" to_rm: "<<to_rm_norm.size() <<std::endl;
//    }
    removePattern_latt(to_rm_norm);
    //std::cout << "Removal done";


    removeIdx(to_rm);
    inserted = true;
    //std::cout << "OUT FUN";
}



//Finding the idx of subgraph of a graph this list of score is given as.
std::vector<patternIdx> subgraph_container::findSubgraph(frag_pattern& pat){
    //We start by getting all possible subgraph from the pattern
    std::vector<Edgep> econt = pat.enumerate_possible_subpattern();
    //if(econt.size()>0) std::cout <<"poney"<<std::endl;
    auto & g = pat.get_g();

    //Results
    std::vector<patternIdx> res;

    patternKey key = pat.get_key();

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
                //std::cout << "iso_test";
                if(is_isomorphic(*bp,eset)){
                    //std::cout << "!OK!";
                    res.push_back(std::make_pair(lab,idx));
                }
                //std::cout << std::endl;
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
//    of<<"bsampling: "; //DEBUB
    int nedges = boost::num_edges(g);
    if(nedges<1){
        return res;
    }
//    of<<"sampling: ";

    int num = 0;
    //The two first vertices are always useless.
    graphTraitsp::vertex_iterator bv,ev;

    //Modificiation picking some random edges.
    graphTraitsp::edge_iterator be,ee;
    boost::tie(be,ee) = boost::edges(g);

    while(num<nrand){
        //We get a random summit
        auto rit = std::next(be, std::rand()%(nedges));
        //of << "cnode_samp: "<<g[(*rit)].lab<<"_"<<num<<"/"<<nrand;
        //We pick an edge at random.
        //graphTraitsp::out_edge_iterator bo,eo;
        //boost::tie(bo,eo) = boost::out_edges(*rit,g);
        //Now a random edges is taken
        //if(bo!=eo){
            //We pick the first edge
            short elab = g[*rit].lab;
            //of<< "rlab: "<< elab <<" size: "<< res.size()<< std::endl;
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
        //}
        num++;
    }
    return res;
}

//Function which find the supergraph of a graph.
std::vector<patternIdx> subgraph_container::findSupergraph(frag_pattern& pat,std::ostream& of){
    //TODO implement an index vector to reduce the number of ismoprhism.
    //std::cout << "findSuperGraph";
    std::vector<patternIdx> to_return;

    //Necessity to reduce the possible subset of this graph.
    auto pid = getInitialIdxSupermotifs(pat,of);
    //of<<"idx_"<<pid.size()<<"_"<<std::endl;

//     std::copy_if (pid.begin(), pid.end(), std::back_inserter(to_return),
//                    [pat,this](patternIdx pid){
//
//                   return is_subgraph(pat,this->get(pid));
//                   } );


    for(auto it = pid.begin();it!=pid.end();it++){
        //of << "t" <<it->first<<"_"<<it->second <<"_";
        if(is_subgraph(pat,get(*it))){
            to_return.push_back(*it);
            //of << "iso" << std::endl;
        }
        //of <<" ";
    }
    //std::cout <<"DONE!!";
    return to_return;

//    for(auto it = pmap.begin();it != pmap.end();it++){
//
//    for(auto it = pmap.begin();it != pmap.end();it++){
//        std::vector<frag_pattern>& cvec = it->second;
//        for(int itt=0;itt<cvec.size();itt++){
//            if(is_subgraph(pat,cvec[itt])){
//                to_return.push_back(std::make_pair(it->first,itt));
//            }
//        }
//    }
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


//std::map<short,std::unordered_set<patternIdx> > lab_idx;

void subgraph_container::insertPatternIndex(frag_pattern& pat,patternIdx& idx,std::ostream& of){
    //We insert the pattern in the idx file if necessary.
    graphp& g = pat.get_g();

    if(boost::num_vertices(g)<=2){
        return;
    }
//    of << "idx_inserting ";


    graphTraitsp::out_edge_iterator bo,eo;
    graphTraitsp::vertex_iterator bv,ev;

    //Vertexp root = pat.get_root();
    //Each edge which don't originate form the first compounds is inserted.
    for(boost::tie(bv,ev)=boost::vertices(g);bv != ev;bv++){
        //if((*bv)==root) continue;
//        of << "i";
        //All the out edges are considered
        boost::tie(bo,eo) = boost::out_edges(*bv,g);
        for(;bo!=eo;bo++){

            //the label is extracted
            short elab = (g[(*bo)]).lab;
//            of << elab << " ";

            //We add it to the index.
            (lab_idx[elab]).push_back(idx);

        }
    }
}


void subgraph_container::removePatternIndex(frag_pattern& pat,patternIdx& idx){
    graphp& g = pat.get_g();

    //We insert the pattern in the idx file if necessary.
    if(boost::num_vertices(g)<=2){
        return;
    }


    graphTraitsp::out_edge_iterator bo,eo;
    //Each edge which don't originate form the first compounds is inserted.
    for(int i=0;i<boost::num_vertices(g);i++){
        //All the out edges are considered
        boost::tie(bo,eo) = boost::out_edges(i,g);
        for(;bo!=eo;bo++){
            //the label is extracted
            short elab = g[*bo].lab;
            auto & tmap = lab_idx[elab];

            auto itpos = std::find(tmap.begin(),tmap.end(),idx);
            //We add it to the index.

            tmap.erase(itpos);
        }
    }
}




        //Remove a set of patterns from the container
void subgraph_container::removeIdx(std::vector<patternIdx>& to_remove){
    if(to_remove.size()==0) return;

    //We sort the contained
    std::sort(to_remove.begin(),to_remove.end(),
     [this](const patternIdx & a, const patternIdx & b) -> bool
        {
            if((a.first)<(b.first)){
                return true;
            }
            if(((a.first)==(b.first))&&((a.second)<(b.second))){
                return true;
            }
            return false;
        });

    //Now that the extensions are sorted we remove them.
    patternKey current_key = 0;
    int shift = 0;
    auto it=to_remove.end();
    for(it--;;it--){
        //We directly remove the elements.
        frag_pattern& tpat = *((pmap[it->first].begin())+it->second);
        //std::cout<<"Erasing:"<<tpat.get_norm()<<"eq: "<<it->first<<" "<<it->second<<std::endl;

        //We remove them form the index.
        removePatternIndex(tpat,*it);
        //pmap[it->first].erase((pmap[it->first].begin())+it->second);
        (pmap[it->first].begin()+it->second)->clearPatternFull();
        (*(pmap[it->first].begin()+it->second)).null = true;

        if(it==to_remove.begin()) break;

    }
    //std::cout <<"OUTTT";
}

void subgraph_container::idx_to_string(std::ostream& of){
//    of <<"IdxMap :" << std::endl;
    //std::map<short,std::vector<patternIdx> > lab_idx;
    for(auto it=lab_idx.begin();it != lab_idx.end(); it++){
        of << it->first <<" : ";
        for(auto itt = it->second.begin(); itt != it->second.end() ; itt++){
            of <<" "<<(*itt).first << "_"<<(*itt).second;
        }
        of << std::endl;
    }
//    of <<"EndIdxMap." << std::endl;
}



///LATTICE PART



void subgraph_container::removePattern_latt(std::string& key){
    //The position of the pattern is found.
    auto it = idxmap.find(key);
    Vertexl pos = (*it).second;
    //std::cout <<"Removing_"<<key<<"_"<<pos<<"  ";
    boost::clear_vertex(pos,lat);
    boost::remove_vertex(pos,lat);

    //TODO eventually remove the value from the map
    idxmap.erase(it);
}


//Same butt for a vector of pattern
void subgraph_container::removePattern_latt(std::vector< std::string>& to_rm){
    for(auto it = to_rm.begin();it!=to_rm.end();it++){
       // std::cout << "rm: "<<(*it)<<" ";
        removePattern_latt(*it);
    }
}


int subgraph_container::numPatterns(){
    return boost::num_vertices(lat);
}

void subgraph_container::addPattern_latt(std::string& key,std::vector<std::string>& subgraphs,std::vector<std::string>&supergraphs){
    //The node is added
    Vertexl node = boost::add_vertex(lat);
    lat[node].key = key;
    lat[node].item = false;

    //It is added to the vertex map
    idxmap.emplace(key,node);
    //if(key=="11.95") std::cout << "k:"<<key<<"n:"<<node;

    //We add the parents
    for(auto it = supergraphs.begin();it!=supergraphs.end();it++){
        auto posp = idxmap[*it];

        //An edge is added
        boost::add_edge(node,posp,lat);
    }

    //We add the descendant
    for(auto it = subgraphs.begin();it!=subgraphs.end();it++){
        auto posp = idxmap[*it];

        //An edge is added
        boost::add_edge(posp,node,lat);
    }
}

void subgraph_container::addRoot(){
    //Root is created
    auto root = boost::add_vertex(lat);

    graphTraitl::vertex_iterator bv,ev;

    boost::tie(bv,ev) = boost::vertices(lat);

    for(;bv!=ev;bv++){
        if(boost::in_degree(*bv,lat)==0&&(*bv)!=root){
            boost::add_edge(root,*bv,lat);
        }
    }
}


void subgraph_container::resetIdxMap(){
    idxmap.clear();
    graphTraitl::vertex_iterator bv,ev;
    boost::tie(bv,ev) = boost::vertices(lat);
    for(;bv!=ev;bv++){
        //if(lat[*bv].key=="11.95") std::cout<< "fffound: "<<*bv;
        idxmap[lat[*bv].key] = (*bv);
    }
}
//Create a mapping between the label and the pattern notation
std::map< Vertexl, patternIdx> subgraph_container::getMapping(){

    std::map< Vertexl, patternIdx> res;
    //of <<"map: "<<std::endl;
    //Associate a value to a node at each step.
    for(auto it=pmap.begin();it != pmap.end(); it++){
                //Associate a value to a node at each step.
        for(int i=0; i < it->second.size();i++){
            patternIdx cpair = std::make_pair(it->first,i);

            //We get the label of the current node
            std::string cnorm = it->second[i].get_norm();
            //if(cnorm=="2.180.238")
            //    std::cout << ""<<cnorm<<"_"<<idxmap[cnorm]<<"_"<<lat[idxmap[cnorm]].key<<"_"<<cpair.first<<"_"<<cpair.second <<std::endl;
            res.insert( std::make_pair(idxmap[cnorm],cpair));
//            if(lat[idxmap[cnorm]].key =="11.95"){
//                std::cout <<"norm: "<<cnorm<<"2ndnorm: "<<lat[idxmap[cnorm]].key <<" cpair: "<<cpair.first << "_"<<cpair.second<<std::endl;
//            }
        }
    }

    return res;
}



//Return the pattern with 0 sons.
std::vector<Vertexl> subgraph_container::getLeafs(){
    graphTraitl::vertex_iterator bv,ev;
    std::vector<Vertexl> res;
    boost::tie(bv,ev) = boost::vertices(lat);
    for(bv;bv!=ev;bv++){
        //We check the number of incoming edge.
        if(boost::in_degree(*bv,lat)==0){
            res.push_back(*bv);
        }

    }
    return res;
}

void removeNull(std::map<short,std::vector <frag_pattern> >& pmap){
    for(auto pit = pmap.begin();pit!=pmap.end();pit++){
        int maxSize = pit->second.size()-1;
        while(maxSize>=0){
            if(maxSize<pit->second.size()&&pit->second[maxSize].null){
                pit->second.erase(pit->second.begin()+maxSize);
            }else{
                maxSize--;
            }
        }
    }
}


void subgraph_container::postProcessing(int num_graphs,std::ostream& of){
    //This part add the objects to the graph.

    //The pattenr list is cleared.
    of<<"Cleaning pattern list...";

    removeNull(pmap);

    of << " remaining "<<numPatterns()<<std::endl;
    std::vector< Vertexl > mg_nodes;
    of << "Post-processing of the lattice with "<<boost::num_vertices(lat) <<
    " vertices and "<< boost::num_edges(lat)<<" edges: "<<std::endl;
    for(int i = 0;i<num_graphs;i++){
        Vertexl temp = boost::add_vertex(lat);
        lat[temp].item = true;
        mg_nodes.push_back(temp);
    }


    //HYPOTHESIS: AN ITEM MAY ONLY BE ADDED TO A LEAF.
    //auto leafs = getLeafs();
    of <<"mapping...";
    resetIdxMap();
    std::map< Vertexl, patternIdx> map_pat = getMapping();
    of <<"done..."<<std::endl;
    //An occurences is bad if it is not found in any of the successors of the data, except if the root is a leaf.
    int counter = 0;
    int current_percent = 0;
    int max_val = boost::num_vertices(lat);

    lab_idx.clear();


    graphTraitl::vertex_iterator bv,ev;
    boost::tie(bv,ev) = boost::vertices(lat);

    for(;bv!=ev;bv++){
        //of << std::endl << std::endl;
        Vertexl cnode = *bv;
        //of << "there"<<lat[cnode].key<<" "<<lat[cnode].item <<" addr: "<<cnode;
        if(lat[cnode].item) continue;
        auto pos = map_pat.find(cnode);
        if(pos == map_pat.end()){
            //std::cout << "NEXT !";
            //std::cout << std::endl;
            continue;
        }

        counter++;
        int ppercent = (counter*10)/max_val;
        if(ppercent!=current_percent){
            of<< ppercent*10 << " ";
            current_percent = ppercent;
        }


        //We get a list of the successors.
        graphTraitl::adjacency_iterator ba,ea;
        boost::tie(ba,ea) = boost::adjacent_vertices(*bv,lat);
        auto & occs = get(map_pat[cnode]).get_occs();

        //If it is a leaf
        if(ba==ea){
            for(auto io = occs.begin(); io != occs.end();io++){
                boost::add_edge(cnode,mg_nodes[(*io).gid],lat);
            }
        }else{
            //In this case we check the successors for occurences
            std::set<short> temp_set;

            for(auto io = occs.begin(); io != occs.end();io++){
                temp_set.insert((*io).gid);
            }
            bool bempty = false;
            //We check the succesors
            do{
                if(temp_set.size()==0){
                    bempty = true;
                    break;
                }
                auto & occs_succ = get(map_pat[*ba]).get_occs();
                for(auto io = occs.begin(); io != occs.end();io++){
                    auto te = temp_set.find((*io).gid);
                    if(te!=temp_set.end()){
                        temp_set.erase(te);
                    }
                }

                ba++;
            }while(ba!=ea);

            bempty = (bempty|ba==ea);
            if(bempty){
                continue;
            }

            //In this case we add an edge between the remaining nodes.
            for(auto its=temp_set.begin();its != temp_set.end();its++){
                boost::add_edge(cnode,mg_nodes[*its],lat);
            }
        }
        //of << "done !";

    }
    of <<"100"<< std::endl;
    of << "Post-processing finished: "<<boost::num_vertices(lat) <<" vertices and "<< boost::num_edges(lat)<<" edges."<<std::endl;
}



//####TO_DEBUG ONLY
void subgraph_container::printPatternsLL(std::ostream& of){
    of << "Patterns linked list:" << std::endl;
    //The normal form is extracted for all the pattern
    std::vector < std::string> nform;

    for(auto it=pmap.begin();it != pmap.end();it++){
            of << "key: "<<it->first <<" num_val: "<<it->second.size() << std::endl;
                int co=0;
                for(auto itt=it->second.begin();itt != it->second.end();itt++){
                    of<< itt->get_norm() <<"_"<<co<<"  ";
                    co++;
                }
                of << std::endl;
    }

}


void subgraph_container::printPatternsLA(std::ostream& of){

    of << "Patterns lattice:" << std::endl;
    //The normal form is extracted for all the pattern
    std::vector < std::string> nform;
    graphTraitl::vertex_iterator bv,ev;
    boost::tie(bv,ev) = boost::vertices(lat);

    for(;bv != ev; bv++){
        nform.push_back(lat[*bv].key);
    }

    std::sort(nform.begin(),nform.end());

    //We print all the labels
    for(auto it = nform.begin();it != nform.end();it++){
        of << *it << std::endl;
    }
}


void subgraph_container::printPatternsB(std::ostream& of){

    of << "Patterns lattice / linked list:" << std::endl;
    //The normal form is extracted for all the pattern
    std::vector < std::string> nform;
    //std::vector < Vertexl> addrnform;
    graphTraitl::vertex_iterator bv,ev;
    boost::tie(bv,ev) = boost::vertices(lat);



    for(;bv != ev; bv++){
        nform.push_back(lat[*bv].key);
        //addrnform.push_back(&(*bv));
    }

    std::sort(nform.begin(),nform.end());

    //The normal form is extracted for all the pattern
    std::vector < std::string> nform2;
    //std::vector < frag_pattern*> addrnform2;
    for(auto it=pmap.begin();it != pmap.end();it++){
        for(auto itt = it->second.begin();itt != it->second.end();itt++){
            nform2.push_back(itt->get_norm());
            //addrnform2.push_back(&(*itt));
        }
    }

    std::sort(nform2.begin(),nform2.end());

    //We print all the labels
    auto it = nform.begin();
//    auto itadr = addrnform.begin();
//    auto itadr2 = addrnform2.begin();
    for(auto it2 = nform2.begin();it!=nform.end()&it2 != nform2.end();it2++,it++){
        of << *it <<"   "<<*it2 <<std::endl;
    }
}


void subgraph_container::printPatternsB(){

    std::cout << "Patterns lattice / linked list:" << std::endl;
    //The normal form is extracted for all the pattern
    std::vector < std::string> nform;
    //std::vector < Vertexl> addrnform;
    graphTraitl::vertex_iterator bv,ev;
    boost::tie(bv,ev) = boost::vertices(lat);



    for(;bv != ev; bv++){
        nform.push_back(lat[*bv].key);
        //addrnform.push_back(&(*bv));
    }

    std::sort(nform.begin(),nform.end());

    //The normal form is extracted for all the pattern
    std::vector < std::string> nform2;
    //std::vector < frag_pattern*> addrnform2;
    for(auto it=pmap.begin();it != pmap.end();it++){
        for(auto itt = it->second.begin();itt != it->second.end();itt++){
            if(!itt->null)
            nform2.push_back(itt->get_norm());
            //addrnform2.push_back(&(*itt));
        }
    }

    std::sort(nform2.begin(),nform2.end());

    //We print all the labels
    auto it = nform.begin();
//    auto itadr = addrnform.begin();
//    auto itadr2 = addrnform2.begin();
    for(auto it2 = nform2.begin();it!=nform.end()&it2 != nform2.end();it2++,it++){
        std::cout << *it <<"__"<<*it2 <<"  ";
    }
}

