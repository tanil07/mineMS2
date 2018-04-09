#include<iostream>
#include <map>

#include "triangles_list.h"
#include "mass_graph.h"

#define INIT_SIZE 512

triangles_list::triangles_list()
{
    //triplets.reserve(INIT_SIZE);
    vtriplets.reserve(INIT_SIZE);
}

triangles_list::~triangles_list()
{
    //dtor
}

        //Adding a graph to the triangle list.
void triangles_list::add_mass_graph(mass_graph& mg){
    //We frst order all the vertices of the data
    graph g = mg.get_g();
    graphTraits::vertex_iterator bv,ev;
    boost::tie(bv,ev)=boost::vertices(g);

    //A copy of the vertices of the graph
    std::vector<Vertex> vvecs(bv,ev);

    //Soritng on the masses.
    std::sort(vvecs.begin(),vvecs.end(),[g](const Vertex a,
                                            const Vertex b) -> bool {return g[a].mz>g[b].mz;});

    //Now we get all the triangles
    for(int ia=0;ia<(vvecs.size()-2);ia++){
        graphTraits::adjacency_iterator bva,eva;
        boost::tie(bva,eva) = boost::adjacent_vertices(vvecs[ia],g);
        std::vector<Vertex> vadj_a(bva,eva);
        //std::cout << "a succs : "<<vadj_a.size()<<" "<<g[vvecs[ia]].mz<<std::endl;
        if(vadj_a.size()==0){ continue;}
        for(int ib=0;ib<(vadj_a.size()-1);ib++){
            //std::cout << "  b : "<<g[vadj_a[ib]].lab<<std::endl;
            graphTraits::adjacency_iterator bvb,evb;
            boost::tie(bvb,evb) = boost::adjacent_vertices(vadj_a[ib],g);
            std::vector<Vertex> vadj_b(bvb,evb);
            if(vadj_b.size()==0){ continue;}
            //Intersection.
            std::vector<Vertex> v_intersection;
            std::set_intersection(vadj_a.begin(),vadj_a.end(),
                                  vadj_b.begin(),vadj_b.end(),
                                  std::back_inserter(v_intersection));

            //Checking if the triplet should be added.
            if(v_intersection.size()>0){
                for(auto ita=v_intersection.begin();ita!=v_intersection.end();ita++){

                    graphTraits::edge_descriptor ea,eb,ec;
                    ea = boost::edge(vvecs[ia],vadj_a[ib],g).first;
                    eb = boost::edge(vadj_a[ib],*ita,g).first;
                    ec = boost::edge(vvecs[ia],*ita,g).first;
                    boost::array<short,3> temp_tri = {g[ea].lab,g[eb].lab,g[ec].lab};

                    bool inserted=true;
                    std::set<boost::array<short, 3> >::iterator it_tri;
                    boost::tie(it_tri,inserted)=triplets.insert(temp_tri);

                    //We copy
                    if(inserted){
                        vtriplets.push_back(temp_tri);
                    }

                }
            }
        }
    }
}


//All the index are constructed there.
void triangles_list::construct_mapping(){

    for(int i=0;i<vtriplets.size();i++){
        short alab = std::get<0>(vtriplets[i]);
        short blab = std::get<1>(vtriplets[i]);
        short clab = std::get<2>(vtriplets[i]);
        std::pair<short,short> ablab = std::make_pair(alab,blab);
        std::pair<short,short> aclab = std::make_pair(alab,clab);

        if(alab==44&clab==116) std::cout << std::endl << alab << "|" << blab << "|" << clab << std::endl;

        //Construction of a index.
        auto ita = idxa.find(alab);
        if(ita==idxa.end()){
            //In this case we create the vector
            std::vector<short> tempa;
            tempa.push_back(i);
            idxa.emplace(alab,tempa);
        }else{
            ((*ita).second).push_back(i);
        }

        //b
        auto itb = idxb.find(blab);
        if(itb==idxb.end()){
            //In this case we create the vector
            std::vector<short> tempb;
            tempb.push_back(i);
            idxb.emplace(blab,tempb);
        }else{
            ((*itb).second).push_back(i);
        }

        //c
        auto itc = idxc.find(clab);
        if(itc==idxc.end()){
            //In this case we create the vector
            std::vector<short> tempc;
            tempc.push_back(i);
            idxc.emplace(clab,tempc);
        }else{
            ((*itc).second).push_back(i);
        }

        //ab
        auto itab = idxab.find(ablab);

        if(itab==idxab.end()){
            idxab.insert(std::make_pair(ablab,clab));
        }else{
            ((*itab).second)=clab;
        }


        //ac
        auto itac = idxac.find(aclab);

        if(itac==idxac.end()){
            idxac.insert(std::make_pair(aclab,blab));
        }else{
            ((*itac).second)=blab;
        }
    }
}

//
////Current implmentation use intersection.
//short triangles_list::get_ab(short a,short b) const {
//    auto iida = idxa.find(a);
//    if(iida==idxa.end()){
//        return NULL_LABEL;
//    }
//    auto iidb = idxb.find(b);
//    if(iidb==idxb.end()){
//        return NULL_LABEL;
//    }
//    auto ida = ((*iida).second);
//    //As there is jsut one possiblity we loop thtought the idx of a until we find b
//    for(auto ia = ida.begin();ia != ida.end();ia++){
//        if(std::get<1>(vtriplets[*ia])==b){
//            return std::get<2>(vtriplets[*ia]);
//        }
//    }
//
//    return NULL_LABEL;
//}

//Version using ab index.
short triangles_list::get_ab(short a,short b) const {
    auto ab = std::make_pair(a,b);
    auto iidab = idxab.find(ab);
    if(iidab==idxab.end()){
        return NULL_LABEL;
    }else{
        return (*iidab).second;
    }
}


//
////return the b term.
//short triangles_list::get_ac(short a,short c) const {
//    //std::cout << "a "<<a<<" c "<<c << std::endl; //DEBUG
//    auto iida = idxa.find(a);
//    if(iida==idxa.end()){
//        //std::cout << "nonA "; //DEBUG
//        return NULL_LABEL;
//    }
//    auto iidc = idxc.find(c);
//    if(iidc==idxc.end()){
//        //std::cout << "nonC "; //DEBUG
//        return NULL_LABEL;
//    }
//    auto ida = ((*iida).second);
//    //As there is jsut one possiblity we loop thtought the idx of a until we find b
//    for(auto ia = ida.begin();ia != ida.end();ia++){
//        if(std::get<2>(vtriplets[*ia])==c){
//            //std::cout << "A&B "; //DEBUG
//            return std::get<1>(vtriplets[*ia]);
//        }
//    }
//    return NULL_LABEL;
//}

short triangles_list::get_ac(short a,short c) const {
    auto ac = std::make_pair(a,c);
    auto iidac = idxac.find(ac);
    if(iidac==idxac.end()){
        return NULL_LABEL;
    }else{
        return (*iidac).second;
    }
}


        //There could be multiples accessors for these values.
std::vector<short> triangles_list::get_a_c(short a) const {

    //We return the vector of index
    auto ita = idxa.find(a);
    if(ita==idxa.end()){
        std::vector<short> temp;
        return temp;
    }
    std::vector<short> ra= (*ita).second;

    //We find the associated C values.
    std::transform(ra.begin(),ra.end(),ra.begin(),
                [this](int i){

                  return std::get<2>(this->vtriplets[i]);
                  });

    return ra;
}

std::vector<short> triangles_list::get_a_b(short a) const {

    //We return the vector of index
    auto ita = idxa.find(a);
    if(ita==idxa.end()){
        std::vector<short> temp;
        return temp;
    }
    std::vector<short> ra= (*ita).second;

    //We find the associated C values.
    std::transform(ra.begin(),ra.end(),ra.begin(),
                [this](int i){
                  return std::get<1>(this->vtriplets[i]);
                  });

    return ra;
}

//std::vector<short> triangles_list::get_c(short c){
//
//}


void triangles_list::to_string(){
    std::cout << "Triangles : "<<triplets.size()<<std::endl;
    if(triplets.size()<10){
        for(auto it_arr=triplets.begin();it_arr!=triplets.end();it_arr++){
            for(int j=0;j<(*it_arr).size();j++){
                std::cout<<(*it_arr)[j]<<"_";
            }
            std::cout << std::endl;
        }
    }
}
//
//void triangles_list::resize_triplets(){
//    //We double the SIze of the vectoR;
//    if(current_size>=triplets.size()){
//            triplets.resize(2*triplets.size())
//    }
//}
