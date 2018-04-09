#include <iostream>
#include <tuple>
#include <string>
#include <iomanip>


#include "subgraph_container.h"


/**
The structure defined in this files are used to create the definitive lattice.


**/






//THis function add a node to exlat for each node of the lattice and return a dictionnary of nodes.
std::map<Vertexl,Vertexel> subgraph_container::latticeMapping(latticeexport& exlat){

    std::map<Vertexl,Vertexel> res;
    //We pass on the lattice
    graphTraitl::vertex_iterator bv,ev;
    boost::tie(bv,ev) = boost::vertices(lat);

    //A node is added for each node of the old lattice.
    for(;bv != ev; bv++){
        Vertexel temp = boost::add_vertex(exlat);
        //if(lat[*bv].key=="2.180.238") std::cout << *bv <<"_" << lat[*bv].key;
        res.insert(std::make_pair(*bv,temp));
    }
    return res;
}


//This function define a mapping between the node of the initial lattice and the node of the new lattices.

std::string int_to_string(int val){
    std::ostringstream convert;
    convert << val;
    return convert.str();
}

std::map<Vertexl,Vertexel> subgraph_container::createExportLattice(latticeexport& elat){

    resetIdxMap();

    //Initialisation
    std::map<Vertexl,Vertexel> vmap = latticeMapping(elat);
    //std::cout << "mapping done" << std::endl; //DEBUG

    //We generate the ids, on for the datasets, one for the nodes.
    int count_M = 0;
    int count_I = 0;

    std::map< Vertexl, patternIdx> map_l_pat = getMapping();
    //std::cout << "2nd mapping done" << std::endl; //DEBUG
    //A mapping which link node is on the original

    for(auto it = vmap.begin();it != vmap.end();it++){
        //The new value is generated
        std::string cid = "";
        //std::cout << "npat: "<<lat[it->first].key << " addr "<<it->first; //DEBUG
        if(lat[it->first].item){
            //std::cout << "item";
            count_I++;
            cid = "G"+int_to_string(count_I);
            elat[it->second].item = 1;
        }else{
            //std::cout << "pattern";
            count_M++;
            cid = "M"+int_to_string(count_M);
            elat[it->second].item = 0;
        }
        elat[it->second].id = cid;

        //Determination of the other node information.
        if(lat[it->first].item){
            elat[it->second].noccs  = 1;
            elat[it->second].num_losses  = 0;
            elat[it->second].noccs_unique = 1;
        }else{
            patternIdx patId = map_l_pat[it->first];
            elat[it->second].noccs  = get(patId).numOccs();
            //Correcting for the loss.
            elat[it->second].num_losses  = get(patId).sizeGraph()-1;
            elat[it->second].noccs_unique  = get(patId).numUniqueOccs().size();
        }
        //std::cout <<"attr" << " "; //DEBUG

        //Currently the socre is initialized to 0:
        elat[it->second].score = 0;

        //In this part we add the corresponding edge.
        graphTraitl::out_edge_iterator bo,eo;
        boost::tie(bo,eo) = boost::out_edges(it->first,lat);

        //The majority of these shit is okay.
        for(;bo!=eo;bo++){
            auto nsource = vmap[boost::source(*bo,lat)];
            auto ntarget = vmap[boost::target(*bo,lat)];

            //The edge is added
            boost::add_edge(nsource,ntarget,elat);
        }
    }

    return vmap;

}

//path_lattice is a file name indicating all the values of the dataset.
//path_path_patterns is a directory of output ofthe data
void subgraph_container::exportResultingLattice(std::string path_lattice,std::string path_patterns,
                                                std::ostream& of){
    //The lattice is created.
    latticeexport elat;
    std::map<Vertexl,Vertexel> vmap = createExportLattice(elat);

    std::map< Vertexl, patternIdx> map_l_pat = getMapping();

    //Node properties are defined.
    boost::dynamic_properties dp(boost::ignore_other_properties);
    dp.property("id", boost::get(&export_node_info::id, elat));
    dp.property("noccs", boost::get(&export_node_info::noccs, elat));
    dp.property("num_losses", boost::get(&export_node_info::num_losses, elat));
    dp.property("score", boost::get(&export_node_info::score, elat));
    dp.property("item", boost::get(&export_node_info::item, elat));
    dp.property("i", boost::get(&export_node_info::item, elat));
    //It's not exported, the number of unique occurences is calculated in the graph.


    std::ofstream ofile(path_lattice);
    boost::write_graphml(ofile, elat, dp, true);
    ofile.close();


    //The pattern are exported using this method.
    graphTraitl::vertex_iterator bv,ev;
    //std::cout << "Writing patterns." << std::endl;
    boost::tie(bv,ev) = boost::vertices(lat);
    for( ;bv != ev ; bv++){
        if(lat[*bv].item) continue;
        std::string idname = elat[vmap[*bv]].id;

        std::string path_file = (path_patterns+"/"+idname+".graphml");
        //std::cout << "WWriting: "<<path_file <<std::endl;
        frag_pattern& fp = get(map_l_pat[*bv]);
        //The node is exported.
        fp.write_graphml(path_file);
        //std::cout << "Written: "<<path_file <<std::endl;
    }
    //of << "Patterns written." << std::endl;
}



