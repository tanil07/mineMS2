#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

#include <tuple>

#define BOOST_NO_AUTO_PTR

#ifdef __APPLE__
    _Pragma("clang diagnostic push")
    _Pragma("clang diagnostic ignored \"-Wnonnull\"")
    _Pragma("clang diagnostic ignored \"-Wparentheses\"")
    _Pragma("clang diagnostic ignored \"-Wmaybe-uninitialized\"")

    #include <boost/graph/graph_traits.hpp>
    #include <boost/graph/adjacency_list.hpp>

    _Pragma("clang diagnostic pop")
#else
    _Pragma("GCC diagnostic push")
    _Pragma("GCC diagnostic ignored \"-Wnonnull\"")
    _Pragma("GCC diagnostic ignored \"-Wparentheses\"")
    _Pragma("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")

    #include <boost/graph/graph_traits.hpp>
    #include <boost/graph/adjacency_list.hpp>

    _Pragma("GCC diagnostic pop")
#endif 


#include <boost/functional/hash.hpp>


#define K_PATH_ROOT -1
#define NULL_LABEL -1
#define REMOVE_EXTENSION true
#define PATH_TREE_ROOT -1
#define ITEM_ID = 0;



//Properties of mass_graph
struct node_info{
    //TODO change it to short value when graphml bug is fixed.
    short lab;
    double mz;
    double intensity;
};

struct edge_info{
    //TODO change it to short value when graphml bug is fixed.
    short lab;
};

//properties of k-tree
struct k_tree_info{
    short lab=0;
    short dist=0;
};



//THis information is filled only when exporting.
struct pgraph_info{
    std::string occs;
    int unique_occs;
    std::string norm;
};

//properties of pattern to avoid copy issue a constructor is added
struct pnode_info{
    pnode_info(short depth=0,short lab=0): depth(depth),lab(lab) {}
    short depth; //depth is short because it not supposed ot be bigger than 3.
    short lab; //correspond to the distance form the precursor.
};

struct pedge_info{
    pedge_info(bool span = false,short lab=-1): span(span),lab(lab) {}
    bool span; //is the edge a spanning edge.
    short lab; //the edge lab.
};

//We need to add na operation for the structre to be hashable.
struct occ{
    short gid;
    short idx;
};
inline bool operator<(const occ& lhs, const occ& rhs)
{
    if(lhs.gid < rhs.gid) return true;
    if(lhs.gid > rhs.gid) return false;
    if(lhs.idx < rhs.idx) return true;
    return false;
}

//Typedef of the k-path tree.
typedef boost::adjacency_list<
    boost::setS, boost::setS, boost::bidirectionalS,
    k_tree_info, boost::no_property> ktree;
typedef boost::graph_traits<ktree> graphTraitst;
typedef boost::graph_traits<ktree>::vertex_descriptor Vertext;

//Occurrences list
typedef std::map<Vertext, std::set<occ> > MapOccurrences;





//adj structure info
struct adjnode_info{
    adjnode_info(short lab=0,graphTraitst::vertex_descriptor v=(graphTraitst::null_vertex())): lab(lab),v(v){}
    short lab; //correspond to the distance form the precursor.
    //Store the corresponding node in the k-path tree for commodity purpose
    graphTraitst::vertex_descriptor v;
};



//Typedef mass_graph
typedef boost::adjacency_list<
    boost::vecS, boost::vecS, boost::bidirectionalS,
    node_info, edge_info> graph;
typedef boost::graph_traits<graph> graphTraits;
typedef boost::graph_traits<graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<graph>::edge_descriptor Edge;
typedef std::map<Edge,short> VisitMap;

//Property map
typedef boost::property_map<graph,boost::vertex_index_t>::type IndexMap;


//Typedef of frag_pattern
//!!!!!!!!!! LET THIS STRUCTURE IN VEC AT ALL COST !!!!!!!!!
typedef boost::adjacency_list<
    boost::setS, boost::vecS, boost::bidirectionalS,
    pnode_info, pedge_info, pgraph_info> graphp;
typedef boost::graph_traits<graphp> graphTraitsp;
typedef boost::graph_traits<graphp>::vertex_descriptor Vertexp;
typedef boost::graph_traits<graphp>::edge_descriptor Edgep;
//We define the extension vector

//Extension, extended node in the pattern, corresponding node in in $T$
typedef std::tuple<Vertexp,Vertext,short> Extension;


//Typedef of the adjacency list
typedef boost::adjacency_list<
    boost::setS, boost::setS, boost::undirectedS,
    adjnode_info, boost::no_property> graphadj;
typedef boost::graph_traits<graphadj> graphTraitsAdj;
typedef boost::graph_traits<graphadj>::vertex_descriptor Vertexadj;



struct lat_info{
    std::string key;
	int id;
    //indicate if the node is an object or None.
    bool item;
};


//Typedef used in the lattice with subgraph relationship
//We need to find the the node quickly, and the edge set
typedef boost::adjacency_list<
    boost::setS, boost::setS, boost::bidirectionalS,
    lat_info, boost::no_property> latticesub;
typedef boost::graph_traits<latticesub> graphTraitl;
typedef boost::graph_traits<graphadj>::vertex_descriptor Vertexl;
typedef boost::graph_traits<graphadj>::edge_descriptor Edgel;

//We store a map between the values.

//TODO check if map is more fficient
typedef boost::unordered_map<std::string,Vertexl> subgraph_idx;




//Typedef used for the pattern container
typedef short patternKey;
typedef std::pair<patternKey,short> patternIdx;

//Operator
struct less_than
{
    inline bool operator() (const Extension& ext1, const Extension& ext2)
    {
        return ((std::get<0>(ext1) < std::get<0> (ext2)) ||
                (std::get<2> (ext1) < std::get<2> (ext2)));
    }
};

#endif // COMMON_H_INCLUDED
