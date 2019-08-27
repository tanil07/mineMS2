#include<iostream>

#include "adjacencyGraph.h"
#include "k_path_tree.h"

adjacencyGraph::adjacencyGraph()
{
    //ctor
}


//Initialize the graph with n nodes correspoding to the n labels
adjacencyGraph::adjacencyGraph(int n){
    for(int i=0;i<n;i++){
        Vertexadj nadj = boost::add_vertex(g);
        g[nadj].lab=i;
        lab_node.insert(std::make_pair(i,nadj));
    }
}

graphadj& adjacencyGraph::get_g(){
    return g;
}



Vertexadj adjacencyGraph::getNode(short lab){
        auto it = lab_node.find(lab);
        if(it==lab_node.end()){
            Vertexadj nadj = boost::add_vertex(g);
            g[nadj].lab=lab;
            lab_node.insert(std::make_pair(lab,nadj));
            return nadj;
        }else{
            return (*it).second;
        }
}


//Function of update which use a node and a graph.
void adjacencyGraph::add_adj(Vertex v,graph& mg){
    graphTraits::out_edge_iterator b1,e1,b2,e2;
    boost::tie(b1,e1) = boost::out_edges(v,mg);
    //boost::tie(b2,e2) = boost::out_edges(v,mg);
    for(;b1!=e1;b1++){
        //We check if the node exists.
        short lab1 = mg[*b1].lab;
        Vertexadj v1 = getNode(lab1);
        for(b2=b1,b2++;b2!=e1;b2++){
            short lab2 = mg[*b2].lab;
            Vertexadj v2 = getNode(lab2);
            boost::add_edge(v1,v2,g);
        }
    }
};


void adjacencyGraph::add_graph(graph& mg){
    graphTraits::vertex_iterator bv,ev;
    boost::tie(bv,ev)=boost::vertices(mg);
    for(;bv!=ev;bv++){
        add_adj(*bv,mg);
    }
}

void adjacencyGraph::addKTreeVertices(k_path_tree& kt){
    const Vertext rt = kt.get_root();
    const ktree& t = kt.get_t();
    //std::cout << "addr1 t "<<&t<< std::endl;//DEBUG
    //For every successors node we add the correct vertex
    graphTraitst::adjacency_iterator bv,ev;
    for(boost::tie(bv,ev)=boost::adjacent_vertices(rt,t);bv!=ev;bv++){
        //We get the correspoding label
        short tlab = t[*bv].lab;
        g[lab_node[tlab]].v = (*bv);
        //std::cout << "added " << tlab << " to " << g[lab_node[tlab]].lab << "updated " <<g[lab_node[tlab]].v<<"vs" <<*bv<<std::endl;//DEBUG
    }
}

//Return the neighbouring nodes.
//To do implement the lower pattern filter.
std::vector<std::pair <Vertext,short> > adjacencyGraph::neighbours(short val){
    graphTraitsAdj::adjacency_iterator bv,ev;

    //We get the set of adjacents vertices.
    boost::tie(bv,ev)=boost::adjacent_vertices(lab_node[val],g);

    //We now comppute the set of label
    std::vector<std::pair <Vertext,short> > res;
    res.reserve(100);

    std::transform(bv,ev,std::back_inserter(res),
                  [this](Vertexadj v)-> std::pair <Vertext,short> {
                  return std::make_pair((this->g[v]).v,(this->g[v]).lab);
                  ;});

    return res;
}

void adjacencyGraph::remove_nodes(std::vector<short>& labs){
//Nodes are removed from the graph if they are not frequent.
    for(auto it = labs.begin();it!=labs.end();it++){
        //We remove the corresponding node
        Vertexadj vadj = lab_node[*it];
        boost::clear_vertex(vadj,g);
        boost::remove_vertex(vadj,g);
    }
}

void adjacencyGraph::keep_nodes(std::vector<short>& labs){
    std::vector<Vertexadj> to_keep;
    std::transform(labs.begin(),labs.end(),std::back_inserter(to_keep),
                   [this](short i)->Vertexadj{
                        return this->lab_node[i];
                   });
    //We remove the values cotained in the differences set.
    std::unordered_set<short> temp_set(labs.begin(),labs.end());

    //We find the difference of the 2 sets.
    std::vector<Vertexadj> to_rm;
    for(auto it= lab_node.begin(); it!= lab_node.end();it++){
        short clab = (*it).first;
        if(temp_set.find(clab)==temp_set.end()){
            to_rm.push_back((*it).second);
        }
    }
    //Nodes are removed from the graph if they are not frequent.
    for(auto it = to_rm.begin();it!=to_rm.end();it++){
        //We remove the corresponding node
        boost::clear_vertex(*it,g);
        boost::remove_vertex(*it,g);
    }
}

adjacencyGraph::~adjacencyGraph()
{
    //dtor
}
