#include <iostream>
#include <algorithm>

#include "k_path_tree.h"
#include "common.h"
#include "utility_functions.h"
#include "lattice_node.h"

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/array.hpp>


k_path_tree::k_path_tree(int k): k(k)
{
    //an empty graph containing only a root is created.
    //ctor

    std::vector<std::vector<Vertext> > tvec(k);
    for(int i =0; i<k; i++)
    {
        std::vector<Vertext> temp;
        tvec[i]=temp;
    }

    root = boost::add_vertex(t);
    t[root].lab = 0;
    pos = tvec;

}

k_path_tree::~k_path_tree()
{
    //dtor
}


int k_path_tree::get_k()
{
    return k;
}

ktree& k_path_tree::get_t()
{
    return t;
}

MapOccurrences& k_path_tree::get_occs()
{
    return moccs;
}

adjacencyGraph& k_path_tree::get_adj()
{
    return adj;
}

Vertext k_path_tree::get_node(Vertext origin,int lab, int dist)
{
    //We iterate on the successors of origin and fin the node if it exsits
    graphTraitst::adjacency_iterator bv,ev;
    for(boost::tie(bv,ev) = boost::adjacent_vertices(origin,t); bv!=ev ; bv++)
    {
        if(t[*bv].lab==lab)
        {
            return *bv;
        }
    }
    //No found vertex. We create the vertex
    Vertext nnode = boost::add_vertex(t);
    boost::add_edge(origin,nnode,t);
    t[nnode].lab = lab;
    t[nnode].dist = dist;
    return nnode;
}


//Modification to add only the frequent path.

//Utility function used to add a node.
//We add the the occurrences values
void k_path_tree::update_pos_adv(int lab, std::vector<Vertex>& pfr,std::vector<int>& plabs,
                                 int lpath, IndexMap& idx_vertex, VisitMap& vm,int gid, graph& G)
{
    bool pres;
    Edge e;
    short elab;

    //We consider that pi is the size of the added path.
    for(int pi=k;pi>2;pi--)
    {
        if(pi>lpath) break;

        //We check that this occurrences needs to be added.
        boost::tie(e,pres)=boost::edge(pfr[lpath-pi],pfr[lpath],G);
        if(pres){
            elab = G[e].lab;
        }else{//TODO check if it is possible to return directly.
            continue;
        }

        if(pos[pi-1].size()==0)
        {
            if(pos[pi-2].size()!=0)
            {
                Vertext vn = this->get_node(pos[pi-2].back(),plabs[lpath-1],elab);
                pos[pi-1].push_back(vn);
                occ oc = {short(gid),short(idx_vertex[pfr[lpath-pi]])};
                addOccs(moccs,vn,oc);
            }
            continue;
        }

        //Case where the last element is added.
        if(pi!=k)
        {
            Vertext old_vertice = pos[pi-2].back();
            //String obtained by adding the value.
            
            Vertext vn = this->get_node(old_vertice,plabs[lpath-1],elab);
            pos[pi-1].push_back(vn);

            //We create the occurrence :
            occ oc = {short(gid),short(idx_vertex[pfr[lpath-pi]])};
            addOccs(moccs,vn,oc);
        }
    }
        occ oc = {short(gid),short(idx_vertex[pfr[lpath-1]])};
    Vertext vn = this->get_node(root,plabs[lpath-1],plabs[lpath-1]);
    pos[0].push_back(vn);
    addOccs(moccs,vn,oc);

}

void k_path_tree::update_pos_back()
{
    //We just need to remove the last element of every node if there is one.
    for(int pi=0; pi<int(pos.size()); pi++)
    {
        if(pos[pi].size()>0)
        {
            pos[pi].pop_back();
        }
        else
        {
            return;
        }
    }
}

//Normally we never have to find path.
Vertext k_path_tree::get_root()
{
    return root;
}



short k_path_tree::get_dist_vertex(Vertext c){
    return(t[c].dist);
}



void k_path_tree::add_graph(mass_graph& G,int gid, bool prec_only)
{

    //We find all the possible root nodes.
    graph& g = G.get_g();

    //If the graph have no edge we stop right there.
    if(boost::num_edges(g)==0) return;

    //We get the visitor map on the value.
    VisitMap vm = G.buildVisitMap();

    //Finding the roots of a graph.
    std::vector<Vertex> roots;

    if(!prec_only)
    {
        roots = G.roots();
    }
    else
    {
        roots.push_back(G.get_precursor());
    }

    //Defining a node index map
    IndexMap imap = boost::get(boost::vertex_index,g);

    //Vector storing the path from the root.
    std::vector<Vertex> pathFromRoot(boost::num_vertices(g));

    //Vector stroing the labs from the root (simplicity)
    std::vector<int> edgesLabels(boost::num_vertices(g)-1);
    Vertex nvertex = graphTraits::null_vertex();


    for(int i=0; i<int(roots.size()); i++)
    {
        //Root initialization.
        Vertex root = roots[i];
        pathFromRoot[0] = root;
        int ppath = 0;
        bool backward = false;

        Vertex dnode,cnode;
        dnode = nextNode(vm, g, root);
        cnode = root;
        bool null_temp = false;
        while( (cnode!=root) | (dnode!=nvertex) )
        {

            bool bc_check = false;
            //We check that the node have not already been explored.
            if(ppath>=1){
                int idfirst =(ppath-k)>0 ? (ppath-k):0;
                Edge e;
                boost::tie(e,null_temp)=edge(pathFromRoot[idfirst],pathFromRoot[idfirst+1],g);
                if(vm[e]==2){
                    bc_check = true;
                }
            }
            
            //Case where we have to backtracks the DFS tree.
            if((dnode==nvertex)|
                    (prec_only&(cnode != root))|bc_check)
            {
                //We goes back to the last
                ppath--;
                cnode = pathFromRoot[ppath];
                backward = true;
                this->update_pos_back();
            }
            else   //Case where we go down in the graph.
            {

                cnode = dnode;
                ppath++;
                pathFromRoot[ppath] = cnode;

                edgesLabels[ppath-1] = g[boost::edge(pathFromRoot[ppath-1],cnode,g).first].lab;
                backward = false;
                if(!backward)
                {
                    this->update_pos_adv(edgesLabels[ppath-1], pathFromRoot,
                                         edgesLabels,ppath, imap,vm, gid, g);

                }
            }
            dnode = nextNode(vm ,g, cnode);
        }
    }
    //Finally the triangle list is added to the tl
    adj.add_graph(g);
}
void k_path_tree::post_processing()
{
//We construct the triangle list mappin
    adj.addKTreeVertices(*this);
}


std::vector<Extension> k_path_tree::getExtensions(Vertexp v,Vertext vt)
{

    //Now we create the extension.
    graphTraitst::adjacency_iterator bo,be;
    boost::tie(bo,be)=boost::adjacent_vertices(vt,t);
    //Now we add all the extension if possible.
    std::vector<Extension> temp_exts;
    std::transform(bo,be,std::back_inserter(temp_exts),
                   [v,this](Vertext ve)->Extension
    {
        Extension t_ext = std::make_tuple(v,ve,this->t[ve].lab);
        return t_ext;
    });
    std::sort(temp_exts.begin(),temp_exts.end(),less_than());
    return temp_exts;
}

//UTILITY FUNCTION
Vertext k_path_tree::find_pos(std::vector<short>path)
{
    Vertext cnode = get_root();
    for(int i = 0; i < int(path.size()); i++)
    {
        short clab = path[i];
        graphTraitst::adjacency_iterator vb,ve;
        for(boost::tie(vb,ve)=boost::adjacent_vertices(cnode,t) ; vb!=ve ; vb++)
        {
            if(t[*vb].lab==clab)
            {
                cnode = *vb;
                break;
            }
        }
        if(vb == ve) return(graphTraitst::null_vertex());
    }
    return cnode;
};

//Filter out the non frequent nodes.
void k_path_tree::filter_frequent_nodes(int noccs)
{
    graphTraitst::vertex_iterator vb,ve;
    boost::tie(vb,ve)=boost::vertices(t);

    //We store the node to remove.
    std::vector<Vertext> to_rm;
    to_rm.reserve(1024);

    Vertext root_v = get_root();
    boost::tie(vb,ve)=boost::vertices(t);
    //TODO pass this in one pass possibly if it's slow.
    for(; vb!=ve; vb++)
    {
        //Never remove the root.
        if((*vb)==root_v)
        {
            continue;
        }
        if(int(moccs[*vb].size())<noccs)
        {
            to_rm.push_back(*vb);
        }
    }
    //Now we remove the chosen node.
    for(int ir=0; ir<int(to_rm.size()); ir++)
    {
        Vertext v= to_rm[ir];
        boost::clear_vertex(v,t);
        boost::remove_vertex(v,t);
        auto itpos = moccs.find(v);
        moccs.erase(itpos);
    }

    //We get the remaining 0 exts en update the adjacency graph
    std::vector<short> r0;
    graphTraitst::adjacency_iterator bv,ev;
    boost::tie(bv,ev) = boost::adjacent_vertices(root_v,t);
    std::transform(bv,ev,std::back_inserter(r0),
                   [this](Vertext v) -> short {return this->t[v].lab;});
    adj.keep_nodes(r0);
}



//Return the path from the root of the vertex v
std::vector<Vertext> k_path_tree::find_predecessors(Vertext v)
{
    Vertext cnode = v;
    Vertext dest = get_root();
    std::vector<Vertext> pathRoot;
    pathRoot.push_back(v);
    int num_it = 0;
    while(cnode!=dest)
    {
        if(num_it>=k)
        {
            std::vector<Vertext> tav;
            return tav;
        }
        graphTraitst::in_edge_iterator eb,ee;
        boost::tie(eb,ee)=boost::in_edges(cnode,t);
        cnode = boost::source(*eb,t);
        pathRoot.insert(pathRoot.begin(),cnode);
        num_it++;
    }
    return pathRoot;
}

//to_string function used to debug.
 void k_path_tree::to_string(std::ostream& of){
     graphTraitst::vertex_iterator bv,ev;

     //For each edge we plot the correspoding path.
     for(boost::tie(bv,ev)=boost::vertices(t);bv!=ev;bv++){
         if(*bv==root) continue;
         std::vector<Vertext> cpath = find_predecessors(*bv);
         //We print all the correspoding labes
         of << "path : ";

         std::vector<int> pstr(cpath.size());
         std::transform(cpath.begin(),cpath.end(),pstr.begin(),
                        [this](Vertext v) -> int {return t[v].lab;});
         for(long long unsigned int ist=0;ist<pstr.size();ist++){
             of << pstr[ist];
             if(ist<pstr.size()-1) of << "_";
         }
         //We the number of occurrences
         of << " noccs : " << moccs[*bv].size()<<std::endl;
         of << "occs : ";
         for(auto it=(moccs[*bv]).begin();it!=(moccs[*bv]).end();it++) of<<(*it).gid<<"_"<<(*it).idx<<" ";
         of << std::endl;

         //TODO add the occurrences enveutally.
     }
 }

short k_path_tree::getLab(Vertext v) const{
    return t[v].dist;
}

std::vector<short> k_path_tree::getLabSuccs(Vertext v) const{
    graphTraitst::adjacency_iterator ba,ea;
    boost::tie(ba,ea) = boost::adjacent_vertices(v,t);
    std::vector<short> to_return;
    std::transform(ba,ea,std::back_inserter(to_return),
                   [this](Vertext v) -> short {return this->t[v].dist;});
    return to_return;
}

//Function used to construct one edged values
//Construct the set of one edge fragment
std::vector<lattice_node> k_path_tree::constructOneEdgeGraphs(std::ostream& of,bool sorted)
{
    //We first get a list of all the labels.
    graphTraitst::adjacency_iterator vb,ve;
    boost::tie(vb,ve)=boost::adjacent_vertices(root,t);

    //Counter for commodity purpose.
    int nadj = 0;
    while(vb!=ve)
    {
        vb++;
        nadj++;
    }
    boost::tie(vb,ve)=boost::adjacent_vertices(root,t);
    std::vector<lattice_node> single_e_patterns;

    int num = 0; //DEBUG
    int percent10 = 0;
    for(; vb!=ve; vb++)
    {
        //We create the pattern;
        int dec = ((double)num*10)/nadj;
        if(dec!=percent10)
        {
            percent10 = dec;
        }
        single_e_patterns.push_back(lattice_node((*vb), *this,of));
        num++;
    }

    //Sorting the patterns.
    if(sorted)
    {
        std::sort(single_e_patterns.begin(),single_e_patterns.end(),
                  [](const lattice_node& a,const lattice_node& b)->bool
        {
            return a.get_key()<b.get_key();
        });
    }
    return single_e_patterns;
}
