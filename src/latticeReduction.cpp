#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;

#include <boost/graph/depth_first_search.hpp>

#include "common.h"


//Definition of the graph used for reduction.
struct latt_info{
	//TODO change it to short value when graphml bug is fixed.
	short name;
	double score;
};

//A new kind of graph is defined jsut to perform the reduction
typedef boost::adjacency_list<
	boost::vecS, boost::vecS, boost::bidirectionalS,
	latt_info, boost::no_property> latt;

typedef boost::graph_traits<latt> graphTraitLatt;

typedef graphTraitLatt::vertex_descriptor VertexLatt;


//Update of the best_socre magtrix using the occurences table.
void update_dist(S4 pattern,IntegerMatrix& names,NumericMatrix& best_scores, int pid, double score){
	//We get the occurences
	IntegerMatrix occs = pattern.slot("occurences");

	IntegerVector gids = occs( _ , 0);

	//NumericVector rscore = pattern.slot("scores");
	//double score = as<double>(rscore);

	IntegerVector::iterator i1,i2;
	for(i1=gids.begin();i1!=gids.end();i1++){
		int ii1 = (*i1)-1;
		for(i2=i1,i2++;i2!=gids.end();i2++){
			int ii2 = (*i2)-1;
			//We check if the leemnt sids the biggest
			if(score>best_scores((ii1),(ii2))){
				best_scores((ii1),(ii2)) = score;
				names((ii1),(ii2)) = pid;
			}
		}
	}
}


//Getting the minimum common nodes.
void min_common_nodes(List& patterns,NumericVector& scores_pat,
                      IntegerMatrix& names,NumericMatrix& scores,
                      int num_objects){

	for(size_t i = 0; i < int(patterns.size()); i++){
		//We get the s4 object
		S4 pat = patterns[i];
		update_dist(pat,names,scores,i,scores_pat[i+num_objects]);
	}
// 	return List::create(Named("names")=names,
//                      Named("scores")=scores);
}

//Calculating the inital n and the inital score
// [[Rcpp::export]]
IntegerVector vecInitialisation(List& patterns,NumericVector& scores_pat,int num_objects){
	IntegerMatrix names(num_objects,num_objects) ;
	std::fill(names.begin(),names.end(),IntegerVector::get_na());
	NumericMatrix scores(num_objects,num_objects);
	std::fill(scores.begin(),scores.end(),0);

	min_common_nodes(patterns,scores_pat,
                  names,scores,num_objects);

	//Now we computes the values.
	IntegerVector n(patterns.size(),0);

	for(auto it=names.begin() ; it!=names.end() ; it++){
		n[*it]++;
	}
	return n;
}




//Function handling lower and upp_space.
std::set<VertexLatt> upper_space(latt& lat, std::map<VertexLatt,bool>& present,
                                 VertexLatt v, bool border, bool self){
	std::set<VertexLatt> res;
	std::deque<VertexLatt> queue;
	queue.push_back(v);

	if( (present[v]) & self){
		res.insert(v);
		return res;
	}
	VertexLatt cnode;
	while(true){
		if(queue.empty()){
			return res;
		}
		cnode = queue.front();

		graphTraitLatt::in_edge_iterator bi,ei;
		boost::tie(bi,ei) = boost::in_edges(cnode,lat);

		std::vector<VertexLatt> predecessors;
		for( ; (bi!=ei);bi++){
			predecessors.push_back(boost::source(*bi,lat));
		}
		if(border){
			for(auto ip = predecessors.begin();ip!=predecessors.end();ip++){
				if(present[*ip]){
					res.insert(*ip);

				}else{
					queue.push_front(*ip);
				}
			}
		}else{
			for(auto ip = predecessors.begin();ip!=predecessors.end();ip++){
				if(!present[*ip]) queue.push_front(*ip);
				res.insert(*ip);
			}
		}
		queue.pop_front();
	}
}

std::set<VertexLatt> lower_space(latt& lat, std::map<VertexLatt,bool>& present,
                                 VertexLatt v, bool border, bool self){
	std::set<VertexLatt> res;
	std::deque<VertexLatt> queue;
	queue.push_back(v);

	if( (present[v]) & self){
		res.insert(v);
		return res;
	}
	VertexLatt cnode;
	while(true){
		if(queue.empty()){
			return res;
		}
		cnode = queue.front();

		graphTraitLatt::out_edge_iterator bi,ei;
		boost::tie(bi,ei) = boost::out_edges(cnode,lat);

		std::vector<VertexLatt> successors;
		for( ; (bi!=ei);bi++){
			successors.push_back(boost::target(*bi,lat));
		}
		if(border){
			for(auto ip = successors.begin();ip!=successors.end();ip++){
				if(present[*ip]){
					res.insert(*ip);

				}else{
					queue.push_front(*ip);
				}
			}
		}else{
			for(auto ip = successors.begin();ip!=successors.end();ip++){
				if(!present[*ip]) queue.push_front(*ip);
				res.insert(*ip);
			}
		}
		queue.pop_front();
	}
}

//Update of the score matrix.
void update_score_removal(latt& l, NumericVector& current_scores,
                          NumericVector& raw_scores, NumericVector& n,
                          VertexLatt v,std::map<VertexLatt,bool>& present, int num_ir){

	auto ub = upper_space(l,present,v,true,false);
	int v_pos = v-num_ir;
	unsigned int nir = (unsigned int)num_ir;
	if(ub.size()>0){
		double best_score = 0;
		//VertexLatt best_node;
		int best_pos = 0;

		//We use the fact tha the graph is implemented as vector.
		for(auto it=ub.begin();it != ub.end();it++){
			if(raw_scores[*it]>best_score){
				//best_node = *it;
				best_pos = (*it)-num_ir;
				best_score = raw_scores[best_pos];
			}
		}

		//Updating the score directly affected to this node.
		current_scores[best_pos] = current_scores[best_pos] + n[v_pos] *
			(pow(raw_scores[v_pos],2)-pow(raw_scores[best_pos],2));
	}


	auto lb = upper_space(l,present,v,false,false);
	for(auto it=lb.begin();it != lb.end();it++){
		VertexLatt vl = *it;
		int l_pos = *it-num_ir;
		if((vl)<nir) continue;
		auto ubl = upper_space(l,present,vl,true,false);
		double best_score_b = 0;
		//VertexLatt best_node_b;
		int best_pos_b = 0;
		if(ubl.size()>0){
			for(auto itt=ubl.begin();itt != ubl.end();itt++){
				if(raw_scores[*itt]>best_score_b){
					//best_node_b = (*itt);
					best_pos_b = (*itt)-num_ir;
					best_score_b = raw_scores[best_pos_b];
				}

				if(raw_scores[v_pos]>raw_scores[best_pos_b]){
					//Updating the score to nodes affected to this node by the removals.
					current_scores[l_pos] = current_scores[l_pos] + n[l_pos] *
						(pow(raw_scores[v_pos],2)-pow(raw_scores[best_pos_b],2));
				}
			}
		}else{//The node cost of the new node is reinitialized.
			if(present[l_pos]){
				current_scores[l_pos] = current_scores[l_pos] + n[l_pos] *
					(pow(raw_scores[l_pos],2)-pow(raw_scores[v_pos],2));
			}

		}
	}
}

// Creation of the lattice
void fillLattice(DataFrame& df_nodes,DataFrame& df_edges, latt& l, std::set<VertexLatt>& items){
	NumericVector scores = df_nodes["scores"];

	IntegerVector names = df_nodes["name"];
	LogicalVector is_item = df_nodes["item"];

	//Adding the nodes.
	std::vector<VertexLatt> mapped_vec(scores.size());
	for(int i=0;i<scores.size();i++){
		VertexLatt nvert = boost::add_vertex(l);
		l[nvert].score = scores[i];
		//This gives an ID to identify the object.
		l[nvert].name =names[i]-1;
		mapped_vec[i] = nvert;
		if(is_item[i]) items.insert(nvert);
	}

	//Adding the edge
	if(df_edges.nrows()!=0){
		IntegerVector from = df_edges["from"];
		IntegerVector to = df_edges["to"];
		for(int i=0;i<from.size();i++){
			boost::add_edge(mapped_vec[from[i]-1],
                            mapped_vec[to[i]-1],l).first;
		}
	}
}

NumericVector initialCosts(latt& l,NumericVector& n,NumericVector& scores, int num_items){
	NumericVector icosts(n.size()-num_items,0);

	graphTraitLatt::vertex_iterator bv,ev;
	graphTraitLatt::in_edge_iterator bi,ei;
	unsigned int items = (unsigned int)num_items;
	boost::tie(bv,ev) = boost::vertices(l);
	for( ; bv != ev ; bv++){
		if((*bv)<num_items) continue;
		boost::tie(bi,ei) = boost::in_edges(*bv,l);

		double best_score = 0;
		if(bi!=ei){
			for(;bi != ei;bi++){
				if(scores[(source((*bi),l)-num_items)]>best_score){
					best_score = scores[(source((*bi),l)-num_items)];
				}
			}
			//We update the initial cost
			int cpos = (*bv)-num_items;
			icosts[cpos]=(pow(scores[cpos],2)-pow(best_score,2))*n[cpos];
		}else{
			int cpos = (*bv)-num_items;
			icosts[cpos]=pow(scores[cpos],2)*n[cpos];
		}
	}
	return icosts;

}

//Return the index of the minimum of a vezctor
int which_min(NumericVector& v){
	int best_pos = 0;
	double best_val = v[0];
	for(int i = 0;i<v.size();i++){
		if(NumericVector::is_na(v[i])) continue;
		if(v[i]<best_val){
			best_val = v[i];
			best_pos = i;
		}
	}
	return best_pos;
}


//Perfoming the lattice reduction.k is th enumber o
// [[Rcpp::export]]
Rcpp::IntegerVector reduceLatticeK_greedy(DataFrame df_nodes,DataFrame df_edges, NumericVector scores,
                                          NumericVector n,int k){
	latt l;
	std::set<VertexLatt> items;
	fillLattice(df_nodes,df_edges,l,items);
	NumericVector costs = initialCosts(l,n,scores,items.size());
	int nitems = items.size();
	//Initializing the presence map.
	std::map<VertexLatt,bool> present;

	graphTraitLatt::vertex_iterator bv,ev;
	boost::tie(bv,ev) = boost::vertices(l);
	for(;bv!=ev;bv++){
		present.emplace(*bv,true);
	}


	//We calculate the number of nodes to remove.
	int nrm = boost::num_vertices(l)-nitems-k;
	int crm = 0;

	//res
	Rcpp::IntegerVector res(nrm);


	while(crm<nrm){
		int pmin = which_min(costs);
		VertexLatt nlat = pmin+nitems;

		//We update the score.
		update_score_removal(l, costs, scores, n, nlat,
                           present, nitems);

		//Removing the elements.
		res[crm] = pmin+1;
		costs[pmin] = NumericVector::get_na();
		crm++;
	}
	return res;
}

// [[Rcpp::export]]
IntegerVector complementIdx(int num_object,IntegerVector rm){
	IntegerVector res(num_object-rm.size());
	std::set<int> removed(rm.begin(),rm.end());
	int cr = 0;
	for(int i = 0;i<num_object;i++){
		if(removed.find(i)==removed.end()){
			res[cr] = i;
			cr ++ ;
		}
	}
	return res;
}
