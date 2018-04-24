#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;

#include <boost/graph/depth_first_search.hpp>

#include "common.h"


//Definition of the graph used for reduction.
struct latt_info{
	//TODO change it to short value when graphml bug is fixed.
	short name;
	short sid;
	bool item;
	//double score;
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

	for(size_t i = 0; i < patterns.size(); i++){
		//Rcout<<i<<" ";
		//We get the s4 object
		S4 pat = patterns[i];
		update_dist(pat,names,scores,i,scores_pat[i+num_objects]);
	}
	//Rcout<<"...done...";
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
	//Rcout<<"min_common_nodes: ";
	min_common_nodes(patterns,scores_pat,
                  names,scores,num_objects);

	//Now we computes the values.
	IntegerVector n(patterns.size(),0);
	//Rcout <<"next loop";

	for(auto it=names.begin() ; it!=names.end() ; it++){
		if(IntegerMatrix::is_na(*it)) continue;
		//Rcout << " i:"<<*it<<" ";
		n[*it]++;
	}
	return n;
}




//Function handling lower and upp_space.
std::set<VertexLatt> upper_space(latt& lat, std::map<VertexLatt,bool>& present,
                                 VertexLatt v, bool border, bool self){
	Rcout <<" in upper_space ";
	std::set<VertexLatt> res;
	std::deque<VertexLatt> queue;
	queue.push_back(v);

	if( (present[v]) & self){
		res.insert(v);
		Rcout <<" out upper_space ";
		return res;
	}
	VertexLatt cnode;
	while(true){
		if(queue.empty()){
			Rcout <<" out upper_space ";
			return res;
		}
		cnode = queue.front();
		queue.pop_front();
		Rcout << " current:"<<cnode;

		graphTraitLatt::in_edge_iterator bi,ei;
		boost::tie(bi,ei) = boost::in_edges(cnode,lat);

		std::vector<VertexLatt> predecessors;
		for( ; (bi!=ei);bi++){
			if(!lat[boost::source(*bi,lat)].item)
			predecessors.push_back(boost::source(*bi,lat));
		}
		if(border){
			for(auto ip = predecessors.begin();ip!=predecessors.end();ip++){
				if(present[*ip]){
					res.insert(*ip);

				}else{
					Rcout<< " ins:"<<(*ip);
					queue.push_front(*ip);
				}
			}
		}else{
			for(auto ip = predecessors.begin();ip!=predecessors.end();ip++){
				if(!present[*ip]){
					Rcout<< " insb:"<<(*ip);
					queue.push_front(*ip);
				}
				res.insert(*ip);
			}
		}
		Rcout<<" pop:"<<queue.size()<<" head "<<queue.front()<<std::endl;
	}
}

std::set<VertexLatt> lower_space(latt& lat, std::map<VertexLatt,bool>& present,
                                 VertexLatt v, bool border, bool self){
	Rcout <<" in lower_space ";
	std::set<VertexLatt> res;
	std::deque<VertexLatt> queue;
	queue.push_back(v);

	if( (present[v]) & self){
		Rcout <<" out lower_space ";
		res.insert(v);
		return res;
	}
	VertexLatt cnode;
	while(true){
		if(queue.empty()){
			Rcout <<" out lower_space ";
			return res;
		}
		cnode = queue.front();
		queue.pop_front();

		graphTraitLatt::out_edge_iterator bi,ei;
		boost::tie(bi,ei) = boost::out_edges(cnode,lat);

		std::vector<VertexLatt> successors;
		for( ; (bi!=ei);bi++){
			if(!lat[boost::target(*bi,lat)].item)
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
		Rcout<<" pop:"<<queue.size()<<" head "<<queue.front()<<std::endl;
	}
}

//Update of the score matrix.
void update_score_removal(latt& l, NumericVector& current_scores,
                          NumericVector& raw_scores, NumericVector& n,
                          VertexLatt v,std::map<VertexLatt,bool>& present, int num_ir){

	Rcout << std::endl;
	auto ub = upper_space(l,present,v,true,false);
	int v_pos = v-num_ir;
	Rcout << " vnode:" << v_pos << " ";
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
		Rcout << " up:"<<best_pos <<" old_sc:"<<current_scores[best_pos]<<" ";
		current_scores[best_pos] = current_scores[best_pos] + n[v_pos] *
			(pow(raw_scores[v_pos],2)-pow(raw_scores[best_pos],2));
		Rcout << "ns:"<<current_scores[best_pos]<<"  ";
	}


	auto lb = lower_space(l,present,v,false,false);
	for(auto it=lb.begin();it != lb.end();it++){
		VertexLatt vl = *it;
		int l_pos = *it-num_ir;
		if((vl)<nir) continue;
		auto ubl = upper_space(l,present,vl,true,false);
		double best_score_b = 0;
		//VertexLatt best_node_b;
		Rcout<< " l:"<< l_pos <<" ";
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
					Rcout << " best_u_l:"<<best_pos_b <<" old_sc:"<<current_scores[l_pos]<<" ";
					current_scores[l_pos] = current_scores[l_pos] + n[l_pos] *
						(pow(raw_scores[v_pos],2)-pow(raw_scores[best_pos_b],2));
					Rcout << "ns:"<<current_scores[l_pos]<<"  ";
				}
			}
		}else{//The node cost of the new node is reinitialized.
			if(present[vl]){
				Rcout << " low:"<<l_pos <<" old_sc:"<<current_scores[l_pos]<<" ";
				current_scores[l_pos] = current_scores[l_pos] + n[l_pos] *
					(pow(raw_scores[l_pos],2)-pow(raw_scores[v_pos],2));
				Rcout << "ns:"<<current_scores[l_pos]<<"  ";
			}

		}
	}
	Rcout << std::endl;
}

// Creation of the lattice
std::map<int,VertexLatt> fillLattice(DataFrame& df_nodes,DataFrame& df_edges, latt& l, std::set<VertexLatt>& items){
	//NumericVector scores = df_nodes["scores"];

	IntegerVector names = df_nodes["name"];
	LogicalVector is_item = df_nodes["item"];
	IntegerVector sids = df_nodes["sid"];
	std::map<int,VertexLatt> nodes_ids;
	int cnode=0;

	//Adding the nodes.
	std::vector<VertexLatt> mapped_vec(names.size());
	for(int i=0;i<names.size();i++){
		VertexLatt nvert = boost::add_vertex(l);
		//l[nvert].score = scores[i];
		//This gives an ID to identify the object.
		l[nvert].name =names[i]-1;
		l[nvert].sid =sids[i]-1;
		mapped_vec[i] = nvert;
		if(is_item[i]){
			items.insert(nvert);
			l[nvert].item =true;
		}else{
			nodes_ids.emplace(sids[i],nvert);
			l[nvert].item =false;

		}
	}

	//Adding the edge
	if(df_edges.nrows()!=0){
		IntegerVector from = df_edges["from"];
		IntegerVector to = df_edges["to"];
		for(int i=0;i<from.size();i++){
			boost::add_edge(mapped_vec[from[i]-1],
                            mapped_vec[to[i]-1],l);
		}
	}
	return nodes_ids;
}

//NOW ALL THE VECTORS ARE STORED USING THEIR SIDS.
NumericVector initialCosts(latt& l,NumericVector& n,NumericVector& scores, int num_items){
	NumericVector icosts(n.size(),0);

	graphTraitLatt::vertex_iterator bv,ev;
	graphTraitLatt::in_edge_iterator bi,ei;
	//unsigned int nitems = (unsigned int)num_items;
	boost::tie(bv,ev) = boost::vertices(l);
	for( ; bv != ev ; bv++){
		if(l[*bv].item) continue;

		int sid = l[*bv].sid;
		boost::tie(bi,ei) = boost::in_edges(*bv,l);
		double best_score = 0;
		if(bi!=ei){
			for(;bi != ei;bi++){
				if(scores[l[source((*bi),l)].sid]>best_score){
					best_score = scores[l[source((*bi),l)].sid];
				}
			}
			//We update the initial cost
			icosts[sid]=(pow(scores[sid],2)-pow(best_score,2))*n[sid];
		}else{
			icosts[sid]=pow(scores[sid],2)*n[sid];
		}
	}
	return icosts;

}

//Return the index of the minimum of a vezctor
int which_min(NumericVector& v){
	int best_pos = 0;
	double best_val = 10000;
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
	std::map<int,VertexLatt> ids = fillLattice(df_nodes,df_edges,l,items);
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
		Rcout << "  vec: "<<costs.size()<<" ";
		for(auto it=costs.begin();it!=costs.end();it++){
			if(Rcpp::NumericVector::is_na(*it)){
				Rcout<<"NA_";
			}else{
				Rcout << (*it) <<"_";
			}
		}
		Rcout << std::endl;

		VertexLatt nlat = ids[pmin];
		Rcout<<"crm: "<<crm<<" rm: "<<pmin<<"aa";

		//We update the score.
		update_score_removal(l, costs, scores, n, nlat,
                           present, nitems);

		//Removing the elements.
		res[crm] = pmin+1;
		costs[pmin] = NA_REAL;
		present[ids[pmin]] = false;
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
