#ifndef TRIANGLES_LIST_H
#define TRIANGLES_LIST_H

#include <unordered_set>

#include <boost/array.hpp>


#include "common.h"
#include "mass_graph.h"

//TODO REPLACE BY BOOST MULTI INDEX VALUE.


//Class which store the infoprmation about the triangles
//present in a graph
class triangles_list
{
    public:
        //Emlpty copnstructor is always called.
        triangles_list();
        virtual ~triangles_list();

        //Adding a graph to the triangle list.
        void add_mass_graph(mass_graph& mg);

        void construct_mapping();

        short get_ab(short a,short b) const;
        short get_ac(short a,short c) const;

        //There could be multiples accessors for these values.
        std::vector<short> get_a_c(short a) const;
        std::vector<short> get_a_b(short a) const;
        //std::vector<short> get_c(short c);

        //void to_string();
    protected:
    private:
        std::set<boost::array<short, 3> > triplets;
        std::vector<boost::array<short, 3> > vtriplets;

        //The index used :
        std::unordered_map<short,std::vector < short > > idxa;
        std::unordered_map<short,std::vector < short > > idxb;
        std::unordered_map<short,std::vector < short > > idxc;

        std::map<std::pair<short,short>,short > idxab;
        std::map<std::pair<short,short>,short > idxac;

};

#endif // TRIANGLES_LIST_H
