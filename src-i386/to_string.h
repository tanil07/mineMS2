#ifndef UTIL_F_H_INCLUDED
#define UTIL_F_H_INCLUDED

#include <string>
#include <sstream>

namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}
#endif // UTIL_F_H_INCLUDED
