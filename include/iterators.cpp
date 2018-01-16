
#ifndef _ITERATORS_CPP_
#define _ITERATORS_CPP_

#include<vector>

using namespace std;

template<class Function, class elementT>
Function for_each(std::vector<elementT> &vec, Function fn) {
    for(typename vector<elementT>::iterator it=vec.begin();it!=vec.end();it++)
        for_each(*it,fn);
    return fn;
}
template<class Function, class elementT>
Function for_each(std::vector<elementT> const &vec, Function fn) {
    for(typename vector<elementT>::const_iterator it=vec.begin();it!=vec.end();it++)
        for_each(*it,fn);
    return fn;
}
template<class Function, class elementT>
Function for_each(elementT &vec, Function fn) { 
    fn(vec); 
    return(fn);
}

#endif

