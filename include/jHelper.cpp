#ifndef _J_HELPER_CPP_
#define _J_HELPER_CPP_

typedef double mt;
#define ex(x) << #x << ": " << left << setw(10) << x << "  "

template <typename T> int sgn(T val) {
        return (T(0) < val) - (val < T(0));
}


#endif //_J_HELPER_CPP_
