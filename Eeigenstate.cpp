#ifndef _EEIGENSTATE_CPP_
#define _EEIGENSTATE_CPP_

typedef complex<double> phiRangeType;

class Eeigenstate {
    public:
        fct<phiRangeType> phi;
        mt E;
        mt k;
        int band;
};

#endif // _EEIGENSTATE_CPP_

