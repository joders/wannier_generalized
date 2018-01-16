
#ifndef _STREAMERS_CPP_
#define _STREAMERS_CPP_

#include<jHelper.cpp>
#include<fct_classes/fct_class.cpp>
#include"Eeigenstate.cpp"
#include"iterators.cpp"

using namespace std;

template<typename streamableT>
std::ostream& operator<<(std::ostream& o, vector<streamableT> const &par) {
    for(int c=0;c<par.size();c++)
        o << par[c] << '\t';
    return o;
}
std::ostream& operator<<(std::ostream& o, fct<complex<mt> > const &par) {
    if(par.size()>0) {
        for(unsigned c=0;c<par.size();c++)
            o << par.x[c] << '\t' << par.f[c].real() << endl;
//        mt start=par.x.front();
//        mt end  =par.x.back();
//    //    if(par.second_deriv_for_fcts_spline==NULL)
//    //        cout << "streaming non-spline-interpolated fct" << endl;
//        for(unsigned c=0;c<=1000;c++) {
//            mt x=start+(end-start)*(mt)c/1000;
//            o << x << '\t' << par(x).real() /*<< '\t' << par(x).imag()*/ << endl;
//        }
    }
    return o;
}
std::ostream& operator<<(std::ostream& o, fct<mt> const &par) {
    if(par.size()>0) {
        mt start=par.x.front();
        mt end  =par.x.back();
    //    if(par.second_deriv_for_fcts_spline==NULL)
    //        cout << "streaming non-spline-interpolated fct" << endl;
        for(unsigned c=0;c<=1000;c++) {
            mt x=start+(end-start)*(mt)c/1000;
            o << x << '\t' << par(x) << endl;
        }
    }
    return o;
}
template<typename fctRangeType>
struct outputValue {
    mt x; 
    std::ostream *o;
    void operator()(fct<fctRangeType> const &par) {
//        if(par.second_deriv_for_fcts_spline==NULL)
//            cout << "streaming non-spline-interpolated fct" << endl;
        *o << par(x) << '\t';
    }
    void operator()(Eeigenstate const &par) {
        *o << par.phi(x).real() << '\t';
    }
};
template<>
void outputValue<complex<mt> >::operator()(fct<complex<mt> > const &par) {
//    if(par.second_deriv_for_fcts_spline==NULL)
//        cout << "streaming non-spline-interpolated fct" << endl;
    *o << par(x).real() << '\t' /*<< par(x).imag() << '\t'*/;
}

template<typename fctRangeType>
std::ostream& operator<<(std::ostream& o, vector<fct<fctRangeType> > const &par) {
    if(par[0].size()<1000) {
        for(unsigned c=0;c<par[0].size();c++) {
            mt x=par[0].x[c];
            o << x << '\t';
            outputValue<fctRangeType> ov;
            ov.x=x;
            ov.o=&o;
            for_each(par,ov);
            o << endl;
/*            for(unsigned d=0;d<par.size();d++)
                o << par[d](x) << '\t';
            o << endl;*/
        }
    } else {
        mt start=par.front().x.front();
        mt end  =par.front().x.back();
        for(unsigned c=0;c<=1000;c++) {
            mt x=start+(end-start)*(mt)c/1000;
            o << x << '\t';
            outputValue<fctRangeType> ov;
            ov.x=x;
            ov.o=&o;
            for_each(par,ov);
            o << endl;
        }
    }
    return o;
}
std::ostream& operator<<(std::ostream& o, vector<Eeigenstate> const &Ee) {
    mt start=Ee.front().phi.x.front();
    mt end  =Ee.front().phi.x.back();
    for(unsigned c=0;c<=1000;c++) {
        mt x=start+(end-start)*(mt)c/1000;
        o << x << '\t';
        outputValue<phiRangeType> ov;
        ov.x=x;
        ov.o=&o;
        for_each(Ee,ov);
        o << endl;
    }
    return o;
}
std::ostream& operator<<(std::ostream& o, vector<vector<Eeigenstate> > const &Ee) {
    mt start=Ee.front().front().phi.x.front();
    mt end  =Ee.front().front().phi.x.back();
    for(unsigned c=0;c<=1000;c++) {
        mt x=start+(end-start)*(mt)c/1000;
        o << x << '\t';
        outputValue<phiRangeType> ov;
        ov.x=x;
        ov.o=&o;
        for_each(Ee,ov);
        o << endl;
    }
    return o;
}

#endif // _STREAMERS_CPP_

