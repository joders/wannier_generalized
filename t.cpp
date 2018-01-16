#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<cfloat>
#include<fstream>
#include<iostream>
#include<iomanip>
#include<set>
#include<map>

#include<armadillo>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>

#include<boost/math/special_functions/hermite.hpp>

#include<sys/stat.h>
#include<argp.h>

#include<ODE_integration.cpp>
#include<utils.cpp>
#include<fct_classes/fct_class.cpp>
#include<fct_classes/fct_class_functions.cpp>
#include<jHelper.cpp>
#include<findMultipleBunchedRoots.cpp>
#include<diag_Householder_QL/diag_hermitian_Householder_QL.cpp>
#include<iterators.cpp>

//#include"../Eeigenstate.cpp"
#include"streamer.cpp"

//#define PFSI // whether to call prepare for spline interpolation 
//(in my experience this only makes it slower and not more accurate)


using namespace std;
using namespace arma;

string dirName;

class fileStdoutStream:public ostream {
    public:
        ofstream file;
        fileStdoutStream(){}
        void open(const char* filename){
            file.open(filename);
        }
        ~fileStdoutStream(){
            file.close();
        }
        template<typename T>
        fileStdoutStream& operator<<(const T & arg) {
            std::cout << arg;
            file << arg;
            return *this;
        }
        typedef std::ostream& (*stream_function)(std::ostream&);
        fileStdoutStream& operator<<(stream_function func) {
            func(std::cout);
            func(file);
            return *this;
        }
        fileStdoutStream& flush() {
            cout.flush();
            file.flush();
            return *this;
        }

};
fileStdoutStream fout;

class physical_parameters { 
    public:
        mt hbar;
        mt m;

        mt a;
        mt lambda;
        mt k_L; // lattice moment of laser
        mt E_r;

        mt s_bound; 

        mt E; // this value is solved for
        mt k;
        mt a_0;
        mt a_s;
        mt g;

        mt ld; // lattice depth
        mt inhS;
        mt alpha;
        mt inhShift;

        mt borderTrapHeight;
        mt inhLatRatio;

//        unsigned nrSitesLat;
//        unsigned nrSitesTrap;
} pp;

    const bool findMultipleBunchedRoots::DEBUG::PROPAGATE_CALLS=false;
    const bool findMultipleBunchedRoots::DEBUG::SEARCH_CALLS=false;
    const bool findMultipleBunchedRoots::DEBUG::FOUND_ROOT=false;

unsigned kmax;

//mt Vlat(mt x, physical_parameters &pp) {
//    mt sinres=sin(pp.k_L*x);
//    return(pp.ld*pp.E_r*(-sinres*sinres+.5/*-1.*/));
//}
//complex<mt> Vinh(mt x/*, physical_parameters &pp*/) {
//    return(pp.ts*(x-pp.inhShift)*(x-pp.inhShift));
//}
mt Vlat(mt x, physical_parameters &pp) {
    return(pp.ld*pp.E_r*cos(2*pp.k_L*x)/2);
}
complex<mt> Vinh(mt x/*, physical_parameters &pp*/) {
//    return(pp.inhS*pp.E_r*cos(2*pp.alpha*pp.k_L*(x-pp.inhShift))/2);
    return(pp.inhS*(x-pp.inhShift)*(x-pp.inhShift));
}
mt Vlat(mt x) {
    return Vlat(x,pp);
}
complex<mt> V(mt x, physical_parameters &pp) {
    return Vlat(x,pp)+Vinh(x);
}

complex<mt> testF(mt x) {
    return(exp(-(x+14)*(x+14)/8.));

}

void ODE_derivs(const double x, double *y, double *dydx, physical_parameters &pp)
{    
    dydx[0]=y[1];
    dydx[1]=(2.*pp.m*(Vlat(x,pp)-pp.E)/pp.hbar/pp.hbar + pp.k*pp.k)*y[0]+
                                                           2.*pp.k *y[3];
    dydx[2]=y[3];
    dydx[3]=(2.*pp.m*(Vlat(x,pp)-pp.E)/pp.hbar/pp.hbar + pp.k*pp.k)*y[2]-
                                                           2.*pp.k *y[1];
}

ODE<physical_parameters> *D1=NULL;
ODE<physical_parameters> *D2=NULL;

/*mt prek=44444;
unsigned kc=0;*/

mt odeStepWidth;
mt propStart;
mt propEnd;
void dglSolve(ODE<physical_parameters> *&D , unsigned eqncount, mt const xinit, mt *const &yinit) {

/*    if(prek!=pp.k) {
        fout ex(pp.k) << ' ' ex(kc) << endl;
        prek=pp.k;
        kc=0;
    }
    kc++;*/

    if((propStart>xinit)||(propEnd<xinit)) {
        fout << "xinit must lie between start and end of propagation" << endl;
        abort();
    }

    if(propStart!=xinit) {
        if(D!=NULL)
            delete D;
        D=new ODE<class physical_parameters>(eqncount,1E-9);
        D->dy_dt_fct=&ODE_derivs;
        D->P.dt_save_step=DBL_MAX;
        D->P.max_step_size_for_dense_output=DBL_MAX;
        D->t[D->P.t_length]=xinit; //-pp.s_bound;   CHECK THIS
        D->P.t_end=propStart;

        D->A=pp;

        kmax=100;
        D->set_initial_y_vec(yinit);

        odeint(D);

        yinit[0]=D->y[D->P.t_length-1][0];
        yinit[1]=D->y[D->P.t_length-1][1];
        yinit[2]=D->y[D->P.t_length-1][2];
        yinit[3]=D->y[D->P.t_length-1][3];
    }

    if(D!=NULL)
        delete D;
    D=new ODE<class physical_parameters>(eqncount,1E-9);
    D->dy_dt_fct=&ODE_derivs;
    D->P.dt_save_step=odeStepWidth;
    D->P.max_step_size_for_dense_output=odeStepWidth;
    D->t[D->P.t_length]=propStart; //-pp.s_bound;   CHECK THIS
    D->P.t_end=propEnd;

    D->A=pp;

    kmax=100;
    D->set_initial_y_vec(yinit);

    odeint(D);  
}

fct<complex<mt> > phiEndPoints;
map<mt,mt> calculatedValues;
mt lastk=DBL_MAX;
mt propagateFunction(mt par) {
    pp.E=par;

    try { 
        mt cv=calculatedValues.at(par); 
        if(lastk!=pp.k) {
            fout << "ERROR " << __LINE__ << endl;
            abort();
        }
        //fout << par << " present: " << cv << endl; 
        return(cv);} 
    catch(out_of_range) {

    //fout << par << " not present: "; 
    

    //mt yinit1[]={1.,1E-12,1E-12,1E-12};
    mt xinit1=0.;
    mt yinit1[]={1.,0.,0.,0.};
    dglSolve(D1,4,xinit1,yinit1);

    mt I1 =D1->y[D1->P.t_length-1][2];
    mt R1d=D1->y[D1->P.t_length-1][1];

    //mt yinit2[]={1E-12,1E-12,1E-12,1.};
    mt xinit2=0.;
    mt yinit2[]={0.,0.,0.,1.};
    dglSolve(D2,4,xinit2,yinit2);

    mt I2 =D2->y[D2->P.t_length-1][2];
    mt R2d=D2->y[D2->P.t_length-1][1];

    //fout << pp.k << endl;

    mt res=(I1*R2d-I2*R1d/*-0.001*/);
    calculatedValues[par]=res;
    //fout << res << endl;
//    phiEndPoints.x.push_back(par);
//    phiEndPoints.f.push_back(I1*R2d-I2*R1d-0.001);

//    fout << "propagateFunction[" << par << "]=" << (I1*R2d-I2*R1d-0.001) << endl;

    lastk=pp.k;

    return(res); // correct this
    }
}
mt invPropagateFunction(mt par) {
    return(-propagateFunction(par));
}
mt propagateFunction(mt par, void*) {
    return propagateFunction(par);
}
mt invPropagateFunction(mt par, void*) {
    return(invPropagateFunction(par));
}

/*template<typename t>
bool smallerThan(t const a, t const b) {return a<b;}
template<typename t>
bool biggerThan(t const a, t const b) {return a>b;}
template<typename t>
vector<unsigned> findAbsExtrema(vector<t> const series, bool (*compare)(mt const &a, mt const &b), vector<unsigned> const*const indices=NULL, unsigned start=0) {
    vector<unsigned> extremaSites;
    if(indices==NULL) {
        if(start<=series.size()-2)
            if(compare(fabs(series[start]),fabs(series[start+1])))
                extremaSites.push_back(start);
        for(unsigned c=start+1;c<=series.size()-2;c++)
            if(compare(fabs(series[c]),fabs(series[c-1])) && 
               compare(fabs(series[c]),fabs(series[c+1]))    )
                extremaSites.push_back(c);
        if(start<series.size())
            if(compare(fabs(series[series.size()-1]),fabs(series[series.size()-2])))
                extremaSites.push_back(series.size()-1);
    } else {
        if(start<=indices->size()-2)
            if(compare(fabs(series[indices->at(start)]),fabs(series[indices->at(start+1)])))
                extremaSites.push_back(start);
        for(unsigned c=start+1;c<=indices->size()-2;c++)
            if(compare(fabs(series[indices->at(c)]),fabs(series[indices->at(c-1)])) && 
               compare(fabs(series[indices->at(c)]),fabs(series[indices->at(c+1)]))    )
                extremaSites.push_back(c);
        if(start<indices->size())
            if(compare(fabs(series[indices->at(indices->size()-1)]),fabs(series[indices->at(indices->size()-2)])))
                extremaSites.push_back(indices->size()-1);
    }
    return(extremaSites);
}
void cutOffEndDivergence(Eeigenstate &Ee) {
    vector<int> absMaxima=findAbsExtrema(Ee.phi.f,biggerThan);
    fout ex(absMaxima);
    for(int e=0;e<absMaxima.size();e++)
        fout << Ee.phi.f[absMaxima[e]] << ' ';
    fout << endl;
    vector<int> maxOfMax=findAbsExtrema(Ee.phi.f,biggerThan,&absMaxima);
    fout ex(maxOfMax);
    vector<int> minMaxFollowingMaxOfMax=findAbsExtrema(Ee.phi.f,smallerThan,&absMaxima,maxOfMax[0]+1);
    fout ex(minMaxFollowingMaxOfMax);
    if(minMaxFollowingMaxOfMax.size()>0) {
        vector<int> endPoint=findAbsExtrema(Ee.phi.f,smallerThan,NULL,absMaxima[minMaxFollowingMaxOfMax[0]]+1);
        fout ex(endPoint);
        if(endPoint.size()>0)
            for(int e=endPoint[0];e<Ee.phi.f.size();e++)
                Ee.phi.f[e]=0.;
    }
}*/

void outputPropagatedWavefunctionEndValueOverE(mt const Estart, mt const Eend, mt const step) {
    ofstream file((dirName+"/propagatedWavefunctionEndValueOverE.gp").c_str());
    odeStepWidth=DBL_MAX;
    for(mt d=Estart;d<=Eend;d+=step) {
        pp.k=0;
        file << d << '\t' << propagateFunction(d) << ' ';
        calculatedValues.clear();
        pp.k=PI/2/pp.a;
        file << propagateFunction(d) << ' ';
        calculatedValues.clear();
        pp.k=PI/pp.a;
        file << propagateFunction(d) << endl;
        calculatedValues.clear();
    }
    file.close();
}

unsigned statesPerBand0;
unsigned statesPerBand;
unsigned bandToIndex0(unsigned band, unsigned stateIndexInBand) {
    unsigned ind=statesPerBand0*band+stateIndexInBand;
/*    if(Eeigenstates[ind].band!=band)
        fout << "ERROR" << __FILE__ << __LINE__;*/
    return ind;
}
unsigned bandToIndex(unsigned band, unsigned stateIndexInBand) {
    unsigned ind=statesPerBand*band+stateIndexInBand;
/*    if(Eeigenstates[ind].band!=band)
        fout << "ERROR" << __FILE__ << __LINE__;*/
    return ind;
}

/*struct outputEvalueAndBand {
    ofstream *o;
    void operator()(Eeigenstate const &par) {
        *o << setprecision(20) << "k=" << par.k << "  E=" << par.E << "  band=" << par.band << endl;
    }
};*/
/*void outputEvaluesAndBands(vector<Eeigenstate> const &Es) {
    ofstream file("output/H0eigenvalues");
    file << "H0eigenvalues: " << endl;
    outputEvalueAndBand oeab;
    oeab.o=&file;
    for_each(Es,oeab);
    file.close();
}*/

/*void normalizeFct(fct &f) {  // substitute with fct::normalize when fct::normalize works with romberg
    mt norm=romberg_integral(f*f,-pp.s_bound,+pp.s_bound,1E-9);
    f=f/sqrt(norm);
}*/
/*/void normalizeComplexFct(fct<complex<double> > &F) {  // substitute with fct::normalize when fct::normalize works with romberg
    complex<double> norm=romberg_integral(F*F,-pp.s_bound,+pp.s_bound,1E-9);  // CHECK THIS
    F/=sqrt(fabs(norm));                                                  // CHECK THIS
}*/

void normalizeFct(fct<complex<mt> > &par) {
    par.normalize_with_L2_norm();      // CHECK THIS
}

/*Mat<cx_double> calculateOverlapMatrix(vector<Eeigenstate> const &Es) {
    Mat<cx_double> overlap(Es.size(),Es.size());
    for(unsigned c=0;c<Es.size();c++)
        for(unsigned d=0;d<Es.size();d++)
//            if(d<c)
//                overlap(c,d)=overlap(d,c);      // OPTIMIZATION
//            else 
            {
                //fct norm=Es[0][c].rePsi*Es[0][c].rePsi+Es[0][c].imPsi*Es[0][c].imPsi;
                overlap(c,d)=romberg_integral(Es[c].phi.conj(), Es[d].phi, Es[c].phi.x[0], +Es[c].phi.x[Es[c].phi.size()-1], 1E-6);
            }
    return overlap;
}*/

Mat<cx_double> calculateOverlapMatrix(vector<fct<complex<double> > > const &Es) { // doesn't make sense overlap cannot be calculated using just one site
    Mat<cx_double> overlap(Es.size(),Es.size());                          // look at analyics
    for(unsigned c=0;c<Es.size();c++)
        for(unsigned d=0;d<Es.size();d++)
            /*if(d<c)
                overlap(c,d)=overlap(d,c);      // OPTIMIZATION
            else */{
                //fct norm=Es[0][c].rePsi*Es[0][c].rePsi+Es[0][c].imPsi*Es[0][c].imPsi;
                overlap(c,d)=romberg_integral(Es[c].conj(), Es[d], Es[c].x[0], +Es[c].x[Es[c].size()-1], 1E-6);
            }
    return overlap;
}

Mat<cx_double> calculateEigenstatesBlockwise(Mat<cx_double> const symmetricMatrix, Col<mt> &eigenvals, unsigned statesPerBand) {
    Mat<cx_double> eigenstates(symmetricMatrix.n_rows,symmetricMatrix.n_rows);
    eigenstates.zeros();
    for(unsigned c=0;c<symmetricMatrix.n_rows/statesPerBand;c++) {
        Col<mt> eigenvalsPart;
        Mat<cx_double> eigenstatesPart;
        Mat<cx_double> symmetricMatrixPart=symmetricMatrix(span(c*statesPerBand,((c+1)*statesPerBand)-1),span(c*statesPerBand,((c+1)*statesPerBand)-1));

        eig_sym(eigenvalsPart,eigenstatesPart,symmetricMatrixPart);

        eigenvals  (span(c*statesPerBand,((c+1)*statesPerBand)-1)                                              )=eigenvalsPart;
        eigenstates(span(c*statesPerBand,((c+1)*statesPerBand)-1),span(c*statesPerBand,((c+1)*statesPerBand)-1))=eigenstatesPart;
    }
    return(eigenstates);
}

Mat<mt> createDiagonalMatrixFromVector(Col<mt> const eigenvalues) {
    Mat<mt> M(eigenvalues.n_elem,eigenvalues.n_elem);
    M.zeros(); 
    M.diag()=eigenvalues;
    return(M);
}

void addPhaseToLocalizeOnRealAxis(fct<complex<mt> > &f) {
    mt max=0;
    unsigned ind;
    for(unsigned d=0;d<f.size();d++)
        if(fabs(imag(f.f[d]))>max) {            //fabs?!
            max=fabs(imag(f.f[d]));               // fabs?!
            ind=d;
        }
    //fout << c << ' ' << f.f[ind] << ' ' << ind << endl;
    complex<mt> factor=fabs(f.f[ind])/f.f[ind];
    for(unsigned d=0;d<f.size();d++)
        f.f[d]*=factor;
}
/*mt calculateVariance(fct wan, mt expectation) {
    mt variance=0;
    for(unsigned c=1;c<=wan.size()-2;c++)
        variance+=(wan.x[c]-expectation)*(wan.x[c]-expectation)*wan.f[c]*wan.f[c]*((wan.x[c+1]-wan.x[c-1])/2);
    return(variance);
}

vector<mt> calculateVariances(vector<fct> wans) {
    vector<mt> variances;
    for(unsigned c=0;c<wans.size();c++) {
        variances.push_back(calculateVariance(wans[c], (mt)c-wans.size()/2));
    }
    return(variances);
}*/

/// the same integration procedure for a product of two functions
template<typename fctRangeType, typename Fret>
fctRangeType trapez_int_refined_f(fct<fctRangeType> func1, Fret (*f)(mt), fctDomainType x_lo, fctDomainType x_hi, int n, fctRangeType &s) {
    fctDomainType tnm, del, x;
    fctRangeType sum;
	int it,j;
	if (n == 1)
		return (s=0.5*(x_hi-x_lo)*( (func1(x_lo)*f(x_lo)) +(func1(x_hi)*f(x_hi))  ));
	else {
		for (it=1,j=1;j<n-1;j++) 
            it <<= 1;
		tnm=it;
		del=(x_hi-x_lo)/tnm;
		x=x_lo+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) 
            sum += (func1(x)*f(x));
		s=0.5*(s+(x_hi-x_lo)*sum/tnm);
		return s;
	}
}

/// the same integration procedure for a product of two functions
template<typename fctRangeType, typename Fret>
fctRangeType trapez_int_refined_f(fct<fctRangeType> func1, fct<fctRangeType> func2, Fret (*f)(mt), fctDomainType x_lo, fctDomainType x_hi, int n, fctRangeType &s) {
    fctDomainType tnm, del, x;
    fctRangeType sum;
	int it,j;
	if (n == 1)
		return (s=0.5*(x_hi-x_lo)*( (func1(x_lo)*func2(x_lo)*f(x_lo)) +(func1(x_hi)*func2(x_hi)*f(x_hi))  ));
	else {
		for (it=1,j=1;j<n-1;j++) 
            it <<= 1;
		tnm=it;
		del=(x_hi-x_lo)/tnm;
		x=x_lo+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) 
            sum += (func1(x)*func2(x)*f(x));
		s=0.5*(s+(x_hi-x_lo)*sum/tnm);
		return s;
	}
}

/// Romberg integration routine taken from numerical recipes in C++, the file qromb.c
template<typename fctRangeType, typename Fret>
fctRangeType romberg_integral_f(fct<fctRangeType> func1,Fret (*f)(mt), double romberg_precision=romberg_integration_precision)
{
    fctDomainType x_lo=func1.x[0], x_hi=func1.x[func1.size()-1];

    if(x_lo<x_hi) {
        fctRangeType s_in_trapez_int_refined; ///initially declared as a static variable in trapez_int_refined - pass this ba reference
        const int JMAX=30, JMAXP=JMAX+1;
        const int K=6;      /// determines quality: if too low the method seems to sometime have thought that it has converged, where it hasn't
        const double EPS=romberg_precision;
        //printf("romberg_integration_precision=%g\n",EPS);
        fctRangeType ss,dss;
        fctRangeType s[JMAXP+1];
        double h[JMAXP+1];
        bool lastDssSmall=false;

        int j;
        h[1]=1.0;
        for (j=1;j<=JMAX;j++) {
            s[j]=trapez_int_refined_f(func1, f, x_lo, x_hi,j, s_in_trapez_int_refined);
            //printf("RI: s[%d]=%1.12g\n",j,s[j]);
            if (j >= K) {
    //            fout ex(s[j-4]) ex(s[j-3]) ex(s[j-2]) ex(s[j-1]) ex(s[j]) << endl;
                polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
                if (fabs(dss) < EPS*fabs(ss) ||  fabs(ss) < 1E-9 ) {
                    if(lastDssSmall)
                        return ss;
                    else
                        lastDssSmall=true;
                } else
                    lastDssSmall=false;
            }
            s[j+1]=s[j];
            h[j+1]=0.25*h[j];
    // 		printf ("qromb: ss = %le, dss = %le, EPS = %le, j=%d\n",ss,dss,EPS,j); 
        }
        printf("Error: Too many steps in routine qromb\n\n");
        return 0.0;
    } else
        return 0.;
}

/// Romberg integration routine taken from numerical recipes in C++, the file qromb.c
template<typename fctRangeType, typename Fret>
fctRangeType romberg_integral_f(fct<fctRangeType> func1,fct<fctRangeType> func2, Fret (*f)(mt), double romberg_precision=romberg_integration_precision)
{
    if(func1.size()==0)
        return 0;
    if(func2.size()==0)
        return 0;
    fctDomainType x_lo, x_hi;
    if(func1.x[0]<func2.x[0])
        x_lo=func2.x[0];
    else
        x_lo=func1.x[0];

    if(func1.x[func1.size()-1]<func2.x[func2.size()-1])
        x_hi=func1.x[func1.size()-1];
    else
        x_hi=func2.x[func2.size()-1];

    if(x_lo<x_hi) {
        fctRangeType s_in_trapez_int_refined; ///initially declared as a static variable in trapez_int_refined - pass this ba reference
        const int JMAX=30, JMAXP=JMAX+1;
        const int K=6;      /// determines quality: if too low the method seems to sometime have thought that it has converged, where it hasn't
        const double EPS=romberg_precision;
        //printf("romberg_integration_precision=%g\n",EPS);
        fctRangeType ss,dss;
        fctRangeType s[JMAXP+1];
//        fctRangeType ssHistory[JMAXP+1];
        double h[JMAXP+1];
        bool lastDssSmall=false;

        int j;
        h[1]=1.0;
        for (j=1;j<=JMAX;j++) {
            s[j]=trapez_int_refined_f(func1, func2, f, x_lo, x_hi,j, s_in_trapez_int_refined);
            //printf("RI: s[%d]=%1.12g\n",j,s[j]);
            if (j >= K) {
    //            fout ex(s[j-4]) ex(s[j-3]) ex(s[j-2]) ex(s[j-1]) ex(s[j]) << endl;
                polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
//                ssHistory[j]=ss;
                if (fabs(dss) < EPS*fabs(ss) ||  fabs(ss) < 1E-9 ) {
                    if(lastDssSmall)
                        return ss;
                    else
                        lastDssSmall=true;
                } else
                    lastDssSmall=false;
            }
            s[j+1]=s[j];
            h[j+1]=0.25*h[j];
    // 		printf ("qromb: ss = %le, dss = %le, EPS = %le, j=%d\n",ss,dss,EPS,j); 
        }
        printf("Error: Too many steps in routine qromb\n\n");
        return 0.0;
    } else
        return 0.;
}

/*void handler (const char * reason, const char * file, int line, int gsl_errno) {
    return ;
}*/

struct blochFindMultipleRoots {
    mt intervalStart;
    mt intervalEnd;
//    mt closestRootSpacing;
    mt rootValueAccuracy;
    mt takenZero;
//    mt gapFactorToNewBunch; 
//    mt iterationsBeforeStepIncrementation;
//    mt incrementationFactor;
//    mt stepsBetweenRootsInBunch;  // amount of steps between roots when roots were to be equidistantly spaced
    mt (*rf)(mt);
    mt (*mf)(mt,void*);
    set<mt> allRoots;
    set<mt> allMins;
    mt intervalForMinSearch;

    void searchDownUp() {
        mt s=intervalStart,m=s+intervalForMinSearch/2,e=s+intervalForMinSearch;
        while(e<intervalEnd) {
            findMin(s,e,rootValueAccuracy);
            //fout ex(s) ex(e) << endl;
            s=m;
            m=e;
            e+=intervalForMinSearch/2;
        }
        e=intervalEnd;
        m=(s+e)/2;
        findMin(s,e,rootValueAccuracy);
        //fout ex(s) ex(e) << endl;
        /*for(set<mt>::iterator i=allMins.begin();i!=allMins.end();i++) {
            fout << *i << ' ' << rf(*i) << endl;
        }
        fout << endl;*/
        map<mt,mt> phiEndSamples=calculatedValues;
        map<mt,mt>::iterator iprev=phiEndSamples.begin();
        map<mt,mt>::iterator i=++phiEndSamples.begin();
        map<mt,mt>::iterator ipast=++(++phiEndSamples.begin());
        while(i!=phiEndSamples.end()) {
            //fout << i->first << ' ' << i->second << endl;
            if(sgn(iprev->second)!=sgn(i->second)) {
                mt foundRoot=Ridder_find_val(rf,0.,iprev->first,i->first,rootValueAccuracy);
                allRoots.insert(foundRoot);
            } else
                if(ipast!=phiEndSamples.end())
                    if(sgn(i->second)==sgn(ipast->second))
                        if(fabs(i->second)<fabs(iprev->second))
                            if(fabs(i->second)<fabs(ipast->second))
                                if(fabs(i->second)<takenZero) {
                                    allRoots.insert(i->first*(1.-1E-8));
                                    allRoots.insert(i->first*(1.+1E-8));
                                }

            iprev++; i++; ipast++;
        }
        //fout.flush();
        for(map<mt,mt>::iterator i=++calculatedValues.begin();i!=calculatedValues.end();i++) {
            //fout << i->first << ' ' << i->second << endl;
        }

        bool found=true;
        while(found) {
            found=false;
            set<mt>::iterator iprev=     allRoots.begin();
            set<mt>::iterator i=       ++allRoots.begin();
            set<mt>::iterator ipast;
            if(i==allRoots.end())
                ipast=allRoots.end();
            else
                ipast=++(++allRoots.begin());
            while(ipast!=allRoots.end()) {
                if(fabs(*iprev-*ipast)<intervalForMinSearch/100) {
                    found=true;
                    allRoots.erase(i);
                    break;
                }
                iprev++; i++; ipast++;
            }
        }

        //fout << calculatedValues.size() << endl;
        calculatedValues.clear();
    }

    void findMin(mt start, mt end, mt precision) {
        int status;
        int iter = 0, max_iter = 100;
        const gsl_min_fminimizer_type *T;
        gsl_min_fminimizer *s;
        double m = (start+end)/2.0;
        double a = start, b = end;
        gsl_function F;
      
        F.function = mf;
        F.params = 0;
      
        T = gsl_min_fminimizer_brent;
        s = gsl_min_fminimizer_alloc (T);
        mt setErr=gsl_min_fminimizer_set (s, &F, m, a, b);
        //fout ex(setErr) << endl;
        if(setErr==4)
            return;
      
      /*  printf ("using %s method\n",
                gsl_min_fminimizer_name (s));
      
        printf ("%5s [%9s, %9s] %9s %10s %9s\n",
                "iter", "lower", "upper", "min",
                "err", "err(est)");
      
        printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
                iter, a, b,
                m, m - m_expected, b - a);*/
      
        do {
            iter++;
            status = gsl_min_fminimizer_iterate (s);
      
            m = gsl_min_fminimizer_x_minimum (s);
            a = gsl_min_fminimizer_x_lower (s);
            b = gsl_min_fminimizer_x_upper (s);
      
            status 
              = gsl_min_test_interval (a, b, precision, 0.0);
      
      /*      if (status == GSL_SUCCESS)
              printf ("Converged:\n");
      
            printf ("%5d [%.7f, %.7f] "
                    "%.7f %+.7f %.7f\n",
                    iter, a, b,
                    m, m - m_expected, b - a);*/
        } while (status == GSL_CONTINUE && iter < max_iter);
      
        gsl_min_fminimizer_free (s);

//        if(m<1*prec)
//            m=-1E-6;
        allMins.insert(m);
    }

};

mt constantOp(mt par) {
    return(1.);
}
mt unityOp(mt par) {
    return par;
}
mt absUnityOp(mt par) {
    return fabs(par);
}
mt squareOp(mt par) {
    return par*par;
}

template<typename fctRangeType>
string getFileExt(vector<fct<fctRangeType> > var) { return(".gpf"); }
template<typename fctRangeType>
string getFileExt(fct<fctRangeType> var) { return(".gpf"); }
template<typename matType>
string getFileExt(Mat<matType> var) { return(".gpm"); }
template<typename type>
string getFileExt(type var) { return(".gp"); }

template<typename fctRangeType>
vector<fct<fctRangeType> > outputModifier(vector<fct<fctRangeType> > var) { return(var); }
template<typename fctRangeType>
fct<fctRangeType> outputModifier(fct<fctRangeType> var) { return(var); }
template<typename matType>
Mat<matType> outputModifier(Mat<matType> var) { return(var); }
Mat<mt> outputModifier(Mat<complex<mt> > var) { return(real(var)); }
template<typename type>
type outputModifier(type var) { return(var); }


class inhLatticeSystem {

    #define outputWithHeader(x,header) file.open((dirName+string("/")+string(#x)+getFileExt(x)).c_str()); \
                                       file << header << endl << outputModifier(x) << endl; \
                                       file.close();
    #define outputWO(x) outputWithHeader(x,"");

//    #define outputFct(x,labels) file.open((string("output/")+string(#x)+string(".gp")).c_str()); \
//            file << "#\t" << labels << endl << x << endl; file.close()
//    #define outputFctWithLegend(x,labels,legend) file.open((string("output/")+string(#x)+string(".gp")).c_str()); \
//            file << "#\t" << labels << endl << "#\t" << legend << endl << x << endl; file.close(); 
//    #define outputWithExt(x,file_ext) file.open((string("output/")+string(#x)+string(file_ext)).c_str()); file << x << endl; file.close()
//    #define output(x) outputWithExt(x,"")

//    #define EXTRACHECKS

    public:
    ofstream file;
    ifstream ifile;

    fct<complex<mt> > potential;
    mt sampleWidth;
    mt odeSampleWidth;
    mt wannier0FunctionCutoff;
    mt HinW0cutoff;

    unsigned kSteps;     // also the number of Bloch states used for calculation of W0func
    unsigned nrSites;    // number of sites in the inh system

    blochFindMultipleRoots fR;

    unsigned minBands;
    clock_t tMark;

    Col<mt> H0eigenvalues;
    Col<mt> H0kvalues;
    vector<fct<complex<mt> > > H0eigenfunctions;
    vector<fct<complex<mt> > > centeredWannier0Functions;
    vector<fct<complex<mt> > > wannier0Functions;
    vector<fct<complex<mt> > > wannier0FunctionsPeriod;
    Col<mt> wannier0Eigenvalues;

    Mat<cx_double> HinW0partV;
    Mat<cx_double> HinW0partH0;
    Mat<cx_double> HinW0;

    Col<mt> HeigenvaluesESorted;
    Mat<cx_double> W0HOverlapESorted;
    Mat<complex<mt> > XinW0;
    Col<mt> Hlocalization;
    Col<mt> bandToEindex;
    Col<mt> Heigenvalues;
    Mat<cx_double> W0HOverlap;
    vector<fct<complex<mt> > > Heigenfunctions;
    Mat<cx_double> wannierInH;
    Col<mt> wannierEigenvalues;
    Mat<cx_double> HWannierOverlap;
    Mat<cx_double> HWannierOverlapRephased;
    vector<fct<complex<mt> > > wannierFunctions;
    vector<fct<complex<mt> > > wannierFunctionsRephased;
    vector<fct<complex<mt> > > HeigenfunctionsRephased;
    Mat<cx_double> hopping;

    void outputPropagatedWavefunctionEndValueOverE() {
        propStart=0.;
        propEnd=pp.s_bound;
        pp.k=0*PI/pp.a; odeStepWidth=DBL_MAX;
            ::outputPropagatedWavefunctionEndValueOverE(fR.intervalStart, fR.intervalEnd, (fR.intervalEnd-fR.intervalStart)/1000);
    }

    void outputPotential0() {
            fct<complex<mt> > potential0;
            for(unsigned c=0;c<=1000;c++) {
                mt x=pp.a*(-(mt)kSteps/2+(mt)c*kSteps/1000);
                potential0.x.push_back(x);
                potential0.f.push_back(Vlat(x,pp));
            }
            outputWithHeader(potential0,"#\t\tx\tV(x)/E_r\t\twith lines");
    }

    void outputPotential() {
            //fct<complex<mt> > potential;
            for(unsigned c=0;c<=1000;c++) {
                mt x=pp.a*(-(mt)nrSites/2+(mt)c*nrSites/1000);
                potential.x.push_back(x);
                potential.f.push_back(Vlat(x,pp)+Vinh(x));
            }
            outputWithHeader(potential,"#\t\tx\tV(x)/E_r\t\twith lines");
    }

    void calculateH0eigenvalues() {
        fout << "calculateH0eigenvalues() executing ... "; fout.flush();
        tMark=clock();
        propStart=0.;
        propEnd=pp.s_bound;
        minBands=UINT_MAX;
        vector<fct<mt> > H0dispersion;
        odeStepWidth=DBL_MAX;
        mt kStepSize=2.*PI/kSteps/pp.a;
        pp.k=-PI/pp.a;
        fR.rf=invPropagateFunction;
        fR.mf=invPropagateFunction;
        fR.searchDownUp();
        set<mt> kPIRoots=fR.allRoots;
        fR.allRoots.clear();
        fR.allMins.clear();
        pp.k=0;
        fR.rf=propagateFunction;
        fR.mf=propagateFunction;
        fR.searchDownUp();
        set<mt> k0Roots=fR.allRoots;
        fR.allRoots.clear();
        fR.allMins.clear();
        minBands=min(kPIRoots.size(),k0Roots.size());
        set<mt>::iterator iPI=kPIRoots.begin();
        set<mt>::iterator i0=k0Roots.begin();
        for(unsigned c=0;c<minBands;c++) {
            H0dispersion.push_back(fct<mt>(-PI/pp.a, PI/pp.a-kStepSize, kSteps));
            H0dispersion[c].x[0]=-PI/pp.a;
            H0dispersion[c].f[0]=*iPI;
            H0dispersion[c].x[kSteps/2]=0;
            H0dispersion[c].f[kSteps/2]=*i0;
            iPI++;
            i0++;
        }
        for(unsigned c=1;c<kSteps/2;c++) {
            pp.k=-PI/pp.a+c*kStepSize;
            set<mt>::iterator iPI=kPIRoots.begin();
            set<mt>::iterator i0=k0Roots.begin();
            for(unsigned d=0;d<minBands;d++) {
                mt minE=min(*iPI,*i0);
                mt maxE=max(*iPI,*i0);
                mt intervalStretch=0.0001;
                while( sgn(propagateFunction(minE))==sgn(propagateFunction(maxE)) ) {
                    minE-=intervalStretch;
                    maxE+=intervalStretch;
                    intervalStretch*=2;
                }
                mt foundRoot=Ridder_find_val(propagateFunction,0.,minE,maxE,fR.rootValueAccuracy);
                calculatedValues.clear();
                H0dispersion[d].x[c]       =pp.k;
                H0dispersion[d].f[c]       =foundRoot;
                H0dispersion[d].x[kSteps-c]=-pp.k;
                H0dispersion[d].f[kSteps-c]=foundRoot;
                iPI++;
                i0++;
            }
        }

        //H0dispersion.resize(minBands);
//        outputFct(H0dispersion,/*"dispersion relation*/"\tk/a\tE/E_r");
        outputWithHeader(H0dispersion,"#\t\tk/a\tE/E_r\t");
        for(unsigned c=0;c<H0dispersion.size();c++)
            H0dispersion[c].setSampleWidth(kStepSize);

        H0eigenvalues=Col<mt>(minBands*statesPerBand0);
        H0kvalues=Col<mt>(minBands*statesPerBand0);
        for(unsigned c=0;c<H0dispersion.size();c++)
            for(unsigned d=0;d<H0dispersion[c].size();d++) {
                H0eigenvalues[c*statesPerBand0+d]=H0dispersion[c].f[d];
                H0kvalues[c*statesPerBand0+d]=H0dispersion[c].x[d];
            }
//        output(H0eigenvalues);
        outputWO(H0eigenvalues);
//        output(H0kvalues);
        outputWO(H0kvalues);

        //Col<mt> H0eigenvaluesSorted=sort(H0eigenvalues);
        //output(H0eigenvaluesSorted);
        //outputWO(H0eigenvaluesSorted);
        fout << (clock()-(float)tMark)/CLOCKS_PER_SEC << 's' << endl;
    }

    void calculateH0eigenfunctions() {
        if(H0eigenvalues.n_elem==0)
            calculateH0eigenvalues();

        fout << "calculateH0eigenfunctions() executing ... "; fout.flush();
        tMark=clock();
        H0eigenfunctions=vector<fct<complex<mt> > >(minBands*statesPerBand0);
        for(unsigned c=0;c<H0eigenvalues.size();c++) {
            odeStepWidth=DBL_MAX;
            propStart=0.;
            propEnd=pp.s_bound;
            pp.k=H0kvalues[c];
            propagateFunction(H0eigenvalues[c]);
            calculatedValues.clear();
            mt R1d=D1->y[D1->P.t_length-1][1];
            mt R2d=D2->y[D2->P.t_length-1][1];
            mt I1=D1->y[D1->P.t_length-1][2];
            mt I2=D2->y[D2->P.t_length-1][2];

            //fout ex(H0eigenvalues[c]) ex(I1) ex(I2) ex(R1d) ex(R2d);
            mt C1,C2,zero=1E-4;
            if(fabs(I1)>zero) {
                if(fabs(I2)>zero) {
                    C2=sqrt(1./(1.+( I2*I2/I1/I1 )));
                    C1=-I2/I1*C2;
                } else {
                    C1=0;
                    C2=1;
                }
            } else {
                if(fabs(I2)>zero) {
                    C1=1;
                    C2=0;
                } else {
                    if(fabs(R1d)>zero) {
                        if(fabs(R2d)>zero) {
                            C2=sqrt(1./(1.+( R2d*R2d/R1d/R1d )));
                            C1=-R2d/R1d*C2;
                        } else {
                            C1=0;
                            C2=1;
                        }
                    } else {
                        if(fabs(R2d)>zero) {
                            C1=1;
                            C2=0;
                        } else {
                            //fout << "degenerate";
                            /*fout << c << endl;
                            fout ex(I1) ex(I2) ex(R1d) ex(R2d) << endl;
                            fout << "Error: all propagated values close to 0!!! break on line: " << __LINE__ << endl;*/

                            if(fabs(H0kvalues[c])>PI/pp.a/2) {
                                if((c/statesPerBand0)%2==0) {
                                    C1=0; C2=1;
                                } else {
                                    C1=1; C2=0;
                                }
                            } else {
                                if((c/statesPerBand0)%2==1) {
                                    C1=0; C2=1;
                                } else {
                                    C1=1; C2=0;
                                }
                            }

                            //fout << H0eigenvalues[c] << ' '; fout.flush();
                            //fout << H0eigenvalues[c-statesPerBand0] << ' '; fout.flush();

                            bool C2set=true;
                            // C1=0; C2=1;
                            if(((signed)c-statesPerBand0)>=0) {
                                if(fabs(H0eigenvalues[c-statesPerBand0]-H0eigenvalues[c])<1E-1) {
                                    C2set=false;
                                    // C1=0; C2=1;
                                }
                            }

                            if(C2set && (C2!=1)) {
                                fout << "failure assigning degenerate eigenstates " << __LINE__ << endl;
                                abort();
                            }
                        }
                    }
                }
            }
    //        fout << endl;


    //      extractEquidistantComplexFctFromDGLsolver() {
                fct<complex<mt> > variableWidthSamples1;
                fct<complex<mt> > variableWidthSamples2;

                odeStepWidth=odeSampleWidth;
                propStart=0.-10*odeStepWidth;
                propEnd=pp.s_bound+10*odeStepWidth;
                propagateFunction(H0eigenvalues[c]);
                calculatedValues.clear();

    //          copySamplesFromDGLsolverToFctClass() {
                    for(unsigned e=0;e<D1->P.t_length;e++) {
                        variableWidthSamples1.x.push_back(D1->t[e]);
                        variableWidthSamples1.f.push_back(complex<double>(D1->y[e][0],D1->y[e][2]));
                    }
    //          }

    //          copySamplesFromDGLsolverToFctClass() {
                    for(unsigned e=0;e<D2->P.t_length;e++) {
                        variableWidthSamples2.x.push_back(D2->t[e]);
                        variableWidthSamples2.f.push_back(complex<double>(D2->y[e][0],D2->y[e][2]));
                    }
    //          }

                fct<complex<double> > tmpF1(0.,+pp.s_bound,sampleWidth);
                tmpF1.interp_spline_from(variableWidthSamples1);

                fct<complex<double> > tmpF2(0.,+pp.s_bound,sampleWidth);
                tmpF2.interp_spline_from(variableWidthSamples2);
    //      }

            fct<complex<double> > U(tmpF1*C1+tmpF2*C2);


            fct<complex<mt> > Ee;
            Ee.setSampleWidth(sampleWidth);

    //      calculatePeriodicPhiFromU() {
                for(unsigned e=0;e<kSteps;e++) {
                    for(unsigned f=0;f<=U.size()-2;f++) {
                        mt x= ((mt)e-(mt)kSteps/2)    *pp.a +U.x[f];
                        Ee.x.push_back(x);
                        Ee.f.push_back(U.f[f]*(cos(H0kvalues[c]*x)+complex<double>(0,1)*sin(H0kvalues[c]*x)));
                    }
                    for(unsigned f=0;f<=U.size()-2;f++) {
                        mt x= ((mt)e-(mt)kSteps/2+1.)*pp.a -U.x[U.size()-f-1];
                        Ee.x.push_back(x);
                        Ee.f.push_back(conj(U.f[U.size()-f-1])*(cos(H0kvalues[c]*x)+complex<double>(0,1.)*sin(H0kvalues[c]*x)));
                    }
                }
    //      }
            Ee.x.push_back(-Ee.x[0]);
            Ee.f.push_back(Ee.f[0]);

            Ee.normalize_with_L2_norm();
#ifdef PFSI
            Ee.prepare_for_spline_interpolation();
#endif
            H0eigenfunctions[c]=Ee;
        }

//        outputWithHeader(H0eigenfunctions,"#\t\tx\t{/Symbol=\\302}({/Symbol=\\152}(x))\t\twith lines");

        stringstream kvalsHeader;
        kvalsHeader << "#\t\tx\t{/Symbol=\\302}({/Symbol=\\152}@_k^{(0)}(x))\t\twith lines" << endl;
        kvalsHeader << "#\t";
        for(unsigned c=0;c<H0eigenfunctions.size();c++)
            kvalsHeader << "k=" << setprecision(3) << H0kvalues[c] << '\t';
        outputWithHeader(H0eigenfunctions,kvalsHeader.str());
        fout << (clock()-(float)tMark)/CLOCKS_PER_SEC << 's' << endl;
    }

    void calculateOverlapWithH0eigenfunctionsFromFile() {
        if(H0eigenfunctions.size()==0)
            calculateH0eigenfunctions();

//      getWannierStatesFromFile() {
            unsigned EDsamples=800;
            unsigned EDbands=401;
            vector<fct<complex<mt> > > EDeigenfunctions(EDbands);
            for(unsigned c=0;c<EDbands;c++)
                EDeigenfunctions[c]=fct<complex<mt> >(-1.,1.-(2./800),(unsigned)EDsamples);// H0eigenfunctions[0]*0.;
            ifile.open("input/psi_400grid_400msize_2sites");
            for(unsigned c=0;c<EDsamples;c++) {
                mt x=-1.+(double)((c+EDsamples/4)%EDsamples)/400;
                for(unsigned d=0;d<EDbands/*minBands*/;d++) { // matrix_size=30 -> d<31 possible
                    mt ReF; ifile>>ReF;
                    mt ImF; ifile>>ImF;
                    complex<mt> f(ReF,ImF);
                    EDeigenfunctions[d].x[(c+EDsamples/4)%EDsamples]=x;
                    EDeigenfunctions[d].f[(c+EDsamples/4)%EDsamples]=f;
                }
            }
            ifile.close();
//      }

//        outputFct(EDeigenfunctions,/*"Wannier functions of lattice-only system not rephased*/"\tx\t{/Symbol=\\152}@_l^{(0)}(x)\t\twith lines");
        outputWithHeader(EDeigenfunctions,"#\t\tx\t{/Symbol=\\152}@_l^{(0)}(x)\t\twith lines");

        for_each(H0eigenfunctions,normalizeFct);
#ifdef PFSI
        for(unsigned c=0;c<H0eigenfunctions.size();c++)
            H0eigenfunctions[c].prepare_for_spline_interpolation();
#endif
        for_each(EDeigenfunctions,normalizeFct);
#ifdef PFSI
        for(unsigned c=0;c<EDeigenfunctions.size();c++)
            EDeigenfunctions[c].prepare_for_spline_interpolation();
#endif

        for(unsigned c=0;c<minBands;c++) {
            fout << c << ' ' << fabs(romberg_integral(H0eigenfunctions[2*c+1].conj(),EDeigenfunctions[c],
                                                 H0eigenfunctions[2*c+1].x[0],H0eigenfunctions[2*c+1].x[H0eigenfunctions[2*c+1].x.size()-1],1E-8)) << endl;
        }
    }

    void calculateSchroedinger0Error() {
        if(H0eigenfunctions.size()==0)
            calculateH0eigenfunctions();

        fout << "calculateSchroedinger0Error() executing ... "; fout.flush();
        tMark=clock();
        vector<fct<complex<mt> > > Schroedinger0(H0eigenfunctions.size());
        Col<mt> Schroedinger0Error(H0eigenfunctions.size());
        Col<mt> Schroedinger0BandError(minBands);
        for(unsigned c=0;c<H0eigenfunctions.size();c++) {
    //            fout ex(H0eigenvalues[c]) ex(H0kvalues[c]) << endl;
            Schroedinger0[c]=H0eigenfunctions[c]*0;
            free(Schroedinger0[c].second_deriv_for_fcts_spline);
            Schroedinger0[c].second_deriv_for_fcts_spline=NULL;
            pp.k=H0kvalues[c];
            pp.E=H0eigenvalues[c];
            fct<complex<mt> > firstDerivative=H0eigenfunctions[c].deriv();
            fct<complex<mt> > secondDerivative=firstDerivative.deriv();
    /*            for(unsigned d=0;d<10;d++)
                fout << Vlat(H0eigenfunctions[c].x[d],pp)-pp.hbar*pp.hbar/2./pp.m*(secondDerivative[d]/H0eigenfunctions[c][d]) << ' ' << H0eigenfunctions[c].f[d] << ' ' << secondDerivative[d] << endl;*/
            for(unsigned d=4;d<H0eigenfunctions[c].size()-4;d++)
                Schroedinger0[c].f[d]=(Vlat(H0eigenfunctions[c].x[d],pp)-pp.E)*H0eigenfunctions[c][d]-pp.hbar*pp.hbar/2./pp.m*secondDerivative[d];
            Schroedinger0Error[c]=Schroedinger0[c].L2_norm();
        }
        for(unsigned c=0;c<minBands;c++) {
            mt sum=0;
            for(unsigned d=0;d<statesPerBand0;d++)
                sum+=Schroedinger0Error[c*statesPerBand0+d];
            Schroedinger0BandError[c]=sum/statesPerBand0;
        }
        fout << (clock()-(float)tMark)/CLOCKS_PER_SEC << 's' << endl;
//        output(Schroedinger0Error);
        outputWO(Schroedinger0Error);
//        output(Schroedinger0BandError);
        outputWO(Schroedinger0BandError);
        for(unsigned c=0;c<minBands;c++)
            fout << Schroedinger0BandError[c] << ' ';
        fout << endl;
//        outputFct(Schroedinger0,"\tx/a\tdE*phi\t\twith lines");
        outputWithHeader(Schroedinger0,"#\t\tx/a\tdE*phi\t\twith lines");
    }

    void H0H0Overlap() {
        if(H0eigenfunctions.size()==0)
            calculateH0eigenfunctions();

        Mat<cx_double> H0H0Overlap=calculateOverlapMatrix(H0eigenfunctions);
//        output(H0H0Overlap);
        outputWO(H0H0Overlap);
    }

    void hopping0OverDl() {
        if(H0eigenvalues.n_elem==0)
             calculateH0eigenvalues();

        vector<fct<complex<mt> > > hopping0OverDl(minBands);
        for(unsigned c=0;c<minBands;c++) {
            hopping0OverDl[c]=fct<complex<mt> >(1,statesPerBand0/2,statesPerBand0/2)*0;
            for(unsigned d=0;d<statesPerBand0/2;d++) {
                complex<mt> sum=0;
                for(unsigned e=0;e<statesPerBand0;e++)
                    sum+=H0eigenvalues[statesPerBand0*c+e]*(cos(H0kvalues[statesPerBand0*c+e]*d)+complex<mt>(0.,1.)*sin(H0kvalues[statesPerBand0*c+e]*d));
                hopping0OverDl[c][d]=fabs(sum);
            }
        }
        stringstream bandnumberHeader;
        bandnumberHeader << "#\t\tDelta l\tJ" << endl;
        bandnumberHeader << "#\t";
        for(unsigned c=0;c<minBands;c++)
            bandnumberHeader << "band " << c+1 << '\t';
//        outputFctWithLegend(hopping0OverDl,"\tDelta l\tJ",bandnumberHeader.str());
        outputWithHeader(hopping0OverDl,bandnumberHeader);
    }

    void calculateCenteredWannier0Functions() {
        if(H0eigenfunctions.size()==0)
            calculateH0eigenfunctions();

        fout << "calculateCenteredWannier0Functions() executing ... "; fout.flush();
        tMark=clock();
        centeredWannier0Functions=vector<fct<complex<mt> > >(minBands);
        for(unsigned c=0;c<minBands;c++) {
            fct<complex<mt> > ncwf=H0eigenfunctions[0]*0.;
    //      makePhaseRealAround.5() {
                if(c%2==0)
                    for(unsigned d=0;d<statesPerBand0;d++) {
                        complex<mt> phaseFactor=fabs(H0eigenfunctions[d+c*statesPerBand0](.5))/H0eigenfunctions[d+c*statesPerBand0](.5);     // possibly optimizable
                        ncwf+=H0eigenfunctions[d+c*statesPerBand0]*phaseFactor;
                    }
                else
                    for(unsigned d=0;d<statesPerBand0;d++) {
                        complex<mt> deriv=H0eigenfunctions[d+c*statesPerBand0].deriv()(.5);                                      // possibly optimizable
    //                    complex<mt> deriv=(H0eigenfunctions[d+c*statesPerBand0](.5001)-H0eigenfunctions[d+c*statesPerBand0](.4999))/0.0002;                                      // possibly optimizable
    /*                    fct<complex<mt> > derivF=H0eigenfunctions[d+c*statesPerBand0].deriv();
#ifdef PFSI     
                        derivF.prepare_for_spline_interpolation();
#endif // PFSI     
                        complex<mt> deriv=derivF(.5);                                      // possibly optimizable*/
                        complex<mt> phaseFactor=fabs(deriv)/deriv;
                        ncwf+=H0eigenfunctions[d+c*statesPerBand0]*phaseFactor;
                    }
    //      }

            for(unsigned d=0;d<ncwf.size();d++) {
                ncwf.x[d]-=0.5;
        /*        if(ncwf.x[d]==0)
                    Wmiddle=d;*/
            }
            unsigned fI; // first Inside
            for(fI=1;ncwf.x[fI]<-(mt)nrSites/2*pp.a-1E-6;fI++) {
                ncwf.x.push_back(ncwf.x[fI]+(mt)nrSites);
                ncwf.f.push_back(ncwf.f[fI]);
            }
            ncwf.x.erase(ncwf.x.begin(),ncwf.x.begin()+fI);
            ncwf.f.erase(ncwf.f.begin(),ncwf.f.begin()+fI);

            unsigned start, end;
            for(start=0;fabs(ncwf.f[start])<wannier0FunctionCutoff;start++)
                if(start==ncwf.size()-1) {
                    fout << "error in " << __FILE__ << ':' << __LINE__ << endl;
                    abort();
                }
            for(end=ncwf.size()-1;fabs(ncwf.f[end])<wannier0FunctionCutoff;end--)
                ;

/*            while(fabs(ncwf.x[start])+1E-6<fabs(ncwf.x[end])) {      probably erase this
                fout << "cutting in " << c << endl;
                end--;
            }
            while(fabs(ncwf.x[end])+1E-6<fabs(ncwf.x[start])) {
                fout << "cutting in " << c << endl;
                start++;
            }*/

            if(end<ncwf.size()-1) {
                unsigned tmps=ncwf.size();
                ncwf.x.erase(ncwf.x.begin()+end+1,ncwf.x.begin()+tmps);
                ncwf.f.erase(ncwf.f.begin()+end+1,ncwf.f.begin()+tmps);
            }
            if(start>0) {
                ncwf.x.erase(ncwf.x.begin()+0,ncwf.x.begin()+start);
                ncwf.f.erase(ncwf.f.begin()+0,ncwf.f.begin()+start);
            }

            centeredWannier0Functions[c]=ncwf;

        }
        H0eigenfunctions.clear();

        for(unsigned c=0;c<centeredWannier0Functions.size();c++) {
            centeredWannier0Functions[c].normalize_with_L2_norm();
#ifdef PFSI
            centeredWannier0Functions[c].prepare_for_spline_interpolation();
#endif //PFSI
        }

//        outputFct(centeredWannier0Functions,/*"centered Wannier functions of lattice-only*/"\tx\t{/Symbol=\\152}@_l^{(0)}(x)\t\twith lines");
        outputWithHeader(centeredWannier0Functions,"#\t\tx\t{/Symbol=\\152}@_l^{(0)}(x)\t\twith lines");
//        outputWithHeader(centeredWannier0Functions[0],"#\t\tx\t{/Symbol=\\152}@_l^{(0)}(x)\t\twith lines");
//        outputWithHeader(centeredWannier0Functions[1],"#\t\tx\t{/Symbol=\\152}@_l^{(0)}(x)\t\twith lines");
//        outputWithHeader(centeredWannier0Functions[2],"#\t\tx\t{/Symbol=\\152}@_l^{(0)}(x)\t\twith lines");
        fout << (clock()-(float)tMark)/CLOCKS_PER_SEC << 's' << endl;
    }

    void calculateWannier0Functions() {
        if(centeredWannier0Functions.size()==0)
            calculateCenteredWannier0Functions();

        fout << "calculateWannier0Functions() executing ... "; fout.flush();
        tMark=clock();
        wannier0Functions=vector<fct<complex<mt> > >(minBands*nrSites);
        wannier0FunctionsPeriod=vector<fct<complex<mt> > >(minBands*nrSites);
        for(unsigned c=0;c<minBands;c++) {
            for(unsigned d=0;d<nrSites;d++) {
                mt xLocOfCenter=(-(mt)nrSites/2+.5+d)*pp.a;
                wannier0Functions[c*nrSites+d].setSampleWidth(sampleWidth);
                for(unsigned e=0;e<centeredWannier0Functions[c].size();e++) {
                    mt xLoc=xLocOfCenter+centeredWannier0Functions[c].x[e];
                    if((-(mt)nrSites/2*pp.a-1E-6 < xLoc) && (xLoc < +(mt)nrSites/2*pp.a+1E-6)) {
                        //if((-(mt)nrSites/2*pp.a-1E-6 < centeredWannier0Functions[c].x[e]) && (centeredWannier0Functions[c].x[e] < +(mt)nrSites/2*pp.a+1E-6)) { // was needed for nrSites>kValues, now nrSites has to be = kValues
                            wannier0Functions[c*nrSites+d].x.push_back(xLoc);
                            wannier0Functions[c*nrSites+d].f.push_back(centeredWannier0Functions[c].f[e]);
                        //}
                    }
                }
#ifdef PFSI
                wannier0Functions[c*nrSites+d].prepare_for_spline_interpolation();
#endif //PFSI

                mt xLocOfCenterPeriod;
                if(xLocOfCenter<0)
                    xLocOfCenterPeriod=(+(mt)nrSites/2+.5+d)*pp.a;
                else
                    xLocOfCenterPeriod=(-(mt)3*nrSites/2+.5+d)*pp.a;
                wannier0FunctionsPeriod[c*nrSites+d].setSampleWidth(sampleWidth);
                for(unsigned e=0;e<centeredWannier0Functions[c].size();e++) {
                    mt xLoc=xLocOfCenterPeriod+centeredWannier0Functions[c].x[e];
                    if((-(mt)nrSites/2*pp.a-1E-6 < xLoc) && (xLoc < +(mt)nrSites/2*pp.a+1E-6)) {
                        //if((-(mt)nrSites/2*pp.a-1E-6 < centeredWannier0Functions[c].x[e]) && (centeredWannier0Functions[c].x[e] < +(mt)nrSites/2*pp.a+1E-6)) {
                            wannier0FunctionsPeriod[c*nrSites+d].x.push_back(xLoc);
                            wannier0FunctionsPeriod[c*nrSites+d].f.push_back(centeredWannier0Functions[c].f[e]);
                        //}
                    }
                }
#ifdef PFSI
                wannier0FunctionsPeriod[c*nrSites+d].prepare_for_spline_interpolation();
#endif //PFSI
            }
        }

//        outputFct(wannier0Functions,/*"centered Wannier functions of lattice-only*/"\tx\t{/Symbol=\\152}@_l^{(0)}(x)\t\twith lines");
        outputWithHeader(wannier0Functions,"#\t\tx\t{/Symbol=\\152}@_l^{(0)}(x)\t\twith lines");
//        outputWithHeader(wannier0Functions[0],"#\t\tx\t{/Symbol=\\152}@_l^{(0)}(x)\t\twith lines");

//        outputFct(wannier0FunctionsPeriod,/*"centered Wannier functions of lattice-only*/"\tx\t{/Symbol=\\152}@_l^{(0)}(x)\t\twith lines");
        outputWithHeader(wannier0FunctionsPeriod,"#\t\tx\t{/Symbol=\\152}@_l^{(0)}(x)\t\twith lines");
//        outputWithHeader(wannier0FunctionsPeriod[0],"#\t\tx\t{/Symbol=\\152}@_l^{(0)}(x)\t\twith lines");
        fout << (clock()-(float)tMark)/CLOCKS_PER_SEC << 's' << endl;
    }

    void calculateWannier0Overlap() {
        if(wannier0Functions.size()==0)
            calculateWannier0Functions();

        fout << "calculateWannier0Overlap() executing ... "; fout.flush();
        tMark=clock();
        Mat<cx_double> wannier0Overlap(wannier0Functions.size(),wannier0Functions.size());
        Mat<cx_double> wannier0withPeriodOverlap(wannier0FunctionsPeriod.size(),wannier0FunctionsPeriod.size());
        for(unsigned c=0;c<wannier0Functions.size();c++)
            for(unsigned d=0;d<wannier0Functions.size();d++)
                wannier0Overlap(c,d)=romberg_integral_f(wannier0Functions[c].conj(),wannier0Functions[d],constantOp,1E-6);
        for(unsigned c=0;c<wannier0FunctionsPeriod.size();c++)
            for(unsigned d=0;d<wannier0FunctionsPeriod.size();d++) {
                wannier0withPeriodOverlap(c,d) =romberg_integral_f(wannier0Functions[c].conj(),      wannier0Functions[d],      constantOp,1E-6);
                wannier0withPeriodOverlap(c,d)+=romberg_integral_f(wannier0FunctionsPeriod[c].conj(),wannier0FunctionsPeriod[d],constantOp,1E-6);
                wannier0withPeriodOverlap(c,d)+=romberg_integral_f(wannier0FunctionsPeriod[c].conj(),wannier0Functions[d],      constantOp,1E-6);
                wannier0withPeriodOverlap(c,d)+=romberg_integral_f(wannier0Functions[c].conj(),      wannier0FunctionsPeriod[d],constantOp,1E-6);
            }
        file.open((dirName+"/wannier0Overlap.gpf").c_str()); 
        for(unsigned c=0;c<wannier0Overlap.n_rows;c++) {
            for(unsigned d=0;d<wannier0Overlap.n_rows;d++)
                if(abs(wannier0Overlap(c,d))<.01)
                    file << "0.      " << ' ';
                else
                    file << left << setw(8) << setprecision(1) << abs(wannier0Overlap(c,d)) << ' ';
            file << endl;
        }
        file.close();
        file.open((dirName+"/wannier0withPeriodOverlap.gpf").c_str());
        for(unsigned c=0;c<wannier0withPeriodOverlap.n_rows;c++) {
            for(unsigned d=0;d<wannier0withPeriodOverlap.n_rows;d++)
                if(abs(wannier0withPeriodOverlap(c,d))<.01)
                    file << "0.      " << ' ';
                else
                    file << left << setw(8) << setprecision(1) << abs(wannier0withPeriodOverlap(c,d)) << ' ';
            file << endl;
        }
        file.close();
        file.precision(4);
        file.width(1);
        fout << (clock()-(float)tMark)/CLOCKS_PER_SEC << 's' << endl;
    }

    void calculateGaussianInWannier0Basis() {
        if(wannier0Functions.size()==0)
            calculateWannier0Functions();

        fout << "calculateWannier0Overlap() executing ... "; fout.flush();
        tMark=clock();
        fct<complex<mt> > testFct=fct<complex<mt> >(-((double)nrSites/2)*pp.a,+((double)nrSites/2)*pp.a,sampleWidth)*0;
        testFct.assign_fct(testF);
        testFct.normalize_with_L2_norm();
        //outputFct(testFct,/*"test function*/"\tx\tf(x)\t\twith lines");
        outputWithHeader(testFct,"#\t\tx\tf(x)\t\twith lines");
        Col<complex<mt> > W0Foverlap(wannier0Functions.size());
        for(unsigned c=0;c<wannier0Functions.size();c++)
            W0Foverlap[c]=romberg_integral_f(wannier0Functions[c],testFct,constantOp,1E-6);
        fct<complex<mt> > testFctRecalc=fct<complex<mt> >(-((double)nrSites/2)*pp.a,+((double)nrSites/2)*pp.a,sampleWidth)*0;
        for(unsigned c=0;c<wannier0Functions.size();c++)
            testFctRecalc+=wannier0Functions[c]*W0Foverlap[c];
//        testFctRecalc.normalize_with_L2_norm();
        //outputFct(testFctRecalc,/*"test function*/"\tx\tf(x)\t\twith lines");
        outputWithHeader(testFctRecalc,"#\t\tx\tf(x)\t\twith lines");
        fout << "testFctOverlap = " << romberg_integral_f(testFct,testFctRecalc,constantOp,1E-6) << endl;

        fout << (clock()-(float)tMark)/CLOCKS_PER_SEC << 's' << endl;
    }

    void calculateOverlapWithWannierStatesFromFile() {
        if(centeredWannier0Functions.size()==0)
            calculateCenteredWannier0Functions();

        fout << "calculateOverlapWithWannierStatesFromFile() executing ... "; fout.flush();
        tMark=clock();
    //  getWannierStatesFromFile() {
            vector<fct<complex<mt> > > EDcenteredWannier0Functions(31);
            ifstream ifile;
            ifile.open("input/EDcenteredWannier0Functions");
            for(unsigned c=0;c<256*40+1;c++) {
                mt x=-128+(double)c/40;
                for(unsigned d=0;d<31/*minBands*/;d++) { // matrix_size=30 -> d<31 possible
                    mt ReF; ifile>>ReF;
                    mt ImF; ifile>>ImF;
                    complex<mt> f(ReF,ImF);
                    EDcenteredWannier0Functions[d].x.push_back(x);
                    EDcenteredWannier0Functions[d].f.push_back(f);
                }
            }
            ifile.close();
    //  }

        //outputFct(EDcenteredWannier0Functions,/*"Wannier functions of lattice-only system not rephased*/"\tx\t{/Symbol=\\152}@_l^{(0)}(x)\t\twith lines");
        outputWithHeader(EDcenteredWannier0Functions,"#\t\tx\t{/Symbol=\\152}@_l^{(0)}(x)\t\twith lines");

        for_each(centeredWannier0Functions,normalizeFct);
#ifdef PFSI
        for(unsigned c=0;c<centeredWannier0Functions.size();c++)
            centeredWannier0Functions[c].prepare_for_spline_interpolation();
#endif //PFSI
        for_each(EDcenteredWannier0Functions,normalizeFct);
#ifdef PFSI
        for(unsigned c=0;c<EDcenteredWannier0Functions.size();c++)
            EDcenteredWannier0Functions[c].prepare_for_spline_interpolation();
#endif //PFSI

        for(unsigned c=0;c<minBands;c++) {
            fout << c << ' ' << fabs(romberg_integral(centeredWannier0Functions[c].conj(),EDcenteredWannier0Functions[c],
                                                     centeredWannier0Functions[c].x[0],centeredWannier0Functions[c].x[centeredWannier0Functions[c].size()-1],1E-6)) << endl;
        }
        for(unsigned c=0;c<30;c++)
            fout << real(romberg_integral_f(EDcenteredWannier0Functions[c].conj(),EDcenteredWannier0Functions[c],squareOp,1E-6)) << ' ';
        fout << endl;
        fout << (clock()-(float)tMark)/CLOCKS_PER_SEC << 's' << endl;
    }

    void hopping0viaInt() {
        if(wannier0Functions.size()==0)
            calculateWannier0Functions();

        Mat<cx_double> hopping0viaInt(wannier0Functions.size(),wannier0Functions.size());
        hopping0viaInt.zeros();
        for(unsigned c=0;c<wannier0Functions.size();c++) {
    //        Schroedinger[c]=Heigenfunctions[c]*0;
    //        pp.E=Heigenvalues[c];
            fct<complex<mt> > firstDerivative=wannier0Functions[c].deriv();
            fct<complex<mt> > secondDerivative=firstDerivative.deriv();
            for(unsigned d=0;d<wannier0Functions.size();d++) {
                fout ex(c) ex(d) << endl;
                if((c/statesPerBand)==(d/statesPerBand)) {
                    if(d>=c) {                                  // this optimization doesn't yield exactly the same result, probably due to numerical inaccuracies
                        hopping0viaInt(d,c)=romberg_integral_f(wannier0Functions[d].conj(),wannier0Functions[c],Vlat,1E-6);
                        hopping0viaInt(d,c)-=pp.hbar*pp.hbar/2/pp.m*
                                     romberg_integral_f(wannier0Functions[d].conj(),secondDerivative,constantOp,1E-6);
                    } else
                        hopping0viaInt(d,c)=conj(hopping0viaInt(c,d));
                }
            }
        }
        //output(hopping0viaInt);
        outputWO(hopping0viaInt)

        vector<fct<mt> > hopping0viaIntOverBandindex(statesPerBand/2);
        for(unsigned c=0;c<statesPerBand/2;c++) {
            hopping0viaIntOverBandindex[c]=fct<mt>(1,minBands,minBands)*0.;
            for(unsigned d=0;d<minBands;d++)
                hopping0viaIntOverBandindex[c].f[d]=fabs(hopping0viaInt(d*statesPerBand+statesPerBand/2-1-c,d*statesPerBand+statesPerBand/2));
        }
        //outputFct(hopping0viaIntOverBandindex,"\tbandindex\tJ");
        outputWithHeader(hopping0viaIntOverBandindex,"#\t\tbandindex\tJ");

        vector<fct<mt> > hopping0viaIntOverDl(minBands);
        for(unsigned c=0;c<minBands;c++) {
            hopping0viaIntOverDl[c]=fct<mt>(1,statesPerBand/2,statesPerBand/2)*0.;
            for(unsigned d=0;d<statesPerBand/2;d++)
                hopping0viaIntOverDl[c].f[d]=fabs(hopping0viaInt(c*statesPerBand+statesPerBand/2-1-d,c*statesPerBand+statesPerBand/2));
        }
        //outputFct(hopping0viaIntOverDl,"\tDelta l\tJ");
        outputWithHeader(hopping0viaIntOverDl,"#\t\tDelta l\tJ");
    }

    bool calculateHinW0(unsigned n, unsigned m, unsigned &count) {
        bool breakLoop=false;
        //fout << n << ' ' << m << endl;
        if(m>=n) {
            HinW0partV(n,m) =romberg_integral_f(wannier0Functions[n].conj()      , wannier0Functions[m]      , Vinh, 1E-6);
            HinW0partV(n,m)+=romberg_integral_f(wannier0FunctionsPeriod[n].conj(), wannier0FunctionsPeriod[m], Vinh, 1E-6);
            HinW0partV(n,m)+=romberg_integral_f(wannier0FunctionsPeriod[n].conj(), wannier0Functions[m]      , Vinh, 1E-6);
            HinW0partV(n,m)+=romberg_integral_f(wannier0Functions[n].conj()      , wannier0FunctionsPeriod[m], Vinh, 1E-6);
//                        HinW0partV(n,m)=romberg_integral(W0conjTimesTrap, wannier0Functions[m], W0conjTimesTrap.x[0], W0conjTimesTrap.x[W0conjTimesTrap.x.size()-1], 1E-6);
            count++;

            if((n/statesPerBand)==(m/statesPerBand)) {
                complex<mt> sum=0;
                for(unsigned f=0;f<statesPerBand0;f++) {
//                        mt angle=H0kvalues[f]*(wannier0Eigenvalues[m]-wannier0Eigenvalues[n]); // this would be correct if the system was infinite
//                                                                                                        // then the eigenvalue difference would be integer
                    mt angle=H0kvalues[bandToIndex0(n/statesPerBand,f)]*                      // for the finite system the rounded quantities i.m. the indices
                                       pp.a*((signed)(m%statesPerBand)-(signed)(n%statesPerBand));
                    sum+=H0eigenvalues[bandToIndex0(n/statesPerBand,f)]*(cos(angle)+complex<mt>(0,1)*sin(angle));
                }
                sum/=(statesPerBand0*pp.a);
                HinW0partH0(n,m)=sum;
            }

            HinW0(n,m)=HinW0partV(n,m)+HinW0partH0(n,m);
            if(fabs(HinW0(n,m))<HinW0cutoff*fabs(HinW0(n,n)))
                breakLoop=true;
        }
        return(breakLoop);
    }

    void calculateHinW0() {
        if(wannier0Functions.size()==0)
            calculateWannier0Functions();
        if(H0eigenvalues.n_elem==0)
            calculateH0eigenvalues();

        fout << "calculateHinW0() executing ... "; fout.flush();
        tMark=clock();
        unsigned count=0;
        HinW0partH0=Mat<cx_double>(wannier0Functions.size(),wannier0Functions.size());
        HinW0partV=Mat<cx_double>(wannier0Functions.size(),wannier0Functions.size());
        HinW0=Mat<cx_double>(wannier0Functions.size(),wannier0Functions.size());
        HinW0partH0.zeros();
        HinW0partV.zeros();
        HinW0.zeros();
        for(unsigned c=0;c<wannier0Functions.size();c++) {
            for(unsigned d=0;d<minBands;d++) {
                for(unsigned e=0;e<statesPerBand/2;e++) {
                    unsigned f=(c%statesPerBand+statesPerBand-e-1)%statesPerBand;
                    f+=d*statesPerBand;
                    if(calculateHinW0(c,f,count))
                        break;
                }
                for(unsigned e=0;e<statesPerBand/2;e++) {
                    unsigned f=(c%statesPerBand+statesPerBand+e)%statesPerBand;
                    f+=d*statesPerBand;
                    if(calculateHinW0(c,f,count))
                        break;
                }
            }
        }
        for(unsigned c=0;c<wannier0Functions.size();c++) {
            for(unsigned d=0;d<c;d++) {
                HinW0partH0(c,d)=conj(HinW0partH0(d,c));
                HinW0partV(c,d) =conj(HinW0partV(d,c));
                HinW0(c,d)      =conj(HinW0(d,c));
            }
        }
//        output(HinW0);

        stringstream s; s<<"#\t\tl_0'+"<<statesPerBand<<"*{/Symbol=\\141}'\tl_0+"<<statesPerBand<<"*{/Symbol=\\141}\t<l_0,{/Symbol=\\141}|H_0/E_r|l_0',{/Symbol=\\141}'>\t\tset yrange [*:*] reverse;";
        outputWithHeader(HinW0partH0,s.str());
        stringstream s2; s2<<"#\t\tl_0'+"<<statesPerBand<<"*{/Symbol=\\141}'\tl_0+"<<statesPerBand<<"*{/Symbol=\\141}\t<l_0,{/Symbol=\\141}|V_{inh}/E_r|l_0',{/Symbol=\\141}'>\t\tset yrange [*:*] reverse;";
        outputWithHeader(HinW0partV,s2.str());
        stringstream s3; s3<<"#\t\tl_0'+"<<statesPerBand<<"*{/Symbol=\\141}'\tl_0+"<<statesPerBand<<"*{/Symbol=\\141}\t<l_0,{/Symbol=\\141}|H/E_r|l_0',{/Symbol=\\141}'>\t\tset yrange [*:*] reverse;";
        outputWithHeader(HinW0,s3.str());
        fout << (clock()-(float)tMark)/CLOCKS_PER_SEC << "s " ex(count) << endl;
    }

//    void calculateH0inW0fromW0inH0() {
//        #ifdef EXTRACHECKS // calculate H0inW0 from W0inH0  not valid for multiband
//            Mat<cx_double> altH0inW0(wannier0Functions.size(),wannier0Functions.size());
//            for(unsigned c=0;c<statesPerBand;c++)
//                for(unsigned d=0;d<statesPerBand;d++) {
//                    complex<mt> sum=0;
//                    for(unsigned e=0;e<statesPerBand0;e++)
//                        for(unsigned f=0;f<statesPerBand0;f++)
//                            sum+=(H0eigenvalues[e]*conj(H0Wannier0Overlap(e,c))*H0Wannier0Overlap(f,d)/pp.a/wannier0Eigenvalues(d))*wannier0InH0(e,f);
//                    altH0inW0(c,d)=sum;
//                }
//    
////            output(altH0inW0);
//            outputWO(altH0inW0);
//        #endif
//    }

    void calculateHeigenvaluesESorted() {
        if(HinW0.n_elem==0)
            calculateHinW0();

        fout << "calculateHeigenvaluesESorted() executing ... "; fout.flush();
        tMark=clock();
        HeigenvaluesESorted=Col<mt>(wannier0Functions.size());
        W0HOverlapESorted=Mat<cx_double>(wannier0Functions.size(),wannier0Functions.size());
        tMark=clock();
        eig_sym(HeigenvaluesESorted,W0HOverlapESorted,HinW0);

        //output(HeigenvaluesESorted);
        outputWO(HeigenvaluesESorted);
        //output(W0HOverlapESorted);
        stringstream s; s<<"#\t\tn\tl_0+"<<statesPerBand<<"*{/Symbol=\\141}\t<l_0,{/Symbol=\\141}|n>\t\tset yrange [*:*] reverse;";
        outputWithHeader(W0HOverlapESorted,s.str());
        fout << (clock()-(float)tMark)/CLOCKS_PER_SEC << 's' << endl;
    }

    void calculateBand0ProjectorsExpectationInHeigenstates() {
        if(W0HOverlapESorted.n_elem==0)
            calculateHeigenvaluesESorted();

        fout << "calculateBand0ProjectorsExpectationInHeigenstates() executing ... "; fout.flush();
        tMark=clock();
        Mat<mt> band0ProjectorsExpectationInHeigenstates(minBands,statesPerBand*minBands);
        for(unsigned c=0;c<minBands;c++) {
            for(unsigned d=0;d<statesPerBand*minBands;d++) {
                mt sum=0;
                for(unsigned e=0;e<statesPerBand;e++)
                    sum+=norm(W0HOverlapESorted(e+c*statesPerBand,d));
                band0ProjectorsExpectationInHeigenstates(c,d)=sum;
            }
        }
//        for(unsigned c=0;c<statesPerBand*minBands;c++) {
//            mt sum=0.;
//            for(unsigned d=statesPerBand/*0*/;d<2*statesPerBand;d++)
//                sum+=norm(W0HOverlapESorted(d,c));
//            band0ProjectorExpectationInHeigenstates(c)=sum;
//        }
//        output(band0ProjectorsExpectationInHeigenstates);
        outputWithHeader(band0ProjectorsExpectationInHeigenstates,"#\t\tn\tP_i\texpectation\t\tset yrange [*:*] reverse;");

//        vector<map<mt,unsigned> > 
        bandToEindex=Col<mt>(nrSites*minBands);
        vector<vector<mt> > bandEtoEindex(minBands);
        for(unsigned c=0;c<statesPerBand*minBands;c++) {
            unsigned maxOvIndex=0;
            mt maxOv=-DBL_MAX;
            for(unsigned d=0;d<minBands;d++) {
                if(band0ProjectorsExpectationInHeigenstates(d,c)>maxOv) {
                    maxOv=band0ProjectorsExpectationInHeigenstates(d,c);
                    maxOvIndex=d;
                }
            }
            bandEtoEindex[maxOvIndex].push_back(c);
        }
        for(unsigned c=0;c<minBands;c++) {
            if(bandEtoEindex[c].size()!=statesPerBand)
                fout << "error assigning states to band " << c << " in " << __FILE__ << ':' << __LINE__ << endl;
        }
        for(unsigned c=0;c<minBands;c++)
            for(unsigned d=0;d<statesPerBand;d++)
                bandToEindex[d+c*statesPerBand]=bandEtoEindex[c][d];

        //for(unsigned c=0;c<nrSites*minBands;c++)
        //    bandToEindex[c]=c;
//        output(bandToEindex);
        outputWO(bandToEindex);
        fout << (clock()-(float)tMark)/CLOCKS_PER_SEC << 's' << endl;
    }

    void calculateHeigenvaluesOverTrapStrength() {
        if(HinW0partH0.n_elem==0)
            calculateHinW0();
        if(W0HOverlapESorted.n_elem==0)
            calculateHeigenvaluesESorted();

        fout << "calculateHeigenvaluesOverTrapStrength() executing ... "; fout.flush();
        tMark=clock();
        Mat<cx_double> HinW0(wannier0Functions.size(),wannier0Functions.size());
        Col<mt> HeigenvaluesESorted(wannier0Functions.size());
        Mat<cx_double> W0HOverlapESorted(wannier0Functions.size(),wannier0Functions.size());
        vector<fct<complex<mt> > > HeigenvaluesOverTrapStrength1(wannier0Functions.size());
        vector<fct<complex<mt> > > HeigenvaluesOverTrapStrength1Deriv(wannier0Functions.size());
        for(double inhFac=0.;inhFac<pp.inhS;inhFac+=(pp.inhS-0)/150) {
            HinW0=HinW0partH0+HinW0partV/pp.inhS*inhFac;
            eig_sym(HeigenvaluesESorted,W0HOverlapESorted,HinW0);
            for(unsigned c=0;c<wannier0Functions.size();c++) {
                HeigenvaluesOverTrapStrength1[c].x.push_back(inhFac);
                HeigenvaluesOverTrapStrength1[c].f.push_back(HeigenvaluesESorted(c));
                HeigenvaluesOverTrapStrength1Deriv[c].x.push_back(inhFac);
                HeigenvaluesOverTrapStrength1Deriv[c].f.push_back(dot(W0HOverlapESorted.col(c).t()*HinW0partV/pp.inhS,W0HOverlapESorted.col(c)));
            }
        }
        //outputFct(HeigenvaluesOverTrapStrength,/*"energy eigenvalues as function of inh strength*/"\tinh strength/E_r\tE/E_r");
        outputWithHeader(HeigenvaluesOverTrapStrength1,"#\t\tinh strength/E_r\tE/E_r\t\twith lines");
        //outputFct(HeigenvaluesOverTrapStrengthDeriv,/*"energy eigenvalues as function of inh strength*/"\tinh strength deriv/E_r\tE/E_r");
        outputWithHeader(HeigenvaluesOverTrapStrength1Deriv,"#\t\tinh strength deriv/E_r\tE/E_r\t\twith lines");
        fout << (clock()-(float)tMark)/CLOCKS_PER_SEC << "s" << endl;
    }

    void calculateXinW0() {
        if(wannier0Functions.size()==0)
            calculateWannier0Functions();

        fout << "calculateXinW0() executing ... "; fout.flush();
        tMark=clock();
        XinW0=Mat<complex<mt> >(wannier0Functions.size(),wannier0Functions.size());
        XinW0.zeros();
        for(unsigned c=0;c<wannier0Functions.size();c++)
            for(unsigned d=0;d<wannier0Functions.size();d++)
                if(d>=c)
                    XinW0(c,d)=romberg_integral_f(wannier0Functions[c].conj(), wannier0Functions[d], unityOp, 1E-6);
                else
                    XinW0(c,d)=conj(XinW0(d,c));
        fout << (clock()-(float)tMark)/CLOCKS_PER_SEC << "s " << endl;
    }

    void calculateHlocalization() {
        if(XinW0.n_elem==0)
            calculateXinW0();
        if(W0HOverlapESorted.n_elem==0)
            calculateHeigenvaluesESorted();

        Hlocalization=Col<mt>(wannier0Functions.size());
        for(unsigned c=0;c<wannier0Functions.size();c++)
            Hlocalization[c]=real(dot(W0HOverlapESorted.col(c).t()*XinW0,W0HOverlapESorted.col(c)));
    }

    void assignBandOrderEqualHeigenvalueOrder() {
        bandToEindex=Col<mt>(nrSites*minBands);
        for(unsigned c=0;c<nrSites*minBands;c++)
            bandToEindex[c]=c;
    }

    void assignBandOrderUsingLocalization() {
        if(Hlocalization.n_elem==0)
            calculateHlocalization();
        if(HeigenvaluesESorted.n_elem==0)
            calculateHeigenvaluesESorted();

        bandToEindex=Col<mt>(nrSites*minBands);
        map<mt,unsigned> locSorter;
        for(unsigned c=0;c<wannier0Functions.size();c++)
            locSorter.insert(pair<mt,unsigned>(Hlocalization[c],c));
    
        map<mt,unsigned>::iterator it=locSorter.begin();
        for(unsigned c=0;c<nrSites;c++) {
            map<mt,unsigned> bandSorter;
            for(unsigned d=0;d<minBands;d++) {
                bandSorter.insert(pair<mt,unsigned>(HeigenvaluesESorted[it->second],it->second));
                it++;
            }
            unsigned band=0;
            for(map<mt,unsigned>::iterator bandIt=bandSorter.begin();bandIt!=bandSorter.end();bandIt++) {
                bandToEindex[c+band*nrSites]=bandIt->second;
                band++;
            }
        }
//        output(bandToEindex);
        outputWO(bandToEindex);
    //  checkLocalization() {
            for(unsigned c=0;c<nrSites*minBands;c++)
                if(round(Hlocalization[bandToEindex[c]]-.5+nrSites/2)!=c%nrSites)
                    fout << "state " << bandToEindex[c] << " is not well localized at " << (mt)(c%nrSites)-nrSites/2+.5 << endl;
    //  }

        fct<mt> EOverX;
        for(unsigned c=0;c<minBands;c++)
            for(unsigned d=0;d<nrSites;d++) {
                EOverX.x.push_back(Hlocalization[bandToEindex[c*nrSites+d]]);
                EOverX.f.push_back(HeigenvaluesESorted[ bandToEindex[c*nrSites+d]]);
            }
    
        file.open((dirName+"/EOverX.gpf").c_str());
        file << "#\t" << '\t' << "x/a\t" << "E/E_r" << endl;
        for(unsigned c=0;c<wannier0Functions.size();c++)
            file << EOverX.x[c] << '\t' << EOverX.f[c] << endl;
        file.close();
    }
    
    void calculateWannier0Eigenvalues() {
        wannier0Eigenvalues=Col<mt>(wannier0Functions.size());
        for(unsigned c=0;c<wannier0Functions.size();c++)
            wannier0Eigenvalues[c]=((mt)(c%nrSites)-nrSites/2+.5)*pp.a;
//        output(wannier0Eigenvalues);
        outputWO(wannier0Eigenvalues);
    }

struct FsOverTS {
    Col<mt> Heigenvalues;
    Col<mt> HeigenvaluesDeriv;
    Col<mt> PositionOverTrapStrength;
    Col<mt> BpPositionOverTrapStrength;
    FsOverTS(unsigned n) : Heigenvalues(n), HeigenvaluesDeriv(n), PositionOverTrapStrength(n), BpPositionOverTrapStrength(n) {
    }
};

FsOverTS calcFsOverTS(mt const inhStrength) {
    FsOverTS F(wannier0Functions.size());
    Mat<cx_double> HinW0(wannier0Functions.size(),wannier0Functions.size());
    Mat<cx_double> W0HOverlap(wannier0Functions.size(),wannier0Functions.size());
    HinW0=HinW0partH0+HinW0partV/pp.inhS*inhStrength;
    eig_sym(F.Heigenvalues,W0HOverlap,HinW0);

    for(unsigned c=0;c<wannier0Functions.size();c++) {
        F.HeigenvaluesDeriv(c)=fabs(dot(W0HOverlap.col(c).t()*HinW0partV/pp.inhS,W0HOverlap.col(c)));  // CHECK FABS
    }
    for(unsigned c=0;c<wannier0Functions.size();c++) {
        F.PositionOverTrapStrength(c)=real(dot(W0HOverlap.col(c).t()*XinW0,W0HOverlap.col(c)));
    }
    for(unsigned c=0;c<wannier0Functions.size();c++) {
        F.BpPositionOverTrapStrength(c)=dot(square(abs(W0HOverlap.col(c))),wannier0Eigenvalues);
    }

    return F;
}

void saveFsOverTS(mt const inhStrength,FsOverTS const &F) {
    for(unsigned c=0;c<wannier0Functions.size();c++) {
        HeigenvaluesOverTrapStrength[c].x.push_back(inhStrength);
        HeigenvaluesOverTrapStrength[c].f.push_back(F.Heigenvalues[c]);
    }
    for(unsigned c=0;c<wannier0Functions.size();c++) {
        HeigenvaluesOverTrapStrengthDeriv[c].x.push_back(inhStrength);
        HeigenvaluesOverTrapStrengthDeriv[c].f.push_back(F.HeigenvaluesDeriv[c]);
    }
    for(unsigned c=0;c<wannier0Functions.size();c++) {
        PositionOverTrapStrength[c].x.push_back(inhStrength);
        PositionOverTrapStrength[c].f.push_back(F.PositionOverTrapStrength[c]);
    }
    for(unsigned c=0;c<wannier0Functions.size();c++) {
        BpPositionOverTrapStrength[c].x.push_back(inhStrength);
        BpPositionOverTrapStrength[c].f.push_back(F.BpPositionOverTrapStrength[c]);
    }
}

unsigned findClosest(mt const &value,Col<mt> const &vec,mt &error) {
    unsigned ind; error=DBL_MAX;
    for(unsigned c=0;c<vec.n_elem;c++) {
        if(abs(value-vec[c])<error) {
            ind=c;
            error=abs(value-vec[c]);
        }
    }
    return ind;
}

// calculate eigenvalues in x_end and assign them to eigenvalues at x_st using linear interpolation
void linAssign(mt const start,         mt const end, 
        FsOverTS const &startF, FsOverTS       &endF) {
    /*const unsigned displ=128;
    for(unsigned c=0;                          c<round((start/pp.inhS)*displ); c++)
        fout << '|';
    for(unsigned c=round((start/pp.inhS)*displ); c<round((  end/pp.inhS)*displ); c++)
        fout << '-';
    for(unsigned c=round((  end/pp.inhS)*displ); c<displ;                      c++)
        fout << '|';
    fout << endl;*/
    Col<mt> bandSortedEndHeigenvalues(wannier0Functions.size());
    Col<mt> bandSortedEndHeigenvaluesDeriv(wannier0Functions.size());
    Col<mt> tmpEndHeigenvalues=endF.Heigenvalues;
    Col<mt> tmpEndHeigenvaluesDeriv=endF.HeigenvaluesDeriv;
    mt accuError=0;
    mt error;
    bool hit=true;

//  checkIfLinearInterpolationHitsEndValuesAndSort() {
        for(unsigned c=0;c<wannier0Functions.size();c++) {
            mt extrapolatedValue=startF.Heigenvalues[c]+(end-start)*startF.HeigenvaluesDeriv[c];
            unsigned endInd=findClosest(extrapolatedValue,tmpEndHeigenvalues,error);
            if(error>linAssignErrThreshold) {
                hit=false;
                break;
            }
            accuError+=error;
            bandSortedEndHeigenvalues[c]     =tmpEndHeigenvalues[endInd];
            bandSortedEndHeigenvaluesDeriv[c]=tmpEndHeigenvaluesDeriv[endInd];
            tmpEndHeigenvalues[endInd]=DBL_MAX;
        }
//  }

    if(hit) {
        endF.Heigenvalues     =bandSortedEndHeigenvalues;
        endF.HeigenvaluesDeriv=bandSortedEndHeigenvaluesDeriv;
    } else {
//      createMiddlePoints() {
            FsOverTS middleF=calcFsOverTS((start+end)/2);
            linAssign(start,(start+end)/2,startF,middleF);
            saveFsOverTS((start+end)/2,middleF);

            linAssign((start+end)/2,end,middleF,endF);
//      }
    }
}

    vector<fct<complex<mt> > > HeigenvaluesOverTrapStrength;
    vector<fct<complex<mt> > > HeigenvaluesOverTrapStrengthDeriv;
    vector<fct<complex<mt> > > PositionOverTrapStrength;
    vector<fct<complex<mt> > > BpPositionOverTrapStrength;
    mt linAssignErrThreshold;

    void calculateFOverTrapStrength() {
        if(HinW0partH0.n_elem==0)
            calculateHinW0();
        if(wannier0Eigenvalues.n_elem==0)
            calculateWannier0Eigenvalues();
        if(XinW0.n_elem==0)
            calculateXinW0();

        fout << "calculateFOverTrapStrength() executing ... "; fout.flush();
        tMark=clock();
        HeigenvaluesOverTrapStrength=vector<fct<complex<mt> > >(wannier0Functions.size());
        HeigenvaluesOverTrapStrengthDeriv=vector<fct<complex<mt> > >(wannier0Functions.size());
        PositionOverTrapStrength=vector<fct<complex<mt> > >(wannier0Functions.size());
        BpPositionOverTrapStrength=vector<fct<complex<mt> > >(wannier0Functions.size());
    
        mt inhStrengthStart=0.;
        mt inhStrengthEnd=pp.inhS;
        mt inhStrengthMaxErr=.01;
    
        linAssignErrThreshold=inhStrengthMaxErr;
    
        FsOverTS startF=calcFsOverTS(inhStrengthStart);
        saveFsOverTS(inhStrengthStart,startF);
    
        FsOverTS endF=calcFsOverTS(inhStrengthEnd);
        linAssign(inhStrengthStart,inhStrengthEnd,startF,endF);
        saveFsOverTS(inhStrengthEnd,endF);
    
        for(unsigned c=0;c<HeigenvaluesOverTrapStrength.size();c++) {
            HeigenvaluesOverTrapStrength[c].sort();
        }
        for(unsigned c=0;c<HeigenvaluesOverTrapStrengthDeriv.size();c++) {
            HeigenvaluesOverTrapStrengthDeriv[c].sort();
        }
    
/*    //  subtractLocalTrapPotential() {
            mt tmpTS=pp.inhS;
            for(unsigned c=0;c<HeigenvaluesOverTrapStrength.size();c++)
                for(unsigned d=0;d<HeigenvaluesOverTrapStrength[c].size();d++) {
                    pp.inhS=HeigenvaluesOverTrapStrength[c].x[d];
                    HeigenvaluesOverTrapStrength[c].f[d]-=Vinh(real(BpPositionOverTrapStrength[c].f[d]));
                }
            pp.inhS=tmpTS;
    //  }*/
    
        
//        outputFct(HeigenvaluesOverTrapStrength,/*"energy eigenvalues as function of inh strength*/"\tinh strength/E_r\tE/E_r");
        outputWithHeader(HeigenvaluesOverTrapStrength,"#\t\tinh pot strength/E_r\tE/E_r");
//        outputFct(HeigenvaluesOverTrapStrengthDeriv,/*"energy eigenvalues as function of inh strength*/"\tinh strength/E_r\tE'/E_r");
        outputWithHeader(HeigenvaluesOverTrapStrengthDeriv,"#\t\tinh strength/E_r\tE'/E_r");
//        outputFct(PositionOverTrapStrength,/*"energy eigenvalues as function of inh strength*/"\tinh strength/E_r\tx/a");
        outputWithHeader(PositionOverTrapStrength,"#\t\tinh strength/E_r\tx/a");
//        outputFct(BpPositionOverTrapStrength,/*"energy eigenvalues as function of inh strength*/"\tinh strength/E_r\tx_{{/Symbol=\\141}}/a");
        outputWithHeader(BpPositionOverTrapStrength,"#\t\tinh strength/E_r\tx_{{/Symbol=\\141}}/a");
        fout << (clock()-(float)tMark)/CLOCKS_PER_SEC << "s " << endl;
    }

    void calculateChangeInLocalizations() {
        if(HinW0partH0.n_elem==0)
            calculateHinW0();

        Mat<cx_double> HinW0(wannier0Functions.size(),wannier0Functions.size());
        Col<mt> HeigenvaluesESorted(wannier0Functions.size());
        Mat<cx_double> W0HOverlapESorted(wannier0Functions.size(),wannier0Functions.size());
        vector<fct<complex<mt> > > HeigenvaluesOverTrapStrength(wannier0Functions.size());
        vector<fct<mt> > PositionOverTrapStrength(wannier0Functions.size());
        fct<mt> QuaddedPositionOverTrapStrengthDeriv;
        mt inhFacStep=.0002;
        for(unsigned c=0;c<pp.inhS/inhFacStep;c++) {
            mt inhFac=c*inhFacStep;
            HinW0=HinW0partH0+HinW0partV/pp.inhS*inhFac;
            eig_sym(HeigenvaluesESorted,W0HOverlapESorted,HinW0);
            for(unsigned d=0;d<wannier0Functions.size();d++) {
                HeigenvaluesOverTrapStrength[d].x.push_back(inhFac);
                HeigenvaluesOverTrapStrength[d].f.push_back(HeigenvaluesESorted(d));
                PositionOverTrapStrength[d].x.push_back(inhFac);
                PositionOverTrapStrength[d].f.push_back(dot(square(abs(W0HOverlapESorted.col(d))),wannier0Eigenvalues));
            }
            if(c!=0) {
                mt sum=0;
                for(unsigned d=0;d<wannier0Functions.size();d++) {
                    mt diff=PositionOverTrapStrength[d].f[c]-PositionOverTrapStrength[d].f[c-1];
                    mt quo=diff/inhFacStep;
                    sum+=quo*quo;
                }
                QuaddedPositionOverTrapStrengthDeriv.x.push_back(inhFac);
                QuaddedPositionOverTrapStrengthDeriv.f.push_back(sum);
            }
        }
//        outputFct(HeigenvaluesOverTrapStrength,/*"energy eigenvalues as function of inh strength*/"\tinh strength/E_r\tE/E_r");
        outputWithHeader(HeigenvaluesOverTrapStrength,"#\t\tinh strength/E_r\tE/E_r");
//        outputFct(PositionOverTrapStrength,/*"energy eigenvalues as function of inh strength*/"\tinh strength/E_r\tx/a");
        outputWithHeader(PositionOverTrapStrength,"#\t\tinh strength/E_r\tx/a");
//        outputFct(QuaddedPositionOverTrapStrengthDeriv,/*"energy eigenvalues as function of inh strength*/"\tinh strength/E_r\tquadded position derivative/(a/E_r)^2");
        outputWithHeader(QuaddedPositionOverTrapStrengthDeriv,"#\t\tinh strength/E_r\tquadded position derivative/(a/E_r)^2");
    }

    void calculateLeftRightSideOfHinW0EigenvalueEquation() { // calculates left and right side of eigenvalue equation
        Mat<mt> HeigenvaluesOnDiagonal=createDiagonalMatrixFromVector(HeigenvaluesESorted);
        Mat<cx_double> HeigenvalueTimesEigenvector=W0HOverlapESorted*HeigenvaluesOnDiagonal;
        Mat<cx_double> HtimesEigenvector          =HinW0                      *W0HOverlapESorted;
//        output(HeigenvalueTimesEigenvector);
        outputWO(HeigenvalueTimesEigenvector);
//        output(HtimesEigenvector);
        outputWO(HtimesEigenvector);
    }

    void calculateProductOfHWannier0overlapMatrixWithItselvesAdjoint() { // ought to be unit matrix
        Mat<complex<mt> > productOfHWannier0overlapMatrixWithItselvesAdjoint=trans(W0HOverlapESorted)*W0HOverlapESorted;
//        output(productOfHWannier0overlapMatrixWithItselvesAdjoint);
        outputWO(productOfHWannier0overlapMatrixWithItselvesAdjoint);
    //        complex<mt> sum;
        for(unsigned e=0;e<W0HOverlapESorted.n_rows;e++) {
            for(unsigned d=0;d<W0HOverlapESorted.n_rows;d++) {
                complex<mt> sum=0.;
                for(unsigned c=0;c<W0HOverlapESorted.n_rows;c++) 
                    sum+=W0HOverlapESorted(c,d)*conj(W0HOverlapESorted(c,e));
                if(fabs(sum)<.1)                                             // fabs?!
                    sum=0.;
                fout << sum << ' ';
            }
            fout << endl;
        }
    }

    void calculateBandOrderedVariables() {
        if(bandToEindex.n_elem==0) {
            fout << __FILE__ << ':' << __LINE__ << ' ' << "bandToEindex not initialized" << endl;
            abort();
        }
        if(W0HOverlapESorted.n_elem==0)
            calculateHeigenvaluesESorted();
        if(HeigenvaluesESorted.n_elem==0)
            calculateHeigenvaluesESorted();

        Heigenvalues=Col<mt>(HeigenvaluesESorted.size());
        for(unsigned c=0;c<HeigenvaluesESorted.size();c++) {
            Heigenvalues[c]=HeigenvaluesESorted[bandToEindex[c]];
        }
//        output(Heigenvalues);
        outputWO(Heigenvalues);
        W0HOverlap=Mat<cx_double>(wannier0Functions.size(),wannier0Functions.size());
        for(unsigned c=0;c<wannier0Functions.size();c++)
            for(unsigned d=0;d<wannier0Functions.size();d++)
                W0HOverlap(c,d)=W0HOverlapESorted(c,bandToEindex[d]);
//        output(W0HOverlap);
        stringstream s; s<<"#\t\tn\tl_0+"<<statesPerBand<<"*{/Symbol=\\141}\t<l_0,{/Symbol=\\141}|n>\t\tset yrange [*:*] reverse;";
        outputWithHeader(W0HOverlap,s.str());
    }

    void calculateHeigenstates() {
        if(W0HOverlap.n_elem==0)
            calculateBandOrderedVariables();

        fout << "calculateHeigenstates() executing ... "; fout.flush();
        tMark=clock();
        for(unsigned d=0;d<W0HOverlap.n_rows;d++) {
            fct<complex<mt> > H=fct<complex<mt> >(-((double)nrSites/2)*pp.a,+((double)nrSites/2)*pp.a,sampleWidth)*0;
            for(unsigned c=0;c<wannier0Functions.size();c++) {
                H+=wannier0Functions[c]      *W0HOverlap(c,d);
                H+=wannier0FunctionsPeriod[c]*W0HOverlap(c,d);
            }
            H.normalize_with_L2_norm();
#ifdef PFSI
            H.prepare_for_spline_interpolation();
#endif
            Heigenfunctions.push_back(H);
        }

//for(unsigned c=0;c<wannier0Functions.size();c++)
#ifdef PFSI
///    Heigenfunctions[c].phi.prepare_for_spline_interpolation();
#endif

/*        for(unsigned d=0;d<wannier0Functions.size();d++) {
            Eeigenstate H;
            H.E=Heigenvalues[d];
            mt from=wannier0Functions[d].x[0];
            mt to  =wannier0Functions[d].x[wannier0Functions[d].x.size()-1];
            for(mt x=from;x<=to;x+=(to-from)/1000) {
                complex<mt> Res=0;
                for(unsigned c=0;c<wannier0Functions.size();c++)
                    Res+=W0HOverlap(c,d)*wannier0Functions[c](x);
                H.phi.x.push_back(x);
                H.phi.f.push_back(Res);
            }
            Heigenfunctions.push_back(H);
        }*/

//        outputFct(Heigenfunctions,/*"energy eigenfunctions of the lattice + inh system*/"\tx\t{/Symbol=\\152}(x)\t\twith lines");
        outputWithHeader(Heigenfunctions,"#\t\tx\t{/Symbol=\\152}(x)\t\twith lines");
        fout << (clock()-(float)tMark)/CLOCKS_PER_SEC << 's' << endl;
    }

    void outputSelectedFunctions() {
        if(Heigenfunctions.size()==0)
            calculateHeigenstates();

        vector<fct<complex<double> > > selectedFunctions;
        selectedFunctions.push_back(wannier0Functions[bandToIndex(0,       3)]);
        selectedFunctions.push_back(wannier0Functions[bandToIndex(0,       statesPerBand-1)]);
        selectedFunctions.push_back(wannier0Functions[bandToIndex(1,       3)]);
        selectedFunctions.push_back(wannier0Functions[bandToIndex(1,       statesPerBand-1)]);
        selectedFunctions.push_back(wannier0Functions[bandToIndex(2,       3)]);
        selectedFunctions.push_back(wannier0Functions[bandToIndex(2,       statesPerBand-1)]);
        selectedFunctions.push_back(Heigenfunctions[0]);
        selectedFunctions.push_back(Heigenfunctions[5]);
        selectedFunctions.push_back(Heigenfunctions[10]);
        selectedFunctions.push_back(Heigenfunctions[15]);
        selectedFunctions.push_back(Heigenfunctions[20]);

        stringstream Evals;
        Evals << "#\t\tx\t{/Symbol=\\152}(x)\t\twith lines" << endl;
        Evals << "E=" << setprecision(3) << Heigenvalues[ 0] << "\t" <<
                 "E=" << setprecision(3) << Heigenvalues[ 5] << "\t" <<
                 "E=" << setprecision(3) << Heigenvalues[10] << "\t" <<
                 "E=" << setprecision(3) << Heigenvalues[15] << "\t" <<
                 "E=" << setprecision(3) << Heigenvalues[20];
    
        selectedFunctions.push_back(Heigenfunctions[bandToIndex(1,statesPerBand/*/4*/-1)]);
//        outputFctWithLegend(selectedFunctions,/*"Heigenfunctions for k~-PI,-PI/2,0*/"\tx\t{/Symbol=\\152}(x)",Evals.str());
        outputWithHeader(selectedFunctions,Evals.str());
    }

    void compareWithHeigenstatesFromFile() {
        if(Heigenfunctions.size()==0)
            calculateHeigenstates();

        Col<mt> EDHeigenvalues(90);
        ifstream ifile;
        ifile.open("input/30sites_40res_20msize_V1E-16_s10Heigenvalues");
        for(unsigned c=0;c<90;c++) {
            mt val; ifile>>val;
            EDHeigenvalues[c]=val;
        }
        ifile.close();
    //  calculateOverlapWithHeigenfunctionsFromFile() {
    //      getWannierStatesFromFile() {
                unsigned EDHsamples=29*40+1;
                unsigned EDHstates=40;  //  581   1161
                vector<fct<complex<mt> > > EDHeigenfunctions(EDHstates);
                for(unsigned c=0;c<EDHstates;c++)
                    EDHeigenfunctions[c]=fct<complex<mt> >(-14.5,+14.5,(unsigned)EDHsamples);// Heigenfunctions[0]*0.;
                //ifstream ifile;
                ifile.open("input/30sites_40res_20msize_V1E-16_s10Heigenfunction");
                for(unsigned c=0;c<EDHsamples;c++) {
                    mt x=-14.5+(double)c*(14.5-(-14.5))/(EDHsamples-1);
                    for(unsigned d=0;d<EDHstates;d++) { // matrix_size=30 -> d<31 possible
                        mt ReF; ifile>>ReF;
    //                    mt ImF; ifile>>ImF;
    //                    complex<mt> f(ReF,ImF);
                        EDHeigenfunctions[d].x[c]=x;
                        EDHeigenfunctions[d].f[c]=ReF;
                    }
                }
                ifile.close();
    //      }

//            outputFct(EDHeigenfunctions,/*"Wannier functions of lattice-only system not rephased*/"\tx\t{/Symbol=\\152}@_l^{(0)}(x)\t\twith lines");
            outputWithHeader(EDHeigenfunctions,"#\t\tx\t{/Symbol=\\152}@_l^{(0)}(x)\t\twith lines");
        
            for_each(Heigenfunctions,normalizeFct);
#ifdef PFSI
            for(unsigned c=0;c<Heigenfunctions.size();c++)
                Heigenfunctions[c].prepare_for_spline_interpolation();
#endif //PFSI
            for_each(EDHeigenfunctions,normalizeFct);
#ifdef PFSI
            for(unsigned c=0;c<EDHeigenfunctions.size();c++)
                EDHeigenfunctions[c].prepare_for_spline_interpolation();
#endif //PFSI
    
            for(unsigned c=0;c<minBands;c++) {
                fout << c << ' ' << fabs(romberg_integral(Heigenfunctions[2*c+1].conj(),EDHeigenfunctions[c],
                                                     Heigenfunctions[2*c+1].x[0],Heigenfunctions[2*c+1].x[Heigenfunctions[2*c+1].x.size()-1],1E-8)) << endl;
            }
    //  }
        fout << "calculating Schroedinger equation deviation for EDHeigenfunction from file... "; fout.flush();
        tMark=clock();
        vector<fct<complex<mt> > > EDSchroedinger(EDHeigenfunctions.size());
        for(unsigned c=0;c<EDHeigenfunctions.size();c++) {
            EDSchroedinger[c]=EDHeigenfunctions[c]*0;
            fct<complex<mt> > firstDerivative=EDHeigenfunctions[c].deriv();
            fct<complex<mt> > secondDerivative=firstDerivative.deriv();
            for(unsigned d=0;d<EDHeigenfunctions[c].size();d++)
                EDSchroedinger[c][d]=(V(EDHeigenfunctions[c].x[d],pp)-EDHeigenvalues[c])*EDHeigenfunctions[c][d]-pp.hbar*pp.hbar/2./pp.m*secondDerivative[d];
        }
        fout << (clock()-(float)tMark)/CLOCKS_PER_SEC << 's' << endl;
//        outputFct(EDSchroedinger,"\tx/a\tdE*phi");
        outputWithHeader(EDSchroedinger,"#\t\tx/a\tdE*phi\t\twith lines");
        exit(0);

    }

    void calculateSchroedingerError() {
        if(Heigenfunctions.size()==0)
            calculateHeigenstates();

        fout << "calculateSchroedingerError() executing ... "; fout.flush();
        tMark=clock();
        vector<fct<complex<mt> > > Schroedinger(Heigenfunctions.size());
        Col<mt> SchroedingerError0(Heigenfunctions.size());
        //Col<mt> SchroedingerError1(Heigenfunctions.size());
        //Col<mt> SchroedingerError2(Heigenfunctions.size());
        Col<mt> SchroedingerBandError0(minBands);
        //Col<mt> SchroedingerBandError1(minBands);
        //Col<mt> SchroedingerBandError2(minBands);
        for(unsigned c=0;c<Heigenfunctions.size();c++) {
            Schroedinger[c]=Heigenfunctions[c]*0;
            free(Schroedinger[c].second_deriv_for_fcts_spline);
            Schroedinger[c].second_deriv_for_fcts_spline=NULL;
            fct<complex<mt> > firstDerivative=Heigenfunctions[c].deriv();
            fct<complex<mt> > secondDerivative=firstDerivative.deriv();
            for(unsigned d=4;d<Heigenfunctions[c].size()-4;d++)
                Schroedinger[c][d]=(V(Heigenfunctions[c].x[d],pp)-Heigenvalues[c])*Heigenfunctions[c][d]-pp.hbar*pp.hbar/2./pp.m*secondDerivative[d];
            SchroedingerError0[c]=sqrt(romberg_integral(Schroedinger[c].L2_norm_squared_elementwise(),Schroedinger[c].x[4],Schroedinger[c].x[Schroedinger[c].x.size()-5],1E-6));
//            SchroedingerError1[c]=sqrt(romberg_integral(Schroedinger[c].L2_norm_squared_elementwise(),-((double)nrSites/2-1)*pp.a,+((double)nrSites/2-1)*pp.a,1E-6));
//            SchroedingerError2[c]=sqrt(romberg_integral(Schroedinger[c].L2_norm_squared_elementwise(),-((double)nrSites/2-2)*pp.a,+((double)nrSites/2-2)*pp.a,1E-6));
        }
        for(unsigned c=0;c<minBands;c++) {
            mt sum=0;
            for(unsigned d=0;d<statesPerBand;d++)
                sum+=SchroedingerError0[c*statesPerBand+d];
            SchroedingerBandError0[c]=sum/statesPerBand;
        }
/*        for(unsigned c=0;c<minBands;c++) {
            mt sum=0;
            for(unsigned d=0;d<statesPerBand;d++)
                sum+=SchroedingerError1[c*statesPerBand+d];
            SchroedingerBandError1[c]=sum/statesPerBand;
        }
        for(unsigned c=0;c<minBands;c++) {
            mt sum=0;
            for(unsigned d=0;d<statesPerBand;d++)
                sum+=SchroedingerError2[c*statesPerBand+d];
            SchroedingerBandError2[c]=sum/statesPerBand;
        }*/
        fout << (clock()-(float)tMark)/CLOCKS_PER_SEC << "s " << endl;
//        output(SchroedingerError0);
        outputWO(SchroedingerError0);
////        output(SchroedingerError1);
//        outputWO(SchroedingerError1);
////        output(SchroedingerError2);
//        outputWO(SchroedingerError2);
//        output(SchroedingerBandError0);
        outputWO(SchroedingerBandError0);
////        output(SchroedingerBandError1);
//        outputWO(SchroedingerBandError1);
////        output(SchroedingerBandError2);
//        outputWO(SchroedingerBandError2);
        for(unsigned c=0;c<minBands;c++)
            fout << SchroedingerBandError0[c] << ' ';
        fout << endl;
/*        for(unsigned c=0;c<minBands;c++)
            fout << SchroedingerBandError1[c] << ' ';
        fout << endl;
        for(unsigned c=0;c<minBands;c++)
            fout << SchroedingerBandError2[c] << ' ';
        fout << endl;*/
//        outputFct(Schroedinger,"\tx/a\tdE*phi");
        outputWithHeader(Schroedinger,"#\t\tx/a\tdE*phi\t\twith lines");
    }

    void calculateOverlapWithHarmonicOscillator() {
        if(Heigenfunctions.size()==0)
            calculateHeigenstates();

        unsigned nrSt=10;
        vector<fct<complex<mt> > > hO(nrSt);
        for(unsigned c=0;c<nrSt;c++) {
            hO[c]=Heigenfunctions[0];
            //hO[c]=fct<complex<mt> >(-5,5,.1);

            mt omega=sqrt(2*pp.inhS/pp.m);
            mt nfak=1;
            for(unsigned d=2;d<=c;d++)
                nfak*=d;
            mt firstFac=1/sqrt(pow(2,c)*nfak); // oneOverSqrt2ToTheNTimesNfak
            mt secondFac=sqrt(sqrt(pp.m*omega/PI/pp.hbar));
            for(unsigned d=0;d<hO[c].size();d++) {
                mt x=hO[c].x[d];
                mt thirdFac=exp(-pp.m*omega*x*x/2/pp.hbar);
                mt hermArg=sqrt(pp.m*omega/pp.hbar)*x;
                hO[c][d]=firstFac*secondFac*thirdFac*boost::math::hermite(c,hermArg);
            }
            hO[c].normalize_with_L2_norm();
#ifdef PFSI
            hO[c].prepare_for_spline_interpolation();
#endif //PFSI
        }
//        outputFct(hO,"\tx/a\tf(x)");
        outputWithHeader(hO,"#\t\tx/a\tf(x)\t\twith lines");

        Mat<mt> hOoverlap(nrSt,nrSt);
        for(unsigned c=0;c<nrSt;c++)
            for(unsigned d=0;d<nrSt;d++)
                hOoverlap(c,d)=norm(romberg_integral_f(Heigenfunctions[c].conj(),hO[d],constantOp,1E-6));
//        output(hOoverlap);
        outputWO(hOoverlap);
    }

    void calculateWinE() {
        if(Heigenfunctions.size()==0)
            calculateHeigenstates();

        fout << "calculateWinE() executing ... "; fout.flush();
        tMark=clock();
        wannierInH=Mat<cx_double>(Heigenfunctions.size(),Heigenfunctions.size());
        for(unsigned c=0;c<Heigenfunctions.size();c++) {
            fct<complex<mt> > phiXconj( (Heigenfunctions[c] *
                                         Heigenfunctions[c].x_fct()).conj() );
#ifdef PFSI
            phiXconj.prepare_for_spline_interpolation();
#endif //PFSI
            for(unsigned d=0;d<Heigenfunctions.size();d++) {
                //fout << c << ' ' << d << endl;
                //if(d<c)
                //    wannierInH(c,d)=wannierInH(d,c);                  // OPTIMIZATION 
                //else
                if(c/statesPerBand==d/statesPerBand)
                    wannierInH(c,d)=romberg_integral(phiXconj,
                                              Heigenfunctions[d],
                                              Heigenfunctions[c].x[0],
                                              Heigenfunctions[c].x[Heigenfunctions[c].size()-1], 1E-7);
                else
                    wannierInH(c,d)=0;
            }
        }
        fout << (clock()-(float)tMark)/CLOCKS_PER_SEC << 's' << endl;

//        output(wannierInH);
        outputWithHeader(wannierInH,"#\t\tn'\tn\t<n|X_{{/Symbol=\\141}}|n'>\t\tset yrange [*:*] reverse;");
    }

    void calculateWannierEigenvalues() {
        if(wannierInH.n_elem==0)
            calculateWinE();

        wannierEigenvalues=Col<mt>(Heigenfunctions.size());
        fout << "calculateWannierEigenvalues() executing ... "; fout.flush();
        tMark=clock();
        HWannierOverlap=calculateEigenstatesBlockwise(wannierInH, wannierEigenvalues, statesPerBand);
        fout << (clock()-(float)tMark)/CLOCKS_PER_SEC << 's' << endl;

//        output(wannierEigenvalues);
        outputWO(wannierEigenvalues);
//        output(HWannierOverlap);
        stringstream s; s<<"#\t\tl+"<<statesPerBand<<"*{/Symbol=\\141}\tn\t<n|l,{/Symbol=\\141}>\t\tset yrange [*:*] reverse;";
        outputWithHeader(HWannierOverlap,s.str());
    }

    void calculateLeftRightSideOfWannierInHEigenvalueEquation() { // calculates left and right side of eigenvalue equation
        Mat<mt> wannierEigenvaluesOnDiagonal=createDiagonalMatrixFromVector(wannierEigenvalues);
        Mat<cx_double> wannierEigenvaluesTimesEigenvectors=HWannierOverlap*wannierEigenvaluesOnDiagonal;
        Mat<cx_double> wannierTimesEigenvectors        =wannierInH           *HWannierOverlap;
//        output(wannierEigenvaluesTimesEigenvectors);
        outputWO(wannierEigenvaluesTimesEigenvectors);
//        output(wannierTimesEigenvectors);
        outputWO(wannierTimesEigenvectors);
    }

    void calculateProductOfWannierHoverlapMatrixWithItselvesAdjoint() {
        Mat<complex<mt> > productOfWannierHoverlapMatrixWithItselvesAdjoint=trans(HWannierOverlap)*HWannierOverlap;
//        output(productOfWannierHoverlapMatrixWithItselvesAdjoint);
        outputWO(productOfWannierHoverlapMatrixWithItselvesAdjoint);
        //complex<mt> sum;
        for(unsigned e=0;e<HWannierOverlap.n_rows;e++) {
            for(unsigned d=0;d<HWannierOverlap.n_rows;d++) {
                complex<mt> sum=0.;
                for(unsigned c=0;c<HWannierOverlap.n_rows;c++)
                    sum+=HWannierOverlap(c,d)*conj(HWannierOverlap(c,e));
                if(fabs(sum)<.1)                                                   //fabs?!
                    sum=0.;
                fout << sum << ' ';
            }
            fout << endl;
        }
    }

    void calculateWannierFunctions() {
        if(Heigenfunctions.size()==0)
            calculateHeigenstates();
        if(HWannierOverlap.n_elem==0)
            calculateWannierEigenvalues();

        fout << "calculateWannierFunctions() executing ... "; fout.flush();
        tMark=clock();
        for(unsigned c=0;c<Heigenfunctions.size();c++) {
            fct<complex<mt> > wan=Heigenfunctions[(c/statesPerBand)*statesPerBand]*HWannierOverlap((c/statesPerBand)*statesPerBand,c);
            for(unsigned d=1;d<statesPerBand;d++)
                wan+=Heigenfunctions[(c/statesPerBand)*statesPerBand+d]*HWannierOverlap((c/statesPerBand)*statesPerBand+d,c);
#ifdef PFSI
            wan.prepare_for_spline_interpolation();
#endif
            wannierFunctions.push_back(wan);
        }
        fout << (clock()-(float)tMark)/CLOCKS_PER_SEC << 's' << endl;
//        outputFct(wannierFunctions,/*"Wannier functions of the lattice + inh system not rephased*/"\tx\t{/Symbol=\\152}_l(x)");
        outputWithHeader(wannierFunctions,"#\t\tx\t{/Symbol=\\152}_l(x)\t\twith lines");
    }

    void calculateFunctionsRephased() {
        if(wannierFunctions.size()==0)
            calculateWannierFunctions();
        if(HWannierOverlap.n_elem==0)
            calculateWannierEigenvalues();

        wannierFunctionsRephased=wannierFunctions;
        HWannierOverlapRephased=Mat<cx_double>(Heigenfunctions.size(),wannierFunctions.size());
        HeigenfunctionsRephased=Heigenfunctions;

        for(unsigned c=0;c<wannierFunctionsRephased.size();c++) {
            complex<mt> factor;
            mt max=0;
            unsigned ind;
            if((c/statesPerBand)%2==0) {
                for(unsigned d=0;d<wannierFunctionsRephased[c].size();d++)
                    if(fabs(imag(wannierFunctionsRephased[c].f[d]))>max) {
                        max=fabs(imag(wannierFunctionsRephased[c].f[d]));
                        ind=d;
                    }
                factor=fabs(wannierFunctionsRephased[c].f[ind])/wannierFunctionsRephased[c].f[ind];
            } else {
                fct<complex<mt> > deriv=wannierFunctionsRephased[c].deriv();
                for(unsigned d=0;d<wannierFunctionsRephased[c].size();d++)
                    if(fabs(imag(deriv.f[d]))>max) {
                        max=fabs(imag(deriv.f[d]));
                        ind=d;
                    }
                factor=fabs(deriv.f[ind])/deriv.f[ind];
            }
            for(unsigned d=0;d<wannierFunctionsRephased[c].size();d++) {
                wannierFunctionsRephased[c].f[d]*=factor;
            }
            for(unsigned d=0;d<Heigenfunctions.size();d++)
                HWannierOverlapRephased(d,c)=HWannierOverlap(d,c)*factor;
        }
        for(unsigned c=0;c<HeigenfunctionsRephased.size();c++) {
            complex<mt> factor;
            mt max=0;
            unsigned ind;
            if((c/statesPerBand)%2==0) {
                for(unsigned d=0;d<HeigenfunctionsRephased[c].size();d++)
                    if(fabs(HeigenfunctionsRephased[c].f[d])>max) {
                        max=fabs(HeigenfunctionsRephased[c].f[d]);
                        ind=d;
                    }
                factor=fabs(HeigenfunctionsRephased[c].f[ind])/HeigenfunctionsRephased[c].f[ind];
            } else {
                fct<complex<mt> > deriv=HeigenfunctionsRephased[c].deriv();
                for(unsigned d=0;d<HeigenfunctionsRephased[c].size();d++)
                    if(fabs(imag(deriv.f[d]))>max) {
                        max=fabs(imag(deriv.f[d]));
                        ind=d;
                    }
                factor=fabs(deriv.f[ind])/deriv.f[ind];
            }
            for(unsigned d=0;d<HeigenfunctionsRephased[c].size();d++) {
                HeigenfunctionsRephased[c].f[d]*=factor;
            }
        }

//        for_each( wannierFunctionsRephased , addPhaseToLocalizeOnRealAxis );
//        outputFct(wannierFunctionsRephased,/*"Wannier functions of lattice + inh system*/"\tx\t{/Symbol=\\152}@_l(x)");
        outputWithHeader(wannierFunctionsRephased,"#\t\tx\t{/Symbol=\\152}@_l(x)\t\twith lines");
//        output(HWannierOverlapRephased);
        stringstream s; s<<"#\t\tl+"<<statesPerBand<<"*{/Symbol=\\141}\tn\t<n|l,{/Symbol=\\141}>\t\tset yrange [*:*] reverse;";
        outputWithHeader(HWannierOverlapRephased,s.str());
        outputWithHeader(HeigenfunctionsRephased,"#\t\tx\t{/Symbol=\\152}(x)\t\twith lines");
    }

    void calculateHopping() {
        if(HWannierOverlapRephased.n_elem==0)
            calculateFunctionsRephased();
        if(Heigenvalues.n_elem==0)
            calculateBandOrderedVariables();

        hopping=Mat<cx_double>(wannierFunctions.size(),wannierFunctions.size());
        Mat<double> absHopping(wannierFunctions.size(),wannierFunctions.size());
        for(unsigned c=0;c<wannierFunctions.size();c++) {
            for(unsigned d=0;d<wannierFunctions.size();d++) {
                complex<mt> sum=0;
                for(unsigned e=0;e<Heigenfunctions.size();e++) {
                    sum+=Heigenvalues[e]*conj(HWannierOverlapRephased(e,c))*HWannierOverlapRephased(e,d);
                }
                hopping(c,d)=sum;
                absHopping(c,d)=abs(sum);
            }
        }
//        output(hopping);
        stringstream s;  s << "#\t\tl'+" << statesPerBand << "*{/Symbol=\\141}'\tl+"<<statesPerBand<<"*{/Symbol=\\141}\t<l,{/Symbol=\\141}|H|l',{/Symbol=\\141}'>\t\tset yrange [*:*] reverse;";
        outputWithHeader(hopping,s.str());
//        output(absHopping);
        stringstream s2; s2<< "#\t\tl'+" << statesPerBand << "*{/Symbol=\\141}'\tl+"<<statesPerBand<<"*{/Symbol=\\141}\t|<l,{/Symbol=\\141}|H|l',{/Symbol=\\141}'>|\t\tset yrange [*:*] reverse;";
        outputWithHeader(absHopping,s2.str());

        unsigned neighbors=4;
        unsigned bands=2;
        vector<fct<complex<mt> > > hoppingF(neighbors*bands);
        vector<fct<complex<mt> > > absHoppingF(neighbors*bands);
        for(unsigned b=0;b<bands;b++) {
            for(unsigned c=0;c<neighbors;c++) {
                for(unsigned d=2;d<statesPerBand-2-c;d++) {
                       hoppingF[c+b*neighbors].x.push_back((mt)d-statesPerBand/2+.5);
                       hoppingF[c+b*neighbors].f.push_back(   hopping(d+b*statesPerBand,(d+c)%statesPerBand+b*statesPerBand));
                    absHoppingF[c+b*neighbors].x.push_back((mt)d-statesPerBand/2+.5);
                    absHoppingF[c+b*neighbors].f.push_back(absHopping(d+b*statesPerBand,(d+c)%statesPerBand+b*statesPerBand));
                }
            }
        }
//        outputFct(hoppingF[0],"\tsite index\thopping/E_r");
        outputWithHeader(hoppingF[0],"#\t\tsite index\thopping/E_r");
//        outputFct(absHoppingF[0],"\tsite index\thopping/E_r");
        outputWithHeader(absHoppingF[0],"#\t\tsite index\thopping/E_r");
//        outputFct(hoppingF[1],"\tsite index\thopping/E_r");
        outputWithHeader(hoppingF[1],"#\t\tsite index\thopping/E_r");
//        outputFct(absHoppingF[1],"\tsite index\thopping/E_r");
        outputWithHeader(absHoppingF[1],"#\t\tsite index\thopping/E_r");
//        outputFct(hoppingF[2],"\tsite index\thopping/E_r");
        outputWithHeader(hoppingF[2],"#\t\tsite index\thopping/E_r");
//        outputFct(absHoppingF[2],"\tsite index\thopping/E_r");
        outputWithHeader(absHoppingF[2],"#\t\tsite index\thopping/E_r");
//        outputFct(hoppingF[3],"\tsite index\thopping/E_r");
        outputWithHeader(hoppingF[3],"#\t\tsite index\thopping/E_r");
//        outputFct(absHoppingF[3],"\tsite index\thopping/E_r");
        outputWithHeader(absHoppingF[3],"#\t\tsite index\thopping/E_r");
//        outputFct(hoppingF[4],"\tsite index\thopping/E_r");
        outputWithHeader(hoppingF[4],"#\t\tsite index\thopping/E_r");
//        outputFct(absHoppingF[4],"\tsite index\thopping/E_r");
        outputWithHeader(absHoppingF[4],"#\t\tsite index\thopping/E_r");
//        outputFct(hoppingF[5],"\tsite index\thopping/E_r");
        outputWithHeader(hoppingF[5],"#\t\tsite index\thopping/E_r");
//        outputFct(absHoppingF[5],"\tsite index\thopping/E_r");
        outputWithHeader(absHoppingF[5],"#\t\tsite index\thopping/E_r");
//        outputFct(hoppingF[6],"\tsite index\thopping/E_r");
        outputWithHeader(hoppingF[6],"#\t\tsite index\thopping/E_r");
//        outputFct(absHoppingF[6],"\tsite index\thopping/E_r");
        outputWithHeader(absHoppingF[6],"#\t\tsite index\thopping/E_r");
//        outputFct(hoppingF[7],"\tsite index\thopping/E_r");
        outputWithHeader(hoppingF[7],"#\t\tsite index\thopping/E_r");
//        outputFct(absHoppingF[7],"\tsite index\thopping/E_r");
        outputWithHeader(absHoppingF[7],"#\t\tsite index\thopping/E_r");
    }

    void calculateU() {
        if(wannierFunctions.size()==0)
            calculateWannierFunctions();

        fct<complex<mt> > U;
        for(unsigned c=0;c<statesPerBand;c++) {
            fct<mt> wanSquared=wannierFunctions[c].L2_norm_squared_elementwise();
#ifdef PFSI
            wanSquared.prepare_for_spline_interpolation();
#endif //PFSI
            U.x.push_back((mt)c-statesPerBand/2+.5);
            U.f.push_back(pp.g/2*romberg_integral_f(wanSquared,wanSquared,constantOp, 1E-6));
        }
//        outputFct(U,"\tsite index\tU/E_r");
        outputWithHeader(U,"#\t\tl\tU/E_r");
    }

    void calculateHeigenstatesFromWannier() { // recalculated former Eigenstates from Wannier states
        vector<fct<complex<mt> > > recalcHeigenfunctions=wannierFunctions;
        for(unsigned c=0;c<wannierFunctions.size();c++) {
            recalcHeigenfunctions[c]=wannierFunctions[(c/statesPerBand)*statesPerBand]*conj(HWannierOverlap(c,(c/statesPerBand)*statesPerBand));
            for(unsigned d=1;d<statesPerBand;d++)
                recalcHeigenfunctions[c]+=wannierFunctions[(c/statesPerBand)*statesPerBand+d]*conj(HWannierOverlap(c,(c/statesPerBand)*statesPerBand+d));
        }
    
//        outputFct(recalcHeigenfunctions,/*"recalculated energy eigenfunctions from Wannier functions*/"\tx\t{/Symbol=\\152}(x)\t\twith lines");
        outputWithHeader(recalcHeigenfunctions,"#\t\tx\t{/Symbol=\\152}(x)\t\twith lines");
    }

    void calculateWannierFunctionValueUsingHeigenfunctions() { // wannier eigenvalue check: this function should output the same two numbers
        for(unsigned WEClIndex=0;WEClIndex<statesPerBand;WEClIndex++) {
            unsigned WECband=0;
            mt WECpos=-2.062;
            Mat<complex<mt> > xWeightedHstateOverlap(statesPerBand,statesPerBand);
            for(unsigned c=0;c<statesPerBand;c++) {
                fct<complex<mt> > blochXconj=Heigenfunctions[WECband*statesPerBand+c].conj();
                blochXconj*=blochXconj.x_fct();
                for(unsigned d=0;d<statesPerBand;d++)
                    xWeightedHstateOverlap(c,d)=romberg_integral(blochXconj, Heigenfunctions[WECband*statesPerBand+d], blochXconj.x[0], blochXconj.x[blochXconj.x.size()-1], 1E-6);
            }
    
        //    fout << xWeightedHstateOverlap << endl;
            complex<mt> WECsum=0;
            for(unsigned c=0;c<statesPerBand;c++)
                for(unsigned d=0;d<statesPerBand;d++)
                    WECsum+=xWeightedHstateOverlap(c,d)*Heigenfunctions[WECband*statesPerBand+c](WECpos)*HWannierOverlap(WECband*statesPerBand+d,WECband*statesPerBand+WEClIndex);
            WECsum/=(wannierEigenvalues[WECband*statesPerBand+WEClIndex]*pp.a);
            fout << WECsum << ' ' << wannierFunctions[WECband*statesPerBand+WEClIndex](WECpos) << endl;
        }
    }

    void calculateWannierOverlap() { // wannierFunction overlap
        if(wannierFunctions.size()==0)
            calculateWannierFunctions();

        fout << "calculateWannierOverlap() ... "; fout.flush();
        tMark=clock();
        Mat<cx_double> wannierStateOverlapMatrix(wannierFunctions.size(),wannierFunctions.size());
        wannierStateOverlapMatrix.zeros();
        for(unsigned c=0;c<wannierFunctions.size();c++)
            wannierStateOverlapMatrix(1,c)=abs(romberg_integral(wannierFunctions[1].conj(), wannierFunctions[c], wannierFunctions[1].x[0], +wannierFunctions[1].x[wannierFunctions[1].size()-1], 1E-6));
        for(unsigned c=0;c<wannierFunctions.size();c++)
            wannierStateOverlapMatrix(2,c)=abs(romberg_integral(wannierFunctions[2].conj(), wannierFunctions[c], wannierFunctions[2].x[0], +wannierFunctions[2].x[wannierFunctions[2].size()-1], 1E-6));
        //Mat<cx_double> wannierStateOverlapMatrix=calculateOverlapMatrix(wannierFunctions);
        fout << (clock()-(float)tMark)/CLOCKS_PER_SEC << 's' << endl;
//        output(wannierStateOverlapMatrix);
        outputWithHeader(wannierStateOverlapMatrix,"#\t\tl'\tl\t<l|l'>\t\tset yrange [*:*] reverse;");
        outputWO(wannierInH);
    }

    void calculateBorderEffects() {
        vector<fct<complex<mt> > > borderEffects;
        for(unsigned c=0;c<wannierFunctions.size();c++) {
            fct<complex<mt> > bE(wannierFunctions[c]);
            for(unsigned d=0;d<bE.size();d++) {
                bE.f[d]=fabs(bE.f[d]);
                mt x=bE.x[d];
                int localizedAt=.5+c-(wannierFunctions.size()/2);
                mt equivSpotOnWannierFunctionInLatticeCenter=x-localizedAt;
                if(fabs(equivSpotOnWannierFunctionInLatticeCenter)<pp.s_bound)
                    bE.f[d]-=fabs(wannierFunctions[wannierFunctions.size()/2](equivSpotOnWannierFunctionInLatticeCenter));
            }
            borderEffects.push_back(bE);
        }
//        outputPlot(borderEffects,"\tx/a");
        outputWithHeader(borderEffects,"#\t\tx/a");
    }



    void calc() {
        bool exCheck=false;

        //outputPropagatedWavefunctionEndValueOverE();
        outputPotential0();
        outputPotential();

        statesPerBand0=kSteps;

        calculateH0eigenvalues();
        calculateH0eigenfunctions();
        //calculateOverlapWithH0eigenfunctionsFromFile();
        calculateSchroedinger0Error();
        //if(exCheck) H0H0Overlap();
        //hopping0OverDl();

        calculateCenteredWannier0Functions();
        calculateWannier0Functions();
        //calculateWannier0Overlap();
        //calculateGaussianInWannier0Basis();
        //calculateOverlapWithWannierStatesFromFile();

        statesPerBand=nrSites;

        //hopping0viaInt();
        calculateHinW0();

        //calculateH0inW0fromW0inH0();

        calculateHeigenvaluesESorted();
        calculateHeigenvaluesOverTrapStrength();
        //calculateXinW0();
        //calculateHlocalization();

        // execute either of the following to calculate bandToEindex
        calculateBand0ProjectorsExpectationInHeigenstates();
        assignBandOrderUsingLocalization();
   
        //calculateWannier0Eigenvalues();
        calculateFOverTrapStrength();
        //calculateChangeInLocalizations();
        assignBandOrderEqualHeigenvalueOrder();
        
        if(exCheck) calculateLeftRightSideOfHinW0EigenvalueEquation();
        if(exCheck) calculateProductOfHWannier0overlapMatrixWithItselvesAdjoint();

        calculateBandOrderedVariables();
        calculateHeigenstates();
        //outputSelectedFunctions();

        //compareWithHeigenstatesFromFile();
        calculateSchroedingerError();
        //calculateOverlapWithHarmonicOscillator();

        calculateWinE();
        calculateWannierEigenvalues();

        if(exCheck) calculateLeftRightSideOfWannierInHEigenvalueEquation();
        if(exCheck) calculateProductOfWannierHoverlapMatrixWithItselvesAdjoint();

        calculateWannierFunctions();
        calculateFunctionsRephased();
        calculateHopping();
        calculateU();

//        if(exCheck) calculateHeigenstatesFromWannier();
//        if(exCheck) calculateWannierFunctionValueUsingHeigenfunctions();
//        //calculateWannierOverlap();
//
//        //if(true) calculateBorderEffects();
//        //vector<mt> variances=calculateVariances(wannierFunctions);
//        //output(variances);*/
//
//        vector<fct<complex<double> > > selectedWannierFunctions;
//        selectedWannierFunctions.push_back(potential);
//        //for(unsigned c=0;c<6;c++)
//        //    selectedWannierFunctions.push_back(Heigenfunctions[c]+Heigenvalues[c]);
//        for(unsigned c=0;c<18;c++)
//            selectedWannierFunctions.push_back(wannierFunctionsRephased[c]+hopping(c,c));
////        outputFctWithLegend(selectedWannierFunctions,/*"Heigenfunctions for k~-PI,-PI/2,0*/"\tx\t{/Symbol=\\152}(x)",Evals.str());
//        outputWithHeader(selectedWannierFunctions,"#\t\tx\tV(x)/E_r\t\twith lines");
//
//        vector<fct<complex<double> > > selectedHeigenFunctions;
//        selectedHeigenFunctions.push_back(potential);
//        for(unsigned c=0;c<18;c++)
//            selectedHeigenFunctions.push_back(HeigenfunctionsRephased[c]+Heigenvalues[c]);
//        //for(unsigned c=0;c<6;c++)
//        //    selectedHeigenFunctions.push_back(wannierFunctions[c]+hopping(c,c));
////        outputFctWithLegend(selectedHeigenFunctions,/*"Heigenfunctions for k~-PI,-PI/2,0*/"\tx\t{/Symbol=\\152}(x)",Evals.str());
//        outputWithHeader(selectedHeigenFunctions,"#\t\tx\tV(x)/E_r\t\twith lines");
    }
};


const char *argp_program_version =
"inWan 1.0";
const char *argp_program_bug_address =
"<jonathanenders@gmx.de>";
     
/* Program documentation. */
static char doc[] =
    "inWan -- calculates Wannier states for inhomogeneous lattices";
 
/* A description of the arguments we accept. */
static char args_doc[] = "";
 
/* The options we understand. */
static struct argp_option options[] =  {
    {"sites",        's', "SIZE",      0,  "amount of lattice sites" },
    {"evIntervalEnd",'b', "NUMBER",    0,  "eigenvalue interval end, corresponds to amount of bands" },
    {"latticeDepth", 'l', "NUMBER",    0,  "depth of the lattice in E_r" },
    {"inhStrength",  'i', "NUMBER",    0,  "strength of the inhomogenic potential in E_r" },
    {"inhShift",     'h', "NUMBER",    0,  "shift between lattice and inhomogenic potential in sites" },
    {"periodRatio",  'r', "NUMBER",    0,  "ratio between the lattice periodicities" },
    {"outPath",      'o', "PATH",      0,  "relative path where to store output" },
    { 0 }
};
 
/* Used by main to communicate with parse_opt. */
struct arguments {
    unsigned sites;
    double latticeDepth, inhStrength, inhShift, periodRatio, evIntervalEnd;
    char* outPath;
};
 
 /* Parse a single option. */
static error_t
parse_opt (int key, char *arg, struct argp_state *state) {
   /* Get the input argument from argp_parse, which we
      know is a pointer to our arguments structure. */
    struct arguments *arguments = (struct arguments*)state->input;
 
    switch (key) {
        case 's':
            sscanf(arg,"%u",&(arguments->sites));
            break;
        case 'b':
            sscanf(arg,"%lf",&(arguments->evIntervalEnd));
            break;
        case 'l':
            sscanf(arg,"%lf",&(arguments->latticeDepth));
            break;
        case 'i':
            sscanf(arg,"%lf",&(arguments->inhStrength));
            break;
        case 'h':
            sscanf(arg,"%lf",&(arguments->inhShift));
            break;
        case 'r':
            sscanf(arg,"%lf",&(arguments->periodRatio));
            break;
        case 'o':
            arguments->outPath=arg;
            break;
 
//     case ARGP_KEY_ARG:
//       if (state->arg_num >= 3)
//         /* Too many arguments. */
//         argp_usage (state);
// 
//       break;
// 
//     case ARGP_KEY_END:
//       if (state->arg_num < 2)
//         /* Not enough arguments. */
//         argp_usage (state);
//       break;
 
     default:
       return ARGP_ERR_UNKNOWN;
     }
   return 0;
}
 
/* Our argp parser. */
static struct argp argp = { options, parse_opt, args_doc, doc };
 

int main(int argc, char** argv) {
    gsl_set_error_handler_off();

    struct arguments arguments;
    arguments.sites=6;
    arguments.evIntervalEnd=6.;  //(finds only first band states)   //15.(finds second band states)
    arguments.latticeDepth=10;
    arguments.inhStrength=2.0;//.6355;//0.05;//0.0240019;
    arguments.inhShift=.3;
    arguments.periodRatio=(sqrt(5)-1)/2;
    arguments.outPath=new char[8];
    strcpy(arguments.outPath,"output");
    argp_parse (&argp, argc, argv, 0, 0, &arguments);
    dirName=string(arguments.outPath);
    fout ex(arguments.sites) ex(arguments.latticeDepth) ex(arguments.inhStrength) ex(arguments.inhShift) ex(arguments.periodRatio) ex(arguments.evIntervalEnd) ex(arguments.outPath) << endl;

    if(mkdir(dirName.c_str(),0770)!=0) {
        fout << "inWan: couldn't create folder \"" << dirName << '\"' << endl;
        return(1);
    }

    // copy code into output folder
    // makes more sense to do that in Makefile because code isn't necessarily available during execution time
    /*ifstream codeFileIn("inWan.cpp", fstream::binary);
    ofstream codeFileOut((dirName+"/code.cpp").c_str(), fstream::trunc|fstream::binary);
    codeFileOut << codeFileIn.rdbuf();*/

    fout.open((dirName+"/stdout").c_str());

    inhLatticeSystem tls;

    pp.hbar=1.;
    pp.E_r=1.;
    pp.a=1.; //=1064nm;
    pp.lambda=2.*pp.a;
    pp.k_L=2.*PI/pp.lambda;
    pp.m=pp.hbar*pp.hbar*pp.k_L*pp.k_L/(2*pp.E_r);
    pp.s_bound=.5*pp.a;
    pp.k=0;
    pp.ld=arguments.latticeDepth;    // lattice depth higher 30 needs a lower intervalStart
    pp.inhShift=arguments.inhShift;
    pp.a_0=4.97180451127819548872e-5; //5.29e-11m/1064e-9m;
    pp.a_s=100.40*pp.a_0;
    pp.g=4*PI*pp.hbar*pp.hbar*pp.a_s/pp.m;

    pp.alpha=arguments.periodRatio;
    
    tls.sampleWidth=1./160./pp.a;//1./100/pp.a;
    tls.odeSampleWidth=1./160./pp.a;//2.*pp.s_bound/10000;
    tls.wannier0FunctionCutoff=1E-3;
    tls.HinW0cutoff=1E-3;

    tls.kSteps=arguments.sites;     // also the number of Bloch states used for calculation of W0func
    tls.nrSites=tls.kSteps;    // number of sites in the inhomogenic system

/*  pp.inhLatRatio=.02;
    pp.borderTrapHeight=pp.inhLatRatio*pp.ld;
    pp.inhS=pp.borderTrapHeight/(tls.nrSites/2*tls.nrSites/2);*/

/*  pp.borderTrapHeight=20.;
    pp.inhS=pp.borderTrapHeight/sqrt(tls.nrSites/2);
    pp.inhLapRatio=pp.borderTrapHeight/pp.ld;*/

    pp.inhS=arguments.inhStrength;//0.0240019;//.0;//.032566;//0.2;//0.22;//0.25;//.5;//1E-16;//.0001;
    pp.borderTrapHeight=pp.inhS*(tls.nrSites/2)*(tls.nrSites/2);
    pp.inhLatRatio=pp.borderTrapHeight/pp.ld;

    fout ex(pp.ld) ex(pp.inhS) ex(pp.inhLatRatio) << endl;

    tls.fR.intervalStart=-15.;//-1;//-20.;
    tls.fR.intervalEnd=arguments.evIntervalEnd;//25.;//10.;
    //tls.fR.closestRootSpacing=0.0001;
    tls.fR.rootValueAccuracy=1E-8;
    tls.fR.takenZero=1E-5;
    //tls.fR.gapFactorToNewBunch=3.; 
    //tls.fR.iterationsBeforeStepIncrementation=16;
    //tls.fR.incrementationFactor=1.9;
    //tls.fR.stepsBetweenRootsInBunch=8;  // amount of steps between roots when roots were to be equidistantly spaced
    tls.fR.intervalForMinSearch=5.;
    tls.fR.rf=propagateFunction;
    tls.fR.mf=propagateFunction;
    //tls.fR.f=invPropagateFunction;

    tls.calc();

    fout << "program finished" << endl;
}


// TODO
// check for case C2==0 what happens when C1 gets sign*1 instead 1
//
//
// check if trans(overlapMatrix)=overlapMatrix
// A^t=A
// (A^T)^*=A
// A^T=A^*
// WEC for multiband
//
//
// wannier functions at the border aren't normlized -> onsite int must to erroneous
//
//calculatedValues.clear();
//
//
//   replace wannier0Functions with wannier0Functions+wannier0FunctionsPeriod
