//collection of routines from Numerical Recipes for C for efficient solution of ordinary ODE's
//here: using the Bulirsch-Stoer method, which is generally higher precision than 4th order Runge-Kutta
//modified to object-orientated structure for C++

//create object structure for passing all parameters and lists of function and variable values between different solvers and drivers


#ifndef screen_output
#define screen_output 0
#endif



#ifndef ODE_INTEGRATION_DEFINED
#define ODE_INTEGRATION_DEFINED

#include"nr/nr.h"

class ODE_parameters{
public:
	unsigned N_eqns;	///the number of (real) equations
	unsigned t_length;	///actual number of t-values and function values tabulated
	unsigned max_t_length;  ///maximum number of values of t to save evaluated solution at (arrays have to be reallocated if this becomes longer)
	unsigned max_steps;      /// the maximum number of allowed Ridder steps
	double t_end;	/// final value of time for next step of calling odeint
	double dt_save_step;  /// approximate time steps for which solution is saved
	double precision;	///precision with which to solve ODE
	double guessed_stepsize;   ///step size of first step taken
	double min_stepsize;       ///the routine is stopped and gives an error if the step size is smaller than this. I don't see the point of this, but its given in "numerical recipes"
	double max_step_size_for_dense_output;  /// if a solution y(t) is required later, this sets the upper allowed maximum spacing between two points. Make this very large if only the solution at a given end point is required such that the Ridder algorithm can automatically choose the best step size.
	unsigned ok_steps;		///keep track while solving how many good steps were taken, accumulates if odeint is executed multiple times
	unsigned bad_steps;		///keep track while solving how many bad steps were taken, accumulates if odeint is executed multiple times

    ODE_parameters(){N_eqns=1;t_length=0;max_t_length=int(1E6);max_steps=int(1E6);t_end=0;dt_save_step=1;precision=1E-7;guessed_stepsize=0.1;min_stepsize=0;max_step_size_for_dense_output=1E9;ok_steps=0;bad_steps=0;}   ///default constructor to set the values to some reasonable values, such that not all of them have to be initialized every time an object is created
	void disp(){printf("\n**********************************\nthis ODE_parameter object:\nN_eqns=%d\nt_length=%d\nmax_t_length=%d\nmax_steps=%d\nt_end=%f\ndt_save_step=%g\nprecision=%g\nguessed_stepsize=%f\nmin_stepsize=%f\nok_steps=%d\nbad_steps=%d\n*********************************\n", N_eqns,t_length,max_t_length,max_steps,t_end,dt_save_step,precision,guessed_stepsize,min_stepsize,ok_steps,bad_steps);}
};

template<class physical_parameters_type>
class ODE{
public:
	physical_parameters_type A;	///an object that contains a collection of all phyiscal parameters and arrays used for computation in the function where the derivative in the differential equation is computed
	ODE_parameters P;   /// all parameters for this ODE are combined in the object P, which allows for an easier passing of these

	///each time a new step is determined, a new N_eqns-dimensional vector is added to y and to t. Only initialized then
	double* t;	/// time (independent variable)
	double** y;	///	the set of functions to be determined (dependent variables) y_i(t) listed as y[t][i], where i runs from [0,N_eqns-1]
	void(*dy_dt_fct)(double t, double* y_in, double* dy_dt_out, physical_parameters_type &A);

	ODE(unsigned N_eqns_tmp=1,double precision_tmp=1E-6) { ///constructor with default arguments if no others are passed
        P.max_t_length=unsigned(1E6); ///the temporal length of saved function values may never exceed max_t_length. t_length counts how many values have been initialized
        P.max_steps=1E6;
        P.t_length=0;
        P.N_eqns=N_eqns_tmp;
        P.precision=precision_tmp;
        P.min_stepsize=0.0;
        P.max_step_size_for_dense_output=1E20;
        P.guessed_stepsize=0.1;
        P.ok_steps=0;
        P.bad_steps=0;
        P.dt_save_step=1.0;		//some default value - should be overwritten
        P.t_end=0.0;

        t=(double*) malloc(P.max_t_length*sizeof(double));
        y=(double**) malloc(P.max_t_length*sizeof(double*));
    }

	ODE(ODE_parameters P_in) { ///constructor to initialize object with this ODE parameter object
        P=P_in;
        t=(double*) malloc(P.max_t_length*sizeof(double));
        y=(double**) malloc(P.max_t_length*sizeof(double*));
    }

    ~ODE() {
        free(t);
        for(unsigned ct=0; ct<P.t_length;ct++) 
            free(y[ct]);
        free(y);
    }

    void set_initial_y_vec(double const*const start_vec) {	///copies the elements of the N_eqns-dimensional vector to the first index position 0 (or one ahead of the rest if previous values exist) into the internal matrix y
        y[P.t_length]=(double*) malloc(P.N_eqns*sizeof(double));
        for(unsigned ct=0;ct<P.N_eqns; ct++) { 
            y[P.t_length][ct]=start_vec[ct];
        }
        P.t_length++;
    }

    void write_sol_to_Matlab_file(const std::string &fileName) {
        FILE *F=fopen(fileName.c_str(),"w");
        fprintf(F,"tmp=[");
        for (unsigned i=0;i<P.t_length;i++) {fprintf(F,"%g %g %g\n",t[i],y[i][0],y[i][1]);}
        fprintf(F,"];\n");
        fprintf(F,"r=tmp(:,1);\n");
        fprintf(F,"u=tmp(:,2);\n");
        fprintf(F,"du_dr=tmp(:,3);\n");
        fprintf(F,"figure;\nplotbold(r,u);\nhold on;\nplotbold(r,du_dr);\nlegend({'u(r)','du(r)/dr'});\n");
        fclose(F);
    }
};



//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//imported from pzextr.cpp

extern Mat_DP *d_p;

void pzextr(const unsigned nv, const unsigned iest, const double xest, double *yest, double *yz, double *dy, double *x_p_not_global)
/**
Use polynomial extrapolation to evaluate nv functions at x = 0 by fitting a polynomial to a
sequence of estimates with progressively smaller values x = xest, and corresponding function
vectors yest[1..nv]. This call is number iest in the sequence of calls. Extrapolated function
values are output as yz[1..nv], and their estimated error is output as dy[1..nv].
*/
{	unsigned j,k1;
	double q,f2,f1,delta;
	double *c=(double*) malloc(nv*sizeof(double));//new double[nv];
	double *x=x_p_not_global;
	Mat_DP &d=*d_p;
	x[iest]=xest; /// Save current independent variable.
	for (j=0;j<nv;j++) 
        dy[j]=yz[j]=yest[j];
	if (iest == 0) { /// Store first estimate in first column.
		for (j=0;j<nv;j++) 
            d[j][0]=yest[j];
	} else {
		for (j=0;j<nv;j++) 
            c[j]=yest[j];
		for (k1=0;k1<iest;k1++) {
			delta=1.0/(x[iest-k1-1]-xest);
			f1=xest*delta;
			f2=x[iest-k1-1]*delta;
			for (j=0;j<nv;j++) { /// Propagate tableau 1 diagonal more.
				q=d[j][k1];
				d[j][k1]=dy[j];
				delta=c[j]-q;
				dy[j]=f1*delta;
				c[j]=f2*delta;
				yz[j] += dy[j];
			}
		}
		for (j=0;j<nv;j++) 
            d[j][iest]=dy[j];
	}
	free(c);
}

template<class physical_parameters_type>
void mmid(unsigned nvar, double *y, double *dydx, const double xs, const double htot, const unsigned nstep, double *yout,	void derivs(const double, double*, double*, physical_parameters_type &), physical_parameters_type &A)
{	unsigned i,n;
	double x,swap,h2,h;
	double *ym=(double*) malloc(nvar*sizeof(double));
	double *yn=(double*) malloc(nvar*sizeof(double));
	h=htot/nstep;
	for (i=0;i<nvar;i++) {
		ym[i]=y[i];
		yn[i]=y[i]+h*dydx[i];
	}
	x=xs+h;
	derivs(x,yn,yout, A);
	h2=2.0*h;
	for (n=1;n<nstep;n++) {
		for (i=0;i<nvar;i++) {
			swap=ym[i]+h2*yout[i];
			ym[i]=yn[i];
			yn[i]=swap;
		}
		x += h;
		derivs(x,yn,yout,A);
	}
	for (i=0;i<nvar;i++)
	yout[i]=0.5*(ym[i]+yn[i]+h*yout[i]);
	free(ym);
	free(yn);
}



// imported from bsstep.cpp
Mat_DP *d_p;

template<class physical_parameters_type>
void bulirsch_stoer_step(const unsigned nv, double *y, double *dydx, double &xx, const double htry,
                         const double eps, double *yscal, double &hdid, double &hnext,
                         void derivs(const double, double*, double*, physical_parameters_type &), physical_parameters_type &A)
/** Bulirsch-Stoer step with monitoring of local truncation error to
ensure accuracy and adjust stepsize. Input are the dependent variable
vector y[1..nv] and its derivative dydx[1..nv] at the starting value of
the independent variable xx. Also input are the stepsize to be attempted
htry, the required accuracy eps, and the vector yscal[1..nv] against
which the error is scaled. On output, y and xx are replaced by their new
values, hdid is the stepsize that was actually accomplished, and hnext
is the estimated next stepsize. derivs is the user-supplied routine that
computes the right-hand side derivatives. Be sure to set htry on
successive steps to the value of hnext returned from the previous step, as is the case if the routine is called
by odeint.*/
{
	const unsigned KMAXX=8;  /// Maximum row number used in the extrapolapolation
	const unsigned IMAXX=(KMAXX+1);
	const double SAFE1=0.25; /// Safety factors - set in numerical recipes
	const double SAFE2=0.7;
	const double REDMAX=1.0e-5; /// Maximum factor for stepsize reduction.
	const double REDMIN=0.7; /// Minimum factor for stepsize reduction.
	const double TINY=1.0e-30;/// Prevents division by zero.
	const double SCALMX=0.1; /// 1/SCALMX is the maximum factor by which a stepsize can be increased.
	static const int nseq_d[IMAXX]={2,4,6,8,10,12,14,16,18};
	static unsigned first=1,kmax,kopt;
	static double epsold = -1.0,xnew;

	//double *a=new double[IMAXX];
	static Vec_DP a(IMAXX);
	static Mat_DP alf(KMAXX,KMAXX);
	bool exitflag=false;
	unsigned i,iq,k,kk,km,reduct;
	double eps1,errmax,fact,h,red,scale,work,wrkmin,xest;
	red=-1E200; //initialize to some value that the compiler does not give a warning
	scale=-1E200;
	Vec_INT nseq(nseq_d,IMAXX);
	Vec_DP err(KMAXX);
	double *yerr=(double*) malloc(nv*sizeof(double));
	double *ysav=(double*) malloc(nv*sizeof(double));
	double *yseq=(double*) malloc(nv*sizeof(double));
	double *x_p_not_global=(double*) malloc(KMAXX*sizeof(double));

	d_p=new Mat_DP(nv,KMAXX);
	if (eps != epsold) {   ///A new tolerance, so reinitialize.
		hnext = xnew = -1.0e29; ///“Impossible” values.
		eps1=SAFE1*eps;
		a[0]=nseq[0]+1;
		for (k=0;k<KMAXX;k++) 
            a[k+1]=a[k]+nseq[k+1];
		for (iq=1;iq<KMAXX;iq++) { ///Compute a(k, q).
			for (k=0;k<iq;k++)
			alf[k][iq]=pow(eps1,(a[k+1]-a[iq+1])/
			((a[iq+1]-a[0]+1.0)*(2*k+3)));
		}
		epsold=eps;
		for (kopt=1;kopt<KMAXX-1;kopt++) ///Determine optimal row number for convergence.
		if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
		kmax=kopt;
	}
	h=htry;
	for (i=0;i<nv;i++) 
        ysav[i]=y[i];  /// Save the starting values.
	if (xx != xnew || h != hnext) { /// A new stepsize or a new integration: re-establish the order window.
		first=1;
		kopt=kmax;
	}
	reduct=0;
	for (;;) {
		for (k=0;k<=kmax;k++) {  /// Evaluate the sequence of modified midpoint integrations
			xnew=xx+h;
			if (xnew == xx)
			{
				printf("step size underflow in bulirsch_stoer_step");
				exit(0);
			}
			mmid(nv,ysav,dydx,xx,h,nseq[k],yseq,derivs,A);
			xest=SQR(h/nseq[k]);
			pzextr(nv,k,xest,yseq,y,yerr,x_p_not_global);  /// Perform extrapolation.
			if (k != 0) { /// Compute normalized error estimate epsilon(k).
				errmax=TINY;
				for (i=0;i<nv;i++) 
                    errmax=MAX(errmax,fabs(yerr[i]/yscal[i]));
				errmax /= eps; /// Scale error relative to tolerance.
				km=k-1;
				err[km]=pow(errmax/SAFE1,1.0/(2*km+3));
			}
			if (k != 0 && (k >= kopt-1 || first)) { /// In order window.
				if (errmax < 1.0) { ///Converged.
					exitflag=true;
					break;
				}
				if (k == kmax || k == kopt+1) {  /// Check for possible stepsize reduction.
					red=SAFE2/err[km];
					break;
				}
				else if (k == kopt && alf[kopt-1][kopt] < err[km]) {
					red=1.0/err[km];
					break;
				}
				else if (kopt == kmax && alf[km][kmax-1] < err[km]) {
					red=alf[km][kmax-1]*SAFE2/err[km];
					break;
				}
				else if (alf[km][kopt] < err[km]) {
					red=alf[km][kopt-1]/err[km];
					break;
				}
			}
		}
		if (exitflag) break;
		red=MIN(red,REDMIN); /// Reduce stepsize by at least REDMIN
		red=MAX(red,REDMAX); /// and at most REDMAX.
		h *= red;
		reduct=1;
	} /// Try again.
	xx=xnew; /// Successful step taken.
	hdid=h;
	first=0;
	wrkmin=1.0e35; /// Compute optimal row for convergence
	for (kk=0;kk<=km;kk++) { /// and corresponding stepsize.
		fact=MAX(err[kk],SCALMX);
		work=fact*a[kk+1];
		if (work < wrkmin) {
			scale=fact;
			wrkmin=work;
			kopt=kk+1;
		}
	}
	hnext=h/scale;
	if (kopt >= k && kopt != kmax && !reduct) {  /// Check for possible order increase, but not if stepsize was just reduced.
		fact=MAX(scale/alf[kopt-1][kopt],SCALMX);
		if (a[kopt+1]*fact <= wrkmin) {
			hnext=h/fact;
			kopt++;
		}
	}

	delete d_p;
	free(x_p_not_global);
	free(yerr);
	free(ysav);
	free(yseq);
}



extern unsigned kmax;

template<class physical_parameters_type>
void odeint(ODE<physical_parameters_type> *D)
/**	Runge-Kutta driver with adaptive stepsize control. Integrate starting values D->y[t_length-1][1..nvar]
from x1 to x2 with accuracy eps, storing intermediate results in global variables. h1 should
be set as a guessed first stepsize, P.min_stepsize as the minimum allowed stepsize (can be zero). On
output nok and nbad are the number of good and bad (but retried and fixed) steps taken, and
ystart is replaced by values at the end of the integration interval. derivs is the user-supplied
routine for calculating the right-hand side derivative, while bulirsch_stoer_step is the name of the stepper
routine to be used. */
{
	unsigned start_index=D->P.t_length-1;
	double x1=D->t[start_index];
	double x2=D->P.t_end;
	double dxsav=D->P.dt_save_step;
	if(screen_output) {
        printf("start_t=%g    end_t=%g\n",x1,x2);
    }
	unsigned nvar=D->P.N_eqns;
	double h1=D->P.guessed_stepsize;
	double eps=D->P.precision;
	const unsigned MAXSTP=D->P.max_steps;
	const double TINY=1.0e-30;
	unsigned i,nstp;
	double tsav,t_here,hnext,hdid,h;
	tsav=-1E200;  //initialize to some value that the compiler does not give a warning
	double *yscal=(double*) malloc(nvar*sizeof(double));
	double *y_here=(double*) malloc(nvar*sizeof(double));
	double *dydx=(double*) malloc(nvar*sizeof(double));
	t_here=x1;
    h1=min(h1,D->P.max_step_size_for_dense_output);
	h=SIGN(h1,x2-x1);
	unsigned kount = 0;
	for (i=0;i<nvar;i++) 
        y_here[i]=D->y[start_index][i];
	if (kmax > 0) 
        tsav=t_here-dxsav*2.0;   /// Assures storage of first step.
	for (nstp=0;nstp<MAXSTP;nstp++) {  ///Take at most MAXSTP steps.
		D->dy_dt_fct(t_here,y_here,dydx,D->A);  ///here the value of dydx is written for this specific t_here and y_here
		for (i=0;i<nvar;i++) /// Scaling used to monitor accuracy. This general-purpose choice can be modified if need be.
            yscal[i]=fabs(y_here[i])+fabs(dydx[i]*h)+TINY;
		if (kmax > 0 && kount < kmax-1 && fabs(t_here-tsav) > fabs(dxsav)) {
			tsav=t_here;
			/// Store intermediate results.
			if(t_here>D->t[start_index])
			{
				D->t[D->P.t_length]=t_here;
				D->y[D->P.t_length]=(double*) malloc(D->P.N_eqns*sizeof(double));
				if(screen_output) {
                    cout<<"y created at t_length="<<D->P.t_length<<"  t_here="<<t_here<<endl;
                }
				for (i=0;i<D->P.N_eqns;i++) 
                    D->y[D->P.t_length][i]=y_here[i];
				D->P.t_length++;
			}
		}

		if ((t_here+h-x2)*(t_here+h-x1) > 0.0) {
            h=x2-t_here;   /// If stepsize can overshoot, decrease.
		}
		bulirsch_stoer_step(nvar,y_here,dydx,t_here,h,eps,yscal,hdid,hnext,D->dy_dt_fct,D->A);
		if (hdid == h) 
            ++(D->P.ok_steps); 
        else 
            ++(D->P.bad_steps);
		if ((t_here-x2)*(x2-x1) >= 0.0) {  ///Are we done?
			if (kmax != 0) {
				///Save final step.
				D->t[D->P.t_length]=t_here;
				D->y[D->P.t_length]=(double*) malloc(D->P.N_eqns*sizeof(double));
				if(screen_output) {
                    cout<<"y created at t_length="<<D->P.t_length<<endl;
                }
				for (i=0;i<D->P.N_eqns;i++) 
                    D->y[D->P.t_length][i]=y_here[i];
				D->P.t_length++;
			}
			free(yscal);
			free(y_here);
			free(dydx);
			return;  /// Normal exit.
		}
		if (fabs(hnext) <= D->P.min_stepsize) {
			printf("Step size too small in odeint");
			exit(0);
		}
		hnext=min(hnext, D->P.max_step_size_for_dense_output);
		h=hnext;

	}
	printf("Too many steps in routine odeint\n\n");
	exit(0);
}


#endif // ODE_INTEGRATION_DEFINED
