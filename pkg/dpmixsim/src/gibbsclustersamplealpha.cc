//$Id: gibbsclustersamplealpha.cc,v 1.1 2010/04/19 17:22:17 afs Exp $
//
// Updating alpha: a0, b0, prior for alpha if alpha is updated
//
// modifications from 1st version: MX1 and MX2
//
//   

#define RECORD

#include <Rmath.h>
#include <iostream>
#include <map>
#include <iterator>
#include <set>
#include <valarray>
#include <numeric>
#include <unistd.h>
#include <stdint.h>
// #include <sys/times.h>
#include <assert.h>
#include <stdio.h>

extern "C" {

using namespace std;
//
//
typedef map<int, set<int> >::const_iterator map_citr ;
typedef set<int>::const_iterator set_citr ;

extern void testprint(map<int, set<int> >& table);

//------------------------------------------------------------------------------

void gibbsdpm(const double* xin, const int* pn,
		const double* M, const double* a, const double* b,
		const int* upalpha, const double* a0, const double* b0,
		const int* maxiter, int *recno, const int *kmax,
		int* krec, double* w, double* phirec, double* sigrec)
{
	//--------------------- step 0: initialization ---------------------//
	// start with x. Assume no repeated value!
	double alpha = *M;
	const int n											= *pn;
	const double shape							= *a + n/2.0; // shape parameter 
	const valarray<double> x(xin,n);
	const valarray<double> xsq 			= pow(x,2.0);
	int* nj = new int[n];
	fill(nj,nj+n,1);
	int* s = new int[n];
	for(int i=0; i < n; i++) 	s[i]	= i;
	double sig											= 1.0;				// initialize to 1
	//
	// record indices
	if (*recno > *maxiter || *recno < 0 ) *recno	= *maxiter;
	const int startrec  													=  *maxiter - *recno;
	// initialize table of n clusters with one element each
  map<int, double> theta;
  map<int, set<int> > stb;
  for (int i=0; i < n; i++) {
		set<int> ms; 
		ms.insert(i);
		stb.insert	(pair<int, set<int> >(i, ms)); 
		theta.insert(pair<int, double   >(i, xin[i])); 
	}
	//
	//============================ iterate... ============================== //
	//
	cout << "Processing " << *maxiter << " iterations:\t '.' = 100 iterations" << endl; 
	//
	for (int iter=0; iter < *maxiter; iter++) {
		const double sig2						= sig*2.0;
		const double sigsqrt				= sqrt(sig);
		const double sigplus				= sig+1.0;
		const double sigplussqrt		= sqrt(sigplus);
		const double sigplus2				=	sigplus*2.0;
		// ------------------- step 1: s_i | theta_(-i), sig, x ---------------//
		// calculate probability
		// reducing the number of weight calculations to the number of clusters
		for(int i=0; i < n; i++) {
			const int clsizep					= stb.size()+1; // accounting for new cluster
			valarray<double> probs(0.0,clsizep);
			int kp 										= 0;
		  for (map_citr it=stb.begin(); it != stb.end(); ++it, kp++) {
				const int key 					= (*it).first;
				const int mdim 					= (*it).second.size();
				// begin MX1 
				double mw=1.0;
				double xth = x[i]-theta[key];
				if(abs(xth) > 1.0e-7) // evaluate only if difference is not negligible, 1.0 otherwise
					mw = exp(-pow(xth,2)/sig2)/sigsqrt;
				// end MX1 
				// previous code: 
				// const double mw					= exp(-(pow((x[i]-theta[key]),2))/sig2)/sigsqrt;
				// 
				if(key != s[i]) 
					probs[kp]							= mdim*mw;
				else
					probs[kp]							= (mdim-1)*mw;
			}
			const double bnew 				= alpha*exp(-xsq[i]/ sigplus2) / sigplussqrt;
			const double bden					= probs.sum() + bnew;
			probs[clsizep-1]					= bnew;
			probs 								 		/= bden;
			// sample probabilities: 1st sampling option
 			partial_sum(&probs[0], &probs[0]+clsizep, &probs[0]);
			const double* rprob				= find_if(&probs[0], &probs[0]+clsizep,
																	bind2nd(greater<double>(), runif(0,1)));
			const int cl							= rprob - (&probs[0]);
			int newcl, choice;
			if(cl < (clsizep-1)) {
			 	newcl 									= 0;
			 	map_citr it							= stb.begin();
				// for(int k=0; k < cl ; k++) it++;
				advance(it, cl); 
			 	const int key		 				= (*it).first;
				//?? random schoice within cluster really required ?
				// begin MX2 
				set_citr itset				= stb[key].begin();
				// end MX2 
				//
				/*
				 * begin previous code:   
				const int mdim 					= (*it).second.size();
				const int schoice				= int(runif(0,1) * mdim);
				// set_citr itset				= stb[key].begin();
				// for(int k=0; k < schoice ; k++) itset++; // time bottleneck for large datasets !!!
				set_citr itset;
				if(schoice < (mdim >> 1)) { // ! a bit better running time
					itset 								= stb[key].begin();
					advance(itset, schoice); 
				}
				else {
					itset = stb[key].end(); // one past last element
					advance(itset, -(mdim-schoice)); // schoice from 0:(mdim-1)

				}
				* end previous code
				*/
				choice 									=  *itset;

			}
			else {
				newcl										= 1;
				choice									= -1; // don't care
			}
			// ------------------------------------------------
			int currentkey;
			if (!newcl) {		  						 // put s[i] into cluster s[choice]
				if (s[i] != s[choice]) {
			 		const int siold 			= s[i];
					// DEATH: decrease cluster size 
					nj[s[i]]--;
					if (nj[s[i]] < 0){	
						cout << "error: negative cluster size \n"; exit(1); }
					// update stb and stb.size() iteratively
					if (!nj[s[i]]) {
						stb.erase(s[i]);
						theta.erase(s[i]);
					}
					nj[s[choice]]++;
					s[i]									= s[choice];
					// remove old element from previous position in set
					map_citr it 					= stb.find(siold);
					if(it != stb.end()) stb[siold].erase(i);
					// insert in appropriate set of clusters;  first nonzero place
					stb[s[choice]].insert(i);
					currentkey						= s[choice];
				}
				else 
					continue;
					// currentkey						= s[choice];
			} else {												// new cluster from G=N(0, 1)
				const int siold 				= s[i];
				// DEATH: decrease cluster size 
				nj[s[i]]--;
				if (nj[s[i]] < 0){	
					cout << "error: negative cluster size \n"; exit(1); }
				// update stb and stb.size() iteratively
				if (!nj[s[i]]) {
					stb.erase(s[i]);
					theta.erase(s[i]);
				}
				// BIRTH of cluster : we need to update table(s) 
				const int newcluster 		= find(nj, nj+n, 0) - nj;
				s[i]										= newcluster;
				nj[newcluster]					= 1;
				// insert a new cluster 
				set<int> ms; 
				ms.insert(i);
				stb.insert(pair<int, set<int> >(newcluster, ms));  
				// remove old element from previous position in set
				map_citr it 						= stb.find(siold);
				if(it != stb.end()) stb[siold].erase(i);
				// !! the previous statement may have erased the inserted element 
				// !! find another solution for this double insertion
				stb[newcluster].insert(i);
				currentkey							= newcluster;
			}
		}
		// ------------------- step 2: phi | s, sig, x ----------------------------//
		// ~ N(sum_(x_i=j) x_i/nj(1+sig), sig/nj(1+sig))
		// sj = sum x_i, s_i = j..... for cluster 1,...,n (with empty clusters)
		for (map_citr it=stb.begin(); it != stb.end(); ++it) {
			const int key 						= (*it).first;
			const int mdim 						= (*it).second.size();
			set_citr t1								= stb[key].begin();
			set_citr t2								= stb[key].end();
			double sjkey							= 0.0;
			for(; t1 != t2; t1++)
				sjkey										+= x[*t1];
			const double den 					= (mdim*sigplus);
			// theta[currentkey]			= rnorm(sjkey/den, sig/den);
			theta[key]								= rnorm(sjkey/den, sqrt(sig/den));
		}
		// ------------------- step 3: sig | s, theta, x ---------------------------//
		// ! This is a costly operation, so we do it once each iteration only
		// ! Alternative: use one sig per cluster
		// ! rgamma.c uses shape and scale
		// sig  =  1/rgamma(1, a+n/2, b+ sum((x-theta)^2)/2 ) // original code
		// rate =  b + sum((x-theta)^2)/2 
		/*
		const double sumx 					= (pow(x-theta, 2.0)).sum();
		const double scale					= 1.0 / (*b + sumx/2.0);
		sig 												= 1.0 / rgamma(shape, scale); 
		*/
		double sumx 								= 0.0;
		for (map_citr it=stb.begin(); it != stb.end(); ++it) {
			const int key 						= (*it).first;
			const int mdim 						= (*it).second.size();
			set_citr t1								= stb[key].begin();
			set_citr t2								= stb[key].end();
			valarray<size_t>index(mdim);
			copy(t1, t2, &index[0]);
			valarray<double>aux 			= x[index];
			sumx 											+= (pow(aux-theta[key], 2.0)).sum();
		}
		const double scale					= 1.0 / (*b + sumx/2.0); 
		sig 												= 1.0 / rgamma(shape, scale); 
		// ------------------------ updating alpha ----------------------------------//
	  if(*upalpha) {
			double dtemp, dtemp1, dtemp2;
			const int ncl 							= stb.size();
	    dtemp1											= alpha + 1.0;
	    dtemp2											= (double) n;
	    dtemp 											= *b0-log(rbeta(dtemp1, dtemp2));
	    dtemp1											= (*a0+ncl-1.0) / (n*dtemp);
	    if(unif_rand() < dtemp1) {
	      dtemp2										= *a0+ncl;
	      alpha											= rgamma(dtemp2, 1.0/dtemp);
	    }
	    else {
	      dtemp2										= *a0+ncl-1.0;
	      alpha											= rgamma(dtemp2, 1.0/dtemp);
	    }
	  }
		// ------------------------ record -----------------------------------------//
#ifdef RECORD
		if(iter >= startrec){ 
			const int k 							= iter-startrec;
			const int ncl 						= stb.size(); // counting the number of clusters
			assert(ncl <= *kmax);
			krec[k] 									= ncl;
			int kn=k*(*kmax);
			for (map_citr it=stb.begin(); it != stb.end(); ++it, ++kn) {
				phirec[kn]	= theta[(*it).first];
				w[kn] 			= double((*it).second.size())/n; // weights
				// w[kn] 			= double((*it).second.size()); // weights
			}
			sigrec[k]									= sig;
		}
		if(iter%100 == 0)
			cout.flush() << ".";
#endif
		// -----------------------------------------------------------------//
	}
	delete[] s;
	delete[] nj;
}

} // extern "C"

