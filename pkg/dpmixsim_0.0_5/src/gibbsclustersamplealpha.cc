//$Id: gibbsclustersamplealpha.cc,v 1.1 2010/07/25 17:29:47 afs Exp afs $
//
// Updating alpha: a0, b0, prior for alpha if alpha is updated
//
// modifications from 1st version: MX1 and MX2
//

// #define DEBUG 
#define RECORD

#include <Rmath.h>
#include <iostream>
#include <map>
#include <vector>
#include <iterator>
#include <set>
#include <valarray>
#include <numeric>
#include <unistd.h>
#include <stdint.h>
// #include <sys/times.h>
#include <assert.h>
#include <stdio.h>

using namespace std;

typedef map<int, set<int> >::const_iterator map_citr ;
typedef set<int>::const_iterator set_citr ;
typedef map<int, double >::const_iterator map_ctab ;
typedef map<int, double >::const_reverse_iterator map_crtab ;


extern "C" {

void gibbsdpm(const double* xin, const int* ps, int* njinit, int* pnjlen, const int* pn,
    const double* M, const double* a, const double* b,
    const double* a0, const double* b0, const double* minvar,
    const int* upalpha, const int* maxiter, int *recno, const int *kmax,
    int* krec, double* w, double* phirec, double* varrec);

extern void testprint(map<int, set<int> >& table);

}

//--------------------------------------------------------------------//

void gibbsdpm(const double* xin, const int* ps, int* njinit, int* pnjlen, const int* pn,
    const double* M, const double* a, const double* b,
    const double* a0, const double* b0, const double* minvar,
    const int* upalpha, const int* maxiter, int *recno, const int *kmax,
    int* krec, double* w, double* phirec, double* varrec)
{
  //--------------------- step 0: initialization ---------------------//
  double alpha                     = *M;
  const int n                      = *pn;
  const int njlen                  = *pnjlen;
  const double shape               = *a + n/2.0; // shape parameter 
  const valarray<double> x(xin,n);
  const valarray<double> xsq       = pow(x,2.0);
  // record indices
  if (*recno > *maxiter || *recno < 0 ) *recno  = *maxiter;
  const int startrec               =  *maxiter - *recno;
  //
  // Customized initialization
  vector<int> s(ps,ps+n);
  vector<int> nj(n,0);
  for(int i=0; i < njlen; ++i)
  copy(njinit, njinit+njlen, nj.begin());
  map<int, double> theta;
  map<int, double> sigmap;
  map<int, set<int> > stb;
  for (int i=0; i < njlen; i++) {
    set<int> ms; 
    stb.insert(pair<int, set<int> >(i, ms)); 
  }
  for (int i=0; i < n; i++) {
    stb[s[i]].insert(i);
    theta.insert(pair<int, double>(i, xin[i])); 
    sigmap.insert(pair<int, double>(i, 0.5)); 
  }
#ifdef DEBUG
  testprint(stb);
#endif
  //
  //============================ iterate... ============================//
  //
  cout << "Processing " << *maxiter << " iterations:\t '.' = 100 iterations" << endl; 
  //
  for (register int iter=0; iter < *maxiter; ++iter) {
    //
    // ----------step 3: sig | s, theta, x ----- AND --- step 2: phi | s, sig, x  ----//
    // Sampling step 2 and step 3 first
    // ! using one sig per cluster
    // ! rgamma.c uses shape and scale
    // sig  =  1/rgamma(1, a+n/2, b+ sum((x-theta)^2)/2 ) 
    // rate =  b + sum((x-theta)^2)/2 
    // ~ N(sum_(x_i=j) x_i/nj(1+sig), sig/nj(1+sig))
    for (map_citr it=stb.begin(); it != stb.end(); ++it) {
      const int key              = (*it).first;
      const int mdim             = (*it).second.size();
      set_citr t1                = stb[key].begin();
      set_citr t2                = stb[key].end();
      valarray<size_t>index(mdim);
      copy(t1, t2, &index[0]);
      valarray<double>aux        = x[index];
      const double sumx          = (pow(aux-theta[key], 2.0)).sum();
      const double scale         = 1.0 / (*b + sumx/2.0); 
      double sigmaval            = 1.0 / rgamma(shape, scale); 
      if(sigmaval < *minvar) sigmaval = *minvar;
      sigmap[key]                = sigmaval;
			const double sjkey         = aux.sum();
      const double den           = mdim*(sigmaval+1.0);
      theta[key]                 = rnorm(sjkey/den, sqrt(sigmaval/den));
    }
    //
    // ------------------- step 1: s_i | theta_(-i), sig, x -------------//
    // calculate probability
    // reducing the number of weight calculations to the number of clusters
    //---------------------
    // initialize valarrays 
    size_t clsize               = stb.size();
    valarray<double> thetaval(clsize);
    valarray<double> sigval(clsize);
    valarray<double> mdimval(clsize); // ! need double for valarray ops
    map_citr it                 = stb.begin();
    map_ctab ith                = theta.begin();
    map_ctab its                = sigmap.begin();
    for (size_t j=0; j < stb.size(); ++j, ++ith, ++its, ++it) {
      thetaval[j]               = (*ith).second;
      sigval[j]                 = (*its).second;
      mdimval[j]                = (*it).second.size();
    }
    //---------------------
    for(register int i=0; i < n; ++i) {
      if(clsize != stb.size()) { // resize and update
        clsize                  = stb.size();
        thetaval.resize(clsize);  
        sigval.resize(clsize);  
        mdimval.resize(clsize);  
        map_citr it             = stb.begin();
        map_ctab ith            = theta.begin();
        map_ctab its            = sigmap.begin();
        for (size_t j=0; j < stb.size(); ++j, ++ith, ++its, ++it) {
          thetaval[j]           = (*ith).second;
          sigval[j]             = (*its).second;
          mdimval[j]            = (*it).second.size();
        }
      }
      // else : use updated mdimval
      //------------------
      valarray<double> prob(clsize);
      valarray<double> mw(clsize);
      mw                        = exp(-(pow((x[i]-thetaval),2))/(sigval * 2.0))/(sqrt(sigval));
      prob                      = mdimval * mw;
      map_citr it0              = stb.begin();
      map_citr itsi             = stb.find(s[i]);
      const int kx              = distance(it0, itsi);
      prob[kx]                  = (mdimval[kx] - 1) * mw[kx]; 
      const double sigplus      = sigval.max()+1.0;
      const double bnew         = alpha*exp(-xsq[i]/(sigplus*2.0))/(sqrt(sigplus));
      const double bden         = prob.sum() + bnew;
      const int clsizep         = stb.size()+1;
      valarray<double> probs(clsizep);
      probs                     = prob[slice(0,clsize,1)];
      probs[clsizep-1]          = bnew;
      probs                    /= bden;
      // sample probabilities: 1st sampling option
      partial_sum(&probs[0], &probs[0]+clsizep, &probs[0]);
      const double* rprob       = find_if(&probs[0], &probs[0]+clsizep,
                                  bind2nd(greater<double>(), runif(0,1)));
      const int cl              = rprob - (&probs[0]);
      //-----------------------------------------------
      int newcl, choice;
      if(cl < (clsizep-1)) {
        newcl                   = 0;
        map_citr it             = stb.begin();
        advance(it, cl); 
        const int key           = (*it).first;
        // begin MX2 //?? random schoice within cluster really required ?
        set_citr itset           = stb[key].begin();
        // end MX2 
        choice                   = *itset;
      }
      else {
        newcl                    = 1;
        choice                   = -1; // don't care
      }
      // ------------------------------------------------
      int currentkey;
      if (!newcl) {                   // put s[i] into cluster s[choice]
        if (s[i] != s[choice]) {
           const int siold       = s[i];
          // DEATH: decrease cluster size 
          --nj[s[i]];
          if (nj[s[i]] < 0){  
            cout << "error: negative cluster size \n"; exit(1); }
          // update stb and stb.size() iteratively
          if (!nj[s[i]]) {
            stb.erase(s[i]);
            theta.erase(s[i]);
            sigmap.erase(s[i]);
          }
          ++nj[s[choice]];
          s[i]                   = s[choice];
          /*
          // remove old element from previous position in set
          map_citr it            = stb.find(siold);
          if(it != stb.end()) stb[siold].erase(i);
          // insert in appropriate set of clusters;  first nonzero place
          stb[s[choice]].insert(i);
          currentkey             = s[choice];
          */
          // remove old element from previous position in set
          // and  update mdimval for old and new position
          map_citr it            = stb.begin();
          map_citr itold         = stb.find(siold);
          if(itold != stb.end()) {
            stb[siold].erase(i);
            size_t ixup          = distance(it, itold); 
            mdimval[ixup]        = (*itold).second.size(); // just update mdimval
          }
          // insert in appropriate set of clusters;  first nonzero place
          stb[s[choice]].insert(i);
          currentkey             = s[choice];
          map_citr itchoice      = stb.find(s[choice]);
          size_t ixup            = distance(it, itchoice);
          mdimval[ixup]          = (*itchoice).second.size(); // just update mdimval
        }
        else 
          continue;
          // currentkey          = s[choice];
      } else {                        // new cluster from G=N(0, 1)
        const int siold          = s[i];
        // DEATH: decrease cluster size 
        --nj[siold];
        if (nj[siold] < 0){ 
          cout << "error: negative cluster size \n"; exit(1); }
        // update stb and stb.size() iteratively
        if (!nj[siold]) {
          stb.erase(siold);
          theta.erase(siold);
          sigmap.erase(siold);
        }
        // BIRTH of cluster : we need to update table(s) 
        // const int newcluster  = find(nj, nj+n, 0) - nj;
        const int newcluster     = find(nj.begin(), nj.end(), 0) - nj.begin();
        s[i]                     = newcluster;
        nj[newcluster]           = 1;
        // insert a new cluster 
        set<int> ms; 
        ms.insert(i);
        stb.insert(pair<int, set<int> >(newcluster, ms));  
        // remove old element from previous position in set
        map_citr it              = stb.find(siold);
        if(it != stb.end()) stb[siold].erase(i);
        // !! the previous statement may have erased the inserted element 
        // !! find another solution for this double insertion
        stb[newcluster].insert(i);
        currentkey               = newcluster;
        // theta[currentkey]     = rnorm(x[i]/sigplus, sqrt(sig/sigplus));
        theta[currentkey]        = rnorm(x[i]/2.0, sqrt(0.5));
        sigmap[currentkey]       = 1.0;
        // map_crtab ti          = sigmap.rbegin();
        // sigmap[currentkey]    = (*ti).second; // use max cluster value
      }
    }
    // ------------------------ updating alpha --------------------------//
    if(*upalpha) {
      const int ncl              = stb.size();
      double dtemp1              = alpha + 1.0;
      double dtemp2              = (double) n;
      const double dtemp         = *b0-log(rbeta(dtemp1, dtemp2));
      dtemp1                     = (*a0+ncl-1.0) / (n*dtemp);
      if(unif_rand() < dtemp1) {
        alpha                    = rgamma(*a0+ncl, 1.0/dtemp);
      }
      else {
        alpha                    = rgamma(*a0+ncl-1.0, 1.0/dtemp);
      }
    }
    // ------------------------ record -----------------------------------------//
#ifdef RECORD
    if(iter >= startrec){ 
      const int k                = iter-startrec;
      const int ncl              = stb.size(); // counting the number of clusters
      assert(ncl <= *kmax);
      krec[k]                    = ncl;
      int kn=k*(*kmax);
      for (map_citr it=stb.begin(); it != stb.end(); ++it, ++kn) {
        phirec[kn]               = theta[(*it).first];
        w[kn]                    = double((*it).second.size())/n; // weights
        // w[kn]                 = double((*it).second.size()); // weights
        varrec[kn]               = sigmap[(*it).first];
      }
    }
    if(iter%100 == 0)
      cout.flush() << ".";
#endif
    // -----------------------------------------------------------------//
  }
}

//======================================================================//
void testprint(map<int, set<int> >& table)
{
  cout << "Table stb: " << endl; 
  map_citr it   = table.begin();
  for (; it != table.end(); ++it) {
    int key     =  (*it).first;
    set<int> S  = (*it).second;
    // cout << "size " << S.size() << "\tkey " << key << ":\t ";
    cout << "key " << key << ":\t ";
    copy(S.begin(), S.end(), ostream_iterator<int>(cout, " "));
    cout << endl;
  }
}
//======================================================================//

