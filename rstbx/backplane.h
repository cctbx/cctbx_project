#ifndef RSTBX_BACKPLANE_H
#define RSTBX_BACKPLANE_H

namespace rstbx {

struct backplane {
 public:
  int boxnbg;
  double boxmean, boxvar, boxstd;
  double Sum_x,Sum_x2;
  backplane(){
  //    cout<<"base constructor"<<endl;
clear(); }
  virtual void accumulate (const int&, const int&, const int& px){
    //cout<<"base accumulate"<<endl;
    Sum_x += px;
    Sum_x2 += (double)(px)*px;
    boxnbg++;
  }
  virtual void clear(){
    boxnbg = 0;
    boxmean = boxvar = boxstd = 0.;
    Sum_x = Sum_x2 = 0.;
  }
  virtual void reinitialize(const int&, const int&){
    //cout<<"base reinitialize"<<endl;
    clear();
  }
  virtual void finish(){
    //cout<<"base finish"<<endl;
    boxmean = Sum_x/boxnbg;
    boxvar = Sum_x2/boxnbg - boxmean*boxmean;
    boxstd = std::sqrt(boxvar);
  }
  virtual inline double localmean (const int&, const int&) {
      //cout<<"base mean "<<boxmean<<endl;
    return boxmean; }

  // Trivial calculation of detector gain; re-check derivation before commenting in
  //inline double gain () const { //Suggested by Mosflm manual: digitized value = GAIN * Equiv # of photons
  //  double sample_average = Sum_x/boxnbg;
  //  double sample_variance = (1./boxnbg)*(Sum_x2 - sample_average*sample_average);
  //  return std::sqrt(sample_variance) / sample_average; // is mean_squared / mean
  //}

  virtual ~backplane(){}
};

struct corrected_backplane: public backplane {
 private:
  int Sum_p2,Sum_pq,Sum_p,Sum_q2,Sum_q;
  double Sum_xp,Sum_xq;
  int xstart,ystart;
  double a,b,c;
  std::vector<int> rho_cache;
  std::vector<int> p_cache;
  std::vector<int> q_cache;
  double rmsd;
  double p,q; //temporary values
 public:
  corrected_backplane(const int& xst, const int& yst):
    xstart(xst),ystart(yst) {
    //cout<<"corrected constructor"<<endl;
    clear();
  }
  void reinitialize(const int& xst, const int& yst){
    //cout<<"corrected reinitialize"<<endl;
    xstart = xst; ystart = yst; clear();
  }
  inline void clear(){
    backplane::clear();
    Sum_p2=0;Sum_pq=0;Sum_p=0;Sum_q2=0;Sum_q=0;Sum_xp=0;Sum_xq=0;
    rho_cache.clear();
    p_cache.clear();
    q_cache.clear();
    rmsd=0.;
  }
  inline void accumulate(const int& x, const int& y, const int& px){
    //cout<<"corrected accumulate"<<endl;
    backplane::accumulate(x,y,px);
    int p = x-xstart;
    int q = y-ystart;
    Sum_p2+=p*p;
    Sum_pq+=p*q;
    Sum_p+=p;
    Sum_q2+=q*q;
    Sum_q+=q;
    Sum_xp+=px*p;
    Sum_xq+=px*q;
    rho_cache.push_back(px);p_cache.push_back(p);q_cache.push_back(q);
  }
  inline void finish(){
    //cout<<"corrected finish"<<endl;
    scitbx::mat3<double> rossmann(Sum_p2,Sum_pq,Sum_p,
                          Sum_pq,Sum_q2,Sum_q,
                          Sum_p,Sum_q,boxnbg);
    scitbx::vec3<double> obs(Sum_xp,Sum_xq,Sum_x);
    scitbx::mat3<double> rinv = rossmann.inverse();
    //scitbx::vec3<double> abc = rossmann.inverse()*obs;
    //a=abc[0]; b= abc[1]; c=abc[2];
    a = rinv[0]*Sum_xp + rinv[1]*Sum_xq +rinv[2]*Sum_x;
    b = rinv[3]*Sum_xp + rinv[4]*Sum_xq +rinv[5]*Sum_x;
    c = rinv[6]*Sum_xp + rinv[7]*Sum_xq +rinv[8]*Sum_x;
    for (int v=0; v<boxnbg; ++v){
      double bgobs_bgplane = rho_cache[v] - a*p_cache[v] - b*q_cache[v] -c;
      rmsd +=  bgobs_bgplane*bgobs_bgplane;
    }
    rmsd = std::sqrt(rmsd/boxnbg); //box standard deviation
    boxstd = rmsd;
  }
  inline double localmean(const int&x, const int&y){
    return a*(x-xstart)+b*(y-ystart)+c;
  }
};
} //namespace
#endif// RSTBX_BACKPLANE_H
