#include <scitbx/lbfgsb/raw.h>
#include <scitbx/array_family/shared.h>
#include <iostream>

namespace {

  using scitbx::af::shared;
  using scitbx::fn::absolute;
  using scitbx::fn::pow2;
  using namespace scitbx::lbfgsb::raw;

  template <typename ElementType>
  ref1<ElementType>
  make_ref1(shared<ElementType>& a)
  {
    return ref1<ElementType>(a.begin(), a.size());
  }

  void
  driver3()
  {
    std::string task, csave;
    shared<bool> lsave_(4);
    ref1<bool> lsave = make_ref1(lsave_);
    int n=1000;
    int m=10;
    int iprint;
    shared<int> nbd_(n);
    ref1<int> nbd = make_ref1(nbd_);
    shared<int> iwa_(3*n);
    ref1<int> iwa = make_ref1(iwa_);
    shared<int> isave_(44);
    ref1<int> isave = make_ref1(isave_);
    double f, factr, pgtol;
    shared<double> x_(n);
    ref1<double> x = make_ref1(x_);
    shared<double> l_(n);
    ref1<double> l = make_ref1(l_);
    shared<double> u_(n);
    ref1<double> u = make_ref1(u_);
    shared<double> g_(n);
    ref1<double> g = make_ref1(g_);
    shared<double> dsave_(29);
    ref1<double> dsave = make_ref1(dsave_);
    shared<double> wa_(2*m*n+4*n+12*m*m+12*m);
    ref1<double> wa = make_ref1(wa_);
    double t1, t2, time1, time2, tlimit;
    tlimit = 0.2;
    iprint = -1;
    factr=0.0e0;
    pgtol=0.0e0;
    for(int i=1;i<=n;i+=2) {
      nbd(i)=2;
      l(i)=1.0e0;
      u(i)=1.0e2;
    }
    for(int i=2;i<=n;i+=2) {
      nbd(i)=2;
      l(i)=-1.0e2;
      u(i)=1.0e2;
    }
    for(int i=1;i<=n;i++) {
      x(i)=3.0e0;
    }
    std::cout << std::endl;
    std::cout << "     Solving sample problem." << std::endl;
    std::cout << "      (f = 0.0 at the optimal solution.)" << std::endl;
    std::cout << std::endl;
    task = "START";
    timer(time1);
    lbl_111:
    setulb(
      n,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint,
      csave,lsave,isave,dsave);
    if (task.substr(0,2) == "FG") {
      timer(time2);
      if (time2-time1 > tlimit) {
        task="STOP: CPU EXCEEDING THE TIME LIMIT.";
//      write (6,*) task;
//      int j = 3*n+2*m*n+12*m*m;
//      write (6,*) 'latest iterate x =';
//      write (6,'((1x,1p, 6(1x,d11.4)))') (wa(i),i = j+1,j+n);
//      write (6,'(a,1p,d12.5,4x,a,1p,d12.5)');
//         'at latest iterate   f =',dsave(2),'|proj g| =',dsave(13);
      }
      else {
        f=.25e0*pow2(x(1)-1.e0);
        for(int i=2;i<=n;i++) {
          f=f+pow2(x(i)-pow2(x(i-1)));
        }
        f=4.e0*f;
        t1=x(2)-pow2(x(1));
        g(1)=2.e0*(x(1)-1.e0)-1.6e1*x(1)*t1;
        for(int i=2;i<=n-1;i++) {
          t2=t1;
          t1=x(i+1)-pow2(x(i));
          g(i)=8.e0*t2-1.6e1*x(i)*t1;
        }
        g(n)=8.e0*t1;
      }
      goto lbl_111;
    }
    if (task.substr(0,5) == "NEW_X") {
      if (isave(34) >= 900) {
        task="STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT";
      }
      if (dsave(13) <= 1.e-10*(1.0e0 + absolute(f))) {
        task="STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL";
      }
//    write (6,'(2(a,i5,4x),a,1p,d12.5,4x,a,1p,d12.5)') 'iterate';
//       ,isave(30),'nfg =',isave(34),'f =',f,'|proj g| =',dsave(13);
      if (task.substr(0,4) == "STOP") {
//      write (6,*) task;
//      write (6,*) 'final x=';
//      write (6,'((1x,1p, 6(1x,d11.4)))') (x(i),i = 1,n);
      }
      goto lbl_111;
    }
  }

} // namespace <anonymous>

int main()
{
  try {
    driver3();
  }
  catch (std::exception const& e) {
    std::cout << e.what() << std::endl;
  }
  return 0;
}
