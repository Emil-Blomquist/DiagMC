#include "DiagrammaticMonteCarlo.h"

template<typename T>
inline void fftshift1D (T& vec) {
  unsigned int pivot = (vec.size() % 2 == 0) ? (vec.size() / 2) : ((vec.size() + 1) / 2);

  T
    front{vec.begin(), vec.begin() + pivot},
    back{vec.begin() + pivot, vec.end()};

  vec.clear();

  vec.insert(vec.end(), back.begin(), back.end());
  vec.insert(vec.end(), front.begin(), front.end());

  // cannot get this to work...
  // vec.insert(vec.begin(), vec.begin() + pivot, vec.end());
  // vec.erase(vec.begin() + 2*pivot, vec.end());
}

template<typename T>
inline void ifftshift1D (T& vec) {
  unsigned int pivot = (vec.size()%2 == 0 ? vec.size()/2 : (vec.size() - 1)/2);

  T
    front{vec.begin(), vec.begin() + pivot},
    back{vec.begin() + pivot, vec.end()};

  vec.clear();

  vec.insert(vec.end(), back.begin(), back.end());
  vec.insert(vec.end(), front.begin(), front.end());
}

void DiagrammaticMonteCarlo::doDyson (
  ArrayXXd& SE_t_all,
  ArrayXXd& dG
) {
  clock_t tStart = clock();

  dG = ArrayXXd::Zero(this->hist.rows(), this->hist.cols());

  // define complex number
  const complex<double> I(0.0, 1.0); 

  // real Fourier transform
  FFT<double> fft;

  // we want to go to at least t=500
  double L = (this->maxLength > 500 ? this->maxLength : 500);
  int
    Nt = L/this->dt,
    dNt = Nt - this->maxLength/this->dt;

  // times 
  ArrayXd t = ArrayXd::LinSpaced(Nt, 0, L - this->dt) + 0.5*this->dt;

  // frequencies
  ArrayXd w;
  double scaleFactor = 2*M_PI/(Nt*this->dt);
  if (Nt%2 == 0) {
    w = ArrayXd::LinSpaced(Nt, -Nt/2, (Nt-1)/2)*scaleFactor;
  } else {
    w = ArrayXd::LinSpaced(Nt, -(Nt-1)/2, (Nt-1)/2)*scaleFactor;
  }

  // singular part of the self energy ('param' must be larger than zero)
  double param = 1;
  ArrayXd SE_t_sing = this->alpha*exp(-param*t)/sqrt(M_PI*t);
  ArrayXcd SE_w_sing = this->alpha/sqrt(param + I*w);

  for (unsigned int ip = 0; ip != this->hist.rows(); ip++) {
    // bare propagator
    double
      p = (ip + 0.5)*this->dp,
      EofP = 0.5*pow(p, 2.0);

    ArrayXd G0_t = exp((this->mu - EofP)*t);
    ArrayXcd G0_w = 1/(I*w + (EofP - this->mu));

    // whole part of self energy corresponding to our selected momenta
    // ArrayXd SE_t = SE_t_all.row(ip);
    ArrayXd SE_t{Nt};
    SE_t << ArrayXd{SE_t_all.row(ip)}, ArrayXd::Zero(dNt);

    // regular part of the self energy
    ArrayXd SE_t_reg = SE_t - SE_t_sing;
    
    // convert from eigen::Array to std::vector before the Fourier transform 
    vector<double> SE_t_reg_vec{SE_t_reg.data(), SE_t_reg.data() + SE_t_reg.size()};

    // forward Fourier transform and convert back to eigen::Array
    vector<complex<double> > SE_w_reg_vec;
    fft.fwd(SE_w_reg_vec, SE_t_reg_vec);
    fftshift1D(SE_w_reg_vec);
    Map<ArrayXcd> SE_w_reg(SE_w_reg_vec.data(), Nt);
    SE_w_reg *= (t[Nt - 1] - t[0])/Nt;

    // total self energy
    ArrayXcd SE_w = SE_w_reg + SE_w_sing;

    // propagator difference (to avoid numerical problems)
    ArrayXcd dG_w = 1/(1/G0_w - SE_w) - G0_w;

    // convert from eigen::Array to std::vector before the Fourier transform
    vector<complex<double> > dG_w_vec{dG_w.data(), dG_w.data() + dG_w.size()};

    // inverse Fourier transform and convert back to eigen::Array
    ifftshift1D(dG_w_vec);
    vector<double> dG_t_vec;
    fft.inv(dG_t_vec, dG_w_vec);
    Map<ArrayXd> dG_t(dG_t_vec.data(), Nt);
    dG_t *= Nt/(t[Nt - 1] - t[0]);

    // append this as row to dG
    dG.row(ip) = dG_t.head(Nt - dNt);
  }

  printf("[Dyson in %.2fs]\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
}