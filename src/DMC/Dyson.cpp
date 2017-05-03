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

void DiagrammaticMonteCarlo::Dyson () {
  // define complex number
  const complex<double> I(0.0, 1.0); 

  // real Fourier transform
  FFT<double> fft;

  // which row of the histogram
  unsigned int ip = 0;

  // whole part of the self energy
  Array<double, Dynamic, Dynamic> SE_t_all = this->normalizeHistogram();

  double
    p = (ip + 0.5)*this->dp,
    EofP = 0.5*pow(p, 2.0);

  int Nt = this->maxLength/this->dt;

  // times 
  ArrayXd t = ArrayXd::LinSpaced(Nt, 0, this->maxLength - this->dt) + 0.5*this->dt;


  // frequencies
  ArrayXd w;
  double scaleFactor = 2*M_PI/(Nt*this->dt);
  if (Nt%2 == 0) {
    w = ArrayXd::LinSpaced(Nt, -Nt/2, (Nt-1)/2)*scaleFactor;
  } else {
    w = ArrayXd::LinSpaced(Nt, -(Nt-1)/2, (Nt-1)/2)*scaleFactor;
  }
  
  // bare propagator
  ArrayXd G0_t = exp((this->mu - EofP)*t);
  Array<complex<double>, Dynamic, 1> G0_w = 1/(I*w + (EofP - this->mu));

  // singular part of the self energy
  // must bel arger than one
  double param = 1;
  ArrayXd SE_t_sing = this->alpha*exp(-param*t)/sqrt(M_PI*t);
  ArrayXcd SE_w_sing = this->alpha/sqrt(param + I*w);

  // whole part of self energy corresponding to our selected momenta
  ArrayXd SE_t = SE_t_all.row(ip);

  // regular part of the self energy
  ArrayXd SE_t_reg = SE_t - SE_t_sing;
  
  // convert from eigen::Array to std::vector before the Fourier transform 
  vector<double> SE_t_reg_vec{SE_t_reg.data(), SE_t_reg.data() + SE_t_reg.size()};
  vector<complex<double> > SE_w_reg_vec;

  // forward fourier transform
  fft.fwd(SE_w_reg_vec, SE_t_reg_vec);
  fftshift1D(SE_w_reg_vec);
  Map<ArrayXcd> SE_w_reg(SE_w_reg_vec.data(), Nt);
  SE_t_reg *= this->dt;


  // vector<double> SE_t{0.54402111, 0.36459873, 0.17034683, -0.03083368, -0.23076008, -0.42130064, -0.59470541, -0.74392141, -0.86287948, -0.94674118, -0.99209556, -0.99709789, -0.96154471, -0.8868821, -0.77614685, -0.63384295, -0.46575841, -0.27872982, -0.0803643, 0.12126992, 0.31797166, 0.50174037, 0.66510151, 0.80141062, 0.90512352, 0.97202182, 0.99938456, 0.98609877, 0.93270486, 0.84137452, 0.7158225, 0.56115544, 0.38366419, 0.19056796, -0.01027934, -0.21070855, -0.40256749, -0.57805259, -0.73002623, -0.85230712, -0.93992165, -0.98930624, -0.99845223, -0.96698762, -0.8961922, -0.78894546, -0.64960951, -0.48385164, -0.2984138, -0.10083842, 0.10083842, 0.2984138, 0.48385164, 0.64960951, 0.78894546, 0.8961922, 0.96698762, 0.99845223, 0.98930624, 0.93992165, 0.85230712, 0.73002623, 0.57805259, 0.40256749, 0.21070855, 0.01027934, -0.19056796, -0.38366419, -0.56115544, -0.7158225, -0.84137452, -0.93270486, -0.98609877, -0.99938456, -0.97202182, -0.90512352, -0.80141062, -0.66510151, -0.50174037, -0.31797166, -0.12126992, 0.0803643, 0.27872982, 0.46575841, 0.63384295, 0.77614685, 0.8868821, 0.96154471, 0.99709789, 0.99209556, 0.94674118, 0.86287948, 0.74392141, 0.59470541, 0.42130064, 0.23076008, 0.03083368, -0.17034683, -0.36459873, -0.54402111};
  // vector<complex<double> > SE_w;

  // fft.fwd(SE_w, t);
  // fft.inv(SE_, freqvec);






  // vector<double> T{t.data(), t.data() + t.size()};
  // vector<double> T{t.data(), t.data() + 4};

  // vector<complex<double> > SE_w{Nt}, SE{Nt};

  // fft.fwd(SE_w, T);



}