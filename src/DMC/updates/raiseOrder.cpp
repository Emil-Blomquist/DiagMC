#include "../DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::raiseOrder (double param) {
  // temp
  if (this->maxDiagramOrder && this->FD.Ds.size() == maxDiagramOrder) {
    return;
  }

  double oldVal = 0;
  if (this->debug) {
    // old configuration value
    oldVal = this->evaluateDiagram();
  }

  unsigned int i1 = 0, i2 = 0;
  double q, theta, phi, wInvQ, t1, t2, wInvt1, wInvt2, wInvG1, boldContribution = 1;
  Vector3d Q, meanP;
  shared_ptr<Electron> g1, g2;

  if ( ! this->externalLegs && this->FD.Ds.size() == 0) {
    // attach the phonon to the existing two end vertices

    t1 = 0;
    t2 = this->FD.length;

    wInvt1 = 1;
    wInvt2 = 1;
    wInvG1 = 1;

    double std = (t2 - t1 < pow(10.0, -10.0) ? 100000 : 1/sqrt(t2 - t1)); 
    normal_distribution<double> normal(0.0, std);

    q = abs(normal(this->pcg));
    theta = acos(1 - 2*this->Udouble(0, 1));
    phi = this->Udouble(0, 2*M_PI);
    wInvQ = 2*M_PI * 2 * sqrt(0.5*M_PI*pow(std, 2.0)) * exp(0.5*pow(q/std, 2.0));

    // the Z-axis defines theta
    Q << q*sin(theta)*cos(phi), q*sin(theta)*sin(phi), q*cos(theta);

    // contribution from boldification
    if (this->bold && this->boldIteration > 0) {
      boldContribution = exp(
                           this->additionalPhase(this->FD.Gs[0]->momentum - Q, t2)
                           - this->additionalPhase(this->FD.Gs[0])
                         );
    }

    meanP = this->FD.Gs[0]->momentum;
  } else {
    // select first electron line to split
    i1 = this->Uint(0, this->FD.Gs.size() - 1);
    wInvG1 = this->FD.Gs.size();
    g1 = this->FD.Gs[i1];

    // first vertex position
    double
      t1low = g1->start->position,
      t1up = g1->end->position;

    t1 = this->Udouble(t1low, t1up);
    wInvt1 = t1up - t1low;

    // second vertex position
    double
      t2low = t1,
      t2up = this->FD.length,
      l = 1,
      r = this->Udouble(0, 1),
      dt2 = -log(1 - r + r*exp(-l*(t2up - t2low)))/l;

    t2 = t2low + dt2;
    wInvt2 = exp(l*dt2)*(1 - exp(-l*(t2up - t2low)))/l;

    // find second electron line to split
    unsigned int
      lowerBound = i1,
      upperBound = this->FD.Gs.size() - 1;

    // binary search
    do {
      i2 = 0.5*(upperBound + lowerBound);

      if (this->FD.Gs[i2]->start->position > t2 && this->FD.Gs[i2]->end->position > t2) {
        upperBound = i2 - 1;
      } else if (this->FD.Gs[i2]->start->position < t2 && this->FD.Gs[i2]->end->position < t2) {
        lowerBound = i2 + 1;
      } else {
        break;
      }
    } while (i2 != upperBound);

    g2 = this->FD.Gs[i2];

    // phonon momentum
    // if dt = 0
    double std = (t2 - t1 < pow(10.0, -10.0) ? 100000 : 1/sqrt(t2 - t1)); 
    normal_distribution<double> normal(0.0, std);

    q = abs(normal(this->pcg));
    theta = acos(1 - 2*this->Udouble(0, 1));
    phi = this->Udouble(0, 2*M_PI);
    wInvQ = 2*M_PI * 2 * sqrt(0.5*M_PI*pow(std, 2.0)) * exp(0.5*pow(q/std, 2.0));

    // the Z-axis defines theta
    Q << q*sin(theta)*cos(phi), q*sin(theta)*sin(phi), q*cos(theta);

    // calculate meanP
    if (i1 == i2) {
      meanP = g1->momentum;
    } else {
      meanP = g1->momentum*(g1->end->position - t1) + g2->momentum*(t2 - g2->start->position);
      meanP += this->calculateMeanP(g1->end, g2->start)*(g2->start->position - g1->end->position);
      meanP /= (t2 - t1);
    }

    // contribution from boldification
    if (this->bold && this->boldIteration > 0) {
      if (i1 == i2) {
        boldContribution = exp(
                             this->additionalPhase(g1->p, t1 - t1low)
                             + this->additionalPhase(g1->momentum - Q, t2 - t1)
                             + this->additionalPhase(g1->p, t1up - t2)
                             - this->additionalPhase(g1)
                           );
      } else {

        // contribution from lower electron which is to be split
        boldContribution = this->additionalPhase(g1->p, t1 - t1low)
                         + this->additionalPhase(g1->momentum - Q, t1up - t1)
                         - this->additionalPhase(g1);

        // contribution from upper electron which is to be split
        boldContribution += this->additionalPhase(g2->momentum - Q, t2 - g2->start->position)
                          + this->additionalPhase(g2->p, g2->end->position - t2)
                          - this->additionalPhase(g2);

        // rest of the electrons under the phonon
        shared_ptr<Vertex> v = g1->end;
        while (v != g2->start) {
          boldContribution += this->additionalPhase(v->G[1]->momentum - Q, v->G[1]->end->position - v->position)
                            - this->additionalPhase(v->G[1]);
          
          v = v->G[1]->end;
        }

        boldContribution = exp(boldContribution);
      }
    }
  }

  double
    wInvd = this->FD.Ds.size() + 1,
    sinTheta = sin(theta),
    alpha = this->alpha,
    dt = t2 - t1,
    exponential = exp(-dt*(this->FD.phononEnergy(q) + 0.5*q*q - Q.dot(meanP)));

  double a;
  if ( ! isfinite(exponential)) {
    a = 1;
  } else if ( ! isfinite(boldContribution)) {
    a = 1;
  } else if (sinTheta == 0 || wInvt2 == 0 || wInvQ == 0) {
    a = 0;
  } else {
    a = boldContribution * exponential * 0.5*sqrt(M_PI)/alpha * alpha/(sqrt(8)*M_PI*M_PI) * (wInvG1*wInvt1*wInvt2*wInvQ)/wInvd;
  }

  // accept or reject update
  bool accepted = false;
  if (a > this->Udouble(0, 1)) {
    if ( ! this->externalLegs && this->FD.Ds.size() == 0) {
      // add phonon
      shared_ptr<Phonon> d = this->FD.addInternalPhonon(this->FD.start, this->FD.end, Q, q, theta, phi);
    } else {
      // split first electron line
      shared_ptr<Vertex> v1 = this->FD.insertVertex(i1, t1);

      // split second electron line
      shared_ptr<Vertex> v2 = this->FD.insertVertex(i2 + 1, t2);

      // add phonon
      shared_ptr<Phonon> d = this->FD.addInternalPhonon(v1, v2, Q, q, theta, phi);
    }

    this->FD.setNewStructure();
    accepted = true;
  }


  // cout << "--------------" << endl
  //      << "raiseOrder" << endl
  //      << "\tn=" << this->FD.Ds.size() << endl
  //      << "\tt=" << this->FD.length << endl
  //      << "\tp=" << this->FD.externalMomentum << endl
  //      << "\ta=" << a << endl
  //      << "\taccepted=" << accepted << endl
  //      << "--------------" << endl;



  if (this->debug) {
    double val = this->evaluateDiagram();

    if (accepted) {
      this->checkAcceptanceRatio(boldContribution * exponential * alpha*sinTheta/(sqrt(8)*M_PI*M_PI)/(val/oldVal), "raiseOrder");
    }

    if (a < 0 || ! isfinite(a)) {
      cout << "--------------------------------------------------------------------" << endl
           << "overflow at DMC::raiseOrder " << endl
           << "accepted=" << accepted << endl
           << "a=" << a << endl
           << "a_diag=" << val/oldVal * (wInvG1*wInvt1*wInvt2*wInvQ)/wInvd << endl
           << "ratio=" << a/(val/oldVal * (wInvG1*wInvt1*wInvt2*wInvQ)/wInvd) << endl
           << "order=" << this->FD.Ds.size() << endl
           << "val=" << val << endl
           << "oldVal=" << oldVal << endl
           << "wInvQ=" << wInvQ << endl
           << "wInvt2=" << wInvt2 << endl
           << "--------------------------------------------------------------------" << endl;
    } else if (this->loud) {
      cout << "raiseOrder: " << accepted << " " << a << " " << boldContribution * exponential * alpha*sinTheta/(sqrt(8)*M_PI*M_PI)/(val/oldVal) << endl;
    }
  }
}