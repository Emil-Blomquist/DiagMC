#include "../DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::BOLDlowerOrder (double param) {
  // requirement to lower: must be at least of order 1
  if (this->FD.Ds.size() == 0) {
    return;
  }

  // old configuration value
  double oldVal = 0;
  if (this->debug) {
    oldVal = this->evaluateDiagram();
  }

  // phonon
  int phononIndex;
  shared_ptr<Phonon> d;

  // vertices
  shared_ptr<Vertex> v1, v2;

  double wInvQ, wInvt1, wInvt2, wInvG1, boldContribution = 0;

  if ( ! this->externalLegs && this->FD.Ds.size() == 1) {
    // go down to zeroth order

    phononIndex = 0;
    d = this->FD.Ds[0];
    v1 = d->start;
    v2 = d->end;

    wInvt1 = 1;
    wInvt2 = 1;
    wInvG1 = 1;

    double std = (v2->position - v1->position < pow(10.0, -10.0) ? 100000 : 1/sqrt(v2->position - v1->position)); 
    normal_distribution<double> normal(0.0, std);

    wInvQ = 2*pow(M_PI, 2.0) * sqrt(0.5*M_PI*pow(std, 2.0)) * exp(0.5*pow(d->q/std, 2.0));

    // contribution from boldification
    boldContribution = this->additionalPhase(this->FD.Gs[0]->momentum + d->momentum, this->FD.length)
                     - this->additionalPhase(this->FD.Gs[0]);

  } else {
    // pick a phonon line
    phononIndex = this->Uint(0, this->FD.Ds.size() - 1);
    d = this->FD.Ds[phononIndex];

    // vertices
    v1 = d->start;
    v2 = d->end;

    if ( ! this->externalLegs && (v1 == this->FD.start || v2 == this->FD.end)) {
      // invalid diagram in order to lower
      return;
    }

    // calculate vertex probability
    wInvG1 = this->FD.Gs.size() - 2;

    // inverse probabilities for internal parameters
    // must be the same probabilities as for raising order!
    wInvt1 = (v1->G[1]->end == v2 ? v2->G[1]->end->position : v1->G[1]->end->position) - v1->G[0]->start->position;

    double
      t2low = v1->position,
      t2up = this->FD.length;

    wInvt2 = t2up - t2low;

    double std = (v2->position - v1->position < pow(10.0, -10.0) ? 100000 : 1/sqrt(v2->position - v1->position));

    wInvQ = 2*pow(M_PI, 2.0) * sqrt(0.5*M_PI*pow(std, 2.0)) * exp(0.5*pow(d->q/std, 2.0));

    // contribution from boldification
    if (v1->G[1]->end == v2) {
      boldContribution = this->additionalPhase(v1->G[0]->p, v2->G[1]->end->position - v1->G[0]->start->position)
                       - this->additionalPhase(v1->G[0])
                       - this->additionalPhase(v1->G[1])
                       - this->additionalPhase(v2->G[1]);
    } else {
      // contribution from lower electron which is to be split
      boldContribution = this->additionalPhase(v1->G[0]->p, v1->G[1]->end->position - v1->G[0]->start->position)
                       - this->additionalPhase(v1->G[0])
                       - this->additionalPhase(v1->G[1]);

      // contribution from upper electron which is to be split
      boldContribution += this->additionalPhase(v2->G[1]->p, v2->G[1]->end->position - v2->G[0]->start->position)
                        - this->additionalPhase(v2->G[0])
                        - this->additionalPhase(v2->G[1]);

      // rest of the electrons under the phonon
      shared_ptr<Vertex> v = v1->G[1]->end;
      while (v != v2->G[0]->start) {
        boldContribution += this->additionalPhase(v->G[1]->momentum + d->momentum, v->G[1]->end->position - v->position)
                          - this->additionalPhase(v->G[1]);
        
        v = v->G[1]->end;
      }
    }
  }

  // to calculate the acceptance ratio
  Vector3d
    Q = d->momentum,
    P0 = this->calculateP0(d);

  double
    wInvd = this->FD.Ds.size(),
    sinTheta = sin(d->theta),
    alpha = this->alpha,
    q2 = d->q*d->q,
    dt = v2->position - v1->position,
    exponential = exp(boldContribution + dt*(this->FD.phononEnergy(d->q) + 0.5*q2 - Q.dot(P0)));

  double a;
  if (wInvQ == 0 || sinTheta == 0 || wInvt2 == 0 || ! isfinite(exponential)) {
    a = 1;
  } else {
    a = exponential * sqrt(8)*M_PI*M_PI/(alpha*sinTheta) * wInvd/(wInvG1*wInvt1*wInvt2*wInvQ);
  }


  // accept or reject update
  bool accepted = false;
  if (a > this->Udouble(0, 1)) {
    // remove phonon
    this->FD.removeInternalPhonon(phononIndex);

    if (
      this->externalLegs ||
      ( ! this->externalLegs && this->FD.Ds.size() > 0))
    {
      // remove vertices
      this->FD.removeVertex(v1);
      this->FD.removeVertex(v2);
    }

    this->FD.setNewStructure();
    accepted = true;
  }

  if (this->debug) {
    double val = this->evaluateDiagram();

    if (accepted) {
      this->checkAcceptanceRatio(exponential * sqrt(8)*M_PI*M_PI/(alpha*sinTheta) /(val/oldVal), "BOLDlowerOrder");
    }

    if (a < 0 || ! isfinite(a)) {
      cout << "--------------------------------------------------------------------" << endl
           << "overflow at DMC::BOLDlowerOrder " << endl
           << "accepted=" << accepted << endl
           << "a=" << a << endl
           << "a_diag=" << val/oldVal * wInvd/(wInvG1*wInvt1*wInvt2*wInvQ) << endl
           << "ratio=" << a/(val/oldVal * wInvd/(wInvG1*wInvt1*wInvt2*wInvQ)) << endl
           << "order=" << this->FD.Ds.size() << endl
           << "val=" << val << endl
           << "oldVal=" << oldVal << endl
           << "wInvQ=" << wInvQ << endl
           << "wInvt2=" << wInvt2 << endl
           << "--------------------------------------------------------------------" << endl;
    } else if (this->loud) {
      cout << "BOLDlowerOrder: " << accepted << " " << a << " " << exponential * sqrt(8)*M_PI*M_PI/(alpha*sinTheta) /(val/oldVal) << endl;
    }
  }
}