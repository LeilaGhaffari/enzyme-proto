adouble MatDetAM1Symmetric(adouble A_sym[6]) {
  return A_sym[0] * (A_sym[1] * A_sym[2] - A_sym[3] * A_sym[3]) +
         A_sym[5] * (A_sym[3] * A_sym[4] - A_sym[5] * A_sym[2]) +
         A_sym[4] * (A_sym[5] * A_sym[3] - A_sym[4] * A_sym[1]) +
         A_sym[0] + A_sym[1] + A_sym[2] +
         A_sym[0] * A_sym[1] + A_sym[0] * A_sym[2] + A_sym[1] * A_sym[2] -
         A_sym[5] * A_sym[5] - A_sym[4] * A_sym[4] - A_sym[3] * A_sym[3];
};

adouble MatTraceSymmetric(adouble A_sym[6]) { return A_sym[0] + A_sym[1] + A_sym[2]; }

adouble VoigtDetAM1(const adouble V[6]) {
    return V[0] * (V[1] * V[2] - V[3] * V[3]) +
           V[5] * (V[3] * V[4] - V[5] * V[2]) +
           V[4] * (V[5] * V[3] - V[4] * V[1]) +
           V[0] + V[1] + V[2] +
           V[0] * V[1] + V[0] * V[2] + V[1] * V[2] -
           V[5] * V[5] - V[4] * V[4] - V[3] * V[3];
};

adouble Log1pSeries(adouble x) {
    adouble sum = 0;
    adouble y = x / (2. + x);
    adouble y2 = y*y;
    sum += y;
    for (int i=0; i<5; i++) {
      y *= y2;
      sum += y / (2*i + 3);
    }
    return 2 * sum;
};

adouble VoigtTrace(adouble V[6]) { return V[0] + V[1] + V[2]; }

adouble StrainEnergy(adouble e_sym[6], const double lambda, const double mu) {
  adouble e2_sym[6];

  // J and log(J)
  for (int i = 0; i < 6; i++) e2_sym[i] = 2 * e_sym[i];
  adouble detbm1 = MatDetAM1Symmetric(e2_sym);
  adouble J      = sqrt(detbm1 + 1);
  adouble logJ   = Log1pSeries(detbm1) / 2.;

  // trace(e)
  adouble trace_e = MatTraceSymmetric(e_sym);

  return lambda * (J * J - 1) / 4 - lambda * logJ / 2 + mu * (-logJ + trace_e);
}

void ComputeGradPsi(double grad[6], double Xp[6], const double lambda, const double mu) {
  // Active section for AD
  int tag = 1;
  auto Ea = new adouble[n];
  auto Fa = new adouble[m];
  auto Fp = new double[m];
  trace_on(tag); // Start tracing floating point operations
  for (int i=0; i<n; i++) Ea[i] <<= Xp[i]; // Assign indXpendent variables
  Fa[0] = StrainEnergy(Ea, lambda, mu); // Evaluate the body of the differentiated code
  Fa[0] >>= Fp[0]; // Assign dXpendent variables
  trace_off();    // End of the active section

  // Compute the gradient
  gradient(tag, n, Xp, grad);
  for (int i=0; i<n; i++) if (i>2) grad[i] /= 2.;
};

void ComputeHessianPsi(double hess[6][6], double Xp[6], const double lambda, const double mu) {
    // Active section for AD
    int tag = 1;
    auto Fp = new double[m];
    adouble* Xa = new adouble[n];
    adouble* Fa = new adouble[m];
    trace_on(tag);
    for (int i=0; i<n; i++) Xa[i] <<= Xp[i];
    Fa[0] = StrainEnergy(Xa, lambda, mu);
    Fa[0] >>= Fp[0];
    trace_off();

    // Allocate data array for the lower half of the hessian matrix
    double **H = (double **)malloc(n * sizeof(double *));
    for(int i=0; i<n; i++) H[i] = (double *)malloc((i+1) * sizeof(double));

    // Compute the hessian matrix
    hessian(tag, n, Xp, H);

    // Populate hess
    for (int i=0; i<n; i++) {
      for (int j=0; j<i+1; j++) {
        hess[i][j] = H[i][j];
        if (i != j) hess[j][i] = hess[i][j];
      }
    }
    for (int i=0; i<n; i++) for (int j=0; j<n; j++) if (i > 2) hess[i][j] /= 2.;
};
