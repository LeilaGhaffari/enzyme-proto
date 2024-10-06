adouble MatDetAM1Symmetric(adouble A_sym[6]) {
  return A_sym[0] * (A_sym[1] * A_sym[2] - A_sym[3] * A_sym[3]) +
         A_sym[5] * (A_sym[3] * A_sym[4] - A_sym[5] * A_sym[2]) +
         A_sym[4] * (A_sym[5] * A_sym[3] - A_sym[4] * A_sym[1]) +
         A_sym[0] + A_sym[1] + A_sym[2] +
         A_sym[0] * A_sym[1] + A_sym[0] * A_sym[2] + A_sym[1] * A_sym[2] -
         A_sym[5] * A_sym[5] - A_sym[4] * A_sym[4] - A_sym[3] * A_sym[3];
};

adouble MatTraceSymmetric(adouble A_sym[6]) { return A_sym[0] + A_sym[1] + A_sym[2]; }

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

adouble StrainEnergy(adouble E_sym[6], const double lambda, const double mu) {
  adouble E2_sym[6];

  // J and log(J)
  for (int i = 0; i < 6; i++) E2_sym[i] = 2 * E_sym[i];
  adouble detCm1 = MatDetAM1Symmetric(E2_sym);
  adouble J      = sqrt(detCm1 + 1);
  adouble logJ   = Log1pSeries(detCm1) / 2.;

  // trace(E)
  adouble traceE = MatTraceSymmetric(E_sym);

  return lambda * (J * J - 1) / 4 - lambda * logJ / 2 + mu * (-logJ + traceE);
}

void ComputeGradPsi(double grad[6], double Xp[6], const double lambda, const double mu, GradientData *data) {
  // Active section for AD
  int tag = 1;
  trace_on(tag); // Start tracing floating point operations
  for (int i=0; i<6; i++) data->Ea[i] <<= Xp[i]; // Assign independent variables
  data->Fa[0] = StrainEnergy(data->Ea, lambda, mu); // Evaluate the body of the differentiated code
  data->Fa[0] >>= data->Fp[0]; // Assign dependent variables
  trace_off();    // End of the active section

  // Compute the gradient
  gradient(tag, 6, Xp, grad);
  for (int i=0; i<6; i++) if (i>2) grad[i] /= 2.;
};

void ComputeHessianPsi(double hess[6][6], double Xp[6], const double lambda, const double mu, HessianData *data) {
    // Active section for AD
    int tag = 1;
    trace_on(tag);
    for (int i = 0; i < 6; i++) data->Ea[i] <<= Xp[i];
    data->Fa[0] = StrainEnergy(data->Ea, lambda, mu);
    data->Fa[0] >>= data->Fp[0];
    trace_off();

    // Compute the hessian matrix
    hessian(tag, 6, Xp, data->H);

    // Populate hess
    for (int i=0; i<6; i++) {
      for (int j=0; j<i+1; j++) {
        hess[i][j] = data->H[i][j];
        if (i != j) hess[j][i] = hess[i][j];
      }
    }
    for (int i=3; i<6; i++) for (int j=0; j<6; j++) hess[i][j] /= 2.;
}
