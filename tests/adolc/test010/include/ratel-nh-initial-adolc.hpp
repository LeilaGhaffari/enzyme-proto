adouble RatelVoigtDetAM1(const adouble V[6]) {
    return V[0] * (V[1] * V[2] - V[3] * V[3]) +
           V[5] * (V[3] * V[4] - V[5] * V[2]) +
           V[4] * (V[5] * V[3] - V[4] * V[1]) +
           V[0] + V[1] + V[2] +
           V[0] * V[1] + V[0] * V[2] + V[1] * V[2] -
           V[5] * V[5] - V[4] * V[4] - V[3] * V[3];
};

adouble RatelLog1pSeries(adouble x) {
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

adouble RatelVoigtTrace(adouble V[6]) { return V[0] + V[1] + V[2]; }

adouble StrainEnergy(adouble E_Voigt[6], const double lambda, const double mu) {
  // Calculate 2*E
  adouble E2_Voigt[6];
  for(int i = 0; i<6; i++) E2_Voigt[i] = E_Voigt[i] * 2;

  // log(J)
  adouble detCm1 = RatelVoigtDetAM1(E2_Voigt);

  //double J      = sqrt(detCm1 + 1);
  adouble logJ   = RatelLog1pSeries(detCm1) / 2.;

  // trace(E)
  adouble traceE = RatelVoigtTrace(E_Voigt);

  return lambda*logJ*logJ/2  + mu * (-logJ + traceE);
};

double *Stress(double Ep[6], const double lambda, const double mu) {
    int tag = 1;
    auto Ea = new adouble[n];
    auto Fa = new adouble[m];
    auto Fp = new double[m];
    auto S = new double[n];
    double **E = myalloc(n, p);
    double **F = myalloc(m, p);

    // Start tracing floating point operations
    trace_on(tag);
    // Assign independent variables
    for (int i=0; i<n; i++) Ea[i] <<= Ep[i];
    // Evaluate the body of the differentiated code
    Fa[0] = StrainEnergy(Ea, lambda, mu);
    // Assign dependent variables
    Fa[0] >>= Fp[0];
    trace_off();    // End of the active section

    for (int i = 0; i < n; i++) { // Identity matrix
        for (int j = 0; j < p; j++) if (i == j) E[i][j] = 1.;
        else  E[i][j] = 0.00;
    }

    fov_forward(tag, m, n, p, Ep, E, Fp, F);

    for (int i=0; i<n; i++) {
      S[i] = F[0][i];
      if (i>2) S[i] /= 2.;
    }

    return S;
};

double *dStress(double Ep[6], const double lambda, const double mu) {
    int tag = 1;
    auto Fp = new double[m];
    adouble* E = new adouble[n];
    adouble* F = new adouble[m];
    auto dS = new double[n];
    double** I = new double*[n];
    double** dtensor; // size = m x dim

    for (int i=0; i<n; i++) { // Identity matrix
        I[i] = new double[p];
        for (int j=0; j<p; j++) I[i][j] = (i == j) ? 1.0 : 0.0;
    }

    trace_on(tag);
    for (int i=0; i<n; i++) E[i] <<= Ep[i];
    F[0] = StrainEnergy(E, lambda, mu);
    F[0] >>= Fp[0];
    trace_off();

    int dim = binomi(p + d, d);
    dtensor = myalloc2(m, dim);
    tensor_eval(tag, m, n, d, p, Ep, dtensor, I);

    for (int i=0; i<dim; i++) dS[i] = dtensor[0][i];

    return dS;
};
