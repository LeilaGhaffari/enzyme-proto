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

adouble StrainEnergy(adouble E_Voigt[6], const double lambda, const double mu) {
  // Calculate 2*E
  adouble E2_Voigt[6];
  for(int i = 0; i<6; i++) E2_Voigt[i] = E_Voigt[i] * 2;

  // log(J)
  adouble detCm1 = VoigtDetAM1(E2_Voigt);

  //double J      = sqrt(detCm1 + 1);
  adouble logJ   = Log1pSeries(detCm1) / 2.;

  // trace(E)
  adouble traceE = VoigtTrace(E_Voigt);

  return lambda*logJ*logJ/2  + mu * (-logJ + traceE);
};

double *Kirchhofftau_sym_NeoHookean_AD_ADOLC(const double lambda, const double mu, double e_p[6]) {
    int tag = 1;
    auto e_a = new adouble[n];
    auto F_a = new adouble[m];
    auto F_p = new double[m];
    double **e_sym = myalloc(n, p);
    double **F = myalloc(m, p);
    double dPsi_sym[n];

    // Start tracing floating point operations
    trace_on(tag);
    // Assign independent variables
    for (int i=0; i<n; i++) e_a[i] <<= e_p[i];
    // Evaluate the body of the differentiated code
    F_a[0] = StrainEnergy(e_a, lambda, mu);
    // Assign dependent variables
    F_a[0] >>= F_p[0];
    trace_off();    // End of the active section

    // Identity matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) if (i == j) e_sym[i][j] = 1.;
        else  e_sym[i][j] = 0.00;
    }

    // dPsi / de
    fov_forward(tag, m, n, p, e_p, e_sym, F_p, F);
    for (int i=0; i<n; i++) {
      dPsi_sym[i] = F[0][i];
      if (i>2) dPsi_sym[i] /= 2.;
    }

    return tau_from_dPsi(dPsi_sym, e_p);
};

double *dtau_ADOLC(double e_p[6], const double lambda, const double mu) {
    int tag = 1;
    auto F_p = new double[m];
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
    for (int i=0; i<n; i++) E[i] <<= e_p[i];
    F[0] = StrainEnergy(E, lambda, mu);
    F[0] >>= F_p[0];
    trace_off();

    int dim = binomi(p + d, d);
    dtensor = myalloc2(m, dim);
    tensor_eval(tag, m, n, d, p, e_p, dtensor, I);

    for (int i=0; i<dim; i++) dS[i] = dtensor[0][i];

    return dS;
};
