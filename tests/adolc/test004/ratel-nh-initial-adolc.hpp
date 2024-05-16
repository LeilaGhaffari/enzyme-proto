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
