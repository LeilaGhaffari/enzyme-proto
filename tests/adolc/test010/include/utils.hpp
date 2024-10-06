struct GradientData {
    double *Ep;
    double* Fp;
    adouble* Ea;
    adouble* Fa;
};

// Define a struct to hold the variables
struct HessianData {
    double *Ep;
    double* Fp;
    adouble* Ea;
    adouble* Fa;
    double** H;
};

void InitializeGradientData(GradientData *data) {
    data->Ep = new double[6];
    data->Fp = new double[1];
    data->Ea = new adouble[6];
    data->Fa = new adouble[1];
}

// Initialize and allocate memory for the struct
void initializeHessianData(HessianData *data) {
    data->Ep = new double[6];
    data->Fp = new double[1];
    data->Ea = new adouble[6];
    data->Fa = new adouble[1];
    data->H = (double**)malloc(6 * sizeof(double*));
    for (int i = 0; i < 6; i++) {
        data->H[i] = (double*)malloc((i + 1) * sizeof(double));
    }
}

// Free the memory allocated in the struct
void freeHessianData(HessianData *data) {
    delete[] data->Ep;
    delete[] data->Fp;
    delete[] data->Ea;
    delete[] data->Fa;
    for (int i = 0; i < 6; i++) {
        free(data->H[i]);
    }
    free(data->H);
    delete data;
}

// Free the memory allocated in the struct
void freeGradientData(GradientData *data) {
    delete[] data->Ep;
    delete[] data->Fp;
    delete[] data->Ea;
    delete[] data->Fa;
}
