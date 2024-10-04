
// Define a struct to hold the variables
struct HessianData {
    double* Fp;
    adouble* Xa;
    adouble* Fa;
    double** H;
};


struct GradientData {
    adouble* Ea;          // Active strain tensor
    adouble* Fa;          // Active scalar for energy
    double* Fp;           // Pointer to store the result of the energy
    double* grad;         // Pointer to store the gradient
};

// Initialize and allocate memory for the struct
void initializeHessianData(HessianData *data) {
    data->Fp = new double[1];
    data->Xa = new adouble[6];
    data->Fa = new adouble[1];

    data->H = (double**)malloc(6 * sizeof(double*));
    for (int i = 0; i < 6; i++) {
        data->H[i] = (double*)malloc((i + 1) * sizeof(double));
    }
}

void InitializeGradientData(GradientData *data) {
    // Allocate memory for the struct members
    data->Ea = new adouble[6];
    data->Fa = new adouble[1];
    data->Fp = new double[1];
    data->grad = new double[6];
}

// Free the memory allocated in the struct
void freeHessianData(HessianData *data) {
    delete[] data->Fp;
    delete[] data->Xa;
    delete[] data->Fa;

    for (int i = 0; i < 6; i++) {
        free(data->H[i]);
    }
    free(data->H);
    delete data;
}

// Free the memory allocated in the struct
void freeGradientData(GradientData *data) {
    delete[] data->Ea;
    delete[] data->Fa;
    delete[] data->Fp;
    delete[] data->grad;
}
