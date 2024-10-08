// Standard Neo-Hookean formulation in current configuration with Enzyme
#include "../include/nh-enzyme.h"

int main() {
        // Allocate memory for the parameters
        NeoHookeanElasticityEnzymeParams *context = (NeoHookeanElasticityEnzymeParams *)malloc(sizeof(NeoHookeanElasticityEnzymeParams));
        context->mu = 1;
        context->lambda = 1;

        // Residual Evaluation
        const double dXdx_initial[3][3] = {
          {0.0702417, 0.4799115, 0.3991242},
          {0.6756593, 0.0633284, 0.0959267},
          {0.2241923, 0.0281781, 0.0917613}
        };
        const double dudX[3][3] = {
          {0.1425560,  0.115120,  0.551640},
          {0.0591922,  0.123535,  0.166572},
          {0.1617210,  0.478828,  0.646217}
        };
        double dXdx[3][3], e_sym[6], f1[3][3];
        f1_NeoHookeanCurrentAD_Enzyme(context, dXdx_initial, dudX, dXdx,  e_sym, f1);

        // Jacobian Evaluation
        double ddudX[3][3] = {
            {0.23617975,  0.60250516,  0.1717169},
            {0.86615524,  0.3365063,   0.17371375},
            {0.0441905,   0.16762188,  0.45047923}
        };
        double df1[3][3];
        df1_NeoHookeanCurrentAD_Enzyme(context, dXdx, e_sym, ddudX, df1);
        free(context);

        return 0;
}

/*
e =
        0.245019612050
        0.050858346002
        0.307230477361
        0.125320141114
        0.365707271890
        0.126198905768

Strain Energy from e =  0.338798323727

tau =
        0.968541661897
        0.580219129800
        1.092963392519
        0.250640282229
        0.731414543779
        0.252397811536

de =
        0.594073472562
        0.520846110579
        0.197688936164
        0.375403844308
        0.357021134303
        0.356184921415

dtau  =
        2.662047097805
        2.515592373839
        1.869278025008
        0.750807688616
        0.714042268606
        0.712369842831
*/
