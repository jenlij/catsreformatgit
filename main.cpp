//
//  main.cpp
//  CatsReformat
//
//  Created by Jennifer Johnson on 5/2/15.
//  Copyright (c) 2015 Jennifer Johnson. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include "MathConstant.h"
#include "Parameters.h"
#include "DreUtils.h"

#include <cstdlib>
#include <stdio.h>
#include <math.h>
using namespace std;
int main(int argc, char ** argv)
{
    printf("CaTS Power Calculator Command Line Edition\n");
    printf("(c) 2006 Andrew Skol and Goncalo Abecasis\n");
    printf("(*) formerly known as PowDre -- Andrew's Power Calculator\n\n");
    
    
    int ncases = 1000;
     int ncontrols = 1000;
     
     double pisamples = 1;
     double pimarkers = 1;
     
     
     double freq = 0.5;
     double risk = 1.5;
     double prevalence = 0.1;
     double alpha = 0.0000070;
    
    /*int ncases = 1000;
    int ncontrols = 1000;
    
    double pisamples = 1;
    double pimarkers = 1;
    
    
    double freq = 0.5;
    double risk = 1.5;
    double prevalence = 0.1;
    double alpha = 0.05 / 300000;*/
    
    bool   additive = false;
    bool   dominant = false;
    bool   recessive = false;
    bool   multiplicative = true;
    
    ParameterList pl;
    
    BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("Sample Size")
    LONG_INTPARAMETER("cases", &ncases)
    LONG_INTPARAMETER("controls", &ncontrols)
    LONG_PARAMETER_GROUP("Study Design")
    LONG_DOUBLEPARAMETER("pisamples", &pisamples)
    LONG_DOUBLEPARAMETER("pimarkers", &pimarkers)
    LONG_DOUBLEPARAMETER("alpha", &alpha)
    LONG_PARAMETER_GROUP("Disease Allele")
    LONG_DOUBLEPARAMETER("frequency", &freq)
    LONG_DOUBLEPARAMETER("risk", &risk)
    LONG_PARAMETER_GROUP("Disease Model")
    EXCLUSIVE_PARAMETER("additive", &additive)
    EXCLUSIVE_PARAMETER("dominant", &dominant)
    EXCLUSIVE_PARAMETER("recessive", &recessive)
    EXCLUSIVE_PARAMETER("multiplicative", &multiplicative)
    LONG_DOUBLEPARAMETER("prevalence", &prevalence)
    END_LONG_PARAMETERS();
    
    pl.Add(new LongParameters("Additional Parameters", longParameters));
    pl.Read(argc, argv);
    pl.Status();
    
    double C     = - ninv(alpha * 0.5);
    /*double C1    = - ninv(pimarkers * 0.5);
    double Crep  = - ninv(alpha / pimarkers);
    
    // Setup parameter array for our fancy functions ...
    double parameters[6];
    
    parameters[0] = pisamples;
    parameters[1] = parameters[2] = 0.0;  // Expected stage 1 and 2 statistics
    parameters[3] = C1;
    parameters[4] = 0.0;
    
    double Cjoint = solve(pjoint, parameters, alpha, Crep, C, 1e-6);
    
    printf("Recommended Thresholds are:\n"
           "   For a one-stage study                     %7.3f\n"
           "   For first stage in two-stage study:       %7.3f\n"
           "   For second stage in replication analysis: %7.3f\n"
           "   For second stage in a joint analysis:     %7.3f\n\n",
           C, C1, Crep, Cjoint); */
    
    double p[3];      // The three genotype frequencies
    double f[3];      // The three penetrances
    
    // Genotype frequencies
    p[0] = square(freq);
    p[1] = 2 * freq * (1. - freq);
    p[2] = square(1. - freq);
    
    // Relative penetrances for low risk genotype is 1.0
    f[2] = 1.0;
    
    // Check that at least one parameter is enabled
    if (!(additive || dominant || multiplicative || recessive))
    {
        printf("No genetic model selected!\n\n"
               "No command line parameter required for the default multiplicative model.\n"
               "Use the --additive, --recessive or --dominant parameters only once to\n"
               "select a different model.\n\n");
        return 1;
    }
    
    // Relative penetrances for high risk genotypes depend on model
    if (additive)       { f[0] = 2.0 * risk - 1.0;   f[1] = risk; }
    if (dominant)       { f[0] = risk;               f[1] = risk; }
    if (recessive)      { f[0] = risk;               f[1] = 1.0;  }
    if (multiplicative) { f[0] = risk * risk;        f[1] = risk; }
    
    double scale = prevalence / (f[0]*p[0] + f[1]*p[1] + f[2]*p[2]);
    
    // Adjusted penetrances
    f[0] *= scale;
    f[1] *= scale;
    f[2] *= scale;
    
    if (f[0] > 1.0)
    {
        printf("I don't like the genetic model you requested!\n");
        return 1;
    }
    
    printf("The probability of disease as a function of genotype is ...\n");
    printf("  For genotype A/A [frequency %.3f]:   %.3f\n", p[0], f[0]);
    printf("  For genotype A/B [frequency %.3f]:   %.3f\n", p[1], f[1]);
    printf("  For genotype B/B [frequency %.3f]:   %.3f\n\n", p[2], f[2]);
    
    double pcases = (f[0] * p[0] + f[1] * p[1] * 0.5) / prevalence;
    double pcontrols = ( (1. - f[0]) * p[0] + (1. - f[1]) * p[1] * 0.5) / (1. - prevalence);
    
    printf("The expected disease allele frequency is:\n"
           "       %.3f in cases\n"
           "       %.3f in controls\n\n",
           pcases, pcontrols);
    
    double vcases = pcases * (1.0 - pcases);
    double vcontrols = pcontrols * (1.0 - pcontrols);
    
  /*  double ncp1 = (pcases - pcontrols) / sqrt( (vcases / ncases + vcontrols / ncontrols) / (pisamples * 2));
    double ncp2 = (pcases - pcontrols) / sqrt( (vcases / ncases + vcontrols / ncontrols) / ((1.0 - pisamples) * 2)); */
    double ncp  = (pcases - pcontrols) / sqrt( (vcases / ncases + vcontrols / ncontrols) * 0.5 );
  
    double P    = ndist(-C - ncp, false) + ndist(C - ncp, true);
  /*double P1   = ndist(-C1 - ncp1, false) + ndist(C1 - ncp1, true);
    double Prep = ndist(-C1 - ncp1, false) * ndist(-Crep - ncp2, false) +
    ndist(C1 - ncp1, true) * ndist(Crep - ncp2, true);
    
    // Setup the parameters for our weird integral ...
    parameters[1] = ncp1;
    parameters[2] = ncp2;
    parameters[3] = C1;
    parameters[4] = Cjoint;
    
    double upper_from = max(C1, ncp1 - 7);
    double upper_to = max(C1, ncp1) + 5;
    
    double Pjoint = integral(stuff_to_integrate, (void *) parameters, upper_from, upper_to, 1e-7) +
    integral(stuff_to_integrate, (void *) parameters, -C1 - 5, -C1, 1e-7); */
    
    /*printf("Expected Power is:\n"
           "   For a one-stage study                     %7.3f\n"
           "   For first stage in two-stage study:       %7.3f\n"
           "   For second stage in replication analysis: %7.3f\n"
           "   For second stage in a joint analysis:     %7.3f\n\n",
           P, P1, Prep, Pjoint); */
    cout << "expected power for one stage study " << P << endl;
    
    
    //AA freq, AA prob, AB freq, AB prob, BB freq, BB prob, cases, controls, power
    
    cout << "\n" << fixed << setprecision(3) << p[0] << " " << fixed << setprecision(3) << f[0] << " " << fixed << setprecision(3) << p[1] << " " << fixed << setprecision(3) << f[1] << " " << fixed << setprecision(3) << p[2] << " " << fixed << setprecision(3) << f[2] << " " << fixed << setprecision(3) << pcases << " " << fixed << setprecision(3) << pcontrols << " " << fixed << setprecision(3) << P << endl;


}



