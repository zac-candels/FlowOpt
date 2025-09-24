#include <math.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <lbm.hh>

// This script simulates a four component layered poiseuille flow setup.
// You can modify the densities and relaxation times via the 'setDensities' and 'setTaus' functions.
// You can modify the body force magnitude in the setMagnitudeX function

#ifdef THREEDIMENSIONAL
const int ndim = 3;
#else
const int ndim = 2;
#endif

//#define NUMCOMPONENTS 6
#ifndef NUMCOMPONENTS
#define NUMCOMPONENTS 6
#endif

const int numcomponents = NUMCOMPONENTS;

int lx = 400;  // Size of domain in x direction
int ly = 100;  // Size of domain in y direction
int lz = 1;
int timesteps = 5000;     // Number of iterations to perform
int saveInterval = 1000;  // Interval to save global data
std::string datadir = "data/";
double surfacetension = 0.005;
std::vector<double> taus(numcomponents,1);
std::vector<double> dens(numcomponents,1);
int ncomp = 0.0;
double prop1; double prop2; double prop3; double prop4; double prop5; double prop6; double prop7; double prop8; double prop9; double prop10; double prop11; double prop12;

InputParameters params;

#ifdef MPIPARALLEL
using Lattice = LatticePropertiesRuntime<ParallelX<2>, ndim>;
#else
using Lattice = LatticePropertiesRuntime<NoParallel, ndim>;
#endif

int initBoundary(const int k) {
    return 0;
}

void initParams(std::string inputfile) {
    params.addParameter<int>(lx, "lx");
    params.addParameter<int>(ly, "ly");
    params.addParameter<int>(lz, "lz");
    params.addParameter<int>(timesteps, "timesteps");
    params.addParameter<int>(saveInterval, "saveInterval");
    params.addParameter<double>(surfacetension, "surfacetension");
    for (int n = 0; n<numcomponents; n++) {
        params.addParameter<double>(dens[n], "dens"+std::to_string(n));
        params.addParameter<double>(taus[n], "tau"+std::to_string(n));
    }
    params.addParameter<std::string>(datadir, "datadir");
    params.addParameter<int>(ncomp, "ncomp");
    params.addParameter<double>(prop1, "prop1");
    params.addParameter<double>(prop2, "prop2");
    params.addParameter<double>(prop3, "prop3");
    params.addParameter<double>(prop4, "prop4");
    params.addParameter<double>(prop5, "prop5");
    params.addParameter<double>(prop6, "prop6");
    params.addParameter<double>(prop7, "prop7");
    params.addParameter<double>(prop8, "prop8");
    params.addParameter<double>(prop9, "prop9");
    params.addParameter<double>(prop10, "prop10");
    params.addParameter<double>(prop11, "prop11");
    params.addParameter<double>(prop12, "prop12");

    params.readInput(inputfile);

    Lattice::init(lx, ly, lz);
}


auto initPressure() {
    
    using Trait = DefaultTraitPressureWellBalancedN2<numcomponents,Lattice>;
    FlowFieldPressureWellBalanced2<Lattice, Trait> pressure;

    pressure.setTaus(taus);

    pressure.setDensities(dens);

    pressure.setCollideID({0});

    return pressure;
}

template <int TN>
auto initCH() {

    using Trait = typename std::conditional<
        TN == 0, typename DefaultTraitWellBalancedCH<0,numcomponents,Lattice>::template AddProcessor<ChemicalPotentialCalculatorNCompBoyerConstS2<numcomponents>>,
                 DefaultTraitWellBalancedCH<TN,numcomponents,Lattice>>::type;

    WellBalancedCH<TN, numcomponents, Lattice, Trait> ch;

    ch.setCollideID({0});
    
    if constexpr (TN==0) {
        ch.template getProcessor<ChemicalPotentialCalculatorNCompBoyerConstS2<numcomponents>>().setSurfaceTension(surfacetension);
        ch.template getProcessor<ChemicalPotentialCalculatorNCompBoyerConstS2<numcomponents>>().setD(3);
    }

    double mobilityprefactor = 4*surfacetension;
    std::vector<double> mobility(numcomponents,0);

    for (int n = 0; n<numcomponents; n++) {
        if(n==TN) mobility[n] = - (numcomponents - 1) / (numcomponents * surfacetension) * mobilityprefactor;
        else mobility[n] = 1 / (numcomponents * surfacetension) * mobilityprefactor;
    }

    ch.setAij(mobility);

    return ch;
}