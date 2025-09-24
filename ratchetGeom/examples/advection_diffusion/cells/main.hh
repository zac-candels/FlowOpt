#include <math.h>
#include <stdlib.h>

#include <lbm.hh>

// This script simulates a four component layered poiseuille flow setup.
// You can modify the densities and relaxation times via the 'setDensities' and 'setTaus' functions.
// You can modify the body force magnitude in the setMagnitudeX function

int lx = 200;             // Size of domain in x direction
int ly = 200;             // Size of domain in y direction
int lz = 50;
int timesteps = 5000;     // Number of iterations to perform
int saveInterval = 1000;  // Interval to save global data
double radius = 25.0;
double theta = 90.0;
double A = 0.003;
double kappa = A * 9 / 8;
double surfacetension = sqrt(2 * A * kappa) / 6;
double interfacewidth = sqrt(8 * kappa / A);
double offsetx = 0.0;
double offsety = -ly / 2. + 12;

double tau = 1;

double diffusivity = 0.008;
double Csat = 0.3;
double Cwall = 0.0;

std::string datadir = "data/";
double inflowmomentum = 0;
double pDiff = 0;
int equilibriumtimesteps = 0;
int postwidth = lx / 4;
int postheight = 96;
int capwidth = 0;
int capthickness = 0;
int factor = 1;
double nbpost = 2;
double postfrac = 0.25;
int spacex = 106;

double sink = -0.0001;

InputParameters params;

#ifdef MPIPARALLEL
using Lattice = LatticePropertiesRuntime<ParallelX<3>, 2>;
#else
using Lattice = LatticePropertiesRuntime<NoParallel, 2>;
#endif

// Function used to define the solid geometry
// Here we set a solid at the top and bottom, in the conditions that return 1;
int initBoundary(const int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, lz, k);
    int zz = computeZ(ly, lz, k);
    //double rr2 = (xx - (lx-1)/2. - offsetx) * (xx - (lx-1)/2. - offsetx) + (yy - (ly-1)/2. - offsety) * (yy - (ly-1)/2. - offsety);
    if (yy>=postheight && yy<ly - 1*factor && (xx <= 2*factor-1)) return 3;
    if (yy>=postheight && yy<ly - 1*factor && (xx >= lx - 3*factor+1)) return 4;
    if ((xx <= 2*factor-1)) return 33;
    if ((xx >= lx - 3*factor+1)) return 44;
    if (yy < 1*factor) return 1;
    if (yy >= ly - 1*factor) return 1;

    if (((int)(xx+postfrac * (lx-spacex) / nbpost/2.-spacex/2.) % (int)((lx-spacex) / nbpost) < postfrac * (lx-spacex) / nbpost||xx>lx-spacex/2.) && // Take the remainder after dividing the current x coordinate by the lx/nbPost.
                                                           // Postfraction of this should be a post, (1-postfraction) should be fluid. // Same for z
        yy < postheight && zz<=lz/2.) // y less than post height
    {
        return 1;
    }
    if (((int)(xx+postfrac * (lx-spacex) / nbpost/2.-spacex/2.) % (int)((lx-spacex) / nbpost) > postfrac * (lx-spacex) / nbpost+4) && // Take the remainder after dividing the current x coordinate by the lx/nbPost.
    ((int)(xx+postfrac * (lx-spacex) / nbpost/2.-spacex/2.) % (int)((lx-spacex) / nbpost) < postfrac * (lx-spacex) / nbpost+2*postfrac * (lx-spacex) / nbpost-3) &&                                                   // Postfraction of this should be a post, (1-postfraction) should be fluid. // Same for z
        yy < postheight && yy>4 && zz<=lz/2.) // y less than post height
    {
        return 8;
    }
    if((xx<spacex/2.||xx>lx-spacex/2.) && yy < postheight) return 1;
    //if(xx>=0.75 * lx / 2&&yy < postheight) return 1;

    return 0;
}

double initDensity(int k) {
    int x = computeXGlobal<Lattice>(k);
    double t = x / (lx - 1.0);
    double rhoIn = 1.0 + pDiff/2;
    double rhoOut = 1.0 - pDiff/2;
    return rhoIn*(1-t) + rhoOut*t;
}

double initVelocity(int k) {

    return 0;
}

double initConcentration(int k) {
    int xx = computeXGlobal<Lattice>(k);

    if ((xx <= 2*factor-1)) return Csat;
    else return 0;
    /*if (yy < postheight + offsety)
        return Csat;
    else {
        
        double t = (xx-postheight - offsety) / (lx - 4*factor - postheight - offsety);
        double rhoIn = Csat;
        double rhoOut = 0;
        return rhoIn*(1-t) + rhoOut*t;
    }*/
}


void initParams(std::string inputfile) {
    params.addParameter<int>(lx, "lx");
    params.addParameter<int>(ly, "ly");
    params.addParameter<int>(lz, "lz");
    params.addParameter<int>(timesteps, "timesteps");
    params.addParameter<int>(saveInterval, "saveInterval");
    params.addParameter<double>(radius, "radius");
    params.addParameter<double>(theta, "theta");
    params.addParameter<double>(A, "A");
    params.addParameter<double>(kappa, "kappa");
    params.addParameter<double>(surfacetension, "surfacetension");
    params.addParameter<double>(interfacewidth, "interfacewidth");

    params.addParameter<double>(diffusivity, "diffusivity");
    params.addParameter<double>(Csat, "Csat");
    params.addParameter<double>(Cwall, "Cwall");
    params.addParameter<double>(offsetx, "offsetx");
    params.addParameter<double>(offsety, "offsety");
    params.addParameter<double>(tau, "tau");
    params.addParameter<double>(nbpost, "nbpost");
    params.addParameter<double>(postfrac, "postfrac");
    params.addParameter<std::string>(datadir, "datadir");
    params.addParameter<int>(equilibriumtimesteps, "equilibriumtimesteps");
    params.addParameter<int>(postwidth, "postwidth");
    params.addParameter<int>(postheight, "postheight");
    params.addParameter<double>(pDiff,"pDiff");
    params.addParameter<int>(capwidth, "capwidth");
    params.addParameter<int>(capthickness, "capthickness");
    params.addParameter<int>(factor, "factor");
    params.addParameter<int>(spacex, "spacex");
    params.addParameter<double>(sink, "sink");


    params.readInput(inputfile);

    #ifndef KINETIC
        diffusivity*=factor;
    #endif
    postwidth *= factor;
    postheight *= factor;
    capwidth *= factor;
    capthickness *= factor;
    lx *= factor;
    ly *= factor;
    offsetx *= factor;
    offsety *= factor;
    pDiff *= 1.0/((double)factor);

    A = 12 * surfacetension / interfacewidth;
    kappa = pow(interfacewidth, 2) / 8.0 * A;

    Lattice::init(lx, ly, lz);
}



using AdvectionDiffusionTrait =
    typename DefaultTrait<Lattice>::template SetForce<AdvectionDiffusionCorrection<Concentration<>, ConcentrationOld<>>>
    ::template SetBoundary<Dirichlet, Convective, FreeSlip, Neumann<Concentration<>>>;



template <typename TTrait = typename AdvectionDiffusionTrait::SetDataType<DataOldNewEquilibrium>>
AdvectionDiffusionCorrected<Concentration<>,ConcentrationOld<>,Lattice, TTrait> initAD() {
    AdvectionDiffusionCorrected<Concentration<>,ConcentrationOld<>,Lattice, TTrait> concentration;

    concentration.setCollideID({0});

    concentration.template getBoundary<FreeSlip>().setNodeID(1);

    concentration.template getBoundary<Dirichlet>().setNodeID({ 3,33});

    concentration.template getBoundary<Neumann<Concentration<>>>().setNodeID({ 8 });
    concentration.template getBoundary<Neumann<Concentration<>>>().setNormalGradient({ sink });
    concentration.template getBoundary<Neumann<Concentration<>>>().cutoffNegative(true);

    concentration.template getBoundary<Convective>().setNodeID({4,44, 7});

    concentration.setDiffusivity(diffusivity);

    return concentration;
}


using TraitFlowField = typename DefaultTrait<Lattice>::
    template SetProcessor<Swapper<Velocity<>, VelocityOld<>, Lattice::NDIM>>:: template SetBoundary<std::tuple<BounceBack>,std::tuple<ZouHeDensity>,std::tuple<BounceBack>>;;//,Convective.

template <typename TTrait =
              typename TraitFlowField::template SetDataType<DataOldNewEquilibrium>::template SetCollisionOperator<MRT> >
              auto initPressure() {
    FlowField<Lattice, TTrait> flowfield;

    flowfield.setCollideID({0, 3,4});
    //flowfield.template getBoundary<FreeSlip>().setNodeID({2});
    //flowfield.template getBoundary<FreeSlip,2>().setNodeID({2});
    flowfield.template getBoundary<BounceBack>().setNodeID({1,44,33,8});
    flowfield.template getBoundary<BounceBack, 2>().setNodeID({1,44,33,8});
    flowfield.template getBoundary<ZouHeDensity>().setNodeID({3,4});
    flowfield.template setTau(tau);

    return flowfield;
}