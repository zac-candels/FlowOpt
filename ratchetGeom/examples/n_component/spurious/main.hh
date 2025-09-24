#include <math.h>
#include <stdlib.h>

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

int lx = 400;  // Size of domain in x direction
int ly = 100;  // Size of domain in y direction
int lz = 1;
int timesteps = 5000;     // Number of iterations to perform
int saveInterval = 1000;  // Interval to save global data

double radius = 15;

double dens1 = 1;
double dens2 = 1;
double dens3 = 1;
double dens4 = 1;
std::string datadir = "data/";
double s01 = 0.001;
double s02 = 0.001;
double s03 = 0.001;
double s12 = 0.001;
double s13 = 0.001;
double s23 = 0.001;
double tau1 = 1.0;
double tau2 = 1.0;
double tau3 = 1.0;
double tau4 = 1.0;
double Lambda1 = 1.0;
double Lambda2 = 1.0;
double Lambda3 = 1.0;
double Lambda4 = 1.0;

InputParameters params;

#ifdef MPIPARALLEL
using Lattice = LatticePropertiesRuntime<ParallelX<2>, ndim>;
#else
using Lattice = LatticePropertiesRuntime<NoParallel, ndim>;
#endif

// Function used to define the solid geometry
// Here we set a solid at the top and bottom, in the conditions that return 1;


double diagfunc(double xx, double yy, double th) {
    return xx*sin(M_PI * th / 180.)-yy*tan(M_PI * th / 180.)*sin(M_PI * th / 180.)+yy/cos(M_PI * th / 180.);
}

double initFluid1(const int k) {
    double xx = computeXGlobal<Lattice>(k);
    double yy = computeY(ly, lz, k);
    int zz = computeZ(ly, lz, k);
        double rr2 = (xx - 1.5*(lx - 1) / 4.) *
                     (xx - 1.5*(lx - 1) / 4.) +
                     (yy - 1*(ly - 1) / 2.) *
                     (yy - 1*(ly - 1) / 2.);
    double rr22 = (xx - 1.5*(lx - 1) / 4.) *
                     (xx - 1.5*(lx - 1) / 4.) +
                     (yy - 1*(ly - 1) / 2.) *
                     (yy - 1*(ly - 1) / 2.);
    //return (0.5 + 0.5 * tanh(2 * (sqrt(rr2)-radius) / 4.))*(0.5 + 0.5 * tanh(2 * (yy-(ly - 1)/2.0) / 4.))*(0.5 - 0.5 * tanh(2 * (yy-(ly)) / 4.))+(0.5 - 0.5 * tanh(2 * (yy) / 4.));
    //return 0;
    //return (0.5 + 0.5 * tanh(2 * (yy-(ly + 1)/2.0) / 4.))*(0.5 + 0.5 * tanh(2 * (xx-(lx + 1)/2.0) / 4.))*(0.5 - 0.5 * tanh(2 * (yy-(ly)) / 4.))*(0.5 - 0.5 * tanh(2 * (xx-(lx)) / 4.))+(0.5 - 0.5 * tanh(2 * (yy) / 4.))*(0.5 + 0.5 * tanh(2 * (xx-(lx + 1)/2.0) / 4.))+(0.5 - 0.5 * tanh(2 * (xx) / 4.))*(0.5 + 0.5 * tanh(2 * (yy-(ly + 1)/2.0) / 4.))*(0.5 - 0.5 * tanh(2 * (yy-ly) / 4.));
    return (0.5 - 0.5 * tanh(2 * (sqrt(rr2)-radius) / 4.))*(0.5 - 0.5 * tanh(2 * (xx-1.65*(lx - 1)/4.0) / 4.));//*(1-(0.5 - 0.5 * tanh(2 * (sqrt(rr22)-radius) / 4.)));//(0.5 + 0.5 * tanh(2 * (xx-2*(lx - 1)/4.0) / 4.));//*(0.5 - 0.5 * tanh(2 * (xx-2.5*(lx - 1)/4.0) / 4.));
    return 0;//(0.5 - 0.5 * tanh(2 * (sqrt(rr2)-radius) / 4.));

    return (0.5 - 0.5 * tanh(2 * (yy-8) / 4.));//topextend>bottomextend ? ((shape>topextend) ? shape : topextend) : ((shape>bottomextend) ? shape : bottomextend) ;
    /*
    if (rr2 >= channelradius * channelradius) {
        return 1;
    } 
    else {
        return 0;
    }*///(0.5 + 0.5 * tanh(2 * (entrance) / 4.))*(0.5 + 0.5 * tanh(2 * (res0) / 4.))
}

int initBoundary(const int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, lz, k);

    return 0;
    if (yy>=ly-2||yy<2||xx>=lx-2||xx<2) return 1;
    //else if(initFluid1(k)>0.5) {
    //    return 7;
    //}
    else return 0;

}
double initFluid3(const int k) {
    double xx = computeXGlobal<Lattice>(k);
    double yy = computeY(ly, lz, k);
    double rr2 = (xx - 2.5*(lx - 1) / 4.) *
                     (xx - 2.5*(lx - 1) / 4.) +
                     (yy - 1*(ly - 1) / 2.) *
                     (yy - 1*(ly - 1) / 2.);
    return (0.5 - 0.5 * tanh(2 * (sqrt(rr2)-radius) / 4.))*(0.5 + 0.5 * tanh(2 * (xx-2.35*(lx - 1)/4.0) / 4.));//(0.5 - 0.5 * tanh(2 * (sqrt(rr2)-radius) / 4.))*(0.5 + 0.5 * tanh(2 * (xx-2.5*(lx - 1)/4.0) / 4.));//(1-initFluid2(k)-initFluid1(k))*(0.5 - 0.5 * tanh(2 * (sqrt(rr2)-radius) / 4.));//-(1-initFluid2(k)-initFluid1(k))*(0.5 - 0.5 * tanh(2 * (sqrt(rr2)-radius) / 4.)));
    //return (0.5 - 0.5 * tanh(2 * (yy-(ly + 1)/2.0) / 4.))*(0.5 + 0.5 * tanh(2 * (xx-(lx + 1)/2.0) / 4.))*(0.5 + 0.5 * tanh(2 * (yy) / 4.))*(0.5 - 0.5 * tanh(2 * (xx-(lx)) / 4.))+(0.5 + 0.5 * tanh(2 * (yy-(ly)) / 4.))*(0.5 + 0.5 * tanh(2 * (xx-(lx + 1)/2.0) / 4.)) + (0.5 - 0.5 * tanh(2 * (xx) / 4.))*(0.5 - 0.5 * tanh(2 * (yy-(ly + 1)/2.0) / 4.))*(0.5 + 0.5 * tanh(2 * (yy) / 4.));
    //return 1-initFluid2(k)-initFluid1(k);
}
double initFluid2(const int k) {
    double xx = computeXGlobal<Lattice>(k);
    double yy = computeY(ly, lz, k);

    int zz = computeZ(ly, lz, k);
    double rr2 = (xx - 2*(lx - 1) / 4.) *
                     (xx - 2*(lx - 1) / 4.) +
                     (yy - 1*(ly - 1) / 2.) *
                     (yy - 1*(ly - 1) / 2.);
                 //(yy - 8) *
                 //    (yy - 8);
    //return (xx<lx/2.)*(1-(topextend>bottomextend ? ((shape>topextend) ? shape : topextend) : ((shape>bottomextend) ? shape : bottomextend)))*(0.5 + 0.5 * tanh(2 * (hh) / 4.));
    //return (0.5 + 0.5 * tanh(2 * (yy-8) / 4.))*(0.5 - 0.5 * tanh(2 * (sqrt(rr2)-radius) / 4.));
    //return ((0.5 + 0.5 * tanh(2 * (xx-(lx - 1)/2.0) / 4.))+(0.5 - 0.5 * tanh(2 * (xx-(lx)) / 4.))+(0.5 - 0.5 * tanh(2 * (xx) / 4.)));

    //return (0.5 - 0.5 * tanh(2 * (sqrt(rr2)-radius) / 4.));//*(0.5 - 0.5 * tanh(2 * (xx-2*(lx - 1)/4.0) / 4.));//*(0.5 - 0.5 * tanh(2 * (xx-(lx - 1)/2.0) / 4.));
    return (1-initFluid1(k)-initFluid3(k))*(0.5 - 0.5 * tanh(2 * (sqrt(rr2)-radius) / 4.));// (0.5 + 0.5 * tanh(2 * (yy-(ly - 1)/2.0) / 4.))* (0.5 - 0.5 * tanh(2 * (yy-(ly)) / 4.))+(0.5 - 0.5 * tanh(2 * (yy) / 4.));
    
    //return (0.5 - 0.5 * tanh(2 * (yy-(ly + 1)/2.0) / 4.))*(0.5 - 0.5 * tanh(2 * (xx-(lx + 1)/2.0) / 4.))*(0.5 + 0.5 * tanh(2 * (yy) / 4.))*(0.5 + 0.5 * tanh(2 * (xx) / 4.))+(0.5 + 0.5 * tanh(2 * (yy-(ly)) / 4.))*(0.5 - 0.5 * tanh(2 * (xx-(lx + 1)/2.0) / 4.))+(0.5 + 0.5 * tanh(2 * (xx-(lx)) / 4.))*(0.5 - 0.5 * tanh(2 * (yy-(ly + 1)/2.0) / 4.))*(0.5 + 0.5 * tanh(2 * (yy) / 4.));
    //return (1.0-initFluid1(k));/*(hh<0);

}




void initParams(std::string inputfile) {
    params.addParameter<int>(lx, "lx");
    params.addParameter<int>(ly, "ly");
    params.addParameter<int>(lz, "lz");
    params.addParameter<int>(timesteps, "timesteps");
    params.addParameter<int>(saveInterval, "saveInterval");
    params.addParameter<double>(s01, "s01");
    params.addParameter<double>(s02, "s02");
    params.addParameter<double>(s03, "s03");
    params.addParameter<double>(s12, "s12");
    params.addParameter<double>(s13, "s13");
    params.addParameter<double>(s23, "s23");
    params.addParameter<double>(dens1, "dens1");
    params.addParameter<double>(dens2, "dens2");
    params.addParameter<double>(dens3, "dens3");
    params.addParameter<double>(dens4, "dens4");
    params.addParameter<std::string>(datadir, "datadir");
    params.addParameter<double>(tau1, "tau1");
    params.addParameter<double>(tau2, "tau2");
    params.addParameter<double>(tau3, "tau3");
    params.addParameter<double>(tau4, "tau4");
    params.addParameter<double>(Lambda1, "Lambda1");
    params.addParameter<double>(Lambda2, "Lambda2");
    params.addParameter<double>(Lambda3, "Lambda3");
    params.addParameter<double>(Lambda4, "Lambda4");
    params.addParameter<double>(radius, "radius");

    params.readInput(inputfile);

    Lattice::init(lx, ly, lz);
}


template <
    typename TTrait = DefaultTraitPressureWellBalancedN2<4,Lattice>>//::template SetBoundary<BounceBack>>//
auto initPressure() {
    
    FlowFieldPressureWellBalanced2<Lattice, TTrait> pressure;

    //pressure.template getBoundary<BounceBack>().setNodeID({1});

    pressure.setTaus(tau1,tau2,tau3,tau4);

    pressure.setDensities(dens1,dens2, dens3, dens4);

    pressure.setCollideID({0});
    return pressure;
}

// template<typename TTrait = typename DefaultTraitBinaryLeeHumidity<Lattice>:: template
// SetDataType<DataOldNewEquilibrium>>
template <int N, typename modeltype>
void initCH(modeltype& model) {
    model.setCollideID({0});
}


template <
    typename TTrait = typename DefaultTraitWellBalancedCH<0,4,Lattice>::template SetDataType<
        DataOldNewEquilibrium>::template AddProcessor<ChemicalPotentialCalculator4ComponentBoyer>>//::
        //template SetBoundary<BounceBack>>//:: template AddProcessor<std::tuple<GradientsMultiStencil<ChemicalPotential<0>,CentralQBounceBack,CentralXYZBounceBack>,GradientsMultiStencil<ChemicalPotential<1>,CentralQBounceBack,CentralXYZBounceBack>,GradientsMultiStencil<ChemicalPotential<2>,CentralQBounceBack,CentralXYZBounceBack>,GradientsMultiStencil<ChemicalPotential<3>,CentralQBounceBack,CentralXYZBounceBack>>>>  
auto initCH1() {
    WellBalancedCH<0, 4, Lattice, TTrait> ternary;

    initCH<0>(ternary);
    ternary.template getProcessor<ChemicalPotentialCalculator4ComponentBoyer>().setSurfaceTension(0,1,s01);
    ternary.template getProcessor<ChemicalPotentialCalculator4ComponentBoyer>().setSurfaceTension(0,2,s02);
    ternary.template getProcessor<ChemicalPotentialCalculator4ComponentBoyer>().setSurfaceTension(1,2,s12);
    ternary.template getProcessor<ChemicalPotentialCalculator4ComponentBoyer>().setSurfaceTension(0,3,s03);
    ternary.template getProcessor<ChemicalPotentialCalculator4ComponentBoyer>().setSurfaceTension(1,3,s13);
    ternary.template getProcessor<ChemicalPotentialCalculator4ComponentBoyer>().setSurfaceTension(2,3,s23);
    ternary.template getProcessor<ChemicalPotentialCalculator4ComponentBoyer>().setLambda(Lambda1,Lambda2,Lambda3,Lambda4);
    ternary.template getProcessor<ChemicalPotentialCalculator4ComponentBoyer>().setD(4);

    //ternary.template getForce<WellBalancedCHSource<0>>().setMij({-148000./907./50.,66000./907./50.,92000./907./50.,-10000./907./50.});
    //ternary.template getForce<WellBalancedCHSource<0>>().setMij({0,0,0,0});
    //ternary.template getForce<WellBalancedCHSource<0>>().setMij({-3,1,1,1});

    ternary.setMij({148000./907./50./(1./3.)/3.,-66000./907./50./(1./3.)/3.,-92000./907./50./(1./3.)/3.,10000./907./50./(1./3.)/3.});
    //ternary.setMij({3,-1,-1,-1});
    //ternary.setMij({0,0,0,0});

    ternary.template getBoundary<BounceBack>().setNodeID({1});

    return ternary;
}

template <
    typename TTrait = DefaultTraitWellBalancedCH<1,4,Lattice>>  
auto initCH2() {
    WellBalancedCH<1, 4, Lattice, TTrait> ternary;
    initCH<1>(ternary);

    //ternary.template getBoundary<BounceBack>().setNodeID({1});
    //ternary.template getForce<WellBalancedCHSource<1>>().setMij({66000./907./50., -152000./907./50.,    8000./907./50.,   78000./907./50.});
    //ternary.template getForce<WellBalancedCHSource<1>>().setMij({0,0,0,0});
    //ternary.template getForce<WellBalancedCHSource<1>>().setMij({1,-3,1,1});

    ternary.setMij({-66000./907./50./(1./3.)/3., 152000./907./50./(1./3.)/3.,    -8000./907./50./(1./3.)/3.,   -78000./907./50./(1./3.)/3.});
    //ternary.setMij({-1,3,-1,-1});
    //ternary.setMij({0,0,0,0});

    return ternary;
}

template <
    typename TTrait = DefaultTraitWellBalancedCH<2,4,Lattice>>  
auto initCH3() {
    WellBalancedCH<2, 4, Lattice, TTrait> ternary;
    initCH<2>(ternary);
    
    //ternary.template getBoundary<BounceBack>().setNodeID({1});
    //ternary.template getForce<WellBalancedCHSource<2>>().setMij({92000./907./50.,    8000./907./50., -167500./907./50.,   67500./907./50.});
    //ternary.template getForce<WellBalancedCHSource<2>>().setMij({0,0,0,0});
    //ternary.template getForce<WellBalancedCHSource<2>>().setMij({1,1,-3,1});

    ternary.setMij({-92000./907./50./(1./3.)/3.,    -8000./907./50./(1./3.)/3., 167500./907./50./(1./3.)/3.,   -67500./907./50./(1./3.)/3.});
    //ternary.setMij({-92000./907./50./(1./3.),    -8000./907./50./(1./3.), 167500./907./50./(1./3.),   -67500./907./50./(1./3.)});
    //ternary.setMij({-1,-1,3,-1});
    //ternary.setMij({0,0,0,0});


    return ternary;
}
