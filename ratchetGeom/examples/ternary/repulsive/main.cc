#include <chrono>
#include <cstdlib>
#include <thread>
#include <math.h>
#include <stdlib.h>

#include <lbm.hh>

// This script simulates a four component layered poiseuille flow setup.
// You can modify the densities and relaxation times via the 'setDensities' and 'setTaus' functions.
// You can modify the body force magnitude in the setMagnitudeX function

int lx = 400;  // Size of domain in x direction
int ly = 100;  // Size of domain in y direction
int timesteps = 5000;     // Number of iterations to perform
int saveInterval = 1000;  // Interval to save global data
double radius1 = 30;
double radius2 = 30;
double spacing = 70;
std::string datadir = "data/";
double s01 = 0.001;
double s02 = 0.001;
double s12 = 0.001;
double gamma0 = s01 + s02 - s12;
double gamma1 = s01 + s12 - s02;
double gamma2 = s02 + s12 - s01;
double interfacewidth = 4;
double tau0 = 1.0;
double tau1 = 1.0;
double tau2 = 1.0;
std::vector<double> gammas = {};

InputParameters params;

using Lattice = LatticePropertiesRuntime<NoParallel, 2>;

int initSolid(const int k) {
    int y = computeY(ly, 1, k);
    //return 0;
    if (y <= 1) return 1;
    if (y >= ly-2) return 1;
    else return 0;
}

double initFluid1(const int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    double rr2 = (xx - (lx - 1) / 2. - spacing/2.) * (xx - (lx - 1) / 2. - spacing/2.) + (yy - (ly - 1) / 2.) * (yy - (ly - 1) / 2.);
    if (rr2 < radius1 * radius1) {
        return 1;
    } else
        return 0;
}

double initFluid2(const int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    double rr2 = (xx - (lx - 1)/2. + spacing/2.) * (xx - (lx - 1)/2. + spacing/2.) + (yy - (ly - 1)/2.) * (yy - (ly - 1)/2.);

    if (rr2 < radius2 * radius2) {
        return 1;
    } else
        return 0;
}

void initParams(std::string inputfile) {
    params.addParameter<int>(lx, "lx");
    params.addParameter<int>(ly, "ly");
    params.addParameter<int>(timesteps, "timesteps");
    params.addParameter<int>(saveInterval, "saveInterval");
    params.addParameter<double>(radius1, "radius1");
    params.addParameter<double>(spacing, "spacing");
    params.addParameter<double>(s01, "surfacetension01");
    params.addParameter<double>(s02, "surfacetension02");
    params.addParameter<double>(s12, "surfacetension12");
    params.addParameter<double>(radius2, "radius2");
    params.addParameter<std::string>(datadir, "datadir");
    params.addParameter<double>(tau0, "tau0");
    params.addParameter<double>(tau1, "tau1");
    params.addParameter<double>(tau2, "tau2");
    params.addParameter<double>(interfacewidth, "interfacewidth");

    params.readInput(inputfile);

    gamma0 = s01 + s02 - s12;
    gamma1 = s01 + s12 - s02;
    gamma2 = s02 + s12 - s01;
    gammas = {gamma0, gamma1, gamma2};
    Lattice::init(lx, ly, 1);
}

//
///////// Simulation details
//


// Function to create the pressure LBM model (solves navier stokes)
auto initPressure() {

    // PressureLee is the model class, and we give the lattice and traits as template parameters.
    FlowFieldPressureWellBalanced<Lattice, DefaultTraitPressureWellBalancedN<3,Lattice>::template AddForce<RepulsiveForce<0,1>,RepulsiveForce<1,0>>> pressure;

    // Boundary ids to apply the LBM model
    pressure.setCollideID({0});
    pressure.template getBoundary<BounceBack>().setNodeID({1});

    // Set relaxation times of each component
    pressure.setTaus(tau0,tau1,tau2);
    pressure.setDensities(1.0,1.0,1.0);

    pressure.template getForce<RepulsiveForce<0,1>>().setMagnitude(0.00007);
    pressure.template getForce<RepulsiveForce<1,0>>().setMagnitude(0.00007);

    // Return the model so it can be used in main.cc
    return pressure;
}

// Function to create the order parameter LBM model (solves cahn hilliard)
template<int N>
auto initTernary() {

    using Trait = typename std::conditional<
        N == 0, typename DefaultTraitWellBalancedCH<N,3,Lattice>:: template SetProcessor<std::tuple<GradientsMultiStencil<OrderParameter<N>,CentralXYZBounceBack,LaplacianCentralWetting/*,BiasedQBounceBack,BiasedXYZBounceBack*/,CentralQBounceBack>>,std::tuple<ChemicalPotentialCalculatorTernaryLee>,std::tuple<GradientsMultiStencil<ChemicalPotential<0>,CentralXYZBounceBack,CentralQBounceBack>,GradientsMultiStencil<ChemicalPotential<1>,CentralXYZBounceBack,CentralQBounceBack>,GradientsMultiStencil<ChemicalPotential<2>,CentralXYZBounceBack,CentralQBounceBack>>>,
                DefaultTraitWellBalancedCH<N,3,Lattice>>::type;

    // WellBalancedCH is the model class, and we give the lattice and traits as template parameters.
    WellBalancedCH<N, 3, Lattice, Trait> ternary;

    // Boundary ids to apply the LBM model
    ternary.setCollideID({0});
    ternary.template getBoundary<BounceBack>().setNodeID({1});

    std::vector<double> mobility = {0,0,0};
    mobility[N]=-1.;

    ternary.setAij(mobility);

    // Chemical potential calculation needs the free energy parameters
    if constexpr (N == 0) {
        ternary.template getProcessor<ChemicalPotentialCalculatorTernaryLee>().setSurfaceTension(s01, s02, s12);
        ternary.template getProcessor<ChemicalPotentialCalculatorTernaryLee>().setInterfaceWidth(interfacewidth);
        ternary.template getProcessor<ChemicalPotentialCalculatorTernaryLee>().setLambda(0.01);
    }

    return ternary;
}

int main(int argc, char **argv) {

    initParams("input.txt");

    auto ternary1 = initTernary<0>();
    auto ternary2 = initTernary<1>();
    auto pressure = initPressure();

    SaveHandler<Lattice> saver(datadir);
    Geometry<Lattice>::initialiseBoundaries(initSolid,{0});
    
    OrderParameter<1>::set<Lattice, 1>(initFluid2);
    OrderParameter<1>::smooth<Lattice, 1>(4);
    OrderParameter<0>::set<Lattice, 1>(initFluid1);
    OrderParameter<0>::smooth<Lattice, 1>(4);

    Algorithm lbm(ternary2, ternary1, pressure);

    saver.saveHeader(timesteps, saveInterval);

    char fdump[256];
    sprintf(fdump, "%s/Header.mat", datadir.c_str());
    std::ofstream fs(fdump, std::ios::out | std::ios::binary | std::ios::app);
    fs.write((char*)(&s01), sizeof(double));
    fs.write((char*)(&s02), sizeof(double));
    fs.write((char*)(&s12), sizeof(double));
    fs.close();

    for (int timestep = 0; timestep <= timesteps; timestep++) {
        TIME = timestep;
        
        // Save the desired parameters, producing a binary file for each.
        if (timestep % saveInterval == 0) {
            if (mpi.rank == 0) std::cout << "Saving at timestep " << timestep << "." << std::endl;

            saver.saveBoundaries(timestep);
            saver.saveParameter<ChemicalPotential<0>>(timestep, true);
            saver.saveParameter<ChemicalPotential<1>>(timestep, true);
            saver.saveParameter<ChemicalPotential<2>>(timestep, true);
            saver.saveParameter<LaplacianChemicalPotential<0>>(timestep, true);
            saver.saveParameter<LaplacianChemicalPotential<1>>(timestep, true);
            saver.saveParameter<LaplacianChemicalPotential<2>>(timestep, true);
            saver.saveParameter<Pressure<>>(timestep);
            saver.saveParameter<OrderParameter<0>>(timestep, true);
            saver.saveParameter<OrderParameter<1>>(timestep, true);
            saver.saveParameter<Velocity<>, Lattice::NDIM>(timestep);
            saver.saveParameter<ForceRepulsive<>, Lattice::NDIM>(timestep);
        }

        lbm.evolve();
        if (ternary1.isNan()) {
            std::cout << "NaN encountered." << std::endl;
            exit(1);
        }
    }
}
