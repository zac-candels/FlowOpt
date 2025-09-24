#include <math.h>
#include <stdlib.h>

#include <lbm.hh>

int lx = 200;             // Simulation size in x direction
int ly = 200;             // Simulation size in y direction
int lz = 1;               // Simulation size in z direction
int timesteps = 5000;     // Total number of iterations
int saveInterval = 1000;  // How frequently to save order parameter, velocity etc
double radius = 25.0;     // Droplet radius
double theta = 90.0;      // Contact angle
double A = 0.003;         // Parameter for free energy functional (depends on surface tension and interface width, you can
                          // Ignore for now)
double kappa = A * 9 / 8; // Parameter for free energy functional (depends on surface tension and interface width, you can ignore for now)
double surfacetension = sqrt(2 * A * kappa) / 6;  // Surface tension in lattice units
double interfacewidth = sqrt(8 * kappa / A);      // Width of the diffuse interface between fluids
double posx = 0.0;        // Droplet position in the x direction
double posy = 0.0;        // Droplet position in the y direction
//double posz = 0.0;        // Droplet position in the z direction, uncomment for full 3D
double dens1 = 1;         // Droplet density
double dens2 = 1;         // Air density
double bodyforcex = 0;    // Magnitude of body force acting on the droplet in the x direction
std::string datadir = "data/"; // Data directory
int equilibriumtimesteps = 0; // Number of timesteps until body force is applied
double tau1 = 1;          // Relaxation time of droplet
double tau2 = 1;          // Relaxation time of droplet
int solidthickness = 1; // Thickness of the solid boundary on the bottom of the domain

// Class that will handle input parameters
InputParameters params;

// Lattice classes, contain information about the size of the simulation and the parallelisation.
#ifdef MPIPARALLEL //If you have -DM
using Lattice = LatticePropertiesRuntime<ParallelX<2> /*Class to handle parallelisation*/, 3 /*Number of physical dimensions*/>;
#else
using Lattice = LatticePropertiesRuntime<NoParallel, 3>;
#endif

// Function used to define the solid geometry
int initBoundary(const int k) {
    // x coordinate of the lattice node k
    //int xx = computeXGlobal<Lattice>(k);
    // y coordinate of the lattice node k
    int yy = computeY(ly, lz, k);
    // z coordinate of the lattice node k

    // Layer of nodes with id 1 at the bottom of the domain, we will later apply bounce back to these nodes
    if (yy <= solidthickness) return 1;

    // Layer of nodes with id 2 at the top of the domain, we will later apply the mirror boundary condition to these nodes
    if (yy >= ly - 2) return 1;

    // If none of these conditions apply, return 0 to signify a fluid node
    return 0;
}

//Initialises the fluid with a bulk value of 1 for the droplet and a bulk value of 0 for air
double initFluid(const int k) {
    // x coordinate of the lattice node k
    int xx = computeXGlobal<Lattice>(k);
    // y coordinate of the lattice node k
    int yy = computeY(ly, lz, k);
    // z coordinate of the lattice node k
    //int zz = computeZ(ly, lz, k); // Uncomment for 3D

    // Radial distance from the centre of the droplet at (posx, posy)
    double rr2 = (xx - posx) * (xx - posx) + (yy - posy) * (yy - posy);
    // Double rr2 = (xx - posx) * (xx - posx) + (yy - posy) * (yy - posy) + (zz - posz) * (zz - posz); // Switch these for 3D

    // Smooth droplet
    return (0.5 - 0.5 * tanh(2 * (sqrt(rr2) - radius) / (4.0))); // Cut off just above the posts at postheight
}

// Function to initialise the parameters from an input file
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
    params.addParameter<double>(posx, "posx");
    params.addParameter<double>(posy, "posy");
    //params.addParameter<double>(posz, "posz"); // Uncomment for 3D
    params.addParameter<double>(dens1, "dens1");
    params.addParameter<double>(dens2, "dens2");
    params.addParameter<double>(tau1, "tau1");
    params.addParameter<double>(tau2, "tau2");
    params.addParameter<std::string>(datadir, "datadir");
    params.addParameter<double>(bodyforcex, "bodyforcex");
    params.addParameter<int>(equilibriumtimesteps, "equilibriumtimesteps");

    /*
    If you want to add a parameter here, follow the format above
    params.addParameter< *parameter type* >( *parameter name in this file*, *parameter name in the input file* );
    e.g. I have added:*/
    params.addParameter<int>( solidthickness, "solidthickness" );

    // Read the input file and initialise the parameters
    params.readInput(inputfile);
    posy = solidthickness;

    // Initialise free energy parameters from the surface tension and interface width
    A = 12 * surfacetension / interfacewidth;
    kappa = pow(interfacewidth, 2) / 8.0 * A;

    // Initialise the lattice class with the simulation size
    Lattice::init(lx, ly, lz);

}

//
///////// Simulation details
//

// Processor that will calculate the gradient in density with the given stencils.
// E.g. CentralXYZBounceBack is central gradients in the x, y and z directions with bounce back reflection applied to the density if the gradient would go into a solid node.
using densitygradients =
    GradientsMultiStencil<Density<>, CentralXYZBounceBack, CentralQBounceBack, MixedXYZBounceBack, MixedQBounceBack>;
using velocitygradients =
    GradientsDirectional<Velocity<>, CentralXYZBounceBackDirectional>;

// This trait class contains anything that is not inherent to the lbm model that we want to simulate. We will give this to the model.
using NewTraitPressureLee =
    typename DefaultTrait<Lattice, 2> // Start from default trait class, will be D2Q9(2D) or D3Q19(3D) and BGK collision operator
    ::template SetBoundary<BounceBack> // Here we specify we want bounceback boundary conditions
    ::template SetForce<PressureLeeForce, PressureLeeForce2, BodyForce<>> // PressureLeeForce, PressureLeeForce2 are necessary for the model (can ignore this).
                                                                          // Bodyforce can be used to apply an external force, you can replace this with the LRP force once you have it.
    ::template SetProcessor<std::tuple<densitygradients,velocitygradients>>; // Calculate the gradient in density as above.
                                               // When we implement the viscous dissipation calculator, we will create a processor and add it here.

// Function to create the pressure LBM model (solves navier stokes)
auto initPressure() {

    // PressureLee is the model class, and we give the lattice and traits as template parameters.
    PressureLee<Lattice, NewTraitPressureLee> pressure;

    // Boundary ids to apply the LBM model
    pressure.setCollideID({0});

    // Apply the bounce back boundary condition on all nodes with id 1
    pressure.template getBoundary<BounceBack>().setNodeID(1);

    // Set magnitude of bodyforce and component to apply it to, will be zero until AfterEquilibration is called
    pressure.template getForce<BodyForce<>>().setForce({0, 0}, // 0 in the x direction, 0 in the y direction
                                                        0); // Act on 0 component (droplet)

    // Set relaxation times of each component
    pressure.setTau1(tau1);
    pressure.setTau2(tau2);

    // Set solid ids for the density gradient processor
    pressure.template getProcessor<densitygradients>().setBoundaryID({1});
    pressure.template getProcessor<velocitygradients>().setBoundaryID({1});

    // Return the model so it can be used in main.cc
    return pressure;
}

// Processor that will calculate the gradient in the order parameter with the given stencils.
// E.g. CentralXYZBounceBack is central gradients in the x, y and z directions with bounce back reflection applied to the order parameter if the gradient would go into a solid node.
// Note that here we have included a stencil to calculate the laplacian with the wetting boundary condition.
using orderparamgradients = GradientsMultiStencil<OrderParameter<>, CentralXYZBounceBack, CentralQBounceBack,
                                                  MixedXYZBounceBack, MixedQBounceBack, LaplacianCentralWetting>;

// Gradients in pressure
using pressuregradients =
    GradientsMultiStencil<Pressure<>, CentralXYZBounceBack, CentralQBounceBack, MixedXYZBounceBack, MixedQBounceBack>;

using NewTraitBinaryLee =
    typename DefaultTrait<Lattice, 2>
    ::template SetBoundary<BounceBack> // Again, bounce back and mirror boundary conditions
    ::template SetProcessor<std::tuple<orderparamgradients, pressuregradients, ChemicalPotentialCalculatorBinaryLee>> // Note the inclusion of a class to calculate the chemical potential
    ::template AddProcessor<std::tuple<GradientsMultiStencil<ChemicalPotential<>, LaplacianCentralBounceBack>>> // Model specific (although the std::tuple<...> means this will be calculated after the processors on the row above have been calculated for every lattice point)
    ::template SetForce<LeeBinarySource, MuSourceLocal, MuSourceNonLocal>; // Model specific

// Function to create the order parameter LBM model (solves cahn hilliard)
auto initBinary() {

    // BinaryLee is the model class, and we give the lattice and traits as template parameters.
    BinaryLee<Lattice, NewTraitBinaryLee> binary;

    // Boundary ids to apply the LBM model
    binary.setCollideID({0});

    // Set densities of each component
    binary.setDensity1(dens1);
    binary.setDensity2(dens2);

    // Set relaxation times of each component (This model also has to know this)
    binary.setTau1(tau1);
    binary.setTau2(tau2);

    // Set solid ids for the gradient processors
    binary.template getProcessor<GradientsMultiStencil<ChemicalPotential<>, LaplacianCentralBounceBack>>()
        .setBoundaryID({1});
    binary.template getProcessor<pressuregradients>().setBoundaryID({1});
    binary.template getProcessor<orderparamgradients>().setBoundaryID({1});

    // Model specific, needs the free energy parameters
    binary.template getForce<MuSourceLocal>().setBeta(A);
    binary.template getForce<MuSourceNonLocal>().setBeta(A);

    // Chemical potential calculation needs the free energy parameters
    binary.template getProcessor<ChemicalPotentialCalculatorBinaryLee>().setA(A);
    binary.template getProcessor<ChemicalPotentialCalculatorBinaryLee>().setKappa(kappa);
    // Stabilisation parameter to stop the order parameter going below 0
    binary.template getProcessor<ChemicalPotentialCalculatorBinaryLee>().setOmega(0.0005);

    // Apply the bounce back boundary condition on all nodes with id 1
    binary.template getBoundary<BounceBack>().setNodeID(1);

    // Pre factor for wetting that depends on the contact angle
    double wettingprefactor = -2 * cos(theta * M_PI / 180.0) * sqrt(2 * A / kappa);
    // Give this to the gradient calculator
    binary.template getProcessor<orderparamgradients>().setWettingPrefactor(wettingprefactor);

    return binary;
}

//Will apply the body force after the equilibrium timesteps
template <typename T>
void AfterEquilibration(int eqsteps, T& model) {
    if (eqsteps == equilibriumtimesteps) model.template getForce<BodyForce<>>().setForce({bodyforcex, 0}, // bodyforcex in the x direction, 0 in the y direction
                                                                                         0); // Act on 0 component (droplet)
}