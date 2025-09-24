#include <math.h>
#include <stdlib.h>

#include <lbm.hh>

// This script simulates a four component layered poiseuille flow setup.
// You can modify the densities and relaxation times via the 'setDensities' and 'setTaus' functions.
// You can modify the body force magnitude in the setMagnitudeX function

int lx = 200;             // Size of domain in x direction
int ly = 200;             // Size of domain in y direction
int outflowoffset = 0;
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
double dens1 = 1;
double dens2 = 1;
#ifdef KINETIC
double evaporationrate = 0.0001;
#else
double diffusivity = 0.008;
double Hsat = 0.3;
double Hwall = 0.0;
#endif
std::string datadir = "data/";
double inflowmomentum = 0;
double pDiff = 0;
int equilibriumtimesteps = 0;
int postwidth = lx / 4;
int postheight = 96;
int capwidth = 0;
int capthickness = 0;
int factor = 1;
double tau1 = 1;
double tau2 = 1;

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
    int yy = computeY(ly, 1, k);
    double rr2 = (xx - (lx - 1) / 2. - offsetx) * (xx - (lx - 1) / 2. - offsetx) +
                 (yy - (ly - 1) / 2. - offsety) * (yy - (ly - 1) / 2. - offsety);

    if (yy <= 1*factor) return 1;
    if (yy > ly - 4*factor) return 2;
    if (yy == ly - 4*factor) return 22;
    if (yy >= ly - 4*factor - outflowoffset) return 123;

    if ((xx < postwidth / 2 || xx > lx - postwidth / 2) && yy < postheight) return 1;
    if ((xx < postwidth / 2 + capwidth / 2 || xx > lx - postwidth / 2 - capwidth / 2) && yy < postheight && yy >= postheight - capthickness) return 1;
    if (yy < postheight + offsety) return 5;

    return 0;
}

double initFluid(const int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    double rr2 = (xx - (lx - 1) / 2. - offsetx) * (xx - (lx - 1) / 2. - offsetx) +
                 (yy - postheight - offsety +40) * (yy - postheight - offsety + 40);

    return 0.5 - 0.5 * tanh(2 * ((yy - postheight - offsety)) / (sqrt(8 * kappa / A)));
    //return 0.5 - 0.5 * tanh(2 * ((sqrt(rr2)-40)) / (sqrt(8 * kappa / A)));
}

#ifndef KINETIC
double initHumidity(int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    double rr2 = (xx - (lx - 1) / 2. - offsetx) * (xx - (lx - 1) / 2. - offsetx) +
                 (yy - (ly - 1) / 2. - offsety) * (yy - (ly - 1) / 2. - offsety);
    if (yy < postheight + offsety)
        return Hsat;
    else {
        
        double t = (yy-postheight - offsety) / (ly - 4*factor - outflowoffset - postheight - offsety);
        double rhoIn = Hsat;
        double rhoOut = 0;
        return rhoIn*(1-t) + rhoOut*t;
    }
}
#endif

bool interfaceCondition(const double& val, int k) { return val < 0.5; }

bool interfaceConditionK(int k) { return OrderParameter<>::get<Lattice>(k) < 0.5; }

// Function used to define the fluid
// Here we set a tanh transition in the y direction. This should match the equilibrium profile.

double distancefunc(int k, int idx) {
    using Stencil = std::conditional_t<Lattice::NDIM == 1, D1Q3,
                           std::conditional_t<Lattice::NDIM == 2, D2Q9, D3Q19>>;;

    double normaldist = -0.5 * sqrt(8 * kappa / A) * atanh((OrderParameter<>::get<Lattice>(k) - 0.5) / 0.5);
    double* gradorderparam = GradientOrderParameter<>::getAddress<Lattice, Lattice::NDIM>(k, 0);

    std::vector<double> normal;
    if constexpr (Lattice::NDIM == 1)
        normal = {-gradorderparam[0] / sqrt(pow(gradorderparam[0], 2))};
    else if constexpr (Lattice::NDIM == 2)
        normal = {-gradorderparam[0] / sqrt(pow(gradorderparam[0], 2) + pow(gradorderparam[1], 2)),
                  -gradorderparam[1] / sqrt(pow(gradorderparam[0], 2) + pow(gradorderparam[1], 2))};
    else
        normal = {-gradorderparam[0] /
                      sqrt(pow(gradorderparam[0], 2) + pow(gradorderparam[1], 2) + pow(gradorderparam[2], 2)),
                  -gradorderparam[1] /
                      sqrt(pow(gradorderparam[0], 2) + pow(gradorderparam[1], 2) + pow(gradorderparam[2], 2)),
                  -gradorderparam[2] /
                      sqrt(pow(gradorderparam[0], 2) + pow(gradorderparam[1], 2) + pow(gradorderparam[2], 2))};
    double normdotci = 0;
    double magci = 0;
    for (int xyz = 0; xyz < Lattice::NDIM; xyz++) {
        normdotci += normal[xyz] * Stencil::Ci_xyz(xyz)[idx];
        magci += Stencil::Ci_xyz(xyz)[idx] * Stencil::Ci_xyz(xyz)[idx];
    }
    double dist = fabs(normaldist / normdotci);

    if (idx == 0 || std::isnan(dist) || std::isinf(dist)) {
        return 0.5;
    } else if (dist <= 0) {
        return 0.5;
    }
    return dist;
}

void initParams(std::string inputfile) {
    params.addParameter<int>(lx, "lx");
    params.addParameter<int>(ly, "ly");
    params.addParameter<int>(timesteps, "timesteps");
    params.addParameter<int>(saveInterval, "saveInterval");
    params.addParameter<double>(radius, "radius");
    params.addParameter<double>(theta, "theta");
    params.addParameter<double>(A, "A");
    params.addParameter<double>(kappa, "kappa");
    params.addParameter<double>(surfacetension, "surfacetension");
    params.addParameter<double>(interfacewidth, "interfacewidth");
    params.addParameter<int>(factor, "factor");
#ifdef KINETIC
    params.addParameter<double>(evaporationrate, "evaporationrate");
#else
    params.addParameter<double>(diffusivity, "diffusivity");
    params.addParameter<double>(Hsat, "Hsat");
    params.addParameter<double>(Hwall, "Hwall");
    diffusivity*=factor;
#endif
    params.addParameter<double>(offsetx, "offsetx");
    params.addParameter<int>(outflowoffset, "outflowoffset");
    params.addParameter<double>(offsety, "offsety");
    params.addParameter<double>(dens1, "dens1");
    params.addParameter<double>(dens2, "dens2");
    params.addParameter<double>(tau1, "tau1");
    params.addParameter<double>(tau2, "tau2");
    params.addParameter<std::string>(datadir, "datadir");
    params.addParameter<int>(equilibriumtimesteps, "equilibriumtimesteps");
    params.addParameter<int>(postwidth, "postwidth");
    params.addParameter<int>(postheight, "postheight");
    params.addParameter<int>(capwidth, "capwidth");
    params.addParameter<int>(capthickness, "capthickness");
    
    

    params.readInput(inputfile);
#ifndef KINETIC
    diffusivity*=factor*factor;
#endif
    postwidth *= factor;
    postheight *= factor;
    capwidth *= factor;
    capthickness *= factor;
    lx *= factor;
    ly *= factor;
    offsetx *= factor;
    offsety *= factor;
    
    A = 12 * surfacetension / interfacewidth;
    kappa = pow(interfacewidth, 2) / 8.0 * A;

    Lattice::init(lx, ly, 1);
}

#ifndef KINETIC
using traithumid =
    DefaultTraitHumidity<Lattice>::template SetBoundary<InterpolatedDirichlet, DirichletTemp, Convective,
                                                                          Refill<Humidity<>>, FreeSlip>::SetCollisionOperator<TRT>;;

template <typename TTrait = typename traithumid::SetDataType<DataOldNewEquilibrium>>
EvaporationHumidity<Lattice, TTrait> initHumidity() {
    EvaporationHumidity<Lattice, TTrait> humidity;

    using dbtype = InterpolatedDirichlet;

    humidity.setCollideID({0, 6, 8, 9});

    humidity.template getBoundary<FreeSlip>().setNodeID(1);

    humidity.template getBoundary<InterpolatedDirichlet>().setInterfaceDistanceFunction(distancefunc);
    humidity.template getBoundary<dbtype>().setNodeID(5);
    humidity.template getBoundary<dbtype>().setInterfaceVal(Hsat);

    humidity.template getBoundary<DirichletTemp>().setNodeID({2, 22, 3, 123});
    humidity.template getBoundary<DirichletTemp>().setInterfaceVal(Hwall);

    humidity.template getBoundary<Convective>().setNodeID({4, 7});

    humidity.setDiffusivity(diffusivity);

    humidity.template getProcessor<HumidityBoundaryLabels, 1>().setInterfaceCondition(interfaceCondition);
    humidity.template getProcessor<SetHumidityLiquid, 1>().setInterfaceVal(Hsat);

    humidity.template getProcessor<Gradients<Humidity<>, CentralXYZInterfaceMirrorSolid>>().setInterfaceDistance(
        distancefunc);
    humidity.template getProcessor<Gradients<Humidity<>, CentralXYZInterfaceMirrorSolid>>().setInterfaceVal(Hsat);
    humidity.template getProcessor<Gradients<Humidity<>, CentralXYZInterfaceMirrorSolid>>().setBoundaryID({1, 3, 4, 5});

    humidity.template getForce<EvaporationHumiditySource<EvaporationSourceMethod>>().setInterfaceHumidity(Hsat);
    humidity.template getForce<EvaporationHumiditySource<EvaporationSourceMethod>>().setGasDensity(dens2);

    return humidity;
}
#endif

template <class TLattice>
using DefaultTraitPressureHumidityInflow = typename DefaultTraitPressureWellBalanced<TLattice>::template SetBoundary<BounceBack,PressureOutflow2>::template AddForce<EvaporationPressureSource<EvaporationSourceMethod>>::SetCollisionOperator<MRT>;//, PressureOutflow

template <typename TTrait =
              typename DefaultTraitPressureHumidityInflow<Lattice>::template SetDataType<DataOldNewEquilibrium> >
              auto initPressure() {
    //PressureLeeHumidity<Lattice, TTrait> pressure;
    FlowFieldPressureWellBalanced3Evap<Lattice, TTrait> pressure;

    pressure.setCollideID({0, 8, 9, 5, 6, 123});
#ifdef KINETIC
    //pressure.template getForce<EvaporationPressureSource<EvaporationSourceMethod>>().setInterfaceHumidity(0);
    //pressure.template getForce<EvaporationPressureSource<EvaporationSourceMethod>>().setGasDensity(1);
#else
    //pressure.template getForce<EvaporationPressureSource<EvaporationSourceMethod>>().setInterfaceHumidity(Hsat);
    //pressure.template getForce<EvaporationPressureSource<EvaporationSourceMethod>>().setGasDensity(dens2);
#endif
    
    pressure.template getBoundary<BounceBack>().setNodeID({1});
    pressure.template getBoundary<PressureOutflow2>().setNodeID({2,22});

    pressure.setTaus(tau1,tau2);
    pressure.setDensities(dens1,dens2);
    pressure
        .template getProcessor<GradientsMultiStencil<Density<>, CentralXYZBounceBack, CentralQBounceBack>>()
        .setBoundaryID({1, 3, 4});

    pressure
        .template getProcessor<GradientsDirectional<Velocity<>, CentralXYZBounceBackDirectional>>()
        .setBoundaryID({1, 3, 4});

    return pressure;
}

using orderparamgradients = GradientsMultiStencil<OrderParameter<0>,CentralXYZBounceBack,LaplacianCentralWetting/*,BiasedQBounceBack,BiasedXYZBounceBack*/,CentralQBounceBack>;
/*template <class TLattice>
using DefaultTraitBinaryLeeHumidityInflow = typename DefaultTraitBinaryLee<TLattice>::
    template SetProcessor<std::tuple<ConstantGradientBoundary<OrderParameter<>>,ConstantGradientBoundary<Pressure<>>,ConstantGradientBoundary<Density<>>,ConstantGradientBoundary<ChemicalPotential<>>>,std::tuple<orderparamgradients, ChemicalPotentialCalculatorBinaryLee>,std::tuple<
        GradientsMultiStencil<ChemicalPotential<>, LaplacianCentralBounceBack>,
        GradientsMultiStencil<Pressure<>, CentralXYZBounceBack, CentralQBounceBack, MixedXYZBounceBack,
                              MixedQBounceBack>,
#ifdef KINETIC
        SimpleMassLossCalculator>>
#else
        MassLossCalculatorInterpolated>>
#endif
        ::template AddForce<EvaporationPhaseSource<EvaporationSourceMethod>>::
        template SetBoundary<BounceBack,Convective>;*/

/*template <class TZerothMoment, class TZerothMomentOld, class TMethod = AdvectionDiffusionSourceMethod>
class CHCorrection : public ForceBase<TMethod> {
   public:
    template <class TTraits>
    inline double computeXYZ(int xyz, int k) const {
        return (OrderParameter::template get<typename TTraits::Lattice>(k) *
                    Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz) -
                OrderParameterOld::template get<typename TTraits::Lattice>(k) *
                    VelocityOld<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz)) /
               TTraits::Lattice::DT;
    }
};*/

template<int N, int Nmax, class TLattice>
using TraitWellBalancedCHHumidityInflow = typename DefaultTrait<TLattice,Nmax> :: template SetBoundary<BounceBack>
/* Change this so only done for N=0 except for order param*/   :: template SetProcessor<std::tuple<ConstantGradientBoundary<OrderParameter<>>,ConstantGradientBoundary<Pressure<>>,ConstantGradientBoundary<Density<>>,ConstantGradientBoundary<ChemicalPotential<>>,MirrorBoundary<Velocity<>,Lattice::NDIM>>,std::tuple<orderparamgradients,ChemicalPotentialCalculatorBinaryLee>,std::tuple<MassLossCalculatorInterpolated,GradientsMultiStencil<ChemicalPotential<>, CentralQBounceBack, CentralXYZBounceBack>,Gradients<ChemicalPotential2<>, CentralQFourthOrderBounceBack>>>
                                                               :: template SetForce< Lin2CompCHSource,EvaporationPhaseSource<EvaporationSourceMethod>>::template SetBoundary<BounceBack,Convective>;

template <typename TTrait =
              typename TraitWellBalancedCHHumidityInflow<0,2,Lattice>::template SetDataType<DataOldNewEquilibrium>>
auto initBinary() {
    //BinaryLeeHumidity<Lattice, TTrait> binary;
    WellBalancedCH<0,2,Lattice,TTrait> binary;
    binary.setCollideID({0, 5, 6, 8, 9, 123});

    binary.template getProcessor<ConstantGradientBoundary<OrderParameter<>>>().setNodeID({2,22});
    binary.template getProcessor<ConstantGradientBoundary<Density<>>>().setNodeID({2,22});
    binary.template getProcessor<ConstantGradientBoundary<Pressure<>>>().setNodeID({2,22});
    binary.template getProcessor<MirrorBoundary<Velocity<>,Lattice::NDIM>>().setNodeID({2,22});
    binary.template getProcessor<ConstantGradientBoundary<ChemicalPotential<>>>().setNodeID({2,22});

    binary.template getProcessor<orderparamgradients>().setBoundaryID({1,2,22,3, 4});
    binary.template getProcessor<GradientsMultiStencil<ChemicalPotential<>, CentralQBounceBack,CentralXYZBounceBack>>().setBoundaryID({1,2,22,3, 4});
    binary.template getProcessor<Gradients<ChemicalPotential2<>, CentralQFourthOrderBounceBack>>().setBoundaryID({1,2,22,3, 4});

    binary.setAij({-0.0,0});
    //binary.setTau(0.8);

    //binary.template getForce<MuSourceLocal>().setBeta(A);
    //binary.template getForce<MuSourceNonLocal>().setBeta(A);

    binary.template getProcessor<ChemicalPotentialCalculatorBinaryLee>().setA(A);
    binary.template getProcessor<ChemicalPotentialCalculatorBinaryLee>().setKappa(kappa);
    binary.template getProcessor<ChemicalPotentialCalculatorBinaryLee>().setOmega(2.5);
#ifdef KINETIC
    binary.template getProcessor<SimpleMassLossCalculator>().setEvaporationRate(0.0);
#else
    binary.template getProcessor<MassLossCalculatorInterpolated>().toggleCalculate(false);
    binary.template getProcessor<MassLossCalculatorInterpolated>().setInterfaceHumidity(Hsat);
    binary.template getProcessor<MassLossCalculatorInterpolated>().setDiffusivity(diffusivity);
    binary.template getProcessor<MassLossCalculatorInterpolated>().setInterfaceWidth(sqrt(8 * kappa / A));
    binary.template getProcessor<MassLossCalculatorInterpolated>().setPhiGasLiquid(0, 1);
    binary.template getProcessor<MassLossCalculatorInterpolated>().setGasDensity(dens2);
#endif
    binary.template getBoundary<BounceBack>().setNodeID({1, 2,22, 3,4});
    //binary.template getBoundary<Dirichlet>().setNodeID({44});
    //binary.template getBoundary<Dirichlet>().setInterfaceVal(0);
    binary.template getBoundary<Convective>().setNodeID({222});
    double wettingprefactor = -2 * cos(theta * M_PI / 180.0) * sqrt(2 * A / kappa);

    binary.template getProcessor<orderparamgradients>().setWettingPrefactor(wettingprefactor);

    return binary;
}

template <typename T>
void AfterEquilibration(int eqsteps, T& model) {
#ifdef KINETIC
    if (eqsteps == equilibriumtimesteps)
        model.template getProcessor<SimpleMassLossCalculator>().setEvaporationRate(evaporationrate);
#else
    if (eqsteps == equilibriumtimesteps)
        model.template getProcessor<MassLossCalculatorInterpolated>().toggleCalculate(true);
#endif
}
