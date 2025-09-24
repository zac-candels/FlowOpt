#include "mainTernary.hh"
#include <cstdlib>
#include <chrono>
#include <thread>

int main(int argc, char **argv){

    using namespace std::this_thread; // sleep_for, sleep_until
    using namespace std::chrono; // nanoseconds, system_clock, seconds

    int seed;
    std::string seedstr;
    #ifdef MPIPARALLEL
    mpi.init();

    MPI_Type_create_resized(MPI_INT, 0L, sizeof(Boundary<Lattice::NDIM>), &mMPIBoundary);
    MPI_Type_commit(&mMPIBoundary);

    int blocklengths[3] = {1,1,Lattice::NDIM};
    MPI_Datatype types[3] = {MPI_INT,MPI_CXX_BOOL,MPI_INT8_T};
    std::array<int8_t,Lattice::NDIM> a;
    MPI_Aint offsets[3];
    offsets[0] = offsetof(Boundary<Lattice::NDIM>, Id);
    offsets[1] = offsetof(Boundary<Lattice::NDIM>, IsCorner);
    offsets[2] = offsetof(Boundary<Lattice::NDIM>, NormalDirection) + (size_t)((((char *)(&a)-(char *)(&a[0]))));// + offsetof(std::array<int8_t,Lattice::NDIM>, NormalDirection);
    
    MPI_Type_create_struct(3, blocklengths, offsets, types, &mMPIBoundary);
    MPI_Type_commit(&mMPIBoundary);
    #endif
    if (argc>1){
        seed=std::atoi(argv[1]);
        seedstr=std::to_string(seed);
        initParams("input/input"+seedstr+".txt");
    }
    else{
        seedstr="";
        initParams("input.txt");
    }

    if(mpi.rank==0&&argc==0){
        int ret;
        std::string tmp="rm input/input"+seedstr+".txt";
        const char *array = tmp.c_str();
        ret = std::system(array);
    }
    Lattice l;

    SaveHandler<Lattice> saver(datadir);

    /*OrderParameter<2>::set<Lattice>(initFluid1);
    OrderParameter<1>::set<Lattice>(initFluid2);
    OrderParameter<0>::set<Lattice>(initFluid3);*/
    OrderParameter<2>::set<Lattice>(initFluid3);
    OrderParameter<1>::set<Lattice>(initFluid1);
    OrderParameter<0>::set<Lattice>(initFluid2);
    
    Geometry<Lattice>::initialiseBoundaries(initBoundary,{0,5,6,7});

    auto pressure = initPressure<>();

    auto ch1 = initCH1<>();
    auto ch2 = initCH2<>();
    //auto ch3 = initCH3<>();
    
    Algorithm lbm(ch2,ch1,pressure);// ternary3,ch3,
    //Algorithm init();
    saver.saveHeader(timesteps, saveInterval);
    for (int timestep=0; timestep<=timesteps; timestep++) {
        TIME=timestep;
        // Save the desired parameters, producing a binary file for each.
        if (timestep%saveInterval==0) {
            if(mpi.rank==0)std::cout<<"Saving at timestep "<<timestep<<"."<<std::endl;
            
            saver.saveBoundaries(timestep);
            saver.saveParameter<ChemicalPotential<0>>(timestep,true);
            saver.saveParameter<ChemicalPotential<1>>(timestep,true);
            saver.saveParameter<ChemicalPotential<2>>(timestep,true);
            saver.saveParameter<ChemicalPotential<3>>(timestep,true);
            saver.saveParameter<LaplacianOrderParameter<>>(timestep);
            saver.saveParameter<LaplacianChemicalPotential<>>(timestep);
            saver.saveParameter<Pressure<>>(timestep);
            saver.saveParameter<Density<>>(timestep);
            saver.saveParameter<OrderParameter<0>>(timestep,true);
            saver.saveParameter<OrderParameter<1>>(timestep,true);
            saver.saveParameter<OrderParameter<2>>(timestep,true);
            saver.saveParameter<Velocity<>,Lattice::NDIM>(timestep);
            
        }

        // Evolve by one timestep
        lbm.evolve();
    }
    
}
