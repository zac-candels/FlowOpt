#include "main.hh"
#include <cstdlib>
#include <chrono>
#include <thread>


template <template <int> typename TParameter, typename TLattice, int NDIM, typename T, T I, T... Is>
void initialiseInstance(int instance,typename TParameter<0>::ParamType val, int k, std::index_sequence<I, Is...>) {
    if (instance == I) {
        return TParameter<I>::template initialise<TLattice>(val, k);
    } else if constexpr (sizeof...(Is) > 0) {
        return getInstance<TParameter, TLattice, NDIM>(instance, std::index_sequence<Is...>{});
    } else {
        throw std::invalid_argument("Requested instance not in index list.");
    }
}

template <template <int> typename TParameter, int N, typename TLattice, int NDIM = 1>
void initialiseInstance(int instance,typename TParameter<0>::ParamType val, int k) {
    if (instance >= N) throw std::invalid_argument("Requested instance must be less than N.");
    return getInstance<TParameter, TLattice, NDIM>(instance, val, k, std::make_index_sequence<N>{});
}

template <int counter, typename... TModels>
struct getAlgorithmType {
    using type = typename getAlgorithmType<counter-1,decltype(initCH<counter>()),TModels...>::type;
};

template <typename... TModels>
struct getAlgorithmType<0,TModels...> {
    using type = Algorithm<TModels...,decltype(initCH<0>()),decltype(initPressure())>;
};

template <int counter,typename... TModels>
struct getModelTuple{
    using type = typename getModelTuple<counter-1,decltype(initCH<counter>()),TModels...>::type;

};

template <typename... TModels>
struct getModelTuple<0,TModels...>{
    using type = std::tuple<TModels...,decltype(initCH<0>()),decltype(initPressure())>;
};

struct getAlgorithm{
    typename getModelTuple<numcomponents-1>::type mt_Models;

    template<int counter,typename T,typename... Ts>
    typename getAlgorithmType<numcomponents-1>::type get(T pressure,Ts... models) {
        if constexpr (sizeof...(models) < numcomponents-1) {
            std::get<decltype(initCH<numcomponents-counter>())>(mt_Models) = initCH<numcomponents-counter>();
            return get<counter-1>(pressure,models...,initCH<numcomponents-counter>());
        } else {
            std::get<decltype(initCH<0>())>(mt_Models) = initCH<0>();
            return std::apply([](auto&... models) { return Algorithm(models...); }, mt_Models);
        }
    }

    typename getAlgorithmType<numcomponents-1>::type get() {
        std::get<decltype(initPressure())>(mt_Models) = initPressure();
        return get<numcomponents-1>(std::get<decltype(initPressure())>(mt_Models));
    }
};

template<template<int> typename TParameter, int N,typename TSaver>
void saverecursive(TSaver& saver, int timestep) {
    saver.template saveParameter<TParameter<N>>(timestep,true);
    if constexpr (N>0) {
        saverecursive<TParameter,N-1,TSaver>(saver,timestep);
    }
}

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

    std::vector<double> prop={prop1,prop2,prop3,prop4,prop5,prop6,prop7,prop8,prop9,prop10,prop11,prop12};
    for (int k : RangeK<Lattice>()) {
        std::vector<double> randomval(ncomp,0);
        double sum = 0;
        for (int n = 0; n<ncomp; n++) {
            double rand = prop[n]*ncomp*static_cast <float> (std::rand()) / static_cast <float> (RAND_MAX);
            sum+=rand;
            randomval[n] = rand;
        }
        for (int n = 0; n<ncomp-1; n++) {
            initialiseInstance<OrderParameter, numcomponents - 1, Lattice>(n,k,randomval[n]/sum);
        }
        for (int n = ncomp-1; n<numcomponents-1; n++) {
            initialiseInstance<OrderParameter, numcomponents - 1, Lattice>(n,k,0.0);
        }

    }

    /*saver.loadParameter<OrderParameter<>>(datadir+"/OrderParameter0_t1000000.mat","bin");
    saver.loadParameter<OrderParameter<1>>(datadir+"/OrderParameter1_t1000000.mat","bin");
    saver.loadParameter<OrderParameter<2>>(datadir+"/OrderParameter2_t1000000.mat","bin");
    saver.loadParameter<OrderParameter<3>>(datadir+"/OrderParameter3_t1000000.mat","bin");
    saver.loadParameter<OrderParameter<4>>(datadir+"/OrderParameter4_t1000000.mat","bin");
    //saver.loadParameter<OrderParameter<5>>(datadir+"/OrderParameter5_t310000.mat","bin");
    //saver.loadParameter<OrderParameter<6>>(datadir+"/OrderParameter6_t310000.mat","bin");
    //saver.loadParameter<OrderParameter<7>>(datadir+"/OrderParameter7_t120000.mat","bin");
    //saver.loadParameter<OrderParameter<8>>(datadir+"/OrderParameter8_t120000.mat","bin");
    //saver.loadParameter<OrderParameter<9>>(datadir+"/OrderParameter9_t120000.mat","bin");
    //saver.loadParameter<OrderParameter<10>>(datadir+"/OrderParameter10_t120000.mat","bin");
    //saver.loadParameter<OrderParameter<11>>(datadir+"/OrderParameter11_t60000.mat","bin");
    saver.loadParameter<Pressure<>>(datadir+"/Pressure_t1000000.mat","bin");
    saver.loadParameter<Velocity<>,Lattice::NDIM>(datadir+"/Velocity_t1000000.mat","bin");*/
    
    Geometry<Lattice>::initialiseBoundaries(initBoundary,{0});
    
    getAlgorithm getter;
    auto lbm = getter.get();

    saver.saveHeader(timesteps, saveInterval);
    for (int timestep=0; timestep<=timesteps; timestep++) {
        TIME=timestep;
        // Save the desired parameters, producing a binary file for each.
        if (timestep%saveInterval==0||(timestep<100000&&timestep%5000==0)) {
            if(mpi.rank==0)std::cout<<"Saving at timestep "<<timestep<<"."<<std::endl;
            
            saver.saveBoundaries(timestep);
            /*saver.saveParameter<ChemicalPotential<0>>(timestep,true);
            saver.saveParameter<ChemicalPotential<1>>(timestep,true);
            saver.saveParameter<ChemicalPotential<2>>(timestep,true);
            saver.saveParameter<ChemicalPotential<3>>(timestep,true);
            saver.saveParameter<ChemicalPotential<4>>(timestep,true);
            saver.saveParameter<ChemicalPotential<5>>(timestep,true);*/
            //saverecursive<ChemicalPotential, numcomponents>(saver, timestep);
            //saver.saveParameter<LaplacianOrderParameter<>>(timestep);
            //saver.saveParameter<LaplacianChemicalPotential<>>(timestep);
            saver.saveParameter<Pressure<>>(timestep);
            //saver.saveParameter<Density<>>(timestep);
            /*saver.saveParameter<OrderParameter<0>>(timestep,true);
            saver.saveParameter<OrderParameter<1>>(timestep,true);
            saver.saveParameter<OrderParameter<2>>(timestep,true);
            saver.saveParameter<OrderParameter<3>>(timestep,true);
            saver.saveParameter<OrderParameter<4>>(timestep,true);*/
            saverecursive<OrderParameter, numcomponents-1>(saver, timestep);
            saver.saveParameter<Velocity<>,Lattice::NDIM>(timestep);
            
        }

        // Evolve by one timestep
        lbm.evolve();
        /*if (timestep>=10000){
            #pragma omp for schedule(guided)
            for (int k = Lattice::HaloSize; k <Lattice::N - Lattice::HaloSize; k++) { //Loop over k
                Velocity<>::get<Lattice,Lattice::NDIM>(k,0)=0;
                Velocity<>::get<Lattice,Lattice::NDIM>(k,1)=0;
            }
        }*/
        
    }
    
}
