#include "main.hh"

#include <chrono>
#include <cstdlib>

int main(int argc, char **argv) {
    int seed;
    std::string seedstr;
#ifdef MPIPARALLEL
    mpi.init();

    MPI_Type_create_resized(MPI_INT, 0L, sizeof(Boundary<Lattice::NDIM>), &mMPIBoundary);
    MPI_Type_commit(&mMPIBoundary);

    int blocklengths[3] = {1, 1, Lattice::NDIM};
    MPI_Datatype types[3] = {MPI_INT, MPI_CXX_BOOL, MPI_INT8_T};
    std::array<int8_t, Lattice::NDIM> a;
    MPI_Aint offsets[3];
    offsets[0] = offsetof(Boundary<Lattice::NDIM>, Id);
    offsets[1] = offsetof(Boundary<Lattice::NDIM>, IsCorner);
    offsets[2] = offsetof(Boundary<Lattice::NDIM>, NormalDirection) + (size_t)((((char *)(&a) - (char *)(&a[0]))));

    MPI_Type_create_struct(3, blocklengths, offsets, types, &mMPIBoundary);
    MPI_Type_commit(&mMPIBoundary);
#endif
    if (argc > 1) {
        seed = std::atoi(argv[1]);
        seedstr = std::to_string(seed);
        initParams("input/input" + seedstr + ".txt");
    } else {
        seedstr = "";
        initParams("input.txt");
    }

    if (mpi.rank == 0 && argc == 0) {
        int ret;
        std::string tmp = "rm input/input" + seedstr + ".txt";
        const char *array = tmp.c_str();
        ret = std::system(array);
    }

    auto binary = initBinary<>();
    auto pressure = initPressure<>();
#ifndef KINETIC
    auto humidity = initHumidity<>();
#endif

    Geometry<Lattice>::initialiseBoundaries(initBoundary, {0, 22, 5, 6, 8, 9,123});
    OrderParameter<>::set<Lattice>(initFluid);
    OrderParameterOld<>::set<Lattice>(initFluid);
#ifndef KINETIC
    Humidity<>::set<Lattice>(initHumidity);
#endif

    /*Gradient<Pressure<>>::set<Lattice,Lattice::NDIM,0>(0.0);
    Gradient<Pressure<>>::set<Lattice,Lattice::NDIM,1>(0.0);
    GradientMixed<Pressure<>>::set<Lattice,Lattice::NDIM,0>(0.0);
    GradientMixed<Pressure<>>::set<Lattice,Lattice::NDIM,1>(0.0);
    Gradient<Pressure<>>::set<Lattice,9,0>(0.0);
    Gradient<Pressure<>>::set<Lattice,9,1>(0.0);
    Gradient<Pressure<>>::set<Lattice,9,2>(0.0);
    Gradient<Pressure<>>::set<Lattice,9,3>(0.0);
    Gradient<Pressure<>>::set<Lattice,9,4>(0.0);
    Gradient<Pressure<>>::set<Lattice,9,5>(0.0);
    Gradient<Pressure<>>::set<Lattice,9,6>(0.0);
    Gradient<Pressure<>>::set<Lattice,9,7>(0.0);
    Gradient<Pressure<>>::set<Lattice,9,8>(0.0);
    GradientMixed<Pressure<>>::set<Lattice,9,0>(0.0);
    GradientMixed<Pressure<>>::set<Lattice,9,1>(0.0);
    GradientMixed<Pressure<>>::set<Lattice,9,2>(0.0);
    GradientMixed<Pressure<>>::set<Lattice,9,3>(0.0);
    GradientMixed<Pressure<>>::set<Lattice,9,4>(0.0);
    GradientMixed<Pressure<>>::set<Lattice,9,5>(0.0);
    GradientMixed<Pressure<>>::set<Lattice,9,6>(0.0);
    GradientMixed<Pressure<>>::set<Lattice,9,7>(0.0);
    GradientMixed<Pressure<>>::set<Lattice,9,8>(0.0);
    ChemicalPotential<>::set<Lattice>(0.0);
    Pressure<>::set<Lattice>(1);
    Data_Base<Lattice, typename DefaultTraitPressureTernaryLee<Lattice>::Stencil>::getInstance().communicate(
        Pressure<>::getInstance<Lattice>());*/


#ifdef KINETIC
    Algorithm lbm(binary, pressure);
#else
    Algorithm lbm(humidity,binary, pressure);
    //Algorithm lbm(binary, pressure);
#endif

    using std::chrono::duration;
    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::milliseconds;

    SaveHandler<Lattice> saver(datadir);
    saver.saveHeader(timesteps, saveInterval);
    auto t1 = std::chrono::high_resolution_clock::now();
    for (int timestep = 0; timestep <= timesteps; timestep++) {
        TIME = timestep;
        // Save the desired parameters, producing a binary file for each.
        if (timestep % saveInterval == 0) {
            if (mpi.rank == 0) std::cout << "Saving at timestep " << timestep << "." << std::endl;

            saver.saveBoundaries(timestep);
#ifndef KINETIC
            saver.saveParameter<GradientHumidity<>, Lattice::NDIM>(timestep);
                        saver.saveParameter<Humidity<>>(timestep);
#endif
            saver.saveParameter<LaplacianChemicalPotential<>>(timestep);
            saver.saveParameter<ChemicalPotential<>>(timestep);
            saver.saveParameter<Density<>>(timestep);
            saver.saveParameter<Pressure<>>(timestep);
            saver.saveParameter<OrderParameter<>>(timestep);
            saver.saveParameter<MassSink<>>(timestep);
            //saver.saveParameter<MassSink2<>>(timestep);
            saver.saveParameter<Velocity<>, Lattice::NDIM>(timestep);
            
        }
        AfterEquilibration(timestep, binary);

        // Evolve by one timestep
#ifdef MPIPARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        lbm.evolve();
    }
    auto t2 = std::chrono::high_resolution_clock::now();

    /* Getting number of milliseconds as an integer. */
    auto ms_int = duration_cast<milliseconds>(t2 - t1);

    /* Getting number of milliseconds as a double. */
    duration<double, std::milli> ms_double = t2 - t1;

    std::cout << ms_int.count() << "ms\n";
    std::cout << ms_double.count() << "ms\n";
}
