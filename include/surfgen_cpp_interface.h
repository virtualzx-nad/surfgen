#ifndef __SURFGEN_CPP_H__
#define __SURFGEN_CPP_H__
#include <vector>

using namespace std;

class SurfgenPotential{
  protected:
    int _natoms;
    int _nstates;

  public:
    SurfgenPotential();
    ~SurfgenPotential();

    int getNumAtoms();
    int getNumStates();

    void startTrajectory(int traj_index);

    void enableParsing(); 
    void disableParsing(); 
    void pauseParsing();
    void resumeParsing();

    void evaluate(double*geom, double&energy, double*grads, double*Hd_matrix, double*dcgrads);
    void getEnergy(double* geom, double*energy);

};
#endif
