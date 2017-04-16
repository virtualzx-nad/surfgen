#include "surfgen_cpp_interface.h"

extern "C"{
  void finalizesurfgen_();
  void opentrajfile_(long long *itraj);
  void initpotential_();
  void evaluatesurfgen_(double*cgeom,double*energy,double*cgrads,double*hmat,double*dcgrads);
  void getenergy_(double* cgeom, double*energy);
  void pauseparsing_();
  void resumeparsing_();
  void enableparsing_();
  void disableparsing_();
  void getinfo_(long long *natoms, long long *nstates);
}

SurfgenPotential::SurfgenPotential(){
  initpotential_();
  long long natoms, nstates;
  getinfo_(&natoms, &nstates);
  _natoms = natoms;
  _nstates = nstates;
}

SurfgenPotential::~SurfgenPotential(){
  finalizesurfgen_();
}

int SurfgenPotential::getNumAtoms(){
  return _natoms;
}
int SurfgenPotential::getNumStates(){
  return _nstates;
}

void SurfgenPotential::startTrajectory(int traj_index){
  long long trajid64 = traj_index;
  opentrajfile_(&trajid64);
}

void SurfgenPotential::enableParsing(){
  enableparsing_();
}

void SurfgenPotential::disableParsing(){
  disableparsing_();
}

void SurfgenPotential::pauseParsing(){
  pauseparsing_();
}
void SurfgenPotential::resumeParsing(){
  resumeparsing_();
}

void SurfgenPotential::evaluate(double*geom, double&energy, double*grads, double*Hd_matrix, double*dcgrads){
  evaluatesurfgen_(geom, &energy, grads, Hd_matrix, dcgrads);
}

void SurfgenPotential::getEnergy(double* geom, double* energy){
  getenergy_(geom, energy);
}
