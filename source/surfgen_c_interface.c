#include "surfgen_c_interface.h"
// C interface for Surfgen evaluation libraries.  These simply translate function names, adjusts to use
// more common `int` types instead of the 64 bit int, and reorder arrays into C order
extern void finalizesurfgen_();
extern void opentrajfile_(long long *itraj);
extern void initpotential_();
extern void evaluatesurfgen_(double*cgeom,double*energy,double*cgrads,double*hmat,double*dcgrads);
extern void getenergy_(double* cgeom, double*energy);
extern void pauseparsing_();
extern void resumeparsing_();
extern void enableparsing_();
extern void disableparsing_();
extern void getinfo_(long long *natoms, long long *nstates);

void InitPotential(){
  initpotential_();
}

void FinalizeSurfgen(){
    finalizesurfgen_();
}

void OpenTrajFile(int traj_id){
  long long trajid64 = traj_id;
  opentrajfile_(&trajid64);
}

void EnableParsing(){
  enableparsing_();
}

void DisableParsing(){
  disableparsing_();
}

void PauseParsing(){
  pauseparsing_();
}

void ResumeParsing(){
  resumeparsing_();
}

void GetInfo(int *num_atoms, int *num_states){
  long long natoms64, nstates64;
  getinfo_(&natoms64, &nstates64);
  *num_atoms = natoms64;
  *num_states = nstates64;
} 

int GetNumAtoms(){
  int num_atoms, num_states;
  GetInfo(&num_atoms, &num_states);
  return num_atoms;
}

int GetNumStates(){
  int num_atoms, num_states;
  GetInfo(&num_atoms, &num_states);
  return num_states;
}

void EvaluateSurfgen(double*geom, double *energy, double*grads, double*hd_matrix, double*diab_grads){
  evaluatesurfgen_(geom, energy, grads, hd_matrix, diab_grads);
}

void GetEnergy(double* geom, double*energy){
  getenergy_(geom, energy);
}


