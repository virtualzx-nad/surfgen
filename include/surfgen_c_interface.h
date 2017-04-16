// C interface for Surfgen evaluation libraries.  These simply translate function names, adjusts to use
// more common `int` types instead of the 64 bit int, and reorder arrays into C order

#ifndef __SURFGEN_C_H__
#define __SURFGEN_C_H__
void InitPotential();
void FinalizeSurfgen();
void OpenTrajFile(int traj_id);
void EnableParsing();
void DisableParsing();
void PauseParsing();
void ResumeParsing();
void GetInfo(int *num_atoms, int *num_states);
int GetNumAtoms();
int GetNumStates();
void EvaluateSurfgen(double*geom, double *energy, double*cgrads, double*hd_matrix, double*diab_grads);
void GetEnergy(double* geom, double* energy);

#endif
