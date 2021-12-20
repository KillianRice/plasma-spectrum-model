#include <md.h>
// This function calculates the length of the simulation box based on the number of particles in the simulation.
// The output of this function is the box length in units of the Wigner-Seitz radius (a).
double Md::getBoxSizeForMIC(int N0)
{
  return pow(N0 * 4. * cts::pi / 3.,1./3.);
}
