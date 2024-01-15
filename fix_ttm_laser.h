/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(ttm/laser,FixTTMLaser);
// clang-format on
#else

#ifndef LMP_FIX_TTM_LASER_H
#define LMP_FIX_TTM_LASER_H

#include "fix.h"
#include <exception>

namespace LAMMPS_NS {

class FixTTMLaser : public Fix {
 public:
  FixTTMLaser(class LAMMPS *, int, char **);
  virtual ~FixTTMLaser();
  virtual void post_constructor();
  int setmask();
  virtual void init();
  void setup(int);
  void post_force_setup(int);
  virtual void post_force(int);
  void post_force_respa_setup(int, int, int);
  void post_force_respa(int, int, int);
  virtual void end_of_step();
  void reset_dt();
  void grow_arrays(int);
  virtual void write_restart(FILE *);
  virtual void restart(char *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();
  virtual double compute_vector(int);
  virtual double memory_usage();
  virtual double electronic_specific_heat(double);
  virtual double gamma_p(double);
  std::vector<double> electronic_specific_heat_array;
  std::vector<double> gamma_p_array;

 protected:
  int nlevels_respa;
  int seed;
  int nxgrid, nygrid, nzgrid;    // size of global grid
  int ngridtotal;                // total size of global grid
  int total_esh_array;
  int deallocate_flag;
  int outflag, outevery, pulseduration;
  double shift, tinit, tinit2;
  double v_esh, n_esh;
  double e_energy, transfer_energy;
  char *infile, *outfile, *beamfile;
  char *Cfile, *Gfile;

  class RanMars *random;
  double electronic_density;
  double electronic_thermal_conductivity;
  double T_dependant_electronic_thermal_conductivity;
  double &electronic_thermal_conductivity_scalar = T_dependant_electronic_thermal_conductivity;
  double electronic_specific_heat_300k;
  double gamma_s, v_0, v_0_sq;

  double *gfactor1, *gfactor2, *ratio, **flangevin;
  double ***T_electron, ***T_electron_old, ***T_beam;
  double ***net_energy_transfer, ***net_energy_transfer_all;
  double ***T_atomic;
  int ***nsum, ***nsum_all;
  double ***sum_vsq, ***sum_vsq_all;
  double ***sum_mass_vsq, ***sum_mass_vsq_all;

  virtual void allocate_grid();
  virtual void deallocate_grid();
  virtual void read_electron_temperatures(const std::string &);
  virtual void read_beam_temperatures(const std::string &);
  virtual void read_electronic_specific_heat(const std::string &);
  virtual void read_gamma_p(const std::string &);
  virtual void write_electron_temperatures(const std::string &);

  class parser_error : public std::exception {
    std::string message;

   public:
    parser_error(const std::string &mesg) { message = mesg; }
    const char *what() const noexcept { return message.c_str(); }
  };
};

}    // namespace LAMMPS_NS

#endif
#endif