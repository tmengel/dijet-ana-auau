#ifndef MYANA_RHOETACALIBLOOKUP_H
#define MYANA_RHOETACALIBLOOKUP_H

// ============================================================================
// RhoEtaCalibLookup.h -- pick the rho eta-shape calibration file for a run.
//
//   #include <myana/RhoEtaCalibLookup.h>
//   const std::string calib = RhoEtaCalibLookup::GetCalibPath( run_number );
//   subrho -> set_etaCalib_directPath( calib );
//   dtb    -> SetEtaCalib_DirectPath( calib );
//
// Placeholder, file-based lookup (no CDB). The dataset is inferred from the run
// number; the run's own calibration is returned if it exists, otherwise the
// dataset default:
//
//   dataset    run numbers        per-run file                    default
//   sim        run < 1000         sim/rho_calib_<run>.root        sim/rho_calib_default.root (= run 30)
//   run2auau   50000 - 59999      run2auau/rho_calib_<run>.root   run2auau/rho_calib_default.root (= run 54590)
//
// Anything else has no calibration: an empty string is returned with a
// warning, and SubtractTowersRhov1 falls back to a flat rho profile (w = 1).
//
// variant "scaled" selects the UE-scaled HIJING calibration
// (sim/rho_calib_<run>_scaled.root, falling back to sim/rho_calib_30_scaled.root).
//
// The base directory can be overridden with the RHO_ETA_CALIB_DIR environment
// variable. Calibrations and how they are made: see README.md in that directory.
// ============================================================================

#include <TString.h>
#include <TSystem.h>

#include <cstdlib>
#include <iostream>
#include <string>

namespace RhoEtaCalibLookup
{
  inline const std::string kDefaultBaseDir = "/sphenix/user/tmengel/dijet-ana-auau/calibrations/rho_eta";

  enum class Dataset { Unknown, Sim, Run2AuAu };

  // run whose calibration each dataset's rho_calib_default.root was copied from
  inline const int kSimDefaultRun = 30;
  inline const int kRun2AuAuDefaultRun = 54590;

  inline std::string BaseDir()
  {
    const char *env = std::getenv("RHO_ETA_CALIB_DIR");
    return (env && *env) ? std::string(env) : kDefaultBaseDir;
  }

  inline Dataset DatasetForRun(const int run_number)
  {
    if (run_number >= 0 && run_number < 1000) { return Dataset::Sim; }
    if (run_number >= 50000 && run_number < 60000) { return Dataset::Run2AuAu; }
    return Dataset::Unknown;
  }

  inline std::string DatasetName(const Dataset d)
  {
    switch (d)
    {
    case Dataset::Sim: return "sim";
    case Dataset::Run2AuAu: return "run2auau";
    default: return "unknown";
    }
  }

  inline bool Exists(const std::string &path) { return !gSystem->AccessPathName(path.c_str()); }

  /// Calibration file for run_number: its own if present, else the dataset default.
  /// Returns "" (no calibration) when the run belongs to no known dataset.
  inline std::string GetCalibPath(const int run_number, const std::string &variant = "", const bool verbose = true)
  {
    const Dataset d = DatasetForRun(run_number);
    if (d == Dataset::Unknown)
    {
      std::cout << "RhoEtaCalibLookup::GetCalibPath - WARNING: run " << run_number
                << " is not in a dataset with a rho eta-shape calibration; returning none "
                << "(flat rho, w = 1)" << std::endl;
      return "";
    }
    if (!variant.empty() && !(d == Dataset::Sim && variant == "scaled"))
    {
      std::cout << "RhoEtaCalibLookup::GetCalibPath - WARNING: unknown variant \"" << variant << "\" for "
                << DatasetName(d) << "; returning none (flat rho, w = 1)" << std::endl;
      return "";
    }

    const std::string dir = BaseDir() + "/" + DatasetName(d);
    const std::string suffix = variant.empty() ? "" : "_" + variant;
    const std::string own = std::string(Form("%s/rho_calib_%d%s.root", dir.c_str(), run_number, suffix.c_str()));
    if (Exists(own))
    {
      if (verbose) { std::cout << "RhoEtaCalibLookup::GetCalibPath - run " << run_number << ": " << own << std::endl; }
      return own;
    }

    // default: plain copy for the nominal calibration; variants fall back to the
    // default run's variant file
    const std::string def = variant.empty()
                                ? dir + "/rho_calib_default.root"
                                : std::string(Form("%s/rho_calib_%d%s.root", dir.c_str(), kSimDefaultRun, suffix.c_str()));
    if (!Exists(def))
    {
      std::cout << "RhoEtaCalibLookup::GetCalibPath - ERROR: no calibration for run " << run_number
                << " and default " << def << " is missing; returning none (flat rho, w = 1)" << std::endl;
      return "";
    }
    if (verbose)
    {
      std::cout << "RhoEtaCalibLookup::GetCalibPath - no calibration for run " << run_number << " ("
                << DatasetName(d) << "), using default: " << def << std::endl;
    }
    return def;
  }
}  // namespace RhoEtaCalibLookup

#endif
