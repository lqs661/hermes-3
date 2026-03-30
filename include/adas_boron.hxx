#pragma once
#ifndef ADAS_BORON_H
#define ADAS_BORON_H

#include "adas_reaction.hxx"

#include <array>
#include <initializer_list>

/// Ionisation energies in eV
/// from https://www.webelements.com/boron/atoms.html
/// Conversion 1 kJ mol^-1 = 1.0364e-2 eV
/// These are added (removed) from the electron energy during recombination (ionisation)
constexpr std::array<BoutReal, 5> boron_ionisation_energy{
    8.30, 25.16, 37.93, 259.37, 340.22};

/// The name of the species.
///
/// b, b+, b+2, b+3, ...
///
/// Special cases for level=0 and 1
///
/// @tparam level  The ionisation level: 0 is neutral, 5 is fully stripped.
template <int level>
constexpr std::initializer_list<char> boron_species_name{'b', '+', '0' + level};

template <>
constexpr std::initializer_list<char> boron_species_name<1>{'b', '+'};

template <>
constexpr std::initializer_list<char> boron_species_name<0>{'b'};

/// ADAS effective ionisation (ADF11)
///
/// @tparam level  The ionisation level of the ion on the left of the reaction
template <int level>
struct ADASBoronIonisation : public OpenADAS {
  ADASBoronIonisation(std::string, Options& alloptions, Solver*)
      : OpenADAS(alloptions, "scd89_b.json", "plt89_b.json",
                 boron_species_name<level>, boron_species_name<level + 1>, level,
                 -boron_ionisation_energy[level]) {}

private:
  void transform_impl(GuardedOptions& state) override {
    calculate_rates(
        state["species"]["e"],                          // Electrons
        state["species"][boron_species_name<level>],    // From this ionisation state
        state["species"][boron_species_name<level + 1>] // To this state
    );
  }
};

/// ADAS effective recombination coefficients (ADF11)
///
/// @tparam level  The ionisation level of the ion on the right of the reaction
template <int level>
struct ADASBoronRecombination : public OpenADAS {
  ADASBoronRecombination(std::string, Options& alloptions, Solver*)
      : OpenADAS(alloptions, "acd89_b.json", "prb89_b.json",
                 boron_species_name<level + 1>, boron_species_name<level>, level,
                 boron_ionisation_energy[level]) {}

private:
  void transform_impl(GuardedOptions& state) override {
    calculate_rates(
        state["species"]["e"],                           // Electrons
        state["species"][boron_species_name<level + 1>], // From this ionisation state
        state["species"][boron_species_name<level>]      // To this state
    );
  }
};

/// @tparam level     The ionisation level of the ion on the right of the reaction
/// @tparam Hisotope  The hydrogen isotope ('h', 'd' or 't')
template <int level, char Hisotope>
struct ADASBoronCX : public OpenADASChargeExchange {
  ADASBoronCX(std::string, Options& alloptions, Solver*)
      : OpenADASChargeExchange(alloptions, "ccd89_b.json",
                               boron_species_name<level + 1>, {Hisotope},
                               boron_species_name<level>, {Hisotope, '+'}, level) {}

private:
  void transform_impl(GuardedOptions& state) override {
    GuardedOptions species = state["species"];
    calculate_rates(
        species["e"],                             // Electrons
        species[boron_species_name<level + 1>],  // From this ionisation state
        species[{Hisotope}],                      // and this neutral hydrogen isotope
        species[boron_species_name<level>],      // To this state
        species[{Hisotope, '+'}]                  // and this hydrogen ion
    );
  }
};

namespace {
// Ionisation by electron-impact
RegisterComponent<ADASBoronIonisation<0>> register_ionisation_b0("b + e -> b+ + 2e");
RegisterComponent<ADASBoronIonisation<1>> register_ionisation_b1("b+ + e -> b+2 + 2e");
RegisterComponent<ADASBoronIonisation<2>> register_ionisation_b2("b+2 + e -> b+3 + 2e");
RegisterComponent<ADASBoronIonisation<3>> register_ionisation_b3("b+3 + e -> b+4 + 2e");
RegisterComponent<ADASBoronIonisation<4>> register_ionisation_b4("b+4 + e -> b+5 + 2e");

// Recombination
RegisterComponent<ADASBoronRecombination<0>> register_recombination_b0("b+ + e -> b");
RegisterComponent<ADASBoronRecombination<1>> register_recombination_b1("b+2 + e -> b+");
RegisterComponent<ADASBoronRecombination<2>> register_recombination_b2("b+3 + e -> b+2");
RegisterComponent<ADASBoronRecombination<3>> register_recombination_b3("b+4 + e -> b+3");
RegisterComponent<ADASBoronRecombination<4>> register_recombination_b4("b+5 + e -> b+4");

// Charge exchange
RegisterComponent<ADASBoronCX<0, 'h'>> register_cx_b0h("b+ + h -> b + h+");
RegisterComponent<ADASBoronCX<1, 'h'>> register_cx_b1h("b+2 + h -> b+ + h+");
RegisterComponent<ADASBoronCX<2, 'h'>> register_cx_b2h("b+3 + h -> b+2 + h+");
RegisterComponent<ADASBoronCX<3, 'h'>> register_cx_b3h("b+4 + h -> b+3 + h+");
RegisterComponent<ADASBoronCX<4, 'h'>> register_cx_b4h("b+5 + h -> b+4 + h+");

RegisterComponent<ADASBoronCX<0, 'd'>> register_cx_b0d("b+ + d -> b + d+");
RegisterComponent<ADASBoronCX<1, 'd'>> register_cx_b1d("b+2 + d -> b+ + d+");
RegisterComponent<ADASBoronCX<2, 'd'>> register_cx_b2d("b+3 + d -> b+2 + d+");
RegisterComponent<ADASBoronCX<3, 'd'>> register_cx_b3d("b+4 + d -> b+3 + d+");
RegisterComponent<ADASBoronCX<4, 'd'>> register_cx_b4d("b+5 + d -> b+4 + d+");

RegisterComponent<ADASBoronCX<0, 't'>> register_cx_b0t("b+ + t -> b + t+");
RegisterComponent<ADASBoronCX<1, 't'>> register_cx_b1t("b+2 + t -> b+ + t+");
RegisterComponent<ADASBoronCX<2, 't'>> register_cx_b2t("b+3 + t -> b+2 + t+");
RegisterComponent<ADASBoronCX<3, 't'>> register_cx_b3t("b+4 + t -> b+3 + t+");
RegisterComponent<ADASBoronCX<4, 't'>> register_cx_b4t("b+5 + t -> b+4 + t+");
} // namespace

#endif // ADAS_BORON_H
