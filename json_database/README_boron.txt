Boron ADAS data files required by Hermes-3
=========================================

To enable boron reaction components (see `include/adas_boron.hxx`),
place the following JSON files in this directory:

- scd89_b.json   (ionisation rates)
- plt89_b.json   (ionisation radiation losses)
- acd89_b.json   (recombination rates)
- prb89_b.json   (recombination radiation losses)
- ccd89_b.json   (charge exchange rates)

Expected format matches other OpenADAS JSON files here (Hermes reads
`log_temperature` = log10(Te/eV), `log_density` = log10(ne/m^-3),
and `log_coeff` with the usual -6 offset from the .dat file).

ADF11 .dat grids (89 and 96): first N1 numbers after the header are
log10(ne/cm^-3); the next N2 numbers are log10(Te/eV). For plt89_b,
that is roughly 10…15 on the density axis (10^10…10^15 cm^-3) and
0…4.7 on the temperature axis (1 eV … ~5e4 eV).

Regenerate JSON from .dat with `python adas_dat_to_json.py` in this folder.

If these files are missing, boron reactions will throw at runtime when
the corresponding component is created.

OpenADAS JSON paths are resolved like Amjuel: by default the directory
json_database next to the Hermes repository (compile-time __FILE__ path),
not the shell current working directory. Override with the root option
json_database_dir if you install binaries without the full source tree.
