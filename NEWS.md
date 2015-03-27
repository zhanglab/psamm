
v0.3
----

- Require reaction files to be explicitly listed in `model.yaml`.
- Add support for TSV reaction file format.
- Change format of YAML reactions (see [README](README.md) for details).
- Add tables of recognized compounds and reaction properties to
  [README](README.md).
- `masscheck` command: Automatically exclude biomass reaction from check.

v0.2
----

- Allow compounds to be specified using YAML, TSV or ModelSEED format. This
  changes the format of the `compounds` key in `model.yaml` (see
  [README](README.md) for more information).
- Allow specifying biomass reaction in `model.yaml` using the `biomass` key.
  The biomass reaction will be used by default for FBA, FVA, etc.
- Allow explicit definition of media. This can be defined using a table format
  or YAML format. See [README](README.md) for more information.
- `chargecheck`/`formulacheck` commands: Only check reactions where all
  compounds have charge/formula specified. The number of skipped reactions is
  reported separately.
- `chargecheck` command: Use charge information from model definition instead
  of requiring a separate charge table file.

v0.1
----

- Initial release.
