# Chi-squared Hypothesis Comparison Tool

The GlueX experiment features nearly full angular coverage and can detect both neutral and charged particles. This enables reliable kinematic fitting for exclusive reactions. Since some particles decay into γ before identification, the kinematic fitter can evaluate alternative hypotheses for γ's possible progenitors. This tool helps identify events that are shared between different reconstruction hypotheses.

## Features

- Compare chi-squared values between N different reconstruction hypotheses
- Support for both best overall combo and best-per-beam-ID matching modes
- Support for two output modes; one conserving disk-space and one preserving all primary tree events
- Output includes matched events with corresponding chi-square values
- Built on ROOT's RDataFrame for efficient data processing

## Requirements

- ROOT environment (tested on ROOT 6.24/04)
- C++ compiler (tested on GCC 4.8.5)
- clang-format for automatic formatting (with the included .clang-format)
- Input files must contain a flat ROOT tree with the following required branches:
  - Event ID
  - Run Number
  - Beam ID
  - ChiSq
  - NDF (Number of Degrees of Freedom)

## Installation

If you have installed gluex_root_analysis from the central Makefile, `chisq_hypothesis_comparison` should have been automatically installed to your $PATH.

If you have made modifications to the code, you will need to rebuild using the Makefile found in the chisq_hypothesis_comparison directory (where this README is located):
```bash
make install
```

## Usage

### Pass the configuration file as the only argument

```bash
./compare_hypotheses config.ini
```


### Configuration file (see the example config.ini for structure and syntax)

- `file`: The path of each tree's .root file
- `tree`: Name of each tree

### Optional Flags

- `outfile`: Custom output filename (default: `<tree2>_hypothesesMatched.root`)
- `best_per_beam`: Match by best combo per beam ID (default: match by best overall combo)
- `preserve_combos`: Preserve all primary tree entries, rather than the default behavior of removing non-unique combos by χ²

## Output Format

By default, the program generates a ROOT file containing:
- A copy of the primary tree with improbable, non-unique combos removed (see the preserve output mode above if non-unique combos need to be preserved)
- New branch with matched secondary combos' χ²/NDF values (the branch is named in the format: [secondary_tree_name]_chisq_ndf)

## Matching Criteria

Two matching modes are available:
1. Best Overall Combo: Selects the combo with lowest χ² per event
2. Best Per Beam ID: Selects the combo with lowest χ² for each unique event+beam ID pair; to be used alongside accidental subtraction

## Performance

The performance of the tool is limited by the ROOT library's I/O efficiency. In a test run comparing two trees with around 300,000 events each, the matching process takes 2 seconds for best overall mode and 3 seconds for best per beam matching mode. Writing the output to file using the ROOT library takes 30 seconds.

## Future Development

- [X] Benchmarking support
- [X] Output only unique combos
- [X] Rewrite implementation of matching logic
- [X] Use configuration file instead of CLI flags
- [X] Make the logging of matching info optional
- [X] Support simultaneous comparison of 3+ hypotheses
- [X] Post warning if user passes an empty tree (likely due to forgetting to fill flat tree)

## Authors

Jake Serwe - Florida State University

## Credit

This tool is intended as a successor to Alex Barnes' [mergeTrees.C](https://github.com/JeffersonLab/halld_recon/tree/master/src/programs/Utilities/mergeTrees)
