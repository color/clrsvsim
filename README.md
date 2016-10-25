# clrsvsim
Color Genomics Structural Variant Simulator

# Background
Because structural variants are relatively rare, a simulation is sometimes used to tune and validate tools aimed at detecting and analyzing them.

Structural variant simulators such as `SVGen` or `SVSim` work by generating a genome that has structural variants (in `fasta` format), then applying simulators like `wgsim` to synthesize reads from the genome, with a noise model that emulates a particular sequencing technology. This is useful when trying to simulate data in the absence of real human samples. However, it has some drawbacks:
- The simulated reads do not model the assay capture process - depth will typically be stable over the entire sample, rather than concentrated around the probes.
- Any synthetic noise model will not accurately match the particular chemistry and setup used in a given lab.
