# Changelog

## v1.3-r1 (code structure, resume and threads) (2026-11-25)

### Added

- `--version` for `EVG`, `graphvcf` and `fastAQ`.
- `-t/--threads` for every `graphvcf` and `fastAQ` subcommand (default 4); results do not depend on the thread count.
- `fastAQ sample -o` writes the reads to a file.

### Changed

- `--restart` is replaced by `--resume`: a command is redone if its options, inputs or outputs changed.
- Clearer code layout and one log format for all three programs.
- htslib is bundled; `.gz` outputs are BGZF and the merged VCFs are indexed with tabix.
- Sample names may contain only letters, digits, `_`, `.`, `+` and `-`, and must not start with `.` or `-`.

### Fixed

- `graphvcf filter` (and `convert --maf/--geno`) dropped every SNP at its defaults.
- Samples with the same read file names overwrote each other's downsampled reads.
- A failed job or a failed part of a pipe did not stop the run.
- Crashes on empty input, on GT not being the first FORMAT field, and on haploid `GT=0` in `split -f`.
- Damaged input and write errors (such as a full disk) were ignored.
- Version checks compared strings (`vg` 1.100 counted as older than 1.63).
- Build with GCC 13.
