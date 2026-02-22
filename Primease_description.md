# Primease — Detailed Technical Summary & Module Relationships

> Purpose: automated primer design pipeline that generates primer candidates with Primer3, screens them by alignment with Bowtie (target + optional host genomes), filters off-target primers, and produces visualizations and downloadable results.

This document explains **how the tool works**, how each module interacts, the runtime data flow, important files...

---

## High-level pipeline (single run)

1. **User input**: upload/paste a DNA template (FASTA / GenBank / SnapGene) and choose a workflow:
   - Single sequence of interest
   - Two overlapping fragments
   - Annotation-based selection from GenBank
2. **Preprocessing**:
   - `fix_genbank_locus()` (if needed) repairs LOCUS line format.
   - `Target` loads the first record and writes `target.fasta`.
3. **Annotation / Region selection**:
   - `SequenceAnnotator` finds a pasted sequence in the template and annotates it.
   - `OverlappingFragmentAnnotator` finds overlap between two fragments and annotates that region.
   - `GenbankfileHandling` extracts regions directly from annotated genbank features and resolves overlapping fragments.
   - These produce `annotated_sequence_expected_structure.gb`.
4. **Selection → ROI objects**:
   - `Selection` subclasses produce a list of `ROI` objects (each ROI encapsulates name, start, end, and expected_structure_sequence).  
5. **Primer design**:
   - `Primer3Interface` writes a `settings` file (default `settings.bak` or user-provided `settings_spec.bak`), appends the ROI block, sets `PRIMER_PRODUCT_SIZE_RANGE`, runs Primer3, parses the text output into a DataFrame and writes `potential_primers.fasta`.
6. **Build Bowtie index & align**:
   - `BowtieInterface.create_index()` moves `target.fasta` into a task-specific `bowtie_index/<uuid>/target/` and runs `bowtie-build` to create an index for the target.
   - Host genomes (selected by UI) are copied into `bowtie_index/<uuid>/host/` and their indices are discovered.
   - `BowtieInterface.run_bowtie()` runs `bowtie -x <index> -a -f potential_primers.fasta -v 3` for the target and each host.
   - Alignments are written to `bt_target.csv` and a concatenated (in case of severa genomes) `bt_host_*.csv` → `bt_host.csv`.
7. **Parse alignments**:
   - `BowtieResult` loads Bowtie outputs, ensures columns (`name`, `strand`, `reference`, `start`, `sequence`, `quality`, `instances`), and derives:
     - `id` (primer identifier) from `name`
     - `orientation` (fwd/rev) from `name`
8. **Off-target filtering**:
   - `OfftargetChecker` collects:
     - **Sponges**: primer `id` that binds ≥ `config.sponge_value` times (grouped by `(id, orientation)`).
     - **Off-target amplicons**: for an `(id, reference)` group, if both strands are present and at least one forward–reverse pair is within `config.offtarget_size_cutoff`, mark the `id` as offtarget.
     - Remove any primer pairs that align anywhere in the host (by exact name).
   - Keep top `config.top` primers per position.
   - Output: `qTagGer_Output.csv` and `final_primers` DataFrame.
9. **Visualization**:
   - `Visualization` maps primer sequences back onto the full sequence (search forward primers directly, reverse primers via reverse complement), draws arrows and connecting lines, highlights SOI region, and saves `primers_plot.png`.
   - `VisualizeGenbank` draws a GenBank-level feature map and saves `Genbank_vizualisation.png`.
10. **Cleanup**:
    - `cleanup()` moves output files to a timestamped directory, removes temp files (e.g., `target.fasta`, `annotated_sequence_expected_structure.gb`), and deletes the `bowtie_index/<uuid>` host/task dir.

---

## Core modules and responsibilities

### `main` (Streamlit app)
- Presents UI (upload/paste, workflow selection, option to upload custom Primer3 settings, host genome choice).
- Accepts sequence(s), triggers the chosen pipeline:
  - `main_sequence_of_interest`
  - `main_fragments` (overlapping fragments)
  - `main_annotated_genbank` (aleary correctly annotated genome)
- Handles session state, progress bars and download buttons.
- Orchestrates `copy_reference_genome()` to assemble host references in a unique task dir.

### `Target`
- Loads first SeqRecord from file (GenBank/SnapGene/FASTA).
- Writes `target.fasta` for Bowtie indexing.
- Exposes `record` and uppercased `seq`.

### `Selection` & subclasses
- `Selection` (abstract) defines `selection()` → list of `ROI`.
- `SequenceOfInterest` scans `target.record.features`, finds `fragment` features or `misc_feature` with `locus_tag='fragment'`, builds `ROI`s.
- `GenbankfileHandling` loads an expected-structure genbank, extracts features, detects overlaps between feature sequences (largest suffix/prefix overlap > `min_overlap`), merges overlapping features into SOIs, and builds `ROI`s and `list_sequence_of_interest`.
- `SequenceOfInterest` and `GenbankfileHandling` are how the pipeline determines the genomic coordinates and sequences that Primer3 should target.

### `ROI` (referenced)
- Expected to contain: `name`, `target_start`, `target_end`, `expected_structure_sequence`.

### `Primer3Interface`
- Writes a Primer3 `settings` file (from default or user-provided), appends the ROI block and ensures `PRIMER_PRODUCT_SIZE_RANGE`.
- Calls Primer3 binary (`<primer3_path> settings > primer3_result`).
- Parses Primer3 text output with regex to extract primer pairs and metrics into a DataFrame (`primer_sites`), names primers as `<SEQUENCE_ID>_<index>`.
- Writes `potential_primers.fasta` with entries `>name_fwd` and `>name_rev` for Bowtie.

### `BowtieInterface`
- Creates per-run task directory under `bowtie_index/<uuid>/`.
- Moves `target.fasta` and builds Bowtie index for the target (`bowtie-build`).
- Detects host FASTA files in `task_dir/host` and runs Bowtie against each host index if present.
- Runs Bowtie with `-a -f` to report all alignments for primers in FASTA, allowing up to 3 mismatches (`-v 3`).
- Aggregates host results into a single host output file.
- Produces `BowtieResult` objects for both target and host outputs.

### `BowtieResult`
- Robustly parses Bowtie tabular output (handles empty/missing files).
- Ensures consistent DataFrame columns and types.
- Extracts `id` (primer id) and `orientation` from `name` tokens.
- Exposes `.result` DataFrame consumed by `OfftargetChecker`.

### `OfftargetChecker`
- Accepts primer candidate DataFrame, `BowtieResult` for target and host, and `Config`.
- Computes:
  - **Sponges**: IDs that occur ≥ `sponge_value` times for a given `(id, orientation)`.
  - **Off-target amplicons**: in a given `(id, reference)`, if forward and reverse alignments exist with proper strand orientation and distance ≤ `offtarget_size_cutoff`, flag the id.
  - Removes primer pairs found aligning to the host (exact names).
- Outputs `final_primers` and writes `qTagGer_Output.csv`.

### `Visualization` & `VisualizeGenbank`
- `Visualization`: maps primers to positions using string matching (forward) and reverse complement (reverse), draws arrow patches for primers, connects forward/reverse pairs, highlights SOI region, centers plot on SOI (±500 bp), and saves `primers_plot.png`.
- `VisualizeGenbank`: uses `dna_features_viewer.BiopythonTranslator` to plot all features and label SOIs; saves `Genbank_vizualisation.png`.

### Utilities
- `handle_file_upload()` saves uploaded files locally.
- `check_file_structure()` validates FASTA/GenBank parseability with `Bio.SeqIO`.
- `fix_genbank_locus()` repairs LOCUS line formatting when necessary.
- `initial_cleanup()` and `cleanup()` manage working files and archival.

---

## File-level dataflow (key files created / consumed)

- **Inputs**
  - `uploaded` FASTA / GenBank / pasted sequence
  - `settings_spec.bak` (optional custom Primer3 settings)

- **Intermediate / generated**
  - `pasted_sequence.fasta` (when user pastes)
  - `annotated_sequence_expected_structure.gb` (after annotation step)
  - `target.fasta` (written by `Target`)
  - `settings` (Primer3 settings file, generated from `settings.bak` or `settings_spec.bak`)
  - `primer3_result` (text output from Primer3)
  - `potential_primers.fasta` (FASTA of forward and reverse primer sequences)
  - `bt_target.csv` (Bowtie alignments vs. target)
  - `bt_host.csv` (concatenated Bowtie alignments vs. host genomes)

- **Outputs**
  - `qTagGer_Output.csv` (final primer table)
  - `primers_plot.png` (visualization centered on SOI)
  - `Genbank_vizualisation.png` (whole-genome/genbank view)
  - Archive folder: `<jobname>_YYYYMMDD_HHMM/` containing outputs

---

## Key configuration & parameters

- `Config.primer3_path` — path to Primer3 binary used to run designs.
- `Config.top` — number of primer candidates to keep per ROI position.
- `Config.sponge_value` — threshold number of alignments that mark a primer as a sponge (repetitive).
- `Config.offtarget_size_cutoff` — maximum allowed distance between forward/reverse alignments to be considered an off-target amplicon (bp).
- Host genome selection is managed by the UI and `copy_reference_genome()` which copies reference FASTA files into the per-task host dir.

---

## Error handling & notable failure conditions

- **Missing files**: `Target.load_record()` raises if input extension is unknown. Bowtie index creation raises `FileNotFoundError` if `target.fasta` missing.
- **Empty or malformed outputs**: `BowtieResult` robustly returns an empty DataFrame for missing/empty Bowtie outputs.
- **No primers from Primer3**: `Primer3Interface.parse_results()` raises `ValueError` explaining possible causes.
- **All primers removed by host alignment**: `OfftargetChecker.remove_host_aligned_pairs()` raises `RuntimeError` if no primer remains.
- **Sequence not found**: `SequenceAnnotator.run()` raises `ValueError` if pasted SOI not found in template.
- **Overlap not found / too short**: `OverlappingFragmentAnnotator` raises `ValueError` when overlaps are missing or < 10 bp.

