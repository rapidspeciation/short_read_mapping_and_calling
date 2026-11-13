# Snakefile
# NOTE — Start fresh screen first. [> screen -S <screen-name>]
# NOTE — Run iinit before running and enter PAM password. [> iinit]
# NOTE — Activate conda environment set up for snakemake. [> module load conda; > conda activate snakemake_7.8.0]
# NOTE — Navigate to the run working directory. [> cd run_directory]
# NOTE — Submit with command: > snakemake -s ../scripts/reads_to_groupBCF.v2.5.snakefile --profile ~/.config/snakemake/lsf-sanger-1 --restart-times 2 --keep-going --rerun-incomplete --resources iRODS_jobs=10 -j 999 -n
# NOTE — Get rule graph: > snakemake -s ../scripts/reads_to_groupBCF.v2.5.snakefile --profile ~/.config/snakemake/lsf-sanger-1 --rerun-incomplete -j 999 --rulegraph | dot -Tsvg > rulegraph.svg
# NOTE — Get run graph: > snakemake -s ../scripts/reads_to_groupBCF.v2.5.snakefile --profile ~/.config/snakemake/lsf-sanger-1 --rerun-incomplete --resources iRODS_jobs=10 -j 999 --dag | dot -Tsvg > dag.svg

# Set the base directory under which all read, CRAM and BCF output folders are written.
BASE_DIR = "/lustre/scratch125/tol/teams/meier/users/jw42/maratus_popgen_extension"

# Set paths to executables.
samtools = "/software/badger/module-builds/samtools/1.22.1/bin/samtools"
fastp = "/software/treeoflife/shpc/0.1.26/wrapper/quay.io/biocontainers/fastp/0.23.4--hadf994f_3/bin/fastp"
minibwa = "/software/team347/jw42/minibwa/minibwa"
bcftools = "/software/badger/module-builds/bcftools/1.24/bin/bcftools"

# Import python modules
from collections import defaultdict
import pandas as pd
import os
import yaml
import re


### LOAD REFERENCE INFORMATION
# Read the reference metadata table.
references_table = pd.read_csv('../scripts/reads_to_groupBCF.references_metadata.tsv', sep="\t", keep_default_na=False).set_index("reference_id")
reference_path_for = references_table["reference_path"].to_dict()
masked_rm2_bed_for = references_table["masked_rm2_bed"].to_dict()
busco_exons_helixer_bed_for = references_table["busco_exons_helixer_bed"].to_dict()
busco_exons_braker3_bed_for = references_table["busco_exons_braker3_bed"].to_dict()
busco_exons_metaeuk_bed_for = references_table["busco_exons_metaeuk_bed"].to_dict()
# List every reference_id in the reference metadata table.
REFERENCE_IDS = list(references_table.index)
# Read the reference intervals table, which splits each reference into named intervals for parallel variant calling.
intervals_table = pd.read_csv('../scripts/reads_to_groupBCF.references_intervals.tsv', sep="\t")
# Make a lookup from (reference_id, interval) to that interval's coordinates.
interval_coordinates = {(row.reference_id, str(row.interval)): row.coordinates for row in intervals_table.itertuples()}
# Make a lookup from reference_id to the list of its interval names.
intervals_by_reference = intervals_table.groupby("reference_id")["interval"].apply(lambda s: [str(v) for v in s]).to_dict()
# Set the ploidy file, which gives per-chromosome ploidy for each sex for bcftools call.
ploidy_file = "../scripts/reads_to_groupBCF.references_ploidy.tsv"


### LOAD SPECIMEN SEX INFORMATION
# Set the specimen sex file, which gives the sex of each specimen for bcftools call.
specimen_sex_file = "../scripts/reads_to_groupBCF.specimen_sex.tsv"


### LOAD BCF TARGETS AND PROCESS
# Set the YAML file listing the BCF targets to produce, read from the directory snakemake is run in.
BCF_TARGETS_FILE = "BCF_targets.yaml"
# Name the four fields that together identify one BCF target.
BCF_KEY_FIELDS = ("reference_id", "dataset_prefix", "bcf_subset", "bcf_callgroup")
# Define a function giving the callgroup file for a bcf_callgroup, which is named after it.
def callgroup_file_path(bcf_callgroup):
    return f"../scripts/reads_to_groupBCF.callgroup.{bcf_callgroup}.tsv"
# Name the columns every run metadata file must have.
RUN_METADATA_COLUMNS = ("sequence_id", "specimen_id", "library_id", "path")
# Define a function giving the run metadata file for a dataset_prefix, which is named after it.
def run_metadata_path(dataset_prefix):
    return f"../scripts/reads_to_groupBCF.run_metadata.{dataset_prefix}.tsv"
# Define a function that loads and checks the BCF targets file, collecting every problem before stopping, and returns the specimen lists, callgroup files, and run metadata tables.
def load_bcf_targets(path):
    # Read the YAML, treating an empty file as an empty list of blocks.
    with open(path) as fh:
        blocks = yaml.safe_load(fh) or []
    # Check the file is a list of blocks, each starting with '- dataset_prefix:'.
    if not isinstance(blocks, list):
        raise ValueError(f"{path}: expected a list of BCF target blocks, each starting with '- dataset_prefix: ...'")
    # Start empty lookups for specimens and callgroup files per target, run tables per dataset_prefix, and an empty list of problems found.
    group_to_specimens, callgroup_file_for, run_tables, errors = {}, {}, {}, []
    # Go through the blocks one at a time, numbering them from 1 for error messages.
    for n, block in enumerate(blocks, start=1):
        # Describe the block by its number until its key fields have been read.
        where = f"block {n}"
        # Check the block is a set of 'key: value' lines.
        if not isinstance(block, dict):
            errors.append(f"{where}: not a 'key: value' block")
            continue
        # List any required fields that are missing or empty.
        missing = [f for f in BCF_KEY_FIELDS + ("specimens",) if block.get(f) in (None, "", [])]
        if missing:
            errors.append(f"{where}: missing {', '.join(missing)}")
            continue
        # List any key fields whose value YAML did not read as text, such as a purely numeric value.
        not_text = [f for f in BCF_KEY_FIELDS if not isinstance(block[f], str)]
        if not_text:
            errors.append(f"{where}: {', '.join(not_text)} must be text; put the value in quotes")
            continue
        # Build the target's key from its four key fields, with surrounding spaces removed.
        key = tuple(block[f].strip() for f in BCF_KEY_FIELDS)
        # Unpack the key into its four separate values.
        reference_id_value, dataset_prefix_value, bcf_subset_value, bcf_callgroup_value = key
        # Describe the block by its number and target name in later error messages.
        where = f"block {n} ({dataset_prefix_value}.{reference_id_value}.{bcf_subset_value}.{bcf_callgroup_value})"
        # List any key fields containing a dot, slash, space or tab, which would break the output file names.
        bad_chars = [f for f, v in zip(BCF_KEY_FIELDS, key) if any(c in v for c in "./ \t")]
        if bad_chars:
            errors.append(f"{where}: {', '.join(bad_chars)} must not contain '.', '/' or spaces (they are used in file names)")
        # Check whether an earlier block already defined this same target.
        if key in group_to_specimens:
            errors.append(f"{where}: duplicate target; another block has the same dataset_prefix, reference_id, bcf_subset and bcf_callgroup")
            continue
        # Check the reference_id is listed in the reference metadata table.
        if reference_id_value not in references_table.index:
            errors.append(f"{where}: reference_id not in references_metadata.tsv")
        # Otherwise check the reference_id has intervals defined.
        elif reference_id_value not in intervals_by_reference:
            errors.append(f"{where}: reference_id not in references_intervals.tsv")
        # Take the block's specimens entry.
        specimens = block["specimens"]
        # Split the specimens on spaces and line breaks if they were written as one piece of text rather than a list.
        specimens = specimens.split() if isinstance(specimens, str) else specimens
        # Check the specimens are a list of text IDs.
        if not isinstance(specimens, list) or not all(isinstance(x, str) for x in specimens):
            errors.append(f"{where}: specimens must be text IDs, one per line; quote any purely numeric ID")
            continue
        # Remove surrounding spaces from each specimen ID and drop any empty ones.
        specimens = [x.strip() for x in specimens if x.strip()]
        # Find any specimen IDs listed more than once in this block.
        repeated = sorted({x for x in specimens if specimens.count(x) > 1})
        # Check whether any specimen IDs are repeated.
        if repeated:
            errors.append(f"{where}: specimens listed more than once: {', '.join(repeated)}")
        # Check whether this dataset_prefix's run metadata file has not been read yet.
        if dataset_prefix_value not in run_tables:
            # Work out the run metadata file path for this dataset_prefix.
            run_path = run_metadata_path(dataset_prefix_value)
            # Check whether the run metadata file exists.
            if not os.path.isfile(run_path):
                run_tables[dataset_prefix_value] = None
                errors.append(f"{where}: run metadata file not found: {run_path}")
            # Otherwise the run metadata file exists.
            else:
                # Read the run metadata file, keeping every value as text.
                table = pd.read_csv(run_path, sep="\t", dtype=str)
                # List any required columns the file lacks.
                absent = [c for c in RUN_METADATA_COLUMNS if c not in table.columns]
                # Store the table, or nothing if it lacks required columns.
                run_tables[dataset_prefix_value] = None if absent else table
                # Check whether any required columns are missing.
                if absent:
                    errors.append(f"{where}: {run_path} is missing column(s): {', '.join(absent)}")
        # Check the specimens against the run metadata, if a usable run table exists for this dataset_prefix.
        if run_tables[dataset_prefix_value] is not None:
            known_specimens = set(run_tables[dataset_prefix_value]["specimen_id"])
            unknown = [x for x in specimens if x not in known_specimens]
            # Check whether any specimens are unknown.
            if unknown:
                errors.append(f"{where}: specimens not in {run_metadata_path(dataset_prefix_value)}: {', '.join(unknown)}")
        # Work out the callgroup file path for this target's bcf_callgroup.
        callgroup_path = callgroup_file_path(bcf_callgroup_value)
        # Check whether the callgroup file exists.
        if not os.path.isfile(callgroup_path):
            errors.append(f"{where}: callgroup file not found: {callgroup_path}")
        # Store the target's specimen list.
        group_to_specimens[key] = specimens
        # Store the target's callgroup file.
        callgroup_file_for[key] = callgroup_path
    # Check whether any problems were found across all blocks.
    if errors:
        raise ValueError(f"Problems in {path}:\n  " + "\n  ".join(errors))
    # Return the lookups of specimens and callgroup files per target, and run tables per dataset_prefix.
    return group_to_specimens, callgroup_file_for, run_tables
# Load the BCF targets, giving lookups from (reference_id, dataset_prefix, bcf_subset, bcf_callgroup) to its specimens and callgroup file, and from dataset_prefix to its run metadata table.
bcf_group_to_specimens, bcf_callgroup_file_for, run_metadata_for_prefix = load_bcf_targets(BCF_TARGETS_FILE)
# Combine the run metadata tables of all dataset_prefixes into one table of distinct runs, for the per-run lookups used by A1 and B1.
run_metadata_table = pd.concat(list(run_metadata_for_prefix.values()), ignore_index=True)[list(RUN_METADATA_COLUMNS)].drop_duplicates()
# Find any sequence_id still listed more than once, meaning different files give it a different specimen_id, library_id, or path.
conflicting_runs = sorted(run_metadata_table.loc[run_metadata_table["sequence_id"].duplicated(), "sequence_id"].unique())
# Check whether any sequence_id has conflicting details, which matters because each run is downloaded once to one file.
if conflicting_runs:
    raise ValueError("sequence_id(s) with different specimen_id, library_id, or path in different run metadata files: " + ", ".join(conflicting_runs))
# Import the run metadata file.
run_metadata_table = run_metadata_table.set_index("sequence_id")
sequence_id = run_metadata_table.index
specimen_id = run_metadata_table["specimen_id"]
library_id = run_metadata_table["library_id"]
iRODS_cram_path = run_metadata_table["path"]
# Sort the BCF target keys into a fixed order.
BCF_GROUPS = sorted(bcf_group_to_specimens.keys())
# List each distinct (specimen_id, reference_id, dataset_prefix) combination that needs a processed CRAM, in first-seen order.
SPECIMEN_CRAM_KEYS = list(dict.fromkeys((specimen_id_value, key[0], key[1]) for key, specimen_ids in bcf_group_to_specimens.items() for specimen_id_value in specimen_ids))
# List the distinct dataset_prefixes used in the BCF targets, for wildcard constraints.
dataset_prefixes = sorted({key[1] for key in bcf_group_to_specimens})
# Make a lookup from (dataset_prefix, specimen_id) to the list of that specimen's sequence_ids in that dataset_prefix's run metadata file.
specimen_sequences = {(dataset_prefix_value, specimen_id_value): list(sequence_ids)
                      # Go through each dataset_prefix and its run metadata table.
                      for dataset_prefix_value, table in run_metadata_for_prefix.items()
                      # Group each table's sequence_ids by specimen_id, keeping the file's order.
                      for specimen_id_value, sequence_ids in table.groupby("specimen_id", sort=False)["sequence_id"]}
# List the distinct specimen_ids that need processed CRAMs, for wildcard constraints.
SAMPLES = sorted({specimen_id_value for specimen_id_value, _, _ in SPECIMEN_CRAM_KEYS})


### DETERMINE OUTPUT TARGETS
# Define a function giving the path of the mapped CRAM for one run, produced by B1.
def cram_output_path(specimen_id, sequence_id, reference_id, dataset_prefix):
    return f"{BASE_DIR}/02_reads/B_reads_aligned/{dataset_prefix}/{reference_id}/cram/{dataset_prefix}.{specimen_id}.{sequence_id}.{reference_id}.cram"

# Define a function giving the path of the merged, duplicate-removed CRAM for one specimen, produced by C1.
def cram_rm_output_path(specimen_id, reference_id, dataset_prefix):
    return f"{BASE_DIR}/02_reads/B_reads_aligned/{dataset_prefix}/{reference_id}/cram_rm/{dataset_prefix}.{specimen_id}.{reference_id}.rm.cram"

# Define a function giving the path of the BCF for one target and interval, produced by D1.
def bcf_output_path(reference_id, dataset_prefix, bcf_subset, bcf_callgroup, interval):
    return f"{BASE_DIR}/03_bcf/{dataset_prefix}_group/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}.{interval}.bcf.gz"
# Start an empty lookup from (specimen_id, reference_id, dataset_prefix) to its list of per-run CRAM paths.
cram_files = defaultdict(list)
# Go through each specimen, reference and dataset_prefix combination that needs a processed CRAM.
for specimen_id_value, reference_id_value, dataset_prefix_value in SPECIMEN_CRAM_KEYS:
    # Go through each of that specimen's runs in that dataset_prefix.
    for sequence_id_value in specimen_sequences[(dataset_prefix_value, specimen_id_value)]:
        # Add the run's mapped CRAM path to the combination's list.
        cram_files[(specimen_id_value, reference_id_value, dataset_prefix_value)].append(cram_output_path(specimen_id_value, sequence_id_value, reference_id_value, dataset_prefix_value))
# List the processed CRAM paths for every combination, used when the pipeline is set to stop at processed CRAMs.
cram_rm_targets = [cram_rm_output_path(s, r, d) for (s, r, d) in sorted(cram_files)]
# Start an empty list of final BCF paths, one per target and interval.
bcf_targets = []
# Go through each BCF target.
for (reference_id_value, dataset_prefix_value, bcf_subset_value, bcf_callgroup_value) in BCF_GROUPS:
    # Go through each interval of the target's reference.
    for interval_value in intervals_by_reference[reference_id_value]:
        # Add the BCF path for this target and interval.
        bcf_targets.append(bcf_output_path(reference_id_value, dataset_prefix_value, bcf_subset_value, bcf_callgroup_value, interval_value))
# Start an empty list of final BCF stats paths.
bcf_stats_targets = []
# Go through each BCF target.
for (reference_id_value, dataset_prefix_value, bcf_subset_value, bcf_callgroup_value) in BCF_GROUPS:
    # Go through each interval of the target's reference.
    for interval_value in intervals_by_reference[reference_id_value]:
        # Build the path prefix for all-sites stats thinned to one site per 50 kb, which are always produced (rule E1).
        stats_prefix_allsites = f"{BASE_DIR}/03_bcf/{dataset_prefix_value}_group/{dataset_prefix_value}.{reference_id_value}.{bcf_subset_value}.{bcf_callgroup_value}/stats_allsites/{dataset_prefix_value}.{reference_id_value}.{bcf_subset_value}.{bcf_callgroup_value}.thin50000.{interval_value}"
        # Add the all-sites stats file.
        bcf_stats_targets.append(stats_prefix_allsites + ".stats")
        # Build the path prefix for all-sites stats outside RepeatModeler2 masked regions, which are always produced (rule E2).
        stats_prefix_allsites_maskedrm2 = f"{BASE_DIR}/03_bcf/{dataset_prefix_value}_group/{dataset_prefix_value}.{reference_id_value}.{bcf_subset_value}.{bcf_callgroup_value}/stats_allsites_maskedrm2/{dataset_prefix_value}.{reference_id_value}.{bcf_subset_value}.{bcf_callgroup_value}.thin50000_maskedrm2.{interval_value}"
        # Add the masked all-sites stats file.
        bcf_stats_targets.append(stats_prefix_allsites_maskedrm2 + ".stats")
        # Check whether the reference has a MetaEuk BUSCO exons BED file.
        if busco_exons_metaeuk_bed_for.get(reference_id_value) not in [None, "NULL", ""]:
            # Build the path prefix for stats on MetaEuk BUSCO exons thinned to one site per 20 bp (rule E5).
            stats_prefix_busco_exons_metaeuk = f"{BASE_DIR}/03_bcf/{dataset_prefix_value}_group/{dataset_prefix_value}.{reference_id_value}.{bcf_subset_value}.{bcf_callgroup_value}/stats_busco_exons_metaeuk/{dataset_prefix_value}.{reference_id_value}.{bcf_subset_value}.{bcf_callgroup_value}.thin20_busco_exons_metaeuk.{interval_value}"
            # Add the MetaEuk BUSCO exons stats file.
            bcf_stats_targets.append(stats_prefix_busco_exons_metaeuk + ".stats")
        # Check whether the reference has a Helixer BUSCO exons BED file.
        if busco_exons_helixer_bed_for.get(reference_id_value) not in [None, "NULL", ""]:
            # Build the path prefix for stats on Helixer BUSCO exons thinned to one site per 20 bp (rule E3).
            stats_prefix_busco_exons_helixer = f"{BASE_DIR}/03_bcf/{dataset_prefix_value}_group/{dataset_prefix_value}.{reference_id_value}.{bcf_subset_value}.{bcf_callgroup_value}/stats_busco_exons_helixer/{dataset_prefix_value}.{reference_id_value}.{bcf_subset_value}.{bcf_callgroup_value}.thin20_busco_exons_helixer.{interval_value}"
            # Add the Helixer BUSCO exons stats file.
            bcf_stats_targets.append(stats_prefix_busco_exons_helixer + ".stats")
        # Check whether the reference has a Braker3 BUSCO exons BED file.
        if busco_exons_braker3_bed_for.get(reference_id_value) not in [None, "NULL", ""]:
            # Build the path prefix for stats on Braker3 BUSCO exons thinned to one site per 20 bp (rule E4).
            stats_prefix_busco_exons_braker3 = f"{BASE_DIR}/03_bcf/{dataset_prefix_value}_group/{dataset_prefix_value}.{reference_id_value}.{bcf_subset_value}.{bcf_callgroup_value}/stats_busco_exons_braker3/{dataset_prefix_value}.{reference_id_value}.{bcf_subset_value}.{bcf_callgroup_value}.thin20_busco_exons_braker3.{interval_value}"
            # Add the Braker3 BUSCO exons stats file.
            bcf_stats_targets.append(stats_prefix_busco_exons_braker3 + ".stats")


### LOAD RESOURCES INFORMATION
# Set the YAML file giving memory and queue for B1/C1 per dataset_prefix and for D1/E per dataset_prefix and bcf_subset; rule A1 keeps its own fixed resources.
BCF_RESOURCES_FILE = "../scripts/reads_to_groupBCF.resources.yaml"
# Name the settings each entry in the 'cram' and 'bcf' sections of the resources file must give.
CRAM_RESOURCE_SETTINGS = ("B1_filter_map_memory", "B1_filter_map_queue", "C1_process_cram_memory", "C1_process_cram_queue")
BCF_RESOURCE_SETTINGS = ("D1_call_variants_memory", "D1_call_variants_queue", "E_bcf_stats_memory", "E_bcf_stats_queue")
# Define a function that turns one resource value (a lambda taking wildcards and attempt, a number of MB, or a queue name) into a function of wildcards and attempt.
def compile_resource_setting(value):
    # Check the value is a number or non-empty text, and not true/false.
    if isinstance(value, bool) or not isinstance(value, (int, float, str)) or str(value).strip() == "":
        raise ValueError(f"must be a lambda, a number or a queue name, not {value!r}")
    # Check whether the value is written as a lambda.
    if isinstance(value, str) and value.strip().startswith("lambda"):
        # Turn the lambda text into a real function, able to use anything defined in this snakefile such as n_runs_for.
        func = eval(value.strip(), globals())
        # Check the text really produced a function.
        if not callable(func):
            raise ValueError(f"{value!r} is not a function")
        return func
    # Remove surrounding spaces from a queue name, leaving numbers unchanged.
    value = value.strip() if isinstance(value, str) else value
    return lambda wildcards, attempt: value
# Define a function that loads one section of the resources file and returns a lookup from each key combination to its settings, where any key field may list several values separated by commas.
def load_resource_section(entries, section, key_fields, settings, errors):
    # Check the section is a list of entries.
    if not isinstance(entries, list):
        errors.append(f"section '{section}' missing or not a list of '- dataset_prefix: ...' entries")
        return {}
    table = {}
    # Go through the entries one at a time, numbering them from 1 for error messages.
    for n, entry in enumerate(entries, start=1):
        # Describe the entry by its section and number until its key fields have been read.
        where = f"{section} entry {n}"
        # Check the entry is a set of 'key: value' lines.
        if not isinstance(entry, dict):
            errors.append(f"{where}: not a 'key: value' block")
            continue
        # Start an empty lookup from key field to its listed values, and an empty list of unreadable key fields.
        values_for, bad_keys = {}, []
        # Go through the key fields this section needs.
        for f in key_fields:
            # Take the entry's value for this key field.
            raw = entry.get(f)
            # Split text on commas, keep a YAML list as it is, and mark anything else as unreadable.
            parts = raw.split(",") if isinstance(raw, str) else raw if isinstance(raw, list) else None
            # Check the value was readable, contains only text, and has at least one non-empty value.
            if parts is None or not all(isinstance(x, str) for x in parts) or not [x for x in parts if x.strip()]:
                bad_keys.append(f)
            else:
                # Store the values with surrounding spaces and empty values removed, keeping each value only once in its original order.
                values_for[f] = list(dict.fromkeys(x.strip() for x in parts if x.strip()))
        # Check whether any key fields were unreadable.
        if bad_keys:
            errors.append(f"{where}: missing {', '.join(bad_keys)} (must be text, or several separated by commas; quote any purely numeric value)")
            continue
        # Start the list of key combinations with a single empty combination.
        keys = [()]
        # Go through the key fields in order.
        for f in key_fields:
            # Extend every combination so far with each listed value of this key field.
            keys = [k + (v,) for k in keys for v in values_for[f]]
        # Describe the entry by its section, number and listed key values in later error messages.
        where = f"{section} entry {n} ({'.'.join(', '.join(values_for[f]) for f in key_fields)})"
        # Find any of this entry's key combinations already covered by an earlier entry.
        repeated = [k for k in keys if k in table]
        # Check whether any combinations are already covered.
        if repeated:
            errors.append(f"{where}: duplicate; another entry already covers {'; '.join('.'.join(k) for k in repeated)}")
            continue
        # List any settings in the entry that are neither key fields nor expected settings, such as typos.
        unknown = [str(k) for k in entry if k not in key_fields + settings]
        if unknown:
            errors.append(f"{where}: unknown settings: {', '.join(unknown)} (expected: {', '.join(settings)})")
        # List any expected settings the entry does not give.
        missing = [s for s in settings if s not in entry]
        if missing:
            errors.append(f"{where}: missing settings: {', '.join(missing)}")
        # Start an empty lookup from setting name to its function.
        compiled = {}
        # Go through the expected settings.
        for s in settings:
            # Only handle settings the entry actually gives.
            if s in entry:
                # Try to turn the setting's value into a function.
                try:
                    # Store the function for this setting.
                    compiled[s] = compile_resource_setting(entry[s])
                # Catch any problem with the value, such as a lambda with a syntax error.
                except Exception as e:
                    errors.append(f"{where}: {s}: {e}")
        # Go through every key combination this entry covers.
        for key in keys:
            table[key] = compiled
    # Return the lookup from key combination to settings.
    return table
# Define a function that loads and checks the resources file, making sure it covers every dataset_prefix and (dataset_prefix, bcf_subset) in the BCF targets; extra entries are allowed.
def load_bcf_resources(path):
    # Check the resources file exists.
    if not os.path.isfile(path):
        raise ValueError(f"resources file not found: {path}")
    # Open the resources file.
    with open(path) as fh:
        # Read the YAML, treating an empty file as having no sections.
        sections = yaml.safe_load(fh) or {}
    # Check the file is laid out as named sections.
    if not isinstance(sections, dict):
        raise ValueError(f"{path}: expected two sections, 'cram:' and 'bcf:'")
    # Start an empty list of problems found.
    errors = []
    # List any sections other than 'cram' and 'bcf'.
    unknown_sections = [str(k) for k in sections if k not in ("cram", "bcf")]
    # Check whether any unknown sections were given.
    if unknown_sections:
        errors.append(f"unknown sections: {', '.join(unknown_sections)} (expected: cram, bcf)")
    # Load the 'cram' section, keyed by dataset_prefix alone.
    cram_table = load_resource_section(sections.get("cram"), "cram", ("dataset_prefix",), CRAM_RESOURCE_SETTINGS, errors)
    # Load the 'bcf' section, keyed by dataset_prefix and bcf_subset.
    bcf_table = load_resource_section(sections.get("bcf"), "bcf", ("dataset_prefix", "bcf_subset"), BCF_RESOURCE_SETTINGS, errors)
    # Go through every dataset_prefix used in the BCF targets.
    for dataset_prefix_value in sorted({key[1] for key in bcf_group_to_specimens}):
        # Check the 'cram' section has an entry for it.
        if (dataset_prefix_value,) not in cram_table:
            errors.append(f"no 'cram' entry for dataset_prefix {dataset_prefix_value} (used in {BCF_TARGETS_FILE})")
    # Go through every (dataset_prefix, bcf_subset) combination used in the BCF targets.
    for dataset_prefix_value, bcf_subset_value in sorted({(key[1], key[2]) for key in bcf_group_to_specimens}):
        # Check the 'bcf' section has an entry for it.
        if (dataset_prefix_value, bcf_subset_value) not in bcf_table:
            errors.append(f"no 'bcf' entry for dataset_prefix {dataset_prefix_value} with bcf_subset {bcf_subset_value} (used in {BCF_TARGETS_FILE})")
    # Check whether any problems were found.
    if errors:
        raise ValueError(f"Problems in {path}:\n  " + "\n  ".join(errors))
    return cram_table, bcf_table
# Load the resources file, giving lookups from dataset_prefix to B1/C1 settings and from (dataset_prefix, bcf_subset) to D1/E settings.
cram_resources_for, bcf_resources_for = load_bcf_resources(BCF_RESOURCES_FILE)
# Define a function that gives a rule the resource function for one setting from the resources file.
def yaml_resource(setting):
    # Define the function snakemake calls for each job, with the job's wildcards and attempt number.
    def resource(wildcards, attempt):
        # Check whether the setting belongs to B1 or C1.
        if setting in CRAM_RESOURCE_SETTINGS:
            # Look up the B1/C1 settings for the job's dataset_prefix.
            settings = cram_resources_for[(wildcards.dataset_prefix,)]
        # Otherwise the setting belongs to D1 or E.
        else:
            # Look up the D1/E settings for the job's dataset_prefix and bcf_subset.
            settings = bcf_resources_for[(wildcards.dataset_prefix, wildcards.bcf_subset)]
        # Work out the setting's value for this job and attempt.
        return settings[setting](wildcards, attempt)
    # Return the function for snakemake to use.
    return resource


### MISCELLANEOUS FUNCTIONS
# Define a function that builds a wildcard constraint matching exactly one of the given values, or nothing at all if the list is empty.
def one_of(values):
    # Join the escaped values with '|' inside brackets, or return '(?!)' (a pattern that never matches) when there are no values.
    return "(" + "|".join(re.escape(v) for v in values) + ")" if values else "(?!)"
# Define a function giving the number of sequencing runs for a C1 job's (specimen_id, reference_id, dataset_prefix).
def n_runs_for(w):
    return len(cram_files[(w.specimen_id, w.reference_id, w.dataset_prefix)])





# Tell snakemake what the final file output will be so that it knows the directionality of the pipeline. Hash out all but one of the options below.
rule all:
    # To finish with all processed cram files:
    # input: cram_rm_targets
    # To finish with BCF files:
    # input: bcf_targets
    # To finish with BCF files and stats:
    input: bcf_targets + bcf_stats_targets




# Step A: Download and convert interleaved cram files from iRODS.
rule A1_iRODS_cram_download:
    params:
        iRODS_cram_path = lambda wildcards: iRODS_cram_path[wildcards.sequence_id]
    output:
        cram_reads=temp(BASE_DIR + "/02_reads/A_reads/{specimen_id}.{sequence_id}.cram")
    resources:
        iRODS_jobs=1,
        mem_mb=lambda wildcards, attempt: 1000 * attempt,
        queue=lambda wildcards, attempt: "normal" if attempt == 1 else "normal"
    threads: 1
    shell:
        """
        iget {params.iRODS_cram_path} {output.cram_reads} -f
        """




### Step B: Filter filters and map them to reference genome.
## samtools fastq options:
# -n: no need to add /1 /2 to reads since interleaved already
## fastp filtering options:
# --length_required 50: minimum read length of 50 bp
# --trim_poly_g: trim poly-G tails (because NovaSeq data)
# --detect_adapter_for_pe: detect adapters for paired-end reads
## minibwa mapping options:
# -R: read group header line with sample, library, platform, and ID information
rule B1_filter_map:
    wildcard_constraints:
        dataset_prefix = one_of(dataset_prefixes),
        reference_id = one_of(REFERENCE_IDS)
    input:
        cram=BASE_DIR + "/02_reads/A_reads/{specimen_id}.{sequence_id}.cram"
    output:
        fastp_html=BASE_DIR + "/02_reads/B_reads_aligned/{dataset_prefix}/{reference_id}/cram/{dataset_prefix}.{specimen_id}.{sequence_id}.{reference_id}.fastp.html",
        fastp_json=BASE_DIR + "/02_reads/B_reads_aligned/{dataset_prefix}/{reference_id}/cram/{dataset_prefix}.{specimen_id}.{sequence_id}.{reference_id}.fastp.json",
        cram=temp(BASE_DIR + "/02_reads/B_reads_aligned/{dataset_prefix}/{reference_id}/cram/{dataset_prefix}.{specimen_id}.{sequence_id}.{reference_id}.cram"),
        cram_stats=BASE_DIR + "/02_reads/B_reads_aligned/{dataset_prefix}/{reference_id}/cram/{dataset_prefix}.{specimen_id}.{sequence_id}.{reference_id}.cram.stats",
        cram_flagstat=BASE_DIR + "/02_reads/B_reads_aligned/{dataset_prefix}/{reference_id}/cram/{dataset_prefix}.{specimen_id}.{sequence_id}.{reference_id}.cram.flagstat"
    params:
        sequence_id = "{sequence_id}",
        specimen_id = "{specimen_id}",
        library_id = lambda wildcards: library_id[wildcards.sequence_id],
        reference_path = lambda wildcards: reference_path_for[wildcards.reference_id]
    resources:
        mem_mb=yaml_resource("B1_filter_map_memory"),
        queue=yaml_resource("B1_filter_map_queue")
    threads: 20
    shell:
        """
        {samtools} fastq -@ {threads} {input.cram} -n | \
            {fastp} -w {threads} --stdin --stdout --interleaved_in --length_required 50 --trim_poly_g --detect_adapter_for_pe -j {output.fastp_json} -h {output.fastp_html} | \
            {minibwa} map -t {threads} {params.reference_path} - -R "@RG\\tSM:{wildcards.specimen_id}\\tLB:{params.library_id}\\tPL:ILLUMINA\\tID:{params.specimen_id}.{params.library_id}" | \
            {samtools} view -@ {threads} -o {output.cram} - --reference {params.reference_path}
        {samtools} stats -@ {threads} --reference {params.reference_path} {output.cram} > {output.cram_stats}
        {samtools} flagstat -@ {threads} {output.cram} > {output.cram_flagstat}
        """




# Step C: Remove duplicates (optical and PCR) and unmapped reads from the CRAM files; index and get stats.
rule C1_process_cram:
    wildcard_constraints:
        specimen_id = one_of(SAMPLES),
        reference_id = one_of(REFERENCE_IDS),
        dataset_prefix = one_of(dataset_prefixes)
    input:
        cram_list = lambda w: cram_files[(w.specimen_id, w.reference_id, w.dataset_prefix)]
    output:
        cram_rm=BASE_DIR + "/02_reads/B_reads_aligned/{dataset_prefix}/{reference_id}/cram_rm/{dataset_prefix}.{specimen_id}.{reference_id}.rm.cram",
        index=BASE_DIR + "/02_reads/B_reads_aligned/{dataset_prefix}/{reference_id}/cram_rm/{dataset_prefix}.{specimen_id}.{reference_id}.rm.cram.crai",
        markdup_stats=BASE_DIR + "/02_reads/B_reads_aligned/{dataset_prefix}/{reference_id}/cram_rm/{dataset_prefix}.{specimen_id}.{reference_id}.rm.cram.samtools_markdup_stats",
        stats=BASE_DIR + "/02_reads/B_reads_aligned/{dataset_prefix}/{reference_id}/cram_rm/{dataset_prefix}.{specimen_id}.{reference_id}.rm.cram.stats",
        flagstat=BASE_DIR + "/02_reads/B_reads_aligned/{dataset_prefix}/{reference_id}/cram_rm/{dataset_prefix}.{specimen_id}.{reference_id}.rm.cram.flagstat"
    params:
        reference_path=lambda w: reference_path_for[w.reference_id],
        n_runs=n_runs_for,
        tmp_directory=lambda w: f"{BASE_DIR}/02_reads/B_reads_aligned/{w.dataset_prefix}/{w.reference_id}/cram_rm/tmp/{w.dataset_prefix}.{w.specimen_id}.{w.reference_id}"
    resources:
        mem_mb=yaml_resource("C1_process_cram_memory"),
        queue=yaml_resource("C1_process_cram_queue")
    threads: 10
    shell:
        """
        mkdir -p {params.tmp_directory}
        TMP=$(mktemp -d {params.tmp_directory}/XXXXXX)
        trap 'rm -rf "$TMP"' EXIT

        if [ {params.n_runs} -eq 1 ]; then
            {samtools} view -@ {threads} --reference {params.reference_path} -u -o - {input.cram_list}
        else
            {samtools} merge -@ {threads} --reference {params.reference_path} -u -o - {input.cram_list}
        fi | \
            {samtools} collate -@ {threads} --reference {params.reference_path} -T $TMP/collate -O -u - | \
            {samtools} fixmate -@ {threads} --reference {params.reference_path} -m -u - - | \
            {samtools} sort -@ {threads} --reference {params.reference_path} -T $TMP/sort -u - | \
            {samtools} markdup -@ {threads} --reference {params.reference_path} --use-read-groups -d 2500 -T $TMP/markdup -r -f {output.markdup_stats} --json -u - - | \
            {samtools} view -@ {threads} --reference {params.reference_path} -F 4 -o {output.cram_rm} -

        {samtools} index {output.cram_rm} -o {output.index}
        {samtools} stats -@ {threads} --reference {params.reference_path} {output.cram_rm} > {output.stats}
        {samtools} flagstat -@ {threads} {output.cram_rm} > {output.flagstat}
        """






### Step D: Call variants per chromosome interval.
## bcftools mpileup options:
# -q 20: minimum mapping quality of 20
# -Q 20: minimum base quality of 20
# -d & -L 250000: essentially disables per-file per-position read depth cap so all reads are considered and depth filtering done later
# -a FORMAT/QS,FORMAT/SP,FORMAT/DP,FORMAT/AD,FORMAT/SCR,INFO/AD,INFO/SCR: annotate with various stats
# -f: reference genome
# -r: region to call variants in
# -Ou: output uncompressed BCF to stdout for piping
## bcftools call options:
# -m: use the multiallelic caller
# -f GQ,GP: annotate with genotype quality and genotype posterior probability
# -G: a file containing samples and the 'calling group' they belong to, so that calling is sensitive to population/species structure
# --ploidy-file: per-chromosome per-sex ploidy definitions (handles X/Y haploid calling for males)
# -S: file containing all samples and their sex, used in conjunction with --ploidy-file
# -Ob: output compressed BCF
rule D1_call_variants:
    wildcard_constraints:
        reference_id = "(" + "|".join(REFERENCE_IDS) + ")" if REFERENCE_IDS else ".+"
    input:
        cram=lambda w: [cram_rm_output_path(specimen_id, w.reference_id, w.dataset_prefix) for specimen_id in bcf_group_to_specimens[(w.reference_id, w.dataset_prefix, w.bcf_subset, w.bcf_callgroup)]],
        crai=lambda w: [cram_rm_output_path(specimen_id, w.reference_id, w.dataset_prefix) + ".crai" for specimen_id in bcf_group_to_specimens[(w.reference_id, w.dataset_prefix, w.bcf_subset, w.bcf_callgroup)]]
    output:
        bcf=BASE_DIR + "/03_bcf/{dataset_prefix}_group/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}.{interval}.bcf.gz",
        csi=BASE_DIR + "/03_bcf/{dataset_prefix}_group/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}.{interval}.bcf.gz.csi",
        specimen_sex_subset=temp(BASE_DIR + "/03_bcf/{dataset_prefix}_group/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}.{interval}.sex_subset.tsv")
    params:
        coordinates=lambda w: interval_coordinates[(w.reference_id, w.interval)],
        reference_path=lambda w: reference_path_for[w.reference_id],
        specimens=lambda w: bcf_group_to_specimens[(w.reference_id, w.dataset_prefix, w.bcf_subset, w.bcf_callgroup)],
        callgroup_file=lambda w: bcf_callgroup_file_for[(w.reference_id, w.dataset_prefix, w.bcf_subset, w.bcf_callgroup)],
        specimen_sex_file=specimen_sex_file,
        ploidy_file=ploidy_file
    resources:
        mem_mb=yaml_resource("D1_call_variants_memory"),
        queue=yaml_resource("D1_call_variants_queue")
    threads: 1
    shell:
        """
        # Filter sex file to only include samples in this BCF group
        awk 'NR==FNR{{a[$1]=1;next}} $1 in a{{print}}' <(printf '%s\n' {params.specimens}) {params.specimen_sex_file} > {output.specimen_sex_subset}
        # Run calling
        {bcftools} mpileup -q 20 -Q 20 -d 250000 -L 250000 -a FORMAT/QS,FORMAT/SP,FORMAT/DP,FORMAT/AD,FORMAT/SCR,INFO/AD,INFO/SCR,INFO/FS -f {params.reference_path} -r {params.coordinates} {input.cram} -Ou | \
            {bcftools} call -m -f GQ,GP -G {params.callgroup_file} --ploidy-file {params.ploidy_file} -S {output.specimen_sex_subset} -Ou | \
            {bcftools} +fill-tags -Ob -o {output.bcf} -- -t AC,AN,HWE,ExcHet,F_MISSING
        # Index output
        {bcftools} index --csi {output.bcf}
        """
    



### Step E: Calculate statistics on the base BCF files.

# E1: Stats on all sites (always runs)
rule E1_bcf_stats_allsites:
    wildcard_constraints:
        reference_id = "(" + "|".join(REFERENCE_IDS) + ")" if REFERENCE_IDS else ".+"
    input:
        bcf=BASE_DIR + "/03_bcf/{dataset_prefix}_group/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}.{interval}.bcf.gz",
        csi=BASE_DIR + "/03_bcf/{dataset_prefix}_group/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}.{interval}.bcf.gz.csi"
    output:
        vcf_thinned=temp(BASE_DIR + "/03_bcf/{dataset_prefix}_group/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}/stats_allsites/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}.thin50000.{interval}.vcf.gz"),
        stats=BASE_DIR + "/03_bcf/{dataset_prefix}_group/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}/stats_allsites/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}.thin50000.{interval}.stats"
    resources:
        mem_mb=yaml_resource("E_bcf_stats_memory"),
        queue=yaml_resource("E_bcf_stats_queue")
    threads: 1
    shell:
        """
        {bcftools} view -i 'TYPE!="indel"' -Ou {input.bcf} | {bcftools} +prune -n 1 -w 50000bp -N rand -Oz -o {output.vcf_thinned}
        {bcftools} query -f '%CHROM\t%POS\t%ALT\t%QUAL\t%F_MISSING\t%FS[\t%DP:%GQ:%SP]\n' {output.vcf_thinned} | cat <({bcftools} query -l {output.vcf_thinned} | tr '\n' '\t' | sed 's/^/CHROM\tPOS\tALT\tQUAL\tF_MISSING\tFS\t/;s/\t$//;s/$/\\n/') - > {output.stats}
        """

# E2: Stats on rm2 masked all sites (always runs)
rule E2_bcf_stats_allsites_maskedrm2:
    wildcard_constraints:
        reference_id = "(" + "|".join(REFERENCE_IDS) + ")" if REFERENCE_IDS else ".+"
    input:
        bcf=BASE_DIR + "/03_bcf/{dataset_prefix}_group/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}.{interval}.bcf.gz",
        csi=BASE_DIR + "/03_bcf/{dataset_prefix}_group/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}.{interval}.bcf.gz.csi"
    output:
        vcf_thinned=temp(BASE_DIR + "/03_bcf/{dataset_prefix}_group/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}/stats_allsites_maskedrm2/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}.thin50000_maskedrm2.{interval}.vcf.gz"),
        stats=BASE_DIR + "/03_bcf/{dataset_prefix}_group/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}/stats_allsites_maskedrm2/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}.thin50000_maskedrm2.{interval}.stats"
    params:
        allsites_masked_rm2_bed=lambda w: masked_rm2_bed_for[w.reference_id]
    resources:
        mem_mb=yaml_resource("E_bcf_stats_memory"),
        queue=yaml_resource("E_bcf_stats_queue")
    threads: 1
    shell:
        """
        {bcftools} view -T ^{params.allsites_masked_rm2_bed} -i 'TYPE!="indel"' -Ou {input.bcf} | {bcftools} +prune -n 1 -w 50000bp -N rand -Oz -o {output.vcf_thinned}
        {bcftools} query -f '%CHROM\t%POS\t%ALT\t%QUAL\t%F_MISSING\t%FS[\t%DP:%GQ:%SP]\n' {output.vcf_thinned} | cat <({bcftools} query -l {output.vcf_thinned} | tr '\n' '\t' | sed 's/^/CHROM\tPOS\tALT\tQUAL\tF_MISSING\tFS\t/;s/\t$//;s/$/\\n/') - > {output.stats}
        """

# E3: Stats on BUSCO exons identified by Helixer (only runs if bed is not NULL)
rule E3_bcf_stats_busco_helixer:
    wildcard_constraints:
        reference_id = "(" + "|".join([r for r in REFERENCE_IDS if busco_exons_helixer_bed_for.get(r) not in [None, "NULL", ""]]) + ")" if REFERENCE_IDS else ".+"
    input:
        bcf=BASE_DIR + "/03_bcf/{dataset_prefix}_group/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}.{interval}.bcf.gz",
        csi=BASE_DIR + "/03_bcf/{dataset_prefix}_group/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}.{interval}.bcf.gz.csi"
    output:
        vcf_thinned=temp(BASE_DIR + "/03_bcf/{dataset_prefix}_group/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}/stats_busco_exons_helixer/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}.thin20_busco_exons_helixer.{interval}.vcf.gz"),
        stats=BASE_DIR + "/03_bcf/{dataset_prefix}_group/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}/stats_busco_exons_helixer/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}.thin20_busco_exons_helixer.{interval}.stats"
    params:
        busco_exons_helixer_bed=lambda w: busco_exons_helixer_bed_for[w.reference_id]
    resources:
        mem_mb=yaml_resource("E_bcf_stats_memory"),
        queue=yaml_resource("E_bcf_stats_queue")
    threads: 1
    shell:
        """
        {bcftools} view -R {params.busco_exons_helixer_bed} -i 'TYPE!="indel"' -Ou {input.bcf} | {bcftools} +prune -n 1 -w 20bp -N rand -Oz -o {output.vcf_thinned}
        {bcftools} query -f '%CHROM\t%POS\t%ALT\t%QUAL\t%F_MISSING\t%FS[\t%DP:%GQ:%SP]\n' {output.vcf_thinned} | cat <({bcftools} query -l {output.vcf_thinned} | tr '\n' '\t' | sed 's/^/CHROM\tPOS\tALT\tQUAL\tF_MISSING\tFS\t/;s/\t$//;s/$/\\n/') - > {output.stats}
        """

# E4: Stats on BUSCO exons identified by Braker3 (only runs if bed is not NULL)
rule E4_bcf_stats_busco_braker3:
    wildcard_constraints:
        reference_id = "(" + "|".join([r for r in REFERENCE_IDS if busco_exons_braker3_bed_for.get(r) not in [None, "NULL", ""]]) + ")" if REFERENCE_IDS else ".+"
    input:
        bcf=BASE_DIR + "/03_bcf/{dataset_prefix}_group/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}.{interval}.bcf.gz",
        csi=BASE_DIR + "/03_bcf/{dataset_prefix}_group/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}.{interval}.bcf.gz.csi"
    output:
        vcf_thinned=temp(BASE_DIR + "/03_bcf/{dataset_prefix}_group/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}/stats_busco_exons_braker3/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}.thin20_busco_exons_braker3.{interval}.vcf.gz"),
        stats=BASE_DIR + "/03_bcf/{dataset_prefix}_group/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}/stats_busco_exons_braker3/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}.thin20_busco_exons_braker3.{interval}.stats"
    params:
        busco_exons_braker3_bed=lambda w: busco_exons_braker3_bed_for[w.reference_id]
    resources:
        mem_mb=yaml_resource("E_bcf_stats_memory"),
        queue=yaml_resource("E_bcf_stats_queue")
    threads: 1
    shell:
        """
        {bcftools} view -R {params.busco_exons_braker3_bed} -i 'TYPE!="indel"' -Ou {input.bcf} | {bcftools} +prune -n 1 -w 20bp -N rand -Oz -o {output.vcf_thinned}
        {bcftools} query -f '%CHROM\t%POS\t%ALT\t%QUAL\t%F_MISSING\t%FS[\t%DP:%GQ:%SP]\n' {output.vcf_thinned} | cat <({bcftools} query -l {output.vcf_thinned} | tr '\n' '\t' | sed 's/^/CHROM\tPOS\tALT\tQUAL\tF_MISSING\tFS\t/;s/\t$//;s/$/\\n/') - > {output.stats}
        """

# E5: Stats on BUSCO exons identified by MetaEuk (only runs if bed is not NULL)
rule E5_bcf_stats_busco_metaeuk:
    wildcard_constraints:
        reference_id = "(" + "|".join([r for r in REFERENCE_IDS if busco_exons_metaeuk_bed_for.get(r) not in [None, "NULL", ""]]) + ")" if REFERENCE_IDS else ".+"
    input:
        bcf=BASE_DIR + "/03_bcf/{dataset_prefix}_group/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}.{interval}.bcf.gz",
        csi=BASE_DIR + "/03_bcf/{dataset_prefix}_group/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}.{interval}.bcf.gz.csi"
    output:
        vcf_thinned=temp(BASE_DIR + "/03_bcf/{dataset_prefix}_group/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}/stats_busco_exons_metaeuk/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}.thin20_busco_exons_metaeuk.{interval}.vcf.gz"),
        stats=BASE_DIR + "/03_bcf/{dataset_prefix}_group/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}/stats_busco_exons_metaeuk/{dataset_prefix}.{reference_id}.{bcf_subset}.{bcf_callgroup}.thin20_busco_exons_metaeuk.{interval}.stats"
    params:
        busco_exons_metaeuk_bed=lambda w: busco_exons_metaeuk_bed_for[w.reference_id]
    resources:
        mem_mb=yaml_resource("E_bcf_stats_memory"),
        queue=yaml_resource("E_bcf_stats_queue")
    threads: 1
    shell:
        """
        {bcftools} view -R {params.busco_exons_metaeuk_bed} -i 'TYPE!="indel"' -Ou {input.bcf} | {bcftools} +prune -n 1 -w 20bp -N rand -Oz -o {output.vcf_thinned}
        {bcftools} query -f '%CHROM\t%POS\t%ALT\t%QUAL\t%F_MISSING\t%FS[\t%DP:%GQ:%SP]\n' {output.vcf_thinned} | cat <({bcftools} query -l {output.vcf_thinned} | tr '\n' '\t' | sed 's/^/CHROM\tPOS\tALT\tQUAL\tF_MISSING\tFS\t/;s/\t$//;s/$/\\n/') - > {output.stats}
        """
