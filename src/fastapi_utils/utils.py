# ONLY include: read_sequences_from_filename, create_sequences_data, get_completed_zip_paths, merge_zip_archives

import json
import logging
import os
import tempfile
import zipfile
from pathlib import Path
from typing import Any, Dict, List
from uuid import uuid4

from shared_utils.schemas import SequencesMapping, SequenceData
from shared_utils.file_utils import get_zip_by_idr_dir, get_zip_by_fasta_dir, get_runs_dir, get_run_metadata

logger = logging.getLogger(__name__)

def create_sequences_data(sequence_strings: set[str]) -> SequencesMapping:
    """Build a mapping of sequence string -> metadata by scanning prior runs.

    If a sequence appears in multiple runs, a cached status is preferred over pending.
    Any sequences not found are initialized as pending with a fresh UUID.
    """
    # Copy to avoid mutating the caller's set
    sequence_strings = set(sequence_strings)
    sequences_data: SequencesMapping = {}
    for path in get_runs_dir().iterdir():
        if path.suffix == ".json" and path.is_file():
            try:
                data = get_run_metadata(path.stem)
                sequences = data["sequences"]
                if not sequences:
                    continue
                for seq_str, seq_data in sequences.items():
                    if seq_str not in sequence_strings:
                        continue
                    sequence_strings.remove(seq_str)

                    seq_status = seq_data["status"]
                    if seq_status == "pending" or seq_status == "pending_external":
                        seq_data["status"] = "pending_external"
                    elif seq_status == "cached" or seq_status == "complete":
                        seq_data["status"] = "cached"
                    else:
                        raise ValueError(f"Invalid sequence status: {seq_status}")

                    if seq_str in sequences_data:
                        # Prefer cached over pending
                        if sequences_data[seq_str]["status"] != "cached":
                            sequences_data[seq_str] = seq_data
                    else:
                        sequences_data[seq_str] = seq_data
            except (json.JSONDecodeError, FileNotFoundError) as e:
                logger.warning(f"Failed to read {path}: {e}")
            except Exception as e:
                logger.warning(f"Unexpected error reading {path}: {e}")

    for seq_str in sequence_strings:
        sequences_data[seq_str] = SequenceData(
            sequence_id="TEMP_ID",
            status="pending",
            start_time=None,
            end_time=None,
            seq_uuid=str(uuid4()),
            zip_path=None,
            job_id=None,
        )
    return sequences_data

def get_completed_zip_paths(
    sequences_data: SequencesMapping, require_all_complete: bool = False
) -> List[str]:
    """Return absolute paths for zipfiles that are available for the provided sequences."""
    completed_zip_paths: List[str] = []
    idr_zip_dir = get_zip_by_idr_dir()
    for seq_data in sequences_data.values():
        if seq_data["zip_path"]:
            completed_zip_paths.append(seq_data["zip_path"])
            continue
        seq_uuid = seq_data["seq_uuid"]
        if seq_uuid:
            candidate = idr_zip_dir / f"{seq_uuid}.zip"
            if candidate.exists():
                completed_zip_paths.append(str(candidate))
        elif require_all_complete:
            raise ValueError(f"Sequence {seq_data['sequence_id']} is not complete")
    return completed_zip_paths


# ------------------------------------------------------------------------------- #


# The following functions are copied from nardini.utils.py
# NOTE: See https://github.com/mshinn23/nardini/blob/main/nardini/utils.py
def read_sequences_from_filename(sequence_filename, default_name, verbose=False):
    """This is a helper function to read in sequences from a sequence FASTA file.
    This is a companion function to `read_sequences_from_string_list`.

    @param sequence_filename (str): The filepath of the file containing sequences.
                                    This file can contain FASTA records, or
                                    sequences that are saved on each newline.

    @param default_name (str):      This is the prefix name that's used when the
                                    file is found to only contain raw sequences
                                    on each new line.

    @param verbose (bool):          The extent to which information should be
                                    displayed throughout the analysis. Default:
                                    False.

    @returns seqio_sequences:       A list of sequence strings that were
                                    extracted from the sequence file.
    """
    from Bio import SeqIO
    from io import StringIO

    seqio_sequences = list()
    if sequence_filename is None:
        raise RuntimeError("This parameter cannot be `None`. Exiting.")

    if type(sequence_filename) is not str:
        raise RuntimeError("This parameter can only be a string. Exiting.")

    if os.path.exists(sequence_filename):
        with open(sequence_filename, "r") as seqfile:
            content = seqfile.read()
            if content.count(">") > 0:
                seqio_sequences = list(SeqIO.parse(sequence_filename, "fasta"))
            else:
                # This means that there is a list of sequences that are separated
                # by newlines. We split by any whitespace so as to capture
                # carriage returns and line feeds.
                raw_sequences = [
                    s.strip() for s in content.split() if len(s.strip()) > 0
                ]
                # seqio_sequences = read_sequences_from_string_list(
                #     raw_sequences, default_name
                # )
                def read_sequences_from_string_list(list_of_sequences, default_name, verbose=False):
                    sequences = list()
                    # This means that we have to create a fake record using the sequence content.
                    for index, sequence in enumerate(list_of_sequences, start=1):
                        fasta = f">{default_name}-{index}\n{sequence}"
                        fake_record = SeqIO.read(StringIO(fasta), "fasta")
                        sequences.append(fake_record)

                    if verbose:
                        print("Number of sequences read: {num}".format(num=len(sequences)), end="\n\n")
                    return sequences
                seqio_sequences = read_sequences_from_string_list(
                    raw_sequences, default_name
                )
    else:
        raise RuntimeError(f'Sequence filename: "{sequence_filename}" not found.')
    return seqio_sequences

def merge_zip_archives(zip_list: List[str], destination_filepath: str) -> str:
    """Merge multiple zip files into a single zip.

    - Copies all files from each input zip into the merged archive, resolving
      filename collisions by appending a suffix.
    - Concatenates all per-zip `sequences.tsv` into a `master_sequences.tsv` in
      the merged archive.
    """
    # helper fxn for mergeZips
    def add_to_master_sequences_tsv(master_tsv_path: str, tsv_content: str) -> None:
        """Append rows from a per-zip sequences.tsv into the master TSV, skipping headers."""
        import re

        lines = tsv_content.splitlines()
        for line in lines:
            line = line.strip()
            if not line:
                continue
            columns = re.split(r" {2,}", line)
            if columns and columns[0] != "ID":
                with open(master_tsv_path, "a") as f:
                    f.write("\t".join(columns) + "\n")


    def _get_unique_filename(original_filename, used_filenames, zip_index):
        """Generate a unique filename by adding a suffix if the original is already used"""
        if original_filename not in used_filenames:
            return original_filename

        # Extract filename and extension
        if "." in original_filename:
            name, ext = original_filename.rsplit(".", 1)
            ext = "." + ext
        else:
            name = original_filename
            ext = ""

        # Try adding zip index first
        candidate = f"{name}_zip{zip_index}{ext}"
        if candidate not in used_filenames:
            return candidate

        # If still duplicate, add a counter
        counter = 1
        while True:
            candidate = f"{name}_zip{zip_index}_{counter}{ext}"
            if candidate not in used_filenames:
                return candidate
            counter += 1

    if not zip_list:
        raise ValueError("No zip files provided to merge")

    destination_filepath = str(destination_filepath)
    if not destination_filepath.endswith(".zip"):
        destination_filepath = destination_filepath + ".zip"

    # Use a temporary file for the master TSV to avoid working directory issues
    with tempfile.NamedTemporaryFile(
        mode="w+", suffix=".tsv", delete=False
    ) as temp_tsv:
        master_tsv_path = temp_tsv.name
        header = [
            "ID",
            "original_seq",
            "most_similar_seq",
            "sum_abs_zscore_original_seq",
            "sum_abs_zscore_scrambled_seq",
        ]
        temp_tsv.write("\t".join(header) + "\n")

    try:
        with zipfile.ZipFile(destination_filepath, "w") as merged_zip:
            used_filenames = set()

            for zip_index, zip_path in enumerate(zip_list):
                with zipfile.ZipFile(zip_path, "r") as zip_ref:
                    for file_info in zip_ref.infolist():
                        file_name = file_info.filename
                        if file_name != "sequences.tsv":
                            unique_filename = _get_unique_filename(
                                file_name, used_filenames, zip_index
                            )
                            used_filenames.add(unique_filename)
                            with zip_ref.open(file_info) as source_file:
                                merged_zip.writestr(unique_filename, source_file.read())
                        else:
                            with zip_ref.open(file_info) as seq_tsv_file:
                                seq_tsv_content = seq_tsv_file.read().decode("utf-8")
                                logger.info(f"Processing sequences.tsv from {zip_path}")
                                add_to_master_sequences_tsv(
                                    master_tsv_path, seq_tsv_content
                                )

            with open(master_tsv_path, "r") as master_tsv_file:
                master_tsv_content = master_tsv_file.read()
                merged_zip.writestr("master_sequences.tsv", master_tsv_content)

    finally:
        # Clean up the temporary file
        path_obj = Path(master_tsv_path)
        if path_obj.exists():
            path_obj.unlink()

    if not Path(destination_filepath).exists():
        raise ValueError("Zip creation failed!")
    logger.info(f"Merged zip file created at {destination_filepath}")
    return destination_filepath


def ensure_volume_directories() -> bool:
    """Ensure the required directories exist in the mounted volume."""
    was_created = False
    for directory in [get_zip_by_fasta_dir(), get_zip_by_idr_dir(), get_runs_dir()]:
        if not directory.exists():
            logger.warning(f"Volume directory {directory} does not exist, creating it")
            directory.mkdir(parents=True, exist_ok=True)
            was_created = True
    return was_created


def sanitize_output_filename(filename: str) -> str:
    if not filename:
        raise ValueError("Filename is required")
    if len(filename) > 250:
        filename = filename[:250]
    # Remove any potential .fasta extension
    suffix = filename.split(".")[-1]
    if suffix in ["fasta", "fa", "fas"]:
        suffix_len = len(suffix) + 1
        filename = filename[:-suffix_len]
    if not filename.endswith(".zip"):
        filename = filename + ".zip"
    return filename
