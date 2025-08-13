# ONLY include: _process_one_sequence_to_volume

import os
import shutil
import tempfile
from pathlib import Path

from shared_utils.schemas import SequenceInput
from shared_utils.file_utils import get_zip_by_idr_dir


def _process_one_sequence_to_volume(sequence_input: SequenceInput) -> str:
    """Shared helper to run Nardini for one sequence and write the zip into the volume.

    Returns absolute path to the written zip in `zipfiles/by_idr`.
    Raises on error.
    """
    from nardini.constants import ( # type: ignore
        DEFAULT_RANDOM_SEED,
        NUM_SCRAMBLED_SEQUENCES,
        TYPEALL,
    )
    from nardini.score_and_plot import calculate_zscore_and_plot # type: ignore
    from nardini.utils import set_random_seed # type: ignore

    sequence = sequence_input["sequence"]
    seq_uuid = sequence_input["seq_uuid"]

    amino_acid_groupings = TYPEALL
    num_scrambled_sequences = NUM_SCRAMBLED_SEQUENCES
    random_seed = set_random_seed(DEFAULT_RANDOM_SEED)
    output_dir = get_zip_by_idr_dir()

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        old_cwd = os.getcwd()
        os.chdir(tmp_path)
        try:
            calculate_zscore_and_plot(
                [sequence],
                amino_acid_groupings,
                num_scrambled_sequences,
                random_seed,
            )

            files = list(tmp_path.iterdir())
            zip_files = [f for f in files if f.name.endswith(".zip")]
            if len(zip_files) != 1:
                raise RuntimeError(
                    f"Expected one zip for sequence {sequence.id}, found {len(zip_files)}"
                )

            produced_zip_path = zip_files[0]
            dest_path = Path(output_dir) / f"{seq_uuid}.zip"
            shutil.move(str(produced_zip_path), str(dest_path))
        finally:
            os.chdir(old_cwd)

    if not Path(dest_path).exists():
        raise RuntimeError(f"Output zip file {dest_path} does not exist")
    return str(dest_path)
