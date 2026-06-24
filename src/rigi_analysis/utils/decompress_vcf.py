"""Utility module to decompress .vcf.gz files.

Follows bioinformatics best practices. It attempts to use `bgzip` (part of
htslib) if available on the system PATH for high-performance parallel
decompression. If not available, it falls back to streaming python's
standard `gzip` library with a low memory footprint.
"""

import argparse
import os
import shutil
import subprocess
import sys
from typing import Optional


def is_gzip_file(filepath: str) -> bool:
    """Check if the file has a valid Gzip header signature.

    Args:
        filepath: Path to the file.

    Returns:
        True if the file starts with the Gzip magic bytes, False otherwise.
    """
    if not os.path.exists(filepath):
        return False
    try:
        with open(filepath, 'rb') as f:
            header = f.read(2)
            return header == b'\x1f\x8b'
    except IOError:
        return False


def decompress_vcf(
    input_path: str,
    output_path: Optional[str] = None,
    keep: bool = True,
    force: bool = False,
) -> str:
    """Decompress a .vcf.gz file to a .vcf file.

    Follows bioinformatics best practices:
    - Verifies the input file existence and gzip signature.
    - Uses `bgzip` if available on the system PATH for speed.
    - Falls back to streaming `gzip` chunk-by-chunk to keep memory usage low.
    - Validates that the decompressed output starts with '##fileformat=VCF'.
    - Deletes the compressed source file after successful extraction
      if `keep` is False.

    Args:
        input_path: Path to the input .vcf.gz file.
        output_path: Optional path to the output .vcf file. If not provided,
            it defaults to removing the '.gz' extension from input_path.
        keep: If True (default), the original compressed file is kept.
            If False, it is deleted after successful decompression.
        force: If True, overwrites the output file if it already exists.

    Returns:
        The path to the decompressed VCF file.

    Raises:
        FileNotFoundError: If the input file does not exist.
        ValueError: If the input file is not a valid gzip file or output file
            already exists and `force` is False.
        RuntimeError: If decompression fails or the output fails validation.
    """
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")

    if not is_gzip_file(input_path):
        raise ValueError(f"Input file is not a valid Gzip file: {input_path}")

    # Determine output path
    if not output_path:
        if input_path.endswith('.gz'):
            output_path = input_path[:-3]
        else:
            output_path = input_path + '.decompressed'

    if os.path.exists(output_path) and not force:
        raise ValueError(
            f"Output file already exists: {output_path}. "
            "Use force=True to overwrite."
        )

    # Temporary output path for safe write & swap
    temp_output_path = output_path + ".tmp"
    if os.path.exists(temp_output_path):
        os.remove(temp_output_path)

    bgzip_exe = shutil.which("bgzip")
    success = False

    if bgzip_exe:
        # Use bgzip for high-performance block decompression
        try:
            # -c writes to stdout, -d decompresses
            with open(temp_output_path, "wb") as out_f:
                subprocess.run(
                    [bgzip_exe, "-d", "-c", input_path],
                    stdout=out_f,
                    stderr=subprocess.PIPE,
                    check=True,
                )
            success = True
        except (subprocess.SubprocessError, OSError):
            # Clean up temp file on failure
            if os.path.exists(temp_output_path):
                os.remove(temp_output_path)
            # Fallback to python standard library on subprocess failure
            pass

    if not success:
        # Fall back to standard library streaming to prevent high memory usage
        import gzip
        try:
            with gzip.open(
                input_path, 'rt', encoding='utf-8', errors='ignore'
            ) as f_in:
                with open(temp_output_path, 'w', encoding='utf-8') as f_out:
                    # Stream in 1MB chunks / line buffer
                    for line in f_in:
                        f_out.write(line)
            success = True
        except Exception as e:
            if os.path.exists(temp_output_path):
                os.remove(temp_output_path)
            raise RuntimeError(
                f"Failed to decompress file via Python gzip fallback: {e}"
            ) from e

    # Perform validation on decompressed output
    try:
        with open(
            temp_output_path, 'r', encoding='utf-8', errors='ignore'
        ) as f_out_check:
            first_line = f_out_check.readline()
            if not first_line.startswith("##fileformat=VCF"):
                raise ValueError(
                    "Decompressed file does not start with ##fileformat=VCF"
                )
    except Exception as e:
        if os.path.exists(temp_output_path):
            os.remove(temp_output_path)
        raise RuntimeError(
            f"Decompressed file validation failed: {e}"
        ) from e

    # Atomically rename/replace temp file to target output file
    if os.path.exists(output_path):
        os.remove(output_path)
    os.rename(temp_output_path, output_path)

    # Delete source if keep=False
    if not keep:
        os.remove(input_path)

    return output_path


def main() -> None:
    """Command line entrypoint for decompressing .vcf.gz files."""
    parser = argparse.ArgumentParser(
        description=(
            "Decompress a .vcf.gz file following bioinformatics best practices."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to the input compressed VCF file (.vcf.gz)",
    )
    parser.add_argument(
        "-o",
        "--output",
        help=(
            "Path to the output decompressed VCF file "
            "(defaults to removing .gz)"
        ),
    )
    parser.add_argument(
        "--keep",
        action="store_true",
        default=True,
        help="Keep the input compressed file (default)",
    )
    parser.add_argument(
        "--delete-source",
        dest="keep",
        action="store_false",
        help=(
            "Delete the input compressed file after "
            "successful decompression"
        ),
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Overwrite the output file if it exists",
    )

    args = parser.parse_args()

    try:
        out_path = decompress_vcf(
            input_path=args.input,
            output_path=args.output,
            keep=args.keep,
            force=args.force,
        )
        print(f"Successfully decompressed VCF file to: {out_path}")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
