"""Unit tests for the VCF decompression utility.

Verifies correct behavior for file type checks, python fallback extraction,
bgzip invocation, output validation, and CLI argument parsing.
"""

import gzip
import os
import subprocess
from unittest.mock import MagicMock, patch

import pytest

from rigi_analysis.utils.decompress_vcf import (
    decompress_vcf,
    is_gzip_file,
    main,
)


@pytest.fixture
def temp_vcf_gz(tmp_path):
    """Fixture to generate a temporary valid .vcf.gz file."""
    filepath = tmp_path / "sample.vcf.gz"
    content = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\n"
        "chr1\t100\t.\tA\tT\n"
    )
    with gzip.open(filepath, "wt", encoding="utf-8") as f:
        f.write(content)
    return filepath, content


@pytest.fixture
def temp_invalid_gz(tmp_path):
    """Fixture to generate a temporary invalid .gz file."""
    filepath = tmp_path / "invalid.vcf.gz"
    with open(filepath, "w", encoding="utf-8") as f:
        f.write("not a gzip file")
    return filepath


class TestVCFDecompression:
    """Test suite for VCF decompression logic."""

    def test_is_gzip_file_valid(self, temp_vcf_gz):
        """Test is_gzip_file with a valid gzipped file."""
        filepath, _ = temp_vcf_gz
        assert is_gzip_file(str(filepath)) is True

    def test_is_gzip_file_invalid(self, temp_invalid_gz):
        """Test is_gzip_file with an invalid gzipped file."""
        assert is_gzip_file(str(temp_invalid_gz)) is False

    def test_is_gzip_file_missing(self, tmp_path):
        """Test is_gzip_file with a non-existent file."""
        missing_path = tmp_path / "missing.vcf.gz"
        assert is_gzip_file(str(missing_path)) is False

    def test_decompress_vcf_python_fallback(self, temp_vcf_gz, tmp_path):
        """Test fallback decompression using Python's gzip library."""
        filepath, expected_content = temp_vcf_gz
        output_path = tmp_path / "decompressed.vcf"

        # Mock shutil.which to return None so bgzip is not found
        with patch("shutil.which", return_value=None):
            result_path = decompress_vcf(
                str(filepath),
                output_path=str(output_path),
                keep=True,
                force=False,
            )
            assert result_path == str(output_path)
            assert os.path.exists(output_path)
            with open(output_path, "r", encoding="utf-8") as f:
                content = f.read()
            assert content == expected_content
            # Source should still exist
            assert os.path.exists(filepath)

    def test_decompress_vcf_default_output_path(self, temp_vcf_gz):
        """Test default output path generation (removing .gz extension)."""
        filepath, expected_content = temp_vcf_gz
        expected_output = str(filepath)[:-3]

        with patch("shutil.which", return_value=None):
            result_path = decompress_vcf(str(filepath), keep=True)
            assert result_path == expected_output
            assert os.path.exists(expected_output)
            with open(expected_output, "r", encoding="utf-8") as f:
                content = f.read()
            assert content == expected_content
            # Clean up
            os.remove(expected_output)

    def test_decompress_vcf_delete_source(self, temp_vcf_gz, tmp_path):
        """Test decompress_vcf deletes source file when keep=False."""
        filepath, expected_content = temp_vcf_gz
        output_path = tmp_path / "decompressed_delete.vcf"

        with patch("shutil.which", return_value=None):
            result_path = decompress_vcf(
                str(filepath), output_path=str(output_path), keep=False
            )
            assert result_path == str(output_path)
            assert os.path.exists(output_path)
            # Source file should be deleted
            assert not os.path.exists(filepath)

    def test_decompress_vcf_missing_input(self, tmp_path):
        """Test decompress_vcf raises error for non-existent input."""
        missing_path = tmp_path / "does_not_exist.vcf.gz"
        with pytest.raises(FileNotFoundError):
            decompress_vcf(str(missing_path))

    def test_decompress_vcf_invalid_gzip(self, temp_invalid_gz):
        """Test decompress_vcf raises error for invalid gzip file."""
        with pytest.raises(ValueError, match="not a valid Gzip file"):
            decompress_vcf(str(temp_invalid_gz))

    def test_decompress_vcf_already_exists_error(self, temp_vcf_gz, tmp_path):
        """Test decompress_vcf raises error if output already exists."""
        filepath, _ = temp_vcf_gz
        output_path = tmp_path / "existing.vcf"
        # Pre-create output file
        with open(output_path, "w", encoding="utf-8") as f:
            f.write("existing")

        with pytest.raises(ValueError, match="already exists"):
            decompress_vcf(
                str(filepath), output_path=str(output_path), force=False
            )

    def test_decompress_vcf_already_exists_force(self, temp_vcf_gz, tmp_path):
        """Test decompress_vcf overwrites output file if force=True."""
        filepath, expected_content = temp_vcf_gz
        output_path = tmp_path / "existing.vcf"
        with open(output_path, "w", encoding="utf-8") as f:
            f.write("existing")

        with patch("shutil.which", return_value=None):
            result_path = decompress_vcf(
                str(filepath),
                output_path=str(output_path),
                force=True,
            )
            assert result_path == str(output_path)
            with open(output_path, "r", encoding="utf-8") as f:
                content = f.read()
            assert content == expected_content

    def test_decompress_vcf_invalid_vcf_header(self, tmp_path):
        """Test decompress_vcf fails validation if output is not VCF."""
        # Create a valid gzip file that does NOT contain valid VCF header
        filepath = tmp_path / "not_vcf.gz"
        with gzip.open(filepath, "wt", encoding="utf-8") as f:
            f.write("invalid vcf header content")

        with patch("shutil.which", return_value=None):
            with pytest.raises(
                RuntimeError, match="Decompressed file validation failed"
            ):
                decompress_vcf(str(filepath))

    def test_decompress_vcf_bgzip_success(self, temp_vcf_gz, tmp_path):
        """Test decompression using the bgzip command-line utility."""
        filepath, expected_content = temp_vcf_gz
        output_path = tmp_path / "bgzip_out.vcf"

        # Mock subprocess.run to simulate bgzip decompressing
        def mock_run(cmd, stdout, **kwargs):
            # Write expected VCF content to output stream
            stdout.write(expected_content.encode("utf-8"))
            return MagicMock(returncode=0)

        with patch("shutil.which", return_value="/usr/bin/bgzip"), \
             patch("subprocess.run", side_effect=mock_run) as mock_sub:
            result_path = decompress_vcf(
                str(filepath), output_path=str(output_path)
            )
            assert result_path == str(output_path)
            assert os.path.exists(output_path)
            with open(output_path, "r", encoding="utf-8") as f:
                content = f.read()
            assert content == expected_content
            mock_sub.assert_called_once()

    def test_decompress_vcf_bgzip_failure_fallback(self, temp_vcf_gz, tmp_path):
        """Test subprocess failure triggers fallback to Python's gzip."""
        filepath, expected_content = temp_vcf_gz
        output_path = tmp_path / "fallback_out.vcf"

        # Mock shutil.which to find bgzip, but subprocess.run raises error
        with patch("shutil.which", return_value="/usr/bin/bgzip"), \
             patch(
                 "subprocess.run",
                 side_effect=subprocess.SubprocessError("Failed"),
             ):
            result_path = decompress_vcf(
                str(filepath), output_path=str(output_path)
            )
            assert result_path == str(output_path)
            assert os.path.exists(output_path)
            with open(output_path, "r", encoding="utf-8") as f:
                content = f.read()
            assert content == expected_content


class TestDecompressCLI:
    """Test suite for CLI arguments and main function."""

    @patch("rigi_analysis.utils.decompress_vcf.decompress_vcf")
    def test_main_cli_success(self, mock_decompress):
        """Test CLI main entry point with correct arguments."""
        mock_decompress.return_value = "output.vcf"
        with patch(
            "sys.argv",
            [
                "rigi-analysis-vcf-decompress",
                "-i",
                "input.vcf.gz",
                "-o",
                "output.vcf",
            ],
        ):
            main()
            mock_decompress.assert_called_once_with(
                input_path="input.vcf.gz",
                output_path="output.vcf",
                keep=True,
                force=False,
            )

    @patch("rigi_analysis.utils.decompress_vcf.decompress_vcf")
    def test_main_cli_delete_source(self, mock_decompress):
        """Test CLI main entry point with --delete-source and -f flags."""
        mock_decompress.return_value = "output.vcf"
        with patch(
            "sys.argv",
            [
                "rigi-analysis-vcf-decompress",
                "-i",
                "input.vcf.gz",
                "--delete-source",
                "-f",
            ],
        ):
            main()
            mock_decompress.assert_called_once_with(
                input_path="input.vcf.gz",
                output_path=None,
                keep=False,
                force=True,
            )

    @patch(
        "rigi_analysis.utils.decompress_vcf.decompress_vcf",
        side_effect=ValueError("Some error"),
    )
    def test_main_cli_error(self, mock_decompress):
        """Test CLI main entry point exits with code 1 upon error."""
        with patch(
            "sys.argv", ["rigi-analysis-vcf-decompress", "-i", "input.vcf.gz"]
        ):
            with pytest.raises(SystemExit) as exc_info:
                main()
            assert exc_info.value.code == 1
