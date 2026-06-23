import argparse
import logging
import os
import shutil
import tarfile
from typing import Optional

from SigProfilerExtractor import sigpro as sig
from SigProfilerMatrixGenerator import install as genInstall


class SigProfilerAnalysis:
    """Run SigProfilerExtractor analysis on VCF files with genome installation."""

    def __init__(self, output_dir: str, reference_genome: str = "GRCh38",
                 log_file: Optional[str] = None):
        """Initialize SigProfilerAnalysis.

        Args:
            output_dir: Output directory for results
            reference_genome: Reference genome version (default: GRCh38)
            log_file: Optional log file path
        """
        self.output_dir = output_dir
        self.reference_genome = reference_genome
        self.results = None

        # Setup logging
        os.makedirs(output_dir, exist_ok=True)
        log_path = log_file or os.path.join(output_dir, 'sigprofiler_analysis.log')

        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_path),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)

    def install_reference_genome(self, offline_files_path: Optional[str] = None):
        """Install reference genome if needed.

        Args:
            offline_files_path: Path to a previously-downloaded
                {reference_genome}.tar.gz. If provided, the file is copied into
                SigProfilerMatrixGenerator's install location and extracted,
                bypassing the FTP download. Useful on VMs without outbound FTP.
        """
        try:
            if offline_files_path:
                self.logger.info(
                    f"Installing reference genome from local tarball: {offline_files_path}"
                )
                self._install_from_local_tarball(offline_files_path)
            else:
                self.logger.info(f"Installing reference genome (online): {self.reference_genome}")
                genInstall.install(self.reference_genome)
                self.logger.info(f"Successfully installed {self.reference_genome}")
        except Exception as e:
            self.logger.error(f"Error with reference genome: {str(e)}")
            raise

    def _install_from_local_tarball(self, tarball_path: str):
        """Place a local {reference_genome}.tar.gz into SigProfilerMatrixGenerator's
        expected install directory and extract it, skipping the FTP download.
        """
        if not os.path.isfile(tarball_path):
            raise FileNotFoundError(f"Offline genome tarball not found: {tarball_path}")

        # SigProfilerMatrixGenerator looks for genomes under
        #   <install_dir>/references/chromosomes/tsb/<reference_genome>/
        import SigProfilerMatrixGenerator
        sigprofiler_dir = os.path.dirname(SigProfilerMatrixGenerator.__file__)
        tsb_dir = os.path.join(sigprofiler_dir, "references", "chromosomes", "tsb")
        extracted_dir = os.path.join(tsb_dir, self.reference_genome)

        # Skip if already extracted and non-empty
        if os.path.isdir(extracted_dir) and os.listdir(extracted_dir):
            self.logger.info(f"Reference genome already extracted at: {extracted_dir}")
            return

        os.makedirs(tsb_dir, exist_ok=True)
        target_tarball = os.path.join(tsb_dir, f"{self.reference_genome}.tar.gz")

        # Copy tarball into place unless it's already there
        src_real = os.path.realpath(tarball_path)
        dst_real = os.path.realpath(target_tarball)
        if src_real != dst_real:
            self.logger.info(f"Copying tarball into SigProfiler tree: {target_tarball}")
            shutil.copy2(tarball_path, target_tarball)
        else:
            self.logger.info(f"Tarball already at: {target_tarball}")

        self.logger.info(f"Extracting → {tsb_dir}")
        with tarfile.open(target_tarball, "r:gz") as tar:
            tar.extractall(path=tsb_dir)

        if not os.path.isdir(extracted_dir):
            raise RuntimeError(
                f"Extraction completed but expected directory not found: {extracted_dir}. "
                f"Check that the tarball contains a top-level '{self.reference_genome}/' directory."
            )

        # Some tarballs are packed with full path replay (e.g.
        # GRCh38/chromosomes/tsb/GRCh38/1.txt) so the chromosome .txt files
        # land at an arbitrary depth instead of directly under extracted_dir.
        # Detect and flatten so SigProfiler can find them.
        self._flatten_chromosome_layout(extracted_dir)

        self.logger.info(f"Reference genome installed offline at: {extracted_dir}")

    def _flatten_chromosome_layout(self, extracted_dir: str):
        """Ensure chromosome files ('1.txt', '2.txt', ..., 'MT.txt') sit directly
        under extracted_dir. If they are nested, walk the tree, find the
        deepest directory containing them, and lift its contents up.
        """
        # Quick check: if 1.txt is already at the top level we are done.
        if os.path.isfile(os.path.join(extracted_dir, "1.txt")):
            return

        self.logger.info("Chromosome files not at expected depth — searching tree...")

        # Find the directory that actually holds the chromosome files
        chrom_dir = None
        for root, _dirs, files in os.walk(extracted_dir):
            if "1.txt" in files and "MT.txt" in files:
                chrom_dir = root
                break

        if chrom_dir is None:
            raise RuntimeError(
                f"Could not locate chromosome files (1.txt, MT.txt, ...) "
                f"anywhere under {extracted_dir}. The tarball does not appear "
                f"to be a SigProfilerMatrixGenerator-formatted reference."
            )

        if os.path.realpath(chrom_dir) == os.path.realpath(extracted_dir):
            return  # already flat

        self.logger.info(f"Lifting chromosome files from {chrom_dir} → {extracted_dir}")
        for entry in os.listdir(chrom_dir):
            src = os.path.join(chrom_dir, entry)
            dst = os.path.join(extracted_dir, entry)
            if os.path.exists(dst):
                # very unlikely but be safe: remove pre-existing duplicate
                if os.path.isdir(dst):
                    shutil.rmtree(dst)
                else:
                    os.remove(dst)
            shutil.move(src, dst)

        # Clean up the now-empty wrapper directories. Walk upward from
        # chrom_dir's parent until we hit extracted_dir, removing empties.
        wrapper = os.path.dirname(chrom_dir)
        while os.path.realpath(wrapper) != os.path.realpath(extracted_dir):
            try:
                os.rmdir(wrapper)
            except OSError:
                break  # directory not empty (other artifacts) — leave alone
            wrapper = os.path.dirname(wrapper)

        self.logger.info("Chromosome layout flattened.")

    def run_analysis(self, input_vcf: str, project_name: str = "SigProfilerProject",
                    min_signatures: int = 1, max_signatures: int = 10,
                    nmf_replicates: int = 100, cpu_count: int = -1,
                    gpu: bool = False) -> dict:
        """Run SigProfilerExtractor analysis on VCF files.

        Args:
            input_vcf: Path to input VCF files or directory
            project_name: Project name for analysis
            min_signatures: Minimum number of signatures (default: 1)
            max_signatures: Maximum number of signatures (default: 10)
            nmf_replicates: Number of NMF replicates (default: 100)
            cpu_count: Number of CPUs to use (-1 for all, default: -1)
            gpu: Run NMF on CUDA-enabled GPU via PyTorch (default: False)

        Returns:
            Analysis results dictionary
        """
        if not os.path.exists(input_vcf):
            raise FileNotFoundError(f"Input VCF path not found: {input_vcf}")

        self.logger.info("Starting SigProfilerExtractor analysis")
        self.logger.info(f"Project: {project_name}")
        self.logger.info(f"Input: {input_vcf}")
        self.logger.info(f"Output: {self.output_dir}")
        self.logger.info(f"Reference genome: {self.reference_genome}")
        self.logger.info(f"GPU: {'enabled' if gpu else 'disabled'}")

        if gpu:
            try:
                import torch
                if not torch.cuda.is_available():
                    self.logger.warning(
                        "GPU requested but torch.cuda.is_available() is False. "
                        "Falling back to CPU. Check CUDA driver and PyTorch build."
                    )
                    gpu = False
                else:
                    self.logger.info(
                        f"CUDA OK: {torch.cuda.device_count()} device(s), "
                        f"using device 0 ({torch.cuda.get_device_name(0)})"
                    )
            except ImportError:
                self.logger.warning(
                    "GPU requested but PyTorch not installed. Falling back to CPU."
                )
                gpu = False

        try:
            results = sig.sigProfilerExtractor(
                input_type="vcf",
                output=self.output_dir,
                input_data=input_vcf,
                reference_genome=self.reference_genome,
                minimum_signatures=min_signatures,
                maximum_signatures=max_signatures,
                nmf_replicates=nmf_replicates,
                cpu=cpu_count,
                gpu=gpu
            )

            self.results = results
            self.logger.info("SigProfilerExtractor analysis completed successfully")
            return results

        except Exception as e:
            self.logger.error(f"Error running SigProfilerExtractor: {str(e)}")
            raise


def main():
    """Command-line interface for SigProfilerExtractor analysis."""
    parser = argparse.ArgumentParser(
        description="Run SigProfilerExtractor analysis on VCF data.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python sigprofiler.py -i ./filtered_vcfs -o ./results
  python sigprofiler.py -i ./filtered_vcfs -o ./results -r GRCh37 -M 10
        """
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input directory containing VCF files"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output directory to store results"
    )
    parser.add_argument(
        "-r", "--reference",
        default="GRCh38",
        choices=["GRCh37", "GRCh38"],
        help="Reference genome version (default: GRCh38)"
    )
    parser.add_argument(
        "-p", "--project",
        default="SigProfilerProject",
        help="Project name for analysis"
    )
    parser.add_argument(
        "-m", "--min-signatures",
        type=int,
        default=1,
        help="Minimum number of signatures (default: 1)"
    )
    parser.add_argument(
        "-M", "--max-signatures",
        type=int,
        default=10,
        help="Maximum number of signatures (default: 10)"
    )
    parser.add_argument(
        "-n", "--nmf-replicates",
        type=int,
        default=100,
        help="Number of NMF replicates (default: 100)"
    )
    parser.add_argument(
        "-c", "--cpu",
        type=int,
        default=-1,
        help="Number of CPUs to use (-1 for all, default: -1)"
    )
    parser.add_argument(
        "--gpu",
        action="store_true",
        help="Run NMF on CUDA-enabled GPU via PyTorch (requires CUDA-matched torch build)"
    )
    parser.add_argument(
        "--offline-genome",
        help="Optional path to offline genome reference files"
    )

    args = parser.parse_args()

    try:
        analysis = SigProfilerAnalysis(
            output_dir=args.output,
            reference_genome=args.reference
        )

        analysis.install_reference_genome(offline_files_path=args.offline_genome)

        analysis.run_analysis(
            input_vcf=args.input,
            project_name=args.project,
            min_signatures=args.min_signatures,
            max_signatures=args.max_signatures,
            nmf_replicates=args.nmf_replicates,
            cpu_count=args.cpu,
            gpu=args.gpu
        )

        print("\n✓ SigProfilerExtractor analysis completed successfully!")
        print(f"Results saved to: {args.output}")

    except Exception as e:
        print(f"\n✗ Analysis failed: {str(e)}")
        exit(1)


if __name__ == "__main__":
    main()
