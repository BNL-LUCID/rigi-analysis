    #!/usr/bin/env python
    """
    Unified mutation pattern analysis: compute both binary (0/1) and categorical (C/T/B/0) patterns.
    Replaces mutation_analysis.py + pattern computation from sankey_visualization.py
    """

    import pandas as pd
    import numpy as np
    import os
    import argparse
    from glob import glob
    from collections import defaultdict
    import re
    from concurrent.futures import ProcessPoolExecutor
    import time


    class MutationPatternAnalyzer:
        def __init__(self, input_dir, output_dir):
            self.input_dir = input_dir
            self.output_dir = output_dir
            os.makedirs(output_dir, exist_ok=True)
            
            # Control pattern for dose identification
            self.control_pattern = r'(?:d0|D0|[Cc]ontrol)'
        
        def load_annotated_mutations(self, chromosome=None, mutation_type='SNV'):
            """
            Load annotated mutation files for specified chromosome and type.
            
            Args:
                chromosome: Specific chromosome (e.g., 'chr1') or None for all
                mutation_type: SNV, DBS, MNS, or ID
            
            Returns:
                DataFrame with all mutations
            """
            if chromosome:
                pattern = os.path.join(self.input_dir, f"{chromosome}_{mutation_type}_annotated.pkl")
            else:
                pattern = os.path.join(self.input_dir, f"*_{mutation_type}_annotated.pkl")
            
            files = glob(pattern)
            
            if not files:
                print(f"No files found matching: {pattern}")
                return None
            
            print(f"Loading {len(files)} file(s) for {mutation_type}...")
            
            all_data = []
            for file_path in files:
                try:
                    df = pd.read_pickle(file_path)
                    # Only derive chromosome from filename if the DataFrame doesn't
                    # already carry real per-row chromosomes. Combined-file inputs
                    # (e.g. all_DBS_annotated.pkl) keep their existing Chromosome column.
                    if 'Chromosome' not in df.columns or df['Chromosome'].nunique() <= 1:
                        chrom = os.path.basename(file_path).split('_')[0]
                        df['Chromosome'] = chrom
                        label = chrom
                    else:
                        label = f"{file_path} ({df['Chromosome'].nunique()} chroms)"
                    all_data.append(df)
                    print(f"  Loaded {len(df)} mutations from {label}")
                except Exception as e:
                    print(f"  Error loading {file_path}: {e}")
            
            if not all_data:
                return None
            
            combined = pd.concat(all_data, ignore_index=True)
            print(f"Total mutations loaded: {len(combined)}")
            
            return combined
        
        def create_permanent_mutation_id(self, df):
            """
            Create consistent mutation IDs across timepoints.
            Format: chr_pos_ref_alt
            """
            if 'PermanentMutationID' in df.columns:
                return df
            
            if all(col in df.columns for col in ['Chromosome', 'Start', 'Ref', 'Alt']):
                df['PermanentMutationID'] = (
                    df['Chromosome'].astype(str) + '_' +
                    df['Start'].astype(str) + '_' +
                    df['Ref'] + '_' +
                    df['Alt']
                )
            elif 'MutationID' in df.columns:
                # Try to extract from existing MutationID
                df['PermanentMutationID'] = df['MutationID'].str.replace(
                    r'(.*_.*_.*_.*)_.*_.*', r'\1', regex=True
                )
            else:
                raise ValueError("Cannot create PermanentMutationID - missing required columns")
            
            print(f"Created PermanentMutationID for {len(df)} mutations")
            return df
        
        def compute_binary_patterns(self, df, doses=None):
            """
            Compute binary (0/1) presence/absence patterns across timepoints.
            
            Args:
                df: DataFrame with mutation data
                doses: List of doses to analyze, or None for all
            
            Returns:
                DataFrame with binary pattern columns
            """
            print("\nComputing binary (0/1) patterns...")
            
            # Get unique doses if not specified
            if doses is None:
                doses = sorted(df['Dose'].unique())
            
            # Get timepoints
            timepoints = sorted([tp for tp in df['Timepoint'].unique() if tp.startswith('W')])
            
            # Group by PermanentMutationID
            grouped = df.groupby('PermanentMutationID')
            
            results = []
            for mut_id, group in grouped:
                record = {'MutationID': mut_id}
                
                # Add genomic coordinates
                first_row = group.iloc[0]
                record['Chromosome'] = first_row['Chromosome']
                if 'Start' in group.columns:
                    record['Start'] = first_row['Start']
                if 'Ref' in group.columns:
                    record['Ref'] = first_row['Ref']
                if 'Alt' in group.columns:
                    record['Alt'] = first_row['Alt']
                
                # Binary patterns for each dose-timepoint combination
                for dose in doses:
                    for tp in timepoints:
                        mask = (group['Dose'] == dose) & (group['Timepoint'] == tp)
                        record[f'{dose}_{tp}'] = 1 if mask.any() else 0
                
                results.append(record)
            
            binary_df = pd.DataFrame(results)
            print(f"Computed binary patterns for {len(binary_df)} unique mutations")
            
            return binary_df
        
        def compute_categorical_patterns(self, df):
            """
            Compute categorical (C/T/B/0) patterns distinguishing control vs treatment.
            
            C = Control only
            T = Treatment only  
            B = Both control and treatment
            0 = Absent
            
            Args:
                df: DataFrame with mutation data
            
            Returns:
                DataFrame with categorical pattern columns
            """
            print("\nComputing categorical (C/T/B/0) patterns...")
            
            # Create binary flags for control vs treatment at each timepoint
            df = self._add_control_treatment_flags(df)
            
            # Group by PermanentMutationID and aggregate
            grouped = df.groupby('PermanentMutationID').agg({
                'Chromosome': 'first',
                'Start': 'first' if 'Start' in df.columns else lambda x: None,
                'Ref': 'first' if 'Ref' in df.columns else lambda x: None,
                'Alt': 'first' if 'Alt' in df.columns else lambda x: None,
                'control_w0': 'max',
                'control_w1': 'max',
                'control_w2': 'max',
                'control_w3': 'max',
                'treated_w1': 'max',
                'treated_w2': 'max',
                'treated_w3': 'max'
            }).reset_index()
            
            grouped = grouped.rename(columns={'PermanentMutationID': 'MutationID'})
            
            # Compute categorical states for each week
            grouped['W0'] = grouped.apply(
                lambda row: 'W0_Present' if row['control_w0'] else 'W0_Absent', 
                axis=1
            )
            
            for week in [1, 2, 3]:
                grouped[f'W{week}'] = grouped.apply(
                    lambda row: self._determine_week_state(
                        row[f'control_w{week}'], 
                        row[f'treated_w{week}'],
                        week
                    ),
                    axis=1
                )
            
            # Create pattern string (e.g., "0T00", "CBBB")
            grouped['Pattern'] = grouped.apply(
                lambda row: self._create_pattern_string(row), 
                axis=1
            )
            
            # Categorize patterns
            grouped['Category'] = grouped.apply(
                lambda row: self._categorize_pattern(row), 
                axis=1
            )
            
            # Drop intermediate columns
            cols_to_keep = [
                'MutationID', 'Chromosome', 'Start', 'Ref', 'Alt',
                'W0', 'W1', 'W2', 'W3', 'Pattern', 'Category'
            ]
            categorical_df = grouped[[col for col in cols_to_keep if col in grouped.columns]]
            
            print(f"Computed categorical patterns for {len(categorical_df)} unique mutations")
            
            return categorical_df
        
        def _add_control_treatment_flags(self, df):
            """Add binary flags for control vs treatment presence."""
            control_mask = df['Dose'].str.contains(self.control_pattern, na=False, regex=True)
            
            timepoint_masks = {
                'w0': df['Timepoint'] == 'W0',
                'w1': df['Timepoint'] == 'W1',
                'w2': df['Timepoint'] == 'W2',
                'w3': df['Timepoint'] == 'W3'
            }
            
            for tp, mask in timepoint_masks.items():
                df[f'control_{tp}'] = mask & control_mask
                if tp != 'w0':  # W0 is always control (baseline)
                    df[f'treated_{tp}'] = mask & ~control_mask
            
            # Fill NAs with False
            flag_cols = [col for col in df.columns if col.startswith(('control_', 'treated_'))]
            for col in flag_cols:
                df[col] = df[col].fillna(False)
            
            return df
        
        def _determine_week_state(self, control, treated, week):
            """Determine categorical state for a given week."""
            if control and treated:
                return f'W{week}_Both'
            elif control:
                return f'W{week}_Control'
            elif treated:
                return f'W{week}_Treatment'
            else:
                return f'W{week}_Lost'
        
        def _create_pattern_string(self, row):
            """Create compact pattern string (e.g., 'CBTT', '0T00')."""
            pattern = []
            
            # W0
            pattern.append('C' if 'Present' in row['W0'] else '0')
            
            # W1-W3
            for week in [1, 2, 3]:
                state = row[f'W{week}']
                if 'Both' in state:
                    pattern.append('B')
                elif 'Control' in state:
                    pattern.append('C')
                elif 'Treatment' in state:
                    pattern.append('T')
                else:
                    pattern.append('0')
            
            return ''.join(pattern)
        
        def _categorize_pattern(self, row):
            """Assign category based on pattern."""
            pattern = row['Pattern']
            w0, w1, w2, w3 = pattern
            
            # Control-only patterns
            if w0 == 'C' and w1 == 'C' and w2 == 'C' and w3 == 'C':
                return 'Control_All_Timepoints'
            
            # Treatment-only patterns
            if w0 == '0' and w1 == 'T' and w2 == 'T' and w3 == 'T':
                return 'Treatment_Only_All'
            
            if w0 == '0' and w1 == 'T' and w2 == '0' and w3 == '0':
                return 'Treatment_Only_W1'
            
            if w0 == '0' and w1 == '0' and w2 == 'T' and w3 == '0':
                return 'Treatment_Only_W2'
            
            if w0 == '0' and w1 == '0' and w2 == '0' and w3 == 'T':
                return 'Treatment_Only_W3'
            
            # Both patterns
            if w0 == 'C' and w1 == 'B' and w2 == 'B' and w3 == 'B':
                return 'Present_All_Timepoints_Both'
            
            # W0 only
            if w0 == 'C' and w1 == '0' and w2 == '0' and w3 == '0':
                return 'Present_W0_Only'
            
            # Generic category based on pattern
            return f'Pattern_{pattern}'
        
        def analyze_by_dose(self, df, mutation_type='SNV'):
            """
            Analyze mutations separately for each dose.
            
            Args:
                df: Combined mutation DataFrame
                mutation_type: Type identifier for output files
            
            Returns:
                dict: {dose: (binary_df, categorical_df)}
            """
            doses = sorted([d for d in df['Dose'].unique() 
                        if not re.match(self.control_pattern, d)])
            
            results = {}
            
            for dose in doses:
                print(f"\n{'='*60}")
                print(f"Processing Dose: {dose}")
                print(f"{'='*60}")
                
                # Filter for this dose + control
                dose_data = df[
                    df['Dose'].str.contains(f'{dose}|{self.control_pattern}', 
                                        na=False, regex=True)
                ].copy()
                
                print(f"Mutations for dose {dose} (including control): {len(dose_data)}")
                
                # Compute both pattern types
                binary_df = self.compute_binary_patterns(dose_data, doses=[dose])
                categorical_df = self.compute_categorical_patterns(dose_data)
                
                # Save outputs
                self._save_dose_outputs(binary_df, categorical_df, dose, mutation_type)
                
                results[dose] = (binary_df, categorical_df)
            
            return results
        
        def _save_dose_outputs(self, binary_df, categorical_df, dose, mutation_type):
            """Save binary and categorical patterns for a specific dose."""
            dose_dir = os.path.join(self.output_dir, f'dose_{dose}')
            os.makedirs(dose_dir, exist_ok=True)
            
            # Save binary patterns
            binary_csv = os.path.join(dose_dir, f'all_mutations_dose_{dose}_{mutation_type}.csv')
            binary_pkl = os.path.join(dose_dir, f'all_mutations_dose_{dose}_{mutation_type}.pkl')
            binary_df.to_csv(binary_csv, index=False)
            binary_df.to_pickle(binary_pkl)
            print(f"\nSaved binary patterns:")
            print(f"  CSV: {binary_csv}")
            print(f"  PKL: {binary_pkl}")
            
            # Save categorical patterns
            categorical_csv = os.path.join(dose_dir, f'mutation_annotations_dose_{dose}_{mutation_type}.csv')
            categorical_pkl = os.path.join(dose_dir, f'mutation_annotations_dose_{dose}_{mutation_type}.pkl')
            categorical_df.to_csv(categorical_csv, index=False)
            categorical_df.to_pickle(categorical_pkl)
            print(f"\nSaved categorical patterns:")
            print(f"  CSV: {categorical_csv}")
            print(f"  PKL: {categorical_pkl}")
        
        def analyze_by_chromosome(self, mutation_type='SNV'):
            """
            Analyze mutations chromosome by chromosome.
            
            Args:
                mutation_type: SNV, DBS, MNS, or ID
            """
            # Get all chromosome files
            pattern = os.path.join(self.input_dir, f"*_{mutation_type}_annotated.pkl")
            files = glob(pattern)
            chromosomes = sorted(list(set([os.path.basename(f).split('_')[0] for f in files])))
            
            print(f"Found {len(chromosomes)} chromosomes to process")
            
            for chrom in chromosomes:
                print(f"\n{'='*60}")
                print(f"Processing Chromosome: {chrom}")
                print(f"{'='*60}")
                
                # Load data for this chromosome
                df = self.load_annotated_mutations(chromosome=chrom, mutation_type=mutation_type)
                
                if df is None or len(df) == 0:
                    print(f"No data for {chrom}, skipping...")
                    continue
                
                # Create permanent IDs
                df = self.create_permanent_mutation_id(df)
                
                # Analyze by dose
                self.analyze_by_dose(df, mutation_type=mutation_type)
        
        def print_pattern_summary(self, categorical_df):
            """Print summary statistics for categorical patterns."""
            print("\n" + "="*60)
            print("PATTERN SUMMARY")
            print("="*60)
            
            # Pattern distribution
            pattern_counts = categorical_df['Pattern'].value_counts()
            print("\nTop 20 patterns:")
            for pattern, count in pattern_counts.head(20).items():
                pct = (count / len(categorical_df)) * 100
                print(f"  {pattern}: {count:,} ({pct:.1f}%)")
            
            # Category distribution
            category_counts = categorical_df['Category'].value_counts()
            print("\nCategory distribution:")
            for category, count in category_counts.items():
                pct = (count / len(categorical_df)) * 100
                print(f"  {category}: {count:,} ({pct:.1f}%)")
            
            # State distribution by week
            for week in ['W0', 'W1', 'W2', 'W3']:
                state_counts = categorical_df[week].value_counts()
                print(f"\n{week} states:")
                for state, count in state_counts.items():
                    pct = (count / len(categorical_df)) * 100
                    print(f"  {state}: {count:,} ({pct:.1f}%)")


    def main():
        parser = argparse.ArgumentParser(
            description='Unified mutation pattern analysis (binary and categorical)'
        )
        parser.add_argument(
            '--input-dir', '-i', required=True,
            help='Directory containing annotated mutation PKL files'
        )
        parser.add_argument(
            '--output-dir', '-o', default='pattern_analysis',
            help='Directory to save output files'
        )
        parser.add_argument(
            '--mutation-type', '-m', default='SNV',
            choices=['SNV', 'DBS', 'MNS', 'ID'],
            help='Mutation type to analyze'
        )
        parser.add_argument(
            '--by-chromosome', '-c', action='store_true',
            help='Process each chromosome separately'
        )
        
        args = parser.parse_args()
        
        start_time = time.time()
        
        # Create analyzer
        analyzer = MutationPatternAnalyzer(
            input_dir=args.input_dir,
            output_dir=args.output_dir
        )
        
        if args.by_chromosome:
            # Process chromosome by chromosome
            analyzer.analyze_by_chromosome(mutation_type=args.mutation_type)
        else:
            # Process all chromosomes together
            print(f"Loading all {args.mutation_type} mutations...")
            df = analyzer.load_annotated_mutations(mutation_type=args.mutation_type)
            
            if df is None or len(df) == 0:
                print("No data loaded. Exiting.")
                return
            
            # Create permanent IDs
            df = analyzer.create_permanent_mutation_id(df)
            
            # Analyze by dose
            results = analyzer.analyze_by_dose(df, mutation_type=args.mutation_type)
            
            # Print summary for first dose
            if results:
                first_dose = list(results.keys())[0]
                _, categorical_df = results[first_dose]
                analyzer.print_pattern_summary(categorical_df)
        
        elapsed = time.time() - start_time
        print(f"\n{'='*60}")
        print(f"Analysis complete in {elapsed:.2f} seconds")
        print(f"{'='*60}")


    if __name__ == "__main__":
        main()
