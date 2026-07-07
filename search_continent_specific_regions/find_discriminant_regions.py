#!/usr/bin/env python3
"""
VCF Window Analyzer - Identifies 500nt genomic windows that distinguish geographic groups 
"""

import sys
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Set

class VCFWindowAnalyzer:
    def __init__(self, vcf_file: str, groups_file: str, window_size: int = 500, threshold: float = 0.75):
        self.vcf_file = vcf_file
        self.groups_file = groups_file
        self.window_size = window_size
        self.threshold = threshold
        
        # Data structures
        self.sample_to_group = {}  # sample_id -> group_name
        self.group_to_samples = defaultdict(list)  # group_name -> [sample_ids]
        self.samples = []  # ordered list of samples from VCF
        self.results = []  # final results
        
    def load_groups(self):
        """Load sample-to-group mapping from file"""
        print(f"Loading groups from {self.groups_file}...")
        with open(self.groups_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= 2:
                    sample, group = parts[0], parts[1]
                    self.sample_to_group[sample] = group
                    self.group_to_samples[group].append(sample)
        
        print(f"Loaded {len(self.sample_to_group)} samples in {len(self.group_to_samples)} groups:")
        for group, samples in self.group_to_samples.items():
            print(f"  {group}: {len(samples)} samples")
    
    def parse_genotype(self, gt_field: str) -> str:
        """Extract genotype from GT:VF format (e.g., '0/0:100' -> '0')"""
        gt = gt_field.split(':')[0]
        if gt == './.':
            return '.'
        # Take only the 3rd character (allele after separator)
        if len(gt) >= 3:
            return gt[2]
        return '.'
    
    def analyze_vcf(self):
        """Main analysis: parse VCF and analyze windows"""
        print(f"\nAnalyzing VCF file {self.vcf_file}...")
        
        # Storage for sliding window
        window_variants = []  # List of (chrom, pos, {sample: genotype})
        current_chrom = None
        variant_count = 0
        windows_analyzed = 0
        
        with open(self.vcf_file, 'r') as f:
            for line in f:
                line = line.strip()
                
                # Skip empty lines
                if not line:
                    continue
                
                # Parse header to get sample names
                if line.startswith('#CHROM'):
                    parts = line.split('\t')
                    self.samples = parts[9:]  # samples start at column 9
                    print(f"Found {len(self.samples)} samples in VCF")
                    
                    # Filter samples to only those in groups file
                    self.samples = [s for s in self.samples if s in self.sample_to_group]
                    print(f"Keeping {len(self.samples)} samples that are in groups file")
                    continue
                
                # Skip other header lines
                if line.startswith('#'):
                    continue
                
                # Parse variant line
                parts = line.split('\t')
                chrom = parts[0]
                pos = int(parts[1])
                genotypes = parts[9:]
                
                variant_count += 1
                if variant_count % 1000 == 0:
                    print(f"  Processed {variant_count} variants, analyzed {windows_analyzed} windows...", end='\r')
                
                # If chromosome changed, clear window
                if current_chrom != chrom:
                    current_chrom = chrom
                    window_variants = []
                    print(f"\n  Starting chromosome {chrom}")
                
                # Extract genotypes for samples in our groups
                variant_genotypes = {}
                for i, sample in enumerate(self.samples):
                    if i < len(genotypes):
                        gt = self.parse_genotype(genotypes[i])
                        variant_genotypes[sample] = gt
                
                # Add current variant to window
                window_variants.append((chrom, pos, variant_genotypes))
                
                # Remove variants outside the window (older than pos - window_size + 1)
                window_start = pos - self.window_size + 1
                window_variants = [(c, p, g) for c, p, g in window_variants if p >= window_start]
                
                # Analyze window ending at current position
                if len(window_variants) > 0:
                    self.analyze_window(window_variants)
                    windows_analyzed += 1
        
        print(f"\n\nTotal variants processed: {variant_count}")
        print(f"Total windows analyzed: {windows_analyzed}")
        print(f"Interesting regions found: {len(self.results)}")
    
    def analyze_window(self, window_variants: List[Tuple[str, int, Dict[str, str]]]):
        """Analyze a single window of variants"""
        if not window_variants:
            return
        
        chrom = window_variants[0][0]
        start_pos = window_variants[0][1]
        end_pos = window_variants[-1][1]
        
        # Build profiles for each sample
        profiles = {}
        for sample in self.samples:
            profile_parts = []
            has_missing = False
            
            for _, _, genotypes in window_variants:
                gt = genotypes.get(sample, '.')
                if gt == '.':
                    has_missing = True
                    break
                profile_parts.append(gt)
            
            # Only keep profile if no missing data
            if not has_missing:
                profiles[sample] = ''.join(profile_parts)
        
        # Analyze profiles by group
        group_profiles = defaultdict(list)  # group -> [profiles]
        for sample, profile in profiles.items():
            group = self.sample_to_group[sample]
            group_profiles[group].append(profile)
        
        # Find consensus profile for each group (≥75% threshold)
        group_consensus = {}  # group -> (consensus_profile, percentage)
        for group, group_prof_list in group_profiles.items():
            if not group_prof_list:
                continue
            
            # Count profile frequencies
            profile_counts = Counter(group_prof_list)
            most_common_profile, count = profile_counts.most_common(1)[0]
            percentage = (count / len(group_prof_list)) * 100
            
            # Check if meets threshold
            if percentage >= self.threshold * 100:
                group_consensus[group] = (most_common_profile, percentage)
        
        # Check if ALL groups have consensus AND they are all different
        # We need ALL groups to have consensus, not just some of them
        all_groups = set(self.group_to_samples.keys())
        groups_with_consensus = set(group_consensus.keys())
        
        if groups_with_consensus == all_groups and len(group_consensus) >= 2:
            # All groups have consensus ≥75%
            consensus_profiles = [prof for prof, _ in group_consensus.values()]
            
            # Check if all consensus profiles are unique
            if len(consensus_profiles) == len(set(consensus_profiles)):
                # This is an interesting region!
                self.results.append({
                    'chrom': chrom,
                    'start': start_pos,
                    'end': end_pos,
                    'group_percentages': {g: pct for g, (_, pct) in group_consensus.items()},
                    'group_sample_counts': {g: len(group_prof_list) for g, group_prof_list in group_profiles.items()},
                    'is_interesting': 'Yes'
                })
    
    def write_results(self, output_file: str):
        """Write results to TSV file"""
        print(f"\nWriting results to {output_file}...")
        
        # Get all groups for header
        all_groups = sorted(self.group_to_samples.keys())
        
        with open(output_file, 'w') as f:
            # Header
            header_parts = ['Chromosome', 'Start', 'End']
            for group in all_groups:
                header_parts.append(f'{group}_percentage')
            for group in all_groups:
                header_parts.append(f'{group}_n_samples')
            header_parts.append('Interesting')
            f.write('\t'.join(header_parts) + '\n')
            
            # Data rows
            for result in self.results:
                row_parts = [
                    result['chrom'],
                    str(result['start']),
                    str(result['end'])
                ]
                
                # Add percentages
                for group in all_groups:
                    pct = result['group_percentages'].get(group, None)
                    if pct is not None:
                        row_parts.append(f'{pct:.2f}')
                    else:
                        row_parts.append('N/A')  # Should not happen with fixed logic
                
                # Add sample counts
                for group in all_groups:
                    count = result['group_sample_counts'].get(group, 0)
                    row_parts.append(str(count))
                
                row_parts.append(result['is_interesting'])
                f.write('\t'.join(row_parts) + '\n')
        
        print(f"Results written successfully!")

def main():
    if len(sys.argv) != 4:
        print("Usage: python vcf_window_analyzer.py <vcf_file> <groups_file> <output_file>")
        print("Example: python vcf_window_analyzer.py data.vcf groups.txt results.tsv")
        sys.exit(1)
    
    vcf_file = sys.argv[1]
    groups_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Create analyzer
    analyzer = VCFWindowAnalyzer(vcf_file, groups_file)
    
    # Load groups
    analyzer.load_groups()
    
    # Analyze VCF
    analyzer.analyze_vcf()
    
    # Write results
    analyzer.write_results(output_file)
    
    print("\nAnalysis complete!")

if __name__ == '__main__':
    main()
