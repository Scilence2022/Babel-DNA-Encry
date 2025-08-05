#!/usr/bin/env python3
"""
Benchmark script for comparing different similarity methods.
Helps users choose the optimal similarity method for their Z-DNA data.
"""

import sys
import time
import random
from alignment_utils import (
    AlignmentConfig,
    benchmark_alignment_methods,
    read_fasta,
    calculate_kmer_similarity,
    calculate_sliding_window_similarity,
    calculate_lcs_similarity,
    calculate_local_alignment_similarity_fast,
    calculate_pairwise_alignment_similarity,
    calculate_local_alignment_similarity
)

def benchmark_on_real_data(data_seq_file, index_seq_file, num_samples=10, iterations_per_test=100, num_repeats=3):
    """
    Benchmark similarity methods on real data files.
    """
    print("=== Similarity Method Benchmark ===")
    print(f"Data sequence file: {data_seq_file}")
    print(f"Index sequence file: {index_seq_file}")
    print(f"Number of sample comparisons: {num_samples}")
    print(f"Iterations per test: {iterations_per_test}")
    print(f"Number of repeats: {num_repeats}")
    print()
    
    # Read sequences
    try:
        data_seq = read_fasta(data_seq_file)[0]
        index_seqs = read_fasta(index_seq_file)
    except FileNotFoundError as e:
        print(f"Error: Could not find file {e}")
        return
    except IndexError:
        print(f"Error: No sequences found in {data_seq_file}")
        return
    
    print(f"Data sequence length: {len(data_seq)}")
    print(f"Number of index sequences: {len(index_seqs)}")
    print(f"Average index sequence length: {sum(len(s) for s in index_seqs) / len(index_seqs):.1f}")
    print()
    
    # Configuration
    config = AlignmentConfig(match_score=2, mismatch_score=-1, gap_open=-2, gap_extend=-0.5)
    
    # Sample random pairs for benchmarking
    results = {}
    methods = [
        ('kmer', lambda s1, s2: calculate_kmer_similarity(s1, s2, k=13)),
        ('sliding', lambda s1, s2: calculate_sliding_window_similarity(s1, s2, window_size=50, step=10)),
        ('lcs', lambda s1, s2: calculate_lcs_similarity(s1, s2)),
        ('local_fast', lambda s1, s2: calculate_local_alignment_similarity_fast(s1, s2, config)),
        ('global', lambda s1, s2: calculate_pairwise_alignment_similarity(s1, s2, config)),
        ('local', lambda s1, s2: calculate_local_alignment_similarity(s1, s2, config))
    ]
    
    for method_name, method_func in methods:
        print(f"Benchmarking {method_name}...")
        repeat_times = []
        repeat_scores = []
        
        for repeat in range(num_repeats):
            print(f"  Repeat {repeat + 1}/{num_repeats}...")
            
            # Prepare test sequences for this repeat
            test_pairs = []
            for i in range(num_samples):
                idx = random.randint(0, len(index_seqs) - 1)
                test_pairs.append((data_seq, index_seqs[idx]))
            
            try:
                # Run the method multiple times for accurate timing
                start_time = time.time()
                total_score = 0.0
                
                for _ in range(iterations_per_test):
                    for seq1, seq2 in test_pairs:
                        score = method_func(seq1, seq2)
                        total_score += score
                
                elapsed_time = time.time() - start_time
                avg_score = total_score / (iterations_per_test * num_samples)
                
                repeat_times.append(elapsed_time)
                repeat_scores.append(avg_score)
                
            except Exception as e:
                print(f"    Error in {method_name}: {e}")
                repeat_times.append(float('inf'))
                repeat_scores.append(0.0)
        
        # Calculate statistics
        import statistics
        avg_time = statistics.mean(repeat_times) if repeat_times else float('inf')
        std_time = statistics.stdev(repeat_times) if len(repeat_times) > 1 else 0.0
        avg_score = statistics.mean(repeat_scores) if repeat_scores else 0.0
        min_score = min(repeat_scores) if repeat_scores else 0.0
        max_score = max(repeat_scores) if repeat_scores else 0.0
        
        # Calculate time per single operation
        time_per_op = avg_time / (iterations_per_test * num_samples) if avg_time != float('inf') else float('inf')
        
        results[method_name] = {
            'avg_time': avg_time,
            'std_time': std_time,
            'time_per_op': time_per_op,
            'avg_score': avg_score,
            'min_score': min_score,
            'max_score': max_score,
            'repeat_times': repeat_times,
            'repeat_scores': repeat_scores
        }
    
    # Print results
    print("\n" + "="*100)
    print("BENCHMARK RESULTS")
    print("="*100)
    print(f"{'Method':<12} {'Avg Time (s)':<12} {'Std Dev (s)':<12} {'Time/Op (ms)':<12} {'Avg Score':<10} {'Min Score':<10} {'Max Score':<10} {'Speed Rank':<10}")
    print("-" * 100)
    
    # Sort by average time for ranking
    sorted_methods = sorted(results.items(), key=lambda x: x[1]['avg_time'])
    
    for rank, (method_name, data) in enumerate(sorted_methods, 1):
        time_per_op_ms = data['time_per_op'] * 1000 if data['time_per_op'] != float('inf') else float('inf')
        print(f"{method_name:<12} {data['avg_time']:<12.4f} {data['std_time']:<12.4f} {time_per_op_ms:<12.2f} "
              f"{data['avg_score']:<10.3f} {data['min_score']:<10.3f} {data['max_score']:<10.3f} {rank:<10}")
    
    print("\n" + "="*100)
    print("RECOMMENDATIONS")
    print("="*100)
    
    # Find fastest methods
    fast_methods = [name for name, data in sorted_methods[:3]]
    print(f"Fastest methods: {', '.join(fast_methods)}")
    
    # Find methods with good score range
    good_range_methods = []
    for name, data in results.items():
        if data['max_score'] - data['min_score'] > 0.1:  # Good discrimination
            good_range_methods.append(name)
    
    if good_range_methods:
        print(f"Methods with good discrimination: {', '.join(good_range_methods)}")
    
    # Specific recommendations
    print("\nRecommendations:")
    
    # For speed
    fastest = sorted_methods[0][0]
    print(f"• For maximum speed: {fastest}")
    
    # For balance
    balanced_candidates = []
    for name, data in sorted_methods:
        if data['time_per_op'] < 0.001 and (data['max_score'] - data['min_score']) > 0.05:
            balanced_candidates.append(name)
    
    if balanced_candidates:
        print(f"• For speed/accuracy balance: {balanced_candidates[0]}")
    else:
        print(f"• For speed/accuracy balance: {sorted_methods[1][0] if len(sorted_methods) > 1 else fastest}")
    
    # For accuracy (methods that use full sequence information)
    accuracy_methods = ['global', 'local', 'local_fast']
    best_accuracy = None
    for name in accuracy_methods:
        if name in results:
            best_accuracy = name
            break
    
    if best_accuracy:
        print(f"• For maximum accuracy: {best_accuracy}")
    
    return results

def quick_benchmark(num_tests_per_length=20, iterations_per_test=100, num_repeats=3):
    """
    Quick benchmark using synthetic sequences with detailed results for each sequence length.
    """
    print("=== Quick Synthetic Benchmark ===")
    print(f"Number of tests per sequence length: {num_tests_per_length}")
    print(f"Iterations per test: {iterations_per_test}")
    print(f"Number of repeats: {num_repeats}")
    
    # Generate synthetic sequences
    def generate_seq(length):
        return ''.join(random.choice('ATCG') for _ in range(length))
    
    def mutate_dna(dna_sequence: str, mutation_rate: float = 0.05) -> str:
        # """
        # Introduces random point mutations into a DNA sequence at a given rate.

        # A point mutation replaces a single base with one of the other three bases.
        # For example, 'A' can be mutated to 'T', 'C', or 'G'.

        # Args:
        #     dna_sequence (str): The original DNA sequence, consisting of 'A', 'T', 
        #                         'C', and 'G'.
        #     mutation_rate (float): The probability (between 0.0 and 1.0) that any
        #                            given base will mutate.

        # Returns:
        #     str: A new DNA sequence with mutations introduced.

        # Raises:
        #     ValueError: If the mutation_rate is not between 0.0 and 1.0.
        #     ValueError: If the dna_sequence contains characters other than 'A', 'T', 'C', 'G'.
        # """
        # 1. Input Validation
        if not 0.0 <= mutation_rate <= 1.0:
            raise ValueError("Mutation rate must be between 0.0 and 1.0.")

        valid_bases = set('ATCG')
        if not set(dna_sequence.upper()).issubset(valid_bases):
            raise ValueError("Input DNA sequence contains invalid characters.")

        # 2. Mutation Logic
        bases = ['A', 'T', 'C', 'G']
        mutated_sequence = []

        for base in dna_sequence.upper():
            # Decide if a mutation should occur for this base
            if random.random() < mutation_rate:
                # It's a mutation! Choose a new base.
                # Create a list of possible new bases (all except the current one)
                possible_mutations = [b for b in bases if b != base]
                # Randomly select one of the possible new bases
                new_base = random.choice(possible_mutations)
                mutated_sequence.append(new_base)
            else:
                # No mutation, keep the original base
                mutated_sequence.append(base)

        # 3. Return the new sequence
        return "".join(mutated_sequence)
    
    # Test different sequence lengths
    test_lengths = [100, 300, 600, 1000]
    config = AlignmentConfig()
    
    # Test each method
    methods = [
        ('kmer', lambda s1, s2: calculate_kmer_similarity(s1, s2, k=13)),
        ('sliding', lambda s1, s2: calculate_sliding_window_similarity(s1, s2)),
        ('lcs', lambda s1, s2: calculate_lcs_similarity(s1, s2)),
        ('local_fast', lambda s1, s2: calculate_local_alignment_similarity_fast(s1, s2, config)),
        ('global', lambda s1, s2: calculate_pairwise_alignment_similarity(s1, s2, config)),
        ('local', lambda s1, s2: calculate_local_alignment_similarity(s1, s2, config))
    ]
    
    # Store results by length and method
    results_by_length = {}
    overall_results = {}
    
    for length in test_lengths:
        print(f"\nTesting sequence length: {length} bp")
        results_by_length[length] = {}
        
        for method_name, method_func in methods:
            print(f"  Benchmarking {method_name}...")
            repeat_times = []
            repeat_scores = []
            
            for repeat in range(num_repeats):
                print(f"    Repeat {repeat + 1}/{num_repeats}...")
                
                # Prepare test sequences for this repeat
                test_pairs = []
                for test_num in range(num_tests_per_length):
                    seq1 = generate_seq(length)
                    seq2 = mutate_dna(seq1, mutation_rate=0.05)
                    test_pairs.append((seq1, seq2))
                
                try:
                    # Run the method multiple times for accurate timing
                    start_time = time.time()
                    total_score = 0.0
                    
                    for _ in range(iterations_per_test):
                        for seq1, seq2 in test_pairs:
                            score = method_func(seq1, seq2)
                            total_score += score
                    
                    elapsed_time = time.time() - start_time
                    avg_score = total_score / (iterations_per_test * num_tests_per_length)
                    
                    repeat_times.append(elapsed_time)
                    repeat_scores.append(avg_score)
                    
                except Exception as e:
                    print(f"      Error in {method_name} at length {length}: {e}")
                    repeat_times.append(float('inf'))
                    repeat_scores.append(0.0)
            
            # Calculate statistics for this length
            import statistics
            avg_time = statistics.mean(repeat_times) if repeat_times else float('inf')
            std_time = statistics.stdev(repeat_times) if len(repeat_times) > 1 else 0.0
            avg_score = statistics.mean(repeat_scores) if repeat_scores else 0.0
            min_score = min(repeat_scores) if repeat_scores else 0.0
            max_score = max(repeat_scores) if repeat_scores else 0.0
            
            # Calculate time per single operation
            total_ops = len(repeat_times) * iterations_per_test * num_tests_per_length
            time_per_op = avg_time / total_ops if avg_time != float('inf') else float('inf')
            
            results_by_length[length][method_name] = {
                'avg_time': avg_time,
                'std_time': std_time,
                'time_per_op': time_per_op,
                'avg_score': avg_score,
                'min_score': min_score,
                'max_score': max_score,
                'repeat_times': repeat_times,
                'repeat_scores': repeat_scores
            }
            
            # Add to overall results
            if method_name not in overall_results:
                overall_results[method_name] = {
                    'repeat_times': [],
                    'repeat_scores': []
                }
            overall_results[method_name]['repeat_times'].extend(repeat_times)
            overall_results[method_name]['repeat_scores'].extend(repeat_scores)
    
    # Calculate overall statistics
    for method_name in overall_results:
        repeat_times = overall_results[method_name]['repeat_times']
        repeat_scores = overall_results[method_name]['repeat_scores']
        
        avg_time = statistics.mean(repeat_times) if repeat_times else float('inf')
        std_time = statistics.stdev(repeat_times) if len(repeat_times) > 1 else 0.0
        avg_score = statistics.mean(repeat_scores) if repeat_scores else 0.0
        min_score = min(repeat_scores) if repeat_scores else 0.0
        max_score = max(repeat_scores) if repeat_scores else 0.0
        
        # Calculate time per single operation
        total_ops = len(repeat_times) * iterations_per_test * num_tests_per_length
        time_per_op = avg_time / total_ops if avg_time != float('inf') else float('inf')
        
        overall_results[method_name].update({
            'avg_time': avg_time,
            'std_time': std_time,
            'time_per_op': time_per_op,
            'avg_score': avg_score,
            'min_score': min_score,
            'max_score': max_score
        })
    
    # Print detailed results for each sequence length
    print("\n" + "="*120)
    print("DETAILED RESULTS BY SEQUENCE LENGTH")
    print("="*120)
    
    for length in test_lengths:
        print(f"\nSequence Length: {length} bp")
        print("-" * 110)
        print(f"{'Method':<12} {'Avg Time (s)':<12} {'Std Dev (s)':<12} {'Time/Op (ms)':<12} {'Avg Score':<10} {'Min Score':<10} {'Max Score':<10} {'Speed Rank':<10}")
        print("-" * 110)
        
        # Sort by average time for ranking
        length_results = results_by_length[length]
        sorted_methods = sorted(length_results.items(), key=lambda x: x[1]['avg_time'])
        
        for rank, (method_name, data) in enumerate(sorted_methods, 1):
            time_per_op_ms = data['time_per_op'] * 1000 if data['time_per_op'] != float('inf') else float('inf')
            print(f"{method_name:<12} {data['avg_time']:<12.4f} {data['std_time']:<12.4f} {time_per_op_ms:<12.2f} "
                  f"{data['avg_score']:<10.3f} {data['min_score']:<10.3f} {data['max_score']:<10.3f} {rank:<10}")
        
        # Show fastest method for this length
        fastest_method = sorted_methods[0][0]
        print(f"Fastest for {length}bp: {fastest_method}")
    
    # Print overall summary
    print("\n" + "="*120)
    print("OVERALL SUMMARY")
    print("="*120)
    total_tests = len(test_lengths) * num_tests_per_length * iterations_per_test * num_repeats
    print(f"Total operations performed: {total_tests} per method")
    print()
    
    print(f"{'Method':<12} {'Avg Time (s)':<12} {'Std Dev (s)':<12} {'Time/Op (ms)':<12} {'Avg Score':<10} {'Min Score':<10} {'Max Score':<10} {'Speed Rank':<10}")
    print("-" * 110)
    
    # Sort by average time for ranking
    sorted_overall = sorted(overall_results.items(), key=lambda x: x[1]['avg_time'])
    
    for rank, (method_name, data) in enumerate(sorted_overall, 1):
        time_per_op_ms = data['time_per_op'] * 1000 if data['time_per_op'] != float('inf') else float('inf')
        print(f"{method_name:<12} {data['avg_time']:<12.4f} {data['std_time']:<12.4f} {time_per_op_ms:<12.2f} "
              f"{data['avg_score']:<10.3f} {data['min_score']:<10.3f} {data['max_score']:<10.3f} {rank:<10}")
    
    # Performance analysis by sequence length
    print("\n" + "="*120)
    print("PERFORMANCE ANALYSIS")
    print("="*120)
    
    print("\nFastest method by sequence length:")
    for length in test_lengths:
        length_results = results_by_length[length]
        fastest = min(length_results.items(), key=lambda x: x[1]['avg_time'])
        print(f"  {length:4d}bp: {fastest[0]:<12} ({fastest[1]['avg_time']:.4f}s ± {fastest[1]['std_time']:.4f}s)")
    
    print("\nMethod performance scaling:")
    print(f"{'Method':<12} {'100bp (ms)':<11} {'300bp (ms)':<11} {'600bp (ms)':<11} {'1000bp (ms)':<12} {'Scaling':<10}")
    print("-" * 75)
    
    for method_name in methods:
        method = method_name[0]
        times_by_length = [results_by_length[length][method]['time_per_op'] * 1000 for length in test_lengths]
        scaling = times_by_length[-1] / times_by_length[0] if times_by_length[0] > 0 else float('inf')
        
        print(f"{method:<12} {times_by_length[0]:<11.2f} {times_by_length[1]:<11.2f} "
              f"{times_by_length[2]:<11.2f} {times_by_length[3]:<12.2f} {scaling:<10.1f}x")
    
    print("\n" + "="*120)
    print("RECOMMENDATIONS")
    print("="*120)
    
    # Find fastest methods overall
    fast_methods = [name for name, data in sorted_overall[:3]]
    print(f"Fastest methods overall: {', '.join(fast_methods)}")
    
    # Find methods with good score range
    good_range_methods = []
    for name, data in overall_results.items():
        if data['max_score'] - data['min_score'] > 0.1:  # Good discrimination
            good_range_methods.append(name)
    
    if good_range_methods:
        print(f"Methods with good discrimination: {', '.join(good_range_methods)}")
    
    # Specific recommendations
    print("\nRecommendations:")
    
    # For speed
    fastest = sorted_overall[0][0]
    print(f"• For maximum speed: {fastest}")
    
    # For balance
    balanced_candidates = []
    for name, data in sorted_overall:
        if data['time_per_op'] < 0.001 and (data['max_score'] - data['min_score']) > 0.05:
            balanced_candidates.append(name)
    
    if balanced_candidates:
        print(f"• For speed/accuracy balance: {balanced_candidates[0]}")
    else:
        print(f"• For speed/accuracy balance: {sorted_overall[1][0] if len(sorted_overall) > 1 else fastest}")
    
    # For accuracy (methods that use full sequence information)
    accuracy_methods = ['global', 'local', 'local_fast']
    best_accuracy = None
    for name in accuracy_methods:
        if name in overall_results:
            best_accuracy = name
            break
    
    if best_accuracy:
        print(f"• For maximum accuracy: {best_accuracy}")
    
    # Length-specific recommendations
    print(f"\nLength-specific recommendations:")
    for length in test_lengths:
        length_results = results_by_length[length]
        fastest = min(length_results.items(), key=lambda x: x[1]['avg_time'])
        print(f"• For {length}bp sequences: {fastest[0]}")
    
    return overall_results, results_by_length

def main():
    """
    Main function for running benchmarks.
    """
    if len(sys.argv) < 2:
        print("Usage:")
        print("  python benchmark_similarity.py quick [num_tests] [iterations] [repeats]     # Quick synthetic benchmark")
        print("  python benchmark_similarity.py <data_seq> <index_seqs> [num_samples] [iterations] [repeats]  # Real data benchmark")
        print()
        print("Examples:")
        print("  python benchmark_similarity.py quick                        # Quick benchmark with default settings")
        print("  python benchmark_similarity.py quick 50                     # Quick benchmark with 50 tests per length")
        print("  python benchmark_similarity.py quick 100 1 3               # Quick benchmark with 100 tests, 1 iterations, 3 repeats")
        print("  python benchmark_similarity.py input_files/data-seq.fa input_files/64-bit-index-seqs.fa")
        print("  python benchmark_similarity.py input_files/data-seq.fa input_files/64-bit-index-seqs.fa 100 1 3")
        sys.exit(1)
    
    if sys.argv[1] == 'quick':
        num_tests = int(sys.argv[2]) if len(sys.argv) > 2 else 20
        iterations = int(sys.argv[3]) if len(sys.argv) > 3 else 100
        repeats = int(sys.argv[4]) if len(sys.argv) > 4 else 3
        quick_benchmark(num_tests, iterations, repeats)
    else:
        if len(sys.argv) < 3:
            print("Error: Need both data sequence file and index sequence file")
            sys.exit(1)
        
        data_seq_file = sys.argv[1]
        index_seq_file = sys.argv[2]
        num_samples = int(sys.argv[3]) if len(sys.argv) > 3 else 10
        iterations = int(sys.argv[4]) if len(sys.argv) > 4 else 100
        repeats = int(sys.argv[5]) if len(sys.argv) > 5 else 3
        
        benchmark_on_real_data(data_seq_file, index_seq_file, num_samples, iterations, repeats)

if __name__ == "__main__":
    main() 