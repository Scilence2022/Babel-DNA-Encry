"""
Alignment-based utility functions for Z-DNA deciphering.
This module provides sequence alignment-based alternatives to k-mer Jaccard similarity methods.
"""

import time
from Bio import pairwise2
from Bio.Seq import Seq
import numpy as np
from difflib import SequenceMatcher
import re

class AlignmentConfig:
    """Configuration class for alignment parameters."""
    def __init__(self, 
                 match_score=2, 
                 mismatch_score=-1, 
                 gap_open=-2, 
                 gap_extend=-0.5):
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_open = gap_open
        self.gap_extend = gap_extend

def calculate_kmer_similarity(seq1, seq2, k=6):
    """
    Fast k-mer based similarity calculation.
    
    Parameters:
        seq1 (str): First DNA sequence
        seq2 (str): Second DNA sequence
        k (int): K-mer length
    
    Returns:
        float: Jaccard similarity score [0, 1]
    """
    if not seq1 or not seq2:
        return 0.0
    
    # Generate k-mers
    kmers1 = set()
    kmers2 = set()
    
    for i in range(len(seq1) - k + 1):
        kmers1.add(seq1[i:i+k])
    
    for i in range(len(seq2) - k + 1):
        kmers2.add(seq2[i:i+k])
    
    if not kmers1 or not kmers2:
        return 0.0
    
    # Calculate Jaccard similarity
    intersection = len(kmers1.intersection(kmers2))
    union = len(kmers1.union(kmers2))
    
    return intersection / union if union > 0 else 0.0

def calculate_sliding_window_similarity(seq1, seq2, window_size=50, step=10):
    """
    Fast sliding window similarity calculation.
    
    Parameters:
        seq1 (str): First DNA sequence
        seq2 (str): Second DNA sequence
        window_size (int): Size of sliding window
        step (int): Step size for sliding window
    
    Returns:
        float: Maximum similarity score [0, 1]
    """
    if not seq1 or not seq2:
        return 0.0
    
    max_similarity = 0.0
    
    # Slide window across seq1
    for i in range(0, len(seq1) - window_size + 1, step):
        window1 = seq1[i:i + window_size]
        
        # Slide window across seq2
        for j in range(0, len(seq2) - window_size + 1, step):
            window2 = seq2[j:j + window_size]
            
            # Calculate simple match ratio
            matches = sum(1 for a, b in zip(window1, window2) if a == b)
            similarity = matches / window_size
            
            max_similarity = max(max_similarity, similarity)
    
    return max_similarity

def calculate_lcs_similarity(seq1, seq2):
    """
    Calculate similarity using Longest Common Subsequence (LCS).
    
    Parameters:
        seq1 (str): First DNA sequence
        seq2 (str): Second DNA sequence
    
    Returns:
        float: LCS-based similarity score [0, 1]
    """
    if not seq1 or not seq2:
        return 0.0
    
    # Use SequenceMatcher for fast LCS-like calculation
    matcher = SequenceMatcher(None, seq1, seq2)
    similarity = matcher.ratio()
    
    return similarity

def calculate_edit_distance_similarity(seq1, seq2, max_distance=None):
    """
    Calculate similarity based on edit distance with early termination.
    
    Parameters:
        seq1 (str): First DNA sequence
        seq2 (str): Second DNA sequence
        max_distance (int): Maximum allowed edit distance (for early termination)
    
    Returns:
        float: Edit distance-based similarity score [0, 1]
    """
    if not seq1 or not seq2:
        return 0.0
    
    len1, len2 = len(seq1), len(seq2)
    max_len = max(len1, len2)
    
    # Early termination if length difference is too large
    if max_distance and abs(len1 - len2) > max_distance:
        return 0.0
    
    # Use a simplified edit distance calculation
    if len1 > len2:
        seq1, seq2 = seq2, seq1
        len1, len2 = len2, len1
    
    # Create distance matrix (only current and previous row needed)
    prev_row = list(range(len2 + 1))
    
    for i in range(1, len1 + 1):
        curr_row = [i] + [0] * len2
        
        for j in range(1, len2 + 1):
            if seq1[i-1] == seq2[j-1]:
                curr_row[j] = prev_row[j-1]
            else:
                curr_row[j] = 1 + min(prev_row[j],      # deletion
                                    curr_row[j-1],    # insertion
                                    prev_row[j-1])    # substitution
            
            # Early termination if distance exceeds threshold
            if max_distance and curr_row[j] > max_distance:
                return 0.0
        
        prev_row = curr_row
    
    edit_distance = prev_row[len2]
    similarity = 1.0 - (edit_distance / max_len)
    
    return max(0.0, similarity)

def calculate_local_alignment_similarity_fast(seq1, seq2, config=None):
    """
    Fast local alignment using BioPython with optimized parameters.
    
    Parameters:
        seq1 (str): First DNA sequence
        seq2 (str): Second DNA sequence
        config (AlignmentConfig): Alignment scoring parameters
    
    Returns:
        float: Normalized similarity score [0, 1]
    """
    if config is None:
        config = AlignmentConfig()
    
    if not seq1 or not seq2:
        return 0.0
    
    # Truncate very long sequences for speed
    max_len = 300
    if len(seq1) > max_len:
        seq1 = seq1[:max_len]
    if len(seq2) > max_len:
        seq2 = seq2[:max_len]
    
    # Use local alignment
    alignments = pairwise2.align.localms(
        seq1, seq2,
        config.match_score, config.mismatch_score,
        config.gap_open, config.gap_extend,
        one_alignment_only=True
    )
    
    if not alignments:
        return 0.0
    
    alignment = alignments[0]
    score = alignment.score
    
    # Better normalization for local alignment
    max_possible_score = min(len(seq1), len(seq2)) * config.match_score
    
    # Normalize similarity to [0, 1] range
    similarity = max(0.0, score / max_possible_score) if max_possible_score > 0 else 0.0
    return similarity

def calculate_hybrid_similarity(seq1, seq2, config=None, method='auto'):
    """
    Hybrid similarity calculation using multiple methods.
    
    Parameters:
        seq1 (str): First DNA sequence
        seq2 (str): Second DNA sequence
        config (AlignmentConfig): Alignment scoring parameters
        method (str): 'auto', 'kmer', 'sliding', 'lcs', 'edit', 'local'
    
    Returns:
        float: Similarity score [0, 1]
    """
    if not seq1 or not seq2:
        return 0.0
    
    if method == 'auto':
        # Choose method based on sequence length
        avg_len = (len(seq1) + len(seq2)) / 2
        if avg_len > 1000:
            method = 'kmer'
        elif avg_len > 500:
            method = 'sliding'
        else:
            method = 'lcs'
    
    if method == 'kmer':
        return calculate_kmer_similarity(seq1, seq2, k=6)
    elif method == 'sliding':
        return calculate_sliding_window_similarity(seq1, seq2, window_size=50, step=10)
    elif method == 'lcs':
        return calculate_lcs_similarity(seq1, seq2)
    elif method == 'edit':
        return calculate_edit_distance_similarity(seq1, seq2, max_distance=100)
    elif method == 'local':
        return calculate_local_alignment_similarity_fast(seq1, seq2, config)
    else:
        # Default to k-mer for speed
        return calculate_kmer_similarity(seq1, seq2, k=6)

def calculate_pairwise_alignment_similarity(seq1, seq2, config=None):
    """
    Calculate similarity between two sequences using pairwise global alignment.
    
    Parameters:
        seq1 (str): First DNA sequence
        seq2 (str): Second DNA sequence
        config (AlignmentConfig): Alignment scoring parameters
    
    Returns:
        float: Normalized similarity score [0, 1]
    """
    if config is None:
        config = AlignmentConfig()
    
    if not seq1 or not seq2:
        return 0.0
    
    # Use global alignment with custom scoring
    alignments = pairwise2.align.globalms(
        seq1, seq2,
        config.match_score, config.mismatch_score,
        config.gap_open, config.gap_extend,
        one_alignment_only=True
    )
    
    if not alignments:
        return 0.0
    
    alignment = alignments[0]
    score = alignment.score
    max_possible_score = max(len(seq1), len(seq2)) * config.match_score
    
    # Normalize similarity to [0, 1] range
    similarity = max(0.0, score / max_possible_score) if max_possible_score > 0 else 0.0
    return similarity

def calculate_local_alignment_similarity(seq1, seq2, config=None):
    """
    Calculate similarity between two sequences using local alignment (Smith-Waterman).
    
    Parameters:
        seq1 (str): First DNA sequence
        seq2 (str): Second DNA sequence
        config (AlignmentConfig): Alignment scoring parameters
    
    Returns:
        float: Normalized similarity score [0, 1]
    """
    if config is None:
        config = AlignmentConfig()
    
    if not seq1 or not seq2:
        return 0.0
    
    # Use local alignment
    alignments = pairwise2.align.localms(
        seq1, seq2,
        config.match_score, config.mismatch_score,
        config.gap_open, config.gap_extend,
        one_alignment_only=True
    )
    
    if not alignments:
        return 0.0
    
    alignment = alignments[0]
    score = alignment.score
    max_possible_score = min(len(seq1), len(seq2)) * config.match_score
    
    # Normalize similarity to [0, 1] range
    similarity = max(0.0, score / max_possible_score) if max_possible_score > 0 else 0.0
    return similarity

def find_best_alignment_match(query_seq, reference_seqs, config=None, alignment_type='hybrid'):
    """
    Find the best matching reference sequence for a query sequence using alignment.
    
    Parameters:
        query_seq (str): Query DNA sequence
        reference_seqs (list): List of reference DNA sequences
        config (AlignmentConfig): Alignment scoring parameters
        alignment_type (str): 'global', 'local', 'hybrid', 'kmer', 'sliding', 'lcs', 'edit'
    
    Returns:
        tuple: (best_index, best_score) where best_index is the index of the best match
               and best_score is the similarity score
    """
    if config is None:
        config = AlignmentConfig()
    
    best_index = -1
    best_score = 0.0
    
    # Choose similarity function based on alignment type
    if alignment_type == 'global':
        similarity_func = lambda s1, s2: calculate_pairwise_alignment_similarity(s1, s2, config)
    elif alignment_type == 'local':
        similarity_func = lambda s1, s2: calculate_local_alignment_similarity(s1, s2, config)
    elif alignment_type == 'kmer':
        similarity_func = lambda s1, s2: calculate_kmer_similarity(s1, s2, k=6)
    elif alignment_type == 'sliding':
        similarity_func = lambda s1, s2: calculate_sliding_window_similarity(s1, s2)
    elif alignment_type == 'lcs':
        similarity_func = lambda s1, s2: calculate_lcs_similarity(s1, s2)
    elif alignment_type == 'edit':
        similarity_func = lambda s1, s2: calculate_edit_distance_similarity(s1, s2)
    else:  # hybrid or default
        similarity_func = lambda s1, s2: calculate_hybrid_similarity(s1, s2, config)
    
    for i, ref_seq in enumerate(reference_seqs):
        score = similarity_func(query_seq, ref_seq)
        if score > best_score:
            best_score = score
            best_index = i
    
    return best_index, best_score

def cluster_sequences_by_alignment(sequences, reference_seqs, config=None, 
                                 threshold=0.35, alignment_type='hybrid'):
    """
    Cluster sequences based on their best alignment match to reference sequences.
    
    Parameters:
        sequences (list): List of sequences to cluster
        reference_seqs (list): List of reference sequences (one per cluster)
        config (AlignmentConfig): Alignment scoring parameters
        threshold (float): Minimum similarity threshold for assignment
        alignment_type (str): Type of alignment/similarity method
    
    Returns:
        list: List of lists, where each sublist contains sequences belonging to that cluster
    """
    if config is None:
        config = AlignmentConfig()
    
    clusters = [[] for _ in range(len(reference_seqs))]
    unassigned = []
    
    for seq in sequences:
        best_idx, best_score = find_best_alignment_match(
            seq, reference_seqs, config, alignment_type
        )
        
        if best_score >= threshold:
            clusters[best_idx].append(seq)
        else:
            unassigned.append(seq)
    
    return clusters, unassigned

def calculate_consensus_similarity(sequences, reference_seq, config=None, 
                                 alignment_type='hybrid', method='average'):
    """
    Calculate consensus similarity between a list of sequences and a reference.
    
    Parameters:
        sequences (list): List of DNA sequences
        reference_seq (str): Reference DNA sequence
        config (AlignmentConfig): Alignment scoring parameters
        alignment_type (str): Type of alignment/similarity method
        method (str): 'average', 'median', or 'max' for consensus calculation
    
    Returns:
        float: Consensus similarity score
    """
    if config is None:
        config = AlignmentConfig()
    
    if not sequences:
        return 0.0
    
    # Choose similarity function
    if alignment_type == 'global':
        similarity_func = lambda s1, s2: calculate_pairwise_alignment_similarity(s1, s2, config)
    elif alignment_type == 'local':
        similarity_func = lambda s1, s2: calculate_local_alignment_similarity(s1, s2, config)
    elif alignment_type == 'kmer':
        similarity_func = lambda s1, s2: calculate_kmer_similarity(s1, s2, k=6)
    else:  # hybrid or default
        similarity_func = lambda s1, s2: calculate_hybrid_similarity(s1, s2, config)
    
    similarities = []
    for seq in sequences:
        sim = similarity_func(seq, reference_seq)
        similarities.append(sim)
    
    if method == 'average':
        return np.mean(similarities)
    elif method == 'median':
        return np.median(similarities)
    elif method == 'max':
        return np.max(similarities)
    else:
        raise ValueError("Method must be 'average', 'median', or 'max'")

def filter_sequences_by_alignment(sequences, reference_seq, config=None, 
                                num_sequences=5, alignment_type='hybrid'):
    """
    Filter sequences by selecting those with moderate alignment similarity to reference.
    
    Parameters:
        sequences (list): List of DNA sequences to filter
        reference_seq (str): Reference sequence
        config (AlignmentConfig): Alignment scoring parameters
        num_sequences (int): Number of sequences to return
        alignment_type (str): Type of alignment/similarity method
    
    Returns:
        list: Filtered list of sequences
    """
    if config is None:
        config = AlignmentConfig()
    
    if len(sequences) <= num_sequences:
        return sequences
    
    # Choose similarity function
    if alignment_type == 'global':
        similarity_func = lambda s1, s2: calculate_pairwise_alignment_similarity(s1, s2, config)
    elif alignment_type == 'local':
        similarity_func = lambda s1, s2: calculate_local_alignment_similarity(s1, s2, config)
    elif alignment_type == 'kmer':
        similarity_func = lambda s1, s2: calculate_kmer_similarity(s1, s2, k=6)
    else:  # hybrid or default
        similarity_func = lambda s1, s2: calculate_hybrid_similarity(s1, s2, config)
    
    # Calculate similarities
    seq_similarities = []
    for seq in sequences:
        sim = similarity_func(seq, reference_seq)
        seq_similarities.append((seq, sim))
    
    # Sort by similarity
    seq_similarities.sort(key=lambda x: x[1])
    
    # Select middle sequences (avoiding extreme values)
    total_seqs = len(seq_similarities)
    start_idx = (total_seqs - num_sequences) // 2
    end_idx = start_idx + num_sequences
    
    selected_sequences = [seq for seq, _ in seq_similarities[start_idx:end_idx]]
    return selected_sequences

def benchmark_alignment_methods(seq1, seq2, config=None):
    """
    Benchmark different alignment methods and return timing and scores.
    
    Parameters:
        seq1 (str): First DNA sequence
        seq2 (str): Second DNA sequence
        config (AlignmentConfig): Alignment scoring parameters
    
    Returns:
        dict: Dictionary with timing and score information
    """
    if config is None:
        config = AlignmentConfig()
    
    results = {}
    methods = [
        ('kmer', lambda: calculate_kmer_similarity(seq1, seq2)),
        ('sliding', lambda: calculate_sliding_window_similarity(seq1, seq2)),
        ('lcs', lambda: calculate_lcs_similarity(seq1, seq2)),
        ('edit', lambda: calculate_edit_distance_similarity(seq1, seq2)),
        ('local_fast', lambda: calculate_local_alignment_similarity_fast(seq1, seq2, config)),
        ('hybrid', lambda: calculate_hybrid_similarity(seq1, seq2, config)),
        ('global', lambda: calculate_pairwise_alignment_similarity(seq1, seq2, config)),
        ('local', lambda: calculate_local_alignment_similarity(seq1, seq2, config))
    ]
    
    for method_name, method_func in methods:
        try:
            start_time = time.time()
            score = method_func()
            elapsed_time = time.time() - start_time
            
            results[method_name] = {
                'score': score,
                'time': elapsed_time
            }
        except Exception as e:
            results[method_name] = {
                'score': 0.0,
                'time': float('inf'),
                'error': str(e)
            }
    
    return results

def read_fasta(file_path):
    """
    Read sequences from a FASTA file.
    
    Parameters:
        file_path (str): Path to FASTA file
    
    Returns:
        list: List of DNA sequences
    """
    sequences = []
    with open(file_path, 'r') as f:
        current_seq = ""
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    sequences.append(current_seq)
                    current_seq = ""
            else:
                current_seq += line
        if current_seq:
            sequences.append(current_seq)
    return sequences 