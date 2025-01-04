import argparse
from test_utils import load_variable_from_gzfile, run_mul_retri_bits, collect_bits_based_on_key

def main():
    parser = argparse.ArgumentParser(
        description='Simulate error rates for Z-DNA key decoding.'
    )
    
    parser.add_argument(
        '-i', '--bit_dec_file_path', 
        type=str, 
        default='input_files/32-Bit-5-of-100Seqs-2E5-Decs.gz',
        help='Path to the bit decoding file (default: input_files/32-Bit-5-of-100Seqs-2E5-Decs.gz)'
    )
    
    parser.add_argument(
        '--rep_size', 
        type=int, 
        default=10000,
        help='Number of repetitions for each simulation (default: 10000)'
    )
    
    parser.add_argument(
        '--m_range', 
        type=int, 
        default=13,
        help='Range for m in simulations (default: 13)'
    )
    
    parser.add_argument(
        '--64B', 
        action='store_true',
        help='Use 64-bit key instead of the default 32-bit key'
    )
    
    args = parser.parse_args()
    
    # Define the correct_bit_key based on the --64B flag
    if args.__dict__['64B']:
        correct_bit_key = [
            1, 1, 0, 0, 0, 0, 1, 0, 
            0, 0, 0, 0, 0, 1, 0, 0, 
            0, 1, 1, 1, 0, 0, 1, 1, 
            0, 0, 1, 0, 1, 0, 0, 1, 
            1, 1, 0, 0, 1, 0, 0, 1, 
            1, 0, 0, 0, 1, 0, 0, 0, 
            1, 0, 0, 1, 1, 1, 0, 1, 
            1, 0, 1, 0, 1, 1, 1, 1
        ]
        print("Using 64-Bit key.")
    else:
        correct_bit_key = [
            1, 0, 1, 0, 1, 1, 1, 0, 
            0, 0, 1, 1, 1, 1, 0, 1, 
            1, 0, 0, 0, 0, 1, 0, 1, 
            0, 0, 1, 1, 0, 1, 1, 0
        ]
        print("Using 32-Bit key.")
    
    # Load the bit decoding array from the specified file
    print(f"Loading bit decoding data from: {args.bit_dec_file_path}")
    dec_bit_arr = load_variable_from_gzfile(args.bit_dec_file_path)
    
    # Collect bits based on the correct_bit_key
    bit_dec_values = collect_bits_based_on_key(dec_bit_arr, correct_bit_key)
    
    # Run multiple retrieval simulations
    result = []
    for i in range(3, args.m_range+1, 2):
        aa = run_mul_retri_bits(bit_dec_values, i, args.rep_size)
        result.append(aa)
        print(f"Number of Multi-Retrievals: {i}\t\tNumber of correct decodings: {aa}/{args.rep_size}")

if __name__ == "__main__":
    main()