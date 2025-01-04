import math
import argparse
from scipy.stats import binom

# Function to parse the command line arguments
def parse_args():
    parser = argparse.ArgumentParser(
        description='Calculate binomial cumulative distribution.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False  # Disable the default help
    )
    # Manually add the help option
    parser.add_argument(
        '-h', '--help',
        action='help',
        default=argparse.SUPPRESS,
        help='Show this help message and exit.'
    )
    # Add other arguments
    parser.add_argument(
        '-E',
        type=float,
        default=0.00372166666666667,
        help='The value of E to use in the calculations'
    )
    return parser.parse_args()

def main():
    # Parse the command line arguments
    args = parse_args()
    E = args.E  # Assign the E value passed from the command line or default value

    print("N\tm\tE\tEn")

    for i in range(1, 11):
        N = i * 2 + 1
        m = math.ceil((N + 1) / 2)

        print(N, end='\t')
        print(m, end='\t')
        print(E, end='\t')
        print(1 - binom.cdf(m - 1, N, E))

if __name__ == '__main__':
    main() 