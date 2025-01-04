import math
import argparse
from scipy.stats import binom

# Function to parse the command line argument
def parse_args():
    parser = argparse.ArgumentParser(description='Calculate binomial cumulative distribution.')
    # Set default value for E and allow optional override from command line
    parser.add_argument('-E', type=float, default=0.00372166666666667, help='The value of E to use in the calculations (default: 0.00372166666666667)')
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
