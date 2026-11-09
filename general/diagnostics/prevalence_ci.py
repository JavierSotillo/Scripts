"""
Prevalence and confidence interval (normal/Wald approximation).

For x positives out of n individuals examined, prevalence is p = x/n and the
confidence interval is p +/- z * sqrt(p(1 - p)/n), with z = 1.96 for 95%.
The Wald interval is adequate for moderate sample sizes, but performs poorly
for very small samples or prevalences close to 0% or 100%; in those cases a
Wilson or exact (Clopper-Pearson) interval is preferable.

Usage:
    python prevalence_ci.py 48 208          # positives, total examined
    python prevalence_ci.py 48 208 2.576    # optional z-value (here 99% CI)

Requires Python 3 only (no additional packages). It can also be run in any
online Python interpreter by editing the values at the bottom of the file.
"""

import math
import sys


def prevalence_ci(positives, total, z=1.96):
    """Return prevalence and Wald CI as percentages: (prevalence, lower, upper)."""
    if total <= 0 or not 0 <= positives <= total:
        raise ValueError("positives must be between 0 and total, and total > 0")
    p = positives / total
    se = math.sqrt(p * (1 - p) / total)
    return 100 * p, 100 * (p - z * se), 100 * (p + z * se)


if __name__ == "__main__":
    if len(sys.argv) in (3, 4):
        positives, total = int(sys.argv[1]), int(sys.argv[2])
        z = float(sys.argv[3]) if len(sys.argv) == 4 else 1.96
    else:
        # Edit these values when running without arguments
        positives, total, z = 48, 208, 1.96

    prev, low, high = prevalence_ci(positives, total, z)
    print(f"Positives: {positives}/{total}")
    print(f"Prevalence: {prev:.2f}%")
    print(f"Confidence interval (z = {z}): {low:.2f}% - {high:.2f}%")
