# This file stays at this path because it is cited in:
#   Mutengo MM, et al. (2026) PLoS Negl Trop Dis 20(10): e0014013.
#   https://doi.org/10.1371/journal.pntd.0014013
# The same script, with the article's documentation, is in
# papers/Mutengo2026_Shaematobium-Lusaka/
#
"""
Prevalence and 95% confidence interval (normal/Wald approximation).

Used in:
    Mutengo MM, et al. (2026) Schistosoma haematobium infection among school
    age children in two peri-urban areas of Lusaka, Zambia: A cross-sectional
    study. PLoS Negl Trop Dis 20(10): e0014013.
    https://doi.org/10.1371/journal.pntd.0014013

Usage:
    python prevalence_ci.py 48 208     # positives, total examined
    python prevalence_ci.py            # reproduces Table 2 of the article

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


def report(positives, total, label=""):
    prev, low, high = prevalence_ci(positives, total)
    prefix = f"{label:<20}" if label else ""
    print(f"{prefix}{positives:>4}/{total:<4} {prev:6.2f}% (95% CI: {low:.2f}-{high:.2f})")


# Table 2: positives and children examined per district and diagnostic method
TABLE_2 = [
    ("Microscopy", "Chongwe", 34, 111),
    ("Microscopy", "Kafue", 14, 97),
    ("Microscopy", "Overall", 48, 208),
    ("qPCR", "Chongwe", 82, 111),
    ("qPCR", "Kafue", 27, 97),
    ("qPCR", "Overall", 109, 208),
    ("ELISA", "Chongwe", 51, 111),
    ("ELISA", "Kafue", 49, 97),
    ("ELISA", "Overall", 100, 208),
]

if __name__ == "__main__":
    if len(sys.argv) == 3:
        report(int(sys.argv[1]), int(sys.argv[2]))
    else:
        print("Table 2 - prevalence of S. haematobium infection by district\n")
        for method, district, positives, total in TABLE_2:
            report(positives, total, f"{method} {district}" if district != "Overall" else f"{method} overall")
