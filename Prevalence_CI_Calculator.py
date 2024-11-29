#Go to https://www.online-python.com/

import math

# Given data
positives = 48
total_population = 92
z = 1.96  # Z-value for 95% confidence

# Prevalence
prevalence = positives / total_population

# Standard Error (SE)
se = math.sqrt(prevalence * (1 - prevalence) / total_population)

# Confidence Interval (CI)
ci_lower = prevalence - z * se
ci_upper = prevalence + z * se

# Convert to percentage
prevalence_percentage = prevalence * 100
ci_lower_percentage = ci_lower * 100
ci_upper_percentage = ci_upper * 100

prevalence_percentage, ci_lower_percentage, ci_upper_percentage
# Print the results
print("Prevalence: {:.2f}%".format(prevalence_percentage))
print("95% Confidence Interval: ({:.2f}%, {:.2f}%)".format(ci_lower_percentage, ci_upper_percentage))