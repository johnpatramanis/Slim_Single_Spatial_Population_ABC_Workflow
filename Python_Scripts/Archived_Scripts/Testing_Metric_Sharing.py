import numpy as np


a_n = 0.20
b_n = 0.20

a_h = 1 - a_n
b_h = 1 - b_n



### For Neanderthal Sharing
observed_n = 0.05
expected_n = a_n * b_n
normalize_n = np.sqrt( a_n - a_n**2) * np.sqrt( b_n - b_n**2)
Metric_n = (observed_n - expected_n) / normalize_n


### For Human Sharing

observed_h = 1 - observed_n
expected_h = a_h * b_h
normalize_h = np.sqrt( a_h - a_h**2) * np.sqrt( b_h - b_h**2)
Metric_h = (observed_h - expected_h) / normalize_h


print(Metric_n, Metric_h)