# Numerical tests

Use a fresh `stats::lm.fit()` on the rows of each complete window as the
reference for rolling fits. When testing nearly collinear predictors, verify
full column rank in each reference fit. Check the first window separately
from subsequent windows to distinguish initialization errors from update or
downdate errors.
