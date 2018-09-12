# rollRegres 0.1.1
* handle data sets with gaps
* add the option to require a minimum number of observations in a window
* output with groups have changed such that the first window max `grp` value will have 
at least `width - 1L` distance from `min(grp)` in `roll_regres` and 
`roll_regres.fit`
