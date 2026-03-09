# DisparityBias

### Disparity project

Goals of the project - exploration of the effect of incomplete sampling (spatial vs non-spatial, temporal) on measures of disparity
- exploration of that combined with effect of extinction selectivity

Tree simulation - B,D rates fixed & homogeneous
* constant rates
* ME w/optional selectivity (no, trait-based, geo-based, correlated trait-geo)
* radiation
-> test magnitude of the ME (low/high) + only keep one for the correlated sims
-> magnitude of the correlation (low/high) & magnitude of radiation (low/high)

Trait simulation
=> simulation of traits on the complete tree - low/high rate of evolution
=> keep 2 traits but test one dataset with varying rate of evolution for each trait
=> decide how to simulate geographically-correlated traits -> OU w/different optima ? BM with different speeds ?

Geographical simulation
* joined trees so defined by the tree simulation
* 1 vs 5 migration events
=> to see how this can be adapted for the ME setup - need to sample the attachment points before the ME

Fossil sampling
- sampling variations (constant, temporal, spatial, spatial in 2nd bin only)
=> test high and low sampling rates
=> trait value based on the position of the fossil specimen for each fossil 
=> control dataset is high constant fossil sampling

Output metrics
look at increase/decrease in disparity - relative change between bins (with sign of change)
+ look at change compared to the benchmark (almost complete sampling)
