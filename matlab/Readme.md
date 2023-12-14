"Entropic Wasserstein? That's esoteric crystal healing!"

This folder contains some implementations of optimal transport related algorithms, including

* (approximate) solvers for linear assignment problems (#sources = #sinks):
  **habr.m, hungarian.m, lapmack.m, lapjv.m, ungarisch.m**
* linear programmes for general transport problems:
  **transsimp.m, wasseropt.m** 
* Sinkhorn algorithms from entropic Wasserstein distances:
  **sinkhorn.m, sinkfast.m**  
* OT-based sensitivity measures for 1D and multidimensional outputs:
  **wassersi.m, sinksens.m**
* Bures-Wasserstein sensitivity measure: **bwsi.m**
* some other sensitivity measures to compare with
* and some utility functions
  
Please see bigwassersteintest.m for the all-purpose test routine, sinksens.m for a stripped down implementation, and examplerun.m for a demo.
