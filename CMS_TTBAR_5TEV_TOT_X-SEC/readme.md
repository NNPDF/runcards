These runcards were used to compute top pair production processes for NNPDF4.1 with MATRIX at NNLO except for the 
following changes:
* The top mass was varied up (mt_+ = 175.0 GeV) and down (mt_- = 170.0 GeV) with respect to the central value mt_0 = 172.5 GeV. 
* The renormalisation and factorisation scales were dynamical but evaluated at the central value mt_0. That is, 
  the scales were set to `HT(mt=mt_0, pT)/4`, with `HT = (ET_t + ET_tx)` and `ET_t = sqrt(mt_0 ** 2 + pT_t ** 2)`, 
  `ET_tx = sqrt(mt_0 ** 2 + pT_tx ** 2)`. This choice is reflected in the `paramter.dat` card,
```
: dynamic_scale     =  11
```
where `dynamic_scale = 11` corresponds to the definition written in `specify.scales.cxx`

To reproduce the NNPDF4.1 predictions, these runcards should be modified by setting in `model.dat`
```
 6 172.5 # M_t
```

