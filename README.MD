## Implementation of popular multiple testing and gatekeeping procedures in SAS 

A collection of SAS macros to support a wide variety of popular multiple testing and gatekeeping procedures, including the computation of multiplicity adjusted p-values and confidence intervals: 

* PvalProc macro computes adjusted p-values for popular nonparametric procedures (Bonferroni, Holm, fixed-sequence, fallback, Hommel, and Hochberg procedures) in multiplicity problems with equally or unequally weighted null hypotheses.

* ParProc macro computes adjusted p-values for parametric procedures from the Dunnett family (single-step and step-down Dunnett procedures) in one-sided hypothesis testing problems with a balanced one-way layout and equally weighted null hypotheses. 

* PvalCI macro computes simultaneous confidence intervals for popular nonparametric procedures (Bonferroni, Holm, fixed-sequence, and fallback procedures) in multiplicity problems with equally or unequally weighted null hypotheses.

* ParCI macro computes simultaneous confidence intervals for parametric procedures from the Dunnett family (single-step and step-down Dunnett procedures) in one-sided hypothesis testing problems with a balanced one-way layout and equally weighted null hypotheses.

* TreeGate macro computes adjusted p-values and supports simulation-based power calculations for tree-gatekeeping procedures based on Bonferroni, Holm, Hochberg and Hommel procedures. A detailed user manual is provided.

Examples of SAS programs for computing adjusted p-values and confidence intervals can be found in Dmitrienko and D'Agostino (2013, Section 11).

## Code for online training courses

The R and SAS code used in the online training courses:

* Key Multiplicity Issues in Clinical Trials (Part I).
* Key Multiplicity Issues in Clinical Trials (Part II).

The courses are available on [Mediana's YouTube channel](https://www.youtube.com/channel/UCAQUAsQxxnzje4rEdJ4870Q). 

## References

* Alosh, M., Bretz, F., Huque, M. Advanced multiplicity adjustment methods in clinical trials. Statistics in Medicine. 33:693-713, 2014.

* Dmitrienko, A., Tamhane, A.C., Bretz, F. (editors). Multiple Testing Problems in Pharmaceutical Statistics. Chapman and Hall/CRC Press, New York, 2009.

* Dmitrienko, A., D'Agostino, R.B., Huque, M.F. Key multiplicity issues in clinical drug development. Statistics in Medicine. 32:1079-1111, 2013.

* Dmitrienko, A., D'Agostino, R.B. Traditional multiplicity adjustment methods in clinical trials. Statistics in Medicine. 32:5172-5218, 2013.


