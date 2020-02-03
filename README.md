# diagnostic-accuracy-meta-analysis
R syntax files for performing meta-analysis of diagnostic accuracy studies

1. bivariate MA diagnostic accuracy.R
    Syntax to perform a random-effects meta-analysis of a series of tests of a single diagnostic test. Input data is a data.frame with true and false positives and negatives.
    Main output is an ROC-curve with AUC-value. Analysis is similar to what the mada-package offers. Difference is that this syntax uses the binomial-normal model, uses a Bayesian framework and offers 95% credibility intervals.
    
2. MA diagnostic accuracies 2 tests FE model 1a.R
    MA diagnostic accuracies 2 tests RE model 1b.R
    MA diagnostic accuracies 2 tests FE model 2a.R
    MA diagnostic accuracies 2 tests RE model 2b.R
    Syntaxes to perform fixed-effects (FE) and random-effects (RE) meta-analysis of a series of tests comparing two diagnostic tests. Input data is a data.frame with true and false positives and negatives for both tests. 
    Version "b" also handles studies that report results of both tests in 2-by-2 tables of both tests for cases and controls separately. Main output is a pair of ROC-curves, AUC-values and their difference. The syntax uses binomial-normal models, uses a  Bayesian Framework and offers 95% credibility intervals.
    
3. netwerk MADAS.R
    As above for multiple diagnostic tests.

4. 