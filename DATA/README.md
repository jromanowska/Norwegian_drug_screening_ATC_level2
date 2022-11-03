## DATA folder

### `expose_2_prescr_PD-4prescr-levo-mao-b_age-time-scale_results_atc2level_all.csv`

All the results from the pooled analysis.

 **column name**  |  **description**  |  **format**  |  **notes**
:----------------:|:------------------|:-------------|:----------------
 ATC_code         | level 2 ATC code  | character    | one letter and two digits
 pvals            | p-value of the Cox-regression HR | numeric | 
 name             | description of the ATC code | character string | 
 group            | ATC group         | character    | the first letter of the ATC code
 HR               | hazard risk       | numeric      | 
 HR_l             | lower bound of 95% confidence interval of the HR value | numeric | 
 HR_u             | higher bound of 95% confidence interval of the HR value | numeric | 
 p.adj.FDR        | false discovery rate-adjusted p-value | numeric | 
 N.nonpd          | number of non-PD users of the drugs | integer | 
 N.pd             | number of users that developed PD | integer | 

Similarly, the files `expose_2_prescr_PD-4prescr-levo-mao-b_age-time-scale_results_atc2level_men.csv`
and `expose_2_prescr_PD-4prescr-levo-mao-b_age-time-scale_results_atc2level_women.csv`
contain the same type of data but for the sex-stratified analyses.

### `signif_res_pooled_all.txt`

Significant results (i.e., FDR < 0.05) from the pooled analysis.

 **column name**  |  **description**  |  **format**  |  **notes**
:----------------:|:------------------|:-------------|:----------------
 ATC_code         | level 2 ATC code  | character    | one letter and two digits
 pvals            | p-value of the Cox-regression HR | numeric | 
 name             | description of the ATC code | character string | 
 group            | ATC group         | character    | the first letter of the ATC code
 HR               | hazard risk       | numeric      | 
 HR_l             | lower bound of 95% confidence interval of the HR value | numeric | 
 HR_u             | higher bound of 95% confidence interval of the HR value | numeric | 
 p.adj.FDR        | false discovery rate-adjusted p-value | numeric | 
 N.nonpd          | number of non-PD users of the drugs | integer | 
 N.pd             | number of users that developed PD | integer | 
 direction_change | which direction is the risk change? | character | two levels: 'increase' or 'decrease'

### `signif_res_compare_pool_strat.txt`

Comparing the significant results from the pooled analysis with the results from
the stratified analyses.

