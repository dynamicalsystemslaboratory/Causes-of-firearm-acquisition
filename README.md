
#    codes

guns_ent2.m
- matlab code to compute conditional entropies of binary multivariate time-series.

guns_boot.m
- matlab code to generate surrogate binary time-series following a local permutation scheme.

Study1.m
- matlab code to perform study 1 and reproduce data in Table 2. This code computes conditional transfer entropy on binary time-series, generates surrogate binary time-series following a local permutation scheme, and estimates p-values for conditional transfer entropy.

Study2.m
- matlab code to perform study 2 and reproduce data in Table 3. This code computes conditional transfer entropy on binary time-series, generates surrogate binary time-series following a local permutation scheme, and estimates p-values for conditional transfer entropy.

Study3.m
- matlab code to perform study 3 and reproduce data in Table 4. This code computes conditional transfer entropy on binary time-series, generates surrogate binary time-series following a local permutation scheme, and estimates p-values for conditional transfer entropy.


#    Excel data  files 



Mass shootings_1999_2017_state_level.xlsx:
- Monthly mass shootings for each state from 1999 to 2017 (source: from Washington Post).


Media Coverage 1999-2017.xlsx: 
- Monthly media coverage on firearm control and shootings from 1999 to 2017 for Chicago Tribune, Los Angeles Times, New York Times, Times Picayune, and Wall Street Journal.


#    matlab data files   


BChecks.mat
- Monthly background checks for each state from 1999 to 2017, detrended and seasonally adjusted for each State except for Connecticut and Hawaii.

D.mat
-  Distance between States was measured from the geodesics between their centroids (source: N. Abaid et al. "The effect of geography and citizen behavior on motorvehicle deaths in the United States,"  PLoS ONE 10 (2015) e0123339).

laws.mat
- 133 firearm-related laws in each State from 1991 to 2017 (source: Firearm Laws Project).


MOs_Patterns.mat
-  Monthly media coverage on shootings from 1999 to 2017 for New York Times, Times Picayune, Chicago Tribune, Los Angeles Times, and Wall Street Journal.

MOfc_Patterns.mat
- Monthly media coverage on firearm control and shootings from 1999 to 2017 for  New York Times, Times Picayune, Chicago Tribune, Los Angeles Times, and Wall Street Journal.

Mshooting_Patterns_states.mat
- Binary time-series at a state level on the incidence of mass shootings obtained from the Washington Post database.

MStotal_Patterns.mat
- Binary time-series ata a ntional level on the incidence of mass shootings obtained from the Washington Post database.

permn.m
- matlab code used to compute all possible values of a binary multivariate time-series.

population.mat
- Average populations by state from 1999 to 2017 (source: US Census Bureau).
