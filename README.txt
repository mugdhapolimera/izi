IZI - A Bayesian Metallicity Esitmator
Python Version

Improvements:

1. Change non-informative uniform log(Z) and log(q) prior to informative prior based on mass metallicity relation.

2. Change the pipeline to parallelly process multiple galaxies instead of one at a time. 
OR
Change the full Bayesian calculation to MCMC. 

3. Use better grid with finer spacing, more lines. Also, change the original grid by removing/adding one line to avoid the 35x35 grid.

4. Error in the model is estimated as e_mod = epsilon * f_mod. (f_mod is the Z or q value from the model; epsilon is a constant). 
epsilon can be changed to a nuisance parameter for more physical errors.

5. Make an automatic selection between using scipy method or gpy method for interpolation of grid. 

6. Automatically accept or reject lines depending on the model.

7. Make d.flags a dict