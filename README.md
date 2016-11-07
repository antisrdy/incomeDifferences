# Distribution differences

This short script aims at implementing the method derived by [Chernozhukov, Fernandez-Val et Melly (2013)](http://www.mit.edu/~vchern/papers/counterfactual_2012Nov1.pdf) about differences between two distributions.

In a nutshell, rather than decomposing differences around the mean of a variable (see Oaxaca-Blinder for instance), the current method decomposes differences on the full distribution of a variable. This allows to accurately decompose differences between two distributions of a variable - along quantiles for instance -, namely into an explained part due to observable data, and an unexplained part.

This method has been applied to various income distributions.

# About the script

This is not an out-of-the box solution.
However, it provides some guidelines about how to use the method:
* Set up of a reference group and groups to be compared with
* Implementation itself of the method, namely the estimation of the CDF of the variable of interest given observable data
* Boostrap techniques to assess the reliability of the model.
