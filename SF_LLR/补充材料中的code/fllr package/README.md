# fllr: Functional local linear regression #
## Frédéric Ferraty and Stanislav Nagy ##

R package fllr with an implementation of the **functional local linear regression methods** in the scalar-on-function regression model. Allows for local linear estimation of both the regression operator, and the functional derivative. 
To install the package *fllr* in Windows you need to have R Tools installed on your computer, see 

https://cran.r-project.org/bin/windows/Rtools/

Then, run the following code


```
#!R

install.packages("https://bitbucket.org/StanislavNagy/fllr/raw/ * key * /fllr_0.0.1.tar.gz",type="source")
library(fllr)
help(package=fllr)
help(fllr)
```


The *key* in the link above (after raw/) must be modified to link to the newest version of the fllr tar.gz file. Alternatively, download the latest version of the tar.gz file to the R working directory in your computer (see getwd() in R) and run

```
#!R


install.packages("fllr_0.0.1.tar.gz", repos = NULL, type="source")
```


In case of questions contact 

***Stanislav Nagy***

*nagy@karlin.mff.cuni.cz*