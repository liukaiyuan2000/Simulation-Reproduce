### Hereunder, we provide the generators for the models obeying the null hypohtesis.
### Copying a selected generator and pasting it to the file ''Power_functions.R''
### in a proper place allows one to calculate the Type I errors of the tests 
### investigated in the paper under the given model. A selection of the parameter
### which defines the hypothesis, the sample size n, the parameter $c$,
### and the critical value of the $VGH$ test (cv.VGH in the file ''Power_functions.R'')
### enables one to calculate the Type I errors under different configurations of the parameters. 

###################################
### FGM copula, Generator begins 
###################################

theta = 0.5		### this parameter changes and it occurs in the tables as $theta$

	X = runif(n)
	T = runif(n)

	A = 1 + theta*(1-2*X)
	B = sqrt( A^2 - 4*(A-1)*T )
	Y = 2*T/(A+B)

##################################
### FGM copula, Generator ends 
##################################

#######################################
### Clayton copula, Generator begins 
#######################################

theta = 1		### this parameter changes and it occurs in the tables as $theta$

phi.f = function(t){
	(t^(-theta) - 1)/theta
}

	U = runif(n)
	T = runif(n)

W= exp( log( U^(-theta-1)/T )/(-theta-1) )
V= ifelse( phi.f(W)-phi.f(U) <= phi.f(0), exp( -log( theta*(phi.f(W)-phi.f(U))+1 )/theta ), 0)

	X = U
	Y = V

#####################################
### Clayton copula, Generator ends 
#####################################

#######################################
### Mardia copula, Generator begins 
#######################################

theta = 0.9		### this parameter changes and it occurs in the tables as $theta$

	P1 = runif(n)
	P2 = runif(n)
	
	U.M = runif(n)	
	
	M1 = U.M
	M2 = U.M

	U.W = runif(n)
	
	W1 = U.W
	W2 = 1-U.W

	U = runif(n)

	X = ifelse(U < theta^2*(1-theta)/2, W1, ifelse(U < theta^2*(1-theta)/2 + 1-theta^2, P1, M1) )
	Y = ifelse(U < theta^2*(1-theta)/2, W2, ifelse(U < theta^2*(1-theta)/2 + 1-theta^2, P2, M2) )

#####################################
### Mardia copula, Generator ends 
#####################################

#######################################
### t-distribution, Generator begins 
#######################################

nu = 5		### this parameter changes and it occurs in the tables as $theta$

	Z1 = rnorm(n)
	Z2 = rnorm(n)

	S = rchisq(n,nu)

	X = Z1/(sqrt(S/nu))
	Y = Z2/(sqrt(S/nu))

#####################################
### t-distribution, Generator ends 
#####################################
