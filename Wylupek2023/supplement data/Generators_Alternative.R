### Hereunder, we provide the generators for the alternatives defined in Examples 1-9.
### Copying a selected generator and pasting it to the file ''Power_functions.R''
### in a proper place allows one to calculate the powers of the tests investigated 
### in the paper under the given alternative. A selection of the parameter
### which defines the alternative, the sample size n, the parameter $c$,
### and the critical value of the $VGH$ test (cv.VGH in the file ''Power_functions.R'')
### enables one to calculate the powers under different configurations of the parameters. 

##################################
### Example 1, Generator begins 
##################################

mu_X = 0
Sigma_X = 1

mu_Y = 0.8		### this parameter changes and it occurs in the figures as $mu$
Sigma_Y = 1

rho = -0.8/(Sigma_X * Sigma_Y)

Z1 = rnorm(n)
Z2 = rnorm(n)

X = Sigma_X * Z1 + mu_X
Y = Sigma_Y * (rho * Z1 + sqrt(1-rho^2) * Z2) + mu_Y

#################################
### Example 1, Generator ends 
#################################

##################################
### Example 2, Generator begins 
##################################

mu_X = 0
Sigma_X = 1

mu_Y = 0
Sigma_Y = 2		### this parameter changes and it occurs in the figures as $sigma$

rho = 0.95/(Sigma_X * Sigma_Y)

Z1 = rnorm(n)
Z2 = rnorm(n)

X = Sigma_X * Z1 + mu_X
Y = Sigma_Y * (rho * Z1 + sqrt(1-rho^2) * Z2) + mu_Y

#################################
### Example 2, Generator ends 
#################################

##################################
### Example 3, Generator begins 
##################################

mu_X1 = 1
Sigma_X1 = 1

mu_Y1 = 1
Sigma_Y1 = 1

rho1 = -0.8/(Sigma_X1 * Sigma_Y1)

Z1 = rnorm(n)
Z2 = rnorm(n)

X1 = Sigma_X1 * Z1 + mu_X1
Y1 = Sigma_Y1 * (rho1 * Z1 + sqrt(1-rho1^2) * Z2) + mu_Y1

mu_X2 = 2
Sigma_X2 = 2

mu_Y2 = 2
Sigma_Y2 = 12	### this parameter changes and it occurs in the figures as $sigma$

rho2 = -3.6/(Sigma_X2 * Sigma_Y2)

Z3 = rnorm(n)
Z4 = rnorm(n)

X2 = Sigma_X2 * Z3 + mu_X2
Y2 = Sigma_Y2 * (rho2 * Z3 + sqrt(1-rho2^2) * Z4) + mu_Y2

U = runif(n)

X = ifelse(U < 0.6, X1, X2)
Y = ifelse(U < 0.6, Y1, Y2)

#################################
### Example 3, Generator ends 
#################################

##################################
### Example 4, Generator begins 
##################################

Beta = 0.7		### this parameter changes and it occurs in the figures as $beta$

	X = rnorm(n)
	E = rnorm(n)
	
	Y = Beta*X^2 + E

#################################
### Example 4, Generator ends 
#################################

##################################
### Example 5, Generator begins 
##################################

Beta = 1.4		### this parameter changes and it occurs in the figures as $beta$

	X = rnorm(n)
	E = rnorm(n)
	
	Y = -(Beta*X^3 + E)

#################################
### Example 5, Generator ends 
#################################

##################################
### Example 6, Generator begins 
##################################

Beta = 0.6		### this parameter changes and it occurs in the figures as $beta$

	X = rnorm(n)
	E = rnorm(n)
	
	Y = Beta*X^4 + E

#################################
### Example 6, Generator ends 
#################################

##################################
### Example 7, Generator begins 
##################################

library(pscl)

theta = 0.14	### this parameter changes and it occurs in the figures as $theta$

par = theta

nu = 5
rho = 0.8

gamma1 =  par
gamma2 =  -par

	Z1 = rnorm(n)
	Z2 = rho*Z1 + sqrt(1-rho^2)*rnorm(n)

	W = rigamma(n,nu/2,nu/2)

	X = gamma1*W + sqrt(W)*Z1
	Y = gamma2*W + sqrt(W)*Z2

#################################
### Example 7, Generator ends 
#################################

##################################
### Example 8, Generator begins 
##################################

a2 = 0.4	### this parameter changes and it occurs in the figures as $a$

b2 = 0
a3 = a2	
b3 = b2
Delta = 0.15		

Y1 = rnorm(n,0,sqrt(2))

	W2 = rexp(n)
	F2 = runif(n, -pi/2, pi/2)
	F20 = -pi*b2*(1-abs(1-a2))/(2*a2)

Y2 = ( sin(a2*(F2-F20))/((cos(F2))^(1/a2)) )*( (cos( F2 - a2*(F2-F20) ))/ W2 )^((1-a2)/a2)

	W3 = rexp(n)
	F3 = runif(n, -pi/2, pi/2)
	F30 = -pi*b3*(1-abs(1-a3))/(2*a3)

Y3 = ( sin(a3*(F3-F30))/((cos(F3))^(1/a3)) )*( (cos( F3 - a3*(F3-F30) ))/ W3 )^((1-a3)/a3)

	X =  Y1 + Delta*Y3
	Y =  Y2 + Delta*Y3 

#################################
### Example 8, Generator ends 
#################################

##################################
### Example 9, Generator begins 
##################################

b2 = 0.8	### this parameter changes and it occurs in the figures as $a$

a1 = 0.4	
b1 = 0.1
a2 = b2
a3 = a2	
b3 = b2

Delta = 0.15		

	W1 = rexp(n)
	F1 = runif(n, -pi/2, pi/2)
	F10 = -pi*b1*(1-abs(1-a1))/(2*a1)

Y1 = ( sin(a1*(F1-F10))/((cos(F1))^(1/a1)) )*( (cos( F1 - a1*(F1-F10) ))/ W1 )^((1-a1)/a1)

	W2 = rexp(n)
	F2 = runif(n, -pi/2, pi/2)
	F20 = -pi*b2*(1-abs(1-a2))/(2*a2)

Y2 = ( sin(a2*(F2-F20))/((cos(F2))^(1/a2)) )*( (cos( F2 - a2*(F2-F20) ))/ W2 )^((1-a2)/a2)

	W3 = rexp(n)
	F3 = runif(n, -pi/2, pi/2)
	F30 = -pi*b3*(1-abs(1-a3))/(2*a3)

Y3 = ( sin(a3*(F3-F30))/((cos(F3))^(1/a3)) )*( (cos( F3 - a3*(F3-F30) ))/ W3 )^((1-a3)/a3)

	X =  Y1 + Delta*Y3
	Y =  Y2 + Delta*Y3 

#################################
### Example 9, Generator ends 
#################################
