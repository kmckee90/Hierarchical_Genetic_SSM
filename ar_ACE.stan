// Author: Kevin L. McKee (klmckee@vt.edu)
// Last updated: 3-7-2021
// Mixed-effects autoregression with ACE decomposition of AR(1) parameters.
// Assumes equal number of MZ and DZs
// Includes standardized ACE variance components in the output.

data {
  int<lower=0> N;              // num individuals
  int<lower=0> T;              // num individuals
  matrix[T, N] MZ1;              // outcomes  
  matrix[T, N] MZ2;              // outcomes  
  matrix[T, N] DZ1;              // outcomes  
  matrix[T, N] DZ2;              // outcomes  

}
parameters {
  vector[2] betaMZ_z[N];
  vector[2] betaDZ_z[N];
  
  real betaMu;
  real <lower=0> q;
  
  real<lower=0> a;
  real<lower=0> c;
  real<lower=0> e;
}
transformed parameters{
  vector[2] betaMZ[N];
  vector[2] betaDZ[N];
  for(i in 1:N){
	betaMZ[i] = betaMZ_z[i]+rep_vector(betaMu,2);
	betaDZ[i] = betaDZ_z[i]+rep_vector(betaMu,2);
  }
}
model {
	matrix[2,2] SigMZ;
	matrix[2,2] SigDZ;
	vector[T*N*4] err;
	real nu;
	int u = 0;

	SigMZ[1,1] = a+c+e;
	SigMZ[2,2] = a+c+e;
	SigMZ[1,2] = a+c;
	SigMZ[2,1] = a+c;
		
	SigDZ[1,1] = a+c+e;
	SigDZ[2,2] = a+c+e;
	SigDZ[1,2] = .5*a+c;
	SigDZ[2,1] = .5*a+c;
		
	for(i in 1:N){ 
				nu = 0;
				u+=1;
				err[u] = MZ1[1,i] - nu;
				for (t in 2:T) {
					u+=1;
					nu = (betaMZ[i,1]) * MZ1[t-1,i];
					err[u] = MZ1[t,i] - nu;
				}
				u+=1;
				nu = 0;
				err[u] = MZ2[1,i] - nu;
				for (t in 2:T) {
					u+=1;
					nu = (betaMZ[i,2]) * MZ2[t-1,i];
					err[u] = MZ2[t,i] - nu;
				}
				u+=1;
				nu = 0;
				err[u] = DZ1[1,i] - nu;
				for (t in 2:T) {
					u+=1;
					nu = (betaDZ[i,1]) * DZ1[t-1,i];
					err[u] = DZ1[t,i] - nu;
				}
				u+=1;
				nu = 0;
				err[u] = DZ2[1,i] - nu;
				for (t in 2:T) {
					u+=1;
					nu = (betaDZ[i,2]) * DZ2[t-1,i];
					err[u] = DZ2[t,i] - nu;
				}
			}
	
	err ~ normal(0, q);	
	betaMZ_z ~ multi_normal(rep_vector(0,2), SigMZ);
	betaDZ_z ~ multi_normal(rep_vector(0,2), SigDZ);
	betaMu ~ normal(0,100);
	
	q ~ normal(0, 100);
	a ~ normal(0, 100);
	c ~ normal(0, 100);
	e ~ normal(0, 100);
}
generated quantities {
  real v = a+c+e;
  real A = a/v;
  real C = c/v;
  real E = e/v;
}

