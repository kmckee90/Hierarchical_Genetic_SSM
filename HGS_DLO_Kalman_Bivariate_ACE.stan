//Author: Kevin McKee (klmckee@vt.edu)
//Kalman-Bucy Filter based continuous-time State Space Model for Stan
//With Multivariate Hierarchical Gene-Environment parameter covariance decomposition.

functions{
//Matrix discretization function for continuous time state-space models.
	matrix[] discretize(matrix A, matrix B, matrix Q, real dt) {
		int m = rows(A);
		int n = m;
		matrix[m,n] expA;
		matrix[m,n] intA; 
		matrix[m,n] I = diag_matrix(rep_vector(1.0, m));
		matrix[m,n] output[3];
	// First Block expm for A integral, and expm(A*deltaT)
		matrix[2*m,2*n] M = append_col(rep_matrix(0, 2*m, n), append_row(I, A));
		M = matrix_exp(M*dt);
		intA = block(M, 1,  n+1, m, n);
		expA = block(M, m+1,  n+1, m, n);
	// Second Block expm for discretized Q
		M = append_col( append_row(-1*A', rep_matrix(0, m, n)), append_row(Q, A));
		M = matrix_exp(M*dt);	
		output[1] = expA;
		output[2] = intA*B;
		output[3] = expA' * block(M, 1,  n+1, m, n);
		return output;
	}	
}
data {
  int<lower=0> nMZ;               // num individuals
  int<lower=0> nDZ;               // num individuals
  int<lower=1> p;				// num indicators
  int<lower=1> q;				// num states
  int<lower=0> T;              	// num occasions
  real dt;						// sampling interval
  matrix[T, p] MZ[nMZ,2];              // outcomes  
  matrix[T, p] DZ[nDZ,2];              // outcomes
}
parameters {
  matrix[nMZ, 2] muMZ;
  matrix[nDZ, 2] muDZ;
  vector<lower=0, upper=pi()/(2*dt)>[2] freqMZ[nMZ];
  vector<lower=0, upper=pi()/(2*dt)>[2] freqDZ[nDZ];
  real<lower=0, upper=pi()/(2*dt)> freqMu;
  vector<lower=0, upper=1>[2] dampingRatioMZ[nMZ];
  vector<lower=0, upper=1>[2] dampingRatioDZ[nDZ]; 
  real dampingRatioMu;
  matrix<lower=0>[nMZ, 2] sigmaMZ_r;
  matrix<lower=0>[nDZ, 2] sigmaDZ_r;
  vector<lower=0>[p] epsilon_r;
  vector[p-1] loadings;
  real<lower=0> a_eta;
  real<lower=0> c_eta;
  real<lower=0> e_eta; 
  real<lower=0> a_zeta;
  real<lower=0> c_zeta;
  real<lower=0> e_zeta;
  real rSE;
  real ra;
  real rc;
  vector[2] x0;
  cov_matrix[2] P0;
}
transformed parameters {
	vector[2] etaMZ[nMZ];
	vector[2] etaDZ[nDZ];
	vector[2] zetaMZ[nMZ];
	vector[2] zetaDZ[nDZ];
	matrix[nMZ, 2] sigmaMZ = sigmaMZ_r .* sigmaMZ_r;
	matrix[nDZ, 2] sigmaDZ = sigmaDZ_r .* sigmaDZ_r;
    vector[p] epsilon = epsilon_r .* epsilon_r;
	for(i in 1:nMZ){
		etaMZ[i] = -(freqMZ[i] .* freqMZ[i]);
		zetaMZ[i]= -(dampingRatioMZ[i].*freqMZ[i]*2.0);
	}
	for(i in 1:nDZ){
		etaDZ[i] = -(freqDZ[i] .* freqDZ[i]);
     	zetaDZ[i]= -(dampingRatioDZ[i].*freqDZ[i]*2.0);
	}
}
model {
  	matrix[4,4] SigMZ = rep_matrix(0.0,4,4);
	matrix[4,4] SigDZ = rep_matrix(0.0,4,4);
	matrix[6,6] SigA_MZ = diag_matrix(rep_vector(0.5,6));
	matrix[6,6] SigA_DZ = diag_matrix(rep_vector(0.5,6));
	matrix[6,6] SigC;
	matrix[6,6] SigE = diag_matrix(rep_vector(1.0,6));
	matrix[4,6] LA = rep_matrix(0.0, 4, 6);
	matrix[4,6] LC = rep_matrix(0.0, 4, 6);
	matrix[4,6] LE = rep_matrix(0.0, 4, 6);
	vector[4] out_pars_MZ[nMZ];
	vector[4] out_pars_DZ[nDZ];		
	vector[4] out_Mu;
	matrix[p,q] H = rep_matrix(0.0,p,q);
	matrix[q,q] A = rep_matrix(0.0,q,q);
	matrix[q,q] B = rep_matrix(0.0,q,q);
	matrix[q,q] Q = rep_matrix(0.0,q,q);
	matrix[q,q] discMats[3];
	matrix[q,q] Ad;
	matrix[q,q] Bd;
	matrix[q,q] Qd;
	matrix[q,q] P;
	matrix[q,q] Pp;
	matrix[p,p] S;
	matrix[p,p] R;
	matrix[q,p] K;
	vector[q] x;
	vector[q] xp;
	vector[p] z;
	vector[p] r;
	vector[T] err;
	R = diag_matrix( epsilon );
	H[1,1] = 1.0;
	if(p>1){
		H[2:p,1] = loadings;
	}
//Genetic Model
	SigA_MZ[4,1] = 1.0;
	SigA_MZ[5,2] = 1.0;
	SigA_MZ[6,3] = 1.0;
	SigA_DZ[4,1] = 0.5;
	SigA_DZ[5,2] = 0.5;
	SigA_DZ[6,3] = 0.5;
	SigA_MZ = SigA_MZ + SigA_MZ';
	SigA_DZ = SigA_DZ + SigA_DZ';
	SigC = SigA_MZ;
	LA[1,1] = a_eta;
	LA[1:2,2] = rep_vector(ra,2);
	LA[2,3] = a_zeta;
	LC[1,1] = c_eta;
	LC[1:2,2] = rep_vector(rc,2);
	LC[2,3] = c_zeta;
	LE[1,1] = e_eta;
	LE[1:2,2] = rep_vector(rSE,2);
	LE[2,3] = e_zeta;
	LA[3,4] = a_eta;
	LA[3:4,5] = rep_vector(ra,2);
	LA[4,6] = a_zeta;
	LC[3,4] = c_eta;
	LC[3:4,5] = rep_vector(rc,2);
	LC[4,6] = c_zeta;
	LE[3,4] = e_eta;
	LE[3:4,5] = rep_vector(rSE,2);
	LE[4,6] = e_zeta;
	SigMZ = LA*SigA_MZ*LA' + LC*SigC*LC' + LE*SigE*LE';
	SigDZ = LA*SigA_DZ*LA' + LC*SigC*LC' + LE*SigE*LE';
//Run Kalman-Bucy filter
	out_Mu[1] = freqMu;
	out_Mu[2] = dampingRatioMu;
	out_Mu[3] = freqMu;
	out_Mu[4] = dampingRatioMu;
//MZ twins	
	for(i in 1:nMZ){
		out_pars_MZ[i,1] = freqMZ[i,1];
		out_pars_MZ[i,2] = dampingRatioMZ[i,1];
		out_pars_MZ[i,3] = freqMZ[i,2];
		out_pars_MZ[i,4] = dampingRatioMZ[i,2];
		for(j in 1:2){
			A[2,1] = etaMZ[i,j];
			A[2,2] = zetaMZ[i,j];
			A[1,2] = 1.0;
			Q[2,2] = sigmaMZ[i,j]; 
			discMats = discretize(A, B, Q, dt);
			Ad = discMats[1];
			Bd = discMats[2];
			Qd = discMats[3];
			x = x0;
			P = P0;
			for(t in 1:T){
				z = MZ[i,j][t]';
			//Prediction
				xp = Ad*x;
				Pp = Ad*P*Ad'+Qd;
			//Correction
				r = z - H*xp - muMZ[i,j];
				S = H*Pp*H'+R;
				K= (inverse(S)*(H*Pp))';
				x = xp + K * r;
				P = Pp - K * H * Pp;
				r ~ normal(0.0, sqrt(S[1,1]));
			}
		}
	}		
//DZ twins
	for(i in 1:nDZ){
		out_pars_DZ[i,1] = freqDZ[i,1];
		out_pars_DZ[i,2] = dampingRatioDZ[i,1];
		out_pars_DZ[i,3] = freqDZ[i,2];
		out_pars_DZ[i,4] = dampingRatioDZ[i,2];
		for(j in 1:2){
			A[2,1] = etaDZ[i,j];
			A[2,2] = zetaDZ[i,j];
			A[1,2] = 1.0;
			Q[2,2] = sigmaDZ[i,j]; 
			discMats = discretize(A, B, Q, dt);
			Ad = discMats[1];
			Bd = discMats[2];
			Qd = discMats[3];
			x = x0;
			P = P0;
			for(t in 1:T){
				z = DZ[i,j][t]';
			//Prediction
				xp = Ad*x;
				Pp = Ad*P*Ad'+Qd;
			//Correction
				r = z - H*xp - muDZ[i,j];
				S = H*Pp*H'+R;
				K= (inverse(S)*(H*Pp))';
				x = xp + K * r;
				P = Pp - K * H * Pp;
				r ~ normal(0.0, sqrt(S[1,1]));
			}
		}
	}		
//Hyperpriors	
	out_pars_MZ ~ multi_normal(out_Mu, SigMZ);
	out_pars_DZ ~ multi_normal(out_Mu, SigDZ);
//Priors
	x0 ~ normal(0, 1);
	to_vector(muMZ) ~ normal( 0.0, 1.0);
	to_vector(muDZ) ~ normal( 0.0, 1.0);
	freqMu ~ normal(0.0, 1.0);
	to_vector(sigmaMZ_r) ~ normal( 0.0, 4.0);
	to_vector(sigmaDZ_r) ~ normal( 0.0, 4.0);
	epsilon_r ~ normal( 0.0, 1.0);  
	loadings  ~ normal( 0.0, 1.0);  
	a_eta ~ normal( 0.0, 1.0);  
	c_eta ~ normal( 0.0, 1.0);  
	e_eta ~ normal( 0.1, 1.0);  
	a_zeta ~ normal( 0.0, 1.0);  
	c_zeta ~ normal( 0.0, 1.0);  
	e_zeta ~ normal( 0.1, 1.0); 
}

