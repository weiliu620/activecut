#include "commonalt.h"

twister_base_gen_type mygenerator;

using boost::math::beta_distribution;

int generateVonMise(vnl_vector<float>& randomVector, float lambda)
{
     // By "Computer Generation of Distributions on the m-sphere, Gary Ulrich, 1984

     unsigned short idx = 0;
     float gammaVar1 = 0, gammaVar2 = 0;

     unsigned short m = randomVector.size();
     float U = 0, Z = 0, W = 0, T = 0;

     // Uniform real random generator.
     boost::uniform_01<> uni_dist; // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_01<> > uni(mygenerator, uni_dist);

     // construct uniform on (m-1) dimensional sphere.
     boost::uniform_on_sphere<float> sphere_dist(m-1);
     boost::variate_generator<twister_base_gen_type&, boost::uniform_on_sphere<float> > sphere_generator(mygenerator, sphere_dist);

     // construct gamma distribution (boost random lib)
     boost::gamma_distribution<float> gamma_distr((m-1)/2);
     boost::variate_generator<twister_base_gen_type&, boost::gamma_distribution<float> > gamma_generator(mygenerator, gamma_distr);

     beta_distribution<> beta_distr((m-1)/2, (m-1)/2);
     
     // step 0, initialization. 
     float b = (-2*lambda + sqrt(4*lambda*lambda + (m-1)*(m-1)))/(m-1);
     float a = ((m-1) + 2*lambda + sqrt(4*lambda*lambda + (m-1)*(m-1)))/4;
     float d = (4*a*b) / (1+b) - (m-1) * log(m-1);

     // Generate beta distribution sample
     do {
	  // U = uni(); // in order to generate beta distribution sample.
	  // Z = quantile(beta_distr, U);
	  
	  // Generate two rvs from Gamma to get Beta distribution.
	  gammaVar1 = gamma_generator();
	  gammaVar2 = gamma_generator();
	  Z = gammaVar1/(gammaVar1 + gammaVar2);

	  U = uni(); // this time, U is used for von Mise.
	  
	  W = (1 - (1+b)*Z)/(1 + (1-b)*Z);
	  T = (2 * a * b) / (1 + (1-b)*Z);
     }
     while((m-1) * log(T) - T + d < log(U));

     // Step 3. 
     std::vector<float> uniformDirection = sphere_generator();
     const float scale = sqrt(1 - W*W);
     for (idx = 0; idx < uniformDirection.size(); idx++) {
	  randomVector[idx] = uniformDirection[idx] * scale;
     }
     // last elements is W.
     randomVector[idx] = W;
     
     return (0);
}

int generateVonMiseWood(vnl_vector<float>& randomVector, float kappa)
{
     // Adapted from "Mondeling Data using Directional Distributions, Inderjit S. Dhillon et. al.

     float gammaVar1 = 0, gammaVar2 = 0;
     unsigned short idx = 0;
     unsigned short d = randomVector.size();     

     // construct uniform on (m-1) dimensional sphere.
     boost::uniform_on_sphere<float> sphere_dist(d-1);
     boost::variate_generator<twister_base_gen_type&, boost::uniform_on_sphere<float> > sphere_generator(mygenerator, sphere_dist);

     // construct gamma distribution (boost random lib)
     boost::gamma_distribution<float> gamma_distr((d-1)/2);
     boost::variate_generator<twister_base_gen_type&, boost::gamma_distribution<float> > gamma_generator(mygenerator, gamma_distr);

     // Uniform real random generator.
     boost::uniform_01<> uni_dist; // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_01<> > uni(mygenerator, uni_dist);


     beta_distribution<> beta_distr((d-1)/2, (d-1)/2);
     
     // step 2.
     float t1 = sqrt(4*kappa*kappa + (d-1)*(d-1));
     // Step 3.
     float b = (-2 * kappa + t1)/(d-1);
     // step 4.
     float x0 = (1 - b)/ (1 + b);

     // Step 6.

     // Step 7. 
     float c = kappa * x0 + (d-1) * log(1 - x0*x0);
     
     float t = 0, u = 0, z = 0, w = 0;
     
     t = -1000;
     u = 1;
     
     while(t < log(u)) {
	  gammaVar1 = gamma_generator();
	  gammaVar2 = gamma_generator();
	  z = gammaVar1/(gammaVar1 + gammaVar2);

	  u = uni(); // this time, U is used for von Mise.
	  w = (1 - (1+b)*z)/(1 - (1-b)*z);
	  // the paper "Modeling data using Directional Distributions"
	  // by Dhillon has a mistake here. We should add '-c'
	  // term. Refere to original Woods algorithm.
//	  t = kappa*w + (d-1) * log(1 - x0 * w); 
	  t = kappa*w + (d-1) * log(1 - x0 * w) - c;
     }
	  
     std::vector<float> uniformDirection = sphere_generator();
     const float scale = sqrt(1 - w*w);	  
     for (idx = 0; idx < uniformDirection.size(); idx++) {
	  randomVector[idx] = uniformDirection[idx] * scale;
     }

     randomVector[d-1] = w;

     return (0);
}

int getRotationMat(vnl_matrix<float>& rotationMat,
		   vnl_vector<float>& srcVec, 
		   vnl_vector<float>& destVec)
{
     // rotation matrix to transform srcVec to destVec. Adapted from
     // "Pivotal bootstrap Methods for k-Sample Problems in
     // Directional Statistics and Shape Analysis, G.J.A.Amaral
     // et. al. See page 699 sec 3.2.1.
     float alpha = acos(cos_angle(srcVec, destVec));
     vnl_vector<float> c = srcVec - inner_product(srcVec, destVec) * destVec;
     c.normalize();
     vnl_matrix<float> A = outer_product(destVec, c) - outer_product(c, destVec);
     
     rotationMat.set_identity();
     vnl_matrix<float> B = outer_product(destVec, destVec) + outer_product(c, c);
     rotationMat = rotationMat +  A * sin(alpha)  +   B * (cos(alpha) - 1);
} 
