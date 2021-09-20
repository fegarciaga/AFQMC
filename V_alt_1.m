function [phi, O, w, invO_matrix_up, invO_matrix_dn] = V_alt_1(phi, phi_T, N_up, N_par, O, w, invO_matrix_up, invO_matrix_dn, deltau, U)
% function [phi, O, w, invO_matrix_up, invO_matrix_dn] = V(phi, phi_T, N_up, N_par, O, w, invO_matrix_up, invO_matrix_dn, aux_fld)
% Sample the potential energy term using the continuous
% Hubbard-Stratonovich transformation
% Inputs:
%   phi: a single row in the matrix of a walker (corresponding to the amplitude of all electrons over a single lattice site in that walker)
%   phi_T: the matrix of the trial wave function
%   N_up: the number of spin up electrons
%   N_par: the total number of electrons
%   O: the overlap of the aforementioned walker
%   w: the weight of the aforementioned walker
%   invO_matrix_up: the inverse of the spin up sector of the walker's overlap matrix 
%   invO_matrix_dn: the inverse of the spin down sector of the walker's overlap matrix
% Outputs:
%   phi: the propagated row of the aforementioned walker
%   O: the overlap after propagation of the aforementioned walker
%   w: the weight after propagation of the aforementioned walker
%   invO_matrix_up: the updated inverse of the spin up sector of the walker's overlap matrix 
%   invO_matrix_dn: the updated inverse of the spin down sector of the walker's overlap matrix 
%
% Huy Nguyen, Hao Shi, Jie Xu and Shiwei Zhang
% ©2014 v1.0
% Package homepage: http://cpmc-lab.wm.edu
% Distributed under the <a href="matlab: web('http://cpc.cs.qub.ac.uk/licence/licence.html')">Computer Physics Communications Non-Profit Use License</a>
% Any publications resulting from either applying or building on the present package 
%   should cite the following journal article (in addition to the relevant literature on the method):
% "CPMC-Lab: A Matlab Package for Constrained Path Monte Carlo Calculations" Comput. Phys. Commun. (2014)

%% Pre-allocate matrices:
Gii=zeros(2,1);
RR=ones(2,1);
matone=RR;

%% Calculate the Green's function
temp1_up=phi(1:N_up)*invO_matrix_up;
temp1_dn=phi(N_up+1:N_par)*invO_matrix_dn;
temp2_up=invO_matrix_up*phi_T(1:N_up)';
temp2_dn=invO_matrix_dn*phi_T(N_up+1:N_par)';
Gii(1)=temp1_up*phi_T(1:N_up)';
Gii(2)=temp1_dn*phi_T(N_up+1:N_par)';
%% Perform the sampling and propagate the walker
F=-1i*sqrt(U*deltau/2)*(Gii(1)+Gii(2));
sample=randn+F;  
aux=zeros(1,2);
aux(1)=exp(1i*sqrt(U*deltau/2)*sample);
aux(2)=exp(1i*sqrt(U*deltau/2)*sample);
RR=(aux-matone).*Gii+matone;  
O_rat=RR(1)*RR(2);
x_phaseless=max(0,cos(angle(O_rat)));
if x_phaseless==0
   display(cos(angle(O_rat))); 
end
w=abs(O_rat*exp(sample*F-0.5*F^2))*x_phaseless;
O=O*O_rat;    
% propagates the walker with the chosen auxiliary field
phi(1:N_up)=phi(1:N_up)*aux(1);
phi(N_up+1:N_par)=phi(N_up+1:N_par)*aux(2);    
% Update the overlap using Sherman-Morrison
invO_matrix_up=invO_matrix_up+(1-aux(1))/RR(1)*temp2_up*temp1_up;
invO_matrix_dn=invO_matrix_dn+(1-aux(2))/RR(2)*temp2_dn*temp1_dn;
end