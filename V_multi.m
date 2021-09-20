function [phi, O, O_prov, w, invO_matrix_up, invO_matrix_dn] = V_multi(phi, phi_T, N_up, N_par, N_det, O, O_prov, w, w_T, invO_matrix_up, invO_matrix_dn, aux_fld)
% function [phi, O, w, invO_matrix_up, invO_matrix_dn] = V(phi, phi_T, N_up, N_par, O, w, invO_matrix_up, invO_matrix_dn, aux_fld)
% Sample the auxiliary field over a single lattice site for a single walker and propagate that walker by the potential energy propagator exp(-deltau*V)
% Inputs:
%   phi: a single row in the matrix of a walker (corresponding to the amplitude of all electrons over a single lattice site in that walker)
%   phi_T: the matrix of the trial wave function
%   N_up: the number of spin up electrons
%   N_par: the total number of electrons
%   O: the overlap of the aforementio<ned walker
%   w: the weight of the aforementioned walker
%   invO_matrix_up: the inverse of the spin up sector of the walker's overlap matrix 
%   invO_matrix_dn: the inverse of the spin down sector of the walker's overlap matrix 
%   aux_fld: the 2x2 matrix containing all the possible values of the quantity exp(gamma*s(sigma)*x_i) (used in V.m only)
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
RR=zeros(2,2, N_det);
matone=ones(2,2);
temp1_up=zeros(N_det, N_up);
temp1_dn=zeros(N_det,N_par-N_up);
temp2_up=zeros(N_up,N_det);
temp2_dn=zeros(N_par-N_up,N_det);
OV=zeros(2,1);
%% Calculate the Green's function
for ii=1:N_det
    temp1_up(ii,:)=phi(1:N_up)*invO_matrix_up(:,:,ii);
    temp1_dn(ii,:)=phi(N_up+1:N_par)*invO_matrix_dn(:,:,ii);
    temp2_up(:,ii)=invO_matrix_up(:,:,ii)*phi_T(1,1:N_up,ii)';
    temp2_dn(:,ii)=invO_matrix_dn(:,:,ii)*phi_T(1,N_up+1:N_par,ii)';
    Gii(1)=temp1_up(ii,:)*phi_T(1,1:N_up,ii)';
    Gii(2)=temp1_dn(ii,:)*phi_T(1,N_up+1:N_par,ii)';
    RR(:,:,ii)=(aux_fld-matone).*horzcat(Gii,Gii)+matone;
end
%% Perform the importance sampling and propagate the walker
% compute overlaps
O_prov_ratio_temp=RR(1,:,:).*RR(2,:,:);

for ii=1:N_det
    OV=OV+O_prov_ratio_temp(1,:,ii)*O_prov(ii)*w_T(ii)/O;
end

O_ratio_temp_real=max(real(OV), zeros(1,2));
% the normalization for the importance-sampled pdf
sum_O_ratio_temp_real=O_ratio_temp_real(1)+O_ratio_temp_real(2);

% if both auxiliary fields lead to negative overlap then kill the walker
if sum_O_ratio_temp_real<=0
    w=0;
end
if w>0
    % Otherwise update the weight
    w=w*0.5*sum_O_ratio_temp_real;    
    
    if O_ratio_temp_real(1)/sum_O_ratio_temp_real>=rand
        x_spin=1;
    else
        x_spin=2;
    end    
    % propagates the walker with the chosen auxiliary field
    phi(1:N_up)=phi(1:N_up)*aux_fld(1,x_spin);
    phi(N_up+1:N_par)=phi(N_up+1:N_par)*aux_fld(2,x_spin);    
    % Update the overlap using Sherman-Morrison
    O=O*OV(x_spin);
    for ii=1:N_det
        invO_matrix_up(:,:,ii)=invO_matrix_up(:,:,ii)+(1-aux_fld(1,x_spin))/RR(1,x_spin,ii)*temp2_up(:,ii)*temp1_up(ii,:);
        invO_matrix_dn(:,:,ii)=invO_matrix_dn(:,:,ii)+(1-aux_fld(2,x_spin))/RR(2,x_spin,ii)*temp2_dn(:,ii)*temp1_dn(ii,:);
        O_prov(ii)=O_prov(ii)*real(O_prov_ratio_temp(1,x_spin,ii));
    end
end
end