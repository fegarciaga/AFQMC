% Script to initialize internal quantities
%
% Huy Nguyen, Hao Shi, Jie Xu and Shiwei Zhang
% ©2014 v1.0
% Package homepage: http://cpmc-lab.wm.edu
% Distributed under the <a href="matlab: web('http://cpc.cs.qub.ac.uk/licence/licence.html')">Computer Physics Communications Non-Profit Use License</a>
% Any publications resulting from either applying or building on the present package 
%   should cite the following journal article (in addition to the relevant literature on the method):
% "CPMC-Lab: A Matlab Package for Constrained Path Monte Carlo Calculations" Comput. Phys. Commun. (2014)

%% Check the validity of user inputs
validation;

%% Initialize internal quantities
N_sites=Lx*Ly*Lz;
N_par=N_up+N_dn;
% form the one-body kinetic Hamiltonian
H_k=H_K_alt(Lx, Ly, Lz, kx, ky, kz, tx, ty, tz, tx2, ty2, tz2);
% the matrix of the operator exp(-deltau*K/2)
Proj_k_half = expm(-0.5*deltau*H_k); 
%% Initialize the trial multideterminant wave function and calculate the ensemble's initial energy 
% initialize multi-determinat
Phi_T=zeros(N_sites, N_par, N_det);

% Fill the multi-determinant with randomly filled resonating valence bonds
% wavefunctions so long range correlations are present
for ii=1:N_det
    for jj=1:N_sites/2
        Phi_T=RVB(Phi_T, N_sites, N_up, ii, jj);
    end
end
% initialize weights of multideterminant trial wavefunction
w_T=1/sqrt(N_det)*ones(N_det, 1);
% initialize matrix for initial green function
invO_matrix_up=zeros(N_up, N_up, N_det);
invO_matrix_dn=zeros(N_dn, N_dn, N_det);

% Calculates inverse matrices
for ii=1:N_det
    invO_matrix_up(:,:,ii)=inv(Phi_T(:,1:N_up,ii)'*Phi_T(:,1:N_up,ii));
    invO_matrix_dn(:,:,ii)=inv(Phi_T(:,N_up+1:N_par,ii)'*Phi_T(:,N_up+1:N_par,ii));
end
E_V=0;
E_K=0;
% Calculates kinetic energy and potential energy
for ii=1:N_det
    temp_up=Phi_T(:,1:N_up,ii)*invO_matrix_up(:,:,ii);
    temp_dn=Phi_T(:,N_up+1:N_par,ii)*invO_matrix_dn(:,:,ii);
    G_up=temp_up*Phi_T(:,1:N_up,ii)';
    G_dn=temp_dn*Phi_T(:,N_up+1:N_par,ii)';
    E_K=E_K+w_T(ii)^2*sum(sum(H_k.'.*(G_up+G_dn)));
    n_r_up=diag(Phi_T(:,1:N_up,ii)*(Phi_T(:,1:N_up,ii))');
    n_r_dn=diag(Phi_T(:,N_up+1:N_par,ii)*(Phi_T(:,N_up+1:N_par,ii))');
    E_V=E_V+w_T(ii)^2*U*n_r_up'*n_r_dn;
    display(G_dn);
    display(G_up);
end
% the total energy of the trial wave function = the initial trial energy
E_T = E_K+E_V;
display(Phi_T);
%% Assemble the initial population of walkers
Phi=zeros(N_sites,N_par,N_wlk);
% initiate each walker to be the trial wave function
for i=1:N_wlk
    % Phi(:,:,i) is the ith walker. Each is a matrix of size N_sites by N_par
    % The left N_sites by N_up block is the spin up sector
    % The rest is the spin down sector
    % They are propagated independently and only share the auxiliary field
    Phi(:,:,i)=Phi_T(:,:,mod(i,N_det)+1); 
end
% initiate the weight and overlap of each walker to 1
w=1/sqrt(N_det)*ones(N_wlk,1);
O=1/sqrt(N_det)*ones(N_wlk,1);
O_prov=zeros(N_wlk, N_det);
for ii=1:N_wlk
    O_prov(ii,mod(ii,N_det)+1)=1;
end
% the arrays that store the energy and weight at each block
E_blk=zeros(N_blk,1);
W_blk=zeros(N_blk,1);

%% initialize auxiliary filed constants
% exponent of the prefactor exp(-deltau*(-E_T)) in the ground state projector 
% fac_norm also include -0.5*U*(N_up+N_dn), the exponent of the prefactor in the Hirsch transformation
fac_norm=(real(E_T)-0.5*U*N_par)*deltau; 
%gamma in Hirsch's transformation
gamma=acosh(exp(0.5*deltau*U)); 
% aux_fld is the 2x2 matrix containing all the possible values of the quantity exp(-gamma*s(sigma)*x_i)
aux_fld=zeros(2,2); 
% The first index corresponds to spin up or down
% The second index corresponds to the auxiliary field x_i=1 or x_i=-1
for i=1:2
    for j=1:2
        aux_fld(i,j)=exp(gamma*(-1)^(i+j));
    end
end
%% filename to be saved
savedFileName=strcat(int2str(Lx),'x',int2str(Ly),'x',int2str(Lz),'_',int2str(N_up),'u',int2str(N_dn),'d_U',num2str(U, '%4.2f'),'_kx',num2str(kx,'%+7.4f'),'_ky',num2str(ky,'%+7.4f'),'_kz',num2str(kz,'%+7.4f'),'_Nwlk_',int2str(N_wlk),suffix,'.mat');

%% randomize the random number generator seed based on the current time
rng('default');