function H=H_K_alt(Lx,Ly,Lz,kx,ky,kz,tx1,ty1,tz1,tx2,ty2,tz2)
% function H=H_K(Lx,Ly,Lz,kx,ky,kz,tx,ty,tz)
% Generate the one-body kinetic term of the Hubbard Hamiltonian with the given parameters
% Input:
%   Lx: The number of lattice sites in the x direction.
%   Ly: The number of lattice sites in the y direction.
%   Lz: The number of lattice sites in the z direction.
%   kx: The x component of the twist angle in TABC (twist-averaging boundary conditions)
%   ky: The y component of the twist angle in TABC
%   kz: The z component of the twist angle in TABC
%   tx1: The hopping amplitude between nearest-neighbor sites in the x direction
%   ty1: The hopping amplitude between nearest neighbor sites in the y direction
%   tz1: The hopping amplitude between nearest neighbor sites in the y direction
%   tx2: The hopping amplitude between next nearest neighbor sites in the
%   y-z plane
%   ty2: The hopping amplitude between next nearest neighbor sites in the
%   x-z plane
%   tz2: The hopping amplitude between next nearest neighbor sites in the
%   x-y plane
%
% Output
%   H: The one-body kinetic Hamiltonian in the form of a square matrix of size (Lx*Ly*Lz) 
%
% Adapted from:
% Huy Nguyen, Hao Shi, Jie Xu and Shiwei Zhang
% ©2014 v1.0
% Package homepage: http://cpmc-lab.wm.edu
% "CPMC-Lab: A Matlab Package for Constrained Path Monte Carlo Calculations" Comput. Phys. Commun. (2014)

    r=0;
    N_sites=Lx*Ly*Lz;   
    kx=sqrt(-1)*pi*kx;
    ky=sqrt(-1)*pi*ky;
    kz=sqrt(-1)*pi*kz;
    H=zeros(N_sites,N_sites);

    for mz=1:Lz
        for iy=1:Ly
            for jx=1:Lx
                r=r+1;      % r=(iy-1)*Lx+jx;
                if Lx~=1
                    if jx==1
                        H(r,r+Lx-1)=H(r,r+Lx-1)-tx1*exp(kx);
                        H(r,r+1)=H(r,r+1)-tx1;
                    elseif jx==Lx
                        H(r,r-1)=H(r,r-1)-tx1;
                        H(r,r+1-Lx)=H(r,r+1-Lx)-tx1*exp(-kx);
                    else
                        H(r,r-1)=-tx1;
                        H(r,r+1)=-tx1;
                    end
                end
                
                if Ly~=1
                    if iy==1
                        H(r,r+(Ly-1)*Lx)=H(r,r+(Ly-1)*Lx)-ty1*exp(ky);
                        H(r,r+Lx)=H(r,r+Lx)-ty1;
                    elseif iy==Ly
                        H(r,r-Lx)=H(r,r-Lx)-ty1;
                        H(r,r-(Ly-1)*Lx)=H(r,r-(Ly-1)*Lx)-ty1*exp(-ky);
                    else
                        H(r,r-Lx)=-ty1;
                        H(r,r+Lx)=-ty1;
                    end
                end

                if Lz~=1
                    if mz==1
                        H(r,r+(Lz-1)*Lx*Ly) = H(r,r+(Lz-1)*Lx*Ly) - tz1*exp(kz);
                        H(r,r+Lx*Ly)= H(r,r+Lx*Ly) - tz1;
                    elseif mz==Lz
                        H(r,r-Lx*Ly) = H(r,r-Lx*Ly) - tz1;
                        H(r,r-(Lz-1)*Lx*Ly) = H(r,r-(Lz-1)*Lx*Ly) - tz1*exp(-kz);
                    else
                        H(r,r-Lx*Ly)=-tz1;
                        H(r,r+Lx*Ly)=-tz1;
                    end
                end
                
                if Ly>3
                    if Lz==1
                        if Lx==1
                            if iy-2<1
                                H(r,r+(Ly-2)*Lx)=H(r,r+(Ly-2)*Lx)-ty2*exp(ky);
                                H(r,r+2*Lx)=-ty2;
                        
                            elseif iy+2>Ly
                                H(r,r-2*Lx)=-ty2;
                                H(r,r-(Ly-2)*Lx)=H(r,r-(Ly-2)*Lx)-ty2*exp(-ky);
                            
                            else
                                H(r,r+2*Lx)=-ty2;
                                H(r,r-2*Lx)=-ty2;
                            end
                        end
                    end
                end
                
                if Lx>3
                    if Lz==1
                        if Ly==1
                            if jx-2<1
                                H(r,r+(Lx-2))=H(r,r+(Lx-2))-tx2*exp(kx);
                                H(r,r+2)=-tx2;
                        
                            elseif jx+2>Lx
                                H(r,r-2)=-tx2;
                                H(r,r-(Lx-2))=H(r,r-(Lx-2))-tx2*exp(-kx);
                            
                            else
                                H(r,r+2)=-tx2;
                                H(r,r-2)=-tx2;
                            end
                        end
                    end
                end
                
                if Lz>3
                    if Ly==1
                        if Lx==1
                            if mz-2<1
                                H(r,r+(Lz-2)*Lx*Ly)=H(r,r+(Lz-2)*Lx*Ly)-tz2*exp(kz);
                                H(r,r+2*Lx*Ly)=-tz2;
                        
                            elseif mz+2>Lz
                                H(r,r-2*Lx*Ly)=-tz2;
                                H(r,r-(Lz-2)*Lx*Ly)=H(r,r-(Lz-2)*Lx*Ly)-tz2*exp(-kz);
                            
                            else
                                H(r,r+2*Lx*Ly)=-tz2;
                                H(r,r-2*Lx*Ly)=-tz2;
                            end
                        end
                    end
                end
                
                if Lx>1
                    if Ly>1
                        if jx==1
                           if iy==1
                              H(r,r+Lx-1+(Ly-1)*Lx)=H(r,r+Lx-1+(Ly-1)*Lx)-tz2*exp(kx)*exp(ky);
                              H(r,r+2*Lx-1)=H(r,r+2*Lx-1)-tz2*exp(kx);
                              H(r,r+1+(Ly-1)*Lx)=H(r,r+1+(Ly-1)*Lx)-tz2*exp(ky);
                              H(r,r+1+Lx)=H(r,r+1+Lx)-tz2;
                           
                           elseif iy==Ly
                               H(r,r+(Lx-1)-(Ly-1)*Lx)=H(r,r+(Lx-1)-(Ly-1)*Lx)-tz2*exp(kx)*exp(-ky);
                               H(r,r+1-(Ly-1)*Lx)=H(r,r+1-(Ly-1)*Lx)-tz2*exp(-ky);
                               H(r,r-1)=H(r,r-1)-tz2*exp(kx);
                               H(r,r+1-Lx)=H(r,r+1-Lx)-tz2;
                               
                           else
                               H(r,r+2*Lx-1)=H(r,r+2*Lx-1)-tz2*exp(kx);
                               H(r,r-1)=H(r,r-1)-tz2*exp(kx);
                               H(r,r+1+Lx)= H(r,r+1+Lx)-tz2;
                               H(r,r+1-Lx)= H(r,r+1-Lx)-tz2;
                           end
                        elseif jx==Lx
                            if iy==1
                                H(r,r-(Lx-1)+(Ly-1)*Lx)=H(r,r-(Lx-1)+(Ly-1)*Lx)-tz2*exp(-kx)*exp(ky);
                                H(r,r+1)=H(r,r+1)-tz2*exp(-kx);
                                H(r,r+(Ly-1)*Lx-1)=H(r,r+(Ly-1)*Lx-1)-tz2*exp(ky);
                                H(r,r-1+Lx)=H(r,r-1+Lx)-tz2;
                                
                            elseif iy==Ly
                                H(r,r-(Lx-1)-(Ly-1)*Lx)=H(r,r-(Lx-1)-(Ly-1)*Lx)-tz2*exp(-kx)*exp(-ky);
                                H(r,r-(2*Lx-1))= H(r,r-(2*Lx-1))-tz2*exp(-kx);
                                H(r,r-1-(Ly-1)*Lx)= H(r,r-1-(Ly-1)*Lx)-tz2*exp(-ky);
                                H(r,r-1-Lx)=H(r,r-1-Lx)-tz2;
                            else
                                H(r,r+1)=H(r,r+1)-tz2*exp(-kx);
                                H(r,r-(2*Lx-1))=H(r,r-(2*Lx-1))-tz2*exp(-kx);
                                H(r,r-1-Lx)=H(r,r-1-Lx)-tz2;
                                H(r,r-1+Lx)=H(r,r-1+Lx)-tz2;
                            end
                        
                        elseif iy==1
                            if mod(jx,Lx)>1
                                H(r,r-1+(Ly-1)*Lx)=H(r,r-1+(Ly-1)*Lx)-tz2*exp(ky);
                                H(r,r+1+(Ly-1)*Lx)=H(r,r+1+(Ly-1)*Lx)-tz2*exp(ky);
                                H(r,r-1+Lx)=H(r,r-1+Lx)-tz2;
                                H(r,r+1+Lx)=H(r,r+1+Lx)-tz2;
                            end
                        elseif iy==Ly
                            if mod(jx,Lx)>1
                                H(r,r-1-(Ly-1)*Lx)= H(r,r-1-(Ly-1)*Lx)-tz2*exp(-ky);
                                H(r,r+1-(Ly-1)*Lx)= H(r,r+1-(Ly-1)*Lx)-tz2*exp(-ky);
                                H(r,r-1-Lx)=H(r,r-1-Lx)-tz2;
                                H(r,r+1-Lx)=H(r,r+1-Lx)-tz2;
                            end
                        else
                            H(r,r+1+Lx)=-tz2;
                            H(r,r-1+Lx)=-tz2;
                            H(r,r+1-Lx)=-tz2;
                            H(r,r-1-Lx)=-tz2;
                        end
                    end
                end
            
                if Ly>1
                    if Lz>1
                        if mz==1
                           if iy==1
                              H(r,r+(Lz-1)*Lx*Ly+(Ly-1)*Lx)=H(r,r+(Lz-1)*Lx*Ly+(Ly-1)*Lx)-tx2*exp(ky)*exp(kz);
                              H(r,r+Lx+Lx*Ly*(Lz-1))= H(r,r+Lx+Lx*Ly*(Lz-1))-tx2*exp(kz);
                              H(r,r+Lx*Ly+(Ly-1)*Lx)=H(r,r+Lx*Ly+(Ly-1)*Lx)-tx2*exp(ky);
                              H(r,r+Lx*Ly+Lx)=H(r,r+Lx*Ly+Lx)-tx2;
                           
                           elseif iy==Ly
                               H(r,r+(Lz-1)*Lx*Ly-(Ly-1)*Lx)=H(r,r+(Lz-1)*Lx*Ly-(Ly-1)*Lx)-tx2*exp(kz)*exp(-ky);
                               H(r,r+Lx*Ly-(Ly-1)*Lx)=H(r,r+Lx*Ly-(Ly-1)*Lx)-tx2*exp(-ky);
                               H(r,r-Lx+(Lz-1)*Lx*Ly)=H(r,r-Lx+(Lz-1)*Lx*Ly)-tx2*exp(kz);
                               H(r,r+Lx*Ly-Lx)=H(r,r+Lx*Ly-Lx)-tx2;
                               
                           else
                               H(r,r+Lx+(Lz-1)*Lx*Ly)=H(r,r+Lx+(Lz-1)*Lx*Ly)-tx2*exp(kz);
                               H(r,r+(Lz-1)*Lx*Ly-Lx)=H(r,r+(Lz-1)*Lx*Ly-Lx)-tx2*exp(kz);
                               H(r,r+Lx*Ly+Lx)= H(r,r+Lx*Ly+Lx)-tx2;
                               H(r,r+Lx*Ly-Lx)= H(r,r+Lx*Ly-Lx)-tx2;
                           end
                        elseif mz==Lz
                            if iy==1
                                H(r,r-(Lz-1)*Lx*Ly+(Ly-1)*Lx)=H(r,r-(Lz-1)*Lx*Ly+(Ly-1)*Lx)-tx2*exp(-kz)*exp(ky);
                                H(r,r-(Lz-1)*Lx*Ly+Lx)=H(r,r-(Lz-1)*Lx*Ly+Lx)-tx2*exp(-kz);
                                H(r,r+Lx)=H(r,r+Lx)-tx2*exp(ky);
                                H(r,r-Lx*Ly+Lx)=H(r,r-Lx*Ly+Lx)-tx2;
                                
                            elseif iy==Ly
                                H(r,r-(Lz-1)*Lx*Ly-(Ly-1)*Lx)=H(r,r-(Lx-1)*Lx*Ly-(Ly-1)*Lx)-tx2*exp(-kz)*exp(-ky);
                                H(r,r-(Lz-1)*Lx*Ly-Lx)= H(r,r-(Lz-1)*Lx*Ly-Lx)-tx2*exp(-kz);
                                H(r,r-Lx*Ly-(Ly-1)*Lx)= H(r,r-Lx*Ly-(Ly-1)*Lx)-tx2*exp(-ky);
                                H(r,r-Lx*Ly-Lx)=H(r,r-Lx*Ly-Lx)-tx2;
                            else
                                H(r,r-(Lz-1)*Lx*Ly+Lx)=H(r,r-(Lz-1)*Lx*Ly+Lx)-tx2*exp(-kz);
                                H(r,r-(Lz-1)*Lx*Ly-Lx)=H(r,r-(Lz-1)*Lx*Ly-Lx)-tx2*exp(-kz);
                                H(r,r-Lx*Ly-Lx)=H(r,r-Lx*Ly-Lx)-tx2;
                                H(r,r-Lx*Ly+Lx)=H(r,r-Lx*Ly+Lx)-tx2;
                            end
                        
                        elseif iy==1
                            if mod(mz,Lz)>1
                                H(r,r-Lx*Ly+(Ly-1)*Lx)=H(r,r-Lx*Ly+(Ly-1)*Lx)-tx2*exp(ky);
                                H(r,r+Lx*Ly+(Ly-1)*Lx)=H(r,r+Lx*Ly+(Ly-1)*Lx)-tx2*exp(ky);
                                H(r,r-Lx*Ly+Lx)=H(r,r-Lx*Ly+Lx)-tx2;
                                H(r,r+Lx*Ly+Lx)=H(r,r+Lx*Ly+Lx)-tx2;
                            end
                        elseif iy==Ly
                            if mod(mz,Lz)>1
                                H(r,r-Lx*Ly-(Ly-1)*Lx)= H(r,r-Lx*Ly-(Ly-1)*Lx)-tx2*exp(-ky);
                                H(r,r+Lx*Ly-(Ly-1)*Lx)= H(r,r+lx*Ly-(Ly-1)*Lx)-tx2*exp(-ky);
                                H(r,r-Lx*Ly-Lx)=H(r,r-Lx*Ly-Lx)-tx2;
                                H(r,r+Lx*Ly-Lx)=H(r,r+Lx*Ly-Lx)-tx2;
                            end
                        else
                            H(r,r+Lx*Ly+Lx)=-tx2;
                            H(r,r-Lx*Ly+Lx)=-tx2;
                            H(r,r+Lx*Ly-Lx)=-tx2;
                            H(r,r-Lx*Ly-Lx)=-tx2;
                        end
                    end
                end
                
                if Lx>1
                    if Lz>1
                        if jx==1
                            if mz==1
                                H(r,r+Lx-1+(Lz-1)*Lx*Ly)=H(r,r+Lx-1+(Lz-1)*Lx*Ly)-ty2*exp(kx)*exp(kz);
                                H(r,r+Lx-1+Lx*Ly)=H(r,r+Lx-1+Lx*Ly)-ty2*exp(kx);
                                H(r,r+1+(Lz-1)*Lx*Ly)=H(r,r+1+(Lz-1)*Lx*Ly)-ty2*exp(kz);
                                H(r,r+1+Lx*Ly)=H(r,r+1+Lx*Ly)-ty2;   
                           
                            elseif mz==Lz
                                H(r,r+(Lx-1)-(Lz-1)*Lx*Ly)=H(r,r+(Lx-1)-(Lz-1)*Lx*Ly)-ty2*exp(kx)*exp(-kz);
                                H(r,r+1-(Lz-1)*Lx*Ly)=H(r,r+1-(Lz-1)*Lx*Ly)-ty2*exp(-kz);
                                H(r,r+(Lx-1)-Lx*Ly)=H(r,r+(Lx-1)-Lx*Ly)-ty2*exp(kx);
                                H(r,r+1-Lx*Ly)=H(r,r+1-Lx*Ly)-ty2;
                               
                            else
                                H(r,r+Lx-1+Lx*Ly)=H(r,r+Lx-1+Lx*Ly)-ty2*exp(kx);
                                H(r,r+(Lx-1)-Lx*Ly)=H(r,r+(Lx-1)-Lx*Ly)-ty2*exp(kx);
                                H(r,r+1+Lx*Ly)= H(r,r+1+Lx*Ly)-ty2;
                                H(r,r+1-Lx*Ly)= H(r,r+1-Lx*Ly)-ty2;
                            end
                        elseif jx==Lx
                            if mz==1
                                H(r,r-(Lx-1)+(Lz-1)*Lx*Ly)=H(r,r-(Lx-1)+(Lz-1)*Lx*Ly)-ty2*exp(-kx)*exp(kz);
                                H(r,r-(Lx-1)+Lx*Ly)=H(r,r-(Lx-1)+Lx*Ly)-ty2*exp(-kx);
                                H(r,r+(Lz-1)*Lx*Ly-1)=H(r,r+(Lz-1)*Lx*Ly-1)-ty2*exp(kz);
                                H(r,r-1+Lx*Ly)=H(r,r-1+Lx*Ly)-ty2;
                                
                            elseif mz==Lz
                                H(r,r-(Lx-1)-(Lz-1)*Lx*Ly)=H(r,r-(Lx-1)-(Lz-1)*Lx*Ly)-ty2*exp(-kx)*exp(-kz);
                                H(r,r-(Lx-1)-Lx*Ly)= H(r,r-(Lx-1)-Lx*Ly)-ty2*exp(-kx);
                                H(r,r-1-(Lz-1)*Lx*Ly)= H(r,r-1-(Lz-1)*Lx*Ly)-ty2*exp(-kz);
                                H(r,r-1-Lx*Ly)=H(r,r-1-Lx*Ly)-ty2;
                               
                            else
                                H(r,r+1)=H(r,r+1)-tz2*exp(-kx);
                                H(r,r-(Lx-1)+Lx*Ly)=H(r,r-(Lx-1)+Lx*Ly)-ty2*exp(-kx);
                                H(r,r-1-Lx*Ly)=H(r,r-1-Lx*Ly)-ty2;
                                H(r,r-1+Lx*Ly)=H(r,r-1+Lx*Ly)-ty2;
                            end
                        
                        elseif mz==1
                            if mod(jx,Lx)>1
                                H(r,r-1+(Lz-1)*Lx*Ly)=H(r,r-1+(Lz-1)*Lx*Ly)-ty2*exp(kz);
                                H(r,r+1+(Lz-1)*Lx*Ly)=H(r,r+1+(Lz-1)*Lx*Ly)-ty2*exp(kz);
                                H(r,r-1+Lx*Ly)=H(r,r-1+Lx*Ly)-ty2;
                                H(r,r+1+Lx*Ly)=H(r,r+1+Lx*Ly)-ty2;
                            end
                        elseif mz==Ly
                            if mod(jx,Lx)>1
                                H(r,r-1-(Lz-1)*Lx*Ly)= H(r,r-1-(Lz-1)*Lx*Ly)-ty2*exp(-kz);
                                H(r,r+1-(Lz-1)*Lx*Ly)= H(r,r+1-(Lz-1)*Lx*Ly)-ty2*exp(-kz);
                                H(r,r-1-Lx*Ly)=H(r,r-1-Lx*Ly)-ty2;
                                H(r,r+1-Lx*Ly)=H(r,r+1-Lx*Ly)-ty2;
                            end
                        else
                            H(r,r+1+Lx*Ly)=-ty2;
                            H(r,r-1+Lx*Ly)=-ty2;
                            H(r,r+1-Lx*Ly)=-ty2;
                            H(r,r-1-Lx*Ly)=-ty2;
                        end
                    end
                end
            end
        end
    end