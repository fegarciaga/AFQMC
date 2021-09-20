function H=H_afm(Lx, Ly, H_k)
    % creates array to convert from bigger lattice to sub domain
    H=zeros(Lx*Ly/2, Lx*Ly/2);
    conv=zeros(Lx*Ly/2,1);
    r=0;
    r_alt=0;
    for jy=1:Ly
        for ix=1:Lx
            r=r+1;
            if mod(ix,2)~=0
                if mod(jy,2)~=0
                    r_alt=r_alt+1;  
                    conv(r_alt)=r;
                end
            end
            if mod(ix,2)==0
                if mod(jy,2)==0
                    r_alt=r_alt+1;
                    conv(r_alt)=r;
                end
            end
        end
    end
    H=H_k(conv, conv);
end