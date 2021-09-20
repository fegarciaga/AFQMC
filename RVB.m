function [Phi_T]=RVB(Phi_T, N_sites, N_up, i, j)
    flag=0;
    while flag==0
        i1=randi([1,N_sites]);
        i2=randi([1,N_sites]);
        if i1~=i2
            if sum(Phi_T(i1, :, i))==0
                if sum(Phi_T(i2, :, i))==0
                    Phi_T(i1, j, i)=1/sqrt(2);    
                    Phi_T(i2, j, i)=1/sqrt(2);  
                    Phi_T(i1, j+N_up, i)=1/sqrt(2);    
                    Phi_T(i2, j+N_up, i)=1/sqrt(2); 
                    flag=1;
                end
            end
        end
    end
end