function H_res=H_col(Lx, kx, tx)
    k_x=pi*sqrt(-1)*kx;
    H_res=zeros(Lx, Lx);
    for i=1:Lx
        if i==1
            H_res(1,Lx)=H_res(1,Lx)-tx*exp(k_x);
            H_res(i,i+1)=H_res(i,i+1)-tx;
        elseif i==Lx
            H_res(Lx,1)=H_res(Lx,1)-tx*exp(-k_x);
            H_res(i,i-1)=H_res(i,i-1)-tx;
        else
            H_res(i,i+1)=-tx;
            H_res(i,i-1)=-tx;
        end
    end
end