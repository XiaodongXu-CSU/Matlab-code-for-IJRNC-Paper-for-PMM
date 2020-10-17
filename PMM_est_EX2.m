function  [cost,ryy] =  PMM_est(Ly,ryy_sample,N,P,W,N_out,N_in,Ndel,KC,MPC_M_fsr,MPC_ERO)


         Dnois_sys   =     cell(N_out,N_out); 
                 z   =     tf('z',1);
    Dnois_sys{1,1}   =     (1 - Ndel(11) * z^(-1) )  /(1 - Ndel(13) * z^(-1) );
    Dnois_sys{1,2}   =     0 * z;    
    Dnois_sys{2,1}   =     0 * z;    
    Dnois_sys{2,2}   =     (1 - Ndel(12) * z^(-1) )  /(1 - Ndel(14) * z^(-1) ); 
    
         for i_out   =   1  :  N_out
                                 for j_out  =   1  :  N_out
                                          Dsys   =     Dnois_sys{i_out,j_out}; 
                                       [H_e,t]   =     impulse(Dsys,W);
                                          hemp   =     H_e; 
        Nois_fir(i_out,W*(j_out-1)+1:W*j_out)    =     hemp(2:end);
                                 end
         end
         Nois_M_fir    =    cell(1,W);
 for mk  =   1  : W
                Ntemp = [];
            for j_out  =  1  :  N_out
                Ntemp = [Ntemp, Nois_fir(:,N*(j_out-1)+mk)];
            end    
       Nois_M_fir{mk} = Ntemp;
 end 



        CReal_g_sys    =     cell(N_out, N_in);
                 s     =     tf('s');

 
    CReal_g_sys{1,1}   =     11.77   *  Ndel(1)  * exp(-0.189* s)            /(10.32 * Ndel(5) *s + 1);
    CReal_g_sys{1,2}   =     - 4.77  *  Ndel(2)  * exp(-3.87 * Ndel(9)  * s) /(8.86  * Ndel(6) *s + 1);
    CReal_g_sys{2,1}   =     - 8.12  *  Ndel(3)  * exp(-4.63 * Ndel(10) * s) /(7.31  * Ndel(7) *s + 1);
    CReal_g_sys{2,2}   =     10.04   *  Ndel(4)                              /(10.34 * Ndel(8) *s + 1);      
       
       
       
 for i_out  =  1  :  N_out
            for j_in  =  1  : N_in
                                  dsys  =   c2d(CReal_g_sys{i_out,j_in},1); 
                 Dreal_sys{i_out,j_in}  =   dsys; 
                            [S_Real,t]  =   step(dsys,N); 
                            temp        =   S_Real';
    Real_fsr(i_out,N*(j_in-1)+1:N*j_in) =   temp(2:end) + MPC_ERO{i_out,j_in}'; 
            end
 end
 
           Real_M_fsr = cell(1,N); 
 for mk  =   1  : N
                Stemp = [];
            for j_in  =  1  :  N_in
               Stemp = [Stemp,Real_fsr(:,N*(j_in-1)+mk)];
            end    
      Real_M_fsr{mk} = Stemp; 
 end 
   

 for k   =   1 : N  
       Del_M_fsr{k}   =   Real_M_fsr{k}  - MPC_M_fsr{k};
 end
 
 
  Real_M_fir{1}   =   Real_M_fsr{1};
   MPC_M_fir{1}   =   MPC_M_fsr{1}; 
   Del_M_fir{1}   =   Del_M_fsr{1}; 
 
for k =  2 : N
    Real_M_fir{k}   =    Real_M_fsr{k}  -   Real_M_fsr{k-1}; 
     MPC_M_fir{k}   =     MPC_M_fsr{k}  -   MPC_M_fsr{k-1}; 
     Del_M_fir{k}   =    Real_M_fir{k}  -   MPC_M_fir{k}; 
end



    B   =   zeros(N_in,N_out);
for k   =   1: P
    B   =   B - KC(1:N_in,N_out*(k-1)+1:N_out*k);
end


C_temp = zeros(N_in,N_in);
for k = 1: P
    C_temp = C_temp + KC(1:N_in,N_out*(k-1)+1:N_out*k)* (MPC_M_fsr{k+1} + Del_M_fsr{1});
end
C{1} = eye(N_in) - C_temp;

for m =  2: N-P
    C_temp = zeros(N_in,N_in);
    for k = 1:P
       C_temp = C_temp + KC(1:N_in,N_out*(k-1)+1:N_out*k)* (MPC_M_fir{k+m} + Del_M_fir{m}); 
    end
   C{m} = - C_temp;   
end

for m =   N-P+1:N-1
    C_temp = zeros(N_in,N_in);
    for k = 1:N-m
       C_temp = C_temp + KC(1:N_in,N_out*(k-1)+1:N_out*k)* (MPC_M_fir{k+m} + Del_M_fir{m}); 
    end
    C_temp1 = zeros(N_in,N_in);
    for k = N-m+1:P
       C_temp1 = C_temp1 + KC(1:N_in,N_out*(k-1)+1:N_out*k)* (Del_M_fir{m}); 
    end
    
   C{m} = - C_temp-C_temp1;
    
end

 C_temp1 = zeros(N_in,N_in);
    for k = 1:P
       C_temp1 = C_temp1 + KC(1:N_in,N_out*(k-1)+1:N_out*k)* (Del_M_fir{N}); 
    end
C{N} = -C_temp1;


for i = 0 : N
    if i == 0
        alpha{i+1} = eye(N_in); 
    elseif i==1
        alpha{i+1} = C{i};
    elseif i<=N
        al_temp = zeros(N_in,N_in);
        for k = 0:i-1
          al_temp = al_temp +  alpha{k+1} * C{i-k};
        end
        alpha{i+1} =  al_temp; 
    elseif i>N
        al_temp = zeros(N_in,N_in);
        for k = i-N:i-1
          al_temp = al_temp + alpha{k+1} * C{i-k};  
        end
        alpha{i+1} =  al_temp;  
    end 
end


for  i  =   0 : W-1
    Gam_tem = zeros(N_in,N_out);
    for j =  0:i
        Gam_tem  =  Gam_tem  + alpha{j+1} * B * Nois_M_fir{i-j+1};    
    end 
    Gama{i+1} = Gam_tem;
end

for i =  0 : W-1
    
    if i == 0
        PI{i+1} = Nois_M_fir{i+1};     
    else
           y_eff_temp = zeros(N_out,N_out);
           for j = 1 : i
              y_eff_temp = y_eff_temp + Real_M_fir{j} * Gama{i-j+1};

           end
         PI{i+1} =  Nois_M_fir{i+1} +  y_eff_temp;
    end
end

%Ly  =  10; 
ryy = [];
for pi = 0 : Ly
    Rytemp = zeros(N_out,N_out);
    for i = 0 : W - pi-1
       Rytemp = Rytemp +  PI{i+1} * PI{i+pi+1}'; 
    end
    Ryy{pi+1} = Rytemp;
          ryy = [ryy,reshape(Ryy{pi+1}',N_out*N_out,1)];
end
 
dryy = ryy - ryy_sample;
DR_v = reshape(dryy,N_out*N_out*(Ly+1),1);

cost = dot(DR_v,DR_v);





