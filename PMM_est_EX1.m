function  [cost,ryy] =  PMM_est(Lu,ruu_sample,N,P,W,N_out,N_in,Ndel,KC,MPC_M_fsr)


 c11  =  Ndel(13); d11 = Ndel(14);c22 = Ndel(15); d22 = Ndel(16);c33 = Ndel(17);d33 = Ndel(18);
 dK11=  Ndel(1);  dt11= Ndel(2); dK13=  Ndel(3);dt13=  Ndel(4);dK22=  Ndel(5);
 dt22=  Ndel(6);  dK31=  Ndel(7);dt31=  Ndel(8); 
 dK32=  Ndel(9);dt32=  Ndel(10);dK33=  Ndel(11);dt33=  Ndel(12);


         Dnois_sys   =     cell(N_out,N_out); 
                 z   =     tf('z',1);
    Dnois_sys{1,1}   =     (1 - c11 * z^(-1) )  /(1 - d11 * z^(-1) );
    Dnois_sys{1,2}   =     0 * z;    
    Dnois_sys{1,3}   =     0 * z; 
    Dnois_sys{2,1}   =     0 * z;    
    Dnois_sys{2,2}   =     (c22 - z^(-1) )      /(1 - d22 * z^(-1) ); 
    Dnois_sys{2,3}   =     0 * z; 
    Dnois_sys{3,1}   =     0 * z;    
    Dnois_sys{3,2}   =     (1 - c33 * z^(-1) )  /(1 - d33 * z^(-1) ); 
    Dnois_sys{3,3}   =     0 * z;   
    
   
    
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
 
    CReal_g_sys{1,1}   =      10.36 *dK11 * exp(- s)  /(19.7*dt11 * s + 1);
    CReal_g_sys{1,2}   =      0*exp(-0*s)        /(0 * s + 1);
    CReal_g_sys{1,3}   =      7.24  *dK13        /(10.5 *dt13* s + 1);
    CReal_g_sys{2,1}   =      0*exp(-0*s)        /(0 * s + 1);
    CReal_g_sys{2,2}   =      -20.28 *dK22* exp(-3*s) /(20.72*dt22 * s + 1);
    CReal_g_sys{2,3}   =      0*exp(-0*s)        /(0 * s + 1);
    CReal_g_sys{3,1}   =      2.52 *dK31* exp(-2*s)   /(7.5 *dt31* s + 1);
    CReal_g_sys{3,2}   =      4.6 *dK32* exp(-7*s)    /(6.72 *dt32* s + 1);
    CReal_g_sys{3,3}   =      10.28 *dK33* exp(-3*s)  /(20.72 *dt33* s + 1);                
                      
                
       
 for i_out  =  1  :  N_out
            for j_in  =  1  : N_in
                                  dsys  =   c2d(CReal_g_sys{i_out,j_in},1); 
                 Dreal_sys{i_out,j_in}  =   dsys; 
                            [S_Real,t]  =   step(dsys,N);
                            temp        =   S_Real';
    Real_fsr(i_out,N*(j_in-1)+1:N*j_in) =   temp(2:end);
                                          
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

 C_temp2 = zeros(N_in,N_in);
    for k = 1:P
       C_temp2 = C_temp2 + KC(1:N_in,N_out*(k-1)+1:N_out*k)* (Del_M_fir{N}); 
    end
C{N} = -C_temp2;


for i = 0 :  N
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
 
ruu1   =  [];ruu2   =  [];ruu3   =  [];
for pi = 0 : Lu-1
     Ru1temp = 0;Ru2temp = 0;Ru3temp = 0;
     for i = 0 : W - pi-1
        Ru1temp = Ru1temp +  Gama{i+1}(1,:) * (Gama{i+pi+1}(1,:))'; 
         Ru2temp = Ru2temp +  Gama{i+1}(2,:) * (Gama{i+pi+1}(2,:))'; 
          Ru3temp = Ru3temp +  Gama{i+1}(3,:) * (Gama{i+pi+1}(3,:))'; 
     end
        ruu1 = [ruu1,Ru1temp];ruu2 = [ruu2,Ru2temp];ruu3 = [ruu3,Ru3temp];
 end





druu1 = ruu1 - ruu_sample(1,:);
druu2 = ruu2 - ruu_sample(2,:);
druu3 = ruu3 - ruu_sample(3,:);
cost = dot(druu1,druu1)*1e2 + dot(druu2,druu2)*1e4 + dot(druu2,druu2)*1e4;





