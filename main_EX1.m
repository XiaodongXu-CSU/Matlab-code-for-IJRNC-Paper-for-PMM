%% Performing Unconstrained DMC and PMM Estimate
clear all;close all;clc

    N_out  =  3; N_in =  3;

Maxum_Time   =     20;
    
    N  =  240; P = N-2; M = N/2; W = N;
    CV_order = {...
               'CV1'
               'CV2'
               'CV3'
                };
    MV_order = {...
       'MV1'
       'MV2'
       'MV3'
        };
         CReal_sys   =     cell(N_out, N_in);
                s    =     tf('s');
    CReal_sys{1,1}   =     15.36 * exp(-s)    /(16.7 * s + 1);
    CReal_sys{1,2}   =     0*exp(-0*s)       /(0 * s + 1);
    CReal_sys{1,3}   =     9.24              /(12.7 * s + 1);
    CReal_sys{2,1}   =     0*exp(-0*s)        /(0 * s + 1);
    CReal_sys{2,2}   =     -23.28 * exp(-3*s) /(18.72 * s + 1);
    CReal_sys{2,3}   =     0*exp(-0*s)        /(0 * s + 1);
    CReal_sys{3,1}   =      4.52 * exp(-2*s)  /(6.7 * s + 1);
    CReal_sys{3,2}   =      6.6 * exp(-7*s)  /(8.72 * s + 1);
    CReal_sys{3,3}   =      13.28 * exp(-3*s)  /(18.72 * s + 1);
    
    
    
    K11   =    15.36;   t11  =  16.7;   s11   =   1;
    K13   =    9.24;    t13  =  12.7;  
    K22   =    -23.28;  t22  =  18.72;  s22   =   3;
    K31   =    4.52;    t31  =  6.7;    s31   =   2;
    K32   =    6.6;     t32  =  8.72;   s32   =   7;
    K33   =    13.28;   t33  =  18.72;  s33   =   3;
    
    
    Kh11   =    10.36;   th11  =  19.7;   sh11   =   2;
    Kh13   =    7.24;    th13  =  10.5;  
    Kh22   =    -20.28;  th22  =  20.72;  sh22   =   3;
    Kh31   =    2.52;    th31  =  7.5;    sh31   =   3;
    Kh32   =    4.6;     th32  =  6.72;   sh32   =   5;
    Kh33   =    10.28;   th33  =  20.72;  sh33   =   4;   
    
    
    
    
    
    
    
    dK11   =    K11/Kh11;    dt11  =  t11/th11;    ds11   =   s11/sh11;
    dK13   =    K13/Kh13;    dt13  =  t13/th13;  
    dK22   =    K22/Kh22;    dt22  =  t22/th22;    ds22   =   s22/sh22;
    dK31   =    K31/Kh31;    dt31  =  t31/th31;    ds31   =   s31/sh31;
    dK32   =    K32/Kh32;    dt32  =  t32/th32;    ds32   =   s32/sh32;
    dK33   =    K33/Kh33;    dt33  =  t33/th33;    ds33   =   s33/sh33;
    
    
    

 for i_out  =  1  :  N_out
            for j_in  =  1  : N_in
                                  dsys  =   c2d(CReal_sys{i_out,j_in},1);
                 Dreal_sys{i_out,j_in}  =   dsys; 
                            [S_Real,t]  =   step(dsys,N);
                            temp        =   S_Real';
    Real_fsr(i_out,N*(j_in-1)+1:N*j_in) =   temp(2:end) ;
                          
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
 

         Dnois_sys   =     cell(N_out,N_out); 
                 z   =     tf('z',1);
    Dnois_sys{1,1}   =     (1 - 0.19 * z^(-1) )  /(1 - 0.98 * z^(-1) );
    Dnois_sys{1,2}   =     0 * z;    
    Dnois_sys{1,3}   =     0 * z; 
    Dnois_sys{2,1}   =     0 * z;    
    Dnois_sys{2,2}   =     (0.02 - z^(-1) )  /(1 - 0.9 * z^(-1) ); 
    Dnois_sys{2,3}   =     0 * z; 
    Dnois_sys{3,1}   =     0 * z;    
    Dnois_sys{3,2}   =     (1 - 0.19 * z^(-1) )  /(1 - 0.98 * z^(-1) ); 
    Dnois_sys{3,3}   =     0 * z;    
    
    
    c11    =    0.19;   d11   =   0.98;   c22    =    0.02;   d22   =   0.9;  
    c33    =    0.19;   d33   =   0.98;  
    

    
    
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
 
        CMPC_sys    =     cell(N_out, N_in);
               s    =     tf('s');
    CMPC_sys{1,1}   =     10.36 * exp(-1*s)  /(19.7 * s + 1);
    CMPC_sys{1,2}   =     0*exp(-0*s)        /(0 * s + 1);
    CMPC_sys{1,3}   =     7.24               /(10.5 * s + 1);
    CMPC_sys{2,1}   =     0*exp(-0*s)        /(0 * s + 1);
    CMPC_sys{2,2}   =     -20.28 * exp(-3*s) /(20.72 * s + 1);
    CMPC_sys{2,3}   =     0*exp(-0*s)        /(0 * s + 1);
    CMPC_sys{3,1}   =      2.52 * exp(-2*s)  /(7.5 * s + 1);
    CMPC_sys{3,2}   =      4.6 * exp(-7*s)   /(6.72 * s + 1);
    CMPC_sys{3,3}   =      10.28 * exp(-3*s) /(20.72 * s + 1);      

 for i_out  =  1  :  N_out
            for j_in  =  1  : N_in
                                 dsys  =   c2d(CMPC_sys{i_out,j_in},1);
                 Dmpc_sys{i_out,j_in}  =   dsys; 
                            [S_MPC,t]  =   step(dsys,N);
                            temp       =   S_MPC';
    MPC_fsr(i_out,N*(j_in-1)+1:N*j_in) =   temp(2:end);
            end
 end 
 
            MPC_M_fsr = cell(1,N);
 for mk  =   1  : N
                Stemp = [];
            for j_in  =  1  :  N_in
                Stemp = [Stemp,MPC_fsr(:,N*(j_in-1)+mk)];
            end    
      MPC_M_fsr{mk} = Stemp; 
 end  
 
ZeorM  =  zeros(N_out,N_in);
S_hat  =  [];
for j  =  1: M
            temp  =   [];
           for i  =  1  :  P
              if i <=j-1
              temp = [temp;ZeorM]; 
              else
              temp = [temp;MPC_M_fsr{i-j+1}];     
              end  
           end
              S_hat = [S_hat,temp];       
end
w_q   =   [0.001 0.001 0.001]';
Q     =   diag(reshape(repmat(w_q,1,P),N_out*P,1));
w_r   =   [1 2 1]';
R     =   diag(reshape(repmat(w_r,1,M),N_in*M,1));
KC    =   (S_hat'*Q*S_hat+R)\S_hat'*Q;
 
         S_M_pre  =  []; 
for i = 1 : P
            temp  =  [];
    for j  =  1   :  N-1
        if j+i<=N
            temp  =  [temp,MPC_M_fsr{i+j}]; 
        else
            temp  =  [temp,ZeorM];
        end 
    end
         S_M_pre  =  [S_M_pre;temp];   
end 
 
 T_sim    =    10000;
 
 
past_u    =    zeros(N_in,N);   
past_e    =    zeros(N_out,W+1);
 
 for cnt   =   1  :  T_sim
        [cnt T_sim] 
        
                  e    =   [randn,randn,randn]';  
       ehist(:,cnt)    =   e;
             past_e    =   [e,past_e(:,1:end-1)];
            past_du    =   [-diff(past_u');past_u(:,end)']';        
        
        
              wtemp    =   zeros(N_out,1);
          for k   =   1:W
              wtemp    =   wtemp  +  Nois_M_fir{k} * past_e(:,k); 
          end
          
            yp_temp    =    zeros(N_out,1); ym_temp  =  zeros(N_out,1);
          for  mk  =  1 : N
            yp_temp    =    yp_temp +   Real_M_fsr{mk} *  past_du(:,mk);    
            ym_temp    =    ym_temp +   MPC_M_fsr{mk}  *  past_du(:,mk);  
          end
           yp(:,cnt)   =    yp_temp + wtemp;    
           ym(:,cnt)   =    ym_temp;  
            
             ypredss   =    MPC_M_fsr{N} * flipud(past_u(:,end-P+1:end)')';
              ypred1   =    S_M_pre * reshape(-diff(past_u')',(N-1)*N_in,1) +reshape(ypredss,P*N_out,1);
              ypred2   =    reshape( repmat((yp(:,cnt)-ym(:,cnt))',P,1)',P*N_out,1); 
                 
                  Eo   =    0 - ypred1 - ypred2;
                DelU   =    KC * Eo;
                   U   =    past_u(:,1) + DelU(1:N_in) ; 
        uhist(:,cnt)   =    U;
              past_u   =    [U, past_u(:,1:end-1)];
        
 end
 
 

 
         Lu   =   15; 
 ruu_sample   =   [];
for  pi    =   0   :  Lu  
       temp   =   zeros(N_out,N_out);
     for  cnt =  1:  T_sim -pi 
       temp   =   temp + (uhist(:,cnt)) * (uhist(:,cnt+pi))';
     end
     hat_Ru{pi+1}   =   (T_sim -pi +1)\ temp;
       ruu_sample   =   [ruu_sample,reshape(hat_Ru{pi+1}',N_out*N_out,1)];
end 



Ruu_true  =  zeros(Lu,1);
for  pi  =  0 :  Lu-1
Ruu1_true(pi+1) = dot(uhist(1,1: T_sim-pi), uhist(1,pi+1: T_sim))/(T_sim-pi);  
end
for  pi  =  0 :  Lu-1
Ruu2_true(pi+1) = dot(uhist(2,1: T_sim-pi), uhist(2,pi+1: T_sim))/(T_sim-pi);  
end
for  pi  =  0 :  Lu-1
Ruu3_true(pi+1) = dot(uhist(3,1: T_sim-pi), uhist(3,pi+1: T_sim))/(T_sim-pi);  
end


RUU_tru = [Ruu1_true;Ruu2_true;Ruu3_true];

 
Tru_DEL         =    [ dK11;dt11;dK13;dt13;dK22;dt22;dK31;dt31; ...
                       dK32;dt32;dK33;dt33;c11;d11;c22;d22;c33;d33];

        cost =   PMM_est(Lu,RUU_tru,N,P,W,N_out,N_in,Tru_DEL ,KC,MPC_M_fsr);           
                   
       Ndel0    =    0.9  *  Tru_DEL;  
       cost0    =    PMM_est(Lu,RUU_tru,N,P,W,N_out,N_in,Ndel0 ,KC,MPC_M_fsr);   
          lb    =    0.01 * ones(18,1);
          ub    =    1.5 *  [ones(12,1);0.3*ones(1,1);1*ones(5,1)];
         opts   =    optiset('solver','ipopt','display','iter','maxtime',Maxum_Time*60); 
        fun     =    @(Ndel)  PMM_est(Lu,RUU_tru,N,P,W,N_out,N_in,Ndel,KC,MPC_M_fsr); 
    opt_N_Del   =    opti_fmincon(fun,Ndel0,[],[],[],[],lb,ub,[],opts);
 
       
    Result_plot(opt_N_Del);  
    
    
    
    
    
    
    
    
    




