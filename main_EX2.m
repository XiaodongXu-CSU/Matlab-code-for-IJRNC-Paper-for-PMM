
clear all;close all;clc

    N_out  =  2; N_in =  2;

    

Maxum_Time   =     20;
    

    N  =  240; P = N-2; M = N/2; W = N;
    CV_order = {...
               'CV1'
               'CV2'
                };
    MV_order = {...
       'MV1'
       'MV2'
        };
         load('Model_Error.mat')  
%% formulating continuous real process
         CReal_sys   =     cell(N_out, N_in);
                s    =     tf('s');
    CReal_sys{1,1}   =     15.36 * exp(-s)    /(7.53 * s + 1);
    CReal_sys{1,2}   =     -9.24              /(12.1 * s + 1);
    CReal_sys{2,1}   =     -4.52 * exp(-2*s)  /(6.07 * s + 1);
    CReal_sys{2,2}   =       6.6 * exp(-3*s)  /(8.72 * s + 1);

 for i_out  =  1  :  N_out
            for j_in  =  1  : N_in
                                  dsys  =   c2d(CReal_sys{i_out,j_in},1); 
                 Dreal_sys{i_out,j_in}  =   dsys; 
                            [S_Real,t]  =   step(dsys,N); 
                            temp        =   S_Real';
    Real_fsr(i_out,N*(j_in-1)+1:N*j_in) =   temp(2:end) + Model_Error{i_out,j_in}'; 
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
    Dnois_sys{1,1}   =     (1 - 0.3 * z^(-1) )  /(1 - 0.5 * z^(-1) );
    Dnois_sys{1,2}   =     0 * z;    
    Dnois_sys{2,1}   =     0 * z;    
    Dnois_sys{2,2}   =     (1 - 0.4 * z^(-1) )  /(1 - 0.6 * z^(-1) ); 
    
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
 
 
 %% MPC model
       CMPC_sys   =   cell(N_out, N_in);
                s    =     tf('s');
    CMPC_sys{1,1}   =     10.36 * exp(-2*s)/(9.53*s+1);
    CMPC_sys{1,2}   =     -6.24            /(10.1*s+1);
    CMPC_sys{2,1}   =     -9.52 * exp(-3*s)/(8.07*s+1);
    CMPC_sys{2,2}   =       8.6 * exp(-2*s)/(10.72*s+1);       

 for i_out  =  1  :  N_out
            for j_in  =  1  : N_in
                                 dsys  =   c2d(CMPC_sys{i_out,j_in},1); 
                 Dmpc_sys{i_out,j_in}  =   dsys; 
                            [S_MPC,t]  =   step(dsys,N); 
                            temp       =   S_MPC';
    MPC_fsr(i_out,N*(j_in-1)+1:N*j_in) =   temp(2:end) + Model_Error{i_out,j_in}'; 
            end
 end 
 
[MPC_tf, MPCn_fsr,MPC_ERO] = Model_fit(MPC_fsr, N,N_out,N_in);
 
 
 
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
w_q   =   [0.001 0.001]';
Q     =   diag(reshape(repmat(w_q,1,P),N_out*P,1));
w_r   =   [1 2]';
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
 
 T_sim    =    24000;
 
 
past_u    =    zeros(N_in,N);   
past_e    =    zeros(N_out,W+1); 
 
 for cnt   =   1  :  T_sim
        [cnt T_sim]
    
        
                  e    =   [randn,randn]'; 
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
                  
                  Eo   =    - ypred1 - ypred2;
                DelU   =    KC * Eo;
                   U   =    past_u(:,1) + DelU(1:N_in) ; 
        uhist(:,cnt)   =    U;
              past_u   =    [U, past_u(:,1:end-1)];
        
 end
 
 
         Ly   =   20; 
 ryy_sample   =   [];
for  pi    =   0   :  Ly  
       temp   =   zeros(N_out,N_out);
     for  cnt =  1:  T_sim -pi 
       temp   =   temp + (yp(:,cnt)) * (yp(:,cnt+pi))';
     end
     hat_Ry{pi+1}   =   (T_sim -pi +1)\ temp;
       ryy_sample   =   [ryy_sample,reshape(hat_Ry{pi+1}',N_out*N_out,1)];
end 
 
  Ndel0    =    0.8  *  ones(10,1); 
  NdelR    =    [2*ones(10,1);0.3;0.5;0.4;0.6];


        Ndel    =      2  *  ones(3,1);   
       Ndel0    =    0.9  *  ones(14,1);   
          lb    =    0.01 *  ones(14,1);
          ub    =    [3.5 *  ones(10,1);0.95*ones(4,1)];
         opts        =   optiset('solver','ipopt','display','iter','maxtime',Maxum_Time*60); 
        fun     =   @(Ndel)  PMM_est(Ly,ryy_sample,N,P,W,N_out,N_in,Ndel,KC,MPC_M_fsr,MPC_ERO); 
    opt_N_Del   =   opti_fmincon(fun,Ndel0,[],[],[],[],lb,ub,[],opts);
 
 
FSRplot(N,W,N_out,N_in,opt_N_Del,Real_fsr,MPC_fsr,MPC_ERO);




