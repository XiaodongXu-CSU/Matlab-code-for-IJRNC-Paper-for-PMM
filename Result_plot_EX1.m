function  []  =    Result_plot(Ndel) 
N_out=3; N_in=3;   N  =  240; P = N-2; M = N/2; W = N;
         CReal_sys   =     cell(N_out, N_in);
                s    =     tf('s');
    CReal_sys{1,1}   =     15.36 * exp(-s)    /(16.7 * s + 1);
    CReal_sys{1,2}   =     0*exp(-0*s)        /(0 * s + 1);
    CReal_sys{1,3}   =     9.24               /(12.7 * s + 1);
    CReal_sys{2,1}   =     0*exp(-0*s)        /(0 * s + 1);
    CReal_sys{2,2}   =     -23.28 * exp(-3*s) /(18.72 * s + 1);
    CReal_sys{2,3}   =     0*exp(-0*s)        /(0 * s + 1);
    CReal_sys{3,1}   =      4.52 * exp(-2*s)  /(6.7 * s + 1);
    CReal_sys{3,2}   =      6.6 * exp(-7*s)   /(8.72 * s + 1);
    CReal_sys{3,3}   =      13.28 * exp(-3*s) /(18.72 * s + 1);
    

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
 
 %% Formulating Noise model process
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
     CMPC_sys{1,1}   =     10.36 * exp(-2*s)  /(19.7 * s + 1);
    CMPC_sys{1,2}   =     0*exp(-0*s)        /(0 * s + 1);
    CMPC_sys{1,3}   =     7.24               /(10.5 * s + 1);
    CMPC_sys{2,1}   =     0*exp(-0*s)        /(0 * s + 1);
    CMPC_sys{2,2}   =     -20.28 * exp(-3*s) /(20.72 * s + 1);
    CMPC_sys{2,3}   =     0*exp(-0*s)        /(0 * s + 1);
    CMPC_sys{3,1}   =      2.52 * exp(-3*s)  /(7.5 * s + 1);
    CMPC_sys{3,2}   =      4.6 * exp(-5*s)   /(6.72 * s + 1);
    CMPC_sys{3,3}   =      10.28 * exp(-4*s) /(20.72 * s + 1);         

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
 

 c11  =  Ndel(13); d11 = Ndel(14);c22 = Ndel(15); d22 = Ndel(16);c33 = Ndel(17);d33 = Ndel(18);
 dK11=  Ndel(1);  dt11= Ndel(2); dK13=  Ndel(3);dt13=  Ndel(4);dK22=  Ndel(5);
 dt22=  Ndel(6);  dK31=  Ndel(7);dt31=  Ndel(8); 
 dK32=  Ndel(9);dt32=  Ndel(10);dK33=  Ndel(11);dt33=  Ndel(12);


%% formulating continuous real(guess) process
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
                          
                            [S_Real,t]  =   step(dsys,N); 
                            temp        =   S_Real';
                Real_g_fsr(i_out,N*(j_in-1)+1:N*j_in) =   temp(2:end);
                                           
            end
 end
 
           Real_gM_fsr = cell(1,N);
 for mk  =   1  : N
                Stemp = [];
            for j_in  =  1  :  N_in
               Stemp = [Stemp,Real_g_fsr(:,N*(j_in-1)+mk)];
            end    
      Real_gM_fsr{mk} = Stemp; 
 end 


 T_sim    =    10000;
 
 

past_u    =    zeros(N_in,N);   
past_e    =    zeros(N_out,W+1); 
 
 for cnt   =   1  :  T_sim
              [cnt T_sim]
          if cnt  <  T_sim/2         
             yr(:,cnt)  =  [-10 10 20]';
          else 
            yr(:,cnt)  =   [10 -20 10]';     
          end  
       
                  e    =   [randn,randn,randn]';  
       ehist(:,cnt)    =   e;
             past_e    =   [e,past_e(:,1:end-1)]; 
            past_du    =   [-diff(past_u');past_u(:,end)']';        
        
        
              wtemp    =   zeros(N_out,1);
          for k   =   1:W
              wtemp    =   wtemp  +  Nois_M_fir{k} * past_e(:,k); 
          end
          
            yp_temp    =    zeros(N_out,1); ym_temp  =  zeros(N_out,1); ygm_temp  =  zeros(N_out,1);
          for  mk  =  1 : N
            yp_temp    =    yp_temp +   Real_M_fsr{mk} *  past_du(:,mk);  
            ygm_temp   =   ygm_temp +   Real_gM_fsr{mk} *  past_du(:,mk); 
            ym_temp    =    ym_temp +   MPC_M_fsr{mk}  *  past_du(:,mk);  
          end
           yp(:,cnt)   =    yp_temp + wtemp; 
          ygm(:,cnt)   =    ygm_temp;
           ym(:,cnt)   =    ym_temp;  
            
             ypredss   =    MPC_M_fsr{N} * flipud(past_u(:,end-P+1:end)')';
              ypred1   =    S_M_pre * reshape(-diff(past_u')',(N-1)*N_in,1) +reshape(ypredss,P*N_out,1);
              ypred2   =    reshape( repmat((yp(:,cnt)-ym(:,cnt))',P,1)',P*N_out,1); 
                  Yr   =    reshape((repmat(yr(:,cnt)',P,1))',P*N_out,1);
                  Eo   =    Yr - ypred1 - ypred2;
                DelU   =    KC * Eo;
                   U   =    past_u(:,1) + DelU(1:N_in) ; 
        uhist(:,cnt)   =    U;
              past_u   =    [U, past_u(:,1:end-1)];
        
 end
 
 
 
  figure()
for i_out = 1: N_out
        for j_in =  1: N_in
        subplot(N_out,N_in,(i_out-1)*N_out+j_in);
        plot( Real_fsr(i_out,N*(j_in-1)+1:N*j_in),'linewidth',1.5,'color', 'b')
        hold on 
        plot( Real_g_fsr(i_out,N*(j_in-1)+1:N*j_in),'--r','linewidth',1.5)
        plot( MPC_fsr(i_out,N*(j_in-1)+1:N*j_in),'--k','linewidth',1.5)
        xlabel('t','FontSize',12)
        
        if i_out ==1 && j_in == 1
        legend('Real model','Upd DMC model','DMC model');
        end
        
      
        grid  
        title(['FSR',num2str(i_out),num2str(j_in)])
        axis([0 N -inf inf])
set(gca,'FontName','Times New Roman','FontSize',12)
        end
end   
 
 
 
 
 
 
 
 
 
 % figure_plot
 
  figure()
 for i_out  =  1  :  N_out
 subplot(1,N_out,i_out)
 plot(ym(i_out,:),'-b');
 hold on
 plot(yp(i_out,:),'-k');
  plot(ygm(i_out,:),'-r');
 legend('ym','yp','ygm')
 end