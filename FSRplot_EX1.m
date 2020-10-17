function  [] =  FSRplot(N,W,N_out,N_in,Ndel,Real_fsr,MPC_fsr,MPC_ERO)

 %% Formulating Noise model process
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
       Nois_M_fir{mk} = Ntemp; % each Si is a N_out * N_in matrix
 end 




%% formulating continuous real(guess) process
        CReal_g_sys    =     cell(N_out, N_in);
                 s     =     tf('s');
%{
REAL_MV_PARA{1,1} = 10.36*N_del(1)* exp(-s)/(7.53*s+1);
REAL_MV_PARA{1,2} = -3.24* N_del(3) /(12.1*s+1);

REAL_MV_PARA{2,1} = -9.52*N_del(2) * exp(-2*s)/(6.07*s+1);
REAL_MV_PARA{2,2} = 8.6 * N_del(4)* exp(-3*s)/(8.72*s+1);
%}
 
    CReal_g_sys{1,1}   =     11.77   *  Ndel(1)  * exp(-0.189* s)            /(10.32 * Ndel(5) *s + 1);
    CReal_g_sys{1,2}   =     - 4.77  *  Ndel(2)  * exp(-3.87 * Ndel(9)  * s) /(8.86  * Ndel(6) *s + 1);
    CReal_g_sys{2,1}   =     - 8.12  *  Ndel(3)*0.6  * exp(-4.63 * Ndel(10) * s) /(7.31  * Ndel(7) *s + 1);
    CReal_g_sys{2,2}   =     10.04   *  Ndel(4)*0.75                              /(10.34 * Ndel(8) *s + 1);      
       
       
       
 for i_out  =  1  :  N_out
            for j_in  =  1  : N_in
                                  dsys  =   c2d(CReal_g_sys{i_out,j_in},1); % convert continous to discrete 
                 Dreal_sys{i_out,j_in}  =   dsys; 
                            [S_Real,t]  =   step(dsys,N); % generating N+1 FSR coefficients  
                            temp        =   S_Real';
    Real_fsrup(i_out,N*(j_in-1)+1:N*j_in) =   temp(2:end) + MPC_ERO{i_out,j_in}';  % generating step response coefficients
            end
 end
  
   

figure()
for i_out = 1: N_out
        for j_in =  1: N_in
        subplot(N_out,N_in,(i_out-1)*N_out+j_in);
        plot(Real_fsr(i_out,N*(j_in-1)+1:N*j_in),'linewidth',1.5)
        hold on 
        plot(MPC_fsr(i_out,N*(j_in-1)+1:N*j_in),'linewidth',1.5)
        plot(Real_fsrup(i_out,N*(j_in-1)+1:N*j_in),'linewidth',1.5,'color', 'k')
        xlabel('t','FontSize',12)
        legend('Plant process','Original MPC','Updated MPC')
        %ylabel('$y_r(t)$','interpreter','latex','FontSize',16)
        grid  
        title(['FSR',num2str(i_out),num2str(j_in)])
        axis([0 240 -inf inf])
%set(gca,'fontsize',20)
set(gca,'FontName','Times New Roman','FontSize',12)
        end
end

figure()
for i_out = 1: N_out
        for j_in =  1: N_in
        subplot(N_out,N_in,(i_out-1)*N_out+j_in);
        plot(Real_fsr(i_out,N*(j_in-1)+1:N*j_in),'linewidth',1.5)
        hold on 
        plot(MPC_fsr(i_out,N*(j_in-1)+1:N*j_in),'linewidth',1.5)
      %  plot(Real_fsrup(i_out,N*(j_in-1)+1:N*j_in),'linewidth',1.5,'color', 'k')
        xlabel('t','FontSize',12)
        legend('Plant process','Original MPC')
        %ylabel('$y_r(t)$','interpreter','latex','FontSize',16)
        grid  
        title(['FSR',num2str(i_out),num2str(j_in)])
        axis([0 240 -inf inf])
%set(gca,'fontsize',20)
set(gca,'FontName','Times New Roman','FontSize',12)
        end
end


