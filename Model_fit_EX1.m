function [MPC_tf, MPCn_fsr,MPC_ERO] = Model_fit(MPC_fsr, N,N_out,N_in)

    type     =   'P1D';
 init_sys    =   idproc(type);  

   MPC_tf    =   cell(N_out,N_in);   
  MPCn_fsr    =   cell(N_out,N_in); 
for i_out  =   1  :  N_out
    for  j_in  =  1  : N_in
         U    =    ones(N+1,1);
      data    =    iddata([0;MPC_fsr(i_out,N*(j_in-1)+1:N*j_in)'],U,1);
      Gsys    =    pem(data,init_sys);
      MPC_tf{i_out,j_in}  =  idtf(Gsys);
      mpctf_temp   =  c2d(MPC_tf{i_out,j_in},1);
      mpcfsr_temp  =  step(mpctf_temp,N);
      MPCn_fsr{i_out,j_in}  =  mpcfsr_temp(2:end);
      MPC_ERO{i_out,j_in}   =  MPC_fsr(i_out,N*(j_in-1)+1:N*j_in)' -  MPCn_fsr{i_out,j_in}; 
    end
end



figure()
for i_out = 1: N_out
        for j_in =  1: N_in
        subplot(N_out,N_in,(i_out-1)*N_out+j_in);
        plot(MPC_fsr(i_out,N*(j_in-1)+1:N*j_in),'linewidth',1.5)
        hold on 
        plot(MPCn_fsr{i_out,j_in},'linewidth',1.5,'color', 'k')
        plot(MPC_ERO{i_out,j_in},'linewidth',1.5,'color', 'r')
        xlabel('t','FontSize',12)
        legend('Original MPC','Fit MPC','Model Error')
        %ylabel('$y_r(t)$','interpreter','latex','FontSize',16)
        grid  
        title(['FSR',num2str(i_out),num2str(j_in)])
        axis([0 240 -inf inf])
%set(gca,'fontsize',20)
set(gca,'FontName','Times New Roman','FontSize',12)
        end
end