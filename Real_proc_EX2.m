

clear all;close all;clc
N_out  =  2; N_in =  2; % F_in = 2;


N  =  240; P = N-2; M = N/2; W = N;
 


 REAL_MV_PARA = cell(N_out, N_in);

       s      =  tf('s');
  REAL_MV_PARA = cell(N_out, N_in);
       s      =  tf('s');
 REAL_MV_PARA{1,1} = 15.36 * exp(-s)/(7.53*s+1);
 REAL_MV_PARA{1,2} = -9.24 /(12.1*s+1);
 
 REAL_MV_PARA{2,1} = -4.52 * exp(-2*s)/(6.07*s+1);
 REAL_MV_PARA{2,2} = 6.6 * exp(-3*s)/(8.72*s+1);
 
 REAL_MV_FSR = cell(N_out, N_in);
 
for i_out  =  1  :  N_out
    for j_in  =  1  : N_in
                           sys  =  c2d(REAL_MV_PARA{i_out,j_in},1);
                    [S_Pred,t]  =  step(sys,N); 
                    [H_Real,t] =  impulse(sys,N);
                    temp        =   S_Pred';
                    temh        =   H_Real';
        REAL_MV_FSR{i_out,j_in} =  temp(2:end); 
        REAL_MV_FIR{i_out,j_in} =  temh(2:end);
        Model_Error{i_out,j_in} =  3*rand(N,1);
        REAL_IR_FSR{i_out,j_in} =  REAL_MV_FSR{i_out,j_in} +Model_Error{i_out,j_in}'; 
    end
end

save('Model_Error.mat','Model_Error');


figure()
for i_out = 1: N_out
        for j_in =  1: N_in
        subplot(N_out,N_in,(i_out-1)*N_out+j_in);
        plot(REAL_IR_FSR{i_out,j_in},'linewidth',1.5)
        xlabel('t','FontSize',12)
      
        grid  
        title(['FSR',num2str(i_out),num2str(j_in)])
        axis([0 240 -inf inf])

set(gca,'FontName','Times New Roman','FontSize',12)
        end
end


N_out  =  2; N_in =  2;


N  =  240; P = N-2; M = N/2; W = N;

 

 REAL_MV_PARA = cell(N_out, N_in);

       s      =  tf('s');
  REAL_MV_PARA = cell(N_out, N_in);
       s      =  tf('s');
 REAL_MV_PARA{1,1} = 10.36 * exp(-2*s)/(9.53*s+1);
 REAL_MV_PARA{1,2} = -6.24 /(10.1*s+1);
 
 REAL_MV_PARA{2,1} = -9.52 * exp(-3*s)/(8.07*s+1);
 REAL_MV_PARA{2,2} = 8.6 * exp(-2*s)/(10.72*s+1);
 
 REAL_MV_FSR = cell(N_out, N_in);
 
for i_out  =  1  :  N_out
    for j_in  =  1  : N_in
                           sys  =  c2d(REAL_MV_PARA{i_out,j_in},1);
                    [S_Pred,t]  =  step(sys,N); 
                    [H_Real,t] =  impulse(sys,N);
                    temp        =   S_Pred';
                    temh        =   H_Real';
        REAL_MV_FSR{i_out,j_in} =  temp(2:end);  
        REAL_MV_FIR{i_out,j_in} =  temh(2:end);
        MPC_IR_FSR{i_out,j_in} =  REAL_MV_FSR{i_out,j_in} +Model_Error{i_out,j_in}'; 
    end
end




figure()
for i_out = 1: N_out
        for j_in =  1: N_in
        subplot(N_out,N_in,(i_out-1)*N_out+j_in);
        plot(REAL_IR_FSR{i_out,j_in},'linewidth',1.5)
        hold on
         plot(MPC_IR_FSR{i_out,j_in},'linewidth',1.5)
         legend('Real FSR','MPC FSR')
        xlabel('t','FontSize',12)
      
        grid  
        title(['FSR',num2str(i_out),num2str(j_in)])
        axis([0 240 -inf inf])
set(gca,'FontName','Times New Roman','FontSize',12)
        end
end








%{
plot(T_scale,ys,'linewidth',1.5,'color','b','linestyle','--')
axis([0 12 -6 6])
ylabel('$y_r(t)$, $y_o(t)$ and $y_s(t)$','interpreter','latex','FontSize',16)
xlabel('$t$','interpreter','latex','FontSize',16)
h=legend('$y_r(t)$','$y_o(t)$','$y_s(t)$')
set(h,'interpreter','latex','FontSize',16)
grid
%}


% matrix type FSR models S_i are matrixces-- this is important.
A_MV_FSR = cell(1,N);   
A_MV_FIR = cell(1,N);
for k   =   1 : N
    for i_out = 1: N_out
        for j_in = 1:N_in
            A_MV_FSR{k}(i_out,j_in) = REAL_MV_FSR{i_out,j_in}(k); 
            A_MV_FIR{k}(i_out,j_in) = REAL_MV_FIR{i_out,j_in}(k); 
        end
    end
end

%% Real Noise Models
Nc{1,1} = [1 -0.19];Nc{1,2} = [0.02 1];Nc{1,3} = [1 -0.19];Nc{1,4} = [0.02 1];
Nd{1,1} = [1 -0.98];Nd{1,2} = [1 -0.9];Nd{1,3} = [1 -0.98];Nd{1,4} = [1 -0.9];
Beta    =  cell(1,N_out);
for i_out =  1  :  N_out
            num  = Nc{1,i_out};
            den  = Nd{1,i_out};
           Dsys  =  tf(num,den,1);
         [H_e,t] =  impulse(Dsys,W);
  BH{1,i_out}  =  H_e';
end

for k  =  1 :  W
    Btemp = [];
    for i_out = 1: N_out
        Btemp = [Btemp,BH{i_out}(k)];
    end
    Beta{k} = diag(Btemp);    
end
