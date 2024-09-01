%-- .- -.. . / -... -.-- / -.-- --- --. . ... ....
%Date:19_11_2023
clc;clear all;close all;
N=50;    %number of time steps;
C=[0 1];
Q=[0.04 0;
    0 0.04];
R=0.09;
P0=1.5;
V0=0;
%initialise covariance matrix
p=[0 0;
    0 0];
y=[];
x=[V0,P0]';
w1=Q(1,1)*randn(N,1);
w2=Q(2,2)*randn(N,1);
v=R*randn(N, 1);
%input generation
x1=x(1);
x2=x(2);
%x1_vec=[x1];
%x2_vec=[x2];
for i=1:N
    x1new=(x1*x2)/(1+x1*x1)+w1(i);
    x2new=x2/(1+x1*x1)+w2(i);
    x1=x1new;
    x2=x2new;
%    x1_vec=[x1_vec;x1];
%    x2_vec=[x2_vec;x2];
    y=[y, x2+v(i)];
end
%y(k) to be arranged for each step
t=1:N;
plot(t, y, 'DisplayName',"Generated data");
hold on
clear x1 x2 x1new x2new w1 w2 v;
xt_initial={[0,0]', [20,-20]', [-5,-5]'};
    for f=1:length(xt_initial)
        xt=xt_initial{f};
        %Estimation_begin
        for k=2:N
            xt_ap=[(xt(1,k-1)*xt(2,k-1))/(1+xt(1,k-1)*xt(1,k-1)) xt(2,k-1)/(1+xt(1,k-1)*xt(1,k-1))]'; %state estimate propagation
            A=[((1+3*xt_ap(1)*xt_ap(1))*xt_ap(2))/((1+xt_ap(1)*xt_ap(1))^2) xt_ap(1)/(1+xt_ap(1)*xt_ap(1));
                -2*xt_ap(1)*xt_ap(2)/((1+xt_ap(1)*xt_ap(1))^2)           1/(1+xt_ap(1)*xt_ap(1))];
            %error c0variance propagation
            p_bar=A*p*A'+Q;
        
            %find kalman gain
            K=p_bar*C'*inv(C*p_bar*C'+R);
        
            %update state estimate
            xt(:,k)=xt_ap+K*[y(k)-C*xt_ap];
        
            %update covariance matrix
            p=(eye(2,2)-K*C)*p_bar; 
            p_chol=chol(p);
            p=p_chol*p_chol';
        end
        t=1:N;
        %print_xt=num2str(xt_initial{f})
        plot(t,xt(2,:), 'DisplayName', ['x initial = ', num2str(xt(1,1)),' , ' num2str(xt(2,1))])
        legend; %("Generated data","Predicted data")
        % Clearing all variables except f and xt_initial at the end
        vars = who;
        clearvars('-except', 'f', 'xt_initial', vars{:});
    end