clear;
clc;
%% Load of Linearized Model

load('linsys2.mat')
load('Linearized_Plant_by_MPC.mat')
plant_C = plant;
%%  Wave Disturbance Model

h_13 = 3; %(m) %wave significant height

%Wave disturbance parameters
a = 2*10^(-1); 
b = -5*10^(-1);
c = 2*10^(3);
d = 2*10^(3);

a = a/2;
b = b/2;
c = c/2;
d = d/2;

B_ = 3.11/(h_13^2);
K = 1.11; %Maximum wave peak
wm = (0.8*B_)/4;
syms s;
s = tf('s');

%Disturbance transfer function
G4 = (K*(s/wm)^2/((s/wm)^2 + s/wm + 1)^2);

%% First MPC Controller for step reproduction for linearized system
%% create MPC controller object with sample time
mpc1 = mpc(plant_C, 0.1);
%% specify prediction horizon
mpc1.PredictionHorizon = 15;
%% specify control horizon
mpc1.ControlHorizon = 1;
%% specify nominal values for inputs and outputs
mpc1.Model.Nominal.U = [0;0;0;0;0;0;-2500;100000];
mpc1.Model.Nominal.Y = [0;0;0;0;0;0];
%% specify overall adjustment factor applied to weights
beta = 4.2207;
%% specify weights
mpc1.Weights.MV = [0 0 0 0 0 0]*beta;
mpc1.Weights.MVRate = [0.1 0.1 0.1 0.1 0.1 0.1]/beta;
mpc1.Weights.OV = [100 100 100 100 100 100]*beta;
mpc1.Weights.ECR = 100000;
%% specify simulation options
options = mpcsimopt();
options.RefLookAhead = 'off';
options.MDLookAhead = 'off';
options.Constraints = 'on';
options.OpenLoop = 'off';
%% run simulation
%%sim(mpc1, 101, mpc1_RefSignal, mpc1_MDSignal, options);

%% Decoupled PID
G = linsys1;
figure(1);
stepplot(G);
% Code inspired by the MIMO controller mathworks tutorial
G.InputName = {'pu1','pu2', 'pu3', 'pu4', 'pu5', 'pu6'};
G.OutputName = 'y';

%The use of the decoupler gain is suggested by mathwork doc
D = tunableGain('Decoupler',eye(6));
D.u = 'e';
D.y = {'u1','u2', 'u3', 'u4', 'u5', 'u6'};

%tunable tool, tunable doesn't find the parameter, but creates a parametric
%controller where the parameter can be computed through the systune()
%function
C_1 = tunablePID('C_1','pid'); 
C_1.u = 'u1'; C_1.y = 'pu1';

C_2 = tunablePID('C_2','pid');  
C_2.u = 'u2'; C_2.y = 'pu2';

C_3 = tunablePID('C_3','pid');  
C_3.u = 'u3'; C_3.y = 'pu3';

C_4 = tunablePID('C_4','pid');  
C_4.u = 'u4'; C_4.y = 'pu4';

C_5 = tunablePID('C_5','pid');  
C_5.u = 'u5'; C_5.y = 'pu5';

C_6 = tunablePID('C_6','pid');  
C_6.u = 'u6'; C_6.y = 'pu6';

Sum = sumblk('e = r - y',6);

%Parametric closed loop, the parameters are:
%Gains of decoupler matrix
%PI controllers parameters (i.e. Kp, Ki)
CL_Parametric = connect(G,D,C_1,C_2,C_3,C_4,C_5,C_6,Sum,'r','y');

%Setting step response requiremnts, to avoid problem are used the same
%mathwork doc parameters
Rtrack = TuningGoal.Tracking('r','y',10);

%Systuning, this function finds the best PI parameters to achieve the
%tuning goal
[CL2,fSoft] = systune(CL_Parametric,Rtrack);

%Plotting closed loop plot, if stable and solves tracking it is OK
figure(2);
stepplot(CL2(1:3,1:6),10);
title("Step plot of the PID controller for the position outputs");
legend("Step");
xlabel('$t$','Interpreter','latex','FontSize',12);
grid;
figure(21);
stepplot(CL2(4:6,1:6),10);
title("Step plot of the PID controller for the orientation outputs");
legend("Step");
xlabel('$t$','Interpreter','latex','FontSize',12);
grid;

%To evaluate the performances are now plotted the sensitivity functions
%Ce is the extended controller
Ce = [tf(CL2.Blocks.C_1) 0 0 0 0 0; 0 tf(CL2.Blocks.C_2) 0 0 0 0; 0 0 tf(CL2.Blocks.C_3) 0 0 0; 0 0 0 tf(CL2.Blocks.C_4) 0 0; 0 0 0 0 tf(CL2.Blocks.C_5) 0; 0 0 0 0 0 tf(CL2.Blocks.C_6)];

%For sake of semplicity, the controller is evaluated through the loopsens
%function
loops = loopsens(G,Ce*tf(CL2.Blocks.Decoupler)); 
K_PID = Ce*tf(CL2.Blocks.Decoupler);

%Plot of sensitivity functions, they defines the perfomances to achieve
%with the mixed sensitivity controller, in this case are evaluated output
%function So,To
figure(3);
bodemag(loops.So(1:3,1:6),'r',loops.To(1:3,1:6),'b',loops.Lo(1:3,1:6),'g')
legend('Sensitivity','Complementary Sensitivity','Loop Transfer');
title('Performances of the PID controller for position outputs');
grid;

figure(31);
bodemag(loops.So(4:6,1:6),'r',loops.To(4:6,1:6),'b',loops.Lo(4:6,1:6),'g')
legend('Sensitivity','Complementary Sensitivity','Loop Transfer');
title('Performances of the PID controller for orientation outputs');
grid;

%% Mixed sensitivity problem
G = linsys1;
%Weigths definition, they are choesen to met the performances of the
%decentralized PI
wsb = 10^(2); As = 10^(-4); Ms = 5; %Ms = 1.77;
wt = 1000; At = 10^(-4); Mt = 5; %Mt = 1.12;

w1=((s/Ms)+wsb)/(s+(wsb*As));
w2 = tf(0.001);
w3 = ((s+(wt/Mt)))/(s*At+(wt));

%The weigth are chosen as diagonal 2x2 matrices
W1 = eye(6)*w1; 
W2 = eye(6)*w2;
W3 = eye(6)*w3;

%Weigthd s connections
W1.u = {'e'}; W1.y = {'z1'}; 
W2.u = {'u'}; W2.y = {'z2'};
W3.u = {'y'}; W3.y = {'z3'};

G.u = {'u'};
G.y = 'y';

Sum = sumblk('e = r - y',6);

%Generalized plant model
P = connect(G,W1,W2,W3,Sum,{'r','u'},{'z1','z2','z3','e'});

%Solution of mixed sensitivity
[K_mix,CL2,gamma] = hinfsyn(P,6,6);

%Taking sensitivity functions of the mixed sensitivity controller
loops2 = loopsens(G,K_mix);
CLH_inf = feedback(G*K_mix,eye(6));
figure(4);
stepplot(CLH_inf(1:3,1:6),10);
title("Step plot of the H infinity controller the position outputs");
legend("Step");
xlabel('$t$','Interpreter','latex','FontSize',12);
grid;
figure(41);
stepplot(CLH_inf(4:6,1:6),10);
title("Step plot of the H infinity controller for orientation outputs");
legend("Step");
xlabel('$t$','Interpreter','latex','FontSize',12);
grid;

%Plot of mixed sensitivity performances
figure(5);
bodemag(loops2.So(1:3,1:6),'r',loops2.To(1:3,1:6),'b',loops2.Lo(1:3,1:6),'g');
legend('Sensitivity','Complementary Sensitivity','Loop Transfer');
title('Performances of the H infinity controller for position outputs');
grid;

figure(51);
bodemag(loops2.So(4:6,1:6),'r',loops2.To(4:6,1:6),'b',loops2.Lo(4:6,1:6),'g');
legend('Sensitivity','Complementary Sensitivity','Loop Transfer');
title('Performances of the H infinity controller for orientation outputs');
grid;

%Plot of comparison between mixed sensitivity and PI
figure(6);
bodemag(loops2.So,'r',loops2.To,'b',loops.So,'r--',loops.To,'b--');
legend('S K_mix','T K_mix','S PID','T PID');
title("Comparision between $K_{mix}$ and  PID performances",'interpreter','latex');
ylim([-100 20]);

grid;

%% LQR

G=linsys1;

%We create a vector representing the 12 states in order to evaluate the step response
%in the Simulink "LQRFINALEMODEL.slx" file
x_desired=[1,1,1,1,1,1,0,0,0,0,0,0];

%% Weights and construction of the controller
Q=[10^8 zeros(1,11);
    0 10^8 zeros(1,10);
    zeros(1,2) 10^7 zeros(1,9);
    zeros(1,3) 5*10^10 zeros(1,8);
    zeros(1,4) 10^10 zeros(1,7);
    zeros(1,5) 10^8 zeros(1,6);
    zeros(1,6) 10^6 zeros(1,5);
    zeros(1,7) 1 zeros(1,4);
    zeros(1,8) 1 zeros(1,3);
    zeros(1,9) 1 zeros(1,2);
    zeros(1,10) 1 0;
    zeros(1,11) 1];
R=[5 zeros(1,5);
    0 5 zeros(1,4);
    zeros(1,2) 5 zeros(1,3);
    zeros(1,3) 5 zeros(1,2);
    zeros(1,4) 5 0
    zeros(1,5) 5];

[K_LQR,S,e]=lqr(G,Q,R);


%% Trajectory building

sample_time = 20;

%% Diving at a depth of (-18,-18,-18)  
x_trajectory = [0 -1 -2 -3 -4 -5 -6 -7 -8 -9 -10 -11 -12 -13 -14 -15 -16 -17 -18];
y_trajectory = [0 -1 -2 -3 -4 -5 -6 -7 -8 -9 -10 -11 -12 -13 -14 -15 -16 -17 -18];
z_trajectory = [0 -1 -2 -3 -4 -5 -6 -7 -8 -9 -10 -11 -12 -13 -14 -15 -16 -17 -18];

roll_trajectory = zeros(1,19);
pich_trajectory = zeros(1,19);
yaw_trajectory = zeros(1,19);

% %% Arriving at depth of (-20) avoiding obstacles during the trip  
% x_trajectory = [0 2 4 5 6 7 8 9 10 11 12 11 11 11 10 5 2 1 0];
% y_trajectory = [0 2 4 4 4 4 4 5 6 7 8 7 6 5 4 3 2 1 0];
% z_trajectory = [0 -2 -6 -10 -14 -18 -20 -20 -20 -20 -20 -20 -20 -20 -20 -20 -20 -20 -20]; 
% 
% roll_trajectory = [0 2 8 6 2 0 0 0 0 0 zeros(1,10)];
% pich_trajectory = [zeros(1,5) 30 45 30 zeros(1,5) 30 45 30 zeros(1,3)];
% yaw_trajectory = zeros(1,19);

% %% Discending and Ascending
% x_trajectory = [0 1 2 3 4 5 6 7 8 9 10 11 12 12 13 14 14 16 16];
% y_trajectory = [zeros(1,10) 10 20 30 30 25 30 35 40 40];
% z_trajectory = [0 -2 -3 -4 -5 -6 -7 -8 -9 -10 -9 -8 -7 -6 -5 -4 -3 -2 0];
% 
% roll_trajectory = [zeros(1,8) 10 20 30 20 10 zeros(1,6)];
% pich_trajectory = [0 10 20 30 zeros(1,6) 10 45 20 10 zeros(1,5)];
% yaw_trajectory = [zeros(1,4) 10 45 60 20 10 zeros(1,10)];

% %% Fossa delle Marianne example
% x_trajectory = [zeros(1,5) -2 -3 -2 -3 -6 -5 -3 -2 -1 0 -1 -2 -5 -6];
% y_trajectory = [0 -1 -3 -5 -7 -9 -10 -12 -14 -17 -19 -23 -25 -29 -30 -35 -39 -45 -90];
% z_trajectory = [0 -200 -400 -600 -1000 -2000 -3000 -4000 -5000 -6000 -7000 -8000 -8500 -8600 -9000 -9500 -10000 -10500 -11000];
% 
% roll_trajectory = zeros(1,19);
% pich_trajectory = [zeros(1,8) -3 -4 -5 -2 10 -2 3 2 3 5 0];
% yaw_trajectory = zeros(1,19);

%% Simulink

sim('Final_Version.slx')


%% MPC plot
out_MPC = ans.yout.getElement('out_MPC');

t = out_MPC.Values.Time;
figure(9);
out_MPC_x = out_MPC.Values.Data(:,1);
out_MPC_y = out_MPC.Values.Data(:,2);
out_MPC_z = out_MPC.Values.Data(:,3);
plot3(out_MPC_x,out_MPC_y,out_MPC_z);
hold on
plot3(x_trajectory, y_trajectory, z_trajectory, '-o', 'Color', 'r');
hold off;
xlabel('x')
ylabel('y')
zlabel('z')
legend("Real Trajectory","Reference Trajectory",'Location','northeast')
grid;

%% H-inf plot
figure(10);
out_Hinf = ans.yout.getElement('out_Hinf');

t = out_Hinf.Values.Time;

out_Hinf_x = out_Hinf.Values.Data(:,1);
out_Hinf_y = out_Hinf.Values.Data(:,2);
out_Hinf_z = out_Hinf.Values.Data(:,3);
plot3(out_Hinf_x,out_Hinf_y,out_Hinf_z);
hold on
plot3(x_trajectory, y_trajectory, z_trajectory, '-o', 'Color', 'r');
hold off;
xlabel('x')
ylabel('y')
zlabel('z')
legend("Real Trajectory","Reference Trajectory",'Location','northeast')
grid;

%% PID plot
figure(11);
out_PID = ans.yout.getElement('out_PID');

t = out_PID.Values.Time;

out_PID_x = out_PID.Values.Data(:,1);
out_PID_y = out_PID.Values.Data(:,2);
out_PID_z = out_PID.Values.Data(:,3);
plot3(out_PID_x,out_PID_y,out_PID_z);
hold on
plot3(x_trajectory, y_trajectory, z_trajectory, '-o', 'Color', 'r');
hold off;
xlabel('x')
ylabel('y')
zlabel('z')
legend("Real Trajectory","Reference Trajectory",'Location','northeast')
grid;

%% Simulink

sim('LQRFINALEMODEL.slx')

%% LQR plot
figure(12);
out_LQR = ans.yout.getElement('out_LQR');

t = out_LQR.Values.Time;

out_LQR_x = out_LQR.Values.Data(:,1);
out_LQR_y = out_LQR.Values.Data(:,2);
out_LQR_z = out_LQR.Values.Data(:,3);
plot3(out_LQR_x,out_LQR_y,out_LQR_z);
hold on
plot3(x_trajectory, y_trajectory, z_trajectory, '-o', 'Color', 'r');
hold off;
xlabel('x')
ylabel('y')
zlabel('z')
legend("Real Trajectory","Reference Trajectory",'Location','northeast')
grid;

