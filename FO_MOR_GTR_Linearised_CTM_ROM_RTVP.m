%%% School of Computing and Mathematical Sciences
%%% University of Greenwich, London, United Kingdom
%%% Sheikh Hassan
%%% Email: s.hassan@gre.ac.uk; rhr171@hotmail.com
%%% Thermal-Mechanical - First Order - GTR - Linearised

% [1] Hassan, S., Rajaguru, P., Stoyanov, S., Bailey, C. and Tilford, T., 2024.
%     Coupled thermal-mechanical analysis of power electronic modules with finite
%     element method and parametric model order reduction. Power Electronic Devices
%     and Components, 8, p.100063.
% [2] Hassan, S., Stoyanov, S., Rajaguru, P. and Bailey, C., 2026. Reduced-order
%     modelling for thermal–mechanical analysis of power electronic modules. Finite
%     Elements in Analysis and Design, 253, p.104488.

%% Required Matrices
% Dr, Kr, L1r, L1mLr, L0r, V

%% Reading System Matrices
% M_Pre = mmread('System_Matrices/MassCoupled.mtx');         %Mass Matrix
% D_Pre = mmread('System_Matrices/DampCoupled.mtx');         %Damping Matrix
% K_Pre = mmread('System_Matrices/StiffCoupled.mtx');        %Stiffness Matrix
% L1_Pre = sparse(mmread('System_Matrices/LoadCoupled_L1.mtx')); %Load Vector L1
% L0_Pre = sparse(mmread('System_Matrices/LoadCoupled_L0.mtx')); %Load Vector L0
% L1mL0_Pre = L1_Pre - L0_Pre;

% Method for GTR - Transformation - Second Order to First Order
% [K] = [sparse(np,np) -speye(np);K_Pre sparse(np,np)];
% [L1] = [sparse(np,1);L1_Pre];
% [L0] = [sparse(np,1);L0_Pre];
% % [L] = [sparse(np,1) sparse(np,1);L1_Pre L0_Pre];
% [L1mL0] = [sparse(np,1); L1mL0_Pre];
% [D] = [speye(np) sparse(np,np);D_Pre M_Pre];
%

% Using PRIMA
% [V] = PRIMA(D,K,L,order);

% Using Arnoldi
% [V] = Arnoldi(D,K,L,order);

% Using SOAR
% [V] = SOAR(M,D,K,N,order,L);

%% Obtaining Reduced Order Matrices using SOAR
% Mr_Pre = V'*M_Pre*V; %Reduced Order Mass Matrix
% Dr_Pre = V'*D_Pre*V; %Reduced Order Damping Matrix
% Kr_Pre = V'*K_Pre*V; %Reduced Order Stiffness Matrix
% L1r_Pre = V'*L1_Pre; %Reduced Order Load Vector
% L0r_Pre = V'*L0_Pre;
% L1mL0r_Pre = L1r_Pre - L0r_Pre;

%% Method for GTR - Transformation - Second Order to First Order
% rp = size(Kr_Pre,1);
% Kr = [sparse(rp,rp) -speye(rp); Kr sparse(rp,rp)];
% Dr = [speye(rp,rp) sparse(rp,rp);Dr Mr];
% L1r = [sparse(rp,1); L1r];
% L0r = [sparse(rp,1); L0r];

%% Obtaining Reduced Order Matrices using PRIMA/Arnoldi
% Dr = V'*D*V; %Reduced Order Damping Matrix
% Kr = V'*K*V; %Reduced Order Stiffness Matrix
% L1r = V'*L1; %Reduced Order Load Vector
% L0r = V'*L0;
% L1mL0r = L1r - L0r;

%%
V_ODPRIMA = V; % clear V;

np=size(V_ODPRIMA,1);

%% Defining GTR Modelling Parameters
theta = 1.0; % 

n = np;

%Load and Time Steps
dt = 0.01;            % Time Step Size, e.g., 0.05s
NumbofLoads = 1/dt;     % Required Load Steps to Get 1s, e.g., 1s/0.05s=20 Steps
k_NumbLoadStep = 12*5.5*NumbofLoads;% Total Load Steps

% Convergence Check (Incomplete)
% Based on load/temperature increment
tol_convergance_ROM = 2; % Temperature Change per Load Steps

%Create 'Output_Files_ROM' Folder
if exist("Output_Files_ROM","dir")
    fprintf("'Output_Files_ROM' Folder Exists\n");
else
    fprintf("Creating 'Output_Files_ROM' Folder ... \n");
    mkdir Output_Files_ROM
    fprintf("'Output_Files_ROM' Folder Has Been Created.\n");
end

%%
order = size(V_ODPRIMA,2);

%%
t2 = tic;
%%%Start_ROM_Solution%%%

%% Initialisation
m = size(Kr,2);

k = 1;

% Initialisation
BRed(1:m,1:k_NumbLoadStep+1) = sparse(m,k_NumbLoadStep+1);
% if exist("Output_Files/BRed.txt","file")
%     delete Output_Files/BRed.txt
% end
% temp_BRed = full(BRed(1:m,k));
% writematrix(temp_BRed,'Output_Files/BRed.txt','Delimiter','tab')
URed_initial_condition=(Kr\L0r)';        % Initial Condition - State Variable U (Displacement)
URed(1:k_NumbLoadStep+1,1:m)   = sparse(k_NumbLoadStep+1,m);    % State Variable URed (Displacement)
URed(k,1:m)   = URed_initial_condition;
if exist("Output_Files_ROM/URed.txt","file")
    delete Output_Files_ROM/URed.txt
end
temp_URed = full(URed(k,1:m));
writematrix(temp_URed,'Output_Files_ROM/URed.txt','Delimiter','tab')

UReduced_pre(1:n,1:k_NumbLoadStep+1) = sparse(n,k_NumbLoadStep+1);
UReduced_pre(1:n,k) = V_ODPRIMA*URed(k,1:order/2)';
UReduced(1:k_NumbLoadStep+1,1:n) = sparse(k_NumbLoadStep+1,n);
UReduced(k,1:n) = UReduced_pre(1:n,k)';
if exist("Output_Files_ROM/UReduced.txt","file")
    delete Output_Files_ROM/UReduced.txt
end
temp_UReduced = full(UReduced(k,1:n));
writematrix(temp_UReduced,'Output_Files_ROM/UReduced.txt','Delimiter','tab')

dURed(1:k_NumbLoadStep+1,1:m)  = sparse(k_NumbLoadStep+1,m);    % State Variable dURed (Velocity)
% if exist("Output_Files/dURed.txt","file")
%     delete Output_Files/dURed.txt
% end
% temp_dURed = full(dURed(k,1:m));
% writematrix(temp_dURed,'Output_Files/dURed.txt','Delimiter','tab')

% % Initialisation
temp3(1:m) = sparse(m,1);
temp3 = temp3';

% % Initialisation
k = 1;
Time_ROM(1:k_NumbLoadStep+1,1) = sparse(k_NumbLoadStep+1,1);

% % Initialisation
i = 1;
diff_abs_UReduced(1:k_NumbLoadStep+1,1:n) = sparse(k_NumbLoadStep+1,n);
% if exist("Output_Files/diff_abs_UReduced.txt","file")
%     delete Output_Files/diff_abs_UReduced.txt
% end
% temp_diff_abs_UReduced = full(diff_abs_UReduced(k,1:m));
% writematrix(temp_diff_abs_UReduced,'Output_Files/diff_abs_UReduced.txt','Delimiter','tab')

% Pre-computation of GTR Integration Parameter Matrices
ARed = (((1/(theta*dt))*Dr)+Kr);
A1Red = (1/(theta*dt))*Dr;
% %%

%% ROM Solution Starts
% % Visualising State/Load Increments
x_live_ROM(1:k_NumbLoadStep+1,1) = sparse(k_NumbLoadStep+1,1);
x_live_ROM(k,:)=k;
y_live_ROM_Max(1:k_NumbLoadStep+1,1) = sparse(k_NumbLoadStep+1,1);
y_live_ROM_Max(k,1) = sparse(max(abs(UReduced(k,181856))));
y_live_ROM_Min(1:k_NumbLoadStep+1,1) = sparse(k_NumbLoadStep+1,1);
y_live_ROM_Min(k,1) = sparse(min(abs(UReduced(k,:))));

close all
set(0,'DefaultFigureWindowStyle','docked')
f3=figure;
plot(x_live_ROM(k,1),y_live_ROM_Max(k,1),'o',LineWidth=2);
hold on
title('Maximum Value (Absolute)')
xlabel('Load Steps')
ylabel('ROM States - Max')
f4=figure;
plot(x_live_ROM(k,1),y_live_ROM_Min(k,1),'o',LineWidth=2);
hold on
title('Minimum Value (Absolute)')
xlabel('Load Steps')
ylabel('ROM States - Min')
% linkdata(f3,'on')
% linkdata(f4,'on')
drawnow

for  k = 1:k_NumbLoadStep

    fprintf('ROM Load Step Number: %d\n',k+1)

    Time_ROM(k+1) =Time_ROM(k)+dt;

    temp3(1:m)= ((1/(theta*dt))*URed(k,1:m))+(((1-theta)/theta)*dURed(k,1:m));
    
    switch true
        % Cycle 1
        case k<=0.01*NumbofLoads
            L55r = L0r + 0.0*L1mL0r + (k-0*NumbofLoads)*((1*L1mL0r - 0.0*L1mL0r)/((0.01*NumbofLoads)-0*NumbofLoads));
        case (k>0.01*NumbofLoads) && (k<=2.5*NumbofLoads)
            L55r = L0r + L1mL0r;
        case (k>2.5*NumbofLoads) && (k<=2.51*NumbofLoads)
            L55r = L0r + L1mL0r + (k-2.5*NumbofLoads)*((0.0*L1mL0r - 1.0*L1mL0r)/(2.51*NumbofLoads - 2.5*NumbofLoads));
        case (k>2.51*NumbofLoads) && (k<=5.5*NumbofLoads)
            L55r = L0r + 0.0*L1mL0r;
        % Cycle 2
        case (k>5.5*NumbofLoads) && (k<=(1*5.5+0.01)*NumbofLoads)
            L55r = L0r + 0.0*L1mL0r + (k-(1*5.5+0)*NumbofLoads)*((1*L1mL0r - 0.0*L1mL0r)/(((1*5.5+0.01)*NumbofLoads)-(1*5.5+0)*NumbofLoads));
        case (k>(1*5.5+0.01)*NumbofLoads) && (k<=(1*5.5+2.5)*NumbofLoads)
            L55r = L0r + L1mL0r;
        case (k>(1*5.5+2.5)*NumbofLoads) && (k<=(1*5.5+2.51)*NumbofLoads)
            L55r = L0r + L1mL0r + (k-(1*5.5+2.5)*NumbofLoads)*((0.0*L1mL0r - 1.0*L1mL0r)/((1*5.5+2.51)*NumbofLoads - (1*5.5+2.5)*NumbofLoads));
        case (k>(1*5.5+2.51)*NumbofLoads) && (k<=(1*5.5+5.5)*NumbofLoads)
            L55r = L0r + 0.0*L1mL0r;
        %Cycle 3
        case (k>(1*5.5+5.5)*NumbofLoads) && k<=(2*5.5+0.01)*NumbofLoads
            L55r = L0r + 0.0*L1mL0r + (k-(2*5.5+0)*NumbofLoads)*((1*L1mL0r - 0.0*L1mL0r)/(((2*5.5+0.01)*NumbofLoads)-(2*5.5+0)*NumbofLoads));
        case (k>(2*5.5+0.01)*NumbofLoads) && (k<=(2*5.5+2.5)*NumbofLoads)
            L55r = L0r + L1mL0r;
        case (k>(2*5.5+2.5)*NumbofLoads) && (k<=(2*5.5+2.51)*NumbofLoads)
            L55r = L0r + L1mL0r + (k-(2*5.5+2.5)*NumbofLoads)*((0.0*L1mL0r - 1.0*L1mL0r)/((2*5.5+2.51)*NumbofLoads - (2*5.5+2.5)*NumbofLoads));
        case (k>(2*5.5+2.51)*NumbofLoads) && (k<=(2*5.5+5.5)*NumbofLoads)
            L55r = L0r + 0.0*L1mL0r;
        % Cycle 4
        case (k>(2*5.5+5.5)*NumbofLoads) && (k<=(3*5.5+0.01)*NumbofLoads)
            L55r = L0r + 0.0*L1mL0r + (k-(3*5.5+0)*NumbofLoads)*((1*L1mL0r - 0.0*L1mL0r)/(((3*5.5+0.01)*NumbofLoads)-(3*5.5+0)*NumbofLoads));
        case (k>(3*5.5+0.01)*NumbofLoads) && (k<=(3*5.5+2.5)*NumbofLoads)
            L55r = L0r + L1mL0r;
        case (k>(3*5.5+2.5)*NumbofLoads) && (k<=(3*5.5+2.51)*NumbofLoads)
            L55r = L0r + L1mL0r + (k-(3*5.5+2.5)*NumbofLoads)*((0.0*L1mL0r - 1.0*L1mL0r)/((3*5.5+2.51)*NumbofLoads - (3*5.5+2.5)*NumbofLoads));
        case (k>(3*5.5+2.51)*NumbofLoads) && (k<=(3*5.5+5.5)*NumbofLoads)
            L55r = L0r + 0.0*L1mL0r;
        % Cycle 5
        case (k>(3*5.5+5.5)*NumbofLoads) && (k<=(4*5.5+0.01)*NumbofLoads)
            L55r = L0r + 0.0*L1mL0r + (k-(4*5.5+0)*NumbofLoads)*((1*L1mL0r - 0.0*L1mL0r)/(((4*5.5+0.01)*NumbofLoads)-(4*5.5+0)*NumbofLoads));
        case (k>(4*5.5+0.01)*NumbofLoads) && (k<=(4*5.5+2.5)*NumbofLoads)
            L55r = L0r + L1mL0r;
        case (k>(4*5.5+2.5)*NumbofLoads) && (k<=(4*5.5+2.51)*NumbofLoads)
            L55r = L0r + L1mL0r + (k-(4*5.5+2.5)*NumbofLoads)*((0.0*L1mL0r - 1.0*L1mL0r)/((4*5.5+2.51)*NumbofLoads - (4*5.5+2.5)*NumbofLoads));
        case (k>(4*5.5+2.51)*NumbofLoads) && (k<=(4*5.5+5.5)*NumbofLoads)
            L55r = L0r + 0.0*L1mL0r;
        % Cycle 6
        case (k>(4*5.5+5.5)*NumbofLoads) && (k<=(5*5.5+0.01)*NumbofLoads)
            L55r = L0r + 0.0*L1mL0r + (k-(5*5.5+0)*NumbofLoads)*((1*L1mL0r - 0.0*L1mL0r)/(((5*5.5+0.01)*NumbofLoads)-(5*5.5+0)*NumbofLoads));
        case (k>(5*5.5+0.01)*NumbofLoads) && (k<=(5*5.5+2.5)*NumbofLoads)
            L55r = L0r + L1mL0r;
        case (k>(5*5.5+2.5)*NumbofLoads) && (k<=(5*5.5+2.51)*NumbofLoads)
            L55r = L0r + L1mL0r + (k-(5*5.5+2.5)*NumbofLoads)*((0.0*L1mL0r - 1.0*L1mL0r)/((5*5.5+2.51)*NumbofLoads - (5*5.5+2.5)*NumbofLoads));
        case (k>(5*5.5+2.51)*NumbofLoads) && (k<=(5*5.5+5.5)*NumbofLoads)
            L55r = L0r + 0.0*L1mL0r;
        % Cycle 7
        case (k>(5*5.5+5.5)*NumbofLoads) && (k<=(6*5.5+0.01)*NumbofLoads)
            L55r = L0r + 0.0*L1mL0r + (k-(6*5.5+0)*NumbofLoads)*((1*L1mL0r - 0.0*L1mL0r)/(((6*5.5+0.01)*NumbofLoads)-(6*5.5+0)*NumbofLoads));
        case (k>(6*5.5+0.01)*NumbofLoads) && (k<=(6*5.5+2.5)*NumbofLoads)
            L55r = L0r + L1mL0r;
        case (k>(6*5.5+2.5)*NumbofLoads) && (k<=(6*5.5+2.51)*NumbofLoads)
            L55r = L0r + L1mL0r + (k-(6*5.5+2.5)*NumbofLoads)*((0.0*L1mL0r - 1.0*L1mL0r)/((6*5.5+2.51)*NumbofLoads - (6*5.5+2.5)*NumbofLoads));
        case (k>(6*5.5+2.51)*NumbofLoads) && (k<=(6*5.5+5.5)*NumbofLoads)
            L55r = L0r + 0.0*L1mL0r;
        % Cycle 8
        case (k>(6*5.5+5.5)*NumbofLoads) && (k<=(7*5.5+0.01)*NumbofLoads)
            L55r = L0r + 0.0*L1mL0r + (k-(7*5.5+0)*NumbofLoads)*((1*L1mL0r - 0.0*L1mL0r)/(((7*5.5+0.01)*NumbofLoads)-(7*5.5+0)*NumbofLoads));
        case (k>(7*5.5+0.01)*NumbofLoads) && (k<=(7*5.5+2.5)*NumbofLoads)
            L55r = L0r + L1mL0r;
        case (k>(7*5.5+2.5)*NumbofLoads) && (k<=(7*5.5+2.51)*NumbofLoads)
            L55r = L0r + L1mL0r + (k-(7*5.5+2.5)*NumbofLoads)*((0.0*L1mL0r - 1.0*L1mL0r)/((7*5.5+2.51)*NumbofLoads - (7*5.5+2.5)*NumbofLoads));
        case (k>(7*5.5+2.51)*NumbofLoads) && (k<=(7*5.5+5.5)*NumbofLoads)
            L55r = L0r + 0.0*L1mL0r;
        % Cycle 9
        case (k>(7*5.5+5.5)*NumbofLoads) && (k<=(8*5.5+0.01)*NumbofLoads)
            L55r = L0r + 0.0*L1mL0r + (k-(8*5.5+0)*NumbofLoads)*((1*L1mL0r - 0.0*L1mL0r)/(((8*5.5+0.01)*NumbofLoads)-(8*5.5+0)*NumbofLoads));
        case (k>(8*5.5+0.01)*NumbofLoads) && (k<=(8*5.5+2.5)*NumbofLoads)
            L55r = L0r + L1mL0r;
        case (k>(8*5.5+2.5)*NumbofLoads) && (k<=(8*5.5+2.51)*NumbofLoads)
            L55r = L0r + L1mL0r + (k-(8*5.5+2.5)*NumbofLoads)*((0.0*L1mL0r - 1.0*L1mL0r)/((8*5.5+2.51)*NumbofLoads - (8*5.5+2.5)*NumbofLoads));
        case (k>(8*5.5+2.51)*NumbofLoads) && (k<=(8*5.5+5.5)*NumbofLoads)
            L55r = L0r + 0.0*L1mL0r;
        % Cycle 10
        case (k>(8*5.5+5.5)*NumbofLoads) && (k<=(9*5.5+0.01)*NumbofLoads)
            L55r = L0r + 0.0*L1mL0r + (k-(9*5.5+0)*NumbofLoads)*((1*L1mL0r - 0.0*L1mL0r)/(((9*5.5+0.01)*NumbofLoads)-(9*5.5+0)*NumbofLoads));
        case (k>(9*5.5+0.01)*NumbofLoads) && (k<=(9*5.5+2.5)*NumbofLoads)
            L55r = L0r + L1mL0r;
        case (k>(9*5.5+2.5)*NumbofLoads) && (k<=(9*5.5+2.51)*NumbofLoads)
            L55r = L0r + L1mL0r + (k-(9*5.5+2.5)*NumbofLoads)*((0.0*L1mL0r - 1.0*L1mL0r)/((9*5.5+2.51)*NumbofLoads - (9*5.5+2.5)*NumbofLoads));
        case (k>(9*5.5+2.51)*NumbofLoads) && (k<=(9*5.5+5.5)*NumbofLoads)
            L55r = L0r + 0.0*L1mL0r;
        % Cycle 11
        case (k>(9*5.5+5.5)*NumbofLoads) && (k<=(10*5.5+0.01)*NumbofLoads)
            L55r = L0r + 0.0*L1mL0r + (k-(10*5.5+0)*NumbofLoads)*((1*L1mL0r - 0.0*L1mL0r)/(((10*5.5+0.01)*NumbofLoads)-(10*5.5+0)*NumbofLoads));
        case (k>(10*5.5+0.01)*NumbofLoads) && (k<=(10*5.5+2.5)*NumbofLoads)
            L55r = L0r + L1mL0r;
        case (k>(10*5.5+2.5)*NumbofLoads) && (k<=(10*5.5+2.51)*NumbofLoads)
            L55r = L0r + L1mL0r + (k-(10*5.5+2.5)*NumbofLoads)*((0.0*L1mL0r - 1.0*L1mL0r)/((10*5.5+2.51)*NumbofLoads - (10*5.5+2.5)*NumbofLoads));
        case (k>(10*5.5+2.51)*NumbofLoads) && (k<=(10*5.5+5.5)*NumbofLoads)
            L55r = L0r + 0.0*L1mL0r;
        % Cycle 12
        case (k>(10*5.5+5.5)*NumbofLoads) && (k<=(11*5.5+0.01)*NumbofLoads)
            L55r = L0r + 0.0*L1mL0r + (k-(11*5.5+0)*NumbofLoads)*((1*L1mL0r - 0.0*L1mL0r)/(((11*5.5+0.01)*NumbofLoads)-(11*5.5+0)*NumbofLoads));
        case (k>(11*5.5+0.01)*NumbofLoads) && (k<=(11*5.5+2.5)*NumbofLoads)
            L55r = L0r + L1mL0r;
        case (k>(11*5.5+2.5)*NumbofLoads) && (k<=(11*5.5+2.51)*NumbofLoads)
            L55r = L0r + L1mL0r + (k-(11*5.5+2.5)*NumbofLoads)*((0.0*L1mL0r - 1.0*L1mL0r)/((11*5.5+2.51)*NumbofLoads - (11*5.5+2.5)*NumbofLoads));
        case (k>(11*5.5+2.51)*NumbofLoads) && (k<=(11*5.5+5.5)*NumbofLoads)
            L55r = L0r + 0.0*L1mL0r;
        % Cycle 13
        case (k>(11*5.5+5.5)*NumbofLoads) && (k<=(12*5.5+0.01)*NumbofLoads)
            L55r = L0r + 0.0*L1mL0r + (k-(12*5.5+0)*NumbofLoads)*((1*L1mL0r - 0.0*L1mL0r)/(((12*5.5+0.01)*NumbofLoads)-(12*5.5+0)*NumbofLoads));
        case (k>(12*5.5+0.01)*NumbofLoads) && (k<=(12*5.5+2.5)*NumbofLoads)
            L55r = L0r + L1mL0r;
        case (k>(12*5.5+2.5)*NumbofLoads) && (k<=(12*5.5+2.51)*NumbofLoads)
            L55r = L0r + L1mL0r + (k-(12*5.5+2.5)*NumbofLoads)*((0.0*L1mL0r - 1.0*L1mL0r)/((12*5.5+2.51)*NumbofLoads - (12*5.5+2.5)*NumbofLoads));
        case (k>(12*5.5+2.51)*NumbofLoads) && (k<=(12*5.5+5.5)*NumbofLoads)
            L55r = L0r + 0.0*L1mL0r;
    end

    Loadr = L55r;

    BRed(:,k+1) = Loadr + (Dr*temp3);%     Bred=Load+Dr*temp3';

    URed(k+1,1:m)  = ARed\BRed(:,k+1);
    dURed(k+1,1:m)=(URed(k+1,1:m)-URed(k,1:m)-(1-theta)*dt*dURed(k,1:m))/(theta*dt);

    % Projection Back to Full Dimension
    UReduced_pre(1:n,k+1) = V_ODPRIMA*URed(k+1,1:order/2)';
    UReduced(k+1,1:n) = UReduced_pre(1:n,k+1)';
    
    x_live_ROM(k+1,:)=k;
    y_live_ROM_Max(k+1,:)=max(abs(UReduced(k+1,181856)));
    figure(f3)
    plot(x_live_ROM(k+1,:),y_live_ROM_Max(k+1,:),'o',LineWidth=2);
    y_live_ROM_Min(k+1,:)=min(abs(UReduced(k+1,:)));
    figure(f4)
    plot(x_live_ROM(k+1,:),y_live_ROM_Min(k+1,:),'o',LineWidth=2);
    drawnow limitrate

    fprintf('ROM Load Step Number: %d\n',k)
    diff_abs_UReduced(k+1,:)=abs(UReduced(k+1,:))-abs(UReduced(k,:));
    switch true
        case any(diff_abs_UReduced(k+1,:) > tol_convergance_ROM)
            fprintf('ROM Solution Have Not Converged For Load Step %d\n',k)
        %break;
            return;
        otherwise
            fprintf('ROM Solution Converged For Load Step %d\n',k)
    end

    % % Reduce Computing Memory
    % temp_BRed = full(BRed(:,k+1));
    % writematrix(temp_BRed,'Output_Files/BRed.txt','Delimiter','tab','WriteMode','append')
    temp_URed = full(URed(k+1,:));
    writematrix(temp_URed,'Output_Files_ROM/URed.txt','Delimiter','tab','WriteMode','append')
    temp_UReduced = full(UReduced(k+1,:));
    writematrix(temp_UReduced,'Output_Files_ROM/UReduced.txt','Delimiter','tab','WriteMode','append')
    % temp_dURed = full(dURed(k+1,:));
    % writematrix(temp_dURed,'Output_Files/dURed.txt','Delimiter','tab','WriteMode','append')
    % temp_diff_abs_UReduced = full(diff_abs_UReduced(k+1,:));
    % writematrix(temp_diff_abs_UReduced,'Output_Files/diff_abs_UReduced.txt','Delimiter','tab','WriteMode','append')
    if k>= 2
        BRed(:,k-1) = sparse(m,1);
        URed(k-1,:) = sparse(1,m);
        UReduced_pre(:,k-1) = sparse(n,1);
        UReduced(k-1,:) = sparse(1,n);
        dURed(k-1,:) = sparse(1,m);
        diff_abs_UReduced(k-1,:) = sparse(1,n);
    end
end

drawnow
%%

%%
clear ARed;
clear BRed;
clear URed;
clear UReduced;
clear diff_abs_UReduced;
clear URed_initial_condition;
clear dURed;
clear L55r;
clear Loadr;
clear temp3;
clear tol_convergance_ROM;

clear A1Red;
clear deltaURed;
clear UReduced_pre;
clear temp_BRed temp_diff_abs_UReduced temp_dURed temp_URed temp_UReduced;

%%%%End_ROM_Solution%%%
CompTimeROM = toc(t2);

%% Read Results from Storage
URed = sparse(readmatrix('Output_Files_ROM/URed.txt','Delimiter','\t'));
UReduced = sparse(readmatrix('Output_Files_ROM/UReduced.txt','Delimiter','\t'));