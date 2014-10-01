%%%% Hodgkin Huxley Model %%%%
% Author: Benjamin Ceron

%%% Initialize Time and Constants %%%
dt = .01;
t = 0:dt:100;        %% time in ms
g_Kbar = 36;         %% maximum conductance of potassium mS/cm^2
g_Nabar = 120;       %% maximum conductance of sodium mS/cm^2
g_Lbar = 0.3;        %% maximum conductance of leakage ions mS/cm^2
Ek = -12;            %% Nernst potential of potassium mV
Ena = 115;           %% Nernst potential of sodium mV
El = 10.6;           %% Nernst potential of leakage mV
Vrest = -70;         %% Resting potential mV
Cm = 1.0;            %% Membrane capacitance uF/cm^2

%%% Initialize Initial Conditions %%%
Vm =0;               %% Membrane voltage mV
VmVector = zeros(1,(100/dt)+1);
VmVector(1) = Vm;

%%% Initial Activate Probabilities %%%
AlphaM = 0.1 * ((25 - Vm)/(exp((25 - Vm)/10)-1));
BetaM = 4*exp(-Vm/18);
%%% Initial Inactivate Probabilities %%%
AlphaH = .07*exp(-Vm/20);
BetaH = 1/(exp((30-Vm)/10) + 1);
%%% Initial Deactivate Probabilities %%%
AlphaN = .01 * ((10 - Vm)/(exp((10-Vm)/10) - 1));
BetaN = .125*exp(-Vm/80);   

%% Initial Rates of Changes %%
m = AlphaM/(AlphaM + BetaM);
n = AlphaN/(AlphaN + BetaN);
h = AlphaH/(AlphaH + BetaH);


%% Initialize Injected Current %%
%I_inj_init = 0;                                                  %% Steady-state neuron. USE FOR PART 1 - comment out every other I_inj line except this one if using for part 1
                                                                  %% adjusted to 50 uA/cm^2 to get significant response
%I_inj = [50*ones(1,(.5/dt)+1), zeros(1, 100/dt - .5/dt)];        %%Injected pulse of 5 uA/cm^2 - USE FOR PART 2 - adjusted to 50 uA/cm^2 
I_inj = [50*ones(1,100/dt + 1)];                                  %% Constant current of 5 uA/cm^2 - USE FOR PART 3 - adjusted to 50 uA/cm^2 to get significant response


%% Initialize Rates of Change Vectors %%%
mVector = [m];
nVector = [n];
hVector = [h];


%% Loop through and calculate Vm and conductances for 100 ms %%
for i = 1:1:(100/dt)
    %%% Calculate/Update Currents of Ions %%%
    I_Na = (m^3)*g_Nabar*h*(Vm - Ena);
    I_K = (n^4)*g_Kbar*(Vm - Ek);
    I_L = g_Lbar*(Vm - El);
    I_ion = (I_inj(i+1))- I_K - I_Na - I_L;          %%% in uA/cm^2
    
    %%% Calculate membrane voltage rate of change and update new Vm %%%
    dVdt = I_ion/Cm;        %% in mV/ms
    Vm = Vm + dt*(dVdt);
     
    %%% Store Vm in vector %%%
    VmVector(i+1) = Vm;

    %%% Update variables by differential eqns %%%
      
    %%% Activate Probabilities %%%
    AlphaM = 0.1 * ((25- Vm)/(exp((25 - Vm)/10)-1));
    BetaM = 4*exp(-Vm/18);
    %%% Inactivate Probabilities %%%
    AlphaH = .07*exp(-Vm/20);
    BetaH = 1/(exp((30-Vm)/10) + 1);
    %%% Deactivate Probabilities %%%
    AlphaN = .01 * ((10 - Vm)/(exp((10-Vm)/10) - 1));
    BetaN = .125*exp(-Vm/80);   
    
    %%% Rates of Change Differentials %%%
    dmdt = AlphaM*(1-m)-BetaM*m;
    dndt = AlphaN*(1-n)-BetaN*n;
    dhdt = AlphaH*(1-h)-BetaH*h;
    
    %%% Update Rates of Change %%%
    m = m + dt*(dmdt);
    n = n + dt*(dndt);
    h = h + dt*(dhdt);
    
    %%% Store Updated Rates in Vectors %%%
    mVector = [mVector m];
    nVector = [nVector n];
    hVector = [hVector h];
end
 
%%% Add resting potential to vector of Vm's %%%
VmVector = VmVector + Vrest;

%%% Calculate conductances %%%
mVector = mVector .^ 3;
nVector = nVector .^ 4;

gNa_Vector = g_Nabar .* mVector .* hVector;
gK_Vector = nVector .* g_Kbar;

%%% Plot membrane voltages and conductances over time %%%
figure(1)
plot(t,VmVector)
axis([0 100 -100 40])
xlabel('Time (ms)')
ylabel('Voltage (mV)')
legend('Voltage')
title('Membrane Potential')

figure(2)
hold on
plot(t,gK_Vector, 'g')          %% g_K is green 
plot(t,gNa_Vector, 'r')         %% g_Na is red
axis([0 100 0 40])
xlabel('Time (ms)')
ylabel('Conductance (mS/cm2)')
legend('gK','gNa')
title('gK and gNa')
    
    

