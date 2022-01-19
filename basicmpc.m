%BASICMPC Basic Predictive Control without constraints. (Script file)
%
% Implements unconstrained MPC for stable SISO discrete-time system 
% defined as LTI object.
%
% The following are editable parameters (most have defaults):
% Tref: Time constant of exponential reference trajectory
% Ts:   Sampling interval 
% plant: Plant definition  (discrete-time SISO LTI object)
% model: Model definition  (discrete-time SISO LTI object)
% P: Vector of coincidence points
% M: Control horizon
% tend: Duration of simulation
% setpoint: Setpoint trajectory (column vector) - length must exceed
%           no of steps in simulation by at least max(P).
% umpast, uppast, ympast, yppast: Initial conditions of plant & model.
%
% Assumes Matlab 5.3. Uses Control Toolbox LTI object class. 
% No other toolboxes required.

%% J.M.Maciejowski, 8 March 1999. Revised 22.3.99, 15.12.99.
%% Copyright(C) 1999. All Rights Reserved.
%% Cambridge University Engineering Department.

%%%%%%%%%%%%%%%%% PROBLEM DEFINITION:

% Define time-constant of reference trajectory Tref:
Tref = 6;

% Define sampling interval Ts (default Tref/10):
if Tref == 0,
  Ts = 1;
else
  Ts = Tref/10;
end

% Define plant as SISO discrete-time 'lti' object 'plant'
%%%%%  CHANGE FROM  HERE TO DEFINE NEW PLANT %%%%%
nump=1;
denp=[1,-1.4,0.45];
plant = tf(nump,denp,Ts);
%%%%%  CHANGE UP TO HERE TO DEFINE NEW PLANT  %%%%%
plant = tf(plant);  % Coerce to transfer function form
nump = get(plant,'num'); nump = nump{:}; % Get numerator polynomial
denp = get(plant,'den'); denp = denp{:}; % Get denominator polynomial
nnump = length(nump)-1; % Degree of plant numerator
ndenp = length(denp)-1; % Degree of plant denominator
if nump(1)~=0, error('Plant must be strictly proper'), end;
if any(abs(roots(denp))>1), disp('Warning: Unstable plant'), end

% Define model as SISO discrete-time 'lti' object 'model' 
% (default model=plant):
%%%%%  CHANGE FROM  HERE TO DEFINE NEW MODEL %%%%%
model = plant;
%%%%%  CHANGE UP TO HERE TO DEFINE NEW MODEL  %%%%%
model = tf(model);  % Coerce to transfer function form
numm = get(model,'num'); numm = numm{:}; % Get numerator polynomial
denm = get(model,'den'); denm = denm{:}; % Get denominator polynomial
nnumm = length(numm)-1; % Degree of model numerator
ndenm = length(denm)-1; % Degree of model denominator
if numm(1)~=0, error('Model must be strictly proper'), end;
if any(abs(roots(denm))>1), disp('Warning: Unstable model'), end

nump=[zeros(1,ndenp-nnump-1),nump]; % Pad numerator with leading zeros
numm=[zeros(1,ndenm-nnumm-1),numm]; % Pad numerator with leading zeros

% Define prediction horizon P (steps)(default corresponds to 0.8*Tref):
if Tref == 0,
  P = 5;
else
  P = round(0.8*Tref/Ts);
end

% Define control horizon (default 1):
M = 1;

% Compute model step response values over coincidence horizon:
stepresp = step(model,[0:Ts:max(P)*Ts]);
theta = zeros(length(P),M);
for j=1:length(P),
  theta(j,:) = [stepresp(P(j):-1:max(P(j)-M+1,1))',zeros(1,M-P(j))];
end
S = stepresp(P);

% Compute reference error factor at coincidence points:
% (Exponential approach of reference trajectory to set-point)
if Tref == 0,  % Immediate jump back to set-point trajectory
  errfac = zeros(length(P),1);
else
  errfac = exp(-P*Ts/Tref);
end

%%%%%%%%%%%%%%%% SIMULATION PARAMETERS:

if Tref == 0,
  tend = 100*Ts;
else
  tend = 10*Tref; % Duration of simulation (default 10*Tref)
end
nsteps = floor(tend/Ts); % (Number of steps in simulation).
tvec = (0:nsteps-1)'*Ts; % Column vector of time points (first one 0)

% Define set-point (column) vector (default constant 1):
setpoint = ones(nsteps+max(P),1);

% Define vectors to hold input and output signals, initialised to 0:
uu = zeros(nsteps,1); % Input
yp = zeros(nsteps,1); % Plant Output
ym = zeros(nsteps,1); % Model Output

% Initial conditions:
umpast = zeros(ndenm,1);
uppast = zeros(ndenp,1);
ympast = zeros(ndenm,1);  % For model response
yppast = zeros(ndenp,1);  % For plant response

%%%%%%%%%%%%%%%% SIMULATION:

for k=1:nsteps,

  % Define reference trajectory at coincidence points:
  errornow = setpoint(k)-yp(k);
  reftraj = setpoint(k+P) - errornow*errfac;

  % Free response of model over prediction horizon:
  yfpast = ympast;
  ufpast = umpast;
  for kk=1:max(P), % Prediction horizon
    ymfree(kk) = numm(2:nnumm+1)*ufpast-denm(2:ndenm+1)*yfpast;
    yfpast=[ymfree(kk);yfpast(1:length(yfpast)-1)];
	ufpast=[ufpast(1);ufpast(1:length(ufpast)-1)];
  end
  
  % Compute input signal uu(k):
  if k>1,
    dutraj = theta\(reftraj-ymfree(P)');
    uu(k) = dutraj(1) + uu(k-1);  
  else
	dutraj = theta\(reftraj-ymfree(P)');
	uu(k) = dutraj(1) + umpast(1);
  end
  
  % Simulate plant:
  % Update past plant inputs
  uppast = [uu(k);uppast(1:length(uppast)-1)];
  yp(k+1) = -denp(2:ndenp+1)*yppast+nump(2:nnump+1)*uppast; % Simulation
  % Update past plant outputs
  yppast = [yp(k+1);yppast(1:length(yppast)-1)];


  % Simulate model:
  % Update past model inputs
  umpast = [uu(k);umpast(1:length(umpast)-1)];
  ym(k+1) = -denm(2:ndenm+1)*ympast+numm(2:nnumm+1)*umpast;  % Simulation
  % Update past model outputs
  ympast = [ym(k+1);ympast(1:length(ympast)-1)];

end % of simulation

%%%%%%%%%%%%%%%% PRESENTATION OF RESULTS:

disp('***** Results from script file BASICMPC :')
disp(['Tref = ',num2str(Tref),',  Ts = ',num2str(Ts),...
	  '  P = ',int2str(P'),' (steps),  M = ',int2str(M)])
diffpm = get(plant-model,'num');
if diffpm{:}==0,
  disp('Model = Plant')
else
  disp('Plant-Model mismatch')
end

figure
subplot(211)
% Plot output, solid line and set-point, dashed line:
plot(tvec,yp(1:nsteps),'-',tvec,setpoint(1:nsteps),'--');
grid; title('Plant output (solid) and set-point (dashed)')
xlabel('Time')

subplot(212)
% plot input signal as staircase graph:
stairs(tvec,uu,'-');
grid; title('Input')
xlabel('Time')
