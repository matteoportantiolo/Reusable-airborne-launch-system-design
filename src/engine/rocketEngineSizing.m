clear;
close all;
clc;

%% Set problem parameter
data.pChamber = 150;            % Chamber pressure in bar
data.pLaunch = 0.26;            % Ambient pressure at launch altitude in bar
data.engineNumber = 5;          % Number of engine desired
data.stackMass = 30000;         % Total mass of the stack in kg
data.thrustWeight = 1;          % Set Thrust to Weight ratio

% To set different propellant couple è un parto dp, change all the cea
% callings

%% Plot and export
settings.plotFlag = true;      % True for plotting
settings.exportFlag = true;    % True for exporting file .mat

%% Oxidizer-to-Fuel Ratio (OF) Evaluation
ofSpan = 1:0.01:5;  % Define range of OF ratios
pChamber = data.pChamber;

% Preallocate arrays for specific impulse and temperature
isp = zeros(length(ofSpan), 1);
temp = zeros(length(ofSpan), 1);

% Evaluate Isp and temperature for each OF ratio
for i = 1:length(ofSpan)
    out = CEA('problem', 'rkt', 'frozen', 'nfz', 1, 'o/f', ofSpan(i), 'p(bar)', pChamber, ...
              'reactants', 'fuel', 'RP-1(L)', 'C', 1., 'H', 1.9423, 'wt%', 100., ...
              'h,cal/mol', -5430., 't(k)', 300.0, 'oxid', 'O2(L)', 'O', 2, 'wt%', 100., ...
              'h,cal/mol', -3032., 't(k)', 94.44, 'output', 'joule', 'transport', 'end');
    
    % Extract specific impulse and chamber temperature
    isp(i) = out.output.froz.isp(end);
    temp(i) = out.output.froz.temperature(1);
end

if settings.plotFlag
    % Plot specific impulse and temperature vs OF ratio
    figure;
    plot(ofSpan, isp, 'DisplayName', 'Specific Impulse');
    ylabel('Specific Impulse [s]');
    hold on;
    yyaxis right;
    plot(ofSpan, temp, 'DisplayName', 'Chamber Temperature');
    ylabel('Chamber Temperature [K]');
    xlabel('Oxidizer-to-Fuel Ratio [OF]');
    legend;
    title('Specific Impulse and Temperature vs OF Ratio with only Convergent Nozzle');
    grid minor;
end

%% Optimal OF Choice
[~, idx] = max(isp);
of = ofSpan(idx);  % Choose OF that maximizes Isp

% Run CEA for optimal OF
out = CEA('problem', 'rkt', 'frozen', 'nfz', 1, 'o/f', of, 'p(bar)', pChamber, ...
          'reactants', 'fuel', 'RP-1(L)', 'C', 1., 'H', 1.9423, 'wt%', 100., ...
          'h,cal/mol', -5430., 't(k)', 300.0, 'oxid', 'O2(L)', 'O', 2, 'wt%', 100., ...
          'h,cal/mol', -3032., 't(k)', 94.44, 'output', 'joule', 'transport', 'end');

% Extract gamma values for chamber and throat
gamma = out.output.froz.gamma(1);
gammaThroat = out.output.froz.gamma(2);

%% Expansion Ratio Calculation
pExitDesign = data.pLaunch;  % Exit pressure at 10 km altitude
pRatio = pExitDesign / pChamber;

% Define isentropic pressure relation for Mach calculation
isoPressureRelation = @(x) (1 + (gamma - 1) / 2 * x^2)^(-gamma / (gamma - 1)) - pRatio;
machExit = fzero(isoPressureRelation, 5);  % Exit Mach number

% Calculate expansion ratio
exRatio = 1 / machExit * ((2 / (gammaThroat + 1)) * (1 + (gammaThroat - 1) / 2 * machExit^2))^...
          ((gammaThroat + 1) / (2 * (gammaThroat - 1)));

%% Nozzle CEA Run
out = CEA('problem', 'rkt', 'frozen', 'nfz', 1, 'o/f', of, 'sup,ae/at', exRatio, 'p(bar)', pChamber, ...
          'reactants', 'fuel', 'RP-1(L)', 'C', 1., 'H', 1.9423, 'wt%', 100., ...
          'h,cal/mol', -5430., 't(k)', 300.0, 'oxid', 'O2(L)', 'O', 2, 'wt%', 100., ...
          'h,cal/mol', -3032., 't(k)', 94.44, 'output', 'joule', 'transport', 'end');

% Extract Isp and characteristic velocity
isp = out.output.froz.isp(end);
cStar = out.output.froz.cstar(1);

%% Chamber Sizing Calculations
g0 = 9.81;                                                      % Gravity acceleration (m/s^2)
thrust = data.thrustWeight*data.stackMass*g0/data.engineNumber; % Target thrust in Newtons
dm = thrust / (isp * g0);                                       % Mass flow rate

lStar = 0.9;  % Characteristic length in meters
aThroat = cStar * dm / (pChamber * 1e5);
dThroat = 2 * sqrt(aThroat / pi);

% Calculate contraction ratio using Huzel-Huang relation
cRatio = 8 * (dThroat * 100)^(-0.6) + 1.25;

% Determine chamber area and diameter
aChamber = cRatio * aThroat;
dChamber = 2 * sqrt(aChamber / pi);

% Chamber volume and length
vChamber = aThroat * lStar;
lChamber = vChamber / aChamber;

% Mach in chamber to check Huzel-Huang validity
machCCFun = @(x) (pi * dChamber^2 / 4) / (pi * dThroat^2 / 4) - 1 / x * ...
                ((2 / (gammaThroat + 1)) * (1 + (gammaThroat - 1) / 2 * x^2))^...
                ((gammaThroat + 1) / (2 * (gammaThroat - 1)));
machChamber = fzero(machCCFun, 0.2);

%% Thrust vs Ambient Pressure
pAmbRange = pExitDesign+0.1:-0.00001:0;
pExitReal = out.output.froz.pressure(end);

% Exit velocity calculation
vExit = out.output.froz.mach(end) * sqrt(out.output.froz.gamma(end) * 8314.15 / ...
                                         out.output.froz.mw(end) * out.output.froz.temperature(end));

% Calculate thrust over pressure range
thrustOverPressure = dm * vExit + (pExitReal - pAmbRange) * 1e5 * aThroat * exRatio;

if settings.plotFlag
    % Plot thrust vs ambient pressure
    figure;
    plot(pAmbRange, thrustOverPressure / 1e3, 'DisplayName', 'Thrust vs Ambient Pressure');
    xlabel('Ambient Pressure [bar]');
    ylabel('Thrust [kN]');
    title('Thrust Performance at Various Ambient Pressures');
    legend;
    xlim([pAmbRange(end) pAmbRange(1)]);
    grid minor;
end

if settings.exportFlag
    engine = struct(...
        'chamberPressure', pChamber, ...
        'launchPressure', data.pLaunch, ...
        'engineNumber', data.engineNumber, ...
        'of', of, ...
        'gammaChamber', gamma, ...
        'gammaThroat', gammaThroat, ...
        'expansionRatio', exRatio, ...
        'isp', isp, ...
        'cStar', cStar, ...
        'thrust', thrust, ...
        'massFlowRate', dm, ...
        'throatArea', aThroat, ...
        'throatDiameter', dThroat, ...
        'contractionRatio', cRatio, ...
        'chamberArea', aChamber, ...
        'chamberDiameter', dChamber, ...
        'chamberVolume', vChamber, ...
        'chamberLength', lChamber, ...
        'machChamber', machChamber, ...
        'exitMach', machExit, ...
        'exitVelocity', vExit, ...
        'thrustOverPressure', thrustOverPressure, ...
        'chamberTemperature', temp(idx)...
    );

    engine = orderfields(engine);

    save('engine', "engine");
end