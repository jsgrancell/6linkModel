% Function to find the loads through each suspension arm. Takes input of
% car parameters and g-force (NOT ACCELERATION) data entered in a matrix of 
% corresponding longitudinal (x) and lateral (y) values. The g-force data 
% must be organized in this order to produce a valid return. 

% This function can also take special cases of load table data. This data 
% must be entered in the form [WT lateral, WT longitudinal, Fz% 
% (percent of vehicle weight on wheel of interest), Fx_car (total 
% longitidunal force on the car), Fy_car (total lateral force on the car), 
% Fx (long force on wheel of interest), Fy (like Fx but lat), Fz (like Fx 
% but vertical)]. This is used for special cases like max bump or reverse 
% braking, where the load table algorithms don't apply.

% This model assumes that the lateral force of each tire is proportional
% to its share of the whole vehicle's lateral force (per the suspension
% solver sheet)

% To find ONLY the forces for each load case, run this program with the 
% fe12params structure and an array of Gs in the format:
% [forcesF, ForcesR] = suspensionSolverPlot(fe12params, Gs)

% To find the forces for each load case AND the overall max loads, 
% run this program with the fe12params structure and an array of Gs in the format:
% [forcesF, ForcesR, frontMaxes, rearMaxes] = suspensionSolverPlot(fe12params, Gs)


function [forcesF, forcesR, frontMaxes, rearMaxes, loadsF, loadsR] = suspensionSolverCases(carParams, Gs, specialCasesF, specialCasesR)
    inboardF = carParams.inboardF;
    outboardF = carParams.outboardF;
    locationsF = [carParams.tireContactPtF(1,1), carParams.tireContactPtF(1,2), -8*25.4];
    locationsF = [locationsF; locationsF; locationsF; locationsF; locationsF; locationsF];
    inboardR = carParams.inboardR;
    outboardR = carParams.outboardR;
    locationsR = [carParams.tireContactPtR(1,1), carParams.tireContactPtR(1,2), -8*25.4];
    locationsR = [locationsR; locationsR; locationsR; locationsR; locationsR; locationsR];
    RF = locationsF-[carParams.tireContactPtF; carParams.tireContactPtF];
    RR = locationsR-[carParams.tireContactPtR; carParams.tireContactPtR];
    rF = outboardF - locationsF;
    rR = outboardR - locationsR;
    A = zeros(6);                                                          % Initializes a matrix to store unit vectors and moment vectors
    [loadTableF, loadTableR] = loadCases(carParams, Gs(:, 1:2))           % Creates a table of x, y, and z loads using car parameters and acceleration pairs
    if exist('specialCasesF', 'var')
        for i = 1:length(specialCasesF(1))
            loadTableF(size(loadTableF,1)+i, :) = specialCasesF;
        end
    end
    if exist('specialCasesR', 'var')
        for i = 1:length(specialCasesR(1))
            loadTableR(size(loadTableR,1)+i, :) = specialCasesR;
        end
    end
    forcesF = zeros(size(loadTableF,1),6);                                 % Initializes a matrix to store the forces with a row for each recorded acceleration pair
    forcesR = zeros(size(loadTableR,1),6);
    for i = 1:6                                                            % For each arm
        link = (outboardF(i,:)-inboardF(i,:))';                            % Finds the vector components of each link
        norm(link);
        A(1:3,i) = link/norm(link);                                        % Turns each component into a unit vector and places it within A
    end
    for i = 1:6                                                            % For each arm
        A(4,i) = A(3,i)*rF(i,2)-A(2,i)*rF(i,3);                            % Calculates the moment vector and places said vector in A
        A(5,i) = A(3,i)*rF(i,1)-A(1,i)*rF(i,3);
        A(6,i) = A(2,i)*rF(i,1)-A(1,i)*rF(i,2);
    end
    for i = 1:size(loadTableF,1)                                           % For each lat/long acceleration pair
        B(1:3,1) = loadTableF(i, 3:5)';                                    % Creates a diagonal matrix with x, y, and z forces
        B(4,1) = B(3,1)*RF(1,2)-B(2,1)*RF(1,3);
        B(5,1) = B(1,1)*RF(1,3)-B(3,1)*RF(1,1);
        B(6,1) = B(2,1)*RF(1,1)-B(1,1)*RF(1,2);
        forcesF(i,:) = (A\B)';                                             % Fills the respective row of the forces matrix with the forces through each arm
    end
    frontMaxes = findMaxes(forcesF);
    frontMaxes.Properties.VariableNames = {'Up-Fore','Up-Aft','Low-Fore','Low-Aft','Pushrod','Tie Rod'};
    frontMaxes.Properties.RowNames = {'Max Tension', 'Max Compression'};
    forcesF = array2table(forcesF);
    forcesF.Properties.VariableNames = {'Up-Fore','Up-Aft','Low-Fore','Low-Aft','Pushrod','Tie Rod'};
    forcesF.Properties.RowNames = {'Max Accel', 'Max Braking', 'Max Cornering Left', 'Max Cornering Right', 'Combined Accel/Cornering Left', 'Combined Accel/Cornering Right', 'Combined Braking/Cornering Left', 'Combined Braking/Cornering Right', 'Max Bump'};
    loadsF = findMaxes(loadTableF);
    loadsF = loadsF(:,3:5);
    loadsF.Properties.VariableNames = {'Longitudinal (+ forwards, - backwards)', 'Lateral (+ right, - left)', 'Vertical (+ down, - up)'};
    A = zeros(6);
    for i = 1:6                                                            
        link = (outboardR(i,:)-inboardR(i,:))';                            
        A(1:3,i) = link/norm(link);                                        
    end                                  
    for i = 1:6                                                            
        A(4,i) = A(3,i)*rR(i,2)-A(2,i)*rR(i,3);                            % Calculates the moment vector and places said vector in A
        A(5,i) = A(3,i)*rR(i,1)-A(1,i)*rR(i,3);
        A(6,i) = A(2,i)*rR(i,1)-A(1,i)*rR(i,2);                      
    end
    for i = 1:size(loadTableR,1) 
        B(1:3,1) = loadTableR(i, 3:5)';                                    % Creates a diagonal matrix with x, y, and z forces
        B(4,1) = B(3,1)*RR(1,2)-B(2,1)*RR(1,3);
        B(5,1) = B(1,1)*RR(1,3)-B(3,1)*RR(1,1);
        B(6,1) = B(2,1)*RR(1,1)-B(1,1)*RR(1,2);
        forcesR(i,:) = (A\B)';     
    end
    rearMaxes = findMaxes(forcesR);
    rearMaxes.Properties.VariableNames = {'Up-Fore','Up-Aft','Low-Fore','Low-Aft','Pushrod','Toe Rod'};
    rearMaxes.Properties.RowNames = {'Max Tension', 'Max Compression'};
    forcesR = array2table(forcesR);
    forcesR.Properties.VariableNames = {'Up-Fore','Up-Aft','Low-Fore','Low-Aft','Pushrod','Toe Rod'};
    forcesR.Properties.RowNames = {'Max Accel', 'Max Braking', 'Max Cornering Left', 'Max Cornering Right', 'Combined Accel/Cornering Left', 'Combined Accel/Cornering Right', 'Combined Braking/Cornering Left', 'Combined Braking/Cornering Right', 'Max Bump'};
    loadsR = findMaxes(loadTableR);
    loadsR = loadsR(:,3:5);
    loadsR.Properties.VariableNames = {'Longitudinal (+ forwards, - backwards)', 'Lateral (+ right, - left)', 'Vertical (+ down, - up)'};
end

% Function to find forces and weight transfer from a structure of car
% parameters and an array of lateral and longitudinal acceleration

function [loadTableF, loadTableR] = loadCases(carParams, accelData)
    loadTableF = zeros(size(accelData,1), 5);                              % Initializes a load table matrix to store the outputs of this function
    loadTableR = loadTableF;
    for i = 1:size(accelData,1)                                            % For each provided value of lat and long G's
        [loadTableF(i,1), loadTableR(i,1)] = SampoWeightTransfer(carParams, accelData(i,2)*9.81);
        loadTableF(i,2) = (carParams.m*9.81*accelData(i,1)*carParams.hCG)/carParams.WB;   % Calculates longitudinal WT
        loadTableF(i,5) = -(carParams.m*9.81*carParams.PFront/2+loadTableF(i,1)-loadTableF(i,2)/2);    % Calculates Fz
        loadTableR(i,2) = (carParams.m*9.81*accelData(i,1)*carParams.hCG)/carParams.WB;   % Calculates longitudinal WT
        loadTableR(i,5) = -(carParams.m*9.81*(1-carParams.PFront)/2+loadTableR(i,1)+loadTableR(i,2)/2);    % Calculates Fz
        if accelData(i,1) > 0 && accelData(i,2) == 0
            loadTableF(i,3) = 0;
            loadTableR(i,3) = ((1.25*carParams.m*9.81*carParams.a_s)/carParams.WB)/(1-carParams.hCG/carParams.WB*1.25)/2;
        elseif accelData(i,1) < 0
            loadTableF(i,3) = 1.25*loadTableF(i,5);
            loadTableR(i,3) = 1.25*loadTableR(i,5);
        elseif accelData(i,1) == 0 && accelData(i,2) ~= 0
            loadTableF(i,3) = 0;
            loadTableR(i,3) = 0;
        elseif accelData(i,1) > 0 && accelData(i,2) ~= 0
            loadTableF(i,3) = 0;
            loadTableR(i,3) = 1.25*-loadTableR(i,5);
        end
        loadTableF(i,4) = -loadTableF(i,5)*accelData(i,2);                  % Calculates Fy
        loadTableR(i,4) = -loadTableR(i,5)*accelData(i,2);                  % Calculates Fy
    end
end

% Function to calculate from and rear weight transfer using chassis
% flexibility

function [deltaFzFront, deltaFzRear] = SampoWeightTransfer(carParams, Ay)
    kF = carParams.kF;                                                     % Front roll stiffness
    kR = carParams.kR;                                                     % Rear roll stiffness
    kC = carParams.kC;                                                     % Chassis torsional stiffness
    m_uF = carParams.m_uF;                                                 % Front unsprung mass
    m_uR = carParams.m_uR;                                                 % Rear unsprung mass
    TWf = carParams.TWf;                                                   % Front track width
    TWr = carParams.TWr;                                                   % Rear track width
    h_uF = carParams.h_uF;                                                 % Height of front unsprung center of gravity
    h_uR = carParams.h_uR;                                                 % Height of rear unsprung center of gravity
    zF = carParams.zF;                                                     % Height of roll axis in the front
    zR = carParams.zR;                                                     % Height of roll axis in the rear
    d_sF = carParams.h_sF - zF;                                            % Distance between front unsprung center of mass and roll axis
    d_sR = carParams.h_sR - zR;                                            % Distance between rear unsprung center of mass and roll axis
    m_sF = carParams.m_s*carParams.b_s/carParams.WB;                       % Front sprung mass
    m_sR = carParams.m_s*carParams.a_s/carParams.WB;                       % Rear sprung mass
    deltaFzFront = ((kF*d_sF*m_sF)/(kF+(kR*kC/(kR+kC)))+((kF*kC/(kF+kC))*d_sR*m_sR)/(kF*kC/(kF+kC)+kR)+zF*m_sF+h_uF*m_uF)*Ay/TWf;
    deltaFzRear = ((kR*kC/(kR+kC)*d_sF*m_sF)/(kF+kR*kC/(kR+kC))+(kR*d_sR*m_sR)/(kF*kC/(kF+kC)+kR)+zR*m_sR+h_uR*m_uR)*Ay/TWr;
end

function maxes = findMaxes(values)
    maxVals = zeros(2,size(values,2));
    for i = 1:size(values, 1)
        for j = 1:size(values,2)
            if values(i,j) > maxVals(1,j)
                maxVals(1,j) = values(i,j);
            elseif values(i,j) < maxVals(2,j)
                maxVals(2,j) = values(i,j);
            end
        end
    end
    maxes = array2table(maxVals);
end