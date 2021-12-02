%% SETUP
clear;
close all

% SIMULATION MANAGEMENT
nrunsTotal = 1; % # of repeat simulations in this file
batchNumber = 1; % file number
rng(81368); % initial random seed

% POTENTIAL SWITCH + TURN ON/OFF MFPT MODE
harmonicPotential = 1; % turn harmonic potential on or off
MFPT_on = 1; % if off, will run for fixed duration; if on, will run until lowest ε is crossed

% TURN TESTS ON AND OFF 
meanSquaredDisplacement = 0; % run analysis of mean squared displacement
initializationDistribution = 0; % check whether particles are initialized in proper Boltzmann distribution
harmonicDistribution = 0; % check whether particles exist in proper Boltzmann distribution during simulation

    %% SIMULATION PARAMETERS

% TIME PARAMETERS
dt = 10^-1; % seconds; timestep duration

if MFPT_on==1
    epsilons = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3]./10^9; % meters; reaction boundaries
else
    duration = 30; % seconds; simulation duration
    nsteps = duration/dt; % number of steps
end

% SYSTEM PARAMETERS
dimension = 2; % dimension of system
D_s = 18.1; % nanometers^2/second; translational diffusion coefficient 
T = 298; % Kelvin; temperature

% CONSTANTS
N_a = 6.02214076*10^23; % 1/mol; Avogadro's Number
D = D_s / 10^18; % meters^2/second; Translational diffusion coefficient -- converted from above
k_b = 1.380649*10^-23; % kilogram meters^2 second^-2 Kelvin^-1; Boltzmann Constant
kT = k_b*T; % Joules; Thermal energy of the system 
gamma = kT/D; % Joule second / meters^2; drag coefficient, automatically calculated from the desired diffusion coefficient
numAtoms = 4; % number of atoms in system -- ***2*pairs if MFPT_on==1***

% HARMONIC POTENTIAL
if harmonicPotential==1
    sigma = 3.5 / 10^9; % meters; standard deviation of harmonic potential
    separation = 7 / 10^9; % meters; separation of harmonic centers of the 2 walkers
    k = kT/sigma^2; % Joules / meters^2; stiffness of harmonic potential
end    

%% RUN SIMULATION

if harmonicPotential==1
    centers = zeros(numAtoms,dimension); % meters; position of harmonic centers -- each row is the XY position of each atom
    centers(1:numAtoms/2,1)= 0; % first atom's center is at (0,0)
    centers(numAtoms/2+1:numAtoms,1) = separation; % second atom's center is at (separation,0)
end

nruns = 0; % counter for number of simulations done

FPT_runs = []; % keep track of First Passage Times in each run, used only if MFPT_on==1
firstFrames = []; % save initial positions of atoms in each run
lastFrames = []; % save final positions of atoms in each run

xyInits = []; % used only if initializationDistribution==1
xyDists = []; % used only if harmonicDistribution==1
MSDs = {}; % used only if meanSquaredDisplacement==1

while nruns ~= nrunsTotal
    timeStart = tic; % start timer for this run

    MSD = []; % used only if meanSquaredDisplacement==1
    
    randomSeed = randi(1000000); % generate seed for random number generation
    rng(randomSeed); % produce random number generator from seed
    
    if harmonicPotential==1
        previousFrame1x = normrnd(0,sigma,numAtoms/2,1); % draw initial x of 1st atom from normal distribution of μ=0 and σ=sigma
        previousFrame2x = normrnd(separation,sigma,numAtoms/2,1); % draw initial x of 2nd atom from normal distribution of μ=separation and σ=sigma
        previousFrame1y = normrnd(0,sigma,numAtoms/2,1); % draw initial y of 1st atom from normal distribution of μ=0 and σ=sigma
        previousFrame2y = normrnd(0,sigma,numAtoms/2,1); % draw initial y of 2nd atom from normal distribution of μ=0 and σ=sigma
        previousFrame = [previousFrame1x previousFrame1y; previousFrame2x previousFrame2y]; % make first frame from these initial positions
    elseif harmonicPotential==0
        previousFrame = (-1 + 2.*rand(numAtoms,dimension)).*100/10^9; % randomly distribute numbers from -100 to 100 nm in both directions
    end    
    
    firstFrame = previousFrame;
    firstFrames = [firstFrames; firstFrame]; % save first frame
    
    if initializationDistribution==1
        xyInits = [xyInits; firstFrame-centers];
    end
    
    if MFPT_on==1
        FPTs = ones(numAtoms/2,size(epsilons,2))*-1; % initialize First Passage Time tracker, with one entry for each boundary value ε
    end
    
    % SET UP INTEGRATION
    step=2;
    run=1;
    
    while run==1
        
        if harmonicPotential == 1
            gradU = k.*(previousFrame-centers);
        else
            gradU = 0;
        end
        
        % Brownian Dynamics Equation of Motion; https://en.wikipedia.org/wiki/Brownian_dynamics
        currentFrame = previousFrame - 1/gamma*dt.*gradU + ...
            sqrt(2*gamma*kT*dt)/gamma.*normrnd(zeros([numAtoms dimension]),ones([numAtoms dimension])); 
        
        previousFrame = currentFrame; % update moving frame

        if MFPT_on==1
            % for loop that checks whether the atoms have passed a boundary for the first time
            % marks time as step*dt if it has
            
            %distance = norm(currentFrame(1,:)-currentFrame(2,:)); 
            distances = sqrt(sum((currentFrame(1:numAtoms/2,:)-currentFrame(numAtoms/2+1:numAtoms,:)).^2,2)); % calculate distance between atom pairs
            for pair_index = 1:numAtoms/2
                for epsilon_index = 1:size(epsilons,2)
                   epsilon=epsilons(epsilon_index);
                   if distances(pair_index)<= epsilon && FPTs(epsilon_index)==-1
                      FPTs(pair_index,epsilon_index) = step*dt;
                   end   
                end
            end

            % terminates run if lowest boundary has been passed
            if FPTs(1) ~= -1
                run = 0;
            end
        else
            if step == nsteps
                run = 0;
            end
        end
        
        % update step counter to measure time
        step = step+1;
        
        % gather data for tests
        if meanSquaredDisplacement==1
            MSD = [MSD; step*dt mean(sqrt(sum((currentFrame-firstFrame).^2,2)).^2)];
        end
        
        if harmonicPotential==1
            xyDists = [xyDists; currentFrame-centers];
        end
        
    end
    
    nruns = nruns+1; % update run counter
    
    if MFPT_on==1
        FPT_runs = [FPT_runs; FPTs]; % save First Passage Times of this run
    end
    
    lastFrames = [lastFrames; previousFrame]; % save final positions of atoms

    if meanSquaredDisplacement==1
        MSDs{end+1} = MSD;
    end
    
    timeElapsed = toc(timeStart); % stop timer for this run
    disp(timeElapsed) % display timer output for this run
end

% OUTPUT RESULTS AND SAVED FRAMES OF SIMULATION
writematrix(FPT_runs,strcat(num2str(batchNumber),"FPTs.csv"));
writematrix(firstFrames,strcat(num2str(batchNumber),"first.csv"));
writematrix(lastFrames,strcat(num2str(batchNumber),"last.csv"));

%% RUN TESTS
% set up graphs and data structures for tests
if meanSquaredDisplacement == 1 % create array and graph to track MSD
    figure(1)

    for k = 1:length(MSDs)
       MSD_graph = MSDs{k};
       title("Mean Squared Displacement vs. Time")
       xlabel("Time (s)")
       ylabel("Mean Squared Displacement (meters^2)")

       plot(MSD_graph(:,1),MSD_graph(:,2),"DisplayName",num2str(k))
       hold on

    end
    
    diffusion_law_string = strcat("MSD = ",num2str(2*dimension),"Dt");
    plot(MSD_graph(:,1),2*dimension*D.*MSD_graph(:,1),'DisplayName',diffusion_law_string)
    hold on
    legend("Location","northwest")
    
end

if initializationDistribution == 1 % check whether particles are initialized in proper Boltzmann distribution
    figure(2)
    histogramValues = xyInits(:);
    h = histogram(histogramValues,"Normalization","pdf");
    distanceRange = linspace(min(histogramValues),max(histogramValues),1000);
    ideal =  sqrt(k ./ kT ./ (2*pi)) .* exp(-1*k/kT/2.*(distanceRange.^2));
    hold on
    plot(distanceRange,ideal)
    legend(["Simulation","√(^{k}/_{k_b*T*2π}) * e^{^{-k}/_{k_b*T*2}*x^2}"])
    title("Initialization Distribution")
    xlabel("Distance from Harmonic Centers (m)")
    ylabel("Probability")
end

if harmonicDistribution == 1 % check whether particles exist in proper Boltzmann distribution during simulation
    figure(3)
    histogramValues = xyDists(:);
    h = histogram(histogramValues,"Normalization","pdf");
    distanceRange = linspace(min(histogramValues),max(histogramValues),1000);
    ideal =  sqrt(k ./ kT ./ (2*pi)) .* exp(-1*k/kT/2.*(distanceRange.^2));
    hold on
    plot(distanceRange,ideal)
    legend(["Simulation","√(^{k}/_{k_b*T*2π}) * e^{^{-k}/_{k_b*T*2}*x^2}"])
    title("Harmonic Distribution")
    xlabel("Distance from Harmonic Centers (m)")
    ylabel("Probability")
end

