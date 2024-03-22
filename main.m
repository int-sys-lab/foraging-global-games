clear; close all; clc;

%runtime options
PLOT = FALSE; %set to false to hide real-time plots
SAVE = false; %set to false to not save the data

N = 12; %number of agents
F = 6; %number of food items


h = 5; %5 %sensing horizon (m)
r = 0.25; %robot radius (m)
Ro = 30; %domain size (m)
Ri = 5; %nest size (m)

speed = 1; %3 %bot speed (m/s)
mem_err = 0.35;% 0.2; %memory error scaling factor (m)

J_max = 100; %maximum units of energy (J)
J = 0.5*J_max; % current energy, not 100% to avoid waiting (J)
food_energy = 5; %5; %energy per food item (J)

dJ = 0.1; %0.1; %rate of energy drain (J/sec)
bot_drain = 0.01; %0.01 %0.1; %energy drain from moving (J/s)

%%%% parameters for the utility function
ca = 1;
k = 0.25; %cost constant offset
lambda = 5.5; %cost exponential term


%%%% simulaiton time
tf = 1800; %final time, seconds
nt = 1800; %360; %number of time steps

%%% human recruitment event
tstep_remove = round(nt/2);
N_remove = floor(N/2); %out of 12


%initialize the agents' position and mode
bots = struct([]);
for i = 1:N
    rad = (Ri-2*r);
    th = (i/N)*2*pi;
    bots(i).px = rad*cos(th);
    bots(i).py = rad*sin(th);
    bots(i).mode = 0;
    bots(i).defecit = 0;
    %uniformly spaced thresholds not including 1
    bots(i).thresh = 1 - 0.25*(i/N);
    %bots can store one coordinate in memory
    bots(i).memory = [];
end

%initialize the location of the food tiems
pF = zeros([F,2]);
dth = 2*pi/F; %put each food item in its own slice
for i = 1:F
    rF = sqrt(rand)*(Ro-2*h) + 2*h;
    thF = rand*dth + (i-1)*dth;
    pF(i,:) = [rF*cos(thF), rF*sin(thF)];
end

%plot the initial condition of the system
plot_env(bots, N, Ri, Ro, pF, h, r);


%%%%% simluation!
t = linspace(0, tf, nt);
dt = t(2) - t(1);
J_hist = zeros([1, nt]); %historical energy of nest
D_hist = J_hist; %historical deficit sum
N_hist = J_hist; %historical number of active robots
for i = 1:nt
    %store energy history
    J_hist(i) = J;
    %remove N robots at timestep tstep_remove
    if i == tstep_remove
        N = N - N_remove;
    end
    %reset foraging counter
    n_foraging = 0;
    %count number of bots already foraging
    J_bots = 0; %energy defecit of the bots
    for n = 1:N
        if bots(n).mode > 0
            n_foraging = n_foraging + 1;
        end
        J_bots = J_bots + bots(n).defecit;
    end
    theta = 1 - (J - J_bots)/J_max;
    %calculate utility for all agents (homogeneous)
    m_u = -ca + k + lambda*exp(-(n_foraging + 1))+ theta;
    %select probability of any agent foraging this time step
    if m_u <= 0
        p = 0;
    else
        p = (-log((ca-k-theta) / lambda) - n_foraging - 1)/(N - n_foraging);
    end
    %update bots dynamics based on their current mode
    for n = 1:N
        v = [0; 0];
        switch bots(n).mode
            %%%%%% Here the bot is idling
            case 0 %do nothing
                if rand <= p
                   bots(n).mode = 1;
                end
            %%%%% Here the bot is foraging %%%%%
            case 1 %forage
                N_hist(i) = N_hist(i) + 1;
                %check how close the nearest food is
                df = sqrt(sum((pF - [bots(n).px, bots(n).py]).^2, 2));

                %pick theta to find food
                if any(df <= h) %move toward any we can see
                    idx = find(df == min(df));
                    th = atan2(pF(idx,2) - bots(n).py, pF(idx,1) - bots(n).px);
                    %forget the last fod we saw
                    bots(n).memory = [];
                elseif ~isempty(bots(n).memory) %otherwise, go to memory
                    th = atan2(bots(n).memory(2) - bots(n).py, ...
                               bots(n).memory(1) - bots(n).px);
                else %otherwise, add a new random memory position nearby
                    th = rand*(2*pi);
                    mem = [bots(n).px; bots(n).py] + h*[cos(th); sin(th)];
                    %project it onto the domain
                    if norm(mem) > Ro
                        mem = mem/norm(mem)*Ro;
                    end
                    %don't search your home for food
                    if norm(mem) < Ri
                        mem = mem/norm(mem)*Ri;
                    end
                    %update memory and move
                    bots(n).memory = mem;
                    th = atan2(mem(2), mem(1));
                end
                %move with max_speed at angle theta
                v = speed*[cos(th); sin(th)];

                %if we can pick up any food, move there instead
                if any(df <= r + speed*dt)
                    idx = find(df == min(df));
                    v = (pF(idx,:)' - [bots(n).px; bots(n).py]) / dt;
                    %add a noisy memory of the food location
                    dist = norm(pF(idx,:));
                    noisy_pf = pF(idx,:) + randn([1, 2])*dist*mem_err;
                    if norm(noisy_pf) > Ro
                        noisy_pf = noisy_pf/norm(noisy_pf) * Ro;
                    end
                    bots(n).memory = noisy_pf;
                    %change mode to "return to nest"
                    bots(n).mode = 2;
                else
                    %if we are close to the memory location, go to it
                    if ~isempty(bots(n).memory)
                        dM = [bots(n).memory(1) - bots(n).px; ... 
                              bots(n).memory(2) - bots(n).py];
                        if  norm(dM) <= speed*dt
                            th = atan2(dM(2), dM(1));
                            v = norm(dM)/dt * [cos(th); sin(th)];
                            bots(n).memory = [];
                        end
                    end
                end
            %%%%% end of foraging case! %%%%%
            case 2 %point toward the nest
                N_hist(i) = N_hist(i) + 1;
                th = pi + atan2(bots(n).py, bots(n).px);
                v = speed * [cos(th); sin(th)];
                %drop the food off at the nest if we are close enough
                if sqrt(bots(n).px^2 + bots(n).py^2) <= Ri
                    %energy exchange
                    J = min(J_max, J + food_energy - bots(n).defecit);
                    bots(n).defecit = 0;
                    
                    %pick whether to forage again or to become idle
                    %%% we re-calculate pi and mu using n-1 to not count
                    %%% ourself
                    m_u_self = -ca + k + lambda*exp(-(n_foraging))+ theta;
                    %select probability of any agent foraging this time step
                    if m_u_self <= 0
                        p_self = 0;
                    else
                        p_self = (-log((ca-k-theta) / lambda) - n_foraging - 2)/(N - n_foraging - 1);
                    end
                     if rand <= p_self
                         bots(n).mode = 1;
                     else
                         bots(n).mode = 0;
                         bots(n).memory = [];
                     end
                end
        end
        %enforce safety constraints and move
        v = CBF(bots, n, v, speed, r, Ro);
        bots(n) = move_bot(bots, n, v, dt, bot_drain);
        %update stored variables for each bot
        D_hist(i) = D_hist(i) + bots(n).defecit;
    end

    %energy dynamics
    J = J - dJ*dt;

    %exit sim if total system energy reaches 0
    if J - D_hist(i) <= 0
        break;
    end

    %plot at the end of the time step
    if PLOT
        plot_env(bots, N, Ri, Ro, pF, h, r);
        title("J = " + J)
        %title("Step " + i + ", J = " + J)
        pause(0)
        %saveas(gcf, ['figs/sim_', num2str(i, '%04g'), '.png'])
    end
end

%% plot energy history graph
lw = 2.5;
figure(2); clf; hold on;
plot(t, J_hist-D_hist, '-.k', 'linewidth', lw/2);
plot(t, J_hist, 'r', 'linewidth', lw);

plot([tstep_remove, tstep_remove], [0, 100])


legend({"System Energy", "Nest Energy"}, 'location', 'south west')


xlabel('Time (sec)')
ylabel('Energy (J)')
axis([min(t), max(t), 0, J_max])
grid on; box on;

%saveas(gcf, "energy.png")


%% plot number of active robots
figure(3); clf; hold on;

%expected number of robots
EV = zeros(1, length(t));

TH = 1 - (J_hist - D_hist)/J_max;
bad_idx = find(TH >= ca - k);
good_idx = find(TH < ca-k);

EV(bad_idx) = N + N_remove;

EV(good_idx) = -log((ca-k-TH(good_idx)) / lambda) - 1;

EV(1:tstep_remove) = min(EV(1:tstep_remove), N + N_remove);
EV(tstep_remove:end) = min(EV(tstep_remove:end), N);

%expected number of robots
%plot(t, EV, '-.', 'linewidth', lw);
area(t, EV, 'facecolor', [0.8, 0.8, 1.0])

%number of robots
plot(t, N_hist, '-k', 'linewidth', lw)
plot([0, tstep_remove, tstep_remove, t(end)], [N+N_remove, N+N_remove, N, N], ...
    'r', 'linewidth', lw)

legend({"Expected", "Actual", "Maximum"})

xlabel('Time (sec)')
ylabel('Active Robots')
axis([min(t), max(t), 0, N+N_remove])
grid on; box on;

%saveas(gcf, "num_robots.png")

%% save data
if SAVE
    data = [n, length(t), t, J_hist, D_hist, N_hist];
    writematrix(data, 'sim_results.csv','WriteMode','append' )
end

%% helper function
function plot_env(bots, N, Ri, Ro, pF, h, r)
    figure(1); clf; hold on;
    %draw the environment
    rectangle("Position", [-Ro, -Ro, 2*Ro, 2*Ro], "Curvature", 1, ...
        "FaceColor", [0.8, 1.0, 0.8]);
    rectangle("Position", [-Ri, -Ri, 2*Ri, 2*Ri], "Curvature", 1, ...
        "LineStyle", "--", "FaceColor", 'w');
    %draw the position of the robots
    bp = [bots.px; bots.py];
    plot(bp(1,1:N), bp(2,1:N), 'ok', 'MarkerFaceColor', 'r')
    for n = 1:N %(we may plot fewer than N bots)
        %draw the robot actual size
        rectangle("Position", [bots(n).px-r, bots(n).py-r, 2*r, 2*r], ...
            "Curvature", 1, "FaceColor", "k")
        %%% draw the sensing distance
        rectangle("Position", [bots(n).px-h, bots(n).py-h, 2*h, 2*h], ...
            "Curvature", 1, "LineStyle", ":")
        if ~isempty(bots(n).memory)
            plot(bots(n).memory(1), bots(n).memory(2), 'kx');
            plot([bots(n).memory(1), bots(n).px], ...
                 [bots(n).memory(2), bots(n).py], ':k');
        end
    end
    %draw the position of the food
    plot(pF(:,1), pF(:,2), 'ok', 'MarkerFaceColor', 'b');
    axis([-Ro, Ro, -Ro, Ro]);
    axis square;
end

%% helper function
function v = CBF(bots, n, vstar, max_speed, r, Ro)
    %linear CBF coefficient
    coef = 0.5;
    N = length(bots);
    %construct the CBF
    A = zeros(N, 2); B = zeros([N,1]);
    %collision avoidance
    for k = [1:n-1, n+1:N]
        dp = [bots(k).px - bots(n).px; ...
              bots(k).py - bots(n).py];
        A(k,:) = 2*(dp)';
        B(k) = coef * (-((2*r)^2 - dot(dp, dp)) - norm(dp)*max_speed/2);
    end
    %stay in the arena
    pn = [bots(n).px; bots(n).py];
    A(N,:) = 2*pn';
    B(N) = (Ro-r)^2 - dot(pn,pn);

    %solve the OCP
    H = [1, 0; 0, 1];
    f = -2*vstar';
    
    opts = optimoptions("quadprog", "display", "none");
    v = quadprog(H, f, A, coef*B, [], [], [], [], vstar, opts);

    if norm(v) > max_speed
        v = v/norm(v)*max_speed;
    end
end


%% helper function to move the bot and handle energy updates
function mybot = move_bot(bots, n, v, dt, bot_drain)

    mybot = bots(n);

    mybot.px = mybot.px + v(1)*dt; % speed*cos(th)*dt;
    mybot.py = mybot.py + v(2)*dt; %speed*sin(th)*dt;
    %drain energy while moving
    mybot.defecit = mybot.defecit + bot_drain*norm(v)*dt;
end
