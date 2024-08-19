clc;
clear;
close all;

% Function to clear the screen
if ispc
    system('cls');
else
    system('clear');
end

% Wait for a moment before showing new content
pause(1);

% Network Parameters
xm = 100;%store the x and y coordinates of the centroid of a cluster.
ym = 100;%store the x and y coordinates of the centroid of a cluster.
n = 12;% Number of nodes
sinkx = 55;%
sinky = 55;%store the x and y coordinates of the sink node in a network.
Eo = 0.5;
Eelec = 50 * 10^(-9);
ETx = 50 * 10^(-9);
ERx = 50 * 10^(-9);
Eamp = 100 * 10^(-12);
EDA = 5 * 10^(-9);
Esens = 5 * 10^(-9);
k = 3000;
p = 0.5;
No = p * n;
rnd = 0;
operating_nodes = n;
transmissions = 0;
temp_val = 0;
flag1stdead = 0;
max_distance = 30;
min_energy_threshold = 0.01;
dead_nodes_per_round = [];

% Membaca data lokasi dari file Excel
%filename = 'lokasi_node.xlsx';
%data_lokasi = readtable(filename);

% Pastikan jumlah node sesuai dengan data lokasi di Excel
%n = height(data_lokasi);

% Initialize WSN Nodes
for i = 1:n
    SN(i).id = i;
    SN(i).x = rand() * xm;
    SN(i).y = rand() * ym;
    SN(i).chx = 0;
    SN(i).chy = 0;
    SN(i).E = Eo;
    SN(i).role = 0;
    SN(i).cluster = 0;
    SN(i).cond = 1;
    SN(i).rop = 0;
    SN(i).rleft = 0;
    SN(i).dtch = 0;
    SN(i).dts = 0;
    SN(i).tel = 0;
    SN(i).rn = 0;
    SN(i).chid = 0;
    SN(i).active = true;
    SN(i).received_adv = false;
    SN(i).is_isol = false;
end
   
    
CL = struct([]);
total_energy_per_round = [];
avg_node = [];
%clusterInfoPerRound = struct([]);
someEnergyThreshold = 0.01;
nrg = [];
tr = [];
op = [];
cluster_counts_per_round = [];
first_isolated_round = [];
isolated_nodes_per_round = [];


% Inisialisasi
rnd = 0;
operating_nodes = n;
transmissions = 0;
cluster_id = 1;
dead_nodes=0;
clusterInfoPerRound = {};

% Initialize energy tracking for each node
energy_per_node_per_round = Eo * ones(n, 1);

while operating_nodes > 0
    rnd = rnd + 1;
    t = p / (1 - p * mod(rnd, 1/p));
    tleft = mod(rnd, 1/p);
    CLheads = 0;
    energy = 0;
    clusterInfoThisRound = struct('x', {}, 'y', {}, 'id', {}, 'cluster', {}, 'energy', {}, 'isolated_nodes', {}, 'CHcoord', {});


% Pemilihan Cluster Head
CLheads = 0;
energy = 0;

% Tentukan jumlah maksimum CH yang diizinkan (misalnya 3)
maxClusters = 3;

for i = 1:n
    SN(i).cluster = 0;
    SN(i).role = 0;
    SN(i).chid = 0;
    SN(i).received_adv = false; % Reset received_adv sebelum pemilihan cluster head
    if SN(i).rleft > 0
        SN(i).rleft = SN(i).rleft - 1;
    end
    if SN(i).E > 0 && SN(i).rleft == 0 && CLheads < maxClusters
        generate = rand();
        if generate < t
            SN(i).role = 1;
            SN(i).chx = SN(i).x;
            SN(i).chy = SN(i).y;
            SN(i).rn = rnd;
            SN(i).tel = SN(i).tel + 1;
            SN(i).rleft = 1 / p - tleft;
            SN(i).dts = sqrt((sinkx - SN(i).x)^2 + (sinky - SN(i).y)^2);
            CLheads = CLheads + 1;
            SN(i).cluster = CLheads;
            SN(i).chid = i;
            cluster_info = struct('x', SN(i).x, 'y', SN(i).y, 'id', i, ...
                                  'cluster', CLheads, 'energy', SN(i).E, ...
                                  'isolated_nodes', [], 'CHcoord', [SN(i).x, SN(i).y]);

            CL = [CL, cluster_info];
            clusterInfoThisRound = [clusterInfoThisRound, cluster_info];
        end
    end
end

    % Pengelompokan node ke cluster dan menghitung jarak
    for i = 1:n
        if SN(i).role == 0 && SN(i).E > 0 && CLheads > 0
            d = arrayfun(@(m) sqrt((CL(m).x - SN(i).x)^2 + (CL(m).y - SN(i).y)^2), 1:CLheads);
            [min_d, min_d_index] = min(d);
            SN(i).cluster = min_d_index;
            SN(i).dtch = min_d;
            SN(i).chid = CL(min_d_index).id;
            SN(i).chx = CL(min_d_index).x;
            SN(i).chy = CL(min_d_index).y;
        end
    end

    % Kirim ADV dan tandai node terisolasi
    [SN, distances] = send_adv_to_cm(find([SN.role] == 1), SN, max_distance);
    [SN,isolated_nodes] = mark_isolated_nodes(SN, max_distance, min_energy_threshold);
    isolated_nodes_per_round = [isolated_nodes_per_round, isolated_nodes];

    % Update informasi kluster
    if CLheads > 0
        for j = 1:numel(clusterInfoThisRound)
            cluster_nodes = SN([SN.cluster] == clusterInfoThisRound(j).cluster);
            clusterInfoThisRound(j).isolated_nodes = [cluster_nodes([cluster_nodes.is_isol]).id];
            clusterInfoThisRound(j).memberNodes = [cluster_nodes.id];
            clusterInfoThisRound(j).memberCoords = [[cluster_nodes.x]' [cluster_nodes.y]'];  % Simpan koordinat
            clusterInfoThisRound(j).isIsolated = [cluster_nodes.is_isol];  % Simpan status isolasi
        end
        clusterInfoPerRound{rnd} = clusterInfoThisRound;  % Assign data to cell array
    else
        clusterInfoPerRound{rnd} = struct('x', {}, 'y', {}, 'id', {}, 'cluster', {}, 'energy', {},'isolated_nodes',{}, 'CHcoord', {});
    end
    
    %###### Debuging #########
    % Cek apakah ronde terisolasi pertama sudah teridentifikasi
    if ~isempty(isolated_nodes) && isempty(first_isolated_round)
        first_isolated_round = rnd;
    end

    %%Cetak informasi node terisolasi untuk 10 ronde pertama
    if rnd <= 10 
        fprintf('Round %d: Total CH selected = %d\n', rnd, CLheads);
            if rnd <= length(clusterInfoPerRound{rnd})
                for j = 1:numel(clusterInfoPerRound{rnd})  % Akses langsung ke data ronde saat ini
                    cluster = clusterInfoPerRound{rnd}(j);
                    fprintf('  Cluster %d with CH ID %d at position (%.2f, %.2f)\n', ...
                        cluster.id, cluster.id, cluster.CHcoord(1), cluster.CHcoord(2));


            if ~isempty(cluster.isolated_nodes)
                fprintf('  Isolated nodes in Cluster %d:\n', cluster.id);
                for isolated_id = cluster.isolated_nodes
                    node = SN(isolated_id);
                    fprintf('    Node ID: %d (X: %.2f, Y: %.2f, Isolated: %s)\n', ...
                            node.id, node.x, node.y, mat2str(node.is_isol));
                end
            else
                fprintf('  No isolated nodes in Cluster %d.\n', cluster.id);
            end
                end
            else
                fprintf('Warning: rnd index exceeds the number of recorded rounds.\n');
            end
    end

    
    % Fase Steady-State: Transmisi CM
    for i = 1:n
        if SN(i).cond == 1 && SN(i).role == 0 && CLheads > 0
            ETx = Eelec * k + Eamp * k * SN(i).dtch^2;
            ERx = Eelec * k;
            Econs = ETx + ERx;
            if SN(i).E > Econs
                if SN(i).is_isol
                    SN(i).E = SN(i).E; % No energy consumed for isolated nodes
                    operating_nodes = operating_nodes - 0.01; % Decrease slightly for simulation
                    dead_nodes = dead_nodes + 0.01;
                else
                    SN(i).E = SN(i).E - Econs;
                    energy = energy + Econs;
                end
            else
                dead_nodes = dead_nodes + 1;
                operating_nodes = operating_nodes - 1;
                SN(i).active = false;
                SN(i).cond = 0;
                SN(i).chid = 0;
                SN(i).rop = rnd;
            end
        end
    end

    % Fase Steady-State: Transmisi CH
    for i = 1:n
        if SN(i).cond == 1 && SN(i).role == 1
            ETx = (Eelec + EDA) * k + Eamp * k * SN(i).dts^2;
            ERx_total = sum((Eelec + EDA) * k * ([SN.cluster] == i + 1 & [SN.role] == 0));
            Econs = ETx + ERx_total;
            if SN(i).E > Econs
                SN(i).E = SN(i).E - Econs;
                energy = energy + Econs;
            else
                dead_nodes = dead_nodes + 1;
                operating_nodes = operating_nodes - 1;
                SN(i).active = false;
                SN(i).cond = 0;
                SN(i).rop = rnd;
            end
        end
    end

    % Update metrics
    tr = [tr, operating_nodes];
    op = [op, operating_nodes];
    avg_energy_per_node = energy / max(operating_nodes, 1);
    avg_node = [avg_node, avg_energy_per_node];
    total_energy_consumed_this_round = sum(Eo - [SN.E]);
    total_energy_per_round = [total_energy_per_round, total_energy_consumed_this_round];
    dead_nodes_per_round = [dead_nodes_per_round, dead_nodes];
    if energy > 0
        nrg = [nrg, energy];
        transmissions = transmissions + 1;
    end

    % Check conditions
    if operating_nodes < n && temp_val == 0
        temp_val = 1;
        flag1stdead = rnd;
    end

    if operating_nodes == 0
        break;
    end
end



% Di luar loop simulasi
if flag1stdead > 0
    fprintf('Node pertama mati pada ronde %d.\n', flag1stdead);
else
    fprintf('Tidak ada node yang mati selama simulasi.\n');
end
fprintf('Length of clusterInfoPerRound: %d\n', numel(clusterInfoPerRound));
%============================== Ploting Node Distribution ============
% Pemanggilan fungsi dan plot berdasarkan ronde yang dipilih
isolated_nodes = [];
for idx = 1:numel(SN)
    if SN(idx).is_isol
        isolated_nodes = [isolated_nodes, SN(idx).id];
    end
end

%display_isolated_info(clusterInfoPerRound);
[cluster_info, isolated_nodes, chosen_round_number] = display_isolated_info(clusterInfoPerRound);
if ~isempty(cluster_info)
    plot_cluster_distribution(cluster_info, isolated_nodes, sinkx, sinky, chosen_round_number);
end

% Plot Energy consumed per Transmission
if ~isempty(nrg)  % Pastikan nrg tidak kosong
    figure;
    plot(1:transmissions, nrg, '-b', 'LineWidth', 1);
    xlabel('Transmission');
    ylabel('Energy (Joule)');
    title('Energy Consumed per Transmission');
    grid on;
end
% Plot Average Energy consumed by a Node per Transmission
figure;
plot(1:transmissions, avg_node(1:transmissions), '-r', 'LineWidth', 1);
xlabel('Transmission');
ylabel('Energy (Joule)');
title('Average Energy Consumed by a Node per Transmission');
grid on;

% Plot Total Energy Consumed by the WSN per Transmission
figure;
plot(1:transmissions, total_energy_per_round(1:transmissions), '-g', 'LineWidth', 1);
xlabel('Transmission');
ylabel('Energy (Joule)');
title('Total Energy Consumed by the WSN per Transmission');
grid on;

% Plot Jumlah node mati per round
figure;
plot(1:rnd, op(1:rnd), '-m', 'LineWidth', 1);
xlabel('Round');
ylabel('Number of Operating Nodes');
title('Number of Operating Nodes per Round');
grid on;

% Plot Jumlah Node yang aktif atau beroperasi tiap round
figure;
plot(1:rnd, tr(1:rnd), '-c', 'LineWidth', 1);
xlabel('Round');
ylabel('Number of Active Nodes');
title('Number of Active Nodes per Round');
grid on;

% Mendapatkan jumlah total ronde dari panjang list dead_nodes_per_round
rounds = 1:length(dead_nodes_per_round);

% Membuat plot
figure;
plot(rounds, dead_nodes_per_round, 'red', 'LineStyle', '-', 'LineWidth', 2);
title('Dead Nodes vs Simulation Rounds');
xlabel('Simulation Rounds');
ylabel('Number of Dead Nodes');
grid on;



function [cluster_info, isolated_nodes, round_number] = display_isolated_info(clusterInfoPerRound)
    round_number = input('\nEnter the round number to display isolated node information: ');
    fprintf('Checking round %d\n', round_number);
    isolated_nodes = [];
    cluster_info = [];

    if round_number > 0 && round_number <= numel(clusterInfoPerRound) && ~isempty(clusterInfoPerRound{round_number})
        cluster_info = clusterInfoPerRound{round_number};
        fprintf('Total clusters to check: %d\n', numel(cluster_info));

        isolated_info = sprintf('Round %d - Isolated Node Details:\n', round_number);
        total_isolated = 0;

        for cluster = cluster_info
            fprintf('Cluster %d with CH ID %d at position (%.2f, %.2f):\n', ...
                    cluster.id, cluster.id, cluster.CHcoord(1), cluster.CHcoord(2));

            if ~isempty(cluster.isolated_nodes)
                fprintf('  Isolated nodes in Cluster %d:\n', cluster.id);
                for node_id = cluster.isolated_nodes
                    fprintf('    Node ID: %d (Isolated)\n', node_id);
                    total_isolated = total_isolated + 1;
                end
            else
                fprintf('  No isolated nodes in Cluster %d.\n', cluster.id);
            end
        end

        isolated_info = [isolated_info sprintf('Total Isolated Nodes: %d\n', total_isolated)];
        disp(isolated_info);
    else
        disp('No clusters formed or round number out of range.');
    end

    return;
end


function plot_cluster_distribution(cluster_info, isolated_nodes, sinkx, sinky, round_number)
    fprintf('Debugging node distribution for Round %d:\n', round_number);
    figure;
    hold on;
    grid on;

    % Plot sink with a more distinctive symbol
    scatter(sinkx, sinky, 150, 'k', 'd', 'filled', 'DisplayName', 'Sink'); % Sink marked with a filled black diamond

    colors = ['b', 'g', 'r'];  % Colors for each cluster
    legendEntries = containers.Map();  % To avoid legend duplication

    % Plot CH and CM from clusterInfo data
    for idx = 1:length(cluster_info)
        cluster = cluster_info(idx);
        ch_x = cluster.CHcoord(1);
        ch_y = cluster.CHcoord(2);
        color = colors(mod(idx - 1, length(colors)) + 1);

        % Plotting CH without specific cluster labels
        scatter(ch_x, ch_y, 70, color, '^', 'filled', 'DisplayName', sprintf('CH %d', cluster.id));

        % Check if cluster data includes member nodes locations and isolated status
        if isfield(cluster, 'memberCoords') && isfield(cluster, 'isIsolated')
            for j = 1:size(cluster.memberCoords, 1)
                cm_x = cluster.memberCoords(j, 1);
                cm_y = cluster.memberCoords(j, 2);
                is_isolated = cluster.isIsolated(j);

                % Decide the marker based on the isolation status
                if is_isolated
                    marker = 'x';
                else
                    marker = 'o';
                end
                scatter(cm_x, cm_y, 50, color, marker);
            end
        end
    end

    title(sprintf('Node Distribution at Round %d', round_number));
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    legend('show', 'Location', 'eastoutside');  % Place legend outside the graph on the right side
    hold off;
end

%====================end Distribution=======================


% Define functions for isolated node identification and energy calculations
function [nodes, distances] = send_adv_to_cm(cluster_head_indices, nodes, max_distance)
    distances = [];
    for idx = 1:numel(cluster_head_indices)
        ch_id = cluster_head_indices(idx);
        %disp(['Cluster head index: ', num2str(ch_id)]);
        if nodes(ch_id).active
            nodes(ch_id).received_adv = true;
            for k = 1:numel(nodes)
                if k ~= ch_id && nodes(k).active && nodes(k).cluster == nodes(ch_id).cluster
                    distance = norm([nodes(ch_id).x, nodes(ch_id).y] - [nodes(k).x, nodes(k).y]);
                    %fprintf('Node %d is %f units from CH %d\n', k, distance, ch_id);
                    if distance <= max_distance
                        nodes(k).received_adv = true;
                        %fprintf('Node %d: Received ADV from CH %d\n', k, ch_id);
                    end
                end
            end
        end
    end
end


function [nodes,isolated_nodes] = mark_isolated_nodes(nodes, max_distance, min_energy_threshold)
    isolated_nodes = [];
    for i = 1:numel(nodes)
        if ~nodes(i).received_adv && nodes(i).active && nodes(i).E > min_energy_threshold
            nodes(i).is_isol = true;
            isolated_nodes = [isolated_nodes, nodes(i)];
        else
            nodes(i).is_isol = false;
        end
    end
end