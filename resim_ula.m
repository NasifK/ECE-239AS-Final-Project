%% visualize Faustini Crater --- make sure everything is importing properly
figure(1);
model = createpde;
importGeometry(model,"moon_10.stl");
pdegplot(model);
title("Faustini Crater Visualization");

%% create point cloud
% point cloud
stlData = stlread('moon_10.stl');
point_cloud = stlData.Points;
figure(2);
scatter3(point_cloud(:,1), point_cloud(:,2), point_cloud(:,3), 5, 'filled');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Point Cloud Visualization');

%% begin raytracing
fc = 2.4e9; % transmitter / recieving recieving at 2.9 GHz
lambda = physconst("lightspeed")/fc;
txArray = arrayConfig("Size",[1 1],"ElementSpacing",0.5*lambda);
rxArray = arrayConfig("Size",[1 1],"ElementSpacing",0.5*lambda);

% define tx and rx for our STL
tx = txsite("cartesian", ...
    "Antenna",txArray, ...
    "AntennaPosition",[100; 220; 22], ...
    'TransmitterFrequency',2.4e9);

rx = rxsite("cartesian", ...
    "Antenna",rxArray, ...
    "AntennaPosition",[210; 220; 30]);

% view transmitter reciever pair
siteviewer("SceneModel","moon_10.stl");
show(tx,"ShowAntennaHeight",true)
show(rx,"ShowAntennaHeight",true)

% define propogation model
pm = propagationModel("raytracing", ...
    "CoordinateSystem","cartesian", ...
    "Method","Image", ...
    "AngularSeparation","low", ...
    "MaxNumReflections",2, ...
    "SurfaceMaterial","custom", ...
    "SurfaceMaterialPermittivity",.297, ...
    "SurfaceMaterialConductivity",0.00848);

rays = raytrace(tx,rx, pm);


%% perform raytracing for each point in the point_cloud
cols = size(point_cloud, 2) + 1; % add one more column for path loss
path_losses_2 = zeros(size(point_cloud, 1), cols);

for i = 1:size(point_cloud, 1)
    % define tx and rx for our STL
    tx = txsite("cartesian", ...
        "Antenna",txArray, ...
        "AntennaPosition",[100; 220; 30], ...
        'TransmitterFrequency',2.4e9);
   
    rx = rxsite("cartesian", ...
    "Antenna",rxArray, ...
    "AntennaPosition",[point_cloud(i,1); point_cloud(i,2); point_cloud(i,3)+2]);

    rays = raytrace(tx,rx, pm);

    % view transmitter reciever pair
    %siteviewer("SceneModel","moon_10.stl");
    %show(tx,"ShowAntennaHeight",true)
    %show(rx,"ShowAntennaHeight",true)

    if numel(rays{1}) < 1
        path_loss = NaN;
    else
        path_loss = rays{1}.PathLoss;
    end

    if isvector(path_loss) && numel(path_loss) > 1
        total_path_loss = sum(path_loss);
    else
        total_path_loss = path_loss;
    end
    disp(total_path_loss);

    % add point to the path losses
    path_losses_2(i, 1) = point_cloud(i,1);
    path_losses_2(i,2) = point_cloud(i,2);
    path_losses_2(i,3) = point_cloud(i,3);
    path_losses_2(i,4) = total_path_loss;
end

%% visualize raytracing results

% Extracting the coordinates and path_loss values
x = path_losses_2(:,1);
y = path_losses_2(:,2);
z = path_losses_2(:,3);
path_loss = path_losses_2(:,4);

% Plotting the scatter plot
figure(4);
scatter3(x(:), y(:), z(:), 50, path_loss(:));
colorbar;
h = colorbar;
ylabel(h, 'Path Loss (dBm)');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Scatter Plot with Path Loss');

%% extending this to a 4 element ULA array
% can we extract an H matrix?
N_r = 4;
N_t = 4;
wavelength = 0.125; % for 2.4 GHz wave
path_losses = cell(size(point_cloud, 1), cols);
pm.MaxRelativePathLoss = Inf;

for i = 1:size(point_cloud, 1)
    H = zeros(N_t,N_r);
    for transmit_id = 1:N_t
         tx = txsite("cartesian", ...
        "Antenna",txArray, ...
        "AntennaPosition",[100; 220 + (0.5*lambda*transmit_id - 1.5); 30], ...
        'TransmitterFrequency',2.4e9);

        for recieve_id = 1:N_r
            rx = rxsite("cartesian", ...
            "Antenna",rxArray, ...
            "AntennaPosition",[point_cloud(i,1); point_cloud(i,2) + (0.5*wavelength*recieve_id - 1.5); point_cloud(i,3)+2]);
        
            rays = raytrace(tx,rx, pm);
            path_loss = rays{1}.PathLoss;
            path_loss = db2mag(path_loss);
                if isvector(path_loss) && numel(path_loss) > 1
                    total_path_loss = sum(path_loss) ./ 1000;
                else
                    total_path_loss = path_loss / 1000;
                end
            H(transmit_id, recieve_id) = total_path_loss;
        end
    end

    disp(H);
    % add matrix array to the path losses
    path_losses{i, 1} = point_cloud(i,1);
    path_losses{i,2} = point_cloud(i,2);
    path_losses{i,3} = point_cloud(i,3);
    path_losses{i,4} = H;
end

%% now that we have H arrays, lets see if MRC+MRT or Spatial Multiplexing should be used?

% lets plot in terms of snr
snr_res = 1000;
snr_db = linspace(-30,30, snr_res);
snr = db2mag(snr_db);

% mrc+mrc
mrc_spectrial_efficiencies = zeros(size(point_cloud, 1), cols);

avg_mrc_spectral_effieciencies = zeros(snr_res, 1);

for idx_snr = 1:snr_res
    
    for i = 1:size(point_cloud, 1)
        H = path_losses{i, 4};
        H = H.^(-1);
        [~, sigma, ~] = svd(H);
        P = [2, 0, 0, 0;0,2,0,0;0,0,2,0;0,0,0,2];
        
        sum = 0;
        for sum_idx = 1:N_t/2
            hi = log2(1 + (snr(idx_snr) / N_t) * P(sum_idx, sum_idx) * sigma(sum_idx, sum_idx)^2);
            sum = sum + hi;
        end
        
        mrc_spectrial_efficiencies(i, 1) = point_cloud(i,1);
        mrc_spectrial_efficiencies(i,2) = point_cloud(i,2);
        mrc_spectrial_efficiencies(i,3) = point_cloud(i,3);
        mrc_spectrial_efficiencies(i, 4) = sum;
    end
    
    avg_mrc_spectral_effieciencies(idx_snr) = mean(mrc_spectrial_efficiencies(:, 4));
end

%%
% spatial multiplexing

spatial_spectral_efficiencies = zeros(size(point_cloud, 1), cols);
avg_spatial_spectral_efficiencies = zeros(snr_res,1);

for idx_snr = 1:snr_res
    for i = 1:size(point_cloud, 1)
        H = path_losses{i, 4};
        H = H.^(-1);
        [~, sigma, ~] = svd(H);
        P = [4, 0,0,0;0,0,0,0;0,0,0,0;0,0,0,0];
    
        sum = 0;
        for sum_idx = 1:N_t/2
            hi = log2(1 + (snr(idx_snr) / N_t) * P(sum_idx, sum_idx) * sigma(sum_idx, sum_idx)^2);
            sum = sum + hi;
        end
    
        spatial_spectral_efficiencies(i, 1) = point_cloud(i,1);
        spatial_spectral_efficiencies(i,2) = point_cloud(i,2);
        spatial_spectral_efficiencies(i,3) = point_cloud(i,3);
        spatial_spectral_efficiencies(i, 4) = sum;
    end
    
    avg_spatial_spectral_efficiencies(idx_snr) = mean(spatial_spectral_efficiencies(:, 4));
end

%%
% spatial multiplexing

spatial_spectral_efficiencies_2 = zeros(size(point_cloud, 1), cols);
avg_spatial_spectral_efficiencies_2 = zeros(snr_res,1);

for idx_snr = 1:snr_res
    for i = 1:size(point_cloud, 1)
        H = path_losses{i, 4};
        H = H.^(-1);
        [~, sigma, ~] = svd(H);
        P = [sqrt(8), 0,0,0;0,sqrt(8),0,0;0,0,0,0;0,0,0,0];
    
        sum = 0;
        for sum_idx = 1:N_t/2
            hi = log2(1 + (snr(idx_snr) / N_t) * P(sum_idx, sum_idx) * sigma(sum_idx, sum_idx)^2);
            sum = sum + hi;
        end
    
        spatial_spectral_efficiencies_2(i, 1) = point_cloud(i,1);
        spatial_spectral_efficiencies_2(i,2) = point_cloud(i,2);
        spatial_spectral_efficiencies_2(i,3) = point_cloud(i,3);
        spatial_spectral_efficiencies_2(i, 4) = sum;
    end
    
    avg_spatial_spectral_efficiencies_2(idx_snr) = mean(spatial_spectral_efficiencies_2(:, 4));
end

%% visualize MRC+MRT results
% recalculate

min_snr_db = 20;
% from https://ntrs.nasa.gov/api/citations/19710001826/downloads/19710001826.pdf
min_snr = db2mag(min_snr_db);

for i = 1:size(point_cloud, 1)
    H = path_losses{i, 4};
    H = H.^(-1);
    [~, sigma, ~] = svd(H);
    P = [2, 0,0,0;0,2,0,0;0,0,2,0;0,0,0,2];

    sum = 0;
    for sum_idx = 1:N_t/2
        hi = log2(1 + (min_snr / N_t) * P(sum_idx, sum_idx) * sigma(sum_idx, sum_idx)^2);
        sum = sum + hi;
    end

    mrc_spectrial_efficiencies(i, 1) = point_cloud(i,1);
    mrc_spectrial_efficiencies(i,2) = point_cloud(i,2);
    mrc_spectrial_efficiencies(i,3) = point_cloud(i,3);
    mrc_spectrial_efficiencies(i, 4) = sum;
end

% Extracting the coordinates and path_loss values
x = mrc_spectrial_efficiencies(:,1);
y = mrc_spectrial_efficiencies(:,2);
z = mrc_spectrial_efficiencies(:,3);
mrc = mrc_spectrial_efficiencies(:,4);

% Plotting the scatter plot
figure(9);
scatter3(x(:), y(:), z(:), 50, mrc(:));
colorbar;
h = colorbar;
ylabel(h, 'Spectral Effiency (bits/second/Hz)');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Scatter Plot with Spectral Effiency for MRC+MRT');

%% visualize spatial multiplexing results for strongest channel only
% recalculate
for i = 1:size(point_cloud, 1)
    H = path_losses{i, 4};
    H = H.^(-1);
    [~, sigma, ~] = svd(H);
    P = [4, 0,0,0;0,0,0,0;0,0,0,0;0,0,0,0];

    sum = 0;
    for sum_idx = 1:N_t/2
        hi = log2(1 + (min_snr / N_t) * P(sum_idx, sum_idx) * sigma(sum_idx, sum_idx)^2);
        sum = sum + hi;
    end

    spatial_spectral_efficiencies(i, 1) = point_cloud(i,1);
    spatial_spectral_efficiencies(i,2) = point_cloud(i,2);
    spatial_spectral_efficiencies(i,3) = point_cloud(i,3);
    spatial_spectral_efficiencies(i, 4) = sum;
end

% Extracting the coordinates and path_loss values
x = spatial_spectral_efficiencies(:,1);
y = spatial_spectral_efficiencies(:,2);
z = spatial_spectral_efficiencies(:,3);
spec = spatial_spectral_efficiencies(:,4);

% Plotting the scatter plot
figure(10);
scatter3(x(:), y(:), z(:), 50, spec(:));
colorbar;
h = colorbar;
ylabel(h, 'Spectral Effiency (bits/second/Hz)');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Scatter Plot with Spectral Effiency for Spatial Multiplexing');

%% visualize spatial multiplexing results for strongest and second strongest channel equal
% recalculate
for i = 1:size(point_cloud, 1)
    H = path_losses{i, 4};
    H = H.^(-1);
    [~, sigma, ~] = svd(H);
    P = [sqrt(8), 0,0,0;0,sqrt(8),0,0;0,0,0,0;0,0,0,0];

    sum = 0;
    for sum_idx = 1:N_t/2
        hi = log2(1 + (min_snr / N_t) * P(sum_idx, sum_idx) * sigma(sum_idx, sum_idx)^2);
        sum = sum + hi;
    end

    spatial_spectral_efficiencies(i, 1) = point_cloud(i,1);
    spatial_spectral_efficiencies(i,2) = point_cloud(i,2);
    spatial_spectral_efficiencies(i,3) = point_cloud(i,3);
    spatial_spectral_efficiencies(i, 4) = sum;
end

% Extracting the coordinates and path_loss values
x = spatial_spectral_efficiencies_2(:,1);
y = spatial_spectral_efficiencies_2(:,2);
z = spatial_spectral_efficiencies_2(:,3);
spec_2 = spatial_spectral_efficiencies_2(:,4);

% Plotting the scatter plot
figure(14);
scatter3(x(:), y(:), z(:), 50, spec_2(:));
colorbar;
h = colorbar;
ylabel(h, 'Spectral Effiency (bits/second/Hz)');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Scatter Plot with Spectral Effiency for Spatial Multiplexing');

%% plot difference

figure(12);
hold on;
plot(snr_db, avg_mrc_spectral_effieciencies, 'DisplayName','MRC+MRT');
plot(snr_db, avg_spatial_spectral_efficiencies, 'DisplayName','Spatial Multiplexing');

hold off;
xlabel("SNR (dB)");
ylabel("Average Spectral Effiency Across all Points (SNR)");
title("Average Spectral Effiency for MRC+MRT vs. Spatial Multiplexing");
legend();

%% visualize how much better spatial multiplexing performs at each point

% Extracting the coordinates and path_loss values
x = spatial_spectral_efficiencies(:,1);
y = spatial_spectral_efficiencies(:,2);
z = spatial_spectral_efficiencies(:,3);
div = spatial_spectral_efficiencies(:,4) ./ mrc_spectrial_efficiencies(:,4);

% Plotting the scatter plot
figure(13);
scatter3(x(:), y(:), z(:), 50, div(:));
colorbar;
h = colorbar;
ylabel(h, 'sigma_1 / sigma_2');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Scatter Plot with Showing Spectral Efficiency for Spatial Multiplexing / Spectral Efficiency for MRC+MRT');

%% visualize how much better spatial multiplexing with one channel performs than 2 channels

% Extracting the coordinates and path_loss values
x = spatial_spectral_efficiencies(:,1);
y = spatial_spectral_efficiencies(:,2);
z = spatial_spectral_efficiencies(:,3);
div2 = spatial_spectral_efficiencies(:,4) ./ spatial_spectral_efficiencies_2(:,4);

% Plotting the scatter plot
figure(15);
scatter3(x(:), y(:), z(:), 50, div2(:));
colorbar;
h = colorbar;
ylabel(h, 'sigma_1 / sigma_2');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Scatter Plot with Showing Spectral Efficiency for Strongest Channel / Spectral Efficiency for Two Strongest Channels');