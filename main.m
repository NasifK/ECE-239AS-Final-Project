%% visualize Faustini Crater --- make sure everything is importing properly
figure(1);
model = createpde;
importGeometry(model,"moon_10.stl");
pdegplot(model);
title("Faustini Crater Visualization");

%% create digital elevation model
% point cloud
stlData = stlread('moon_25.stl');
point_cloud = stlData.Points;
figure(2);
scatter3(point_cloud(:,1), point_cloud(:,2), point_cloud(:,3), 5, 'filled');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Point Cloud Visualization');

% create digital elevation model
ptCloud = pointCloud(point_cloud);
digital_elevation_model = pc2dem(ptCloud);
figure(3);
pcshow(ptCloud.Location)
title("Point Cloud")

figure(4);
imagesc(digital_elevation_model)
colormap(gray)
title("Digital Elevation Model")

%% create cubic DEM
% Define cube size
cubeSize = 2; % Adjust as needed

% Define cube color
cubeColor = [0.5, 0.5, 0.5]; % Gray color

% Extract x, y, z coordinates from the point cloud
x = ptCloud.Location(:, 1);
y = ptCloud.Location(:, 2);
z = ptCloud.Location(:, 3);

% Create a new figure for the DEM visualization using cubes
figure;
hold on;
axis equal;

% Loop through each point in the point cloud and create a cube at that point
for i = 1:size(point_cloud, 1)
    cube = [x(i)-cubeSize/2, y(i)-cubeSize/2, z(i)-cubeSize/2, cubeSize, cubeSize, cubeSize];
    drawCube(cube, cubeColor);
end

xlabel('X');
ylabel('Y');
zlabel('Z');
title('Digital Elevation Model (Cubes)');
view([45, 45]);

%% begin raytracing
fc = 2.4e9; % transmitter / recieving recieving at 2.9 GHz
lambda = physconst("lightspeed")/fc;
txArray = arrayConfig("Size",[4 1],"ElementSpacing",2*lambda);
rxArray = arrayConfig("Size",[4 1],"ElementSpacing",lambda);

% define tx and rx for our STL
tx = txsite("cartesian", ...
    "Antenna",txArray, ...
    "AntennaPosition",[210; 220; 22], ...
    'TransmitterFrequency',2.4e9);

rx = rxsite("cartesian", ...
    "Antenna",rxArray, ...
    "AntennaPosition",[100; 200; 30]);

% view transmitter reciever pair
siteviewer("SceneModel","moon_25.stl");
show(tx,"ShowAntennaHeight",true)
show(rx,"ShowAntennaHeight",true)

% define propogation model
pm = propagationModel("raytracing", ...
    "CoordinateSystem","cartesian", ...
    "Method","Image", ...
    "AngularSeparation","low", ...
    "MaxNumReflections",2, ...
    "SurfaceMaterial","custom");

rays = raytrace(tx,rx, pm);
path_loss = rays{1}.PathLoss;
total_path_loss = sum(path_loss);
disp(total_path_loss)

% show point cloud
ptCloud = pointCloud(point_cloud);
figure(3);
pcshow(ptCloud.Location)
title("Point Cloud")

%% perform large-scale tx-rx analysis
% grab point cloud and STL data

% point cloud
stlData = stlread('moon_10.stl');
point_cloud = stlData.Points;
figure(2);
scatter3(point_cloud(:,1), point_cloud(:,2), point_cloud(:,3), 5, 'filled');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Point Cloud Visualization');

% create digital elevation model
ptCloud = pointCloud(point_cloud);
digital_elevation_model = pc2dem(ptCloud);
figure(3);
pcshow(ptCloud.Location)
title("Point Cloud")

% define ONE tx site
tx = txsite("cartesian", ...
    "Antenna",txArray, ...
    "AntennaPosition",[210; 220; 20], ...
    'TransmitterFrequency',2.4e9);


rx = rxsite("cartesian", ...
    "Antenna",rxArray, ...
    "AntennaPosition",[100; 200; 30]);

%test
rays = raytrace(tx,rx, pm);
path_loss = rays{1}.PathLoss;
total_path_loss = sum(path_loss);
disp(total_path_loss)


%%
% point cloud
stlData = stlread('moon_10.stl');
point_cloud = stlData.Points;

% define propogation model
pm = propagationModel("raytracing", ...
    "CoordinateSystem","cartesian", ...
    "Method","Image", ...
    "AngularSeparation","low", ...
    "MaxNumReflections",2, ...
    "SurfaceMaterial","custom", ...
    "SurfaceMaterialPermittivity",.297, ...
    "SurfaceMaterialConductivity",0.00848);

% permittivity:
% https://earth-planets-space.springeropen.com/articles/10.1186/s40623-020-01259-2,
% ε1 = 0.297 S + 0.0107

% conductivity
% https://academic.oup.com/nsr/article/9/11/nwac175/6677425

cols = size(point_cloud, 2) + 1; % add one more column for path loss
path_losses = zeros(size(point_cloud, 1), cols);
disp("Starting raytracing...")

% create digital elevation model
ptCloud = pointCloud(point_cloud);
digital_elevation_model = pc2dem(ptCloud);
figure(3);
pcshow(ptCloud.Location)
title("Point Cloud")

%% help
% define tx and rx for our STL
tx = txsite("cartesian", ...
    "Antenna",txArray, ...
    "AntennaPosition",[190; 320; 40], ...
    'TransmitterFrequency',2.4e9);

%rx = rxsite("cartesian", ...
%    "Antenna",rxArray, ...
%    "AntennaPosition",[100; 200; 30]);


i = 1;
rx = rxsite("cartesian", ...
"Antenna",rxArray, ...
"AntennaPosition",[point_cloud(i,1); point_cloud(i,2); point_cloud(i,3)+2]);

siteviewer("SceneModel","moon_10.stl");
show(tx,"ShowAntennaHeight",true)
show(rx,"ShowAntennaHeight",true)

rays = raytrace(tx,rx, pm);
path_loss = rays{1}.PathLoss;
if isvector(path_loss) && numel(path_loss) > 1
    total_path_loss = sum(path_loss);
else
    total_path_loss = path_loss;
end
%total_path_loss = sum(path_loss);
disp(total_path_loss);
%% perform raytracing for each point in the point_cloud

for i = 1:size(point_cloud, 1)
    rx = rxsite("cartesian", ...
    "Antenna",rxArray, ...
    "AntennaPosition",[point_cloud(i,1); point_cloud(i,2); point_cloud(i,3)+2]);

    rays = raytrace(tx,rx, pm);
    path_loss = rays{1}.PathLoss;
    if isvector(path_loss) && numel(path_loss) > 1
        total_path_loss = sum(path_loss);
    else
        total_path_loss = path_loss;
    end
    disp(total_path_loss);

    % add point to the path losses
    path_losses(i, 1) = point_cloud(i,1);
    path_losses(i,2) = point_cloud(i,2);
    path_losses(i,3) = point_cloud(i,3);
    path_losses(i,4) = total_path_loss;
end

%% visualize raytracing results

% Extracting the coordinates and path_loss values
x = path_losses(:,1);
y = path_losses(:,2);
z = path_losses(:,3);
path_loss = path_losses(:,4);

% Plotting the scatter plot
figure(13);
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
N = 4;
wavelength = 0.125; % for 2.4 GHz wave
path_losses = cell(size(point_cloud, 1), cols);
pm.MaxRelativePathLoss = Inf;

for i = 1:size(point_cloud, 1)
    H = zeros(1,N);
    for antenna_id = 1:N
        rx = rxsite("cartesian", ...
        "Antenna",rxArray, ...
        "AntennaPosition",[point_cloud(i,1); point_cloud(i,2) + (0.5*wavelength*antenna_id - 1.5); point_cloud(i,3)+2]);
    
        rays = raytrace(tx,rx, pm);
        path_loss = rays{1}.PathLoss;
        path_loss = db2mag(path_loss);
            if isvector(path_loss) && numel(path_loss) > 1
                total_path_loss = sum(path_loss) ./ 1000;
            else
                total_path_loss = path_loss / 1000;
            end
        H(antenna_id) = total_path_loss;
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
min_snr_db = 14;
% from https://ntrs.nasa.gov/api/citations/19710001826/downloads/19710001826.pdf
min_snr = mag2db(min_snr_db);

avg_mrc_spectral_effieciencies = zeros(snr_res, 1);

for idx_snr = 1:snr_res
    
    for i = 1:size(point_cloud, 1)
        H = path_losses{i, 4};
        H = H.^(-1);
        H = reshape(H, 2,2);
        [~, sigma, ~] = svd(H);
        P = [1, 0;0,1];
        
        sum = 0;
        for sum_idx = 1:N/2
            hi = log2(1 + (snr(idx_snr) / N) * P(sum_idx, sum_idx) * sigma(sum_idx, sum_idx)^2);
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
        H = reshape(H, 2,2);
        [~, sigma, ~] = svd(H);
        P = [2, 0;0,0];
    
        sum = 0;
        for sum_idx = 1:N/2
            hi = log2(1 + (snr(idx_snr) / N) * P(sum_idx, sum_idx) * sigma(sum_idx, sum_idx)^2);
            sum = sum + hi;
        end
    
        spatial_spectral_efficiencies(i, 1) = point_cloud(i,1);
        spatial_spectral_efficiencies(i,2) = point_cloud(i,2);
        spatial_spectral_efficiencies(i,3) = point_cloud(i,3);
        spatial_spectral_efficiencies(i, 4) = sum;
    end
    
    avg_spatial_spectral_efficiencies(idx_snr) = mean(spatial_spectral_efficiencies(:, 4));
end
%% visualize MRC+MRT results

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

%% visualize spatial multiplexing results

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