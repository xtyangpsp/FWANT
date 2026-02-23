% MATLAB script using Mapping Toolbox
filename = 'station_riverrise_withdata_formatok.txt';

% 1. Read the data
% readtable handles the whitespace-separated columns automatically
opts = delimitedTextImportOptions("NumVariables", 5);
opts.VariableNames = ["Network", "Station", "Longitude", "Latitude", "Elevation"];
opts.VariableTypes = ["string", "string", "double", "double", "double"];
opts.Delimiter = " ";
T = readtable(filename, opts);

% 2. Define the Reference Ellipsoid (WGS84)
% This provides more accuracy than a simple sphere
wgs84 = referenceEllipsoid('wgs84', 'm');

% 3. Compute Pairwise Distances
numStations = size(T, 1);
numPairs = nchoosek(numStations, 2);
distances = zeros(numPairs, 1);

lat = T.Latitude;
lon = T.Longitude;

pairCount = 1;
for i = 1:numStations
    for j = i+1:numStations
        % Built-in distance function from Mapping Toolbox
        % Returns distance in meters when using the wgs84 object
        distances(pairCount) = distance(lat(i), lon(i), lat(j), lon(j), wgs84);
        pairCount = pairCount + 1;
    end
end

% 4. Visualization
figure('Color', 'w', 'Position', [100, 100, 1100, 450]);

% Subplot 1: Histogram
subplot(1, 2, 1);
h = histogram(distances, 30);
h.FaceColor = [0.2 0.6 0.8];
h.EdgeColor = 'w';
title('Distribution of All-Pair Distances');
xlabel('Great Circle Distance (meters)');
ylabel('Frequency');
grid on;

% Subplot 2: Scatter Plot (Geographic Distribution)
subplot(1, 2, 2);
scatter(T.Longitude, T.Latitude, 60, T.Elevation, 'filled', 'MarkerEdgeColor', 'k');
colormap(parula); % or 'viridis' if available in your version
cb = colorbar;
cb.Label.String = 'Elevation (m)';
title('Station Geographic Layout');
xlabel('Longitude');
ylabel('Latitude');
grid on;

% Print Summary
fprintf('Processed %d stations.\n', numStations);
fprintf('Mean separation: %.2f meters.\n', mean(distances));