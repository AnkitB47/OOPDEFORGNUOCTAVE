% Define the output file name
outputFile = 'channel_geometry.txt';

% Run the function to generate the channel geometry
[vertices, edges] = generateChannelGeometry(outputFile);

outputFileName = 'channelgeometery.txt';

generatePolyhedralGeometry(vertices, outputFileName, 2);


