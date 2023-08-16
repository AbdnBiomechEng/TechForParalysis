function data = importfile_mot(filename)

dataLines = [6, Inf];

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 15);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "  ";

% Specify column names and types
opts.VariableNames = ["time", "TH_x", "TH_y", "TH_z", "SC_y", "SC_z", "SC_x", "AC_y", "AC_z", "AC_x", "GH_y", "GH_z", "GH_yy", "EL_x", "PS_y"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
data = readtable(filename, opts);
data{:,2:15} = data{:,2:15}*pi/180; % convert to radians

end