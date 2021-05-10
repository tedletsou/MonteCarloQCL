%close all

WaveFuncData = importdata("WaveFunctions.txt");

Eigenenergies=WaveFuncData(1,:);
WaveFunctions=WaveFuncData(3:end,:);

GridData = importdata("Potential.txt");

Z = GridData(:, 2);
CBE = GridData(:, 1);

figure
plot(Z(1:(end)), CBE(1:(end)), 'k')
hold on
for i1 = 1:length(Eigenenergies)
    plot(Z(1:(end)), Eigenenergies(i1) + 0.15 * max(CBE) * (abs(WaveFunctions(:, i1)) .^ 2  ./ max(abs(WaveFunctions(:, i1)) .^ 2)))
end
Ang = char(197);
xlabel("Position (" + Ang +")")
ylabel("Energy (eV)")
% xlim([min(Z), max(Z)])

NumEigenStates=size(Eigenenergies,2)

%For Importing LOPhonon Scattering Rate Data

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 5);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = " \t";

% Specify column names and types
opts.VariableNames = ["VarName10", "VarName11", "Kix", "Ratex", "VarName15"];
opts.VariableTypes = ["double", "double", "double", "double", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "VarName15", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "VarName15", "EmptyFieldRule", "auto");

% Import the data
tbl = readtable("C:\Users\andre\source\repos\MonteCarloQCL\MonteCarloQCL\ScatteringRateLOEmission.txt", opts);

%% Convert to output type
VarName10 = tbl.VarName10;
VarName11 = tbl.VarName11;
Kix = tbl.Kix;
Ratex = tbl.Ratex;
VarName15 = tbl.VarName15;

size(Kix)

NumK=length(Kix)/(NumEigenStates^2)

for n=1:NumEigenStates
    for m=1:NumEigenStates
        for a=1:NumK
            RateMatSpline(n,m,a)=RateSpline((n-1)*NumK*NumEigenStates+(m-1)*NumK+a);
        end
    end
end
for n=1:NumEigenStates
    for m=1:NumEigenStates
        for a=1:NumK
            KiMatSpline(n,m,a)=KiSpline((n-1)*NumK*NumEigenStates+(m-1)*NumK+a);
        end
    end
end

figure
hold on
for n=1:6
    plot(reshape((KiMatSpline(n,1,:)),1,[]),reshape((RateMatSpline(n,1,:)),1,[]),ColorPick2(n))
end