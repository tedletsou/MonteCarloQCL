% Import WaveFunctions.txt
close all

WaveFuncData = importdata("WaveFunctions.txt");

Eigenenergies=WaveFuncData(1,:);
WaveFunctions=WaveFuncData(3:end,:);

% Import Potential.txt
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

%% Phonon Rates 


% Import the LO Phonon Emission RAte
tbl = readtable("C:\Users\andre\source\repos\MonteCarloQCL\MonteCarloQCL\ScatteringRateLOEmission.txt", opts);

% Import the LO Phonon Absoprtion RAte
tb2 = readtable("C:\Users\andre\source\repos\MonteCarloQCL\MonteCarloQCL\ScatteringRateLOAbsorb.txt", opts);


%% EE Form Factor

clear opts

% Change opts to import EE Form Factor
opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = " \t";

% Specify column names and types
opts.VariableNames = ["i", "f", "g", "H", "qx", "FFx"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";


% Import the EE Form Factor Original UpSampled
tb3 = readtable("C:\Users\andre\source\repos\MonteCarloQCL\MonteCarloQCL\FormFactorEEOriginal.txt", opts);


%% EE Form Factor Upsampled

clear opts

% Change opts to import EE Form Factor
opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = " \t";

% Specify column names and types
opts.VariableNames = ["i", "f", "g", "H", "qx", "FFx"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the EE Form Factor Original UpSampled
tb4 = readtable("C:\Users\andre\source\repos\MonteCarloQCL\MonteCarloQCL\FormFactorEEUpSampled.txt", opts);


%% EE Scattering Rate

clear opts

% Change opts to import EE Form Factor
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = " \t";

% Specify column names and types
opts.VariableNames = ["i", "f", "k", "EEx"];
opts.VariableTypes = ["double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the EE Scattering Rate
tb5 = readtable("C:\Users\andre\source\repos\MonteCarloQCL\MonteCarloQCL\ScatteringRateEE.txt", opts);


%% Convert to output type
VarName10 = tbl.VarName10;
VarName11 = tbl.VarName11;
KiE = tbl.Kix;
RateE = tbl.Ratex;
VarName15 = tbl.VarName15;

KiA = tb2.Kix;
RateA = tb2.Ratex;

NumK=length(KiE)/(NumEigenStates^2);



Qvec = tb3.qx;
FFEEvec = tb3.FFx;

NumQ=length(Qvec)/(NumEigenStates^4)

QvecUp = tb4.qx;
FFEEUpvec = tb4.FFx;

NumQUp=length(QvecUp)/(NumEigenStates^4)

kEE = tb5.k;
EERate = tb5.EEx;

for n=1:NumEigenStates
    for m=1:NumEigenStates
        for a=1:NumK
            RateMatEmit(n,m,a)=RateE((n-1)*NumK*NumEigenStates+(m-1)*NumK+a);
        end
    end
end
for n=1:NumEigenStates
    for m=1:NumEigenStates
        for a=1:NumK
            KiMatEmit(n,m,a)=KiE((n-1)*NumK*NumEigenStates+(m-1)*NumK+a);
        end
    end
end

figure
hold on
for n=1:6
    plot(reshape((KiMatEmit(n,1,:)),1,[]),reshape((RateMatEmit(n,1,:)),1,[]),ColorPick2(n))
end

for n=1:NumEigenStates
    for m=1:NumEigenStates
        for a=1:NumK
            RateMatAbs(n,m,a)=RateA((n-1)*NumK*NumEigenStates+(m-1)*NumK+a);
        end
    end
end
for n=1:NumEigenStates
    for m=1:NumEigenStates
        for a=1:NumK
            KiMatAbs(n,m,a)=KiA((n-1)*NumK*NumEigenStates+(m-1)*NumK+a);
        end
    end
end

figure
hold on
for n=1:6
    plot(reshape((KiMatAbs(n,1,:)),1,[]),reshape((RateMatAbs(n,1,:)),1,[]),ColorPick2(n))
end

for n=1:NumEigenStates
    for m=1:NumEigenStates
        for a=1:NumEigenStates
            for b=1:NumEigenStates
                for c=1:NumQ
                    FFEE(n,m,a,b,c)=FFEEvec((n-1)*NumQ*NumEigenStates^3+(m-1)*NumQ*NumEigenStates^2+(a-1)*NumQ*NumEigenStates+(b-1)*NumQ+c);
                end
            end
        end
    end
end

for n=1:NumEigenStates
    for m=1:NumEigenStates
        for a=1:NumEigenStates
            for b=1:NumEigenStates
                for c=1:NumQ
                    qEE(n,m,a,b,c)=Qvec((n-1)*NumQ*NumEigenStates^3+(m-1)*NumQ*NumEigenStates^2+(a-1)*NumQ*NumEigenStates+(b-1)*NumQ+c);
                end
            end
        end
    end
end

figure
hold on
for n=1:1
    plot(reshape((qEE(n,1,1,2,:)),1,[]),(reshape((FFEE(n,1,1,2,:)),1,[])),strcat(ColorPick2(n),'o'))
end
%axis([1e9 3e9 0 1e-20])


for n=1:NumEigenStates
    for m=1:NumEigenStates
        for a=1:NumEigenStates
            for b=1:NumEigenStates
                for c=1:NumQUp
                    FFEEUp(n,m,a,b,c)=FFEEUpvec((n-1)*NumQUp*NumEigenStates^3+(m-1)*NumQUp*NumEigenStates^2+(a-1)*NumQUp*NumEigenStates+(b-1)*NumQUp+c);
                end
            end
        end
    end
end

for n=1:NumEigenStates
    for m=1:NumEigenStates
        for a=1:NumEigenStates
            for b=1:NumEigenStates
                for c=1:NumQUp
                    qEEUp(n,m,a,b,c)=QvecUp((n-1)*NumQUp*NumEigenStates^3+(m-1)*NumQUp*NumEigenStates^2+(a-1)*NumQUp*NumEigenStates+(b-1)*NumQUp+c);
                end
            end
        end
    end
end

figure
hold on
for n=1:1
    plot(reshape((qEEUp(n,1,1,1,:)),1,[]),(reshape(FFEEUp(n,1,1,1,:),1,[])),strcat(ColorPick2(n),'o'))
end


for n=1:NumEigenStates
    for m=1:NumEigenStates
        for a=1:NumK
            RateMatEE(n,m,a)=EERate((n-1)*NumK*NumEigenStates+(m-1)*NumK+a);
        end
    end
end

for n=1:NumEigenStates
    for m=1:NumEigenStates
        for a=1:NumK
            KiMatEE(n,m,a)=kEE((n-1)*NumK*NumEigenStates+(m-1)*NumK+a);
        end
    end
end


figure
hold on
for n=1:6
    plot(reshape((KiMatEE(n,1,:)),1,[]),reshape((RateMatEE(n,1,:)),1,[]).^2,ColorPick2(n))
end

%axis([1e9 3e9 0 1e-20])