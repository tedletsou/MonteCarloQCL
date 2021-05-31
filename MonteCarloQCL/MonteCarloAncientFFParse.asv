

clear opts

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = " \t";

% Specify column names and types
opts.VariableNames = ["q", "i", "f", "g", "h","FFOld"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the LO Phonon Emission RAte
tbl = readtable("C:\Users\andre\source\repos\MonteCarloQCL\MonteCarloQCL\formfactorsAncient.dat", opts);


NumEigenStates=12;
NumQ=50;

%% Phonon Rates 



FF=tbl.FFOld;
q=tbl.q;

for n=1:NumEigenStates
    for m=1:NumEigenStates
        for a=1:NumEigenStates
            for b=1:NumEigenStates
                for c=1:NumQ
                    FF(n,m,a,b,c)=FFvec((n-1)*NumQ*NumEigenStates^3+(m-1)*NumQ*NumEigenStates^2+(a-1)*NumQ*NumEigenStates+(b-1)*NumQ+c);
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
                    q(n,m,a,b,c)=qvec((n-1)*NumQ*NumEigenStates^3+(m-1)*NumQ*NumEigenStates^2+(a-1)*NumQ*NumEigenStates+(b-1)*NumQ+c);
                end
            end
        end
    end
end

figure
hold on
for n=1:6
    plot(reshape((q(n,1,:)),1,[]),reshape((FF(n,1,:)),1,[]),ColorPick2(n))
end