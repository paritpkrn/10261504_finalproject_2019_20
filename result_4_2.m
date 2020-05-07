initCobraToolbox (false) %initate cobra toolbox

solverName = 'gurobi'; %set solver

solverType = 'LP'; %set type of solver

changeCobraSolver(solverName, solverType);

fileName = 'Srimosus_Oct2019_validated_balanced.mat'; %load the srimosus model (validated and balanced)

if ~exist('modelOri','var')
    
modelOri = readCbModel(fileName);
end

model = modelOri;

ReactionTable = [{'Reaction ID', 'Lower Bound', 'Upper Bound'};...
 
model.rxns, num2cell(model.lb), num2cell(model.ub)];

model=changeRxnBounds(model,'Biomass',1000,'u');
%Phenotypic phase plane analysis
%When performing robustness analysis, one parameter is varied and the network state is calculated.  
%It is also possible to vary two parameters simultaneously and plot the results as a phenotypic phase plane16.  
%These plots can reveal the interactions between two reactions in interesting ways
growthRates = zeros(201);
for i = 0:200
	for j = 0:200
		model = changeRxnBounds(model,'EX_glc(e)',-i,'b');
		model = changeRxnBounds(model,'EX_o2(e)',-j,'b');
		FBAsolution = optimizeCbModel(model,'max');
		growthRates(i+1,j+1) = FBAsolution.f;
	end
end
surfl (growthRates) %create 3D plot
%pcolor (growthRates) %create 2D plot
