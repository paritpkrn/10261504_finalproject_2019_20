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
model = changeRxnBounds(model,'Biomass',1000,'u');
model = changeObjective(model,'Biomass');
FBAsolution_WT = optimizeCbModel(model,'max');
OptimalBiomass = FBAsolution_WT.f;
OTCprod_inWT = FBAsolution_WT.x(model.rxns=="OTCprod"); %หาค่า OTCprod ในสภาวะที่ biomass is optimal
glucose_inWT = FBAsolution_WT.x(model.rxns=="EX_glc(e)");
o2_inWT = FBAsolution_WT.x(model.rxns=="EX_o2(e)");
K = FBAsolution_WT.f;

%robustness analysis + plot the growth rate
%varied glucose
model = changeRxnBounds(model,'EX_o2(e)',o2_inWT,'b');

growthRates = zeros(101,1);

for i = 0:100
	model = changeRxnBounds(model,'EX_glc(e)',-i,'b');
	FBAsolution = optimizeCbModel(model,'max');
	growthRates(i+1) = FBAsolution.f;
end
plot(growthRates)
xlabel('Glucose uptake rate(mmol gDW-1 hr-1)')
ylabel('growth rate (hr-1)')

%varied o2
model = changeRxnBounds(model,'EX_glc(e)', glucose_inWT,'b');

growthRates = zeros(201,1);

for i = 0:200
	model = changeRxnBounds(model,'EX_glc(e)',-i,'b');
	FBAsolution = optimizeCbModel(model,'max');
	growthRates(i+1) = FBAsolution.f;
end
plot(growthRates)
xlabel('oxgen uptake rate(mmol gDW-1 hr-1)')
ylabel('growth rate (hr-1)')