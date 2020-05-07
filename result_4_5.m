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

model = changeRxnBounds(model,'Biomass',1000,'u'); %change biomass to unlimit

%Identify essentail reactions: perform a gene knocked-out analysis.
[grRatio, grRateKO, grRateWT, delRxns,...
hasEffect] = singleGeneDeletion(model, 'FBA', model.genes);

essential_genes = [];

essential_genes_location =[];

non_EG = [];

tol = 1e-6;

for n = 1:length(grRateKO)

if (grRateKO(n)<tol)||(isnan(grRateKO(n)) == 1)

essential_genes = [essential_genes; model.genes(n)];

essential_genes_location = [essential_genes_location; n];
else

non_EG = [non_EG; model.genes(n)];

end

end

%signle find essential reactions
RxnRatio = singleRxnDeletion(model);

RxnRatio(isnan(RxnRatio)) = 0;

pFBAEssentialRxns = model.rxns(RxnRatio < tol);

%Identify non-essentail reactions that can or cannot carry flux:
[minFluxglc, maxFluxglc] = fluxVariability(model, 0);

pFBAnoFluxRxn = [];

pFBAfluxRxn = [];

for i=1:length(model.rxns)

if (abs(minFluxglc(i))<tol)&&(abs(maxFluxglc(i))<tol)

pFBAnoFluxRxn = [pFBAnoFluxRxn i];

else

pFBAfluxRxn = [pFBAfluxRxn i];

end

end

pFBAfluxRxn = pFBAfluxRxn';

ZeroFluxRxns = model.rxns(pFBAnoFluxRxn);