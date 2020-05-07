initCobraToolbox (false) %initate cobra toolbox

solverName = 'gurobi'; %set solver

solverType = 'LP'; %set type of solver

changeCobraSolver(solverName, solverType);

fileName = 'Srimosus_Oct2019_validated_balanced.mat'; %load the srimosus model (validated and balanced)

if ~exist('modelOri','var')
    
modelOri = readCbModel(fileName);
end

model = modelOri;

% Build the rxnGeneMat based on the given models rules field
model = buildRxnGeneMat(model);

ReactionTable = [{'Reaction ID', 'Lower Bound', 'Upper Bound'};...
 
model.rxns, num2cell(model.lb), num2cell(model.ub)];
%cacaluate optimal biomass to find glucose value
model = changeRxnBounds(model,'Biomass',1000,'u'); %change biomass to unlimit

model = changeObjective(model,'Biomass');

FBAsolution = optimizeCbModel(model,'max');

Biomass = FBAsolution.f;

glucose_inWT = FBAsolution.x(model.rxns=="EX_glc(e)");
%set ip for pFBA
model = changeRxnBounds(model, 'EX_glc(e)', 0.5*glucose_inWT, 'l');

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

non_EG = [non_EG; n];

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


%it is necessary to know which genes are associated to the reactions not carrying any flux.
RxnGMat = full(model.rxnGeneMat);
pFBAfluxGenes = non_EG;
pFBAnoFluxGenes = [];
for i = 1:length(pFBAnoFluxRxn)
listGenes = find(RxnGMat(pFBAnoFluxRxn(i),:));
for n = 1:length(listGenes)
pos = find(non_EG==listGenes(n));
if pos
pFBAnoFluxGenes = [pFBAnoFluxGenes; model.genes(non_EG(pos))];
pFBAfluxGenes(pos) = 0;
end
end
end
pFBAnoFluxGenes = unique(pFBAnoFluxGenes);
pFBAfluxGenes(pFBAfluxGenes==0) = [];

%Identify MLE reactions:
% The list of reactions carying flux will be scanned, and the ones that 
% are "turned off" when the system is forced to achieve certain biomass production 
% are MLE reactions. MLE reations will be stored in the vector, |RxnMLE|, and 
% the remaining reactions will be stored in the vector, |restRxn|.
FBAsolution = optimizeCbModel(model, 'max'); 
model = changeRxnBounds(model, 'Biomass', FBAsolution.f, 'l');
%the percentage of optimal solution was set up 
% to 95%. This simulation provides a minimum and a maximum flux balance solution 
% that allows at least a 95% of the optimal solution for the objective function.
[minFlux2, maxFlux2] = fluxVariability(model,95);
RxnMLE = [];
restRxn = [];
for i = 1:length(pFBAfluxRxn)
    if (abs(minFlux2(pFBAfluxRxn(i)))<tol)&&(abs(maxFlux2(pFBAfluxRxn(i)))<tol);    
        RxnMLE = [RxnMLE pFBAfluxRxn(i)];
    else
        restRxn = [restRxn pFBAfluxRxn(i)];
    end
end

RxnMLEname = model.rxns(RxnMLE);

% Identify Optimal and ELE reactions:
FBAsolution = optimizeCbModel(model,'min','one');
model = changeRxnBounds(model, model.rxns, FBAsolution.x, 'b');
%run one last FVA for 100% of the optimal solution. The remaining 
% reactions in the |restRxn| variable were then clasified as Enzymatially Less 
% Eficient Reactions (RxnELE), if the reactions cannot carry any flux, or as Optimal 
% Reactions (RxnOptima), if they can carry flux.

[minFlux3, maxFlux3] = fluxVariability(model, 100);
pFBAopt_Rxns = model.rxns((abs(minFlux3)+abs(maxFlux3))>=tol);
pFBAopt_Rxns = unique(regexprep(pFBAopt_Rxns, '_[f|b]$',''));
pFBAopt_Rxns = setdiff(pFBAopt_Rxns, pFBAEssentialRxns);
ELE_Rxns = model.rxns((abs(minFlux3)+abs(maxFlux3))<=tol);
ELE_Rxns = setdiff(ELE_Rxns, RxnMLEname);
ELE_Rxns = setdiff(ELE_Rxns, ZeroFluxRxns);
RxnELE = findRxnIDs(model, ELE_Rxns);
RxnOptima = findRxnIDs(model, pFBAopt_Rxns);

%Classify the genes:
%genes that are related with each reaction. 
% The main point of this is to classify the genes into the 5 different
% groups  and store them in different vector

% # _*Essential genes:_* metabolic genes necessary for growth in the given media 
% ('essential_genes').
% # _*pFBA optima: _*non-essential genes contributing to the optimal growth 
% rate and minimum gene-associated flux ('OptimaGenes').
% # _*Enzymatically less efficient (ELE): _*genes requiring more flux through 
% enzymatic steps than alternative pathways that meet the same predicted growth 
% rate ('ELEGenes').
% # _*Metabolically less efficient (MLE):_* genes requiring a growth rate reduction 
% if used ('MLEGenes').
% # _*pFBA no-flux:_* genes that are unable to carry flux in the experimental 
% conditions ('pFBAnoFluxGenes').

OptimaGenes = [];
restGenes = pFBAfluxGenes;
for i = 1:length(RxnOptima)
    listGenes = find(RxnGMat(RxnOptima(i), :));
    for n = 1:length(listGenes)
        pos = find(pFBAfluxGenes==listGenes(n));
        if pos 
            OptimaGenes = [OptimaGenes; model.genes(pFBAfluxGenes(pos))];
            restGenes(pos,1) = 0;
        end 
    end
end
OptimaGenes = unique(OptimaGenes);
restGenes(restGenes==0) = [];    

ELEGenes = [];
restGenes2 = restGenes;
for i = 1:length(RxnELE)
    listGenes = find(RxnGMat(RxnELE(i), :));
    for n = 1:length(listGenes)
        pos = find(restGenes==listGenes(n));
        if pos 
            ELEGenes = [ELEGenes; model.genes(restGenes(pos))];
            restGenes2(pos, 1) = 0;
        end 
    end
end
ELEGenes = unique(ELEGenes);
restGenes2(restGenes2==0) = [];

MLEGenes = [];
finalRemainingGenes = restGenes2;
for i = 1:length(RxnMLE)
    listGenes = find(RxnGMat(RxnMLE(i),:));
    for n = 1:length(listGenes)
        pos = find(restGenes2==listGenes(n));
        if pos 
            MLEGenes = [MLEGenes; model.genes(restGenes2(pos))];
            finalRemainingGenes(pos, 1) = 0;
        end 
    end
end
MLEGenes = unique(MLEGenes);
finalRemainingGenes(finalRemainingGenes==0) = []; 

remainingGenes = [];
for n = 1:length(finalRemainingGenes)
        remainingGenes = [remainingGenes; model.genes(finalRemainingGenes(n))];
end
