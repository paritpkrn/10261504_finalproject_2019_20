initCobraToolbox (false) %initate cobra toolbox
solverName = 'gurobi'; %set solver
solverType = 'LP'; %set type of solver
changeCobraSolver(solverName, solverType);
fileName = 'Srimosus_Oct2019_validated_balanced.mat';
if ~exist('modelOri','var')
modelOri = readCbModel(fileName);
end
model = modelOri;
model=changeRxnBounds(model,'Biomass',1000,'u');
model = changeObjective(model,'Biomass');
FBAsolution = optimizeCbModel(model,'max');
biomass=zeros(10,1);
OTC= zeros(10,1);
model = changeObjective(model,'OTCprod');
for i=1:length(OTC)
model= changeRxnBounds(model,'Biomass',i*0.1*FBAsolution.f,'l');
solution=optimizeCbModel(model);
OTC(i)=solution.f;
biomass(i)=solution.x(model.rxns=="Biomass");
end
SlopeBetweenPoints = diff(OTC)./diff(biomass);
SlopeBetweenPoints;
plot(biomass,OTC,'-o')