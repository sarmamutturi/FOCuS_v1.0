function [index]  = indexing(model,biom,rxn_max,bound_glu,bound_O2,bound_ATPM)

%%%%%%%%%%%%% inputs to indexing program %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% model             COBRA stoichiometric model in standard structure form
% biom              Biomass reaction
% rxn_max           Target reaction
% bound_glu         input flux value for glucose
% bound_O2          input flux value for oxygen
% bound_ATPM        input flux value for maintenance ATP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% Input Parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RXN Name %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model_rxns=model.rxns;

rxns_size = size(model.rxns);
rxns_count = rxns_size(1);

for i=1:1:rxns_count
    
    rxn_ind = strcmp(model.rxns(i),rxn_max);
    if rxn_ind > 0
        rxn_max_no = i;
    end
    
end
clear i

disp('Target reaction number =')
disp(rxn_max_no);

%biom='Ec_biomass_iAF1260_core_59p81M'

for i=1:1:rxns_count
    
    biomass_ind = strcmp(model.rxns(i),biom);
    if biomass_ind > 0
        biomass_max_no = i;
    end
    
end

disp('Biomass reaction number =')
disp(biomass_max_no);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Lethal genes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

biomass_lethal=ones(length(model_rxns),1);
ind = 1:length(model.rxns);

model = changeRxnBounds(model,'EX_glc(e)',bound_glu,'l');           %glucose
model = changeRxnBounds(model,'EX_o2(e)',bound_O2,'l');        %oxygen uptake
model = changeRxnBounds(model,'ATPM',bound_ATPM,'l');



for i=1:1:length(model_rxns)
    
    model_rxns{i};
    model_del=changeRxnBounds(model,model_rxns{i},0,'b');
    FBAsolution = optimizeCbModel(model_del,'max');
    if FBAsolution.f==0
        biomass_lethal(i)=0;
    else
        biomass_lethal(i) = FBAsolution.x(biomass_max_no);
    end
    clear model_del FBAsolution
    
end

p_3=find(not(biomass_lethal));
sizp3=size(p_3);

for i=1:sizp3
    rxn = p_3(i);
    ind(rxn) = 0;
end

model.rxns(p_3)
clear rxn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Not to delete following rxns %%%%%%%%%%%%%%%%%%

SpecRxnsRemove = {biom,'ATPM', rxn_max};
length(SpecRxnsRemove);

for i=1:length(SpecRxnsRemove)
    rxn = SpecRxnsRemove{i};
    Rxn_ID = find(strcmp(rxn,model.rxns));
    Rxn_IN = ind==Rxn_ID;
    ind(Rxn_IN) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear Rxn_ID Rxn_IN

k_1=zeros(length(model.grRules),1);
k_2=zeros(1,length(model.grRules));

for i=1:length(model.grRules)
    k_1(i) = strcmp(model.grRules(i),'');
    if k_1(i) == 1
        k_2(i) = 1;
    else
        k_2(i) = 0;
    end
end

Nogenes_ind = find(k_2);

clear k_1 k_2

for i = 1:length(Nogenes_ind)
    Rxn_ID = Nogenes_ind(i);
    Rxn_IN = ind==Rxn_ID;
    ind(Rxn_IN) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind=find(ind);
index=ind';

xlswrite('C:\Users\....\MATLAB\cobra\FOCUS\acetate_index.xlsx',index)
end



