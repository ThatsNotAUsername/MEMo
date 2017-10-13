function [ rxns_index] = get_index_of_rxns(network,varargin)
%get indexes of the reactions in varargin. If varargin is empty the default
%reaction is biomass reaction. I several reactions have the same the first
%one is chosen.
%INPUT : network structure with the field rxnNames
%        varargin names of reactions
rxns_index=[];
if nargin==1
    
cells_biomass=strfind(upper(network.rxnNames),'BIOMASS');

rxns_index=find(~cellfun(@isempty,cells_biomass));

l=length(rxns_index);
    if l==0      
        error(' There is not any biomass reaction')

    elseif  l>1
        error('There is more than one biomass reaction')
    end
    rxns_index=rxns_index(1) ;
    
    else for i=1:(nargin-1)
        rxn_name=varargin{i};

        cells_obj_rxns=strfind(upper(network.rxnNames),upper(rxn_name)); 

        new_index=find(~cellfun(@isempty,cells_obj_rxns));
        l= length(new_index);
        if l==0
            error(' There is not any %s reaction. \n',upper(rxn_name))
          
        elseif length(new_index)>1
            error('There is more than one %s reaction. \n',rxn_name)
        end
        rxns_index=[rxns_index new_index(1)];
    end
end

end

