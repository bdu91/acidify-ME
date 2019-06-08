import pandas as pd
import numpy as np
import cobrame
from cobrame.util import mu
from cobrame.core.component import Metabolite
from math import pi, log

#define constants
phospholipid_mol_fraction = {'pe':0.765, 'pg':0.184, 'clpn':0.051} #the composition between three different kinds of phospholipids, add up to 1
phospholipid_cyto_fraction = {'pe': 0.28, 'pg': 0.538, 'clpn':0.0}#percent of each phospholipid in the cytoplasm (used for inner membrane), rest used for outer membrane
phospholipid_SA_dict = {'pe':50.0,'pg':50.0,'clpn':100.0,'cpe':50.0,'cpg':50.0,'cclpn':100.0,'kdo':121.2}
lipid_overall_wt_percent = 0.093 #0.093 is the weight fraction of total lipid over biomass

#load the fatty acid composition under different conditions, control means that it is grown under pH 7, habituated meaning adapted under pH 5, dehabituated means first under pH 5 then pH 7
fatty_acid_comp_df = pd.read_csv('data/lipid_fatty_acid_composition_under_acid.csv', index_col=0)
fatty_acid_comp_df.index = fatty_acid_comp_df.index.map(unicode)

def get_string_dict(rxn_dict):
    #only for cobra model reaction dict, convert the reaction dict with keys as metabolite object to metabolite name
    rxn_string_dict = {}
    for met, stoich in rxn_dict.iteritems():
        rxn_string_dict[met.id] = stoich
    return rxn_string_dict

def get_lipid_type(lipid_met_id):
    #typically lipid metabolite id are in format letters followed by numbers, so we extract all letters on the left side of the first number occurred
    lipid_type = ''; index = 0
    while (lipid_met_id[index] >= 'a' and lipid_met_id[index] <= 'z'):
        lipid_type += lipid_met_id[index]
        index += 1
    return lipid_type

def get_lipid_composition(me_model, fatty_acid_comp_df, cur_strain_name):
    #get the coefficient for different types of lipid for the ME model lipid demand reactions
    #adjust the fatty acid fraction based on differen strains to obtain the coefficient for different lipids on the membrane
    cur_strain_fa_comp_df = fatty_acid_comp_df[cur_strain_name].loc[(fatty_acid_comp_df[cur_strain_name].notnull()) & ~(fatty_acid_comp_df.index.str.contains('total'))]
    fatty_acid_fraction = (cur_strain_fa_comp_df/cur_strain_fa_comp_df.sum()).to_dict()

    lipid_composition_df = pd.DataFrame(columns=['met','mol_fraction'])
    cur_index = 0
    for cur_lipid, cur_lipid_mol_fraction in phospholipid_mol_fraction.iteritems():
        for cur_fa, cur_fa_mol_fraction in fatty_acid_fraction.iteritems():
            cyto_frac = phospholipid_cyto_fraction[cur_lipid]; peri_frac = 1 - cyto_frac
            if len(cur_fa.split('_')) == 1: #not cyclopropane fatty acid
                cur_phospholipid_name = cur_lipid + cur_fa.split('_')[0]
            elif 'clpn' in cur_lipid: #clpn not have cyclopropane fatty acid attached for now
                continue
            else: 
                cur_phospholipid_name = cur_fa.split('_')[0] + cur_lipid + cur_fa.split('_')[1] #to consider the name of cyclopropane fatty acid
                cyto_frac = 1.0; peri_frac = 0.0 #cyclopropane fatty acid is only present in inner membrane
            if cyto_frac > 0:
                lipid_composition_df.loc[cur_index] = {'met': cur_phospholipid_name + '_c', 'mol_fraction': cur_lipid_mol_fraction*cyto_frac*cur_fa_mol_fraction}
                cur_index += 1
            if peri_frac > 0:
                lipid_composition_df.loc[cur_index] = {'met': cur_phospholipid_name + '_p', 'mol_fraction': cur_lipid_mol_fraction*peri_frac*cur_fa_mol_fraction}
                cur_index += 1
    
    #now calculate the demand coefficient for each lipid molecule in terms of mmol_gDW step by step
    lipid_composition_df['MW'] = [me_model.metabolites.get_by_id(cur_name).formula_weight for cur_name in lipid_composition_df.met.tolist()]
    lipid_composition_df['mg_if_total_1mmol'] = lipid_composition_df['mol_fraction'] * lipid_composition_df['MW']
    lipid_composition_df['lipid_weight_fraction'] = lipid_composition_df['mg_if_total_1mmol']/lipid_composition_df['mg_if_total_1mmol'].sum()
    lipid_composition_df['biomass_weight_fraction'] = lipid_composition_df['lipid_weight_fraction'] * lipid_overall_wt_percent
    lipid_composition_df['mmol_gDW'] = lipid_composition_df['biomass_weight_fraction'] / lipid_composition_df['MW'] * 1000
    return lipid_composition_df

def add_lipid_demand(me_model, lipid_composition_df):
    #add/modify the demand reaction for all the lipids
    for i, row in lipid_composition_df.iterrows():
        component_mass = me_model.metabolites.get_by_id(row['met']).formula_weight / 1000.
        cur_rid = 'Demand_' + row['met']
        if me_model.reactions.has_id(cur_rid): #if the reaction exist then we will modify the stoichiometry
            me_met_obj = me_model.metabolites.get_by_id(row['met'])
            me_lipid_biomass_obj = me_model.metabolites.get_by_id('lipid_biomass')
            me_model.reactions.get_by_id(cur_rid)._metabolites[me_met_obj] = -1 * abs(row['mmol_gDW'])
            me_model.reactions.get_by_id(cur_rid)._metabolites[me_lipid_biomass_obj] = component_mass * abs(row['mmol_gDW'])
        else: #if the reaction does not exist we will add the reaction to the model
            rxn = cobrame.SummaryVariable(cur_rid)
            rxn.add_metabolites({me_model.metabolites.get_by_id(row['met']): -1 * abs(row['mmol_gDW']),
                                 me_model.metabolites.get_by_id('lipid_biomass'): component_mass * abs(row['mmol_gDW'])})
            rxn.lower_bound = mu; rxn.upper_bound = 1000.
            me_model.add_reaction(rxn)

def add_lipid_fraction_constraint(me_model, lipid_composition_df, constraint_type = 'mass'):
    #add the lipid fraction constraint on surface area between different kind of lipids
    clpn_pg_constraint = Metabolite('clpn_pg_constraint'); pe_pg_constraint = Metabolite('pe_pg_constraint')
    for i, row in lipid_composition_df.iterrows():
        me_met_obj = me_model.metabolites.get_by_id(row['met']); met_demand_rxn = 'Demand_' + me_met_obj.id
        if len(set(get_string_dict(me_model.reactions.get_by_id(met_demand_rxn).metabolites).keys()).intersection(set(['clpn_pg_constraint','pe_pg_constraint']))) > 0: continue
        cur_met_coeff = abs(me_model.reactions.get_by_id(met_demand_rxn).get_coefficient(row['met']))
        if 'pe' == me_met_obj.id[:2] or 'cpe' == me_met_obj.id[:3]: #add pe_pg_constraint
            cur_pe_pg_coeff = phospholipid_mol_fraction['pg'] if constraint_type != 'mass' else phospholipid_mol_fraction['pg'] * me_met_obj.formula_weight
            me_model.reactions.get_by_id(met_demand_rxn).add_metabolites({pe_pg_constraint: cur_pe_pg_coeff*cur_met_coeff})
        elif 'pg' == me_met_obj.id[:2] or 'cpg' == me_met_obj.id[:3]: #add pe_pg_constraint and clpn_pg_constraint
            cur_pe_pg_coeff = phospholipid_mol_fraction['pe'] if constraint_type != 'mass' else phospholipid_mol_fraction['pe'] * me_met_obj.formula_weight
            cur_clpn_pg_coeff = phospholipid_mol_fraction['clpn'] if constraint_type != 'mass' else phospholipid_mol_fraction['clpn'] * me_met_obj.formula_weight
            me_model.reactions.get_by_id(met_demand_rxn).add_metabolites({pe_pg_constraint: -cur_pe_pg_coeff*cur_met_coeff, clpn_pg_constraint: -cur_clpn_pg_coeff*cur_met_coeff})
        elif 'clpn' == me_met_obj.id[:4]: #add clpn_pg_constraint
            cur_clpn_pg_coeff = phospholipid_mol_fraction['pg'] if constraint_type != 'mass' else phospholipid_mol_fraction['pg'] * me_met_obj.formula_weight
            me_model.reactions.get_by_id(met_demand_rxn).add_metabolites({clpn_pg_constraint: cur_clpn_pg_coeff*cur_met_coeff})

def add_membrane_protein_total_SA_constraint(me_model, query_name, constraint_name):
    #add the protein total surface area constraint based on the coefficient from iJL1678
    #query_name can be 'Inner' or 'Outer'
    lipo_protein_SA = 50; b4513_SA = 175.896144
    if me_model.metabolites.has_id(constraint_name): return
    cur_constraint = Metabolite(constraint_name)
    unique_genes = np.unique([cur_met.id.split('_')[1] for cur_met in me_model.metabolites.query(query_name)])
    for cur_gene in unique_genes:
        if me_model.metabolites.has_id('protein_' + str(cur_gene) + '_lipoprotein_' + query_name + '_Membrane'):
            protein_name = 'protein_' + str(cur_gene) + '_lipoprotein_' + query_name + '_Membrane'
            cur_coeff = lipo_protein_SA
        elif me_model.metabolites.has_id('protein_' + str(cur_gene) + '_' + query_name + '_Membrane'):
            protein_name = 'protein_' + str(cur_gene) + '_' + query_name + '_Membrane'
            cur_coeff = b4513_SA if cur_gene == 'b4513' else me_model.metabolites.get_by_id(protein_name).formula_weight*1.21/42*2 #protein occupies both sides of the membrane
        else:
            continue

        for rxn in me_model.metabolites.get_by_id(protein_name).reactions:
            if rxn.get_coefficient(protein_name) > 0: #reactions that produce the protein, it will contribute to the membrane surface area
                rxn.add_metabolites({cur_constraint: cur_coeff*abs(rxn.get_coefficient(protein_name))})

def add_membrane_lipid_total_SA_constraint(me_model, constraint_name):
    #add the lipid total surface area constraint based on the coefficient from iJL1678
    if me_model.metabolites.has_id(constraint_name):
        cur_constraint = me_model.metabolites.get_by_id(constraint_name)
    else:
        cur_constraint = Metabolite(constraint_name)
    total_demand_rids = [cur_rxn.id for cur_rxn in me_model.reactions.query('Demand')] 
    if 'inner' in constraint_name: #for inner membrane
        target_demand_rids = [cur_rid for cur_rid in total_demand_rids if '_c' == cur_rid[-2:] or 'clpn' in cur_rid]
    else: #for outer membrane
        target_demand_rids = [cur_rid for cur_rid in total_demand_rids if '_c' != cur_rid[-2:] and 'clpn' not in cur_rid]
    for cur_demand_rid in target_demand_rids:
        cur_lipid_met_name = '_'.join(cur_demand_rid.split('_')[1:])
        cur_rxn = me_model.reactions.get_by_id(cur_demand_rid)
        cur_lipid_type = get_lipid_type(cur_lipid_met_name)
        cur_rxn.add_metabolites({cur_constraint: phospholipid_SA_dict[cur_lipid_type]*abs(cur_rxn.get_coefficient(cur_lipid_met_name))})

def get_cell_SA(growth_rate):
    #get the cell surface area as a function growth rate mu
    l=lambda mu: 1.5*2.6*2**(mu*log(2)/3)
    r=lambda mu: 1.5*0.15204137*2**(mu*log(2)/3)
    sa=lambda mu: 2*pi*r(mu)*(l(mu)-2*r(mu))+4*pi*r(mu)**2
    return sa(growth_rate)

def get_gDW_per_cell(growth_rate):
    #get gDW per cell as a function of growth rate
    l=lambda mu: 1.5*2.6*2**(mu*log(2)/3)
    r=lambda mu: 1.5*0.15204137*2**(mu*log(2)/3)
    v=lambda mu: (l(mu)-2*r(mu))*pi*r(mu)**2 + 4./3*pi*r(mu)**3 #previous version had 4/3 not returning the float number
    gDW_per_vol=(0.3*1.105)*1e-12 #70% water, density is g/mL, to convert to um^3, need to 1/((10^4)^3)
    gDW_per_cell = gDW_per_vol*v(growth_rate)
    return gDW_per_cell

def get_SA_formulation(growth_rate):
    #get the formulation for the total surface area to add to the ME model as constraint
    total_SA = 2*get_cell_SA(growth_rate)*10**8 #surface area has two sides, and we convert unit from um to angstrom
    gDW_per_cell = get_gDW_per_cell(growth_rate)
    avogadro=6.02214e23
    #SUM(flux/1000.*avogadro*gDW_per_cell/mu*SA_individual) = total_SA
    #we extract the coefficient and just have
    #SUM(flux*SA_individual) = total_SA*mu/gDW_per_cell/avogadro*1000
    #(mmol/gDW/hour)*(1/1000)*(/mol)*(g)/(1/hour)*(angstrom^2) = (um^2)*(angstrom^2/um^2)
    return total_SA*growth_rate/gDW_per_cell/avogadro*1000

def get_protein_total_SA(me_model, solution_dict, query_name, constraint_name):
    #get the total surface area contributed from protein based on the numerical solution of ME model
    total_protein_SA = 0.0
    unique_genes = np.unique([cur_met.id.split('_')[1] for cur_met in me_model.metabolites.query(query_name)])
    for cur_gene in unique_genes:
        if me_model.metabolites.has_id('protein_' + str(cur_gene) + '_lipoprotein_' + query_name + '_Membrane'):
            protein_name = 'protein_' + str(cur_gene) + '_lipoprotein_' + query_name + '_Membrane'
        elif me_model.metabolites.has_id('protein_' + str(cur_gene) + '_' + query_name + '_Membrane'):
            protein_name = 'protein_' + str(cur_gene) + '_' + query_name + '_Membrane'
        else:
            continue
        counter = 0
        for rxn in me_model.metabolites.get_by_id(protein_name).reactions:
            if rxn.get_coefficient(protein_name) > 0 and abs(solution_dict[rxn.id]) > 0:
                cur_SA = rxn.get_coefficient(constraint_name) * solution_dict[rxn.id]
                if cur_SA < 0:
                    print 'negative flux encountered'
                else:
                    total_protein_SA += cur_SA
                    counter += 1
        if counter > 1:
            print protein_name, ' contributes to SA in more than two reactions'
    return total_protein_SA

def get_lipid_total_SA(me_model, solution_dict, constraint_name):
    #get the total surface area of contributed from the lipid based on the numerical solution of ME model
    total_lipid_SA = 0.0
    total_demand_rids = [cur_rxn.id for cur_rxn in me_model.reactions.query('Demand')] 
    if 'inner' in constraint_name: #for inner membrane
        target_demand_rids = [cur_rid for cur_rid in total_demand_rids if '_c' == cur_rid[-2:] or 'clpn' in cur_rid]
    else: #for outer membrane
        target_demand_rids = [cur_rid for cur_rid in total_demand_rids if '_c' != cur_rid[-2:] and 'clpn' not in cur_rid]
        
    for cur_demand_rid in target_demand_rids:
        if solution_dict[cur_demand_rid] < 0:
            print 'negative flux encountered'
        elif solution_dict[cur_demand_rid] > 0:
            total_lipid_SA += me_model.reactions.get_by_id(cur_demand_rid).get_coefficient(constraint_name) * solution_dict[cur_demand_rid]
    return total_lipid_SA

def add_membrane_constraint(me_model, cur_strain_name):
    #add lipid demand and lipid fraction constraint
    lipid_composition_df = get_lipid_composition(me_model, fatty_acid_comp_df, cur_strain_name)
    add_lipid_demand(me_model, lipid_composition_df)
    add_lipid_fraction_constraint(me_model, lipid_composition_df, constraint_type='mol')
    #add total surface area constraint on protein and lipid
    add_membrane_protein_total_SA_constraint(me_model, 'Inner', 'inner_membrane_size')
    add_membrane_protein_total_SA_constraint(me_model, 'Outer', 'outer_membrane_size')
    add_membrane_lipid_total_SA_constraint(me_model, 'inner_membrane_size')
    add_membrane_lipid_total_SA_constraint(me_model, 'outer_membrane_size')
    me_model.metabolites.inner_membrane_size._bound = get_SA_formulation(mu)
    me_model.metabolites.outer_membrane_size._bound = get_SA_formulation(mu)