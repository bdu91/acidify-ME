#import modules
from qminospy.me1 import ME_NLP1
from cobrame.io.json import load_json_me_model
from periplasmic_proteome import *
from membrane_lipid_constraint import *
from membrane_protein_activity import *
from proton_influx import *

#load the ME-model
me_model = load_json_me_model('iJL1678b.json')

#add different constraints into the ME-model to simulate acid stress, here we are obtaining a solution under pH 5.5
add_membrane_constraint(me_model, 'MJR_habituated')
add_protein_stability_constraint(me_model, 5.5, constrain_periplasm_biomass=False)
add_proton_leak_rxn(me_model, 5.5)
modify_membrane_protein_activity(me_model, 7.0)

#solve the model
me_nlp = ME_NLP1(me_model, growth_key='mu')
muopt, hs, xopt, cache = me_nlp.bisectmu(precision=1e-6, mumax=1.5)
me_model.solution.f = me_model.solution.x_dict['biomass_dilution']