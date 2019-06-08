from numpy.polynomial.polynomial import polyfit
import numpy as np
from cobrame.core.reaction import MEReaction

def get_pH_in(pH_out):
    #pHi = 7.6 + 0.1*(pHe - 7.6)
    #applicable to pH range between 5 to 9
    return 7.6 + 0.1 * (pH_out - 7.6)
    
def get_membrane_conductance(pH_out, degree = 3):
    #default fitted with a 3rd order polynomial
    #data under aerobic respiration
    measured_pH_out = [4.496389892, 5.496389892, 6.5, 7.5]
    #the conductance is measured in unit umol H+/sec/pH/gDW, the actual flux is mmol/hour/gDW, so need to multiply 3.6 to get mmol H+/hour/pH/gDW
    measured_cond = [1.25648855, 0.842748092, 0.79389313, 0.783206107]
    fit_coef = polyfit(measured_pH_out, measured_cond, degree)
    pH_out_powered = [np.power(pH_out,p) for p in range(degree + 1)]
    return np.dot(fit_coef, pH_out_powered)*3.6

def get_proton_leak_rate(pH_out):
    pH_in = get_pH_in(pH_out)
    delta_pH = pH_in - pH_out
    membrane_C = get_membrane_conductance(pH_out)
    return membrane_C*(pH_in-pH_out)

#add in the proton leak reaction
def add_proton_leak_rxn(me_model, pH_out, ko = False):
    rxn_lb = get_proton_leak_rate(pH_out)
    if len(me_model.reactions.query('H_leak')) == 0: #need to add the proton leak reaction
        H_leak = MEReaction('H_leak')
        H_leak.add_metabolites({me_model.metabolites.h_p:-1, me_model.metabolites.h_c:1})
        me_model.add_reaction(H_leak)
    me_model.reactions.H_leak.lower_bound = rxn_lb
    if ko:
        me_model.reactions.H_leak.knock_out()