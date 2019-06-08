from numpy.polynomial.polynomial import polyfit
from glob import glob
import pandas as pd
import numpy as np
from Bio import SeqIO
import math
import os
import cobrame
from cobrame.core.processdata import PostTranslationData
from cobrame.core.component import ProcessedProtein, Metabolite
from cobrame.core.reaction import PostTranslationReaction, MEReaction
from cobrame.util import mu, building


peptide_rog_df = pd.read_csv('data/peptide_radius_of_gyration.csv')
periplasm_fold_rate_df = pd.read_csv('data/periplasm_protein_fold_rate.csv', index_col=0)
charged_aa_pKa_df = pd.read_csv('data/charged_aa_side_chain_pKa.csv')
DATA_DIR = 'data/'; all_proteins_dir = 'proteins/'

def get_protein_net_charge_table():
	#get the charge of folded proteins at different pHs
	proteins_w_charge_calc = []
	for cur_gene_path in glob(DATA_DIR + all_proteins_dir + '*'): #periplasm_genes
	    if os.path.exists(cur_gene_path + '/sum_crg.out'):
	        proteins_w_charge_calc.append(cur_gene_path.split('/')[-1])
	#define the dataframe with the pH as the index and gene id as the column
	protein_net_charge_df = pd.DataFrame(index=np.arange(0, 14.5, 0.5), columns=proteins_w_charge_calc)

	#now read the charge output file for each gene and fill in the table
	for cur_gene in proteins_w_charge_calc:
	    cur_file = open(DATA_DIR + all_proteins_dir + cur_gene + '/sum_crg.out')
	    for line in cur_file:
	        if 'Net_Charge' not in line: continue
	        protein_net_charge_df[cur_gene] = map(float, line.split()[1:])
	    cur_file.close()
	return protein_net_charge_df

protein_net_charge_df = get_protein_net_charge_table()

def calc_peptide_Rg(residue_num):
    '''
    calculate based on residue length using the best fit line from data
    data as indicated on the source column in the excel sheet
    http://biostat.jhsph.edu/~iruczins/presentations/ruczinski.05.03.rutgers.pdf
    Principles of polymer chemistry, by Flory 1953 is where the scaling exponent 0.588 comes from
    '''
    residual_power = 0.588
    residue_num_list = [np.power(cur_num, residual_power) for cur_num in peptide_rog_df.Corrected_residue.values]
    rog_list = peptide_rog_df.Rg_u.values
    fit_coef = polyfit(residue_num_list, rog_list, 1)
    residue_num_powered = [np.power(residue_num, residual_power*p) for p in range(2)] #only first order
    return np.dot(fit_coef, residue_num_powered)

def get_AA_seq_from_pdb(pdb_file_path):
    #read pdb file and return the one letter code amino acid sequence code for protein
    aa_letter_dict = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',\
    'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
    input_file = open(pdb_file_path)
    seq = ''; prev = '0'
    for line in input_file:
        toks = line.split()
        if len(toks) < 1: continue
        if toks[0] != 'ATOM': continue
        if toks[5] != prev:
            seq += aa_letter_dict[toks[3]]
            prev = toks[5]
    input_file.close()
    return seq

def calc_tot_folding_energy(T, aa_num):
    #the folding energy calculated here is assumed to be pH 7
    #T should be in kelvin
    dH = 4.0*aa_num + 143.0 #in unit kJ/mol
    dS = (13.27*aa_num + 448.0)/1000 #in unit kJ/mol
    dCp = 0.049*aa_num + 0.85 #in unit kJ/mol
    T_h = 373.5; T_s = 385
    dG = dH + dCp*(T - T_h) - T*dS - T*dCp*np.log(T/T_s)
    #dG2 = (4.0*aa_num + 143.0) + (0.049*aa_num + 0.85)*(T-T_h) - T*(0.01327*aa_num + 0.448) - T*(0.049*aa_num + 0.85)*np.log(T/T_s)
    return -dG #G_folded - G_unfolded, in kJ/mol

def calc_tot_folding_energy_from_b_num(T, b_num):
    return calc_tot_folding_energy(T, len(get_AA_seq_from_pdb(get_pdb_file_path(b_num))))

def calc_protein_Rg(pdb_file_path):
    '''
    Calculates the Radius of Gyration (Rg) of a protein given its .pdb 
    structure file. Returns the Rg in Angstrom.
    This function is adapted from https://github.com/sarisabban/Rg/blob/master/Rg.py
    '''
    coord = list()
    mass = list()
    Structure = open(pdb_file_path, 'r')
    for line in Structure:
        try:
            line = line.split()
            x = float(line[6])
            y = float(line[7])
            z = float(line[8])
            coord.append([x, y, z])
            if line[-1] == 'C':
                mass.append(12.0107)
            elif line[-1] == 'O':
                mass.append(15.9994)
            elif line[-1] == 'N':
                mass.append(14.0067)
            elif line[-1] == 'S':
                mass.append(32.065)
        except:
            pass
    xm = [(m*i, m*j, m*k) for (i, j, k), m in zip(coord, mass)]
    tmass = sum(mass)
    rr = sum(mi*i + mj*j + mk*k for (i, j, k), (mi, mj, mk) in zip(coord, xm))
    mm = sum((sum(i) / tmass)**2 for i in zip(*xm))
    rg = math.sqrt(rr / tmass-mm)
    Structure.close()
    return rg

def get_aa_charge(aa_one_letter, pH):
    #get the charge for the individual amino acid, based on individual aa pKa
    if aa_one_letter not in charged_aa_pKa_df.one_letter.values: return 0
    cur_aa_info = charged_aa_pKa_df[charged_aa_pKa_df.one_letter == aa_one_letter]
    cur_aa_pKa = cur_aa_info.pKa.values[0]
    cur_aa_type = cur_aa_info.type.values[0]
    pH_pKa_val = np.power(10, pH - cur_aa_pKa)
    if cur_aa_type == 'acid':
        return -pH_pKa_val/(1 + pH_pKa_val)
    elif cur_aa_type == 'base':
        return 1/(1 + pH_pKa_val)

def get_peptide_charge(aa_one_letter_seq, pH):
    #get the charge for the whole peptide sequence, based on individual aa pKa
    acid_charge = lambda pKa: np.power(10, pH - pKa)/(1 + np.power(10, pH - pKa))
    base_charge = lambda pKa: np.power(10, pKa - pH)/(1 + np.power(10, pKa - pH))
    total_charge = 0
    for aa_one_letter in aa_one_letter_seq:
        if aa_one_letter not in charged_aa_pKa_df.one_letter.values: continue
        cur_aa_info = charged_aa_pKa_df[charged_aa_pKa_df.one_letter == aa_one_letter]
        cur_aa_pKa = cur_aa_info.pKa.values[0]
        cur_aa_type = cur_aa_info.type.values[0]
        if cur_aa_type == 'acid':
            total_charge -= acid_charge(cur_aa_pKa)
        elif cur_aa_type == 'base':
            total_charge += base_charge(cur_aa_pKa)
    return total_charge

def get_protein_charge(bnum, pH_to_predict, degree = 15):
    #use polynomial fit based on the existing pH profile, overfitting is fine here as we are trying to draw a smooth curve
    if type(pH_to_predict) != list:
        pH_to_predict = [pH_to_predict]
    pH_list = protein_net_charge_df.index
    charge_list = protein_net_charge_df[bnum].values
    fit_coef = polyfit(pH_list, charge_list, degree)
    pH_out_powered = [np.power(pH_to_predict,p) for p in range(degree + 1)]
    charge_predicted = np.dot(fit_coef, pH_out_powered)
    if len(charge_predicted) == 1:
        return charge_predicted[0]
    else:
        return charge_predicted

def get_pdb_file_path(b_num):
    return DATA_DIR + all_proteins_dir + '%s/prot.pdb' %b_num

def calc_e_folding_energy(b_num, pH, T=310.15, IS=0.25):
    diele_water = 80.4; R = 8.314
    l_b = 1.39/10**4/diele_water/R/T #in meter
    k = np.sqrt(2*(IS*1000)*l_b) #convert IS from mol/L to mol/m3
    
    cur_bnum_pdb_path = get_pdb_file_path(b_num)
    cur_peptide_seq = get_AA_seq_from_pdb(cur_bnum_pdb_path)
    
    protein_charge = protein_net_charge_df.at[pH, b_num] if pH in protein_net_charge_df.index else get_protein_charge(b_num, pH)
    protein_Rg = calc_protein_Rg(cur_bnum_pdb_path) / 10.0**10 #get angstrom, need to convert to meter
    
    peptide_charge = get_peptide_charge(cur_peptide_seq, pH)
    peptide_Rg = calc_peptide_Rg(len(cur_peptide_seq)) / 10.0**10 #get angstrom, need to convert to meter
    
    dG_e = R*T*(np.power(protein_charge, 2) * l_b / (2 * protein_Rg * (1 + k * protein_Rg)) - \
                np.power(peptide_charge, 2) * l_b / (2 * peptide_Rg * (1 + k * peptide_Rg)))
    return dG_e/1000 #unit in kJ/mol

def calc_tot_folding_energy_at_pH(b_num, pH, T = 310.15):
    dG_folding_pH7 = calc_tot_folding_energy_from_b_num(T, b_num)
    dGe_folding_pH7 = calc_e_folding_energy(b_num, 7.0, T)
    dGe_folding_pH = calc_e_folding_energy(b_num, float(pH), T)
    return dG_folding_pH7 - dGe_folding_pH7 + dGe_folding_pH #in kJ/mol

def get_fold_rate(seq, secstruct):
    """
    This function is obtained from ssbio package, https://github.com/SBRG/ssbio/blob/master/ssbio/protein/sequence/properties/kinetic_folding_rate.py
    Submit sequence and structural class to FOLD-RATE calculator (http://www.iitm.ac.in/bioinfo/fold-rate/)
    to calculate kinetic folding rate.
    Args:
        seq: Amino acid sequence in string format
        secstruct (str): Structural class: `all-alpha``, ``all-beta``, ``mixed``, or ``unknown``
    Returns:
        float: Kinetic folding rate k_f
    """

    url = 'http://www.iitm.ac.in/bioinfo/cgi-bin/fold-rate/foldrateCalculator.pl'

    values = {'sequence': seq, 'eqn': secstruct}
    data = urllib.urlencode(values)
    data = data.encode('ASCII')
    response = urllib.urlopen(url, data)

    result = str(response.read())
    ind = str.find(result, 'The folding rate,')
    result2 = result[ind:ind + 70]
    ind1 = str.find(result2, '=')
    ind2 = str.find(result2, '/sec')
    rate = result2[ind1 + 2:ind2]

    return rate

def merge_list_of_list(list_of_list):
    return [item for sublist in list_of_list for item in sublist]

def get_chaperone_translocation_rxns(complex_name, me_model):
    return list(set([cur_rxn for cur_rxn in merge_list_of_list([met.reactions for met in me_model.metabolites.query(complex_name)]) if 'translocation_' in cur_rxn.id]))

def add_periplasm_protein_folded(me_model):
    #add periplasm proteins into the model and make corresponding changes in reactions
    unique_Tat_complex_translocation_rxns = get_chaperone_translocation_rxns('Tat', me_model)
    unique_Sec_complex_translocation_rxns = get_chaperone_translocation_rxns('Sec', me_model)
    periplasm_proteins_in_Tat_pathway = list(set([cur_protein.id for cur_protein in merge_list_of_list([rxn.products for rxn in unique_Tat_complex_translocation_rxns]) if 'Periplasm' in cur_protein.id]))
    periplasm_proteins_in_Sec_pathway = list(set([cur_protein.id for cur_protein in merge_list_of_list([rxn.products for rxn in unique_Sec_complex_translocation_rxns]) if 'Periplasm' in cur_protein.id]))
    bnum_in_Tat_pathway = [cur_protein.split('_')[1] for cur_protein in periplasm_proteins_in_Tat_pathway]
    bnum_in_Sec_pathway = [cur_protein.split('_')[1] for cur_protein in periplasm_proteins_in_Sec_pathway]
    
    #add folded periplasm protein objects into the model
    for cur_protein in me_model.metabolites.query('Periplasm'):
        cur_folded_pid = cur_protein.id + '_folded'
        folded_protein = ProcessedProtein(cur_folded_pid, cur_protein.id)
        me_model.add_metabolites([folded_protein])
    
    #the proteins translocated by Tat are already folded, so we need to modify change all the reactions they are in as folded protein
    for cur_protein_id in periplasm_proteins_in_Tat_pathway:
        cur_protein_folded = me_model.metabolites.get_by_id(cur_protein_id + '_folded')
        cur_protein_rxns = me_model.metabolites.get_by_id(cur_protein_id).reactions
        for cur_rxn in cur_protein_rxns:
            cur_protein_coeff = cur_rxn.pop(cur_protein_id)
            cur_rxn.add_metabolites({cur_protein_folded: cur_protein_coeff})
            
    #adjust the folding event occuring in cytoplasm for Tat and Sec assisted pathway proteins
    for data in me_model.translation_data:
        #add folding reactions to proteins translocated by Tat
        if data.id in bnum_in_Tat_pathway:
            data.subreactions['GroEL_dependent_folding'] = 1
        #remove folding reactions from proteins translocated by Sec
        if data.id in bnum_in_Sec_pathway:
            for subreaction in list(data.subreactions.keys()):
                if 'folding' in subreaction:
                    data.subreactions.pop(subreaction)

    #for proteins in Sec pathway, we need to modify their complex formation reactions to have folded protein
    #for translocation reaction it will still be unfolded protein
    for cur_protein_id in periplasm_proteins_in_Sec_pathway:
        cur_protein_folded = me_model.metabolites.get_by_id(cur_protein_id + '_folded')
        cur_protein_rxns = me_model.metabolites.get_by_id(cur_protein_id).reactions
        for cur_rxn in cur_protein_rxns:
            if 'formation' in cur_rxn.id:
                cur_protein_coeff = cur_rxn.pop(cur_protein_id)
                cur_rxn.add_metabolites({cur_protein_folded: cur_protein_coeff}) 
                
def add_periplasm_protein_folding_reaction(me_model, pH, T = 310.15):
    #this function allows repeatedly modifying the folding reaction given different pH
    #https://github.com/SBRG/ssbio/tree/master/ssbio/protein/sequence/properties
    R = 8.314
    for cur_protein in me_model.metabolites.query('Periplasm'):
        if 'folded' in cur_protein.id or 'b3509' in cur_protein.id: continue #just handle unfolded protein
        cur_folded_pid = cur_protein.id + '_folded'
        try:
            folded_protein = me_model.metabolites.get_by_id(cur_folded_pid)
        except KeyError:
            folded_protein = ProcessedProtein(cur_folded_pid, cur_protein.id)
            me_model.add_metabolites([folded_protein])
        folding_id = 'folding_' + cur_protein.id + '_folding_spontaneous'
        try:
            folding_rxn = me_model.reactions.query(folding_id)[0]
        except IndexError:
            folding_rxn = PostTranslationReaction(folding_id)
            me_model.add_reaction(folding_rxn)
        folding_rxn.clear_metabolites() #first remove all previous metabolites
        
        cur_bnum = cur_protein.id.split('_')[1]
        if cur_bnum not in protein_net_charge_df.columns:
            cur_folding_energy = -100000 #assume favorable towards folding
        else:
            cur_folding_energy = calc_tot_folding_energy_at_pH(cur_bnum, pH, T) * 1000 #to make the unit J/mol
        keq_folding = 1.0/np.exp(cur_folding_energy/(-R*T)) #unfolded/folded
        
        k_folding = periplasm_fold_rate_df.at[cur_bnum, 'mixed']
        if k_folding < 0:
            k_folding = max([periplasm_fold_rate_df.at[cur_bnum, cur_struct] for cur_struct in periplasm_fold_rate_df.columns]) #get a positive value
        k_folding = k_folding * 3600 #in sec-1 originally, to get to hr-1 need to multiply by 3600
        
        dilution = keq_folding + mu / k_folding
        folding_rxn.add_metabolites({cur_protein: -(dilution + 1.0), folded_protein: 1.0})
        #need to add Keq_folding * HdeB, and the corresponding complex
        folding_rxn.lower_bound = -1000.0 #make it reversible
        folding_rxn.notes = {'keq':keq_folding} #store the current keq for later use


def add_HdeB_complex(me_model):
    # Get the sequence
    b3509_seq = SeqIO.read('data/b3509.fasta', 'fasta')
    nucleo_seq = str(b3509_seq.seq).upper()
    
    # Add transcript to model
    building.create_transcribed_gene(me_model, 'b3509', 'mRNA', nucleo_seq, strand = '-')
    
    # Link transcription data to reaction
    building.add_transcription_reaction(me_model, 'TU_b3509', {'b3509'}, nucleo_seq)
    me_model.process_data.TU_b3509.RNA_polymerase = ''
    me_model.process_data.TU_b3509.subreactions = {}
    me_model.reactions.transcription_TU_b3509.update()
    
    # Add translation reaction
    building.add_translation_reaction(me_model, 'b3509', nucleo_seq)
    me_model.process_data.b3509.add_initiation_subreactions(start_codons=set(), start_subreactions=set())
    me_model.process_data.b3509.add_elongation_subreactions(elongation_subreactions=set())
    me_model.process_data.b3509.add_termination_subreactions(translation_terminator_dict=None)
    me_model.reactions.translation_b3509.update()
    
    #add translocation pathway
    me_model.add_metabolites([ProcessedProtein('protein_b3509_Periplasm', 'protein_b3509')])
    translocation_data = PostTranslationData('translocation_b3509', me_model, 'protein_b3509_Periplasm', 'protein_b3509')
    translocation_rxn = PostTranslationReaction('translocation_b3509')
    me_model.add_reaction(translocation_rxn)
    translocation_rxn.posttranslation_data = translocation_data
    translocation_rxn.update()    
    
    #now add the complex formation reaction, HdeB is a dimer
    complex_data = cobrame.ComplexData('complex_HdeB', me_model)
    complex_data.stoichiometry = {'protein_b3509_Periplasm': 2}
    complex_data.subreactions = {} 
    # complex data has utility to link data and add reaction
    complex_data.create_complex_formation()
    
def add_periplasma_chaperone_protection(me_model):
    '''
    add HdeB protection to the folding reaction, forms HdeB-denature_protein complex, also add a dilution
    reaction for the complex, this function can be reused
    '''
    chaperon_name = 'HdeB'
    for cur_reaction in me_model.reactions.query('Periplasm_folding_spontaneous'):
        cur_bnum = cur_reaction.id.split('_')[2]
        cur_protein_id = 'protein_' + cur_bnum + '_Periplasm'
        cur_complex_id = cur_protein_id + '_' + chaperon_name + '_complex'
        
        try:
            chaperone_complex = me_model.metabolites.get_by_id(cur_complex_id)
        except KeyError:
            chaperone_complex = ProcessedProtein(cur_complex_id, cur_protein_id)
            me_model.add_metabolites([chaperone_complex])
        try:
            cur_reaction.pop('complex_' + chaperon_name)
            cur_reaction.pop(cur_complex_id)
        except KeyError:
            pass
        
        #add the chaperone reaction into the current folding reactions
        cur_reaction.add_metabolites({me_model.metabolites.get_by_id('complex_' + chaperon_name): -cur_reaction.notes['keq'], chaperone_complex: cur_reaction.notes['keq']})
        if len(me_model.reactions.query(cur_complex_id + '_dilution')) == 0:
            complex_dilution_rxn = MEReaction(cur_complex_id + '_dilution')
            complex_dilution_rxn.add_metabolites({me_model.metabolites.get_by_id(cur_complex_id): -1.0})
            me_model.add_reaction(complex_dilution_rxn)
        
def add_periplasm_biomass_constraint(me_model, modeled_protein_fraction, periplasm_protein_fraction, unmodeled_protein_fraction):
    if len(me_model.metabolites.query('periplasm_protein_biomass')) == 0:
        #first change the biomass dilution for periplasm proteins in the translation reactions
        periplasm_protein_genes = list(set([met.id.split('_')[1] for met in me_model.metabolites.query('Periplasm')])) #remove duplicates
        periplasm_protein_biomass = Metabolite('periplasm_protein_biomass')
        for cur_gene in periplasm_protein_genes:
            cur_translation_rxn = me_model.reactions.get_by_id('translation_' + cur_gene)
            cur_protein_biomass_coeff = cur_translation_rxn.pop('protein_biomass')
            cur_translation_rxn.add_metabolites({periplasm_protein_biomass: cur_protein_biomass_coeff})
    
    #now modify the protein biomass and biomass reaction
    protein_biomass_to_biomass_rxn_mets = me_model.reactions.protein_biomass_to_biomass.metabolites
    #pop out the original metabolties
    for cur_met in protein_biomass_to_biomass_rxn_mets:
        me_model.reactions.protein_biomass_to_biomass.pop(cur_met)
    me_model.reactions.protein_biomass_to_biomass.add_metabolites({me_model.metabolites.protein_biomass: -modeled_protein_fraction, \
    me_model.metabolites.periplasm_protein_biomass: -periplasm_protein_fraction, me_model.metabolites.unmodeled_protein_biomass: -unmodeled_protein_fraction, \
    me_model.metabolites.biomass: modeled_protein_fraction + periplasm_protein_fraction + unmodeled_protein_fraction})

def add_protein_stability_constraint(me_model, pH, constrain_periplasm_biomass=True):
    add_periplasm_protein_folded(me_model)
    add_periplasm_protein_folding_reaction(me_model, pH)
    add_HdeB_complex(me_model) #add HdeB later so that b3509 (HdeB gene) itself is not bound by HdeB
    add_periplasma_chaperone_protection(me_model)
    if constrain_periplasm_biomass:
        add_periplasm_biomass_constraint(me_model, 0.65, 0.1, 0.25)