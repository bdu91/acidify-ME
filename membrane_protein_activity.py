def modify_membrane_protein_activity(me_model, fold_change):
    #just need to update ATP synthesis rate here
    ATP_synthase_name = 'ATPS4rpp_FWD_ATPSYN-CPLX_mod_mg2'
    me_model.reactions.get_by_id(ATP_synthase_name).keff = me_model.reactions.get_by_id(ATP_synthase_name).keff * fold_change
    me_model.reactions.get_by_id(ATP_synthase_name).update()