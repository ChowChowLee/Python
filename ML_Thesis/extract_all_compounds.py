import sys
import numpy as np
import pandas as pd
import getpass

import rdkit
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

from enviPath import *
from objects import *

# sys.path.insert(0, "C:\\Users\\leetseng\\enviPath-python\\enviPath_python\\")
# sys.path.insert(1, "C:\\Users\\leetseng\\enviPath-python\\")

sys.path.insert(0, '..//')
file_location = ("C://Users//Lee//PycharmProjects//pythonProjectAlpacaBaby//Python//ML_Thesis")
# input_file_path = file_location+'input/sludge_compounds_final.txt'
output_file_path_full = file_location+'//output//sludge_Rich_raw.tsv'

LEO_SLUDGE = 'https://envipath.org/package/195bc500-f0c6-4bcb-b2fe-f1602b5f20a2'
RICH_SLUDGE = 'https://envipath.org/package/8d3d7ca2-ae4e-4779-a6a2-d3539237c439'
ORIGINAL_SLUDGE = 'https://envipath.org/package/4a3cd0f4-4d2b-4f00-b3e6-a29e721f7038'
INSTANCE_HOST = 'https://envipath.org'
USERNAME = 'leetseng'
PASSWORD = getpass.getpass(prompt='Enter your password:')

def main():

    D = {
         "scenario_id": [], "compound_id": [], "compound_name": [], "smiles": [], "reduced_smiles": [],
         "halflife_raw": [], "halflife_unit": [], "halflife_model_TF": [], "halflife_comment": [],
         "rateconstant": [], "rateconstant_unit": [], "rateconstant_comment": [], "halflife_model": [],
         "acidity": [], "acidity_unit": [],
         "temperature": [], "temperature_unit": [],
         "original_sludge_amount": [], "original_sludge_amount_unit": [],
         "sludge_retention_time": [], "sludge_retention_time_unit": [], "sludge_retention_time_type": [],
         "total_suspended_solids_concentration_start": [], "total_suspended_solids_concentration_end": [], "total_suspended_solids_concentration_unit": [],
         # "volatile_suspended_solids_concentration_start": [], "volatile_suspended_solids_concentration_end": [], "volatile_suspended_solids_concentration_unit": [],
         "addition_of_nutrients": [], "biological_treatment_technology": [],
         "bioreactor_type": [], "bioreactor_value": [], "bioreactor_value_unit": [],
         "nitrogen_content_type": [], "nitrogen_content_influent": [],
         "oxygen_demand_type": [], "oxygen_demand_value": [],
         "oxygen_uptake_rate": [], "oxygen_uptake_rate_unit": [],
         "phosphorus_content": [],
         "redox": [],
         "source_of_liquid_matrix": [],
         "type_of_addition": [],
         "type_of_aeration": [],
         "inoculum_source": [],
         "location": [],
         "purpose_of_wwtp": [],
          }

    i = 0
    for pathway in all_pathways:
        for node in pathway.get_nodes():
            print(node.get_id())  #this is not compound_id, it should be the Node id.
            i += 1
            print("checking node # ", i)
            try:
                scenarios = node.get_scenarios()
            except:
                continue
            else:
                for scenario in scenarios:
                    D = add_scenario_information(D, scenario, node)
    sludge_df = pd.DataFrame.from_dict(D)
    sludge_df.index = np.arange(1, len(sludge_df) + 1)
    sludge_df.index.name = 'index'
    print(sludge_df.describe())
    sludge_df.to_csv(output_file_path_full, mode='w', sep="\t") #a, w, r

def add_scenario_information(D, scenario, node):
    full_scenario = Scenario(eP.requester, id=scenario.get_id())
    add_info = full_scenario.get_additional_information()
    try:
        halflife_object = add_info.get_halflife()
        has_hf = True
    except AttributeError:
        has_hf = False

    try:
        rateconstant_object = add_info.get_rateconstant()
        has_rateconstant = True
    except AttributeError:
        has_rateconstant = False

    if has_hf or has_rateconstant:  ###
        D['compound_id'].append(arcessere_compound_id(node))
        D['compound_name'].append(arcessere_compound_name(node))
        D['smiles'].append(arcessere_smiles(node))
        D['reduced_smiles'].append(canonicalize_smiles(arcessere_smiles(node)))
        # D['reduced_smiles'].append(reduced_smiles)  # cropped_canonical_smiles_no_stereo
        print(full_scenario.get_id())
        D['acidity'].append(arcessere_acidity(add_info))
        D['acidity_unit'].append(arcessere_acidity_unit(add_info))
        D['addition_of_nutrients'].append(arcessere_addition_of_nutrients(add_info))
        D['biological_treatment_technology'].append(arcessere_biological_treatment_technology(add_info))
        D['bioreactor_type'].append(arcessere_bioreactor_type(add_info))
        D['bioreactor_value'].append(arcessere_bioreactor_value(add_info))
        D['bioreactor_value_unit'].append(arcessere_bioreactor_value_unit(add_info))
        # D['confidencelevel'].append(arcessere_confidencelevel(add_info))
        D['halflife_raw'].append(arcessere_halflife(add_info))
        D['halflife_unit'].append(arcessere_halflife_unit(add_info))
        D['halflife_model_TF'].append(arcessere_halflife_model(add_info))
        D['halflife_comment'].append(arcessere_halflife_comment(add_info))
        D['inoculum_source'].append(arcessere_inoculum_source(add_info))
        D['location'].append(arcessere_location(add_info))
        D['nitrogen_content_type'].append(arcessere_nitrogen_content_type(add_info))
        D['nitrogen_content_influent'].append(arcessere_nitrogen_content_influent(add_info))
        D['original_sludge_amount'].append(arcessere_original_sludge_amount(add_info))
        D['original_sludge_amount_unit'].append(arcessere_original_sludge_amount_unit(add_info))
        D['oxygen_demand_type'].append(arcessere_oxygen_demand_type(add_info))
        D['oxygen_demand_value'].append(arcessere_oxygen_demand_value(add_info))
        D['oxygen_uptake_rate_unit'].append(arcessere_oxygen_uptake_rate_unit(add_info))
        D['oxygen_uptake_rate'].append(arcessere_oxygen_uptake_rate(add_info))
        D['phosphorus_content'].append(arcessere_phosphorus_content(add_info))
        D['purpose_of_wwtp'].append(arcessere_purpose_of_wwtp(add_info))
        D['rateconstant'].append(arcessere_rate_constant(add_info))
        D['rateconstant_unit'].append(arcessere_rate_constant_unit(add_info))
        D['halflife_model'].append(arcessere_reaction_order(add_info))
        D['rateconstant_comment'].append(arcessere_rate_constant_comment(add_info))
        D['redox'].append(arcessere_redox(add_info))
        D['scenario_id'].append(scenario.get_id())
        D['sludge_retention_time_type'].append(arcessere_sludge_retention_time_type(add_info))
        D['sludge_retention_time'].append(arcessere_sludge_retention_time(add_info))
        D['sludge_retention_time_unit'].append(arcessere_sludge_retention_time_unit(add_info))
        D['source_of_liquid_matrix'].append(arcessere_source_of_liquid_matrix(add_info))
        D['temperature'].append(arcessere_temperature(add_info))
        D['temperature_unit'].append(arcessere_temperature_unit(add_info))
        D['total_suspended_solids_concentration_start'].append(arcessere_tss_start(add_info))
        D['total_suspended_solids_concentration_end'].append(arcessere_tss_end(add_info))
        D['total_suspended_solids_concentration_unit'].append(arcessere_tss_unit(add_info))
        D['type_of_addition'].append(arcessere_type_of_addition(add_info))
        D['type_of_aeration'].append(arcessere_type_of_aeration(add_info))
        # D['volatile_suspended_solids_concentration_start'].append(arcessere_volatile_ss_start(add_info))
        # D['volatile_suspended_solids_concentration_end'].append(arcessere_volatile_ss_end(add_info))
        # D['volatile_suspended_solids_concentration_unit'].append(arcessere_volatile_ss_unit(add_info))
    return D


def arcessere_compound_id(node):
    try:
        id_from_node = node.get_default_structure().get_id()
    except:
        return ''
    else:
        return id_from_node

def arcessere_compound_name(node):
    try:
        name_from_node = node.get_default_structure().get_name()
    except ValueError:
        return ''
    else:
        return name_from_node             #.split(',')[0]  only pick up one compound name

def arcessere_smiles(node):
    try:
        smiles_from_node = node.get_default_structure().get_smiles()
    except:
        return ''
    else:
        return smiles_from_node

def canonicalize_smiles(smiles_from_node):
    mol = Chem.MolFromSmiles(smiles_from_node) # creates mol object from SMILES
    uncharger = rdMolStandardize.Uncharger() # easier to access
    uncharged = uncharger.uncharge(mol) # protonates or deprotonates the mol object
    new_smiles = rdkit.Chem.rdmolfiles.MolToSmiles(uncharged) # converts mol object to canonical SMILES
    can_smiles = Chem.CanonSmiles(new_smiles)
    return can_smiles

def arcessere_acidity(add_info):
    try:
        raw_pH = add_info.get_acidity().get_value()
    except:
        return np.NaN
    else:
        if ';' in raw_pH:
            if ' - ' in raw_pH:
                pH = range_to_average(raw_pH.split(';')[0])
            else:
                pH = float(raw_pH.split(';')[0])
            return np.round(pH, 1)

def arcessere_acidity_unit(add_info):
    try:
        pH_unit = add_info.get_acidity().get_unit()
    except:
        return ''
    else:
        return pH_unit

def arcessere_addition_of_nutrients(add_info):
    try:
        addition_of_nutrients = add_info.get_additionofnutrients().get_value()
    except:
        return ''
    else:
        return addition_of_nutrients

def arcessere_biological_treatment_technology(add_info):
    try:
        biological_treatment_technology = add_info.get_biologicaltreatmenttechnology().get_value()
    except:
        return ''
    else:
        return biological_treatment_technology

def arcessere_bioreactor_type(add_info):    #########################################
    try:
        bioreactor_type = add_info.get_bioreactor().get_value().split(',')[0]
    except:
        return ''
    else:
        return bioreactor_type


def arcessere_bioreactor_value(add_info):
    try:
        bioreactor = float(add_info.get_bioreactor().get_value().split(',')[1])
    except ValueError:
        return np.NaN
    except:
        return np.NaN
    else:
        return bioreactor

def arcessere_bioreactor_value_unit(add_info):
    try:
        bioreactor_unit = add_info.get_bioreactor().get_unit()
    except:
        return ''
    else:
        return bioreactor_unit

def arcessere_confidencelevel(add_info):
    try:
        confidencelevel = add_info.get_confidencelevel().get_value()
    except:
        return np.NaN
    else:
        return float(confidencelevel)

def arcessere_inoculum_source(add_info):
    try:
        inoculumsource = add_info.get_inoculumsource().get_value()
    except:
        return ''
    else:
        return inoculumsource

def arcessere_location(add_info):
    try:
        location = add_info.get_location().get_value()
    except:
        return ''
    else:
        return location

def arcessere_minormajor(add_info):
    try:
        minormajor = add_info.get_minormajor().get_value()
    except:
        return ''
    else:
        return minormajor

def arcessere_nitrogen_content_type(add_info):
    try:
        nitrogencontent = add_info.get_nitrogencontent().get_value()
    except:
        return ''
    else:
        if '&#8322' in nitrogencontent.split(';')[0]:
            return nitrogencontent.split(';')[0].replace('&#8322', '\u2082')
        elif '&#8323' in nitrogencontent.split(';')[0]:
            return nitrogencontent.split(';')[0].replace('&#8323', '\u2083')
        elif '&#8324' in nitrogencontent.split(';')[0]:
            return nitrogencontent.split(';')[0].replace('&#8324', '\u2084')

def arcessere_nitrogen_content_influent(add_info):
    try:
        nitrogencontent = add_info.get_nitrogencontent().get_value()
    except:
        return np.NaN
    else:
        return float(nitrogencontent.split(';')[1])

def arcessere_original_sludge_amount(add_info):
    try:
        originalsludgeamount = add_info.get_originalsludgeamount().get_value()
    except:
        return np.NaN
    else:
        return originalsludgeamount

def arcessere_original_sludge_amount_unit(add_info):
    try:
        originalsludgeamount_unit = add_info.get_originalsludgeamount().get_unit()
    except:
        return ''
    else:
        return originalsludgeamount_unit

def arcessere_oxygen_demand_type(add_info):
    try:
        oxygendemand = add_info.get_oxygendemand().get_value()
    except:
        return ''
    else:
        return oxygendemand.split(';')[0]

def arcessere_oxygen_demand_value(add_info):
    try:
        oxygendemand = add_info.get_oxygendemand().get_value()
    except:
        return np.NaN
    else:
        return float(oxygendemand.split(';')[1])

def arcessere_oxygen_uptake_rate(add_info):
    try:
        our = add_info.get_oxygenuptakerate.get_value()
    except:
        return np.NaN
    else:
        return range_to_average(our)
def arcessere_oxygen_uptake_rate_unit(add_info):
    try:
        sludgeretentiontime_unit = add_info.get_oxygenuptakerate().get_unit()
    except:
        return ''
    else:
        # return sludgeretentiontime_unit
        if "&#8315&#185" in sludgeretentiontime_unit.split(' ')[1] and "&#8315&#185" in sludgeretentiontime_unit.split(' ')[2]:
            return sludgeretentiontime_unit.split(' ')[0] + '/(L * h)'

def arcessere_phosphorus_content(add_info):
    try:
        phosphoruscontent = add_info.get_phosphoruscontent().get_value()
    except:
        return np.NaN
    else:
        return phosphoruscontent.split(';')[0]

def arcessere_proposed_intermediate(add_info):
    try:
        proposedintermediate = add_info.get_proposedintermediate().get_value()
    except:
        return ''
    else:
        return proposedintermediate

def arcessere_purpose_of_wwtp(add_info):
    try:
        purposeofwwtp = add_info.get_purposeofwwtp().get_value()
    except:
        return ''
    else:
        return purposeofwwtp

def arcessere_redox(add_info):
    try:
        redox = add_info.get_redox().get_value()
    except:
        return ''
    else:
        return redox

def arcessere_sludge_retention_time_type(add_info):
    try:
        sludgeretentiontime = add_info.get_sludgeretentiontime().get_value()
    except:
        return ''
    else:
        return sludgeretentiontime.split(';')[0]

def arcessere_sludge_retention_time(add_info):
    try:
        sludgeretentiontime = add_info.get_sludgeretentiontime().get_value()
    except:
        return np.NaN
    else:
        return float(sludgeretentiontime.split(';')[1])

def arcessere_sludge_retention_time_unit(add_info):
    try:
        sludgeretentiontime_unit = add_info.get_sludgeretentiontime().get_unit()
    except:
        return ''
    else:
        return sludgeretentiontime_unit

def arcessere_source_of_liquid_matrix(add_info):
    try:
        sourceofliquidmatrix = add_info.get_sourceofliquidmatrix().get_value()
    except:
        return ''
    else:
        return sourceofliquidmatrix
def arcessere_tss_start(add_info):
    try:
        tts = add_info.get_tts().get_value()
    except:
        return np.NaN
    else:
        return float(tts.split(' _ ')[0].split(' - ')[0])

def arcessere_tss_end(add_info):
    try:
        tts = add_info.get_tts().get_value()
    except:
        return np.NaN
    else:
        return float(tts.split(' _ ')[0].split(' - ')[1])

def arcessere_tss_unit(add_info):
    try:
        tts_unit = add_info.get_tts().get_unit()
    except:
        return ''
    else:
        return tts_unit

def arcessere_type_of_addition(add_info):
    try:
        typeofaddition = add_info.get_typeofaddition().get_value()
    except:
        return ''
    else:
        return typeofaddition
def arcessere_type_of_aeration(add_info):
    try:
        typeofaeration = add_info.get_typeofaeration().get_value()
    except:
        return ''
    else:
        return typeofaeration

def arcessere_volatile_ss_start(add_info):
    try:
        vss_start = add_info.get_volatiletts().get_value()
    except:
        return np.NaN   #not sure if this return type is correct or not
    else:
        if ' - ' in vss_start:
            return float(vss_start.split(' - ')[0])
        else:
            return float(vss_start)

def arcessere_volatile_ss_end(add_info):
    try:
        vss_end = add_info.get_volatiletts().get_value()
    except:
        return np.NaN
    else:
        if ' - ' in vss_end:
            return float(vss_end.split(' - ')[1])
        else:
            return float(vss_end)

def arcessere_volatile_ss_unit(add_info):
    try:
        vss_unit = add_info.get_volatiletts().get_unit()
    except:
        return ''
    else:
        return vss_unit

def range_to_average(input_string):
    if '-' in input_string:
        min = float(input_string.split(' - ')[0])
        max = float(input_string.split(' - ')[1])
        average = np.average([min, max])
    elif ';' in input_string:
        min = float(input_string.split(';')[0])
        max = float(input_string.split(';')[1])
        average = np.average([min, max])
    else:
        average = input_string
    return average

def arcessere_rate_constant(add_info):
    try:
        rate_constant = add_info.get_rateconstant().get_value()
    except:
        return np.NaN
    else:
        min = rate_constant.split(';')[2].split(' - ')[0]
        max = rate_constant.split(';')[2].split(' - ')[1]
        if min != 'NaN' and max != 'NaN':
            average = np.average([float(min), float(max)])
        elif min != 'NaN' and max == 'NaN':
            average = float(min)
        elif min == 'NaN' and max != 'NaN':
            average = float(max)
        return average

def arcessere_rate_constant_unit(add_info):
    try:
        rate_constant_unit = add_info.get_rateconstant().get_unit()
    except:
        return ''
    else:
        if '&#956;g' in rate_constant_unit:
            return rate_constant_unit.replace('&#956;g', '\u338D')
        else:
            return rate_constant_unit

def arcessere_reaction_order(add_info):
    try:
        rate_constant = add_info.get_rateconstant().get_value()
    except:
        return ''
    else:
        return rate_constant.split(';')[0]

def arcessere_rate_constant_comment(add_info):
    try:
        rate_constant_comment = add_info.get_rateconstant().get_value().split(';')[3]
    except:
        return ''
    else:
        return rate_constant_comment

def arcessere_halflife(add_info):
    try:
        hf = add_info.get_halflife().get_value()
    except:
        return np.NaN
    else:
        # return float(hf.split(';')[3].split(' - ')[0])
        min = float(hf.split(';')[3].split(' - ')[0])
        max = float(hf.split(';')[3].split(' - ')[1])
        average = np.average([min, max])
        return average


def arcessere_halflife_unit(add_info):
    try:
        hf_unit = add_info.get_halflife().get_unit()
    except:
        return ''
    else:
        return hf_unit

def arcessere_halflife_model(add_info):
    try:
        hl = add_info.get_halflife().get_value()
    except:
        return ''
    else:
        return hl.split(';')[0]

def arcessere_halflife_comment(add_info):
    try:
        hl = add_info.get_halflife().get_value()
    except:
        return ''
    else:
        return hl.split(';')[2]

def arcessere_halflife_source(add_info):
    try:
        hl = add_info.get_halflife().get_value()
    except:
        return ''
    else:
        return hl.split(';')[4]

def arcessere_temperature(add_info):
    try:
        temp = add_info.get_temperature().get_value()
    except:
        return np.NaN
    else:
        min = float(temp.split(';')[0])
        max = float(temp.split(';')[1])
        return np.round(np.average([min, max]), 0)

def arcessere_temperature_unit(add_info):
    try:
        temp_unit = add_info.get_temperature().get_unit()
    except:
        return ''
    else:
        return '\u2103'


if __name__ == "__main__":

    try:
        eP = enviPath(INSTANCE_HOST)
        eP.login(USERNAME, PASSWORD)  # getpass.getpass()
        print(eP.who_am_i().get_name())

        package_id = ORIGINAL_SLUDGE
        package = Package(eP.requester, id=package_id)
        all_scenarios = package.get_scenarios()
        all_pathways = package.get_pathways()
        print('Number of scenarios found:', len(all_scenarios))
        print('Number of pathways found:', len(all_pathways))

        main()

    finally:
        eP.logout()