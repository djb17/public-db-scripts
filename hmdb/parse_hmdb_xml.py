"""
Dongjo Ban
Short example on how to traverse through HMDB XML file to extract desired information using
chemical formulas of the small compounds
02/03/2023
"""
import pickle
import xml.etree.ElementTree as ET

def retrieve_via_formula(metabolites, formula, current_dict):
    
    if formula in list(current_dict.keys()):
        return(current_dict)
    
    info_dict = {'name': [], 'synonyms':[], 'accession': [], 'foodb_id':[], 'drugbank_id':[],
                 'average_molecular_weight':[], 'description': [], 'diseases':[],
                 'origin':[], 'pathways': [], 'kingdom':[], 'super_class':[], 'met_class':[],
                 'sub_class':[], 'molecular_framework':[]}
    for metabolite in metabolites:

        # get the chem formula of the current metabolite
        chemical_formula = metabolite.find("{http://www.hmdb.ca}chemical_formula")
        
        # if the formula exists, then check if its the same as the one we're querying
        if chemical_formula is not None and chemical_formula.text == formula:
            name_element = metabolite.find("{http://www.hmdb.ca}name")
            synm_element = metabolite.find("{http://www.hmdb.ca}synonyms")
            
            desc_element = metabolite.find("{http://www.hmdb.ca}description")
            acsn_element = metabolite.find("{http://www.hmdb.ca}accession")
            molw_element = metabolite.find("{http://www.hmdb.ca}average_molecular_weight")
            
            food_element = metabolite.find("{http://www.hmdb.ca}foodb_id")
            drug_element = metabolite.find("{http://www.hmdb.ca}drugbank_id")
            sick_element = metabolite.find("{http://www.hmdb.ca}diseases")
            
            taxonomy = metabolite.find("{http://www.hmdb.ca}taxonomy")
            kingdom = taxonomy.find("{http://www.hmdb.ca}kingdom")
            
            super_class = taxonomy.find("{http://www.hmdb.ca}super_class")
            met_class = taxonomy.find("{http://www.hmdb.ca}class")
            sub_class = taxonomy.find("{http://www.hmdb.ca}sub_class")
            molecular_framework = taxonomy.find("{http://www.hmdb.ca}molecular_framework")

            info_dict['name'].append('' if name_element is None else name_element.text)
            info_dict['description'].append('' if desc_element is None else desc_element.text)
            info_dict['accession'].append('' if acsn_element is None else acsn_element.text)
            info_dict['average_molecular_weight'].append('' if molw_element is None else molw_element.text)

            info_dict['kingdom'].append('' if kingdom is None else kingdom.text)
            info_dict['super_class'].append('' if super_class is None else super_class.text)
            info_dict['met_class'].append('' if met_class is None else met_class.text)
            info_dict['sub_class'].append('' if sub_class is None else sub_class.text)
            info_dict['molecular_framework'].append('' if molecular_framework is None else molecular_framework.text)
            
            info_dict['foodb_id'].append('' if food_element is None else food_element.text)
            info_dict['drugbank_id'].append('' if drug_element is None else drug_element.text)
            
            
            info_dict['diseases'].append(';'.join([i.find("{http://www.hmdb.ca}name").text for i in sick_element.findall('*')]))
            info_dict['synonyms'].append([i.text for i in synm_element.findall('*')])
                      
            # look for origin of the compound under 'Ontology'
            onto_element = metabolite.find("{http://www.hmdb.ca}ontology")
            if len(onto_element) > 0:
                
                onto_desc = onto_element.find(".//{http://www.hmdb.ca}descendant[{http://www.hmdb.ca}definition='Natural or synthetic origin of a chemical.']")
                if onto_desc is not None:
                    onto_descs = onto_desc.find('./{http://www.hmdb.ca}descendants')
                    onto_desc_2 = onto_descs.find('./{http://www.hmdb.ca}descendant')
                    origin = onto_desc_2.find('{http://www.hmdb.ca}term').text
                    info_dict['origin'].append(origin)
                else:
                    info_dict['origin'].append('')
            else:
                info_dict['origin'].append('')
            
            # Look for any pathway information under 'Biological Properties'
            biop_element = metabolite.find("{http://www.hmdb.ca}biological_properties")
            if len(biop_element) > 0:
                biop_paths = biop_element.find("{http://www.hmdb.ca}pathways")
                info_dict['pathways'].append(';'.join([i.find("{http://www.hmdb.ca}name").text for i in biop_paths.findall("*")]))

    current_dict[formula] = info_dict
    return(current_dict)

def retrieve_wrapper(metabolites, formula_dict):
    hmdb_dict = dict()
    for i,f in enumerate(list(formula_dict.values())):
        if i % 500 == 0:
            print(f"Processing {i}th formula")
        retrieve_via_formula(metabolites, f, hmdb_dict)
    return(hmdb_dict)

def main():
    tree = ET.parse('../data/hmdb/nov2021/hmdb_metabolites.xml')
    root = tree.getroot()

    # .// <- find all metabolite elements
    metabolites = root.findall("./{http://www.hmdb.ca}metabolite")
    print(len(metabolites)) # total # of metabolites

    # preliminary look at available tags for the small compounds found in XML file
    for m in metabolites[0]:
        print(m.tag)

    # below is what the input dictionary will need to look like to run the script
    sample_feat_dict = {'RN2': 'C3H6O2', 'RN5': 'C3H4O3',
                        'RN9': 'C6H8O6', 'RN12': 'C3H6O3',
                        'RN19': 'C6H12O6', 'RN20': 'C2H7N2P'}

    sample_hmdb = retrieve_wrapper(metabolites, sample_feat_dict)
    pickle.dump(sample_hmdb, open('./sample_hmdb_byformula.pickle', 'wb'))

if __name__ == “__main__”:
         main()