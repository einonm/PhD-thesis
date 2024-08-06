# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

#
# # Display Results
#
# Some library functions for transforming and displaying results
#

# +
treatment_antihypertensives = {
    "A02AC02",
    "A10BH52",
    "B05XA04",
    "C01DA02",
    "C01DA04",
    "C01DA08",
    "C01DA14",
    "C01DA58",
    "C02AA02",
    "C02AA52",
    "C02AB01",
    "C02AC01",
    "C02AC02",
    "C02AC05",
    "C02AC06",
    "C02BB01",
    "C02CA01",
    "C02CA02",
    "C02CA04",
    "C02CA06",
    "C02CC02",
    "C02DB02",
    "C02DD01",
    "C02KX01",
    "C02LA01",
    "C02LA51",
    "C02LA71",
    "C02LC05",
    "C02LE01",
    "C03AA01",
    "C03AA02",
    "C03AA03",
    "C03AA04",
    "C03AA05",
    "C03AA08",
    "C03AB01",
    "C03AB02",
    "C03AB03",
    "C03AB04",
    "C03AB05",
    "C03AB08",
    "C03AH01",
    "C03AH02",
    "C03AX01",
    "C03BA03",
    "C03BA04",
    "C03BA08",
    "C03BA10",
    "C03BA11",
    "C03BB03",
    "C03CA01",
    "C03CA03",
    "C03CA04",
    "C03CB01",
    "C03CC01",
    "C03DA01",
    "C03DA02",
    "C03DA04",
    "C03DA04",
    "C03DB01",
    "C03DB02",
    "C03EA03",
    "C03EA13",
    "C03EB01",
    "C04AB01",
    "C04AB02",
    "C04AX02",
    "C05AE02",
    "C05AE03",
    "C06BA06",
    "C07AA02",
    "C07AA03",
    "C07AA05",
    "C07AA06",
    "C07AA12",
    "C07AA15",
    "C07AA23",
    "C07AB02",
    "C07AB03",
    "C07AB04",
    "C07AB05",
    "C07AB07",
    "C07AB08",
    "C07AB12",
    "C07AG01",
    "C07AG02",
    "C07BA05",
    "C07BA12",
    "C07BB02",
    "C07BB07",
    "C07BB12",
    "C07BB52",
    "C07BG01",
    "C07CA03",
    "C07CB02",
    "C07CB03",
    "C07CB53",
    "C07CG01",
    "C07DB01",
    "C07FB02",
    "C07FB03",
    "C07FB07",
    "C07FB12",
    "C07FB13",
    "C07FX01",
    "C07FX03",
    "C07FX04",
    "C07FX05",
    "C07FX06",
    "C08CA01",
    "C08CA02",
    "C08CA03",
    "C08CA04",
    "C08CA05",
    "C08CA06",
    "C08CA07",
    "C08CA08",
    "C08CA09",
    "C08CA13",
    "C08CA14",
    "C08CA15",
    "C08CA16",
    "C08CA51",
    "C08CA55",
    "C08DA01",
    "C08DA51",
    "C08DB01",
    "C08GA02",
    "C08GA02",
    "C09AA01",
    "C09AA02",
    "C09AA03",
    "C09AA04",
    "C09AA05",
    "C09AA06",
    "C09AA07",
    "C09AA08",
    "C09AA09",
    "C09AA10",
    "C09AA12",
    "C09AA13",
    "C09AA15",
    "C09BA01",
    "C09BA02",
    "C09BA03",
    "C09BA04",
    "C09BA05",
    "C09BA07",
    "C09BA08",
    "C09BA09",
    "C09BA13",
    "C09BA15",
    "C09BB02",
    "C09BB03",
    "C09BB04",
    "C09BB05",
    "C09BB06",
    "C09BB07",
    "C09BB10",
    "C09BB12",
    "C09BX01",
    "C09BX02",
    "C09BX03",
    "C09BX04",
    "C09BX11",
    "C09CA01",
    "C09CA02",
    "C09CA03",
    "C09CA04",
    "C09CA06",
    "C09CA07",
    "C09CA08",
    "C09CA09",
    "C09CA10",
    "C09DA01",
    "C09DA02",
    "C09DA03",
    "C09DA04",
    "C09DA07",
    "C09DA09",
    "C09DB01",
    "C09DB02",
    "C09DB04",
    "C09DB05",
    "C09DB06",
    "C09DB07",
    "C09DB08",
    "C09DB09",
    "C09DX01",
    "C09DX02",
    "C09DX03",
    "C09DX04",
    "C09DX05",
    "C09DX06",
    "C09DX07",
    "C09XA02",
    "C09XA52",
    "C09XA53",
    "C09XA54",
    "C10AA05",
    "C10AA07",
    "C10BA05",
    "C10BA06",
    "C10BX03",
    "C10BX04",
    "C10BX05",
    "C10BX06",
    "C10BX07",
    "C10BX08",
    "C10BX09",
    "C10BX10",
    "C10BX11",
    "C10BX12",
    "C10BX13",
    "C10BX14",
    "C10BX15",
    "C10BX16",
    "C10BX17",
    "G01AE10",
    "G04CA03",
}

# Using drugbank.ca

# +
drug_hypertensives = {  # https://go.drugbank.com/categories/DBCAT004567
    "R03AC08",
    "R07AB01",
    "C01BD07",
    "R03AC14",
    "N06AG02",
    "N02CC03",
    "R03CC63",
    "R03AC02",
    "R03AK13",
    "R03AK08",
}

drug_bp_side_effects = {
    "C02KX05",
    "R07AX01",
    "C07AA01",
    "A06AD04",
    "C08CA10",
    "C07AA17",
    "C08CA11",
    "D03AX10",
    "P01BC01",
    "A14AB01",
    "A04AA01",
    "N07AB01",
    "L01XE10",
    "N03AB02",
    "L01XX23",
    "L04AA10",
    "L01XE46",
    "M05BB06",
    "J01MA02",
}

# C02KX05 Riociguat - Pulmonary hypertension treatment
# L01XE35 Osimertinib - anticancer drug
# L01XX18 Tiazofurine - anticancer drug
# C01EB18 Ranolazine - stable angina treatment]
# R07AX01 Nitric oxide - imapirment of NO is arisk factor in hypertension
# C07AA01 Alprenolol - used to treat hypertension, but no longer marketed by AZ
# A06AD04 Magnesium sulfate - antihypertensive (user in ER?)
# C08CA10 Nilvadipine - used to treat arterial hypertension
# R06AX26 Fexofenadine - antihistamine (nothing related to BP found?)
# N07XX02 Riluzole - Used to treat ALS
# A01AB09 Miconazole - antifungal
# C07AA17 Bopindolol - used for the management of hypertension
# C08CA11 Manidipine - used for mild to moderate hypertension
# D03AX10 Enoxolone - also a flavouring (?). From liquorice, not recommended for those with hypertension.
# P01BC01 Quinine - decreased blood pressure is a side effect
# A14AB01 Nandrolone - hypertension listed as a side effect (NICE)
# A04AA01 Ondansetron - hypertension listed as an uncommon side effect
# N07AB01 Carbachol - some evidence that it raises blood pressure
# L01XE10 Everolimus - can cause hypertension
# N03AB02 Phenytoin - decreases blood pressure
# L01XX23 Mitotane - hypertension side effect
# M03BB03 Chlorzoxazone - muscle relaxant
# L04AA10 Sirolimus - Hypertension a common side effect
# L01XE46 Encorafenib - Hypertension a common side effect
# M05BB06 Alendronic acid and alfacalcidol, sequential - Hypertension a uncommon side effect of alfacalcidol
# J01MA02 Ciprofloxacin - less common side effects are hyper/hypotension
# R06AX11 Astemizole - antihistamine, some heart rhythm side effects but not BP (?)
# -

treatment_antihypercholesterol = {  # https://go.drugbank.com/drugs/DB04377
    "B01AC06",
    "C07FX04",
    "C10BX04",
    "M01BA03",
    "C10BX02",
    "B01AC56",
    "N02AJ07",
    "N02AJ02",
    "N02BA01",
    "C10BX05",
    "N02BA51",
    "C10BX01",
    "C07FX03",
    "N02AJ18",
    "C10BX12",
    "C10BX06",
    "C07FX02",
    "C10AX14",
    "C10BX15",
    "C10AA05",
    "C10BX03",
    "C10BA05",
    "C10BX11",
    "C10BX08",
    "C10BX06",
    "C10BX12",
    "C10AX15",
    "C10AC04",
    "C10AX13",
    "C10AA04",
    "C10AA02",
    "C10BA01",
    "C10AA03",
    "C10BA03",
    "C10BX02",
    "C10BX10",
    "C10BX05",
    "C10BX07",
    "C10AA07",
    "A10BH52",
    "C10BX13",
    "C10BX09",
    "C10BA09",
    "C10BX16",
    "C10BX17",
    "C10BX14",
    "C10BA06",
    "A10BH51",
    "C10AA01",
    "C10BA02",
    "C10BA04",
    "C10AA06",
    "C10AB01",
    "C10AX09",
    "C10AC01",
    "C10AX02",
    "C10AX11",
    "C10AA08",
    "C10AB09",
    "C10AB11",
    "C10AX05",
}

# TODO
drug_cholesterol = {
    "C10AB02",
    "C10AB04",
}
drug_cl_side_effects = {
    "A10BD17",
    "G03AA09",
    "G03CA01",
    "G03AA04",
    "G03AA01",
    "G03AA05",
    "G03AA13",
    "G03AA12",
    "A14AA05",
    "A02BA03",
    "D10AB01",
    "G03AA07",
    "G03AA06",
    "G03AA14",
    "G03AA08",
}

# ATC_CODE - drugbank indication
# ----------------------------------
indications = {
    'A01AB02':'ear / skin infections',
    'A01AB23':'infections',
    'A02AA01':'heartburn / indegestion',
    'A02BC01':'GERD',
    'A02BC02':'GERD',
    'A02BC03':'gastrointestinal ulcers / GERD',
    'A02BD01':'GERD / bacterial infections',
    'A02BD02':'gastrointestinal ulcers/GERD/bacterial infections',
    'A02BD11':'GERD / bacterial infections',
    'A03AA05':'irritable bowel syndrome',
    'A03AB12':'peptic ulcer',
    'A03BA03':'IBS/heart block/parkinsonism/peptic ulcer/colic/rhinitis',
    'A04AD12':'chemotherapy-induced nausea and vomiting',
    'A06AD04':'constipation/convulsions/hypomagnesemia',      
    'A08AA02':'seizures', 
    'A08AA05':'Duchenne muscular dystrophy / appetite supressor', 
    'A09AB02':'homocystinuria', 
    'A10BB03':'hyperglycemia (type II diabetes)',
    'A11AA01':'vitamin deficiency',
    'A11CC01':'vitamin deficiency / anemia',
    'A11HA07':'nutritional ingredient (under investigation to treat depression/anxiety/psychiatric disorder)',
    'B01AC04':'antiplatelet / vascular & heart disease',
    'B01AC22':'vascular & heart disease',
    'B01AC23':'intermittent claudication',
    'B03AE01':'vitamin deficiency / anemia',
    'B03AE02':'vitamin deficiency / anemia',
    'C01BA03':'ventricular arrhythmia', 
    'C01BD01':'atrial fibrillation / tachycardia',
    'C01BD04':'atrial fibrillation',
    'C01BD07':'atrial fibrillation',
    'C01CA14':'heart failure (wikipedia)',
    'C01CA17':'hypotension',
    'C01EB10':'tachycardia',
    'C01EB18':'angina / arrhythmia',
    'C02AA02':'hypertension',
    'C02BB01':'hypertension',
    'C02CA04':'hypertension',
    'C03DA01':'acne/heart failure/hypertension',
    'C04AB01':'hypertensive crisis / anasthetic reversal',
    'C04AC01':'high cholesterol',
    'C07AA05':'tremor caused by lithium/hypertension/migraine/anxiety',
    'C07AA07':'ventricular arrhytmias',
    'C07AG01':'hypertension',
    'C07FB02':'hypertension',
    'C07FB03':'hypertension',
    'C07FB07':'hypertension',
    'C07FB12':'hypertension',
    'C07FB13':'hypertension / heart failure',
    'C07FX05':'heart failure',
    'C07FX06':'heart failure',
    'C08CA02':'hypertension',                                                                                                                                                           
    'C08CA03':'hypertension',                                                                                                                                                           
    'C08CA04':'hypertension/angina/migraine',
    'C08CA05':'hypertension / angina',
    'C08CA06':'delayed ischemic neurological deficit',                                                                                                                                  
    'C08CA07':'hypertension',                                                                                                                                                           
    'C08CA08':'hypertension',
    'C08CA09':'hypertension',
    'C08CA10':'hypertension arterial',                                                                                                                                                  
    'C08CA11':'hypertension',                                                                                                                                                           
    'C08CA16':'hypertension',                                                                                                                                                           
    'C08CX01':'hypertension',                                                                                                                                                           
    'C08DA01':'hypertension/angina/cluster headache/tachycardia',                                                                                                                       
    'C09AA01':'hypertension / heart failure',
    'C09AA02':'hypertension',
    'C09AA04':'hypertension / heart failure',   
    'C09BB02':'hypertension',     
    'C09BB03':'hypertension',                                                                                                                                                           
    'C09BB04':'hypertension / heart failure',                                                                                                                                           
    'C09BB05':'hypertension',                                                                                                                                                           
    'C09BB06':'hypertension',                                                                                                                                                           
    'C09BB07':'hypertension',                                                                                                                                                           
    'C09BB10':'hypertension',                                                                                                                                                           
    'C09BB12':'hypertension / heart failure',                                                                                                                                           
    'C09BX01':'hypertension',  
    'C09BX02':'hypertension / heart failure',  
    'C09BX03':'hypertension',
    'C09DB01':'hypertension',
    'C09DB02':'hypertension',                                                                                                                                                           
    'C09DB04':'hypertension',           
    'C09DB05':'hypertension',           
    'C09DB06':'hypertension',                                                                                                                                                           
    'C09DB07':'hypertension',
    'C09DX01':'hypertension / heart failure',
    'C09DX03':'hypertension',                                                                                                                                                           
    'C09DX06':'hypertension',                                                                                                                                                           
    'C09XA53':'hypertension',
    'C09XA54':'hypertension',
    'C10AA02':'high cholesterol',
    'C10AX02':'high cholesterol',
    'C10BX03':'high cholesterol / heart failure',
    'C10BX07':'high cholesterol / heart failure',
    'C10BX09':'high cholesterol / heart failure',
    'C10BX11':'high cholesterol / heart failure',
    'C10BX14':'high cholesterol / heart failure',
    'D01AC03':'skin allergies',
    'D01AC08':'fungal infection',
    'D01AC11':'fungal infection',
    'D01AC15':'fungal infection',
    'D04AA10':'allergies / anaphylaxis',
    'D04AX01':'hives',
    'D08AX08':'infections',
    'D10AD03':'acne',
    'G02AB02':'hemorrhage / angina',
    'G02AB03':'hemorrhage / angina',
    'G02CB02':'Parkinson\'s disease',
    'G02CB05':'amenorrhea/diarrhea/migraine',
    'G02CB06':'hyperprolactinemia (wikipedia)',
    'G02CX02':'hypoactive sexual desire disorder',
    'G03XC01':'breast cancer / osteoporosis',
    'G04BD10':'overactive bladder syndrome',
    'G04BE07':'advanced parkinson\'s disease / decreased mobility',
    'G04CA02':'enlarged prostate',
    'G04CA52':'enlarged prostate',
    'G04CA53':'enlarged prostate',
    'G04BD06':'urinary issues',
    'J01CA04':'infection (antibiotic)', 
    'J01DD01':'infection (antibiotic)', 
    'J01DD04':'infection (antibiotic)', 
    'J01DF01':'infection (antibiotic)', 
    'J01EA01':'urinary/respiratory/gastrointestinal tract infection (antibiotic)', 
    'J01EC02':'encephalitis / infection (antibiotic)', 
    'J01ED08':'infection (antibiotic)', 
    'J01FA09':'infection (antibiotic)',
    'J01MA01':'infection (antibiotic)',
    'J01MA02':'infection (antibiotic)',
    'J01MA03':'infection (antibiotic)',
    'J01MA06':'infection (antibiotic)',
    'J01MA09':'infection (antibiotic)',
    'J01MB06':'infection (antibiotic)',
    'J01RA07':'infection (antibiotic)',
    'J01RA10':'infection (antibiotic)',
    'J02AC01':'infection / tubercolosis',
    'J02AC03':'infection (antifungal)', 
    'J04AB02':'infection (antibiotic)',
    'J04AC01':'infection (antibiotic)',
    'J04AM02':'infection (antibiotic)',
    'J04AM05':'infection (antibiotic)',
    'J04AM06':'infection (antibiotic)',
    'J04AM07':'infection (antibiotic)',
    'J04AM08':'infection (antibiotic)', 
    'J04AK01':'infection / tubercolosis',
    'J04AK02':'infection / tubercolosis',
    'J05AE02':'HIV',
    'J05AF01':'HIV',
    'J05AF02':'HIV',
    'J05AF04':'HIV',
    'L01AD01':'Brain tumors/myeloma/lymphoma ',
    'L01BB04':'leukemia/lymphomas/multiple sclerosis',
    'L01CB02':'leukemia/lymphomas/multiple sclerosis',
    'L01CB02':'leukemia',
    'L01CD01':'cancers',
    'L01CX01':'cancers',
    'L01DA01':'cancers (RNA synthesis inhibitor)',
    'L01DB02':'leukemia (DNA replication inhibitor)',
    'L01DB03':'breast/lung/bladder cancer (DNA/RNA inhibitor)',
    'L01DB07':'leukemia / multiple sclerosis (DNA repair inhibitor)',
    'L01XA01':'ovarian/bladder/testicular cancer',
    'L01XE03':'pancreatic / lung cancer',
    'L01XE15':'melanoma / lung cancer',
    'L01XE27':'leukemia / lymphoma',
    'L01XX54':'ovarian cancer',
    'L01XY01':'leukemia',
    'L04AA29':'colitis / arthritis',
    'M01AB02':'osteoarthritis/rheumatoid arthritis/ankylosing spondylitis',
    'M01AC01':'osteoarthritis / rheumatoid arthritis',
    'M01AG01':'osteoarthritis / rheumatoid arthritis',
    'M03AC04':'neuromuscular blocker/muscle relaxer under mechanically ventialted anesthesia',
    'M03BX01':'malignant hyperthermia / spasticity',
    'M03CA01':'alcohol dependency / spasticity',
    'M05BA01':'osteoporosis / ossification',
    'N01AB01':'anesthesia',
    'N01AB04':'anesthesia',
    'N01AB06':'anesthesia',
    'N01AB07':'anesthesia',
    'N01AB08':'anesthesia',
    'N01AF02':'anesthesia',
    'N01AF03':'anesthesia / seizures',
    'N01AX03':'anesthesia',
    'N01AX11':'cataplexy / excessive sleepiness',
    'N01AX14':'depression/MDD/pain',
    'N01BB03':'anesthesia',
    'N01BC01':'anesthesia',
    'N02AF01':'pain relief',  
    'N02AB02':'pain relief',  
    'N02AF02':'pain relief',  
    'N02BG06':'pain relief',  
    'N02BG09':'pain relief',  
    'N02CA01':'cluster headache / migraine',  
    'N02CA02':'cluster headache / migraine',
    'N02CA04':'cluster headache / migraine',
    'N02CC01':'migraine',
    'N02CC02':'migraine',
    'N02CC06':'migraine',
    'N02CX01':'migraine',
    'N03AA02':'seizures/anxiety/insomnia',  
    'N03AA30':'epilepsy', 
    'N03AX11':'seizures / migraine',
    'N03AX12':'seizures / neuropathic pain', 
    'N03AX14':'epilepsies / seizures', 
    'N03AX16':'epilepsies/anxiety/neuropathic pain/seizures',
    'N04AA01':'movement disorder / parkinsonism',
    'N04AA02':'movement disorder / parkinsonism',
    'N04AA04':'movement disorder / parkinsonism',
    'N04BA04':'Parkinson\'s disease',
    'N04BC08':'Parkinson\'s disease',
    'N04BC09':'Parkinson\'s disease / restless leg syndrome',
    'N05AA01':'schizophrenia / hyperactivity',                                                                                                                                            
    'N05AA02':'schizophrenia/bipolar disorder/anxiety/neuralgia',
    'N05AA03':'psychosis / schizophrenia',
    'N05AA05':'psychosis',
    'N05AB03':'schizophrenia/anxiety/depression',                                                                                                                                       
    'N05AC04':'schizophrenia',                                                            
    'N05AD01':'schizophrenia/tourettes/huntingtons/OCD/hyperactivity',
    'N05AD08':'agitation/delerium/nausea',
    'N05AF03':'schizophrenia', 
    'N05AF04':'schizophrenia',                                                                                                                                                          
    'N05AG02':'delusional parasitosis / tourettes tics',  
    'N05AH01':'schizophrenia',
    'N05AH03':'schizophrenia/depression/bipolar/MDD/PTSD',
    'N05AL05':'schizophrenia',
    'N05AX12':'schizophrenia/MDD/bipolar 1/tourettes/agitation/manic depression',                                                                                 
    'N05AX13':'schizoaffective disorders/schizophrenia',                                                                                                                                
    'N05AX14':'schizophrenia',                                                                                                                                                          
    'N05AX15':'schizophrenia/acute depression/bipolar disorder',                                                                                                                        
    'N05AX16':'schizophrenia/MDD', 
    'N05BA10':'anxiety',
    'N05BB01':'anxiety',   
    'N05BC01':'anxiety',                                                                                                                                                                
    'N05CA01':'convulsions / insomnia',                                                                                                                                                 
    'N05CA04':'sedative / hypnotic',                                                                                                                                                      
    'N05CA06':'insomnia',                                                                                                                                                               
    'N05CB01':'insomnia',                                                                                                                                                               
    'N05CD03':'anxiety / insomnia',     
    'N05CF02':'insomnia',     
    'N05CM06':'insomnia',                                                                                                                                                               
    'N05CM08':'insomnia',
    'N06AB03':'bulimia/anxiety/MDD/OCD',
    'N06AB06':'bulimia/anxiety/MDD/OCD/PTSD',
    'N06AB08':'bulimia/MDD/OCD',
    'N06AG02':'MDD',
    'N06BX13':'optic neuropathy/Alzheimer\'s disease',
    'N06AA01':'depression/insomnia/anorexia/bulimia/IBS/pain',
    'N06AA04':'depression / OCD',
    'N06AA12':'anxiety/bipolar disorder/depression/insomnia',
    'N06AA17':'anxiety/depression',
    'N06AA21':'anxiety/bipolar disorder/depression',
    'N06AB03':'anorexia/bulimia/depression/MDD/OCD',
    'N06AB06':'bulimia/depression/MDD/OCD/PTSD',
    'N06AF04':'MDD / depression',
    'N06AX04':'(withdrawn) depression', 
    'N06AX05':'Alzheimer\'s/anxiety/dementia/insomnia/MDD/schizophrenia', 
    'N06AX07':'depression',
    'N06AX11':'depression/anxiety/insomnia/MDD/OCS/pain',
    'N06AX21':'depression/pain/anxiety',
    'N06AX26':'MDD',
    'N06BA01':'ADHD/depression/narcolepsy/obesity/pain',
    'N06BA03':'ADHD / narcolepsy',
    'N06BX13':'Alzheimer\'s disease / optic neuropathy',
    'N06DA02':'Alzheimer\'s disease / dementia',
    'N06DA52':'Alzheimer\'s disease / dementia',
    'N06DX01':'Alzheimer\'s disease / dementia',
    'N07BB03':'alcohol dependency',
    'N07CA02':'balance disorders/vertigo/dizziness/inadequate cerebral circulation',
    'N07CA03':'migraine',   
    'N07XX02':'amyotrophic lateral sclerosis',
    'N07XX06':'tourette\'s/huntington\'s disease/tardive dyskinesia',
    'N07XX07':'multiple sclerosis',
    'P03AC02':'lice and scabies infestation',
    'R01AD04':'asthma',
    'R03AL04':'COPD',
    'R03AL05':'COPD',
    'R03BB05':'COPD',
    'R05CB01':'rhinitis',
    'R06AX02':'allergies / anaphylaxis',     
    'R06AX17':'allergic rhinitis',
    'R06AX22':'allergic rhinitis',
    'R07AB01':'COPD / respiratory depression',
    'S01EX02':'Iatrogenically induced mydriasis',
    'V03AB05':'endocrine/rheumatic/hematologic disorders',
    'V03AB32':'nerve disorders/neuropathy/hangover',
    'V03AF02':'doxorubicin induced cardiomyopathy / drug extravasation',
}

# Drugs with no studies investigating their effect on schizophrenia
sz_novel_drugs = {
    "A06AD04",
    "N07CA02",
    "C09BB04",
    "C08CA07",
    "C08CA11",
    "C09BB12",
    "C01BD07",
    "C08CA08",
    "C08CA16",
    "C09BB05",
    "C09XA53",
    "C08CA02",
    "C09DB07",
    "L01XA01",
    "P03AC02",
    "C09BB03",
    "C09BB07",
    "C09DB06",
    "C09DX06",
    "C08CA04",
    "C08CX01",
    "C09BX01",
    "C09BX03",
    "C09DB02",
    "G04BD06",
    "N06AX07",
    "A03AA05",
    "A11CC01",
    "C01EB18",
    "C08CA09",
    "N02CA01",
    "C01BD01",
    "C07FB13",
    "C07FX06",
    "C09DB01",
    "C09DX01",
    "C09XA54",
    "D01AC03",
    "G02AB03",
    "G04CA02",
    "G04CA52",
    "G04CA53",
    "L01XX54",
    "N04BA04",
    "R03AL04",
    "R03BB05",
    "R06AX22",
    "A01AB02",
    "A03AB12",
    "A03BA03",
    "C01CA14",
    "C04AB01",
    "C07FB02",
    "C07FB07",
    "C07FB12",
    "C07FX05",
    "C09AA01",
    "G04BD10",
    "N02CX01",
    "N04BC08",
    "N05CM06",
    "R03AL05",
    "S01EX02",
}

mdd_novel_drugs = {
    "N05AA05",
    "N05BA10",
    "C01BA03",
    "L01XA01",
    "N02BG06",
    "N01AB07",
    "N02CC02",
    "N02CC06",
    "A10BB03",
    "B01AC04",
    "C08CA07",
    "G02AB03",
    "G03XC01",
    "J01MA01",
    "J04AK02",
    "M05BA01",
    "N01AB04",
    "N02CA01",
    "N05AX13",
    "N06BX13",
}

bp_novel_drugs = {
    "A06AD04",
    "C08CA10",
    "N05CM08",
    "N05CM08",
    "C08CA07",
    "C08CA03",
    "C08CA08",
    "C08CA16",
    "N05CA01",
    "N05CA06",
    "N05CA06",
    "C09BB12",
    "C09BB12",
    "N01AF02",
    "C08CX01",
    "N05CD03",
    "N06AX07",
    "C08CA09",
    "C09BB06",
    "G02CB05",
    "N01AB01",
    "N01AF03",
    "N03AA30",
    "N05CA04",
    "N05CB01",
    "C01BD07",
    "D10AD03",
    "L04AA29",
    "N01AB04",
    "N01AB08",
    "N02CA04",
    "P03AC02",
    "A03AA05",
    "A09AB02",
    "C09AA02",
    "C09BB02",
    "C09BB05",
    "G02AB02",
    "G02AB03",
    "G02CB02",
    "L01CD01",
    "L01CX01",
    "L01XE15",
    "L01XE27",
    "L01XX54",
    "N01AB06",
    "N01AB07",
    "N02AB02",
    "N02CX01",
    "N05AC04",
    "N07XX07",
    "R06AX17",
}

al_novel_drugs = {
    "M03AC04",
    "A02AA01",
}

pd_novel_drugs = {
    "A02BC02",
    "A02BD01",
    "A02BC03",
    "B01AC22",
    "D01AC11",
    "D01AC15",
    "G02CX02",
    "J01CA04",
    "J01DF01",
    "J01ED08",
    "J01RA10",
    "J02AC03",
    "J04AK01",
    "J05AE02",
    "J05AF04",
    "L01AD01",
    "M01AG01",
    "N01BB03",
    "N02AF01",
    "N05BB01",
    "A02BC01",
    "A02BD11",
    "B01AC04",
    "B01AC23",
    "C01BD04",
    "C07AA07",
    "J01MA02",
    "J01RA07",
    "J05AF02",
    "L01XE03",
    "R01AD04",
}

hd_novel_drugs = {
    "L01DB02",
    "A08AA05",
    "A11AA01",
    "J01MA02",
    "J01MA09",
    "J01RA10",
    "L01CB02",
    "L01DB03",
    "L01DB07",
    "V03AF02",
    "C07AG01",
    "J01MA01",
    "J01MA03",
    "J01MA06",
    "J01MB06",
    "L01BB04",
    "R03AL05",
    "R03BB05",
}

# +
import os
import sys
import pandas as pd
import numpy as np

# #!pip install scipy
from scipy import stats

# #!pip install matplotlib
# %matplotlib inline
import matplotlib.pyplot as plt
import seaborn as sns

# #!pip install matplotlib_venn
from matplotlib_venn import venn2

import qvalue
from IPython.display import display, HTML, Latex, display_latex

# +
# #!pip install rpy2
# import rpy2

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects import numpy2ri

from rpy2.robjects.conversion import localconverter

rqvalue = importr("qvalue")
numpy2ri.activate()


# +
# Captions for diagrams and figures
def start_caption(formatting, bold_caption, normal_caption, newpage=False):
    if formatting == "HTML":
        display("Figure: " + bold_caption + " " + normal_caption)
    if formatting == "LaTeX":
        if newpage:
            display(Latex(r"\newpage"))
        display(Latex(r"\begin{figure}[htbp]"))
        display(Latex(r"\centering"))
        display(
            Latex(
                r"\caption{" + bold_caption + r" \normalfont{" + normal_caption + "}}"
            )
        )


def end_caption(formatting, label):
    if formatting == "LaTeX":
        display(Latex(r"\label{" + label + "}"))
        display(Latex(r"\end{figure}"))


# +
# Summary table functions
# Returns ATC data (ATC code & description) in a dataframe, given an ATC table file of the same
def get_atc_file(filename):
    result = pd.DataFrame()
    with open(filename, "r") as file:
        lines = file.readlines()
        for line in lines:
            code = line.split(" ")[0]
            desc = " ".join(line.split(" ")[1:]).strip()
            result = result.append(
                {"atc_code": code, "Description": desc}, ignore_index=True
            )

    return result


# Returns ATC data from all 5 ATC levels in one dataframe
def get_atc_levels(data_path):
    atc_level_files = [
        data_path + "atc/level-" + str(i) + "-atc.txt" for i in range(1, 6)
    ]

    # convert each atc_level_file to a 2-column data frame, and cat all of them together
    atc_levels = pd.concat(
        [get_atc_file(atc_level_file) for atc_level_file in atc_level_files]
    )
    atc_levels.set_index("atc_code", inplace=True)
    return atc_levels


# Read all results from a magma geneset analysis for this gwas and annotation type
def read_results_files(
    gwas,
    annot_type,
    magma_results_dir,
    run_id,
    compsets,
    magma_ver=106,
):
    prefix = os.path.join(magma_results_dir, gwas, run_id)
    if magma_ver >= 108:
        postfix = "-" + compsets + ".gsa.out"
    else:
        postfix = "-" + compsets + ".sets.out"

    file = os.path.join(
        prefix,
        "magma_geneset_result-" + gwas + "-" + annot_type + "-" + run_id + postfix,
    )

    df = pd.read_csv(file, comment="#", delim_whitespace=True)
    if magma_ver >= 108:
        df.drop(["TYPE"], axis=1, inplace=True)

    # Add Q-values
    df.loc[:, "Q " + annot_type + " " + compsets] = qvalue.estimate(np.array(df["P"]))

    df.drop(["BETA", "BETA_STD", "SE"], axis=1, inplace=True)

    if magma_ver >= 108:
        rename_column = "VARIABLE"
    else:
        rename_column = "SET"
    df.rename(
        columns={rename_column: "ATC_CODE", "P": "P " + annot_type + " " + compsets},
        inplace=True,
    )

    return df


# Summarise results by applying significance and number of gene thresholds, and store in a file
# Generates results for both P and Q values.
def summarise_drug_results(
    gwas,
    annot_type,
    atc_levels,
    magma_results_dir,
    summary_results_dir,
    compsets,
    run_id,
    magma_ver=106,
    n_gene_thresh=4,
    signif_thresh=0.05,
):

    result = read_results_files(
        gwas, annot_type, magma_results_dir, run_id, compsets, magma_ver
    )

    # only consider classes/drugs with Qval < QVAL_THRESH and NGENES >= N_GENE_THRESH
    significants = result[result["NGENES"] >= n_gene_thresh]

    if significants.empty:
        print("returning - no significants for " + gwas + " " + annot_type)
        return

    q_significant = significants[
        significants["Q " + annot_type + " " + compsets] < signif_thresh
    ]
    q_final = pd.merge(
        q_significant, atc_levels, right_index=True, left_on="ATC_CODE"
    ).sort_values("Q " + annot_type + " " + compsets)

    q_final.to_csv(
        os.path.join(
            summary_results_dir,
            "drugs_found-"
            + gwas
            + "-"
            + run_id
            + "-"
            + annot_type
            + "-"
            + compsets
            + "_qvals.tsv",
        ),
        sep="\t",
        index=False,
    )

    p_significant = significants[
        significants["P " + annot_type + " " + compsets] < signif_thresh
    ]
    p_final = pd.merge(
        p_significant, atc_levels, right_index=True, left_on="ATC_CODE"
    ).sort_values("P " + annot_type + " " + compsets)
    p_final.to_csv(
        os.path.join(
            summary_results_dir,
            "drugs_found-"
            + gwas
            + "-"
            + run_id
            + "-"
            + annot_type
            + "-"
            + compsets
            + "_pvals.tsv",
        ),
        sep="\t",
        index=False,
    )


def summarise_gopath_results(
    gwas,
    annot_type,
    atc_levels,
    magma_results_dir,
    summary_results_dir,
    run_id,
    magma_ver=106,
    n_gene_thresh=4,
    signif_thresh=0.05,
):
    result = read_results_files(
        gwas, annot_type, magma_results_dir, run_id, "gopaths", magma_ver
    )

    # only consider classes/drugs with Qval < QVAL_THRESH and NGENES >= N_GENE_THRESH
    significants = result[result["NGENES"] >= n_gene_thresh]

    if significants.empty:
        return

    significants.rename(
        columns={"ATC_CODE": "SHORT_NAME", "FULL_NAME": "ATC_CODE"}, inplace=True
    )

    q_significant = significants[
        significants["Q " + annot_type + " gopaths"] < signif_thresh
    ]
    q_significant.to_csv(
        os.path.join(
            summary_results_dir,
            "pathways_found-"
            + gwas
            + "-"
            + run_id
            + "-"
            + annot_type
            + "-go_qvals.tsv",
        ),
        sep="\t",
        index=False,
    )

    p_significant = significants[
        significants["P " + annot_type + " gopaths"] < signif_thresh
    ]
    p_significant.to_csv(
        os.path.join(
            summary_results_dir,
            "pathways_found-"
            + gwas
            + "-"
            + run_id
            + "-"
            + annot_type
            + "-go_pvals.tsv",
        ),
        sep="\t",
        index=False,
    )


def get_annot_results(gwas, annot, magma_results_dir, run_id):
    file = os.path.join(
        magma_results_dir,
        gwas,
        run_id,
        "magma_gene_result-" + gwas + "-" + annot + "-" + run_id + ".genes.out",
    )
    df = None

    if os.path.exists(file):
        df = pd.read_csv(file, delim_whitespace=True)
    else:
        display("File not found: " + file)

    return df


def add_genesets(table, zscores_df, genesets_df, gene_hgnc_df, gene_twas_df):
    table["Ensembl"] = ""
    table["HGNC"] = ""
    table["GWAS ZSTAT"] = ""
    table["Chrom."] = ""

    if len(gene_twas_df):
        table["TWAS ZSTAT"] = ""

    for index, row in table.iterrows():
        aset = genesets_df.loc[index].str.split(" ").tolist()
        setgenes_df = pd.DataFrame(aset).melt()[["value"]]
        setgenes_df.columns = ["Ensembl"]
        if not setgenes_df.empty:
            setgenes2_df = pd.merge(
                setgenes_df, zscores_df, how="left", left_on="Ensembl", right_on="GENE"
            )[["Ensembl", "ZSTAT"]]
            hugo_df = pd.merge(setgenes2_df, gene_hgnc_df, how="left", on="Ensembl")
            if len(gene_twas_df):
                hugo_df = pd.merge(
                    hugo_df, gene_twas_df, how="left", left_on="Ensembl", right_on="ID"
                )

            hugo_df["absZSTAT"] = hugo_df["ZSTAT"].abs()
            hugo_df.sort_values(by=["absZSTAT"], ascending=False, inplace=True)
            hugo_df.drop("absZSTAT", axis=1, inplace=True)
            hugo_df["ZSTAT"] = hugo_df["ZSTAT"].round(3).astype(str)
            hugo_df["chrom"] = hugo_df["chrom"].astype(str)

            # Limit to top 100 for now
            hugo_df = hugo_df.head(100)

            table.at[index, "Ensembl"] = "\n".join(hugo_df["Ensembl"].tolist())
            table.at[index, "HGNC"] = "\n".join(hugo_df["HGNC"].fillna("-?-").tolist())
            table.at[index, "GWAS ZSTAT"] = "\n".join(hugo_df["ZSTAT"].tolist())
            table.at[index, "Chrom."] = "\n".join(hugo_df["chrom"].tolist())

            if len(gene_twas_df):
                hugo_df["TWAS.Z"] = hugo_df["TWAS.Z"].round(3).astype(str)
                table.at[index, "TWAS ZSTAT"] = "\n".join(hugo_df["TWAS.Z"].tolist())

    # Hack for gosets
    if "SHORT_NAME" in table.columns:
        table.set_index("SHORT_NAME", inplace=True)
        table.index.rename("ATC_CODE", inplace=True)


# Display a set of results tables for the given GWAS list, file type and column (P or Q)
def display_tables(
    gwas,
    run_id,
    annot_type,
    file_postfix,
    column,
    zscores_df,
    summary_results_dir,
    genesets_df,
    gene_hgnc_df,
    limit_size=True,
    gene_twas_df=pd.DataFrame(),
):
    pval_files = [
        file
        for file in os.listdir(summary_results_dir)
        if file.endswith(file_postfix)
        and "found-" + gwas + "-" + run_id + "-" + annot_type + "-" in file
    ]
    pval_files.sort()

    for file in pval_files:
        table = pd.read_csv(os.path.join(summary_results_dir, file), sep="\t")
        table.sort_values(column, inplace=True)
        table.set_index("ATC_CODE", inplace=True)

        print("\033[1m\033[4m" + file[:-4] + "\033[0m - " + str(len(table)) + " items")
        if not table.empty:
            if limit_size == True:
                print(
                    "\033[1m\033[4m"
                    + "First 100 without geneset members displayed"
                    + "\033[0m"
                )
            
                if "SHORT_NAME" in table.columns:
                    display(table.drop("SHORT_NAME", axis=1).head(100))
                else:
                    display(table.head(100))
            else:
                if "SHORT_NAME" in table.columns:
                    display(table.drop("SHORT_NAME", axis=1))
                else:
                    display(table)
                
            add_genesets(table, zscores_df, genesets_df, gene_hgnc_df, gene_twas_df)
            if limit_size == True:
                print(
                    "\033[1m\033[4m"
                    + "First 10 with first 50 geneset members displayed"
                    + "\033[0m"
                )
                display(HTML(table.head(10).to_html().replace("\\n", "<br>")))
            else:
                display(HTML(table.to_html().replace("\\n", "<br>")))
        else:
            print("Empty table")


# -


def add_treatment_state(
    gwas,
    table,
):
    if gwas == "GA_I10":
        treatment_drugs = treatment_antihypertensives
        side_effects = drug_hypertensives
        side_effects.update(drug_bp_side_effects)
        affects_col = "Affects BP"
    elif gwas == "GA_E78":
        treatment_drugs = treatment_antihypercholesterol
        side_effects = drug_cholesterol
        side_effects.update(drug_cl_side_effects)
        affects_col = "Affects cholesterol"
    elif gwas == "PGC3_SZ":
        treatment_drugs = []
        side_effects = []
        affects_col = "N/A"       
    
    table["Treatment drug"] = [
            "yes" if x in treatment_drugs else "no" for x in table.index
        ]
    table[affects_col] = ["yes" if x in side_effects else "no" for x in table.index]
        
    for index, row in table.iterrows():
        if row["Treatment drug"] == "yes":
            table.at[index, affects_col] = "yes"
            
    return table


def display_joint_table_thesis(
    formatting,
    gwas,
    run_id,
    annot_type,
    file_postfix1,
    file_postfix2,
    column1,
    column2,
    zscores_df,
    summary_results_dir,
    genesets_df,
    gene_hgnc_df,
    bold_caption,
    normal_caption,
    latex_label,
    with_genes=False,
    n_gene_thresh=4,
    signif_thresh=0.05,
):
    
    pval_files1 = [
        file
        for file in os.listdir(summary_results_dir)
        if file.endswith(file_postfix1)
        and "found-" + gwas + "-" + run_id + "-" + annot_type in file
    ]

    pval_files2 = [
        file
        for file in os.listdir(summary_results_dir)
        if file.endswith(file_postfix2)
        and "found-" + gwas + "-" + run_id + "-" + annot_type in file
    ]

    table1 = pd.read_csv(os.path.join(summary_results_dir, pval_files1[0]), sep="\t")
    table2 = pd.read_csv(os.path.join(summary_results_dir, pval_files2[0]), sep="\t")

    table = pd.merge(table1, table2, how="outer", on="ATC_CODE")

    table["Description"] = table["Description_x"].fillna(table["Description_y"])
    table.drop(["Description_x", "Description_y"], axis=1, inplace=True)

    table = table[table["NGENES_x"] >= n_gene_thresh]
    table = table[
        (table["Q prox duggie"] <= signif_thresh)
        | (table["Q prox stitch"] <= signif_thresh)
    ]

    if not table.empty:
        if "SHORT_NAME" in table.columns:
            table.drop("SHORT_NAME", axis=1, inplace=True)

        table.set_index("ATC_CODE", inplace=True)
        table.index.rename("ATC code", inplace=True)
        table["NGENES_x"] = [str(int(x)) if x > 0 else "-" for x in table["NGENES_x"]]
        table["Q prox duggie"] = [
            x if ~np.isnan(x) else "-" for x in table["Q prox duggie"]
        ]
        table["Q prox duggie"] = [
            x if float(x) <= signif_thresh else "(" + "{:.3e}".format(x) + ")"
            for x in table["Q prox duggie"]
        ]
        table["NGENES_y"] = [str(int(x)) if x > 0 else "-" for x in table["NGENES_y"]]
        table["Q prox stitch"] = [
            x if ~np.isnan(x) else 0 for x in table["Q prox stitch"]
        ]
        table["P prox stitch"] = [
            x if ~np.isnan(x) else 0 for x in table["P prox stitch"]
        ]
        table["Q prox stitch"] = [
            x if float(x) <= signif_thresh else "(" + "{:.3e}".format(x) + ")"
            for x in table["Q prox stitch"]
        ]
        table["Q prox stitch"] = [x if x != 0 else "-" for x in table["Q prox stitch"]]
        table["P prox stitch"] = [x if x != 0 else "-" for x in table["P prox stitch"]]

        table.columns = [
            "Num genes DUGGIE",
            "P DUGGIE",
            "Q DUGGIE",
            "Num genes STITCH",
            "P STITCH",
            "Q STITCH",
            "Description",
        ]
        table = table.reindex(
            [
                "Description",
                "Num genes DUGGIE",
                "P DUGGIE",
                "Q DUGGIE",
                "Num genes STITCH",
                "P STITCH",
                "Q STITCH",
            ],
            axis=1,
        )

        table = add_treatment_state(gwas, table)

        if with_genes:
            table = add_genesets_thesis(table, zscores_df, genesets_df, gene_hgnc_df)
            table.columns = ["ATC code", "Num genes", "Q", "Description", "HGNC", "Z"]
            if formatting == "LaTeX":
                display(Latex(r"\scriptsize"))
                display(
                    Latex(
                        table.to_latex(
                            column_format="p{2cm}p{2cm}p{2cm}p{4cm}p{2cm}p{2cm}p{2cm}",
                            multirow=True,
                            multicolumn=True,
                            index=False,
                            caption=bold_caption
                            + r"\normalfont{"
                            + normal_caption
                            + "}",
                            label=latex_label,
                        )
                    )
                )
                display(Latex(r"\normalsize"))
            else:
                display(bold_caption + " " + normal_caption)
                display(HTML(table.to_html()))

        else:
            table.drop(['P DUGGIE', 'P STITCH'], axis=1, inplace=True)
            if formatting == "LaTeX":
                display(Latex(r"\scriptsize"))
                display(
                    Latex(
                        table.to_latex(
                            column_format="lR{3.4cm}p{1.6cm}R{1.5cm}p{1.6cm}R{1.5cm}p{1.2cm}p{1.2cm}",
                            multirow=True,
                            multicolumn=True,
                            caption=bold_caption
                            + r"\normalfont{"
                            + normal_caption
                            + "}",
                            label=latex_label,
                        )
                    )
                )
                display(Latex(r"\normalsize"))
            else:
                display(bold_caption + " " + normal_caption)
                display(HTML(table.to_html()))


def display_joint_table_thesis_ch4(
    formatting,
    gwas,
    run_id,
    annot_type,
    file_postfix1,
    column1,
    zscores_df,
    summary_results_dir,
    genesets_df,
    gene_hgnc_df,
    bold_caption,
    normal_caption,
    latex_label,
    with_genes=False,
    n_gene_thresh=4,
    signif_thresh=0.05,
):
    
    pval_files1 = [
        file
        for file in os.listdir(summary_results_dir)
        if file.endswith(file_postfix1)
        and "found-" + gwas + "-" + run_id + "-" + annot_type in file
    ]

    table1 = pd.read_csv(os.path.join(summary_results_dir, pval_files1[0]), sep="\t")

    table = table1

    table = table[table["NGENES"] >= n_gene_thresh]
    table = table[
        (table["Q " + column1] <= signif_thresh)
    ]

    if not table.empty:
        if "SHORT_NAME" in table.columns:
            table.drop("SHORT_NAME", axis=1, inplace=True)

        table.set_index("ATC_CODE", inplace=True)
        table.index.rename("ATC code", inplace=True)
        table["NGENES"] = [str(int(x)) if x > 0 else "-" for x in table["NGENES"]]
        table["Q " + column1] = [
            x if ~np.isnan(x) else "-" for x in table["Q " + column1]
        ]
        table["Q " + column1] = [
            x if float(x) <= signif_thresh else "(" + "{:.3e}".format(x) + ")"
            for x in table["Q " + column1]
        ]

        table.columns = [
            "Num genes DUGGIE",
            "P DUGGIE",
            "Q DUGGIE",
            "Description",
        ]
        table = table.reindex(
            [
                "Description",
                "Num genes DUGGIE",
                "P DUGGIE",
                "Q DUGGIE",
            ],
            axis=1,
        )
        
        table = add_treatment_state(gwas, table)

        if with_genes:
            table = add_genesets_thesis(table, zscores_df, genesets_df, gene_hgnc_df)
            table.columns = ["ATC code", "Num genes", "Q", "Description", "HGNC", "Z"]
            if formatting == "LaTeX":
                display(Latex(r"\scriptsize"))
                display(
                    Latex(
                        table.to_latex(
                            column_format="p{2cm}p{2cm}p{2cm}p{4cm}p{2cm}p{2cm}p{2cm}",
                            multirow=True,
                            multicolumn=True,
                            index=False,
                            caption=bold_caption
                            + r"\normalfont{"
                            + normal_caption
                            + "}",
                            label=latex_label,
                        )
                    )
                )
                display(Latex(r"\normalsize"))
            else:
                display(bold_caption + " " + normal_caption)
                display(HTML(table.to_html()))

        else:
            if formatting == "LaTeX":
                display(Latex(r"\scriptsize"))
                display(
                    Latex(
                        table.to_latex(
                            column_format="p{1.3cm}p{3.4cm}p{1.2cm}p{2cm}p{2cm}p{1.2cm}p{1.2cm}",
                            multirow=True,
                            multicolumn=True,
                            caption=bold_caption
                            + r"\normalfont{"
                            + normal_caption
                            + "}",
                            label=latex_label,
                        )
                    )
                )
                display(Latex(r"\normalsize"))
            else:
                display(bold_caption + " " + normal_caption)
                display(HTML(table.to_html()))


def display_dgi_venn_thesis(
    gwas,
    run_id,
    annot_type,
    file_postfix1,
    file_postfix2,
    column1,
    column2,
    zscores_df,
    summary_results_dir,
    genesets_df,
    gene_hgnc_df,
    n_gene_thresh=4,
    signif_thresh=0.05,
):
    if gwas == "GA_I10":
        treatment_drugs = treatment_antihypertensives
        name = "hypertension"
    elif gwas == "GA_E78":
        treatment_drugs = treatment_antihypercholesterol
        name = "hypercholesterolemia"

    pval_files1 = [
        file
        for file in os.listdir(summary_results_dir)
        if file.endswith(file_postfix1)
        and "found-" + gwas + "-" + run_id + "-" + annot_type in file
    ]

    pval_files2 = [
        file
        for file in os.listdir(summary_results_dir)
        if file.endswith(file_postfix2)
        and "found-" + gwas + "-" + run_id + "-" + annot_type in file
    ]

    table1 = pd.read_csv(os.path.join(summary_results_dir, pval_files1[0]), sep="\t")
    table2 = pd.read_csv(os.path.join(summary_results_dir, pval_files2[0]), sep="\t")
    table1["Treatment drug"] = [
        1 if x in treatment_drugs else 0 for x in table1["ATC_CODE"]
    ]
    table2["Treatment drug"] = [
        1 if x in treatment_drugs else 0 for x in table2["ATC_CODE"]
    ]

    table3 = table1[(table1["Q prox duggie"] <= signif_thresh)]
    table4 = table2[(table2["Q prox stitch"] <= signif_thresh)]
    set1 = set(table3["ATC_CODE"])
    set2 = set(table4["ATC_CODE"])
    set3 = set(table3[(table3["Treatment drug"] == 1)]["ATC_CODE"])
    set4 = set(table4[(table4["Treatment drug"] == 1)]["ATC_CODE"])

    plt.figure(figsize=(14, 14)).subplot_mosaic([[1, 2]])
    plt.subplot(2, 2, 1)
    plt.title("(a) all " + name + " drugs")
    v1 = venn2([set1, set2], ("DUGGIE", "STITCH"))
    plt.subplot(2, 2, 2)
    plt.title("(b) " + name + " treatment drugs")
    v2 = venn2([set3, set4], ("DUGGIE", "STITCH"))


def tally_genesets_thesis(table, zscores_df, genesets_df, gene_hgnc_df):
    table2 = pd.DataFrame(columns=["HGNC", "chrom", "start_pos", "count", "Z"])
    table2.set_index("HGNC", inplace=True)
    zscores_df.set_index("GENE", inplace=True)

    for index, row in table.iterrows():
        aset = genesets_df.loc[index].str.split(" ").tolist()
        setgenes_df = pd.DataFrame(aset).melt()[["value"]]
        setgenes_df.columns = ["Ensembl"]
        if not setgenes_df.empty:
            hugo_df = pd.merge(setgenes_df, gene_hgnc_df, how="left", on="Ensembl")

            for index, row in hugo_df.iterrows():
                gene = row["HGNC"]
                if gene in table2.index:
                    table2.at[gene, "count"] = table2.at[gene, "count"] + 1
                else:
                    table2.at[gene, "count"] = 1
                    if row["Ensembl"] in zscores_df.index:
                        table2.at[gene, "Z"] = str(
                            zscores_df.at[row["Ensembl"], "ZSTAT"]
                        )
                        table2.at[gene, "chrom"] = str(
                            zscores_df.at[row["Ensembl"], "CHR"]
                        )
                        table2.at[gene, "start_pos"] = str(
                            zscores_df.at[row["Ensembl"], "START"]
                        )

    return table2


def tally_top_genes_thesis(table, zscores_df, genesets_df, gene_hgnc_df):
    table2 = pd.DataFrame(columns=["HGNC", "count", "Z"])
    table2.set_index("HGNC", inplace=True)
    zscores_df.set_index("GENE", inplace=True)

    for index, row in table.iterrows():
        aset = genesets_df.loc[index].str.split(" ").tolist()
        setgenes_df = pd.DataFrame(aset).melt()[["value"]]
        setgenes_df.columns = ["Ensembl"]
        if not setgenes_df.empty:
            setgenes2_df = pd.merge(
                setgenes_df, zscores_df, how="left", left_on="Ensembl", right_on="GENE"
            )[["Ensembl", "ZSTAT"]]

            hugo_df = pd.merge(setgenes2_df, gene_hgnc_df, how="left", on="Ensembl")
            hugo_df["absZSTAT"] = hugo_df["ZSTAT"].abs()
            hugo_df.sort_values(by=["absZSTAT"], ascending=False, inplace=True)
            hugo_df.drop("absZSTAT", axis=1, inplace=True)

            gene = hugo_df["HGNC"][0]
            row = hugo_df.head(1)
            if gene in table2.index:
                table2.at[gene, "count"] = table2.at[gene, "count"] + 1
            else:
                table2.at[gene, "count"] = 1
                if hugo_df["Ensembl"][0] in zscores_df.index:
                    table2.at[gene, "Z"] = str(
                        zscores_df.at[hugo_df["Ensembl"][0], "ZSTAT"]
                    )

    return table2


# +
def add_genesets_thesis(table, zscores_df, genesets_df, gene_hgnc_df):
    table2 = pd.DataFrame(columns=table.columns)

    for index, row in table.iterrows():
        aset = genesets_df.loc[index].str.split(" ").tolist()
        setgenes_df = pd.DataFrame(aset).melt()[["value"]]
        setgenes_df.columns = ["Ensembl"]
        if not setgenes_df.empty:
            setgenes2_df = pd.merge(
                setgenes_df, zscores_df, how="left", left_on="Ensembl", right_on="GENE"
            )[["Ensembl", "ZSTAT"]]
            hugo_df = pd.merge(setgenes2_df, gene_hgnc_df, how="left", on="Ensembl")
            hugo_df.drop(["Ensembl"], axis=1, inplace=True)

            hugo_df["absZSTAT"] = hugo_df["ZSTAT"].abs()
            hugo_df.sort_values(by=["absZSTAT"], ascending=False, inplace=True)
            hugo_df.drop("absZSTAT", axis=1, inplace=True)
            hugo_df["ZSTAT"] = hugo_df["ZSTAT"].round(5).astype(str)
            hugo_df.reset_index(inplace=True)

            table2 = table2.append(row)
            table2.at[index, "HGNC"] = hugo_df["HGNC"][0]
            table2.at[index, "ZSTAT"] = hugo_df["ZSTAT"][0]
            table2.at[index, "ATC code"] = index
            table2 = table2.append(hugo_df[1:])

    table2.fillna(" ", inplace=True)
    table2 = table2.reindex(
        ["ATC code", "Num genes", "Q", "Description", "HGNC", "ZSTAT"], axis=1
    )

    return table2


def display_table_thesis(
    gwas,
    run_id,
    annot_type,
    file_postfix,
    column,
    zscores_df,
    summary_results_dir,
    genesets_df,
    gene_hgnc_df,
    caption,
    latex_label,
    with_genes=False,
):
    pval_files = [
        file
        for file in os.listdir(summary_results_dir)
        if file.endswith(file_postfix)
        and "found-" + gwas + "-" + run_id + "-" + annot_type in file
    ]
    pval_files.sort()

    for file in pval_files:
        table = pd.read_csv(os.path.join(summary_results_dir, file), sep="\t")
        table.sort_values(column, inplace=True)
        table.set_index("ATC_CODE", inplace=True)

        if not table.empty:
            if "SHORT_NAME" in table.columns:
                table.drop("SHORT_NAME", axis=1, inplace=True)

            table.index.rename("ATC code", inplace=True)
            table.columns = ["Num genes", "P", "Q", "Description"]

            if with_genes:
                table.drop(["P"], axis=1, inplace=True)
                tally_df = tally_top_genes_thesis(
                    table, zscores_df, genesets_df, gene_hgnc_df
                )
                table = add_genesets_thesis(
                    table, zscores_df, genesets_df, gene_hgnc_df
                )
                table.columns = [
                    "ATC code",
                    "Num genes",
                    "Q",
                    "Description",
                    "HGNC",
                    "Z",
                ]
                display(
                    Latex(
                        table.to_latex(
                            column_format="p{2cm}p{2cm}p{2cm}p{4cm}p{2cm}p{2cm}",
                            multirow=True,
                            multicolumn=True,
                            index=False,
                            caption=caption,
                            label=latex_label,
                        )
                    )
                )
                return tally_df
            else:

                display(
                    Latex(
                        table.to_latex(
                            column_format="p{2cm}p{1cm}p{3cm}p{3cm}p{6cm}",
                            multirow=True,
                            multicolumn=True,
                            caption=caption,
                            label=latex_label,
                        )
                    )
                )


# -

# ### p-value histograms and correlations

# +
# Convert a set of p-values to q-values and a pi0 estimate, using the R implementation
def get_r_pi0est(pvalues):
    r_result = rqvalue.pi0est(np.array(pvalues))
    with localconverter(ro.default_converter + pandas2ri.converter):
        result = ro.conversion.py2rpy(r_result)
    return result[0][0]


def get_r_qvalues(pvalues):
    result = rqvalue.qvalue(np.array(pvalues))
    return np.array(result[2])


def get_r_hist_qvalues(pvalues):
    result = rqvalue.hist_qvalue(rqvalue.qvalue(np.array(pvalues)))
    return result


# +
def plot_density_histogram(gwas, annot, analysis, df, pval_type):
    if df is not None:
        sig_range = 3
        # The first (sig_range) bins represent p < 0.05, which will be coloured orange
        N, bins, patches = plt.hist(df[pval_type], bins=20 * sig_range, density=True)
        [x.set_facecolor("orange") for x in patches[0:sig_range]]
        pi0 = get_r_pi0est(df[pval_type])
        plt.axhline(
            pi0, color="red", dashes=(5, 5), label="pi0est = " + str(format(pi0, ".3f"))
        )
        qv = get_r_qvalues(df[pval_type])
        plt.xlabel(pval_type)
        plt.ylabel("Density")
        plt.title(gwas + " " + annot + " " + analysis + " p-value density distribution")
        plt.legend()


def plot_gene_scatter(
    gwas,
    annot1,
    df1,
    annot2,
    df2,
    pval_label,
    axlimit,
    pval_label2="",
):
    if pval_label2 == "":
        pval_label2 = pval_label + "_y"
    else:
        pval_label2 = pval_label2 + "_y"

    pval_label = pval_label + "_x"

    df = pd.merge(df1, df2, how="inner", on="GENE")[["GENE", pval_label, pval_label2]]

    log_label1 = annot1 + " -log(" + pval_label + ")"
    log_label2 = annot2 + " -log(" + pval_label2 + ")"
    df.columns = ["GENE", log_label1, log_label2]

    # take -log10 of Pgene
    df[log_label1] = -np.log10(df[log_label1])
    df[log_label2] = -np.log10(df[log_label2])

    # linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        df[log_label1], df[log_label2]
    )

    # plotting
    diag = np.linspace(0, 160, 500)
    plt.plot(diag, diag, "-g")

    plt.scatter(df[log_label1], df[log_label2])
    plt.xlabel(log_label1)
    plt.ylabel(log_label2)
    plt.title(gwas + " " + annot1 + "/" + annot2 + " correlation")
    plt.xlim(0, axlimit)
    plt.ylim(0, axlimit)
    ax2 = plt.gca()
    plt.text(
        0,
        -0.25,
        "r_value = "
        + "{0:.3f}".format(r_value)
        + ", std err "
        + "{0:.3f}".format(std_err),
        transform=ax2.transAxes,
    )


# -


def group_drugs(
    df,
    gwas,
):
    if gwas == "GA_E78":
        group_labels = ["non-treatment drugs", "treatment drugs"]
        colours = ["pink", "blue"]
        for index, row in df.iterrows():
            if row["ATC_CODE"] in treatment_antihypercholesterol:
                colour_id = 1
            else:
                colour_id = 0

            df.at[index, "colour_id"] = colour_id

    elif gwas == "GA_I10" or gwas == "GA_I25":
        group_labels = ["non-treatment drugs", "treatment drugs"]
        colours = ["pink", "blue"]
        for index, row in df.iterrows():
            if row["ATC_CODE"] in treatment_antihypertensives:
                colour_id = 1
            else:
                colour_id = 0

            df.at[index, "colour_id"] = colour_id

    elif gwas == "SZ" or gwas == "PGC3_SZ":
        group_labels = ["non-treatment drugs", "antipsychotics"]
        colours = ["pink", "blue", "lightblue", "green"]
        for index, row in df.iterrows():
            if row["ATC_CODE"].startswith("N05A"):
                colour_id = 1
            else:
                colour_id = 0

            df.at[index, "colour_id"] = colour_id
    elif gwas == "GA_E10":
        group_labels = ["non-treatment drugs", "diabetes drugs"]
        colours = ["pink", "blue"]
        for index, row in df.iterrows():
            if row["ATC_CODE"].startswith("A10"):
                colour_id = 1
            else:
                colour_id = 0

            df.at[index, "colour_id"] = colour_id
    elif gwas == "GA_K50":
        group_labels = ["non-treatment drugs", "immunosupressants", "antimetabolites"]
        colours = ["pink", "blue", "lightblue"]
        for index, row in df.iterrows():
            if row["ATC_CODE"].startswith("L04A"):
                colour_id = 1
            elif row["ATC_CODE"].startswith("L01BB02"):
                colour_id = 2
            else:
                colour_id = 0

            df.at[index, "colour_id"] = colour_id
    elif gwas == "GA_C50":
        group_labels = ["non-treatment drugs", "antineoplasics", "hormones"]
        colours = ["pink", "blue", "lightblue"]
        for index, row in df.iterrows():
            if row["ATC_CODE"].startswith("L"):
                colour_id = 1
            elif row["ATC_CODE"].startswith("G03"):
                colour_id = 2
            else:
                colour_id = 0

            df.at[index, "colour_id"] = colour_id
    else:
        group_labels = ["non-treatment drugs"]
        colours = ["pink"]
        for index, row in df.iterrows():
            colour_id = 0

            df.at[index, "colour_id"] = colour_id

    return df, group_labels, colours


def plot_drug_scatter(gwas, annot1, df1, annot2, df2, pval_label, axlimit):
    df = pd.merge(df1, df2, how="inner", left_on="ATC_CODE", right_on="ATC_CODE")[
        [
            "ATC_CODE",
            pval_label + " " + annot1,
            pval_label + " " + annot2,
            "Q " + annot1,
            "Q " + annot2,
        ]
    ]

    log_label1 = annot1.upper() + " -log10(" + pval_label + ")"
    log_label2 = annot2.upper() + " -log10(" + pval_label + ")"
    df.columns = ["ATC_CODE", log_label1, log_label2, "Q " + annot1, "Q " + annot2]

    # take -log10 of Pgene
    df[log_label1] = -np.log10(df[log_label1])
    df[log_label2] = -np.log10(df[log_label2])

    # linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        df[log_label1], df[log_label2]
    )

    # Check Q state of both stitch and duggie
    # 0 = Neither significant, 1 = annot1 significant, 2 = annot2 significant, 3 = both significant
    markers = [".", "v", "^", "D"]
    for index, row in df.iterrows():
        state = 0
        if row["Q " + annot1] < 0.05:
            state += 1
        if row["Q " + annot2] < 0.05:
            state += 2

        df.at[index, "Q state"] = state

    df, group_labels, colours = group_drugs(df, gwas)
    groups = df.groupby(["colour_id", "Q state"])

    for name, group in groups:
        plt.scatter(
            group[log_label1],
            group[log_label2],
            label=group_labels[int(name[0])],
            c=colours[int(name[0])],
            marker=markers[int(name[1])],
        )

    diag = np.linspace(0, 160, 500)
    plt.plot(diag, diag, "-g")

    plt.xlabel(log_label1)
    plt.ylabel(log_label2)
    plt.title(gwas + " " + annot1 + "/" + annot2 + " correlation")
    plt.xlim(0, axlimit)
    plt.ylim(0, axlimit)
    ax2 = plt.gca()
    plt.text(
        0,
        -0.25,
        "r_value = "
        + "{0:.3f}".format(r_value)
        + ", std err "
        + "{0:.3f}".format(std_err),
        transform=ax2.transAxes,
    )


def plot_drug_genediff_scatter(gwas, annot1, df1, annot2, df2, pval_label, axlimit):

    df = pd.merge(
        df1, df2, how="inner", on="ATC_CODE", suffixes=[" " + annot1, " " + annot2]
    )
    log_label1 = annot1 + " -log(" + pval_label + ")"
    log_label2 = annot2 + " -log(" + pval_label + ")"

    df["diff"] = df["NGENES " + annot2] - df["NGENES " + annot1]

    # Hack to deal with zero p values in PERMP output
    df[log_label1] = df[pval_label + " " + annot1].replace(0, 0.00001)
    df[log_label2] = df[pval_label + " " + annot2].replace(0, 0.00001)

    # take -log10 of Pgene
    df[log_label1] = -np.log10(df[log_label1])
    df[log_label2] = -np.log10(df[log_label2])

    # linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        df[log_label1], df[log_label2]
    )

    # plotting
    group_labels = ["diff<4", "4<=diff<11", "11<=diff<22", "diff>=22"]
    for index, row in df.iterrows():
        if row["diff"] < 4:
            group_id = 0
        elif row["diff"] < 11:
            group_id = 1
        elif row["diff"] < 22:
            group_id = 2
        else:
            group_id = 3

        df.at[index, "group_id"] = group_id
        df.at[index, "group"] = group_labels[group_id]

    groups = df.groupby("group_id")

    for name, group in groups:
        plt.scatter(group[log_label1], group[log_label2], label=group_labels[int(name)])

    diag = np.linspace(0, 160, 500)
    plt.plot(diag, diag, "-g")

    plt.xlabel(log_label1)
    plt.ylabel(log_label2)
    plt.title(gwas + " " + annot1 + "/" + annot2 + " correlation")
    plt.xlim(0, axlimit)
    plt.ylim(0, axlimit)
    ax2 = plt.gca()
    plt.text(
        0,
        -0.25,
        "r_value = "
        + "{0:.3f}".format(r_value)
        + ", std err "
        + "{0:.3f}".format(std_err),
        transform=ax2.transAxes,
    )

    legend = ax2.legend([])
    plt.legend()

    return df


def show_treatment_drug_stats(dfs, db_col, sig):

    columns = [
        "DTI DB",
        "Annotation",
        "Sig.",
        "Sig. Treatment",
        "Sig. Non-treatment",
        "Non-sig.",
        "Non-sig. Treatment",
        "Non-sig. Non-treatment",
        "Sens.",
        "Spec.",
        "FDR",
        "Mann-Whitney p-value",
        "Fisher p-value",
    ]
    results_df = pd.DataFrame(columns=columns)
    table_fdr = dfs[0].columns[db_col].split()[0]

    for i in range(0, len(dfs)):
        df = dfs[i]

        treatment_df = df[df["colour_id"] != 0]
        TP_df = treatment_df[treatment_df.iloc[:, db_col] < sig]
        FN_df = treatment_df[treatment_df.iloc[:, db_col] >= sig]
        nontreatment_df = df[df["colour_id"] == 0]
        FP_df = nontreatment_df[nontreatment_df.iloc[:, db_col] < sig]
        TN_df = nontreatment_df[nontreatment_df.iloc[:, db_col] >= sig]

        table_dti = treatment_df.columns[db_col].split()[2]
        table_annot = treatment_df.columns[db_col].split()[1]

        try:
            # population difference between treatment and non-treatment drugs
            table_diff_mannwhitneyp = "{:.2e}".format(
                stats.mannwhitneyu(
                    nontreatment_df.iloc[:, db_col], treatment_df.iloc[:, db_col]
                ).pvalue
            )
        except:
            table_diff_mannwhitneyp = "ERROR"

        fisher_contingency = "{:.2e}".format(
            stats.fisher_exact([[len(TP_df), len(FP_df)], [len(FN_df), len(TN_df)]])[1]
        )
        if len(TP_df) or len(FN_df):
            sensitivity = len(TP_df) / (len(TP_df) + len(FN_df))
        else:
            sensitivity = 0

        specificity = len(TN_df) / (len(TN_df) + len(FP_df))

        if  ((len(FP_df) + len(TP_df)) == 0):
            fdr = 1
        else:
            fdr = len(FP_df) / (len(FP_df) + len(TP_df))

        new_row = {
            "DTI DB": table_dti,
            "Annotation": table_annot,
            "Sig.": len(TP_df) + len(FP_df),
            "Sig. Treatment": len(TP_df),
            "Sig. Non-treatment": len(FP_df),
            "Non-sig.": len(TN_df) + len(FN_df),
            "Non-sig. Treatment": len(FN_df),
            "Non-sig. Non-treatment": len(TN_df),
            "Sens.": sensitivity,
            "Spec.": specificity,
            "FDR": fdr,
            "Mann-Whitney p-value": table_diff_mannwhitneyp,
            "Fisher p-value": fisher_contingency,
        }

        results_df = results_df.append(new_row, ignore_index=True)

    return results_df


# define Jaccard similarity function
def jaccard(set1, set2):
    intersection = len(list((set1).intersection(set2)))
    union = (len(list(set1)) + len(list(set2))) - intersection

    if union:
        return float(intersection) / union
    else:
        return 0


def load_drug_sets_ch5(dtis, 
                       annots, 
                       display_annots, 
                       summary_results_dir, 
                       gwas, 
                       run_id, 
                       signif_thresh=0.05):
    
    drug_sets_df = pd.DataFrame()

    for dti in dtis:
        file_postfix = dti + "_qvals.tsv"
        for annot in annots:
            label = dti.upper() + " " + display_annots[annots.index(annot)]
            pval_file = [
                file
                for file in os.listdir(summary_results_dir)
                if file.endswith(file_postfix)
                and "found-" + gwas + "-" + run_id + "-" + annot + "-" in file
            ]
            table = pd.read_csv(os.path.join(summary_results_dir, pval_file[0]), sep="\t")
            table2 = table[table["Q " + annot + " " + dti] < signif_thresh]
            drug_set = table2["ATC_CODE"]
            # Add series as column to get a matrix of boolean columns
            drug_sets_df[label] = False
            for drug in drug_set:
                if drug not in drug_sets_df.index:
                    drug_sets_df = drug_sets_df.append(pd.Series([False], name=drug))
                drug_sets_df.at[drug, label] = True
        
    drug_sets_df = drug_sets_df.fillna(False)
    drug_sets_df = drug_sets_df.drop(0, axis=1)
    
    return drug_sets_df


def get_top_drugs_ch5(drug_sets_df, db_indications_df, ATC_LEVELS, novel_drugs):
    top_drugs_df = pd.DataFrame(drug_sets_df.transpose().sum(), columns=['Count'])
    top_drugs_df = top_drugs_df[top_drugs_df['Count'] > 0]
    top_drugs_df = pd.merge(top_drugs_df, ATC_LEVELS, left_index=True, right_index=True)
    top_drugs_df.index.set_names(["ATC"], inplace=True)
    top_drugs_df.sort_values(by=['Count', 'ATC'], ascending=[False, True], inplace=True)   
    top_drugs_df["Novel"] = np.where(top_drugs_df.index.isin(novel_drugs), "yes", "no")
    top_drugs_df = pd.merge(top_drugs_df, db_indications_df, left_index=True, right_index=True, how='left')                                                                                       
    top_drugs_df.rename(columns={0:'Indications'}, inplace=True)
    top_drugs_df.index.set_names(["ATC"], inplace=True)

    return top_drugs_df
