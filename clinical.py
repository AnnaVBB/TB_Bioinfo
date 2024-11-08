import pandas as pd
import tarfile 

tar_path= r'C:\Users\Usuário\Downloads\Bioinformatica\Projeto\clinical.project-tcga-lgg.2024-11-07.tar.gz'

with tarfile.open(tar_path, 'r:gz') as tar_ref:
    tar_ref.extractall(r'C:\Users\Usuário\Downloads\Bioinformatica\Projeto')

# Carregar cada arquivo TSV em um DataFrame do Pandas
clinical_df = pd.read_csv(r'C:\Users\Usuário\Downloads\Bioinformatica\Projeto\clinical.tsv', sep='\t')
exposure_df = pd.read_csv(r'C:\Users\Usuário\Downloads\Bioinformatica\Projeto\exposure.tsv', sep='\t')
family_history_df = pd.read_csv(r'C:\Users\Usuário\Downloads\Bioinformatica\Projeto\family_history.tsv', sep='\t')
follow_up_df = pd.read_csv(r'C:\Users\Usuário\Downloads\Bioinformatica\Projeto\follow_up.tsv', sep='\t')
pathology_detail_df = pd.read_csv(r'C:\Users\Usuário\Downloads\Bioinformatica\Projeto\pathology_detail.tsv', sep='\t')

# Verifique os primeiros registros de cada DataFrame
pd.set_option('display.max_columns', None)
#print(clinical_df.head())
#print(clinical_df.columns.tolist())

#Remover colunas com muitos valores ausentes:
threshold = len(clinical_df) * 0.5
clinical_df_cleaned = clinical_df.dropna(axis=1, thresh=threshold)
columns_to_remove = ['country_of_birth', 'country_of_residence_at_enrollment', 'education_level', 'marital_status','occupation_duration_years','premature_at_birth',
'vital_status', 'weeks_gestation_at_birth','days_to_birth', 'days_to_death','adrenal_hormone','ann_arbor_b_symptoms', 'ann_arbor_b_symptoms_described', 'ann_arbor_clinical_stage', 'ann_arbor_extranodal_involvement', 'ann_arbor_pathologic_stage', 
'burkitt_lymphoma_clinical_variant', 'calgb_risk_group','child_pugh_classification', 'clark_level', 'cog_liver_stage', 'days_to_best_overall_response', 'days_to_diagnosis', 'days_to_last_follow_up', 'days_to_last_known_disease_status', 'days_to_recurrence',
'diagnosis_is_primary_disease', 'double_expressor_lymphoma', 'double_hit_lymphoma', 'eln_risk_classification', 'enneking_msts_grade', 'enneking_msts_metastasis', 'enneking_msts_stage', 'enneking_msts_tumor_site', 'ensat_clinical_m', 'ensat_pathologic_n',
'ensat_pathologic_stage', 'ensat_pathologic_t', 'esophageal_columnar_dysplasia_degree', 'esophageal_columnar_metaplasia_present', 'fab_morphology_code','figo_stage', 'figo_staging_edition_year', 'first_symptom_longest_duration', 'first_symptom_prior_to_diagnosis', 
'gastric_esophageal_junction_involvement', 'gleason_grade_group', 'gleason_grade_tertiary', 'gleason_patterns_percent', 'gleason_score', 'goblet_cells_columnar_mucosa_present', 'igcccg_stage', 'inpc_grade', 'inpc_histologic_group', 'inrg_stage', 'inss_stage',
'international_prognostic_index', 'irs_group', 'irs_stage', 'ishak_fibrosis_score', 'iss_stage', 'last_known_disease_status', 'laterality', 'masaoka_stage', 'max_tumor_bulk_site','medulloblastoma_molecular_classification', 'melanoma_known_primary',
'method_of_diagnosis', 'micropapillary_features', 'mitosis_karyorrhexis_index', 'mitotic_count', 'morphology', 'ovarian_specimen_status', 'ovarian_surface_involvement', 'papillary_renal_cell_type', 'pediatric_kidney_staging', 'peritoneal_fluid_cytological_status',
'pregnant_at_diagnosis', 'primary_diagnosis', 'primary_disease', 'primary_gleason_grade', 'prior_malignancy','residual_disease', 'satellite_nodule_present', 'secondary_gleason_grade','uicc_clinical_m', 'uicc_clinical_n', 'uicc_clinical_stage', 'uicc_clinical_t', 
'uicc_pathologic_m', 'uicc_pathologic_n', 'uicc_pathologic_stage', 'uicc_pathologic_t', 'uicc_staging_system_edition', 'ulceration_indicator', 'weiss_assessment_findings', 'weiss_assessment_score', 'who_nte_grade', 'wilms_tumor_histologic_subtype', 'breslow_thickness', 
'circumferential_resection_margin', 'gross_tumor_weight', 'largest_extrapelvic_peritoneal_focus', 'lymph_node_involved_site', 'lymph_nodes_positive', 'lymph_nodes_tested', 'lymphatic_invasion_present', 'non_nodal_regional_disease', 'non_nodal_tumor_deposits', 'percent_tumor_invasion', 
'perineural_invasion_present', 'peripancreatic_lymph_nodes_positive', 'peripancreatic_lymph_nodes_tested', 'transglottic_extension',
'vascular_invasion_present', 'vascular_invasion_type', 'clinical_trial_indicator', 'course_number', 'days_to_treatment_end', 'days_to_treatment_start', 'drug_category', 'embolic_agent', 'lesions_treated_number', 'margin_distance.1', 
'margin_status', 'margins_involved_site.1', 'number_of_fractions', 'prescribed_dose', 'prescribed_dose_units', 'pretreatment', 'protocol_identifier', 'radiosensitizing_agent', 'regimen_or_line_of_therapy', 'residual_disease.1', 
'route_of_administration', 'therapeutic_level_achieved', 'therapeutic_levels_achieved', 'therapeutic_target_level','treatment_anatomic_site', 'treatment_anatomic_sites', 'treatment_arm', 'treatment_dose', 'treatment_dose_max', 'treatment_dose_units', 'treatment_duration', 
'treatment_frequency', 'treatment_intent_type', 'treatment_or_therapy', 'treatment_outcome_duration']  # Substitua pelos nomes das colunas reais
clinical_df_cleaned = clinical_df.drop(columns=columns_to_remove)

#print(clinical_df_cleaned.head())
print("Colunas restantes:", clinical_df_cleaned.columns.tolist())
