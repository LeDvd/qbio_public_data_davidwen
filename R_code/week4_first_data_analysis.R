#Installs TCGAbiolinks if not already present
if (!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
if(!require(pacman)) install.packages("pacman")
library(pacman)

#######    Group 2: clinical   ###########
p_load(survival) #What do you need to do before running this line?
p_load(survminer) #What do you need to do before running this line?
p_load(arsenal) #What do you need to do before running this line?
p_load(dplyr)
clin_query <- GDCquery(project = "TCGA-BRCA", data.category="Clinical", file.type="xml")
# GDCdownload( clin_query ) #should only need this command once. This downloads the files onto your system.
clinic <- GDCprepare_clinic(clin_query, clinical.info="patient")
names(clinic)[names(clinic) == "days_to_last_followup"] = "days_to_last_follow_up"

#Add a new column to clinic called "age_category"
#If age_at_initial_pathologic_diagnosis is < 40, define patient as "Young" in new column
#If age_at_initial_pathologic_diagnosis is >= 60, define patient as "Old" in new column
#Other patients (between 40 and 60), define as "Mid"

age_clinical = clinic$age_at_initial_pathologic_diagnosis

clinic$age_category = ifelse(age_clinical < 40, "Young", ifelse(age_clinical >= 60, "Old", "Mid"))

# use tableone to create summary of clinic data
p_load(tableone)
clinic_summary = CreateTableOne(data=clinic)
subtypes <- TCGAquery_subtype(tumor = "BRCA")
age_subs = subtypes$age_at_initial_pathologic_diagnosis
subtypes$age_category = ifelse(age_subs < 40, "Young", ifelse(age_subs >= 60, "Old", "Mid"))

# Add the age_category column in the same way as above
# How many patients are in subtypes vs clinic? Why?

# BELOW CODE WILL NOT RUN. This is an example from another dataset
# Consider what you need to change for it to run on clinic
# Hint: age_category stays the same. "Oncotree.Code" was the name of a column
# Option to use either subtypes or clinic dataframe
names(subtypes) = make.names(names(subtypes))
table_arse <- tableby(age_category ~ (pathologic_stage)+ (mRNA.Clusters)+ (BRCA_Pathology),
	          data = subtypes, numeric.test="kwt", cat.test="fe", 
	          numeric.stats = c("Nmiss", "meansd"), total=FALSE)
df <- as.data.frame( summary(table_arse, text=TRUE, pfootnote=TRUE) 

# error: try increasing the size of the workspace and possibly 'mult'                                                                                                                                                                                             df <- as.data.frame( summary(table_arse, text=TRUE, pfootnote=TRUE) )
write.csv(df, “filepath/filename.csv”, row.names=FALSE)

# Modify below code for different variables of interest
# The "gender" can be changed to any of the column names of the clinic dataframe
# Look at the dataframe via the str(clinic) command
TCGAanalyze_survival( clinic, "gender", filename="survival_gender.pdf")
TCGAanalyze_survival( clinic, "age_category", filename="survival_age_category.pdf")
# use ls to confirm the file was created
# rsync to copy and view on local computer

overall_survival <- as.integer( ifelse( is.na(clinic$days_to_death), clinic$days_to_last_follow_up, clinic$days_to_death) )
clinic$overall_survival <- overall_survival
clinic$death_event <- ifelse(clinic$patient.vital_status == 'alive', 0,1)
colnames(clinic) #use if want to visually see affect of ^^
cox_fit <- coxph(Surv(overall_survival, death_event)~age_at_initial_pathologic_diagnosis, data=clinic)
jpeg("cox_plot_age_continuous.jpg")
ggadjustedcurves(fit, data=lung)
dev.off()
# Use rsync to copy figure onto local system and view
