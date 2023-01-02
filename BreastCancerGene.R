# Breast Cancer Gene Expression Profiles (METABRIC)

# Libraries


# Extracting Data from Sources

#Load Data

getwd()
setwd("C:/Users/erenj/OneDrive/Desktop/R-Code")

gene <- read.csv("./METABRIC_RNA_Mutation.csv")
str(gene)

gene <- read.table("./METABRIC_RNA_Mutation.csv",
                   sep = ",", header = TRUE)

View(gene)

head(gene)

#clinic_data = df.loc[:, df.columns[:31]]
#genomic_data = df.loc[:, df.columns[31:]]


# Detecting Missing Values
#Finding missing data and the percentage of it in each column

#total = genomic_data.isnull().sum().sort_values(ascending = False)
#percent = (genomic_data.isnull().sum() / genomic_data.isnull().count()).sort_values(ascending=False)
#missing_genomic = pd.concat([total, percent], axis = 1, keys = ['total_null', 'percent_null'])
#missing_genomic.head()

#Finding missing data and the percentage of it in each column

#total = clinic_data.isnull().sum().sort_values(ascending = False)
#percent = (clinic_data.isnull().sum() / clinic_data.isnull().count()).sort_values(ascending=False)
#missing_clinic = pd.concat([total, percent], axis = 1, keys = ['total_null', 'percent_null'])
#missing_clinic.head(14)


# Filling the Missing Values
# Visualize the Columns w/ Missing Values




# Drop Columns & Value Counts & Interpreting Data





# Label Encoding (Ordinal)





# One-Hot Encoding (Nominal)


# Add ordinal_df


# Add nominal_df


# Rounding Age of Patients


# Drop Some Columns
#df['cancer_type'].value_counts()





# Find Outliers

#new clinical data
#new_clinic_data = df[df.columns[:25]]
#new_clinic_data.head()



# Correlation Matrix between Numerical Clinical Data



# Update Data Set

#new_clinic_data.shape
#new_df = df.drop(new_clinic_data.columns, axis = 1, inplace = False)
#new_df = new_df.join(new_clinic_data)
#new_df.isnull().sum().sum()




# Treatment Types and Survivals
#treatments = ['chemotherapy', 'hormone_therapy', 'radio_therapy']

#died = new_df[new_df['overall_survival']==0]
#survived = new_df[new_df['overall_survival']==1]




# Statistical Summaries
# Statistical summary for categorical clinical attributes 
#new_df[new_df.columns].describe().T





# Genomic Data
#genomic_data = new_df.loc[:, new_df.dtypes == np.object]
#genomic_data.head()




# MinMax Scaler





# Choose Best PCA




# PCA Implementation




# Model





# Final Scores
#names = ['KNN', 'Logistic', 'Random Forest', 'Decision Tree', 'Extra Trees', 'AdaBoost']




# Machine-Learning-Verfahren
