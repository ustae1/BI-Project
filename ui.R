library(dplyr)
library(tidyverse)
library(purrr)
library(tidyr)
library(magrittr)
library(Hmisc)
library(corrplot)
library(plotly)

getwd()
setwd("C:/Users/erenj/OneDrive/Desktop/kMeansApp")

gene <- read.csv("./METABRIC_RNA_Mutation.csv")

clinic_data <- gene[0:31]
genomic_data <- gene[32:693]

# Counting missing values overall
clinic_data %>%
  summarize(count = sum(is.na(clinic_data)))

# counting missing values in each of the columns
clinic_data %>% purrr::map(~sum(is.na(.)))

# Droping the whole row:
# Drop the missing values
clinic_data2 <- clinic_data %>% drop_na(neoplasm_histologic_grade, 
                                         mutation_count, 
                                         tumor_size)

# Counting missing values overall
clinic_data2 %>%
  summarize(count = sum(is.na(clinic_data2)))

# counting missing values in each of the columns
clinic_data2 %>% purrr::map(~sum(is.na(.)))

# Replace the missing values in tumor_stage
clinic_data3 <- clinic_data2 %>% replace_na(list(tumor_stage = 99))

# Counting missing values overall
clinic_data3 %>%
  summarize(count = sum(is.na(clinic_data3)))


# list the data types for each column
clinic_data3 %>% 
  summarize_all(class) %>% 
  gather(variable, class)


# Datenbereinigung - Part I
clinic_data3[clinic_data3 == ''] <- NA

# counting missing values in each of the columns
clinic_data3 %>% purrr::map(~sum(is.na(.)))

# Droping the whole row:
# Drop the missing values
clinic_data4 <- clinic_data3 %>% drop_na(type_of_breast_surgery, 
                                        cancer_type_detailed, 
                                        cellularity,
                                        er_status_measured_by_ihc,
                                        tumor_other_histologic_subtype,
                                        oncotree_code,
                                        death_from_cancer)

# Counting missing values overall
clinic_data4 %>%
  summarize(count = sum(is.na(clinic_data4)))

# counting missing values in each of the columns
clinic_data4 %>% purrr::map(~sum(is.na(.)))


clinic_data4$primary_tumor_laterality
clinic_data4$X3.gene_classifier_subtype

# Datenbereinigung - Part II
na_cols <- is.na(clinic_data4)
clinic_data5 <- replace(clinic_data4, na_cols, "Unknown")


# Counting missing values overall
clinic_data5 %>%
  summarize(count = sum(is.na(clinic_data5)))

# Datentyp von chr zu factor
clinic_data6 <- clinic_data5 %>%
  mutate_if(sapply(clinic_data5, is.character), as.factor)

# restliche Datentypen 
clinic_data6$chemotherapy <- factor(clinic_data6$chemotherapy, levels = c(0, 1), labels = c("No", "Yes"))
clinic_data6$hormone_therapy <- factor(clinic_data6$hormone_therapy, levels = c(0, 1), labels = c("No", "Yes"))
clinic_data6$overall_survival <- factor(clinic_data6$overall_survival, levels = c(0, 1), labels = c("No", "Yes"))
clinic_data6$radio_therapy <- factor(clinic_data6$radio_therapy, levels = c(0, 1), labels = c("No", "Yes"))
clinic_data6$cohort <- factor(clinic_data6$cohort)
clinic_data6$neoplasm_histologic_grade <- factor(clinic_data6$neoplasm_histologic_grade)
clinic_data6$tumor_stage <- factor(clinic_data6$tumor_stage)


clinic_data6$lymph_nodes_examined_positive
clinic_data6$mutation_count
clinic_data6$nottingham_prognostic_index

# Droping the whole column:
clinic_data7 <- clinic_data6 %>% select(-patient_id)

str(clinic_data7)
sum(is.na(clinic_data7))

# Correlation matrix
# Correlations between multiple variables
correlation <- clinic_data7 %>%
  select(age_at_diagnosis, lymph_nodes_examined_positive, 
         mutation_count, nottingham_prognostic_index, 
         overall_survival_months, tumor_size) %>% 
  cor(use = "pairwise.complete.obs", method = "pearson")

correlation

numeric_data <- clinic_data7 %>% 
  select(age_at_diagnosis,
         mutation_count,
         overall_survival_months, 
         tumor_size)

numeric_data_cor <- cor(numeric_data,
                        method = "pearson", 
                        use='pairwise.complete.obs')

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

corrplot(numeric_data_cor, method = "color", col = col(200),  
         type = "upper", order = "hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
)


# list the data types for each column
genomic_data %>% 
  summarize_all(class) %>% 
  gather(variable, class)

genomic_data2 <- genomic_data[1:489]

# list the data types for each column
genomic_data2 %>% 
  summarize_all(class) %>% 
  gather(variable, class)

no_genomic_data <- genomic_data[490:662]

# list the data types for each column
no_genomic_data %>% 
  summarize_all(class) %>% 
  gather(variable, class)

mutation_data <- no_genomic_data

mutation_data

# counting missing values in each of the columns
mutation_data %>% purrr::map(~sum(is.na(.)))

# Datentyp von chr zu factor
mutation_data2 <- mutation_data %>%
  mutate_if(sapply(mutation_data, is.character), as.factor)

str(mutation_data2)



# ggplotly() !



# k-means only works with numerical variables,
# so don't give the user the option to select
# a categorical variable
vars <- setdiff(names(iris), "Species")

pageWithSidebar(
  headerPanel('Iris k-means clustering'),
  sidebarPanel(
    selectInput('xcol', 'X Variable', vars),
    selectInput('ycol', 'Y Variable', vars, selected = vars[[2]]),
    numericInput('clusters', 'Cluster count', 3, min = 1, max = 9)
  ),
  mainPanel(
    plotOutput('plot1')
  )
)
