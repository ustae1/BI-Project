# Load packages
library(corrplot)
library(VIM)
library(ggplot2)
library(gridExtra)
library(scales)
library(tidyr)
library(dplyr)
library(readr)
library(ggformula)
library(imputeMissings)
library(gapminder)
library(caret)
library(rpart) 
library(rattle)
library(kernlab)
library(randomForest)
library(DMwR2)
library(car)
library(MLeval)
###
library(dplyr)
library(tidyverse)
library(purrr)
library(tidyr)
library(magrittr)
library(Hmisc)
library(corrplot)
library(plotly)
library(ggplot2)
library(corrplot)  
library(ggcorrplot) 
library(broom)
library(ggthemes)
library(leaflet)
library(DT)
library(gridExtra)

getwd()
setwd("C:/Users/erenj/OneDrive/Desktop/R-Code")

gene <- read.csv("./METABRIC_RNA_Mutation.csv")

clinic_data <- gene[0:31]
genomic_data <- gene[32:693]

summary(clinic_data)

# Counting missing values overall
clinic_data %>%
  summarize(count = sum(is.na(clinic_data)))

# counting missing values in each of the columns
datatable(clinic_data %>% purrr::map(~sum(is.na(.))))


# Droping the whole row:
# Drop the missing values
clinic_data2 <- clinic_data %>% drop_na(neoplasm_histologic_grade, 
                                        mutation_count, 
                                        tumor_size)

# Counting missing values overall
clinic_data2 %>%
  dplyr::summarize(count = sum(is.na(clinic_data2)))

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


####
# Load data
df <- read.csv("./METABRIC_RNA_Mutation.csv", header=TRUE, sep=",", stringsAsFactors = FALSE);
dim(df)

df <- gene

## Preparing and Cleaning the Data

#To replace “ ” with  NA
df[df == ""] <- NA            
df_t = df[, c(1:31)] 

# to calculate percent of missing data for each variable
M_Val <- function(x){mm_v=round((sum(is.na(x))/length(x))*100,3)} 
nam_1<-colnames(df_t)
c<-apply(df_t,2,M_Val)

df_p <- data.frame(nam_1,c)
rownames(df_p)<-NULL
df_p1 <- df_p %>%  filter(c > 5 )  %>%  mutate(c1=ifelse(c < 5, c-.4, c-2.1))
gf_boxplot(c~nam_1 , data = df_p) %>%   
  gf_labs(title = "Missing data percentage",
          x = "",  
          y = "Percent (%)")%>%   gf_theme(axis.text.x=element_text(angle=65, hjust=1)) %>%
  gf_label(c1~nam_1, 
           label =~ c, 
           data = df_p1, col = "black", 
           fill =~ nam_1, alpha = .2,
           show.legend = FALSE)



####






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
  dplyr::summarize(count = sum(is.na(clinic_data4)))

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

View(no_genomic_data %>% purrr::map(~sum(is.na(.))))

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

# count of each number of tumor stage
qplot(data = clinic_data7, x = tumor_stage, geom = "bar")

# count of each number of overall survival
qplot(data = clinic_data7, x = overall_survival, geom = "bar")


# count of each number of overall survival
qplot(data = clinic_data7, x = cancer_type, geom = "bar")

# Droping the whole column:
clinic_data8 <- clinic_data7 %>% select(-cancer_type)


# count of each number of overall survival
qplot(data = clinic_data8, x = chemotherapy, geom = "bar")

# count of each number of overall survival
qplot(data = clinic_data8, x = neoplasm_histologic_grade, geom = "bar")

# count of each number of overall survival
qplot(data = clinic_data8, x = hormone_therapy, geom = "bar")

# count of each number of overall survival
qplot(data = clinic_data8, x = integrative_cluster, geom = "bar")

# count of each number of overall survival
qplot(data = clinic_data8, x = radio_therapy, geom = "bar")

# count of each number of overall survival
qplot(data = clinic_data8, x = death_from_cancer, geom = "bar")

# count of each number of overall survival
qplot(data = clinic_data8, x = type_of_breast_surgery, geom = "bar")



### Maybe splitting type of breast surgery in two factor variable 
### As breast conserving and Mastectomy
# ...

str(clinic_data)
clinic_data9 <- clinic_data8
clinic_data9$type_of_breast_surgery <- as.numeric(clinic_data9$type_of_breast_surgery)
clinic_data9$type_of_breast_surgery[clinic_data9$type_of_breast_surgery == '1'] <- '0'
clinic_data9$type_of_breast_surgery[clinic_data9$type_of_breast_surgery == '2'] <- '1'
clinic_data9$type_of_breast_surgery <- as.numeric(clinic_data9$type_of_breast_surgery)

clinic_data9$chemotherapy <- as.numeric(clinic_data9$chemotherapy)
clinic_data9$chemotherapy[clinic_data9$chemotherapy == '1'] <- '0'
clinic_data9$chemotherapy[clinic_data9$chemotherapy == '2'] <- '1'
clinic_data9$chemotherapy <- as.numeric(clinic_data9$chemotherapy)

clinic_data9$hormone_therapy <- as.numeric(clinic_data9$hormone_therapy)
clinic_data9$hormone_therapy[clinic_data9$hormone_therapy == '1'] <- '0'
clinic_data9$hormone_therapy[clinic_data9$hormone_therapy == '2'] <- '1'
clinic_data9$hormone_therapy <- as.numeric(clinic_data9$hormone_therapy)

clinic_data9$radio_therapy <- as.numeric(clinic_data9$radio_therapy)
clinic_data9$radio_therapy[clinic_data9$radio_therapy == '1'] <- '0'
clinic_data9$radio_therapy[clinic_data9$radio_therapy == '2'] <- '1'
clinic_data9$radio_therapy <- as.numeric(clinic_data9$radio_therapy)


data <- data.frame(clinic_data9$type_of_breast_surgery, 
                   clinic_data9$chemotherapy, 
                   clinic_data9$hormone_therapy,
                   clinic_data9$radio_therapy)

cor(data) 
corrplot(cor(data), method = "circle")    # Apply corrplot function


ggcorrplot(cor(data))                     # Apply ggcorrplot function

genomic_data2.cor = cor(genomic_data2)
genomic_data2.cor = cor(genomic_data2, method = c("spearman"))


genomic_data2.rcorr = rcorr(as.matrix(genomic_data2))
genomic_data2.rcorr

genomic_data2.coeff = genomic_data2.rcorr$r
genomic_data2.p = genomic_data2.rcorr$P

corrplot(genomic_data2.cor)

palette = colorRampPalette(c("green", "white", "red")) (20)

heatmap(x = genomic_data2.cor, col = palette, symm = TRUE)


ggplot(clinic_data9, aes(x = age_at_diagnosis, y = overall_survival_months, color = factor(tumor_stage), shape = factor(tumor_stage))) +
  geom_point(size=2) + 
  labs(x = "age_at_diagnosis", 
       y = "overall_survival_months", 
       color = "tumor_stage", 
       shape = "tumor_stage", 
       title = "overall_survival_months by age_at_diagnosis and tumor_stage",
       subtitle = "Data source: Breast Cancer Gene Expression Profiles (METABRIC)") +
  theme_dark() +
  scale_color_brewer(palette = "Set2")

ggplot(clinic_data9, aes(x = overall_survival_months, y = age_at_diagnosis, color = factor(tumor_stage), shape = factor(tumor_stage))) +
  geom_point(size=2) + 
  labs(x = "overall_survival_months", 
       y = "age_at_diagnosis", 
       color = "tumor_stage", 
       shape = "tumor_stage", 
       title = "age_at_diagnosis by overall_survival_months and tumor_stage",
       subtitle = "Data source: Breast Cancer Gene Expression Profiles (METABRIC)") +
  theme_dark() +
  scale_color_brewer(palette = "Set2")


# Droping the whole row:
# Drop the missing values, which were changed to "99"
clinic_data10 <- clinic_data2 %>% drop_na(tumor_stage)

ggplot(clinic_data10, aes(x = age_at_diagnosis, y = overall_survival_months, color = factor(tumor_stage), shape = factor(tumor_stage))) +
  geom_point(size=2) + 
  labs(x = "age_at_diagnosis", 
       y = "overall_survival_months", 
       color = "tumor_stage", 
       shape = "tumor_stage", 
       title = "overall_survival_months by age_at_diagnosis and tumor_stage",
       subtitle = "Data source: Breast Cancer Gene Expression Profiles (METABRIC)") +
  theme_dark() +
  scale_color_brewer(palette = "Set2")

ggplot(clinic_data10, aes(x = overall_survival_months, y = age_at_diagnosis, color = factor(tumor_stage), shape = factor(tumor_stage))) +
  geom_point(size=2) + 
  labs(x = "overall_survival_months", 
       y = "age_at_diagnosis", 
       color = "tumor_stage", 
       shape = "tumor_stage", 
       title = "age_at_diagnosis by overall_survival_months and tumor_stage",
       subtitle = "Data source: Breast Cancer Gene Expression Profiles (METABRIC)") +
  theme_dark() +
  scale_color_brewer(palette = "Set2")


# Convert cyl to factor
clinic_data10 <- clinic_data10 %>% 
  mutate(tumor_stage_factor = as.factor(tumor_stage))


# Plot grouped bar chart
ggplot(data = clinic_data10, 
       aes(x = tumor_stage, 
           fill = tumor_stage_factor)) + 
  geom_bar(position = "dodge")


# Stacked bar chart
ggplot(data = clinic_data10, 
       aes(x = " ", 
           fill = tumor_stage_factor)) + 
  geom_bar(position = "stack")


# Pie chart
ggplot(data = clinic_data10, 
       aes(x = " ", fill = tumor_stage_factor)) + 
  geom_bar(position = "stack") +
  coord_polar(theta = "y")

ggplot(data = clinic_data10, 
       aes(x = " ", fill = tumor_stage_factor)) + 
  geom_bar(position = "stack") +
  coord_polar(theta = "y") + 
  theme_minimal()

ggplot(data = clinic_data10, 
       aes(x = " ", fill = tumor_stage_factor)) + 
  geom_bar(position = "stack") +
  coord_polar(theta = "y") + 
  scale_fill_brewer(palette  = "Dark2")


qplot(numeric_data$tumor_size, numeric_data$mutation_count, data = numeric_data)
qplot(numeric_data$age_at_diagnosis, numeric_data$tumor_size, data = numeric_data)
qplot(numeric_data$overall_survival_months, numeric_data$age_at_diagnosis, data = numeric_data)


ggplot(numeric_data, aes(x = tumor_size, y = mutation_count)) + geom_point(shape = 1)

# numerical to categorical -> This created a new column, cyl_tumor_stag in the dataframe cinic_data
# cinic_data$cyl_tumor_stag <- factor(cinic_data$tumor_stage)

ggplot(clinic_data9, aes(x = age_at_diagnosis, y = mutation_count, shape = tumor_stage)) + geom_point()


ggplot(clinic_data9, aes(x = age_at_diagnosis, y = mutation_count, shape = tumor_stage)) + geom_point(color = "blue")


ggplot(clinic_data9, aes(x = age_at_diagnosis, y = mutation_count, color = tumor_stage)) + geom_point()


ggplot(clinic_data10, aes(x = age_at_diagnosis, y = mutation_count, color = tumor_stage_factor)) + geom_point()


# interactive boxplot with plotly
p <- ggplot(clinic_data10, aes(y = age_at_diagnosis)) + geom_boxplot()
ggplotly(p)

q <- ggplot(clinic_data10, aes(x = tumor_stage_factor, y = age_at_diagnosis)) + geom_boxplot()
ggplotly(q)

k <- ggplot(clinic_data10, aes(x = factor(tumor_stage), y = age_at_diagnosis)) + geom_boxplot()
ggplotly(k)

qplot(tumor_stage_factor, age_at_diagnosis, data = clinic_data10, geom = "boxplot")



# Faceting
ggplot(clinic_data10, aes(x = age_at_diagnosis)) + 
  geom_histogram() + 
  facet_wrap(.~tumor_stage)

ggplot(clinic_data10, aes(x = age_at_diagnosis)) + 
  geom_histogram() + 
  facet_wrap(.~cohort)

ggplot(clinic_data9, aes(x = factor(cohort))) + 
  geom_bar() + 
  facet_wrap(.~overall_survival) + 
  labs(x = "Cohort", title = "Overall survival by cohort (No is death, Yes is living)")


# Map
# Adding markers for map with caption
map <- leaflet() %>% addTiles() %>% 
  addMarkers(lng = 0.135524, lat = 52.176998, popup = 'Cambridge Research Institute')
map

map <- leaflet() %>% addTiles() %>% 
  addMarkers(lng = -123.119223, lat = 49.262668, popup = 'British Columbia Cancer Centre')
map



### For Shiny App -> Make Simple qplots with counting of factor variables
# ...
# ...

# ggplotly() !