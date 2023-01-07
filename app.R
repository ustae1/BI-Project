#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
#library(tidyverse)
library(ggplot2)
library(dplyr)
library(plotly)
library(ggforce)
bc <- read.csv("C:\\Users\\minhm\\Documents\\RWorkSpace\\BI-Project\\METABRIC_RNA_Mutation.csv", stringsAsFactors = FALSE)

# Datenaufbereitung
ds <- bc[0:31]
summary(ds)
ds[ds == ''] <- NA
na_cols <- is.na(ds)
ds <- replace(ds, na_cols, "Unknown")

ds_clean <- ds %>%
  mutate_if(sapply(ds, is.character), as.factor)

ds_clean$chemotherapy <- factor(ds_clean$chemotherapy, levels = c(0, 1), labels = c("No", "Yes"))
ds_clean$hormone_therapy <- factor(ds_clean$hormone_therapy, levels = c(0, 1), labels = c("No", "Yes"))
ds_clean$overall_survival <- factor(ds_clean$overall_survival, levels = c(0, 1), labels = c("No", "Yes"))
ds_clean$radio_therapy <- factor(ds_clean$radio_therapy, levels = c(0, 1), labels = c("No", "Yes"))
ds_clean$cohort <- factor(ds_clean$cohort)
ds_clean$tumor_size <- as.numeric(ds_clean$tumor_size)
ds_clean$mutation_count <- as.numeric(ds_clean$mutation_count)
ds_clean$age_at_diagnosis <- round(ds_clean$age_at_diagnosis)


ds_survivalMonths <- ds_clean %>%
  group_by(age_at_diagnosis, cohort) %>%
  summarise_at(vars(overall_survival_months), list(average_survival_months = mean))

ds_survival <- ds_clean %>% group_by(chemotherapy, radio_therapy, hormone_therapy) %>% count(overall_survival)


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Breast Cancer"),

    # Sidebar with a select input for number of cohort 
    sidebarLayout(
        sidebarPanel(
            
            selectInput(inputId = "cohort1", label = "Choose Cohort", choices = c("All",
                                                                                "Cohort 1",
                                                                                "Cohort 2",
                                                                                "Cohort 3",
                                                                                "Cohort 4",
                                                                                "Cohort 5"),
                        multiple = F),
            
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotlyOutput("distPlot1"),
           plotOutput("distPlot2")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot1 <- renderPlotly({
        if(input$cohort1 == "All"){
          p <- ggplotly(ggplot(ds_survivalMonths, aes(x=age_at_diagnosis, y=average_survival_months)) + geom_line())
        }else if(input$cohort1 == "Cohort 1"){
          p <- ggplotly(ggplot(data = filter(ds_survivalMonths, cohort == "1"), aes(x=age_at_diagnosis, y=average_survival_months)) + geom_line())
        }else if(input$cohort1 == "Cohort 2"){
          p <- ggplotly(ggplot(data = filter(ds_survivalMonths, cohort == "2"), aes(x=age_at_diagnosis, y=average_survival_months)) + geom_line())
        }else if(input$cohort1 == "Cohort 3"){
          p <- ggplotly(ggplot(data = filter(ds_survivalMonths, cohort == "3"), aes(x=age_at_diagnosis, y=average_survival_months)) + geom_line())
        }else if(input$cohort1 == "Cohort 4"){
          p <- ggplotly(ggplot(data = filter(ds_survivalMonths, cohort == "4"), aes(x=age_at_diagnosis, y=average_survival_months)) + geom_line())
        }else if(input$cohort1 == "Cohort 5"){
          p <- ggplotly(ggplot(data = filter(ds_survivalMonths, cohort == "5"), aes(x=age_at_diagnosis, y=average_survival_months)) + geom_line())
        }
        return(p)
    })
    
    output$distPlot2 <- renderPlot({
      data <- data.frame(x = c(0, 1, -1),
                         y = c(-0.5, 1, 1),
                         tx = c(0, 1.5, -1.5),
                         ty = c(-1, 1.3, 1.3),
                         cat = c("Chemo \n 17|28", "Hormone \n 135|270", "Radio \n 119|109"))
      
      ggplot(data, aes(x0 = x, y0 = y, r = 1.5, fill = cat)) +
        geom_circle(alpha = 0.25, size = 1, color = "black", show.legend = F) +
        geom_text(aes(x = tx, y = ty, label = cat)) +
        annotate(geom="text", x=0, y=1.5, label="253|333") + 
        annotate(geom="text", x=-0.9, y=0, label="75|93") + 
        annotate(geom="text", x=0.9, y=0, label="9|19") + 
        annotate(geom="text", x=0, y=0.5, label="83|72") + 
        theme_void()
    })
    
    
}
print(head(p))

# Run the application 
shinyApp(ui = ui, server = server)
