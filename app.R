#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
bc <- read.csv("C:\\Users\\minhm\\Documents\\RWorkSpace\\BI-Project\\METABRIC_RNA_Mutation.csv", stringsAsFactors = FALSE)

# Datenbereinigung
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

str(ds_clean)
sum(is.na(ds_clean))


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Old Faithful Geyser Data"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("bins",
                        "Number of bins:",
                        min = 1,
                        max = 50,
                        value = 30)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white',
             xlab = 'Waiting time to next eruption (in mins)',
             main = 'Histogram of waiting times')
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
