# Load necessary libraries
library(shiny)
library(MASS)
library(ggplot2)

# Define UI for the application
ui <- fluidPage(
  titlePanel("Clinical Trial Design Simulator"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("mean_bp", "Mean Baseline Blood Pressure:", 150),
      numericInput("mean_chol", "Mean Cholesterol:", 200),
      numericInput("mean_age", "Mean Age:", 50),
      
      numericInput("sd_bp", "SD Baseline Blood Pressure:", 12),
      numericInput("sd_chol", "SD Cholesterol:", 20),
      numericInput("sd_age", "SD Age:", 10),
      
      numericInput("corr_bp_chol", "Correlation BP & Cholesterol:", 0.5),
      numericInput("corr_bp_age", "Correlation BP & Age:", 0.3),
      numericInput("corr_chol_age", "Correlation Cholesterol & Age:", 0.2),
      
      numericInput("n", "Number of Participants per Group:", 100),
      
      actionButton("simulate", "Run Simulation")
    ),
    
    mainPanel(
      verbatimTextOutput("summary"),
      verbatimTextOutput("pvalue"),
      plotOutput("plot")
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  simulate_data <- eventReactive(input$simulate, {
    # Define the means vector
    means <- c(input$mean_bp, input$mean_chol, input$mean_age)
    
    # Define the standard deviations vector
    sds <- c(input$sd_bp, input$sd_chol, input$sd_age)
    
    # Define the correlation matrix
    corr_matrix <- matrix(c(
      1, input$corr_bp_chol, input$corr_bp_age,
      input$corr_bp_chol, 1, input$corr_chol_age,
      input$corr_bp_age, input$corr_chol_age, 1
    ), nrow = 3)
    
    # Create the covariance matrix
    cov_matrix <- diag(sds) %*% corr_matrix %*% diag(sds)
    
    # Simulate the data
    group_A_data <- mvrnorm(input$n, mu = means, Sigma = cov_matrix)
    group_B_data <- mvrnorm(input$n, mu = means, Sigma = cov_matrix)
    
    # Combine into a data frame
    data <- data.frame(
      Group = rep(c("A", "B"), each = input$n),
      Baseline_BP = c(group_A_data[, 1], group_B_data[, 1]),
      Cholesterol = c(group_A_data[, 2], group_B_data[, 2]),
      Age = c(group_A_data[, 3], group_B_data[, 3])
    )
    
    # Simulate treatment effects
    effect_A <- rnorm(input$n, mean = -15, sd = 5)
    effect_B <- rnorm(input$n, mean = -10, sd = 5)
    
    # Calculate final blood pressure
    data$Final_BP <- c(data$Baseline_BP[1:input$n] + effect_A,
                       data$Baseline_BP[(input$n+1):(2*input$n)] + effect_B)
    
    return(data)
  })
  
  output$summary <- renderPrint({
    data <- simulate_data()
    summary(data)
  })
  
  output$pvalue <- renderPrint({
    data <- simulate_data()
    model <- lm(Final_BP ~ Group + Baseline_BP + Cholesterol + Age, data = data)
    summary_model <- summary(model)
    
    # Extracting the p-value for the treatment effect (GroupB)
    p_value <- summary_model$coefficients["GroupB", "Pr(>|t|)"]
    
    cat("P-value for the effect of Treatment Group on Final Blood Pressure:", p_value, "\n")
  })
  
  output$plot <- renderPlot({
    data <- simulate_data()
    ggplot(data, aes(x = Group, y = Final_BP, fill = Group)) +
      geom_boxplot() +
      labs(title = "Final Systolic Blood Pressure by Treatment Group",
           x = "Treatment Group",
           y = "Systolic Blood Pressure (mmHg)") +
      theme_minimal()
  })
}

# Run the application
shinyApp(ui = ui, server = server)
