# Load necessary libraries
library(shiny)
library(MASS)
library(ggplot2)
library(shinydashboard)
library(DT)
library(plotly)
library(pwr)

# Define UI for the application
ui <- dashboardPage(
  dashboardHeader(title = "Clinical Trial Design Simulator"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Simulation", tabName = "simulation", icon = icon("dashboard")),
      menuItem("Results", tabName = "results", icon = icon("chart-bar")),
      menuItem("Data", tabName = "data", icon = icon("table")),
      menuItem("Power Analysis", tabName = "power", icon = icon("bolt"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "simulation",
              fluidRow(
                box(
                  title = "Baseline Parameters", status = "primary", solidHeader = TRUE,
                  numericInput("mean_bp", "Mean Baseline Blood Pressure:", 150, min = 0, max = 300),
                  numericInput("mean_chol", "Mean Cholesterol:", 200, min = 0, max = 500),
                  numericInput("mean_age", "Mean Age:", 50, min = 0, max = 120),
                  numericInput("sd_bp", "SD Baseline Blood Pressure:", 12, min = 0),
                  numericInput("sd_chol", "SD Cholesterol:", 20, min = 0),
                  numericInput("sd_age", "SD Age:", 10, min = 0)
                ),
                box(
                  title = "Correlations", status = "primary", solidHeader = TRUE,
                  sliderInput("corr_bp_chol", "Correlation BP & Cholesterol:", -1, 1, 0.5, step = 0.1),
                  sliderInput("corr_bp_age", "Correlation BP & Age:", -1, 1, 0.3, step = 0.1),
                  sliderInput("corr_chol_age", "Correlation Cholesterol & Age:", -1, 1, 0.2, step = 0.1)
                )
              ),
              fluidRow(
                box(
                  title = "Treatment Effects", status = "warning", solidHeader = TRUE,
                  numericInput("effect_A_mean", "Group A Effect Mean:", -15),
                  numericInput("effect_A_sd", "Group A Effect SD:", 5, min = 0),
                  numericInput("effect_B_mean", "Group B Effect Mean:", -10),
                  numericInput("effect_B_sd", "Group B Effect SD:", 5, min = 0)
                ),
                box(
                  title = "Simulation Settings", status = "success", solidHeader = TRUE,
                  numericInput("n", "Number of Participants per Group:", 100, min = 10),
                  numericInput("sim_count", "Number of Simulations:", 1, min = 1, max = 1000),
                  actionButton("simulate", "Run Simulation", class = "btn-lg btn-primary")
                )
              )
      ),
      tabItem(tabName = "results",
              fluidRow(
                box(
                  title = "Summary Statistics", status = "info", solidHeader = TRUE,
                  verbatimTextOutput("summary")
                ),
                box(
                  title = "Statistical Analysis", status = "info", solidHeader = TRUE,
                  verbatimTextOutput("pvalue")
                )
              ),
              fluidRow(
                box(
                  title = "Blood Pressure Distribution", status = "primary", solidHeader = TRUE,
                  plotlyOutput("bp_dist_plot")
                ),
                box(
                  title = "Treatment Effect", status = "primary", solidHeader = TRUE,
                  plotlyOutput("treatment_effect_plot")
                )
              ),
              fluidRow(
                box(
                  title = "Correlation Plot", status = "primary", solidHeader = TRUE, width = 12,
                  plotlyOutput("correlation_plot")
                )
              )
      ),
      tabItem(tabName = "data",
              box(
                title = "Simulated Data", status = "primary", solidHeader = TRUE, width = 12,
                DTOutput("data_table")
              )
      ),
      tabItem(tabName = "power",
              fluidRow(
                box(
                  title = "Power Analysis Parameters", status = "primary", solidHeader = TRUE,
                  numericInput("effect_size", "Expected Effect Size (Cohen's d):", 0.5, min = 0, max = 2, step = 0.1),
                  numericInput("sig_level", "Significance Level (alpha):", 0.05, min = 0.01, max = 0.1, step = 0.01),
                  numericInput("desired_power", "Desired Power:", 0.8, min = 0.5, max = 0.99, step = 0.01),
                  actionButton("run_power", "Run Power Analysis", class = "btn-lg btn-primary")
                ),
                box(
                  title = "Power Analysis Results", status = "info", solidHeader = TRUE,
                  verbatimTextOutput("power_results")
                )
              ),
              fluidRow(
                box(
                  title = "Power Curve", status = "primary", solidHeader = TRUE, width = 12,
                  plotlyOutput("power_curve")
                )
              )
      )
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  simulate_data <- eventReactive(input$simulate, {
    withProgress(message = 'Simulating data', value = 0, {
      sim_results <- lapply(1:input$sim_count, function(i) {
        incProgress(1/input$sim_count, detail = paste("Simulation", i))
        
        means <- c(input$mean_bp, input$mean_chol, input$mean_age)
        sds <- c(input$sd_bp, input$sd_chol, input$sd_age)
        corr_matrix <- matrix(c(
          1, input$corr_bp_chol, input$corr_bp_age,
          input$corr_bp_chol, 1, input$corr_chol_age,
          input$corr_bp_age, input$corr_chol_age, 1
        ), nrow = 3)
        cov_matrix <- diag(sds) %*% corr_matrix %*% diag(sds)
        
        group_A_data <- mvrnorm(input$n, mu = means, Sigma = cov_matrix)
        group_B_data <- mvrnorm(input$n, mu = means, Sigma = cov_matrix)
        
        data <- data.frame(
          Simulation = i,
          Group = rep(c("A", "B"), each = input$n),
          Baseline_BP = c(group_A_data[, 1], group_B_data[, 1]),
          Cholesterol = c(group_A_data[, 2], group_B_data[, 2]),
          Age = c(group_A_data[, 3], group_B_data[, 3])
        )
        
        effect_A <- rnorm(input$n, mean = input$effect_A_mean, sd = input$effect_A_sd)
        effect_B <- rnorm(input$n, mean = input$effect_B_mean, sd = input$effect_B_sd)
        
        data$Final_BP <- c(data$Baseline_BP[1:input$n] + effect_A,
                           data$Baseline_BP[(input$n+1):(2*input$n)] + effect_B)
        
        return(data)
      })
      
      do.call(rbind, sim_results)
    })
  })
  
  output$summary <- renderPrint({
    data <- simulate_data()
    summary(data)
  })
  
  output$pvalue <- renderPrint({
    data <- simulate_data()
    model <- lm(Final_BP ~ Group + Baseline_BP + Cholesterol + Age, data = data)
    summary_model <- summary(model)
    
    p_value <- summary_model$coefficients["GroupB", "Pr(>|t|)"]
    
    cat("P-value for the effect of Treatment Group on Final Blood Pressure:", p_value, "\n")
    cat("Coefficient estimate for Group B:", summary_model$coefficients["GroupB", "Estimate"], "\n")
  })
  
  output$bp_dist_plot <- renderPlotly({
    data <- simulate_data()
    p <- ggplot(data, aes(x = Final_BP, fill = Group)) +
      geom_density(alpha = 0.7) +
      labs(title = "Distribution of Final Blood Pressure by Group",
           x = "Final Blood Pressure (mmHg)",
           y = "Density") +
      theme_minimal()
    ggplotly(p)
  })
  
  output$treatment_effect_plot <- renderPlotly({
    data <- simulate_data()
    data$BP_Change <- data$Final_BP - data$Baseline_BP
    p <- ggplot(data, aes(x = Group, y = BP_Change, fill = Group)) +
      geom_boxplot() +
      labs(title = "Blood Pressure Change by Treatment Group",
           x = "Treatment Group",
           y = "Change in Blood Pressure (mmHg)") +
      theme_minimal()
    ggplotly(p)
  })
  
  output$correlation_plot <- renderPlotly({
    data <- simulate_data()
    p <- plot_ly(data, x = ~Baseline_BP, y = ~Cholesterol, z = ~Age, color = ~Group,
                 type = "scatter3d", mode = "markers") %>%
      layout(scene = list(xaxis = list(title = "Baseline BP"),
                          yaxis = list(title = "Cholesterol"),
                          zaxis = list(title = "Age")))
    p
  })
  
  output$data_table <- renderDT({
    datatable(simulate_data(), options = list(pageLength = 10))
  })
  
  # Power Analysis
  observeEvent(input$run_power, {
    output$power_results <- renderPrint({
      power_result <- pwr.t.test(d = input$effect_size, 
                                 sig.level = input$sig_level, 
                                 power = input$desired_power, 
                                 type = "two.sample", 
                                 alternative = "two.sided")
      
      cat("Required sample size per group:", ceiling(power_result$n), "\n")
      cat("Total sample size:", ceiling(power_result$n) * 2, "\n")
      cat("Actual power achieved:", power_result$power, "\n")
    })
    
    output$power_curve <- renderPlotly({
      n_range <- seq(10, 500, by = 10)
      power_values <- sapply(n_range, function(n) {
        pwr.t.test(n = n, 
                   d = input$effect_size, 
                   sig.level = input$sig_level, 
                   type = "two.sample", 
                   alternative = "two.sided")$power
      })
      
      data <- data.frame(n = n_range, power = power_values)
      
      p <- ggplot(data, aes(x = n, y = power)) +
        geom_line() +
        geom_hline(yintercept = input$desired_power, linetype = "dashed", color = "red") +
        geom_vline(xintercept = ceiling(pwr.t.test(d = input$effect_size, 
                                                   sig.level = input$sig_level, 
                                                   power = input$desired_power, 
                                                   type = "two.sample", 
                                                   alternative = "two.sided")$n), 
                   linetype = "dashed", color = "blue") +
        labs(title = "Power Curve",
             x = "Sample Size per Group",
             y = "Power") +
        theme_minimal()
      
      ggplotly(p)
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server)