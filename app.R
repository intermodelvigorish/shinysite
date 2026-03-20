imv.binary<-function(y, p1, p2, sigma=1e-4) {
    p1<-ifelse(p1<sigma,sigma,p1)
    p2<-ifelse(p2<sigma,sigma,p2)
    p1<-ifelse(p1>1-sigma,1-sigma,p1)
    p2<-ifelse(p2>1-sigma,1-sigma,p2)
    
    ll<-function(x,p) {
        z<-log(p)*x+log(1-p)*(1-x)
        z<-sum(z)/length(x)
        exp(z)
    }    
    loglik1<-ll(y,p1)
    loglik2<-ll(y,p2)
    getcoins<-function(a) {
        f<-function(p,a) abs(p*log(p)+(1-p)*log(1-p)-log(a))
        nlminb(.5,f,lower=0.001,upper=.999,a=a)$par
    }
    c1<-getcoins(loglik1)
    c2<-getcoins(loglik2)
    ew<-function(p1,p0) (p1-p0)/p0
    imv<-ew(c2,c1)
    imv
}

library(shiny)

est <- function(x, th, fixc = FALSE) {
  ll <- function(pars, x, th) {
    a <- pars[1]
    b <- pars[2]
    if (length(pars) == 3) c <- pars[3] else c <- 0
    p <- c + (1 - c) / (1 + exp(-a * (th - b)))
    p <- pmax(pmin(p, 1 - 1e-10), 1e-10)
    -1 * sum(x * log(p) + (1 - x) * log(1 - p))
  }
  if (fixc) pars <- c(1, 0) else pars <- c(1, 0, .1)
  optim(pars, ll, x = x, th = th)$par
}

ui <- fluidPage(
  titlePanel("IRT OM Sensitivity Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Tab 1: Sensitivity to c"),
      sliderInput("n_persons_tab1", "Sample Size:", 
                  min = 1000, max = 10000, value = 5000, step = 1000),
      actionButton("simulate_tab1", "Generate New Data", class = "btn-primary"),
      hr(),
      h4("Tab 2: Sensitivity to Sample Size"),
      sliderInput("c_param_tab2", "Guessing Parameter (c):", 
                  min = 0, max = 0.5, value = 0.2, step = 0.01),
      actionButton("simulate_tab2", "Generate New Data", class = "btn-primary"),
      hr(),
      p("Fixed parameters: a=1.2, b=-0.5")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Sensitivity to c",
                 h3("OM as a Function of c Parameter"),
                 p("100 random c values sampled from Uniform(0, 0.5)"),
                 plotOutput("sensitivity_c_plot", height = "500px"),
                 hr(),
                 h4("Summary Statistics"),
                 verbatimTextOutput("summary_c_output")
        ),
        
        tabPanel("Sensitivity to Sample Size",
                 h3("OM as a Function of Sample Size"),
                 p("Fixed c parameter, varying sample size"),
                 plotOutput("sensitivity_n_plot", height = "500px"),
                 hr(),
                 h4("Summary Statistics"),
                 verbatimTextOutput("summary_n_output")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  a_param <- 1.2
  b_param <- -0.5
  n_samples_c <- 100
  
  sim_data <- reactiveValues(sensitivity_c = NULL, sensitivity_n = NULL)
  
  generate_data_c <- function() {
    withProgress(message = 'Generating c sensitivity data...', value = 0, {
      th <- rnorm(input$n_persons_tab1, mean = 0, sd = 1)
      p_fixed <- 0.2 + (1 - 0.2) / (1 + exp(-a_param * (th - b_param)))
      responses2 <- rbinom(input$n_persons_tab1, 1, p_fixed)
      
      set.seed(NULL)
      c_values <- sort(runif(n_samples_c, 0, 0.5))
      sens_results <- data.frame(c = numeric(), om = numeric(), 
                                 omu = numeric(), omstd = numeric())
      
      for (i in seq_along(c_values)) {
        incProgress(1/n_samples_c, detail = paste("Sample", i, "of", n_samples_c))
        
        c_val <- c_values[i]
        p_temp <- c_val + (1 - c_val) / (1 + exp(-a_param * (th - b_param)))
        resp_temp <- rbinom(input$n_persons_tab1, 1, p_temp)
        
        est3_temp <- est(th = th, x = resp_temp)
        est2_temp <- est(th = th, x = resp_temp, fixc = TRUE)
        
        p3_temp <- est3_temp[3] + (1 - est3_temp[3]) / 
          (1 + exp(-est3_temp[1] * (th - est3_temp[2])))
        p2_temp <- 1 / (1 + exp(-est2_temp[1] * (th - est2_temp[2])))
        
        omv_temp <- imv.binary(responses2, p2_temp, p3_temp)
        
        sens_results <- rbind(sens_results, data.frame(
          c = c_val, om = omv_temp[1], omu = omv_temp[2], omstd = omv_temp[3]
        ))
      }
      sim_data$sensitivity_c <- sens_results
    })
  }
  
  generate_data_n <- function() {
    withProgress(message = 'Generating sample size sensitivity data...', value = 0, {
      n_values <- seq(500, 10000, by = 500)
      sens_results <- data.frame(n = numeric(), om = numeric(), 
                                 omu = numeric(), omstd = numeric())
      
      for (i in seq_along(n_values)) {
        incProgress(1/length(n_values), detail = paste("n =", n_values[i]))
        
        n_val <- n_values[i]
        th <- rnorm(n_val, mean = 0, sd = 1)
        
        p <- input$c_param_tab2 + (1 - input$c_param_tab2) / 
          (1 + exp(-a_param * (th - b_param)))
        responses <- rbinom(n_val, 1, p)
        responses2 <- rbinom(n_val, 1, p)
        
        est3 <- est(th = th, x = responses)
        est2 <- est(th = th, x = responses, fixc = TRUE)
        
        p3 <- est3[3] + (1 - est3[3]) / (1 + exp(-est3[1] * (th - est3[2])))
        p2 <- 1 / (1 + exp(-est2[1] * (th - est2[2])))
        
        omv <- imv.binary(responses2, p2, p3)
        
        sens_results <- rbind(sens_results, data.frame(
          n = n_val, om = omv[1], omu = omv[2], omstd = omv[3]
        ))
      }
      sim_data$sensitivity_n <- sens_results
    })
  }
  
  observe({ 
    generate_data_c() 
    generate_data_n()
  })
  
  observeEvent(input$simulate_tab1, { generate_data_c() })
  observeEvent(input$simulate_tab2, { generate_data_n() })
  
  output$summary_c_output <- renderPrint({
    req(sim_data$sensitivity_c)
    cat("Summary of OM Statistics across 100 c values:\n")
    cat("=============================================\n\n")
    for(stat in c("om", "omu", "omstd")) {
      vals <- sim_data$sensitivity_c[[stat]]
      cat(sprintf("%s:\n  Range: [%.4f, %.4f]\n  Mean: %.4f\n  SD: %.4f\n\n", 
                  stat, min(vals), max(vals), mean(vals), sd(vals)))
    }
  })
  
  output$sensitivity_c_plot <- renderPlot({
    req(sim_data$sensitivity_c)
    cols <- c("red", "blue", "green")
    
    par(mar = c(5, 5, 4, 2))
    plot(sim_data$sensitivity_c$c, sim_data$sensitivity_c$om, 
         type = "n", ylim = c(-0.003, 0.003),
         xlab = "Guessing Parameter (c)", ylab = "OM Value",
         main = paste("Sensitivity to c (n =", input$n_persons_tab1, ")"),
         cex.lab = 1.2)
    grid(col = "lightgray")
    
    for(i in 1:3) {
      stat <- c("om", "omu", "omstd")[i]
      y <- sim_data$sensitivity_c[[stat]]
      points(sim_data$sensitivity_c$c, y, col = adjustcolor(cols[i], 0.4), pch = 16)
      lines(lowess(sim_data$sensitivity_c$c, y, f = 0.3), col = cols[i], lwd = 2)
    }
    legend("topleft", legend = c("om", "omu", "omstd"), 
           col = cols, lwd = 2, pch = 16, bg = "white")
  })
  
  output$summary_n_output <- renderPrint({
    req(sim_data$sensitivity_n)
    cat("Summary of OM Statistics across sample sizes:\n")
    cat("=============================================\n")
    cat(sprintf("c parameter: %.2f\n\n", input$c_param_tab2))
    for(stat in c("om", "omu", "omstd")) {
      vals <- sim_data$sensitivity_n[[stat]]
      cat(sprintf("%s:\n  Range: [%.4f, %.4f]\n  Mean: %.4f\n  SD: %.4f\n\n", 
                  stat, min(vals), max(vals), mean(vals), sd(vals)))
    }
  })
  
  output$sensitivity_n_plot <- renderPlot({
    req(sim_data$sensitivity_n)
    cols <- c("red", "blue", "green")
    
    par(mar = c(5, 5, 4, 2))
    plot(sim_data$sensitivity_n$n, sim_data$sensitivity_n$om, 
         type = "n", ylim = c(-0.003, 0.003),
         xlab = "Sample Size", ylab = "OM Value",
         main = paste("Sensitivity to Sample Size (c =", input$c_param_tab2, ")"),
         cex.lab = 1.2)
    grid(col = "lightgray")
    
    for(i in 1:3) {
      stat <- c("om", "omu", "omstd")[i]
      y <- sim_data$sensitivity_n[[stat]]
      points(sim_data$sensitivity_n$n, y, col = cols[i], pch = 16, cex = 1.5)
      lines(lowess(sim_data$sensitivity_n$n, y, f = 0.3), col = cols[i], lwd = 2)
    }
    legend("topleft", legend = c("om", "omu", "omstd"), 
           col = cols, lwd = 2, pch = 16, bg = "white")
  })
}

shinyApp(ui = ui, server = server)
