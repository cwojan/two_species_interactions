####
#' Two Species Interaction Shiny App
#' EEB3407
#' Initial Code by George Furey
#' Revised for Spring 2023 by Chris Wojan
####

## Load required libraries
library(shiny)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(deSolve)

## Function to run differential equations of two species interaction
two_sp_growth <- function(times, y, parms) {
  with(as.list(c(y, parms)), {
    dN1 <- N1 * (r1 + (s11 * N1) + (s12 * N2)) #basis RSN model
    dN2 <- N2 * (r2 + (s22 * N2) + (s21 * N1)) #basis RSN model
    return(list(c(dN1, dN2)))
  })
}

## Setup the user interface
ui <- fluidPage(
  titlePanel("Exploring Two-Species Interactions"),
  ## The sidebar includes all the editable parameters
  sidebarLayout(
    sidebarPanel(
      sliderInput("r1", HTML(paste0("Species 1 Growth Rate (r", tags$sub('1'), ")")), value = -1, min = -1.5, max = 1.5, step = 0.1),
      sliderInput("r2", HTML(paste0("Species 2 Growth Rate (r", tags$sub('2'), ")")), value = 0.4, min = -1.5, max = 1.5, step = 0.1),
      sliderInput("s11", HTML(paste0("Effect of Species 1 on 1 (s", tags$sub('1,1'), ")")), value = -0.2, min = -1.5, max = 0, step = 0.1),
      sliderInput("s22", HTML(paste0("Effect of Species 2 on 2 (s", tags$sub('2,2'), ")")), value = -0.2, min = -1.5, max = 0, step = 0.1),
      sliderInput("s12", HTML(paste0("Effect of Species 2 on 1 (s", tags$sub('1,2'), ")")), value = 1.5, min = -1.5, max = 1.5, step = 0.1),
      sliderInput("s21", HTML(paste0("Effect of Species 1 on 2 (s", tags$sub('2,1'), ")")), value = -1.3, min = -1.5, max = 1.5, step = 0.1),
      sliderInput("N1", HTML(paste0("Initial Species 1 Abundance (N", tags$sub('1'), ")")), value = 2, min = 0, max = 5, step = 0.5),
      sliderInput("N2", HTML(paste0("Initial Species 2 Abundance (N", tags$sub('2'), ")")), value = 1, min = 0, max = 5, step = 0.5),
      radioButtons("time", "Max Timesteps", choices = c(2, 10, 100), selected = 100),
      withMathJax(helpText("Formula $$\\frac{1}{N_1}\\frac{dN_1}{dt} = r_1 + s_{1,1}N_1 + s_{1,2}N_2$$ \n
                           $$\\frac{1}{N_2}\\frac{dN_2}{dt} = r_2 + s_{2,2}N_1 + s_{2,1}N_1$$"))
          ),
    ## The main display shows tabs, one for each plot, and a comparison
    mainPanel(
      tabsetPanel(
       tabPanel("Pop. Size", plotOutput("pop_size")),
       tabPanel("Phase Plot", plotOutput("phase_plot", width = "600px", height = "600px")),
       tabPanel("Comparison",fluidRow(plotOutput("pop_size2"), plotOutput("phase_plot2")))
      )
    )
  )
)

server <- function(input, output) {
  
  ## Create a reactive conductor to run the differential equations
  ## This can be called by separate outputs without rerunning the code
  current_growth <- reactive({
    ## Run differential equations and save as long data
    ## i.e., species is a column
    long <- data.frame(ode(y = c(N1 = input$N1, N2 = input$N2), 
                          times = seq_len(as.numeric(input$time) * 10) / 10, 
                          func = two_sp_growth,
                          parms = c(
                            r1 = input$r1, r2 = input$r2,
                            s11 = input$s11, s22 = input$s22,
                            s12 = input$s12, s21 = input$s21
                          ),
                          method = "lsoda",
                          maxsteps = 20000
    )) %>%
      pivot_longer(cols = !time,
                   names_to = "species",
                   values_to = "N") %>%
      group_by(species) %>%
      mutate(dN = c(diff(N),NA), ## calculations for growth rates
             per_capita = dN / N)
    
    ## Reformat data as wide for isocline plotting
    wide <- long %>%
      pivot_wider(names_from = species, values_from = c("N", "dN", "per_capita")) %>%
      rename(N1 = N_N1,
             N2 = N_N2)
    
    ## return a list of both data organizations
    return(list(long = long, wide = wide))
  })
  
  ## Create a reactive conductor to store isocline slopes and intercepts
  ## (currently only called by the phase plot output)
  current_iso <- reactive({
    tibble(N1_slope = -input$s11 / input$s12,
           N2_slope = -input$s21 / input$s22,
           N1_intercept = -input$r1 / input$s12,
           N2_intercept = -input$r2 / input$s22,
           N2_xintercept = -input$r2 / input$s21)
  })
  
  ## Plot population over time
  output$pop_size <- renderPlot({
    ggplot(current_growth()$long, aes(x = time, y = N, color = species)) +
      ## draw growth line
      geom_line(aes(linetype = species), linewidth = 2) +
      ## label axes
      labs(x = "Time", y = "Population Size") +
      ## specify legends
      scale_color_manual(name = "Species", values = c("#d95f02", "#7570b3")) +
      scale_linetype_discrete(name = "Species") +
      coord_cartesian(ylim = c(0, max(current_growth()$long$N))) +
      theme_bw() +
      theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16),
            legend.text = element_text(size = 14), legend.title = element_text(size = 16))
  })
  
  ## Population plot with timestamps
  output$pop_size2 <- renderPlot({
    ggplot(current_growth()$long, aes(x = time, y = N, color = species)) +
      ## draw growth line
      geom_line(aes(linetype = species), alpha = 0.3, linewidth = 1) +
      ## print timestamps
      geom_text(data = current_growth()$long %>% 
                  filter(case_when(time <= 12 ~ (time%%1) == 0, 
                                   time <= 24 & time > 12 ~ (time%%2) == 0,
                                   time >= 24 ~ (time%%4) == 0)), 
                aes(label = time)) +
      ## label axes
      labs(x = "Time", y = "Population Size") +
      ## specify legends
      scale_color_manual(name = "Species", values = c("#d95f02", "#7570b3")) +
      scale_linetype_discrete(name = "Species") +
      coord_cartesian(ylim = c(0, max(current_growth()$long$N))) +
      theme_bw() +
      theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16),
            legend.text = element_text(size = 14), legend.title = element_text(size = 16))
  })
  
  ## Plot isoclines on a phase plot
  output$phase_plot <- renderPlot({
    ggplot(current_growth()$wide, aes(x = N1, y = N2)) +
      ## draw isoclines
      with(current_iso(), geom_abline(slope = N1_slope, intercept = N1_intercept, color = "#d95f02", linewidth = 3)) +
      with(current_iso(),
           if(N2_slope %in% c(Inf, -Inf)){ ## if s22 is 0...
             geom_vline(xintercept = N2_xintercept, color = "#7570b3", linewidth = 3)
           } else {
             geom_abline(slope = N2_slope, intercept = N2_intercept, color = "#7570b3", linewidth = 3)
           }
      ) +
      ## draw invisible lines to force a color legend for isoclines
      geom_line(data = current_growth()$long, aes(x = N, y = N, color = species), alpha = 0, linewidth = 3) +
      ## draw points over time
      geom_point(aes(fill = time), shape = 21, stroke = 0, size = 2) +
      ## set scales, axes, etc.
      scale_fill_viridis_c(name = "Timestep") +
      scale_color_manual(name = "Species Isocline",
                         values = c("N1" = "#d95f02", "N2" = "#7570b3"),
                         guide = guide_legend(override.aes = list(alpha = 1))) +
      labs(x = "N1", y = "N2") +
      coord_fixed(xlim = c(0, max(c(max(current_growth()$wide$N1), max(current_growth()$wide$N2)))),
                  ylim = c(0, max(c(max(current_growth()$wide$N1), max(current_growth()$wide$N2))))) +
      theme_bw() +
      theme(legend.position = "bottom", axis.text = element_text(size = 14), axis.title = element_text(size = 16),
            legend.text = element_text(size = 14), legend.title = element_text(size = 16))
  })
  
  ## Isocline plot with timestamps
  output$phase_plot2 <- renderPlot({
    ggplot(current_growth()$wide, 
           aes(x = N1, y = N2)) +
      ## draw isoclines
      with(current_iso(), geom_abline(slope = N1_slope, intercept = N1_intercept, color = "#d95f02", linewidth = 3)) +
      with(current_iso(),
           if(N2_slope %in% c(Inf, -Inf)){ ## if s22 is 0...
             geom_vline(xintercept = N2_xintercept, color = "#7570b3", linewidth = 3)
           } else {
             geom_abline(slope = N2_slope, intercept = N2_intercept, color = "#7570b3", linewidth = 3)
           }
      ) +
      ## draw invisible lines to force color legend for isoclines
      geom_line(data = current_growth()$long, aes(x = N, y = N, color = species), alpha = 0, linewidth = 3) +
      ## draw points over time with timestamps
      geom_point(aes(fill = time), shape = 21, stroke = 0, size = 2, alpha = 0.2) +
      geom_text(data = current_growth()$wide %>% 
                  filter(case_when(time <= 12 ~ (time%%1) == 0, 
                                   time <= 24 & time > 12 ~ (time%%2) == 0,
                                   time >= 24 ~ (time%%4) == 0)),
                aes(label = time)) +
      ## set scales, axes, etc.
      scale_color_manual(name = "Species Isocline",
                         values = c("N1" = "#d95f02", "N2" = "#7570b3"),
                         guide = guide_legend(override.aes = list(alpha = 1))) +
      scale_fill_viridis_c(name = "Timestep") +
      labs(x = "N1", y = "N2") +
      coord_fixed(xlim = c(0, max(c(max(current_growth()$wide$N1), max(current_growth()$wide$N2)))),
                  ylim = c(0, max(c(max(current_growth()$wide$N1), max(current_growth()$wide$N2))))) +
      theme_bw() +
      theme(legend.position = "bottom", axis.text = element_text(size = 14), axis.title = element_text(size = 16),
            legend.text = element_text(size = 14), legend.title = element_text(size = 16))
  })

}

shinyApp(ui, server)


