#' Interactive Timeline Plot Shiny Module UI
#'
#' @param id Character string, namespace identifier for the module
#' @return Shiny UI elements for the timeline plot module
#' @export
#' @importFrom shiny NS fluidRow column selectInput checkboxInput numericInput
#'   dateRangeInput actionButton plotOutput verbatimTextOutput h3 h4 hr
#' @importFrom DT DTOutput
#' @examples
#' \dontrun{
#' # In your Shiny UI
#' timeline_plot_UI("timeline_module")
#' }
timeline_plot_UI <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::fluidPage(
    shiny::titlePanel("MSK-CHORD Interactive Timeline Explorer"),
    
    shiny::fluidRow(
      # Sidebar with controls
      shiny::column(
        width = 4,
        shiny::wellPanel(
          shiny::h4("Plot Configuration"),
          
          # Biomarker selection
          shiny::selectInput(
            ns("biomarker"),
            "Select Biomarker:",
            choices = c("PSA" = "PSA", 
                        "CEA" = "CEA", 
                        "CA 15-3" = "CA_15-3", 
                        "CA 19-9" = "CA_19-9", 
                        "Gleason Score" = "GLEASON"),
            selected = "PSA"
          ),
          
          # Cancer type filter
          shiny::selectInput(
            ns("cancer_type"),
            "Cancer Type:",
            choices = c("All Cancer Types" = "all",
                        "Breast Cancer" = "Breast Cancer",
                        "Colorectal Cancer" = "Colorectal Cancer", 
                        "Non-Small Cell Lung Cancer" = "Non-Small Cell Lung Cancer",
                        "Pancreatic Cancer" = "Pancreatic Cancer",
                        "Prostate Cancer" = "Prostate Cancer"),
            selected = "all"
          ),
          
          shiny::hr(),
          
          # Patient selection method
          shiny::h4("Patient Selection"),
          shiny::radioButtons(
            ns("patient_selection_method"),
            "Selection Method:",
            choices = c("Random Sample" = "random",
                        "Specific Patients" = "specific",
                        "Filter by Criteria" = "filter"),
            selected = "random"
          ),
          
          # Conditional panels for different selection methods
          shiny::conditionalPanel(
            condition = paste0("input['", ns("patient_selection_method"), "'] == 'random'"),
            shiny::numericInput(
              ns("n_patients"),
              "Number of Patients:",
              value = 20,
              min = 1,
              max = 100,
              step = 1
            )
          ),
          
          shiny::conditionalPanel(
            condition = paste0("input['", ns("patient_selection_method"), "'] == 'specific'"),
            shiny::selectizeInput(
              ns("specific_patients"),
              "Search and Select Patients:",
              choices = NULL,
              multiple = TRUE,
              options = list(
                placeholder = "Type patient ID to search...",
                maxItems = 50,
                searchLimit = 100
              )
            )
          ),
          
          shiny::conditionalPanel(
            condition = paste0("input['", ns("patient_selection_method"), "'] == 'filter'"),
            shiny::selectInput(
              ns("filter_gender"),
              "Gender:",
              choices = c("All" = "all", "Female" = "Female", "Male" = "Male"),
              selected = "all"
            ),
            shiny::selectInput(
              ns("filter_stage"),
              "Stage:",
              choices = c("All" = "all"),  # Will be populated dynamically
              selected = "all"
            ),
            shiny::selectInput(
              ns("filter_treatment"),
              "Has Treatment History:",
              choices = c("All" = "all", "Yes" = "yes", "No" = "no"),
              selected = "all"
            ),
            shiny::numericInput(
              ns("filter_age_min"),
              "Minimum Age:",
              value = NA,
              min = 0,
              max = 100
            ),
            shiny::numericInput(
              ns("filter_age_max"),
              "Maximum Age:",
              value = NA,
              min = 0,
              max = 100
            )
          ),
          
          shiny::hr(),
          
          # Plot options
          shiny::h4("Plot Options"),
          shiny::checkboxInput(
            ns("log_scale"),
            "Use Log Scale (for lab values)",
            value = TRUE
          ),
          
          shiny::checkboxInput(
            ns("add_events"),
            "Show Clinical Events",
            value = TRUE
          ),
          
          shiny::selectInput(
            ns("facet_by"),
            "Facet By:",
            choices = c("None" = "none"),
            selected = "none"
          ),
          
          # Date range
          shiny::checkboxInput(
            ns("use_date_range"),
            "Limit Date Range",
            value = FALSE
          ),
          
          shiny::conditionalPanel(
            condition = paste0("input['", ns("use_date_range"), "'] == true"),
            shiny::numericInput(
              ns("date_min"),
              "Start Date (days):",
              value = -1000
            ),
            shiny::numericInput(
              ns("date_max"),
              "End Date (days):",
              value = 1000
            )
          ),
          
          shiny::hr(),
          shiny::actionButton(
            ns("update_plot"),
            "Update Plot",
            class = "btn-primary",
            style = "width: 100%;"
          )
        )
      ),
      
      # Main panel with plot and info
      shiny::column(
        width = 8,
        shiny::tabsetPanel(
          id = ns("main_tabs"),
          
          shiny::tabPanel(
            "Timeline Plot",
            shiny::br(),
            shiny::plotOutput(
              ns("timeline_plot"),
              height = "600px"
            ),
            shiny::br(),
            shiny::h4("Plot Information"),
            shiny::verbatimTextOutput(ns("plot_info"))
          ),
          
          shiny::tabPanel(
            "Interactive Plot",
            shiny::br(),
            plotly::plotlyOutput(
              ns("interactive_plot"),
              height = "600px"
            )
          ),
          
          shiny::tabPanel(
            "Data Summary",
            shiny::br(),
            shiny::h4("Selected Patients Summary"),
            DT::DTOutput(ns("patient_summary")),
            shiny::br(),
            shiny::h4("Biomarker Data Summary"),
            shiny::verbatimTextOutput(ns("data_summary"))
          ),
          
          shiny::tabPanel(
            "Clinical Events",
            shiny::br(),
            shiny::h4("Clinical Events for Selected Patients"),
            DT::DTOutput(ns("events_table"))
          )
        )
      )
    )
  )
}


#' Interactive Timeline Plot Shiny Module Server
#'
#' @param id Character string, namespace identifier for the module
#' @param msk_data Reactive expression containing the MSK-CHORD data
#' @return Server logic for the timeline plot module
#' @export
#' @importFrom shiny moduleServer reactive observe observeEvent updateSelectInput
#'   updateSelectizeInput req renderPlot renderText
#' @importFrom DT renderDT datatable
#' @importFrom plotly renderPlotly ggplotly
#' @examples
#' \dontrun{
#' # In your Shiny server
#' timeline_plot_Server("timeline_module", reactive(msk_data))
#' }
timeline_plot_Server <- function(id, msk_data) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # Reactive values for storing data
    values <- shiny::reactiveValues(
      available_patients = NULL,
      clinical_data = NULL,
      current_plot = NULL,
      selected_patients = NULL
    )
    
    # Initialize data when msk_data is available
    shiny::observe({
      shiny::req(msk_data())
      
      # Get clinical data
      values$clinical_data <- as.data.frame(MultiAssayExperiment::colData(msk_data()))
      
      # Get unique patients
      values$available_patients <- unique(values$clinical_data$PATIENT_ID)
      
      # Update patient selection choices
      shiny::updateSelectizeInput(
        session,
        "specific_patients",
        choices = values$available_patients,
        server = TRUE
      )
      
      # Update stage filter choices
      stages <- unique(values$clinical_data$STAGE_HIGHEST_RECORDED)
      stages <- stages[!is.na(stages) & stages != ""]
      stage_choices <- c("All" = "all")
      if (length(stages) > 0) {
        names(stages) <- stages
        stage_choices <- c(stage_choices, stages)
      }
      
      shiny::updateSelectInput(
        session,
        "filter_stage",
        choices = stage_choices
      )
      
      # Update faceting options based on biomarker
      update_facet_choices()
    })
    
    # Update faceting options when biomarker changes
    shiny::observe({
      update_facet_choices()
    })
    
    update_facet_choices <- function() {
      facet_choices <- c("None" = "none")
      
      if (input$biomarker == "GLEASON") {
        facet_choices <- c(facet_choices, "Gleason Category" = "gleason_category")
      }
      
      facet_choices <- c(facet_choices, 
                         "Gender" = "GENDER",
                         "Cancer Type" = "CANCER_TYPE",
                         "Stage" = "STAGE_HIGHEST_RECORDED")
      
      shiny::updateSelectInput(
        session,
        "facet_by",
        choices = facet_choices
      )
    }
    
    # Get filtered patients based on selection method
    filtered_patients <- shiny::reactive({
      shiny::req(values$clinical_data)
      
      clinical_data <- values$clinical_data
      
      # Filter by cancer type if not "all"
      if (input$cancer_type != "all") {
        clinical_data <- clinical_data[clinical_data$CANCER_TYPE == input$cancer_type, ]
      }
      
      if (input$patient_selection_method == "random") {
        available_patients <- unique(clinical_data$PATIENT_ID)
        n_to_sample <- min(input$n_patients, length(available_patients))
        return(sample(available_patients, n_to_sample))
        
      } else if (input$patient_selection_method == "specific") {
        return(input$specific_patients)
        
      } else if (input$patient_selection_method == "filter") {
        # Apply filters
        filtered_data <- clinical_data
        
        # Gender filter
        if (input$filter_gender != "all") {
          filtered_data <- filtered_data[filtered_data$GENDER == input$filter_gender, ]
        }
        
        # Stage filter
        if (input$filter_stage != "all") {
          filtered_data <- filtered_data[filtered_data$STAGE_HIGHEST_RECORDED == input$filter_stage, ]
        }
        
        # Age filters
        if (!is.na(input$filter_age_min)) {
          filtered_data <- filtered_data[filtered_data$CURRENT_AGE_DEID >= input$filter_age_min, ]
        }
        if (!is.na(input$filter_age_max)) {
          filtered_data <- filtered_data[filtered_data$CURRENT_AGE_DEID <= input$filter_age_max, ]
        }
        
        # Treatment filter
        if (input$filter_treatment != "all") {
          metadata_list <- MultiAssayExperiment::metadata(msk_data())
          if ("timeline_treatment" %in% names(metadata_list)) {
            treatment_data <- as.data.frame(metadata_list[["timeline_treatment"]])
            patients_with_treatment <- unique(treatment_data$PATIENT_ID)
            
            if (input$filter_treatment == "yes") {
              filtered_data <- filtered_data[filtered_data$PATIENT_ID %in% patients_with_treatment, ]
            } else {
              filtered_data <- filtered_data[!filtered_data$PATIENT_ID %in% patients_with_treatment, ]
            }
          }
        }
        
        return(unique(filtered_data$PATIENT_ID))
      }
    })
    
    # Generate plot when update button is clicked
    shiny::observeEvent(input$update_plot, {
      shiny::req(msk_data(), filtered_patients())
      
      patients <- filtered_patients()
      
      if (length(patients) == 0) {
        values$current_plot <- NULL
        return()
      }
      
      values$selected_patients <- patients
      
      # Prepare plot parameters
      cancer_type_param <- if (input$cancer_type == "all") NULL else input$cancer_type
      facet_param <- if (input$facet_by == "none") NULL else input$facet_by
      date_range_param <- if (input$use_date_range) c(input$date_min, input$date_max) else NULL
      
      # Generate the plot
      tryCatch({
        values$current_plot <- plot_biomarker_timeline(
          msk_data = msk_data(),
          biomarker = input$biomarker,
          cancer_type = cancer_type_param,
          patient_ids = patients,
          log_scale = input$log_scale,
          facet_by = facet_param,
          add_events = input$add_events,
          date_range = date_range_param
        )
      }, error = function(e) {
        values$current_plot <- NULL
        showNotification(paste("Error generating plot:", e$message), type = "error")
      })
    })
    
    # Render the main plot
    output$timeline_plot <- shiny::renderPlot({
      if (is.null(values$current_plot)) {
        plot.new()
        text(0.5, 0.5, "Click 'Update Plot' to generate timeline visualization", 
             cex = 1.2, font = 2)
      } else {
        values$current_plot
      }
    })
    
    # Render interactive plot
    output$interactive_plot <- plotly::renderPlotly({
      if (is.null(values$current_plot)) {
        return(NULL)
      }
      
      tryCatch({
        plotly::ggplotly(values$current_plot, tooltip = c("x", "y", "colour"))
      }, error = function(e) {
        return(NULL)
      })
    })
    
    # Render plot information
    output$plot_info <- shiny::renderText({
      if (is.null(values$selected_patients)) {
        return("No patients selected. Configure settings and click 'Update Plot'.")
      }
      
      info <- paste(
        "Biomarker:", input$biomarker,
        "\nNumber of Patients:", length(values$selected_patients),
        "\nCancer Type:", if (input$cancer_type == "all") "All Types" else input$cancer_type,
        "\nLog Scale:", input$log_scale,
        "\nClinical Events:", input$add_events,
        "\nFaceting:", if (input$facet_by == "none") "None" else input$facet_by
      )
      
      if (input$use_date_range) {
        info <- paste(info, "\nDate Range:", input$date_min, "to", input$date_max, "days")
      }
      
      return(info)
    })
    
    # Render patient summary table
    output$patient_summary <- DT::renderDT({
      shiny::req(values$selected_patients, values$clinical_data)
      
      # Get clinical data for selected patients
      patient_data <- values$clinical_data[values$clinical_data$PATIENT_ID %in% values$selected_patients, ]
      
      # Deduplicate by patient (keep first sample per patient)
      patient_data <- patient_data[!duplicated(patient_data$PATIENT_ID), ]
      
      # Select relevant columns
      summary_cols <- c("PATIENT_ID", "CANCER_TYPE", "GENDER", "CURRENT_AGE_DEID", 
                        "STAGE_HIGHEST_RECORDED", "OS_STATUS", "OS_MONTHS")
      available_cols <- intersect(summary_cols, colnames(patient_data))
      
      display_data <- patient_data[, available_cols, drop = FALSE]
      
      DT::datatable(
        display_data,
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          autoWidth = TRUE
        ),
        rownames = FALSE
      )
    })
    
    # Render data summary
    output$data_summary <- shiny::renderText({
      shiny::req(msk_data(), input$biomarker)
      
      # Get biomarker summary
      biomarker_info <- get_timeline_biomarkers(msk_data())
      biomarker_row <- biomarker_info[biomarker_info$biomarker == input$biomarker, ]
      
      if (nrow(biomarker_row) > 0) {
        summary_text <- paste(
          "Biomarker:", biomarker_row$biomarker,
          "\nTotal Measurements:", biomarker_row$n_measurements,
          "\nTotal Patients with Data:", biomarker_row$n_patients,
          "\nDate Range:", biomarker_row$date_range, "days",
          "\nValue Range:", biomarker_row$value_range
        )
      } else {
        summary_text <- "No data available for selected biomarker."
      }
      
      return(summary_text)
    })
    
    # Render clinical events table
    output$events_table <- DT::renderDT({
      shiny::req(values$selected_patients, msk_data())
      
      events <- get_clinical_events(msk_data(), values$selected_patients)
      
      if (nrow(events) > 0) {
        DT::datatable(
          events,
          options = list(
            pageLength = 15,
            scrollX = TRUE,
            autoWidth = TRUE,
            order = list(list(1, 'asc'))  # Sort by START_DATE
          ),
          rownames = FALSE
        )
      } else {
        DT::datatable(
          data.frame(Message = "No clinical events found for selected patients"),
          options = list(pageLength = 5),
          rownames = FALSE
        )
      }
    })
  })
}


#' Launch Shiny Timeline Explorer App
#'
#' @param msk_data A MultiAssayExperiment object from download_data()
#' @param port Integer specifying the port to run the app on (default: random)
#' @param launch.browser Logical indicating whether to launch the app in browser
#' @return Shiny app object
#' @export
#' @importFrom shiny shinyApp fluidPage
#' @examples
#' \dontrun{
#' msk_data <- download_data()
#' launch_timeline_explorer(msk_data)
#' }
launch_timeline_explorer <- function(msk_data, port = NULL, launch.browser = TRUE) {
  
  # Check required packages
  required_packages <- c("shiny", "DT", "plotly")
  missing_packages <- character()
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_packages <- c(missing_packages, pkg)
    }
  }
  
  if (length(missing_packages) > 0) {
    stop("Required packages missing: ", paste(missing_packages, collapse = ", "), 
         "\nInstall with: install.packages(c('", paste(missing_packages, collapse = "', '"), "'))")
  }
  
  ui <- shiny::fluidPage(
    # Add custom CSS for better styling
    shiny::tags$head(
      shiny::tags$style(shiny::HTML("
        .content-wrapper, .right-side {
          background-color: #f4f4f4;
        }
        .well {
          background-color: white;
          border: 1px solid #ddd;
        }
        .btn-primary {
          background-color: #337ab7;
          border-color: #2e6da4;
        }
      "))
    ),
    
    timeline_plot_UI("timeline_explorer")
  )
  
  server <- function(input, output, session) {
    timeline_plot_Server("timeline_explorer", shiny::reactive(msk_data))
  }
  
  app <- shiny::shinyApp(ui = ui, server = server)
  
  if (is.null(port)) {
    shiny::runApp(app, launch.browser = launch.browser)
  } else {
    shiny::runApp(app, port = port, launch.browser = launch.browser)
  }
}