#
# inst/shiny/lsspca_gui/app.R

options(shiny.maxRequestSize = 200 * 1024^2)

library(shiny)
library(bslib)
library(DT)

# ADD read scales
# ADD plots
# ADD variable selection choice

ui <- page_navbar(
  title = "PMlsspca – LSSPCA GUI",
  theme = bs_theme(version = 5),
  nav_panel(
    "Data",
    layout_sidebar(
      sidebar = sidebar(
        fileInput("file", "Upload CSV", accept = ".csv"),
        textInput("sep", "Separator", value = ","),
        checkboxInput("header", "Header", TRUE),
        checkboxInput("stringsAsFactors", "Strings as factors", FALSE),
        hr(),
        uiOutput("vars_ui"),
        hr(),
        checkboxInput("center", "Center", TRUE),
        checkboxInput("scale", "Scale to unit variance", FALSE),
        width = 360
      ),
      card(card_header("Preview"), DTOutput("preview"))
    )
  ),
  nav_panel(
    "Diagnostics",
    layout_sidebar(
      sidebar = sidebar(
        numericInput("nplot", "No. eigenvalues to plot", value = 20, min = 1, step = 1),
        checkboxInput("corr_mat", "Treat as correlation matrix (trace = p)", TRUE),
        numericInput("nfit_line", "Fit line using last k points (0 = none, -r excludes, r includes)", value = 0, min = 0, step = 1),
        width = 360
      ),
      layout_column_wrap(
        width = 1/2,
        card(card_header("Scree plot"), plotOutput("scree", height = 420)),
        card(card_header("Wachter QQ-plot"), plotOutput("wachter", height = 420))
      )
    )
  ),
  nav_panel(
    "Model",
    layout_sidebar(
      sidebar = sidebar(
        selectInput("variant", "Variant", choices = c("pSPCA", "uSPCA", "cSPCA"), selected = "pSPCA"),
        numericInput("ncomp", "Components", value = 4, min = 1, step = 1),
        numericInput("alpha", "Target recovered variance (alpha)", value = 0.95, min = 0.50, max = 0.999, step = 0.01),
        actionButton("run", "Run", class = "btn-primary"),
        br(), br(),
        verbatimTextOutput("run_msg"),
        width = 360
      ),
      card(card_header("Status"), verbatimTextOutput("status"))
    )
  ),
  nav_panel(
    "Results",
    layout_sidebar(
      sidebar = sidebar(
        checkboxInput("plotcontrib", "Plot contributions as % (if supported)", TRUE),
        downloadButton("dl_rds", "Download fit (.rds)"),
        width = 360
      ),
      layout_column_wrap(
        width = 1/2,
        card(card_header("Summary"), DTOutput("sumtbl")),
        card(card_header("Loadings / Contributions plot"), plotOutput("fitplot", height = 420))
      )
    )
  )
)

server <- function(input, output, session) {
  
  # ---------- Data ----------
  dat <- reactive({
    req(input$file)
    sep <- input$sep
    if (!nzchar(sep)) sep <- ","
    read.csv(
      input$file$datapath,
      sep = sep,
      header = isTRUE(input$header),
      stringsAsFactors = isTRUE(input$stringsAsFactors),
      check.names = FALSE
    )
  })
  
  output$vars_ui <- renderUI({
    df <- dat()
    is_num <- vapply(df, function(z) is.numeric(z) || is.integer(z), logical(1))
    num_cols <- names(df)[is_num]
    if (!length(num_cols)) {
      return(tags$div(class = "text-danger", "No numeric columns detected."))
    }
    selectInput("vars", "Variables", choices = num_cols, selected = num_cols, multiple = TRUE)
  })
  
  Xmat <- reactive({
    df <- dat()
    req(input$vars)
    X <- as.matrix(df[, input$vars, drop = FALSE])
    storage.mode(X) <- "double"
    if (isTRUE(input$center)) X <- scale(X, center = TRUE, scale = FALSE)
    if (isTRUE(input$scale))  X <- scale(X, center = FALSE, scale = TRUE)
    X
  })
  
  output$preview <- renderDT({
    df <- dat()
    DT::datatable(head(df, 50), options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # ---------- Diagnostics (eigenvalues) ----------
  eigvals <- reactive({
    X <- Xmat()
    # eigenvalues of sample covariance/correlation as appropriate
    S <- stats::cov(X)
    ev <- sort(eigen(S, symmetric = TRUE, only.values = TRUE)$values, decreasing = TRUE)
    ev
  })
  
  output$scree <- renderPlot({
    req(eigvals())
    if (!requireNamespace("PMlsspca", quietly = TRUE)) {
      plot.new(); text(0.5, 0.5, "Package PMlsspca not available.")
      return()
    }
    ev <- eigvals()
    nplot <- min(length(ev), as.integer(input$nplot))
    pl <- PMlsspca::screeplot(ev, nplot = nplot, perc = TRUE, prn = FALSE, rtn = TRUE)
    print(pl)
  })
  
  output$wachter <- renderPlot({
    req(eigvals())
    if (!requireNamespace("PMlsspca", quietly = TRUE)) {
      plot.new(); text(0.5, 0.5, "Package PMlsspca not available.")
      return()
    }
    ev <- eigvals()
    p <- ncol(Xmat()); n <- nrow(Xmat())
    nplot <- min(length(ev), as.integer(input$nplot))
    nfl <- as.integer(input$nfit_line)
    nfl <- if (is.na(nfl) || nfl <= 0) NULL else nfl
    
    pl <- PMlsspca::wachterqq(
      eigvals = ev, p = p, n = n, gamma = n/p,
      cor = isTRUE(input$corr_mat),
      nplot = nplot,
      nfit_line = nfl,
      prn = FALSE, rtn = TRUE
    )
    print(pl)
  })
  
  # ---------- Fit ----------
  fit <- reactiveVal(NULL)
  run_err <- reactiveVal(NULL)
  
  observeEvent(input$run, {
    run_err(NULL)
    output$status <- renderText("Running...")
    
    if (!requireNamespace("PMlsspca", quietly = TRUE)) {
      run_err("Package 'PMlsspca' not available.")
      fit(NULL)
      output$status <- renderText("ERROR")
      return()
    }
    
    X <- Xmat()
    
    # ===== EDIT THIS BLOCK ONLY (your computation entry-point) =====
    # Replace PMlsspca::lsspca(...) with your actual function/signature.
    obj <- tryCatch({
      PMlsspca::lsspca(
        X = X,
        ncomp = as.integer(input$ncomp),
        alpha = input$alpha,
        variant = input$variant
      )
    }, error = function(e) e)
    # ===============================================================
    
    if (inherits(obj, "error")) {
      run_err(conditionMessage(obj))
      fit(NULL)
      output$status <- renderText(paste("ERROR:", conditionMessage(obj)))
    } else {
      fit(obj)
      output$status <- renderText("Done.")
    }
  })
  
  output$run_msg <- renderText({
    if (!is.null(run_err())) run_err() else ""
  })
  
  # ---------- Results ----------
  output$sumtbl <- renderDT({
    req(fit())
    s <- tryCatch(summary(fit()), error = function(e) NULL)
    validate(need(!is.null(s), "summary(fit) failed."))
    if (is.matrix(s) || is.data.frame(s)) {
      DT::datatable(as.data.frame(s), options = list(pageLength = 10, scrollX = TRUE))
    } else if (is.list(s)) {
      df <- data.frame(
        Item = names(s),
        Value = vapply(s, function(z) {
          if (is.atomic(z) && length(z) == 1) as.character(z) else paste(class(z), collapse = "/")
        }, character(1)),
        row.names = NULL
      )
      DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE))
    } else {
      df <- data.frame(Value = as.character(s))
      DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE))
    }
  })
  
  output$fitplot <- renderPlot({
    req(fit())
    obj <- fit()
    ok <- TRUE
    tryCatch({
      # if your plot.spca supports plotcontributions, it will use it; otherwise ignore
      plot(obj, plotcontributions = isTRUE(input$plotcontrib))
    }, error = function(e) {
      ok <<- FALSE
      plot.new()
      text(0.5, 0.5, paste("Plot failed:", conditionMessage(e)))
    })
    invisible(ok)
  })
  
  output$dl_rds <- downloadHandler(
    filename = function() paste0("PMlsspca_fit_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds"),
    content = function(file) saveRDS(fit(), file = file)
  )
}

shinyApp(ui, server)

# Run the application 
shinyApp(ui = ui, server = server)
