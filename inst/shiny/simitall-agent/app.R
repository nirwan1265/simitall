library(shiny)
library(simitall)

ui <- fluidPage(
  tags$head(
    tags$title("simitall agent"),
    tags$style(HTML("\
      :root { --ink:#18251f; --muted:#617168; --paper:#f7f3e9; --card:#fffdf8; --line:#d8ded4; --forest:#0d5f4c; --forest-dark:#084435; --mint:#dff0e7; }
      body { background:var(--paper); color:var(--ink); font-family:Avenir Next, Avenir, Helvetica Neue, sans-serif; }
      .container-fluid { max-width:1200px; padding:32px 24px 54px; }
      .agent-header { align-items:flex-start; display:flex; gap:24px; justify-content:space-between; margin:0 auto 28px; max-width:980px; }
      .agent-mark { color:var(--forest); font-size:12px; font-weight:800; letter-spacing:.16em; margin-bottom:8px; text-transform:uppercase; }
      h1 { font-family:Iowan Old Style, Palatino, Georgia, serif; font-size:clamp(34px, 5vw, 54px); font-weight:700; letter-spacing:-.04em; line-height:1; margin:0 0 12px; }
      .agent-subtitle { color:var(--muted); font-size:17px; line-height:1.5; margin:0; max-width:680px; }
      .privacy-note { background:var(--mint); border-radius:999px; color:var(--forest-dark); font-size:12px; font-weight:700; padding:8px 12px; white-space:nowrap; }
      .chat-shell { margin:0 auto; max-width:980px; }
      .message { border:1px solid var(--line); border-radius:18px; box-shadow:0 10px 24px rgba(27,52,38,.06); margin-bottom:18px; overflow:hidden; }
      .message-label { color:var(--muted); font-size:11px; font-weight:800; letter-spacing:.12em; padding:14px 20px 0; text-transform:uppercase; }
      .assistant-message { background:var(--card); }
      .assistant-answer { background:var(--card); border:0; color:var(--ink); font-family:ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace; font-size:13px; line-height:1.58; margin:0; max-height:620px; overflow:auto; padding:14px 20px 20px; white-space:pre-wrap; word-break:break-word; }
      .empty-answer { color:var(--muted); font-size:16px; line-height:1.55; padding:14px 20px 20px; }
      .composer { background:var(--card); border:1px solid var(--line); border-radius:18px; box-shadow:0 10px 24px rgba(27,52,38,.06); padding:15px; }
      .composer .form-group { margin:0; }
      .composer-row { align-items:center; display:flex; gap:10px; }
      .composer-row .form-group { flex:1; }
      .composer .form-control { background:#fff; border:1px solid #bdcbbf; border-radius:11px; box-shadow:none; color:var(--ink); font-size:15px; height:48px; padding:12px 14px; }
      .composer .form-control:focus { border-color:var(--forest); box-shadow:0 0 0 3px rgba(13,95,76,.13); }
      .send-button { background:var(--forest); border:0; border-radius:11px; font-weight:700; height:48px; padding:0 20px; }
      .send-button:hover, .send-button:focus { background:var(--forest-dark); }
      .composer-options { align-items:center; display:flex; flex-wrap:wrap; gap:12px 20px; margin:12px 4px 0; }
      .composer-options .checkbox { margin:0; }
      .composer-options label { color:var(--muted); font-size:13px; font-weight:600; }
      .composer-options .radio { display:inline-block; margin:0 10px 0 0; }
      .composer-options .shiny-options-group { margin:0; }
      details { background:rgba(223,240,231,.36); border-radius:12px; color:var(--muted); margin-top:14px; padding:10px 13px; }
      details summary { color:var(--forest-dark); cursor:pointer; font-size:13px; font-weight:700; }
      .advanced-grid { display:grid; gap:10px; grid-template-columns:1fr 1fr; margin-top:12px; }
      .advanced-grid .form-group { margin-bottom:0; }
      .advanced-grid label { font-size:12px; }
      .advanced-grid .form-control { font-size:13px; height:38px; }
      .parameter-panel { margin-top:18px; }
      .parameter-grid { display:grid; gap:12px; grid-template-columns:repeat(3, 1fr); margin-top:10px; }
      .parameter-card { background:var(--card); border:1px solid var(--line); border-radius:10px; padding:12px; }
      .parameter-card strong { color:var(--forest-dark); display:block; font-size:12px; margin-bottom:4px; }
      .parameter-card span { font-size:12px; line-height:1.4; }
      @media (max-width:700px) { .container-fluid { padding:22px 14px 40px; } .agent-header { display:block; } .privacy-note { display:inline-block; margin-top:14px; } .composer-row { align-items:stretch; flex-direction:column; } .send-button { width:100%; } .advanced-grid, .parameter-grid { grid-template-columns:1fr; } }
    ")),
    tags$script(HTML("\
      $(document).on('keydown', '#question', function(event) {
        if (event.key === 'Enter' && !event.shiftKey) {
          event.preventDefault(); $('#ask').click();
        }
      });
    "))
  ),
  div(
    class = "agent-header",
    div(
      div(class = "agent-mark", "Simulation planning assistant"),
      h1("simitall agent"),
      p(class = "agent-subtitle", "Ask one study-design question. The agent retrieves local project knowledge, identifies required inputs, and returns a reviewable plan with verified R code when available.")
    ),
    div(class = "privacy-note", "Local by default - no API key")
  ),
  div(
    class = "chat-shell",
    div(class = "message assistant-message", div(class = "message-label", "simitall"), uiOutput("answer")),
    div(
      class = "composer",
      div(
        class = "composer-row",
        textInput("question", NULL, placeholder = "Describe the simulation or analysis you want to plan..."),
        actionButton("ask", "Send", class = "btn-primary send-button")
      ),
      div(
        class = "composer-options",
        checkboxInput("include_code", "Include verified runnable code", value = TRUE),
        radioButtons("provider", NULL, choices = c("Free local" = "local", "Ollama" = "ollama", "OpenAI API" = "openai"), selected = "local", inline = TRUE)
      ),
      tags$details(
        tags$summary("Optional input details"),
        div(
          class = "advanced-grid",
          selectInput("data_source", "Data source", choices = c("Ask for required data" = "auto", "Generate synthetic inputs" = "synthetic", "I supplied local files" = "provided", "Use package demo only" = "package_demo"), selected = "auto"),
          textInput("input_files", "Local input file paths", placeholder = "Comma-separated paths; files stay on this computer"),
          numericInput("n_context", "Knowledge documents", value = 6, min = 1, max = 12)
        )
      )
    ),
    tags$details(
      class = "parameter-panel",
      tags$summary("Parameters you can tune in generated recipes"),
      div(
        class = "parameter-grid",
        div(class = "parameter-card", tags$strong("Population"), tags$span("Sample size, founder count, number of generations or inbreeding rounds, and family structure.")),
        div(class = "parameter-card", tags$strong("Genome and variation"), tags$span("Chromosome length, marker density, allele-frequency bounds, SNP/indel rates, and ploidy.")),
        div(class = "parameter-card", tags$strong("Recombination and LD"), tags$span("Genetic-map positions, hotspots/cold spots, interference, and designed high-LD intervals.")),
        div(class = "parameter-card", tags$strong("Trait and omics"), tags$span("Heritability, QTL counts and effects, environments, batches, expression programs, and eQTLs.")),
        div(class = "parameter-card", tags$strong("Technical error"), tags$span("Missingness, genotype error, sequencing coverage, read type, library depth, and random seed.")),
        div(class = "parameter-card", tags$strong("Outputs and QC"), tags$span("VCF/FASTA/GFF3/BED/TSV paths, truth files, figures, and pass criteria."))
      )
    )
  )
)

.simitall_shiny_hide_grounding <- function(answer) {
  answer <- gsub("[[:space:]]*\\[SOURCE: [^]]+\\]", "", answer)
  marker <- "\n\nGrounding documents:\n"
  start <- regexpr(marker, answer, fixed = TRUE)[1L]
  if (start < 0L) return(answer)
  before <- substr(answer, 1L, start - 1L)
  remainder <- substr(answer, start + nchar(marker), nchar(answer))
  recipe_start <- regexpr(
    "\n\n(?:Verified |Synthetic human |Input-dependent )", remainder,
    perl = TRUE
  )[1L]
  if (recipe_start < 0L) return(before)
  paste0(before, substr(remainder, recipe_start, nchar(remainder)))
}

server <- function(input, output, session) {
  response <- eventReactive(input$ask, {
    question <- trimws(input$question)
    req(nzchar(question))
    data_source <- input$data_source
    # In the chat interface, an explicit request to simulate a chromosome is
    # consent to a labelled synthetic chromosome. All other chromosome prompts
    # keep the conservative "ask for data" policy unless the user chooses a
    # different source mode themselves.
    explicit_chromosome_simulation <- grepl(
      "\\bsimulate\\s+(?:a\\s+)?(?:synthetic\\s+)?(?:human\\s+)?(?:chromosome\\s*\\d+|chr\\s*\\d+|chr\\d+)",
      tolower(question), perl = TRUE
    )
    if (identical(data_source, "auto") && explicit_chromosome_simulation) {
      data_source <- "synthetic"
    }
    raw_input_files <- input$input_files
    if (is.null(raw_input_files)) raw_input_files <- ""
    input_files <- trimws(unlist(strsplit(raw_input_files, "[,\\n]")))
    input_files <- input_files[nzchar(input_files)]
    tryCatch(
      simitall_ask(
        question,
        provider = input$provider,
        n_context = input$n_context,
        include_code = isTRUE(input$include_code),
        input_files = input_files,
        data_source = data_source
      ),
      error = function(e) structure(
        list(
          answer = paste("Error:", conditionMessage(e)),
          sources = tryCatch(
            search_simitall_knowledge(question, n_results = input$n_context),
            error = function(search_error) data.frame(category = character(), title = character(), path = character(), stringsAsFactors = FALSE)
          )
        ),
        class = "simitall_agent_answer"
      )
    )
  })

  output$answer <- renderUI({
    if (is.null(response())) {
      return(div(class = "empty-answer", "I am ready when you are. Try: 'Plan a synthetic 500-line maize NAM population for chromosome 10 with a high-LD region, GWAS, RNA-seq, and genomic selection.'"))
    }
    tags$pre(class = "assistant-answer", .simitall_shiny_hide_grounding(response()$answer))
  })
}

shinyApp(ui, server)
