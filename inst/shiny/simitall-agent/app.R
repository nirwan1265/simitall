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
      .question-message { background:var(--mint); border-color:#b8d5c5; }
      .question-text { color:var(--ink); font-size:16px; line-height:1.55; padding:12px 20px 20px; }
      .empty-answer { color:var(--muted); font-size:16px; line-height:1.55; padding:14px 20px 20px; }
      .answer-browser { background:var(--card); border:1px solid var(--line); border-radius:18px; box-shadow:0 10px 24px rgba(27,52,38,.06); margin-bottom:18px; padding:16px; }
      .section-nav { display:grid; gap:9px; grid-template-columns:repeat(3, 1fr); }
      .section-card { background:#f5f7f1; border:1px solid var(--line); border-radius:11px; color:var(--forest-dark); cursor:pointer; font-size:13px; font-weight:800; min-height:48px; padding:10px 12px; text-align:left; transition:background .15s ease, border-color .15s ease, transform .15s ease; }
      .section-card:hover { background:var(--mint); border-color:#9abfa9; transform:translateY(-1px); }
      .section-card.is-active { background:var(--forest); border-color:var(--forest); color:#fff; }
      .section-detail { background:#fbfcf8; border:1px solid var(--line); border-radius:12px; margin-top:14px; min-height:185px; overflow:hidden; }
      .section-detail-label { background:rgba(223,240,231,.58); color:var(--forest-dark); font-size:11px; font-weight:800; letter-spacing:.12em; padding:12px 15px; text-transform:uppercase; }
      .section-text { color:var(--ink); font-size:14px; line-height:1.65; padding:16px; white-space:normal; word-break:break-word; }
      .section-code { background:#17231e; border:0; border-radius:0; color:#eef6ef; font-family:ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace; font-size:12px; line-height:1.55; margin:0; max-height:600px; overflow:auto; padding:16px; white-space:pre; }
      .run-button { background:#d97435; border:0; border-radius:12px; font-size:15px; font-weight:850; letter-spacing:.08em; margin:0 0 18px; min-height:52px; width:100%; }
      .run-button:hover, .run-button:focus { background:#b95524; }
      .run-note { color:var(--muted); font-size:12px; margin:-9px 0 16px; text-align:center; }
      .parameter-inputs { display:grid; gap:10px; grid-template-columns:repeat(3, 1fr); padding:0 16px 16px; }
      .parameter-inputs .form-group { margin:0; }
      .parameter-inputs label { color:var(--forest-dark); font-size:12px; font-weight:700; }
      .parameter-inputs .form-control { border-color:#bdcbbf; font-size:13px; height:38px; }
      .upload-area { border-top:1px solid var(--line); margin-top:8px; padding:14px 16px 16px; }
      .upload-area label { color:var(--forest-dark); font-size:13px; font-weight:800; }
      .upload-note { color:var(--muted); font-size:12px; line-height:1.45; margin-top:6px; }
      .run-output { border-top:1px solid var(--line); margin-top:8px; padding:14px 16px 16px; }
      .run-output pre { background:#17231e; border-radius:9px; color:#eef6ef; font-size:12px; line-height:1.45; max-height:240px; overflow:auto; padding:12px; white-space:pre-wrap; }
      .run-output .btn { background:var(--forest); border:0; border-radius:9px; font-weight:700; margin:0 0 12px; }
      .run-status { background:#17231e; border:1px solid #365f4e; border-radius:14px; bottom:18px; box-shadow:0 14px 34px rgba(0,0,0,.22); color:#eff7f0; max-width:min(390px, calc(100vw - 36px)); padding:13px 15px; position:fixed; right:18px; z-index:1200; }
      .run-status-label { color:#a8dfbb; font-size:10px; font-weight:850; letter-spacing:.13em; text-transform:uppercase; }
      .run-status-stage { font-size:14px; font-weight:800; margin:4px 0; }
      .run-status-detail { color:#c9d6cd; font-family:ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace; font-size:11px; line-height:1.4; overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }
      .run-status.is-running { border-color:#e4a05d; }
      .run-status.is-error { border-color:#d97767; }
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
      .knowledge-scale { color:var(--muted); display:flex; font-size:11px; justify-content:space-between; margin-top:-6px; }
      .advanced-grid .form-control { font-size:13px; height:38px; }
      .parameter-panel { margin-top:18px; }
      .parameter-grid { display:grid; gap:12px; grid-template-columns:repeat(3, 1fr); margin-top:10px; }
      .parameter-card { background:var(--card); border:1px solid var(--line); border-radius:10px; padding:12px; }
      .parameter-card strong { color:var(--forest-dark); display:block; font-size:12px; margin-bottom:4px; }
      .parameter-card span { font-size:12px; line-height:1.4; }
      @media (max-width:700px) { .container-fluid { padding:22px 14px 40px; } .agent-header { display:block; } .privacy-note { display:inline-block; margin-top:14px; } .composer-row { align-items:stretch; flex-direction:column; } .send-button { width:100%; } .advanced-grid, .parameter-grid, .section-nav, .parameter-inputs { grid-template-columns:1fr; } }
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
    uiOutput("question_top"),
    div(
      class = "answer-browser",
      uiOutput("section_nav"),
      uiOutput("section_detail")
    ),
    uiOutput("run_button"),
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
          div(
            sliderInput("context_level", "Knowledge depth", min = 1, max = 5,
                        value = 3, step = 1, ticks = TRUE),
            div(class = "knowledge-scale",
                span("Low: ~5"), span("Medium: ~20"), span("High: ~50"))
          )
        )
      )
    )
  ),
  uiOutput("run_status")
)

# A focused evidence set is more reliable than every document at once. These
# levels intentionally expand from a quick answer to a broad literature scan.
.simitall_shiny_context_count <- function(level) {
  levels <- c(5L, 10L, 20L, 35L, 50L)
  level <- suppressWarnings(as.integer(level))
  if (is.na(level) || level < 1L || level > length(levels)) level <- 3L
  levels[[level]]
}

.simitall_shiny_hide_grounding <- function(answer) {
  answer <- gsub("[[:space:]]*\\[SOURCE: [^]]+\\]", "", answer)
  marker <- "\n\nGrounding documents:\n"
  start <- regexpr(marker, answer, fixed = TRUE)[1L]
  if (start < 0L) return(answer)
  before <- substr(answer, 1L, start - 1L)
  remainder <- substr(answer, start + nchar(marker), nchar(answer))
  recipe_start <- regexpr(
    "\n\n(?:Verified |Synthetic human |Input-dependent |Bundled synthetic )", remainder,
    perl = TRUE
  )[1L]
  if (recipe_start < 0L) return(before)
  paste0(before, substr(remainder, recipe_start, nchar(remainder)))
}

.simitall_shiny_extract_section <- function(answer, heading, following) {
  marker <- paste0("\n\n", heading, ":\n")
  start <- regexpr(marker, answer, fixed = TRUE)[1L]
  if (start < 0L) return("")
  content_start <- start + nchar(marker)
  remaining <- substr(answer, content_start, nchar(answer))
  if (!length(following)) return(trimws(remaining))
  endings <- vapply(
    following,
    function(next_heading) {
      match <- regexpr(paste0("\n\n", next_heading, ":\n"), remaining, fixed = TRUE)[1L]
      if (match < 0L) nchar(remaining) + 1L else match
    },
    integer(1)
  )
  trimws(substr(remaining, 1L, min(endings) - 1L))
}

.simitall_shiny_code_blocks <- function(answer) {
  lines <- strsplit(answer, "\n", fixed = TRUE)[[1L]]
  fences <- which(grepl("^```", lines))
  if (length(fences) < 2L) return("")
  starts <- fences[seq(1L, length(fences) - 1L, by = 2L)]
  ends <- fences[seq(2L, length(fences), by = 2L)]
  paste(vapply(
    seq_along(starts),
    function(i) paste(lines[(starts[i] + 1L):(ends[i] - 1L)], collapse = "\n"),
    character(1)
  ), collapse = "\n\n# ---- Next verified recipe block ----\n\n")
}

.simitall_shiny_parameter_text <- function() {
  paste(
    "Population: sample size, founder count, generations/inbreeding rounds, and family structure.",
    "\n\nGenome and variation: chromosome length, marker density, allele-frequency bounds, SNP/indel rates, and ploidy.",
    "\n\nRecombination and LD: map positions, hotspots/coldspots, interference, and designed high-LD intervals.",
    "\n\nTrait and omics: heritability, QTL effects, environments, batches, expression programs, and eQTLs.",
    "\n\nTechnical error and reproducibility: missingness, genotype error, coverage, read type, library depth, and seed.",
    "\n\nOutputs and QC: FASTA/VCF/GFF3/BED/TSV paths, truth files, figures, and explicit pass criteria."
  )
}

.simitall_shiny_default <- function(value, fallback) {
  if (is.null(value) || !length(value) || is.na(value)) fallback else value
}

.simitall_shiny_recipe_defaults <- function(code) {
  first_number <- function(pattern, fallback) {
    hit <- regmatches(code, regexpr(pattern, code, perl = TRUE))
    if (!length(hit) || !nzchar(hit)) return(fallback)
    value <- suppressWarnings(as.numeric(sub(".*?([0-9]+(?:\\.[0-9]+)?).*", "\\1", hit, perl = TRUE)))
    if (is.finite(value)) value else fallback
  }
  list(
    seed = first_number("(?:seed\\s*<-|seed\\s*=)\\s*[0-9]+", 2026),
    n_lines = first_number("(?:n_lines\\s*<-|n_lines\\s*=|n_offspring\\s*=)\\s*[0-9]+", 500),
    n_generations = first_number("(?:n_rounds\\s*<-|n_generations\\s*<-)\\s*[0-9]+", 10),
    n_markers = first_number("(?:n_markers\\s*<-|n_markers\\s*=)\\s*[0-9]+", 500),
    heritability = first_number("heritability\\s*=\\s*[0-9]+(?:\\.[0-9]+)?", 0.6)
  )
}

.simitall_shiny_apply_overrides <- function(code, values) {
  # These substitutions cover the common literal recipe knobs. A recipe that
  # does not expose a knob is left unchanged rather than silently rewritten.
  replace_assignment <- function(text, name, value) {
    gsub(
      paste0("(?m)(^\\s*", name, "\\s*<-\\s*)[^#\\n]+"),
      paste0("\\1", value), text, perl = TRUE
    )
  }
  replace_argument <- function(text, name, value) {
    gsub(
      paste0("(\\b", name, "\\s*=\\s*)[0-9.]+L?"),
      paste0("\\1", value), text, perl = TRUE
    )
  }
  code <- replace_assignment(code, "seed", values$seed)
  code <- replace_argument(code, "seed", values$seed)
  code <- replace_assignment(code, "n_lines", paste0(values$n_lines, "L"))
  code <- replace_argument(code, "n_lines", paste0(values$n_lines, "L"))
  code <- replace_argument(code, "n_offspring", values$n_lines)
  code <- replace_assignment(code, "n_rounds", paste0(values$n_generations, "L"))
  code <- replace_assignment(code, "n_generations", paste0(values$n_generations, "L"))
  code <- replace_assignment(code, "n_markers", paste0(values$n_markers, "L"))
  code <- replace_argument(code, "n_markers", values$n_markers)
  code <- replace_argument(code, "heritability", values$heritability)
  code
}

.simitall_shiny_project_dir <- function() {
  path <- normalizePath(getwd(), mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) return(normalizePath(getwd(), mustWork = TRUE))
    path <- parent
  }
}

.simitall_shiny_result_files <- function(project_dir, started_at) {
  roots <- c(file.path(project_dir, "results"), file.path(project_dir, "analysis", "results"))
  files <- unlist(lapply(roots[file.exists(roots)], function(root) {
    list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  }), use.names = FALSE)
  files <- files[file.info(files)$isdir %in% FALSE]
  files[file.info(files)$mtime >= started_at - 2]
}

.simitall_shiny_sections <- function(answer, question) {
  answer <- .simitall_shiny_hide_grounding(answer)
  headings <- c(
    "Interpreted design", "Proposed stages", "Required inputs",
    "Truth and QC outputs", "Limitations", "Expected outputs"
  )
  stages <- .simitall_shiny_extract_section(answer, "Proposed stages", headings[-2L])
  design <- .simitall_shiny_extract_section(answer, "Interpreted design", headings[-1L])
  inputs <- .simitall_shiny_extract_section(answer, "Required inputs", headings[-3L])
  outputs <- paste(
    .simitall_shiny_extract_section(answer, "Truth and QC outputs", headings[-4L]),
    .simitall_shiny_extract_section(answer, "Expected outputs", character()),
    sep = "\n\n"
  )
  limitations <- .simitall_shiny_extract_section(answer, "Limitations", headings[-5L])
  # A verified recipe is deliberately appended after the planning prose. It
  # belongs only in the R Code panel, never inside the Limitations panel.
  limitations <- sub(
    "(?s)\\n\\n(?:Input-dependent|Verified|Synthetic|Package-verified)[^\\n]*:\\n.*$",
    "", limitations, perl = TRUE
  )
  plan <- if (nzchar(stages)) {
    paste0("A simulation/analysis plan was generated with ", length(grep("^[0-9]+\\.", strsplit(stages, "\n")[[1L]])), " proposed stage(s). Select Design, Stages, Inputs, Outputs, or R Code for the details.")
  } else {
    "The agent returned a planning response. Select the available sections for details."
  }
  list(
    Question = question,
    Plan = plan,
    Design = if (nzchar(design)) design else "No separate interpreted design was returned for this question.",
    Stages = if (nzchar(stages)) stages else "No staged workflow was returned for this question.",
    Inputs = if (nzchar(inputs)) inputs else "No additional inputs were identified.",
    Limitations = if (nzchar(limitations)) limitations else "No additional limitations were returned.",
    `R Code` = .simitall_shiny_code_blocks(answer),
    `Parameter Tunings` = .simitall_shiny_parameter_text(),
    Outputs = if (nzchar(trimws(outputs))) trimws(outputs) else "No output contract was returned."
  )
}

server <- function(input, output, session) {
  selected_section <- reactiveVal("Plan")
  run_state <- reactiveValues(
    log = "No recipe has been run in this session.", files = character(),
    completed = FALSE, running = FALSE, stage = "Ready",
    detail = "No local recipe is running.", process = NULL, started_at = NULL
  )

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
    uploaded_paths <- input$uploaded_files$datapath
    if (is.null(uploaded_paths)) uploaded_paths <- character()
    input_files <- unique(c(input_files, uploaded_paths))
    n_context <- .simitall_shiny_context_count(input$context_level)
    tryCatch(
      simitall_ask(
        question,
        provider = input$provider,
        n_context = n_context,
        include_code = isTRUE(input$include_code),
        input_files = input_files,
        data_source = data_source
      ),
      error = function(e) structure(
        list(
          answer = paste("Error:", conditionMessage(e)),
          sources = tryCatch(
            search_simitall_knowledge(question, n_results = n_context),
            error = function(search_error) data.frame(category = character(), title = character(), path = character(), stringsAsFactors = FALSE)
          )
        ),
        class = "simitall_agent_answer"
      )
    )
  })

  observeEvent(response(), {
    selected_section("Plan")
    if (!isTRUE(run_state$running)) {
      run_state$stage <- "Ready"
      run_state$detail <- "A new recipe is ready to run locally."
      run_state$completed <- FALSE
      run_state$files <- character()
      run_state$log <- "No recipe has been run for this answer."
    }
  }, ignoreInit = TRUE)

  observeEvent(input$answer_section, {
    selected_section(input$answer_section)
  })

  sections <- reactive({
    req(response())
    .simitall_shiny_sections(
      response()$answer,
      trimws(input$question)
    )
  })

  output$question_top <- renderUI({
    if (is.null(response())) {
      return(div(
        class = "message assistant-message",
        div(class = "message-label", "simitall"),
        div(class = "empty-answer", "I am ready when you are. Try: 'Plan a synthetic 500-line maize NAM population for chromosome 10 with a high-LD region, GWAS, RNA-seq, and genomic selection.'")
      ))
    }
    div(
      class = "message question-message",
      div(class = "message-label", "Question"),
      div(class = "question-text", sections()[["Question"]])
    )
  })

  output$section_nav <- renderUI({
    if (is.null(response())) return(NULL)
    choices <- names(sections())
    div(
      class = "section-nav",
      lapply(choices, function(choice) {
        is_active <- identical(choice, selected_section())
        tags$button(
          type = "button",
          class = paste("section-card", if (is_active) "is-active"),
          onclick = sprintf(
            "Shiny.setInputValue('answer_section', '%s', {priority: 'event'});",
            choice
          ),
          choice
        )
      })
    )
  })

  output$section_detail <- renderUI({
    if (is.null(response())) {
      return(div(
        class = "section-detail",
        div(class = "section-detail-label", "Answer sections"),
        div(class = "section-text", "After you send a question, use these cards to inspect the plan, inputs, outputs, runnable R code, and parameters separately.")
      ))
    }
    choice <- selected_section()
    content <- sections()[[choice]]
    if (identical(choice, "R Code")) {
      if (!nzchar(content)) {
        content <- "No verified R code was returned. Enable 'Include verified runnable code', provide any required local inputs, or ask for a synthetic design."
      }
      body <- tags$pre(class = "section-code", content)
    } else if (identical(choice, "Inputs")) {
      upload_summary <- if (is.null(input$uploaded_files)) {
        div(class = "upload-note", "No files uploaded yet. Upload data here, then press Send again so the next plan can use those local files.")
      } else {
        tags$ul(class = "upload-note", lapply(seq_len(nrow(input$uploaded_files)), function(i) {
          paste0(input$uploaded_files$name[i], " (", format(input$uploaded_files$size[i], big.mark = ","), " bytes)")
        }))
      }
      body <- tagList(
        div(class = "section-text", HTML(gsub("\n", "<br>", htmltools::htmlEscape(content), fixed = TRUE))),
        div(
          class = "upload-area",
          fileInput("uploaded_files", "Upload local input files", multiple = TRUE,
                    accept = c(".fa", ".fasta", ".fna", ".vcf", ".vcf.gz", ".bcf", ".tsv", ".csv", ".gff", ".gff3", ".gtf", ".bed", ".map")),
          upload_summary
        )
      )
    } else if (identical(choice, "Parameter Tunings")) {
      defaults <- .simitall_shiny_recipe_defaults(sections()[["R Code"]])
      body <- tagList(
        div(class = "section-text", HTML(gsub("\n", "<br>", htmltools::htmlEscape(content), fixed = TRUE))),
        div(
          class = "parameter-inputs",
          numericInput("recipe_seed", "Random seed", value = defaults$seed, min = 1, step = 1),
          numericInput("recipe_n_lines", "Lines / offspring", value = defaults$n_lines, min = 2, step = 1),
          numericInput("recipe_n_generations", "Generations / rounds", value = defaults$n_generations, min = 1, step = 1),
          numericInput("recipe_n_markers", "Markers", value = defaults$n_markers, min = 10, step = 10),
          numericInput("recipe_heritability", "Heritability", value = defaults$heritability, min = 0, max = 1, step = 0.05)
        )
      )
    } else if (identical(choice, "Outputs")) {
      body <- tagList(
        div(class = "section-text", HTML(gsub("\n", "<br>", htmltools::htmlEscape(content), fixed = TRUE))),
        uiOutput("run_output")
      )
    } else {
      body <- div(class = "section-text", HTML(gsub("\n", "<br>", htmltools::htmlEscape(content), fixed = TRUE)))
    }
    div(
      class = "section-detail",
      div(class = "section-detail-label", choice),
      body
    )
  })

  output$run_button <- renderUI({
    if (is.null(response())) return(NULL)
    has_code <- nzchar(sections()[["R Code"]])
    tagList(
      actionButton("run_recipe", "RUN", class = "btn-primary run-button",
                   disabled = if (!has_code || isTRUE(run_state$running)) "disabled"),
      div(class = "run-note", "Runs the displayed R recipe locally with the current parameter values. Outputs remain on this computer.")
    )
  })

  observeEvent(input$run_recipe, {
    req(response())
    code <- sections()[["R Code"]]
    validate(need(nzchar(code), "No runnable R code is available for this answer."))
    if (grepl("PATH_TO_", code, fixed = TRUE)) {
      run_state$completed <- FALSE
      run_state$running <- FALSE
      run_state$files <- character()
      run_state$stage <- "Needs updated recipe"
      run_state$detail <- "Choose a data source, then press Send to regenerate the recipe."
      run_state$log <- paste(
        "This recipe still contains PATH_TO_* placeholders, so it was not run.",
        "Upload or provide the requested real files, or choose 'Use package demo only' under Optional input details and send the question again.",
        sep = "\n"
      )
      selected_section("Outputs")
      return(invisible(NULL))
    }
    values <- list(
      seed = .simitall_shiny_default(input$recipe_seed, 2026),
      n_lines = .simitall_shiny_default(input$recipe_n_lines, 500),
      n_generations = .simitall_shiny_default(input$recipe_n_generations, 10),
      n_markers = .simitall_shiny_default(input$recipe_n_markers, 500),
      heritability = .simitall_shiny_default(input$recipe_heritability, 0.6)
    )
    code <- .simitall_shiny_apply_overrides(code, values)
    project_dir <- .simitall_shiny_project_dir()
    run_dir <- file.path(project_dir, "results", "simitall_runs", format(Sys.time(), "%Y%m%d_%H%M%S"))
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    script <- file.path(run_dir, "recipe.R")
    # dQuote() can emit typographic quotes on macOS; encodeString() always
    # writes a parser-safe R string literal for the local project directory.
    project_literal <- encodeString(project_dir, quote = "\"")
    bootstrap <- c(
      sprintf("setwd(%s)", project_literal),
      "# RUN starts a fresh R process, so load this development package first.",
      "if (!requireNamespace('simitall', quietly = TRUE)) {",
      "  if (!requireNamespace('devtools', quietly = TRUE)) {",
      "    stop('simitall is not installed and devtools is unavailable to load this project.')",
      "  }",
      sprintf("  devtools::load_all(%s, quiet = TRUE)", project_literal),
      "}"
    )
    writeLines(c(bootstrap, code), script)
    if (!requireNamespace("processx", quietly = TRUE)) {
      run_state$completed <- FALSE
      run_state$stage <- "Cannot start"
      run_state$detail <- "Install the R package processx to run recipes in the background."
      run_state$log <- run_state$detail
      selected_section("Outputs")
      return(invisible(NULL))
    }
    run_state$started_at <- Sys.time()
    run_state$running <- TRUE
    run_state$completed <- FALSE
    run_state$files <- character()
    run_state$stage <- "Running"
    run_state$detail <- paste("Executing", basename(script))
    run_state$log <- paste("Recipe:", script, "\nProcess started in the background.")
    run_state$process <- processx::process$new(
      command = file.path(R.home("bin"), "Rscript"), args = script,
      stdout = "|", stderr = "|", cleanup = FALSE
    )
    selected_section("Outputs")
  })

  observe({
    invalidateLater(600, session)
    if (!isTRUE(run_state$running) || is.null(run_state$process)) return()
    process <- run_state$process
    output_lines <- tryCatch(process$read_output_lines(), error = function(e) character())
    error_lines <- tryCatch(process$read_error_lines(), error = function(e) character())
    new_lines <- c(output_lines, error_lines)
    if (length(new_lines)) {
      run_state$log <- paste(c(run_state$log, new_lines), collapse = "\n")
      run_state$detail <- tail(new_lines, 1L)
    }
    if (isTRUE(process$is_alive())) return()
    status <- process$get_exit_status()
    files <- .simitall_shiny_result_files(.simitall_shiny_project_dir(), run_state$started_at)
    run_state$files <- files
    run_state$running <- FALSE
    run_state$completed <- identical(status, 0L)
    run_state$stage <- if (run_state$completed) "Completed" else "Error"
    run_state$detail <- if (run_state$completed) {
      paste("Finished:", length(files), "output file(s) detected")
    } else {
      paste("Rscript exited with status", status)
    }
    run_state$log <- paste(run_state$log, paste("Process exited with status", status), sep = "\n")
    run_state$process <- NULL
  })

  output$run_output <- renderUI({
    if (!isTRUE(run_state$completed) && identical(run_state$log, "No recipe has been run in this session.")) {
      return(div(class = "run-output", div(class = "upload-note", "Run a recipe to see generated figures, files, and the local execution log here.")))
    }
    div(
      class = "run-output",
      if (length(run_state$files)) downloadButton("download_run", "Download run outputs") else NULL,
      div(class = "upload-note", if (run_state$running) "Recipe is running in the background; the log updates live." else if (run_state$completed) paste("Run completed. Detected", length(run_state$files), "new or updated output file(s), including any validation figures. Use Download run outputs to retrieve them.") else "Run returned an error. Inspect the local log below."),
      tags$pre(run_state$log)
    )
  })

  output$download_run <- downloadHandler(
    filename = function() paste0("simitall_run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip"),
    content = function(file) {
      req(length(run_state$files))
      project_dir <- .simitall_shiny_project_dir()
      files <- normalizePath(run_state$files, mustWork = TRUE)
      # Zip relative paths so extracting never recreates /Users/... folders.
      relative_files <- substring(files, nchar(project_dir) + 2L)
      old_wd <- setwd(project_dir)
      on.exit(setwd(old_wd), add = TRUE)
      utils::zip(file, files = relative_files)
    }
  )

  output$run_status <- renderUI({
    status_class <- if (isTRUE(run_state$running)) "is-running" else if (!isTRUE(run_state$completed) && !identical(run_state$stage, "Ready")) "is-error" else ""
    div(
      class = paste("run-status", status_class),
      div(class = "run-status-label", "Local runner"),
      div(class = "run-status-stage", run_state$stage),
      div(class = "run-status-detail", run_state$detail)
    )
  })
}

shinyApp(ui, server)
