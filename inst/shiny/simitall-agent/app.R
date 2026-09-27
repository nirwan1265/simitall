library(shiny)
library(simitall)

ui <- fluidPage(
  titlePanel("simitall agent"),
  tags$p(
    "Local, retrieval-grounded help for simitall. This first release answers",
    "questions and plans workflows; it never launches simulations."
  ),
  sidebarLayout(
    sidebarPanel(
      textAreaInput(
        "question", "Question", rows = 7,
        placeholder = "Example: Plan a multi-chromosome maize NAM population for GWAS."
      ),
      selectInput(
        "provider", "Answer backend",
        choices = c(
          "Free local evidence briefing" = "local",
          "Free Ollama model on this Mac" = "ollama",
          "OpenAI API (paid)" = "openai"
        ),
        selected = "local"
      ),
      numericInput("n_context", "Knowledge documents", value = 6, min = 1, max = 12),
      actionButton("ask", "Ask simitall", class = "btn-primary"),
      tags$hr(),
      tags$small(
        "The default local mode is free and needs no key. Ollama runs a free model",
        "on this Mac after installation. OpenAI is optional and uses paid API billing."
      )
    ),
    mainPanel(
      h3("Answer"),
      verbatimTextOutput("answer"),
      h3("Grounding documents"),
      tableOutput("sources")
    )
  )
)

server <- function(input, output, session) {
  response <- eventReactive(input$ask, {
    req(nzchar(trimws(input$question)))
    tryCatch(
      simitall_ask(
        input$question,
        provider = input$provider,
        n_context = input$n_context
      ),
      error = function(e) structure(
        list(
          answer = paste("Error:", conditionMessage(e)),
          sources = tryCatch(
            search_simitall_knowledge(input$question, n_results = input$n_context),
            error = function(search_error) data.frame(
              category = character(), title = character(), path = character(),
              stringsAsFactors = FALSE
            )
          )
        ),
        class = "simitall_agent_answer"
      )
    )
  })

  output$answer <- renderText({
    req(response())
    response()$answer
  })
  output$sources <- renderTable({
    req(response())
    sources <- response()$sources
    required <- c("category", "title", "path")
    if (!all(required %in% names(sources))) {
      sources <- data.frame(
        category = character(), title = character(), path = character(),
        stringsAsFactors = FALSE
      )
    }
    sources[, required, drop = FALSE]
  }, striped = TRUE, bordered = TRUE, spacing = "s")
}

shinyApp(ui, server)
