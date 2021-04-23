library(shiny)
library(shinydashboard)
library(shinyjs)

tags$head(
        HTML(
          "
          <script>
          var socket_timeout_interval
          var n = 0
          $(document).on('shiny:connected', function(event) {
          socket_timeout_interval = setInterval(function(){
          Shiny.onInputChange('count', n++)
          }, 15000)
          });
          $(document).on('shiny:disconnected', function(event) {
          clearInterval(socket_timeout_interval)
          });
          </script>
          "
        )
        )

textOutput("keepAlive")

fluidPage(inlineCSS(list(body = "color:DarkBlue")),
          titlePanel("Severe Testing"),
          tabsetPanel(
            tabPanel("Severity curves",
                     tabsetPanel(
                       tabPanel("One sample",
                                sidebarLayout(
                                  sidebarPanel(width=3,
                                               radioButtons("dist", "Distribution", inline=T,
                                                            list("Normal" = "Z", "Student's T" = "t")),
                                               numericInput("obs", "Sample size", min = 2, value = 50),
                                               numericInput("dif", "Observed mean difference", value = 0.3, step=0.1),
                                               numericInput("sd", "Sigma", min = 0, value=1, step=0.1),
                                               numericInput("alpha", "alpha (one-sided)", min=0, max=0.5, value=0.025, step=0.005)),
                                  mainPanel(fluidRow(
                                    splitLayout(style = "border: 1px solid silver:", cellWidths = 625,
                                                column(12, plotOutput("plot1", width = "100%", height = "450px"),
                                                       checkboxInput('adjust', 'Select which x-axis values to display'),
                                                       conditionalPanel(
                                                         condition = "input.adjust" , splitLayout(
                                                           numericInput("from", "From", value=-0.1, step=0.1, width="70%"),
                                                           numericInput("to", "to", value=0.7, step=0.1, width="70%"), numericInput("by", "by", value=0.1, step=0.05, width="70%")))),
                                                sidebarPanel(width=6, textOutput("test"), tableOutput('dat'), tableOutput('dat1'),
                                                             numericInput("severity", "Enter a Severity value", min=0, max=1, step=0.05, value=0.5),
                                                             numericInput("discrepancy", "Discrepancy (gamma)", value = 0),
                                                             numericInput('power', 'Power', value = 0))))))),
                       tabPanel("Paired",
                                sidebarLayout(
                                  sidebarPanel(width=3,
                                               radioButtons("dist2", "Distribution", inline=T,
                                                            list("Normal" = "Z", "Student's T" = "t")),
                                               numericInput("obs2", "Sample size", min = 2, value = 50),
                                               numericInput("dif2", "Observed mean difference", value = 0.3, step=0.1),
                                               selectInput("select", "Choose the input", c('Sigma of the scores differences'='sigmascorediff',
                                                                                           'Correlation, pre and post-test sigmas'='cor')),
                                               conditionalPanel(
                                                 condition = "input.select=='sigmascorediff'" ,
                                                 numericInput("sd2", "Sigma of the score differences", min = 0, value=1, step=0.1)),
                                               conditionalPanel(
                                                 condition = "input.select=='cor'" ,
                                                 numericInput("sda", "Sigma pre-test", min = 0, value=1, step=0.1),
                                                 numericInput("sdb", "Sigma post-test", min = 0, value=1, step=0.1),
                                                 numericInput("correlation", "Correlation", min = 0, max=1, value=0, step=0.1)),
                                               numericInput("alpha2", "alpha (one-sided)", 0.025, min=0, max=0.5, step=0.005)),
                                  mainPanel(fluidRow(
                                    splitLayout(style = "border: 1px solid silver:", cellWidths = 625,
                                                column(12, plotOutput("plot2", width = "100%", height = "450px"),
                                                       checkboxInput('adjust2', 'Select which x-axis values to display'),
                                                       conditionalPanel(
                                                         condition = "input.adjust2" , splitLayout(
                                                           numericInput("from2", "From", value=-0.1, step=0.1, width="70%"),
                                                           numericInput("to2", "to", value=0.7, step=0.1, width="70%"), numericInput("by2", "by", value=0.1, step=0.05, width="70%")))),
                                                sidebarPanel(width=6, textOutput("test2"), tableOutput('dat2'), tableOutput('dat21'),
                                                             numericInput("severity2", "Enter a Severity value", min=0, max=1, step=0.05, value = 0.5),
                                                             numericInput("discrepancy2", "Discrepancy (gamma)", value = 0),
                                                             numericInput('power2', 'Power', value = 0))))))),
                       tabPanel("Two samples",
                                sidebarLayout(
                                  sidebarPanel(width=3,
                                               radioButtons("dist3", "Distribution", inline=T,
                                                            list("Normal" = "Z", "Student's T" = "t")),
                                               numericInput("obs3", "Sample size (control group)", min = 2, value = 50),
                                               numericInput("sd3", "Sigma (control group)", min = 0, value = 1, step=0.1),
                                               numericInput("obs4", "Sample size (Treatment group)", min = 2, value = 50),
                                               numericInput("sd4", "Sigma (Treatment group)", min = 0, value = 1, step=0.1),
                                               numericInput("dif3", "Observed mean difference", value = 0.3, step=0.1),
                                               numericInput("alpha3", "alpha (one-sided)", 0.025, min=0, max=0.5, step=0.005)),
                                  mainPanel(fluidRow(
                                    splitLayout(style = "border: 1px solid silver:", cellWidths = 625,
                                                column(12, plotOutput("plot3", width = "100%", height = "450px"),
                                                       checkboxInput('adjust3', 'Select which x-axis values to display'),
                                                       conditionalPanel(
                                                         condition = "input.adjust3" , splitLayout(
                                                           numericInput("from3", "From", value=-0.2, step=0.1, width="70%"),
                                                           numericInput("to3", "to", value=0.9, step=0.1, width="70%"), numericInput("by3", "by", value=0.1, step=0.05, width="70%")))),
                                                sidebarPanel(width=6, textOutput("test3"), tableOutput('dat3'), tableOutput('dat31'),
                                                             numericInput("severity3", "Enter a Severity value", min=0, max=1, step=0.05, value = 0.5),
                                                             numericInput("discrepancy3", "Discrepancy (gamma)", value = 0),
                                                             numericInput('power3', 'Power', value = 0))))))))),
            tabPanel("Distributions",
                     tabsetPanel(
                       tabPanel("One sample",
                                sidebarLayout(
                                  sidebarPanel(width=3,
                                               radioButtons("dist10", "Distribution", inline=T,
                                                            list("Normal" = "Z", "Student's T" = "t")),
                                               numericInput("obs10", "Sample size", min = 2, value = 50),
                                               numericInput("dif10", "Observed mean difference", value = 0.3, step=0.1),
                                               numericInput("sd10", "Sigma", min = 0, value=1, step=0.1),
                                               numericInput("alpha10", "alpha (one-sided)", min=0, max=0.5, value=0.025, step=0.005),
                                               radioButtons("display", "Display", inline=T,
                                                            list("Severity" = "S", "Power" = "P", "P-value"="P-val"))),
                                  mainPanel(fluidRow(
                                    splitLayout(style = "border: 1px solid silver:", cellWidths = 650,
                                                column(12, plotOutput("plot4", width = "100%", height = "450px"),
                                                       checkboxInput('adjust4', 'Select which x-axis values to display'),
                                                       conditionalPanel(
                                                         condition = "input.adjust4" , splitLayout(
                                                           numericInput("from4", "From", value=-4, step=1, width="70%"),
                                                           numericInput("to4", "to", value=4, step=1, width="70%"), numericInput("by4", "by", value=1, step=0.5, width="70%")))),
                                                sidebarPanel(width=6, textOutput("test10"), tableOutput('dat10'), tableOutput('dat100'),
                                                             numericInput("discrepancy10", "Enter a Discrepancy value", step=0.05, value=0.1),
                                                             numericInput("severity10", "Severity", value = 0),
                                                             numericInput('power10', 'Power', value = 0))))))),
                       tabPanel("Paired",
                                sidebarLayout(
                                  sidebarPanel(width=3,
                                               radioButtons("dist20", "Distribution", inline=T,
                                                            list("Normal" = "Z", "Student's T" = "t")),
                                               numericInput("obs20", "Sample size", min = 2, value = 50),
                                               numericInput("dif20", "Observed mean difference", value = 0.3, step=0.1),
                                               selectInput("select2", "Choose the input", c('Sigma of the scores differences'='sigmascorediff2',
                                                                                            'Correlation, pre and post-test sigmas'='cor2')),
                                               conditionalPanel(
                                                 condition = "input.select2=='sigmascorediff2'" ,
                                                 numericInput("sd20", "Sigma of the score differences", min = 0, value=1, step=0.1)),
                                               conditionalPanel(
                                                 condition = "input.select2=='cor2'" ,
                                                 numericInput("sda2", "Sigma pre-test", min = 0, value=1, step=0.1),
                                                 numericInput("sdb2", "Sigma post-test", min = 0, value=1, step=0.1),
                                                 numericInput("correlation2", "Correlation", min = 0, max=1, value=0, step=0.1)),
                                               numericInput("alpha20", "alpha (one-sided)", 0.025, min=0, max=0.5, step=0.005),
                                               radioButtons("display2", "Display", inline=T,
                                                            list("Severity" = "S", "Power" = "P", "P-value"="P-val"))),
                                  mainPanel(fluidRow(
                                    splitLayout(style = "border: 1px solid silver:", cellWidths = 650,
                                                column(12, plotOutput("plot5", width = "100%", height = "450px"),
                                                       checkboxInput('adjust5', 'Select which x-axis values to display'),
                                                       conditionalPanel(
                                                         condition = "input.adjust5" , splitLayout(
                                                           numericInput("from5", "From", value=-4, step=1, width="70%"),
                                                           numericInput("to5", "to", value=4, step=1, width="70%"), numericInput("by5", "by", value=1, step=0.5, width="70%")))),
                                                sidebarPanel(width=6, textOutput("test20"), tableOutput('dat20'), tableOutput('dat200'),
                                                             numericInput("discrepancy20", "Enter a Discrepancy value", step=0.05, value = 0.1),
                                                             numericInput("severity20", "Severity", value = 0),
                                                             numericInput('power20', 'Power', value = 0))))))),
                       tabPanel("Two samples",
                                sidebarLayout(
                                  sidebarPanel(width=3,
                                               radioButtons("dist30", "Distribution", inline=T,
                                                            list("Normal" = "Z", "Student's T" = "t")),
                                               numericInput("obs30", "Sample size (control group)", min = 2, value = 50),
                                               numericInput("sd30", "Sigma (control group)", min = 0, value = 1, step=0.1),
                                               numericInput("obs40", "Sample size (Treatment group)", min = 2, value = 50),
                                               numericInput("sd40", "Sigma (Treatment group)", min = 0, value = 1, step=0.1),
                                               numericInput("dif30", "Observed mean difference", value = 0.3, step=0.1),
                                               numericInput("alpha30", "alpha (one-sided)", 0.025, min=0, max=0.5, step=0.005),
                                               radioButtons("display3", "Display", inline=T,
                                                            list("Severity" = "S", "Power" = "P", "P-value"="P-val"))),
                                  mainPanel(fluidRow(
                                    splitLayout(style = "border: 1px solid silver:", cellWidths = 650,
                                                column(12, plotOutput("plot6", width = "100%", height = "450px"),
                                                       checkboxInput('adjust6', 'Select which x-axis values to display'),
                                                       conditionalPanel(
                                                         condition = "input.adjust6" , splitLayout(
                                                           numericInput("from6", "From", value=-4, step=1, width="70%"),
                                                           numericInput("to6", "to", value=4, step=1, width="70%"), numericInput("by6", "by", value=1, step=0.5, width="70%")))),
                                                sidebarPanel(width=6, textOutput("test30"), tableOutput('dat30'), tableOutput('dat300'),
                                                             numericInput("discrepancy30", "Enter a Discrepancy value", step=0.05, value = 0.1),
                                                             numericInput("severity30", "Severity", value = 0),
                                                             numericInput('power30', 'Power', value = 0)))))))))))
