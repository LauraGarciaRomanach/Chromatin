library(shiny)
library(rtracklayer)
library(GenomicRanges)
library(Gviz)
library(dplyr)
library(tidyr)
options(ucscChromosomeNames = FALSE) 

gene_index <- readRDS("annotations/gene_index_sel.rds")

sample_names <- c("LD", "SD4W", "SD10W", "LT6W", "WT2W")
colors <- c("#2E8B57",
            "#7BA05B",
            "#BDB76B",
            "#4682B4",
            "#5F9EA0")

#Add a slider to choose how many kb up and down you want to see
ui <- fluidPage(
  titlePanel("PopEpigen"),
  sidebarLayout(
    sidebarPanel(
      selectizeInput(
        inputId = "gene",
        label = "Enter Gene ID",
        choices = NULL,
        selected = character(0),  # Use character(0) instead of NULL
        options = list(
          placeholder = 'Type gene ID',
          maxItems = 1  # Force empty initialization
        )
      ),
      width = 2
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Expression profile and chromatin accessibility",
                 fluidRow(
                   column(6, plotOutput("gene_plot", height = "600px")),
                   column(6, plotOutput("atac_plot", height = "600px")),
                   column(12, style = "margin-top: 30px;",  htmlOutput("rna_atac_foot"))
                   )),
        
        #          plotOutput("gene_plot", height = "800px")),
        # tabPanel("ATAC-Seq",
        #          plotOutput("atac_plot", height = "800px")),
        tabPanel("Chromatin modification",
                 fluidRow(
                   # column(6, plotOutput("gene_plot", height = "500px")),
                   # column(6, plotOutput("atac_plot", height = "500px")),
                   column(6, plotOutput("k4", height = "600px")),
                   column(6, plotOutput("k27me", height = "600px")),
                   column(6, plotOutput("k27ac", height = "600px")),
                   column(12, style = "margin-top: 30px;",  htmlOutput("chr_foot"))
                 ))

        )
      )
    )
  )



server <- function(input, output, session) {

  updateSelectizeInput(
    session,
    'gene',
    choices = gene_index$geneID,
    server = TRUE,
    selected = ""
  )

  gene <- reactive({
    req(input$gene)

    # Lookup chromosome
  gene_index[gene_index$geneID == input$gene,]

  })
    # validate(need(nrow(gene_info) > 0, "Gene not found in index"))
   annot <- reactive({
    # Load only the required chromosome's annotation
    readRDS(file.path("annotations",
                               paste0(gene()$chromosome, ".rds")))
  })

  gene_data <- reactive({
    gene2 <- annot()[
      annot()$type == "gene" &
        annot()$ID == input$gene
    ]
    #
    # # Validate that gene was found
    # if (length(gene) == 0) {
    #   showNotification("Gene not found in annotation", type = "error")
    #   return(NULL)
    # }

    list(
      gene = gene2,
      transcripts = annot()[annot()$type == "mRNA" & grepl(input$gene, annot()$ID)],
      exons = annot()[annot()$type %in% c("exon") & grepl(input$gene, annot()$ID)],
      gene_region = GRanges(
        seqnames = seqnames(gene2),
        ranges = IRanges(
          start = start(gene2) - 5000,
          end = end(gene2) + 5000
        )
      )
    )
  })

  # Function to load just the needed chromosome data
  load_track_data <- function(group_name, sample_names, gene_region) {
    chr <- as.character(seqnames(gene_region)[1])
    lapply(sample_names, function(sample) {
      file <- file.path("RData", group_name, paste0(sample, "_", chr, "_sel.rds"))
      if(file.exists(file)) {
        gr <- readRDS(file)
        subsetByOverlaps(gr, gene_region)
      } else {
        NULL
      }
    })
  }

  gene_track <- reactive({
    data <- gene_data()
    req(data)

    # annot <- readRDS("annotations/chr8.rds")

    gene3 <- annot()[
      annot()$type == "gene" &
        annot()$ID == input$gene
    ]

    transcripts <- annot()[annot()$type == "mRNA" & grepl(input$gene, annot()$ID)]
    # exons <- annot[annot$type == "exon" & grepl("Potra2n8c17315", annot$ID)]
    # mcols(exons)$transcript <- mcols(exons)$Parent
    # utrs <- annot[grepl("UTR", annot$type) & grepl("Potra2n8c17315", annot$ID)]

    features <- annot()[ annot()$type %in% c("exon") & grepl(input$gene, annot()$ID)]

    mcols(features)$transcript <- mcols(features)$Parent

  GeneRegionTrack(
      range = features,
      chromosome = seqnames(data$gene_region),
      name = mcols(gene3)$ID,
      fill = "grey",
      col = "black",
      shape = "smallArrow",
      showId = TRUE,
      transcriptAnnotation = "transcript",
      collapseTranscripts = "longest",  # Changed from "meta"
      just.group = "below",             # Position labels
      geneSymbol = TRUE  ,
      background.title = "transparent",
      cex.group = 1,
      col.group = "black"
    )



    # exons <- data$exons
    # mcols(exons)$transcript <- mcols(exons)$Parent
    #
    # gene_track <-  GeneRegionTrack(
    #   range = exons,
    #   chromosome = seqnames(data$gene_region),
    #   name = mcols(data$gene)$ID,
    #   fill = "grey",
    #   col = "black",
    #   shape = "arrow",
    #   arrowHeadWidth = 12,
    #   arrowHeadMaxWidth = 15,
    #   arrowTailWidth = 2,
    #   showId = TRUE,
    #   transcriptAnnotation = "transcript",
    #   collapseTranscripts = "longest",
    #   just.group = "below",
    #   geneSymbol = TRUE,
    #   background.title = "transparent",
    #   cex.group = 1.2,
    #   col.group = "black"
    # )
  })

  genome_axis <- reactive({
    req(gene_data())

    GenomeAxisTrack(
      col = "darkgrey",
      fontcolor = "black",
      exponent = 3,  # Display as kb (use 6 for Mb)
      labelPos = "below"
    )
  })

  # Plot output
  plot_track <- function(group_name, title){
    data <- gene_data()
    req(data, gene_track(), genome_axis())

    track_data <- load_track_data(
      group_name,
      c("LD", "SD4W", "SD10W", "LT6W", "WT2W"),
      data$gene_region
    )

    # global_max <- max(
    #   sapply(input_file, function(file) {
    #     data <- readRDS(file, which = data$gene_region)
    #     max(data$score, na.rm = TRUE)
    #   }
    #   ))

    global_max <- max(sapply(track_data, function(x) max(x$score, na.rm = TRUE)))

    track <- lapply(1:5, function(i) {
      DataTrack(
        range = track_data[[i]],
        # range = import(input_file[i], which = data$gene_region),
        type = "histogram",
        name = sample_names[i],
        col.histogram = colors[i],
        fill.histogram = colors[i],
        ylim = c(0, global_max),
        window = 150,
        showAxis = TRUE,
        col.axis = "black",
        background.title = "transparent",
        col.title = "black",
        cex.axis = 0.8,
        cex.title = 1.2
      )
    })

    plotTracks(
      c(track, list(gene_track(), genome_axis())),
      from = start(data$gene_region),
      to = end(data$gene_region),
      chromosome = seqnames(data$gene_region),
      sizes = c(rep(0.8, 5), 0.4, 0.5),
      main = title)
  }

  output$gene_plot <- renderPlot({
    req(input$gene)
    plot_track("rna", "RNA-Seq")
  })

  output$atac_plot <- renderPlot({
    req(input$gene)
    plot_track("atac","ATAC-Seq")
  })
  
  output$rna_atac_foot <- renderUI({
    HTML(
    paste("<p style='text-align:justify'>","Tracks show RNA-seq (transcriptional activity) and ATAC-seq (chromatin accessibility) signal intensities across the gene locus. Signal heights correspond to bins per million mapped reads (BPM). Exons are depicted as grey boxes. LD: long day (7 weeks). SD: short day. LT: low temperature. WT: warm temperature.  The number shown is the number of weeks in the treatment.", "</p>"
    ))
  })

  output$k4 <- renderPlot({
    req(input$gene)
    plot_track("k4", "H3K4me3")
  })

  output$k27me <- renderPlot({
    req(input$gene)
    plot_track("k27me", "H3K27me3")
  })

  output$k27ac <- renderPlot({
    req(input$gene)
    plot_track("k27ac", "H3K27ac")
  })
  
  output$chr_foot <-renderUI({
    HTML(
    paste("<p style='text-align:justify'>","Tracks show ChIP-seq enrichment for histone modifications across the gene locus. Signal heights correspond to bins per million mapped reads (BPM). Exons are depicted as grey boxes. LD: long day (7 weeks). SD: short day. LT: low temperature. WT: warm temperature. The number shown is the number of weeks in the treatment.", "</p>"
    ))
  })

  # gene <- annot[
  #   annot$type == "gene" &
  #     annot$ID == "Potra2n14c26778"
  # ]
  #
  # transcripts <- annot[annot$type == "mRNA" & grepl("Potra2n14c26778", annot$ID)]
  # exons <- annot[ annot$type %in% c("exon") & grepl("Potra2n14c26778", annot$ID)]
  # mcols(exons)$transcript <- mcols(exons)$Parent
  #
  # gene_region <- GRanges(
  #   seqnames = seqnames(gene),
  #   ranges = IRanges(
  #     start = start(gene) - 5000,
  #     end = end(gene) + 5000
  #   )
  # )
  #
  # gene_track <- GeneRegionTrack(
  #   range = exons,
  #   chromosome = seqnames(gene_region),
  #   name = mcols(gene)$ID,
  #   fill = "grey",
  #   col = "black",
  #   shape = "smallArrow",
  #   showId = TRUE,
  #   transcriptAnnotation = "transcript",
  #   collapseTranscripts = "longest",
  #   just.group = "below",
  #   geneSymbol = TRUE  ,
  #   background.title = "transparent",
  #   cex.group = 1,
  #   col.group = "black"
  # )
  #
  # genome_axis <- GenomeAxisTrack(col = "darkgrey",
  #                                fontcolor = "black",
  #                                labelPos = "below"
  # )
  #
  # plot_track <- function(input_file){
  #
  #   global_max <- max(
  #     sapply(rna, function(file) {
  #       data <- import(file, which = gene_region)
  #       max(data$score, na.rm = TRUE)
  #     }
  #     )
  #   )
  #
  #
  #   track <- lapply(1:5, function(i) {
  #     DataTrack(
  #       range = import(rna[i], which = gene_region),
  #       type = "histogram",
  #       name = sample_names[i],
  #       col.histogram = colors[i],
  #       fill.histogram = colors[i],
  #       ylim = c(0, global_max),
  #       window = 150,
  #       showAxis = TRUE,
  #       col.axis = "black",
  #       background.title = "transparent",
  #       col.title = "black",
  #       cex.axis = 0.8   ,
  #       cex.title = 1
  #     )
  #   })
  #
  #   plotTracks(
  #     c(track, list(gene_track, genome_axis)),
  #     from = start(gene_region),
  #     to = end(gene_region),
  #     chromosome = seqnames(gene_region),
  #     sizes = c(rep(0.8, 5), 0.5, 0.5)
  #   )
  # }


}

# Run the application
shinyApp(ui = ui, server = server)

