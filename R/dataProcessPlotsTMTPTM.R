#' Visualization for explanatory data analysis - TMT experiment
#'
#' To illustrate the quantitative data and quality control of MS runs,
#' dataProcessPlotsTMT takes the quantitative data from MSstatsTMT converter functions as input
#' and generate two types of figures in pdf files as output :
#' (1) profile plot (specify "ProfilePlot" in option type), to identify the potential sources of variation for each protein;
#' (2) quality control plot (specify "QCPlot" in option type), to evaluate the systematic bias between MS runs.
#'
#' @export
#' @import ggplot2
#' @importFrom graphics axis image legend mtext par plot.new title plot
#' @importFrom grDevices dev.off hcl pdf
#' @importFrom dplyr mutate n
#' @importFrom reshape2 dcast
#' @importFrom MSstatsTMT proteinSummarization
#' @importFrom gridExtra grid.arrange
#' @param data.ptm name of the data with PTM sites in protein name, which can be the output of MSstatsTMT converter functions.
#' @param data.protein name of the data with peptide level, which can be the output of MSstatsTMT converter functions.
#' @param data.ptm.summarization name of the data with ptm sites in protein-level name, which can be the output of
#' the MSstatsTMT \code{\link{proteinSummarization}} function.
#' @param data.protein.summarization name of the data with protein-level, which can be the output of the
#' MSstatsTMT \code{\link{proteinSummarization}} function.
#' @param type choice of visualization. "ProfilePlot" represents profile plot of log intensities across MS runs.
#' "QCPlot" represents box plots of log intensities across channels and MS runs.
#' @param ylimUp upper limit for y-axis in the log scale.
#' FALSE(Default) for Profile Plot and QC Plot uses the upper limit as rounded off maximum of log2(intensities) after normalization + 3..
#' @param ylimDown lower limit for y-axis in the log scale. FALSE(Default) for Profile Plot and QC Plot uses 0..
#' @param x.axis.size size of x-axis labeling for "Run" and "channel in Profile Plot and QC Plot.
#' @param y.axis.size size of y-axis labels. Default is 10.
#' @param text.size size of labels represented each condition at the top of Profile plot and QC plot. Default is 4.
#' @param text.angle angle of labels represented each condition at the top of Profile plot and QC plot. Default is 0.
#' @param legend.size size of legend above Profile plot. Default is 7.
#' @param dot.size.profile size of dots in Profile plot. Default is 2.
#' @param ncol.guide number of columns for legends at the top of plot. Default is 5.
#' @param width width of the saved pdf file. Default is 10.
#' @param height height of the saved pdf file. Default is 10.
#' @param which.Protein Protein list to draw plots. List can be names of Proteins or order numbers of Proteins.
#' Default is "all", which generates all plots for each protein. For QC plot, "allonly" will generate one QC plot with all proteins.
#' @param originalPlot TRUE(default) draws original profile plots, without normalization.
#' @param summaryPlot TRUE(default) draws profile plots with protein summarization for each channel and MS run.
#' @param address the name of folder that will store the results. Default folder is the current working directory.
#' The other assigned folder has to be existed under the current working directory.
#' An output pdf file is automatically created with the default name of "ProfilePlot.pdf" or "QCplot.pdf".
#' The command address can help to specify where to store the file as well as how to modify the beginning of the file name.
#' If address=FALSE, plot will be not saved as pdf file but showed in window.
#' @return plot or pdf
#' @examples
#' data(raw.ptm)
#' data(raw.protein)
#' data(quant.msstats.ptm)
#' data(quant.msstats.protein)
#'
#' ## Profile plot
#' dataProcessPlotsTMTPTM(data.ptm=raw.ptm,
#'                    data.protein=raw.protein,
#'                    data.ptm.summarization=quant.msstats.ptm,
#'                    data.protein.summarization=quant.msstats.protein,
#'                    type='ProfilePlot')
#'
#' ## NottoRun: QC plot
#' # dataProcessPlotsTMTPTM(data.ptm=raw.ptm,
#'                     # data.protein=raw.protein,
#'                     # data.ptm.summarization=quant.msstats.ptm,
#'                     # data.protein.summarization=quant.msstats.protein,
#'                     # type='QCPlot')
dataProcessPlotsTMTPTM <- function(data.ptm,
                                   data.protein,
                                   data.ptm.summarization,
                                   data.protein.summarization,
                                   type,
                                   ylimUp = FALSE,
                                   ylimDown = FALSE,
                                   x.axis.size = 10,
                                   y.axis.size = 10,
                                   text.size = 4,
                                   text.angle = 90,
                                   legend.size = 7,
                                   dot.size.profile = 2,
                                   ncol.guide = 5,
                                   width = 10,
                                   height = 12,
                                   which.Protein = "all",
                                   originalPlot = TRUE,
                                   summaryPlot = TRUE,
                                   address = "") {

  ## save process output in each step
  allfiles <- list.files()
  filenaming <- "msstatstmtptm"

  if (length(grep(filenaming,allfiles)) == 0) {

    finalfile <- "msstatstmtptm.log"
    processout <- NULL

  } else {

    num <- 0
    finalfile <- "msstatstmtptm.log"

    while (is.element(finalfile, allfiles)) {
      num <- num + 1
      lastfilename <- finalfile ## in order to rea
      finalfile <- paste0(paste(filenaming, num, sep="-"), ".log")
    }

    finalfile <- lastfilename
    processout <- as.matrix(read.table(finalfile, header=TRUE, sep="\t"))
  }

  processout <- rbind(processout,
                      as.matrix(c(" ", " ", "MSstatsTMTPTM - dataProcessPlotsTMTPTM function", " "), ncol=1))

  ## Checking for input variables
  type <- toupper(type)

  if (length(setdiff(type, c("PROFILEPLOT", "QCPLOT"))) != 0) {

    processout <- rbind(processout,
                        c(paste0("Input for type=", type,
                                 ". However,'type' should be one of ProfilePlot, QCPlot.")))
    write.table(processout, file=finalfile, row.names=FALSE)

    stop(paste0("Input for type=", type,
                ". However,'type' should be one of ProfilePlot, QCPlot."))
  }


  Condition = Run = xorder = Channel = NULL
  PeptideSequence = PSM = ProteinName = NULL
  GlobalProtein = Protein = NULL
  groupAxis = cumGroupAxis = abundance = analysis = NULL

  datafeature.protein <- data.protein
  datafeature.ptm <- data.ptm
  datarun.protein <- data.protein.summarization
  datarun.ptm <- data.ptm.summarization

  # conditions in feature data
  fea.conds.protein <- as.character(unique(datafeature.protein$Condition))
  fea.conds.ptm <- as.character(unique(datafeature.ptm$Condition))
  # conditions in protein data
  run.conds.protein <- as.character(unique(datarun.protein$Condition))
  run.conds.ptm <- as.character(unique(datarun.ptm$Condition))

  # only keep the overlapped conditions between feature data and protein data
  shared.conds <- Reduce(intersect, list(fea.conds.protein, fea.conds.ptm, run.conds.protein, run.conds.ptm))
  datafeature.protein <- datafeature.protein[datafeature.protein$Condition %in% shared.conds,]
  datafeature.ptm <- datafeature.ptm[datafeature.ptm$Condition %in% shared.conds,]
  datarun.protein <- datarun.protein[datarun.protein$Condition %in% shared.conds,]
  datarun.ptm <- datarun.ptm[datarun.ptm$Condition %in% shared.conds,]

  # make sure condition is factor
  datafeature.protein$Condition <- factor(datafeature.protein$Condition)
  datafeature.ptm$Condition <- factor(datafeature.ptm$Condition)
  datarun.protein$Condition <- factor(datarun.protein$Condition)
  datarun.ptm$Condition <- factor(datarun.ptm$Condition)

  ## Remove Site from protein name
  regex_protein <- '([^-]+)(?:_[^-]+){1}$'
  datafeature.ptm <- datafeature.ptm %>% mutate(GlobalProtein = str_match(ProteinName, regex_protein)[,2])

  colnames(datafeature.protein)[colnames(datafeature.protein) == 'ProteinName'] <- 'Protein'
  colnames(datafeature.ptm)[colnames(datafeature.ptm) == 'ProteinName'] <- 'Protein'

  datafeature.protein$Protein <- factor(datafeature.protein$Protein)
  datafeature.ptm$Protein <- factor(datafeature.ptm$Protein)
  datafeature.ptm$GlobalProtein <- factor(datafeature.ptm$GlobalProtein)

  datarun.protein$Protein <- factor(datarun.protein$Protein)
  datarun.ptm$Protein <- factor(datarun.ptm$Protein)

  ## feature level data : log2 transform
  datafeature.protein$abundance <- log2(datafeature.protein$Intensity)
  datafeature.ptm$abundance <- log2(datafeature.ptm$Intensity)
  datafeature.protein[!is.na(datafeature.protein$Intensity) &
                        datafeature.protein$Intensity < 1, 'abundance'] <- 0
  datarun.ptm[!is.na(datarun.ptm$Intensity) &
                datarun.ptm$Intensity < 1, 'abundance'] <- 0

  if (length(setdiff(toupper(type), c(toupper("ProfilePlot"), toupper("QCPlot")))) != 0) {
    stop(paste0("Input for type=", type,
                ". However,'type' should be one of \"ProfilePlot\", \"QCPlot\"."))
  }

  if (address == FALSE){
    ## here I used == FALSE, instead of !address. Because address can be logical or characters.
    if (which.Protein == 'all') {
      stop('** Cannnot generate all plots in a screen. Please set one protein at a time.')
    } else if (length(which.Protein) > 1) {
      stop('** Cannnot generate multiple plots in a screen. Please set one protein at a time.')
    }
  }

  ## Profile plot ##
  ## ---------------
  if (toupper(type) == "PROFILEPLOT") {

    processout <- rbind(processout,
                        c("ProfilePlot plotting started."))

    ## choose Proteins or not
    if (which.Protein != "all") {
      ## check which.Protein is name of Protein
      if (is.character(which.Protein)) {
        temp.name <- which.Protein

        ## message if name of Protein is wrong.
        if (length(setdiff(temp.name,unique(datafeature.ptm$Protein))) > 0) {
          stop(paste0("Please check protein name. Data set does not have this protein. - ",
                      toString(temp.name)))
        }
      }

      ## check which.Protein is order number of Protein
      if (is.numeric(which.Protein)) {
        temp.name <- levels(datafeature.ptm$Protein)[which.Protein]

        ## message if name of Protein is wrong.
        if (length(levels(datafeature.ptm$Protein)) < max(which.Protein)) {
          stop(paste0("Please check your ion of proteins. There are ",
                      length(levels(datafeature.ptm$Protein))," proteins in this dataset."))
        }
      }

      ## use only assigned proteins
      datafeature.ptm <- datafeature.ptm[which(datafeature.ptm$Protein %in% temp.name), ]
      temp_proteins <- as.character((datafeature.ptm %>% distinct(GlobalProtein))[[1]])
      datafeature.ptm$Protein <- factor(datafeature.ptm$Protein)

      datafeature.protein <- datafeature.protein[which(datafeature.protein$Protein %in% temp_proteins), ]
      datafeature.protein$Protein <- factor(datafeature.protein$Protein)

      datarun.protein <- datarun.protein[which(datarun.protein$Protein %in% temp_proteins), ]
      datarun.ptm <- datarun.ptm[which(datarun.ptm$Protein %in% temp.name), ]
      datarun.protein$Protein <- factor(datarun.protein$Protein)
      datarun.ptm$Protein <- factor(datarun.ptm$Protein)
    }

    ## assign upper or lower limit
    y.limup <- ceiling(max(datafeature.protein$abundance, datafeature.ptm$abundance, na.rm = TRUE) + 5)

    if (is.numeric(ylimUp)) {
      y.limup <- ylimUp
    }

    y.limdown <- 0
    if (is.numeric(ylimDown)) {
      y.limdown <- ylimDown
    }

    datafeature.protein <- datafeature.protein[with(datafeature.protein, order(Run, Condition, Channel)), ]
    datafeature.ptm <- datafeature.ptm[with(datafeature.ptm, order(Run, Condition, Channel)), ]
    datafeature.protein$Run <- factor(datafeature.protein$Run)
    datafeature.ptm$Run <- factor(datafeature.ptm$Run)
    datarun.protein$Run <- factor(datarun.protein$Run)
    datarun.ptm$Run <- factor(datarun.ptm$Run)

    ## !! important: order of x-axis
    ## can be reorder by group and then channel, WITHIN Run
    ## first make new column for x-axis
    datafeature.protein$group.channel <- paste(datafeature.protein$Condition, datafeature.protein$Channel, sep = "_")
    datafeature.ptm$group.channel <- paste(datafeature.ptm$Condition, datafeature.ptm$Channel, sep = "_")

    ## not sure better way for coding
    ## potentially change it.
    datafeature.protein$xorder <- NA
    datafeature.ptm$xorder <- NA

    for (k in seq_along(unique(datafeature.protein$Run))) {

      runid <- unique(datafeature.protein$Run)[k]
      datafeature.protein[datafeature.protein$Run == runid, ]$xorder <- factor(datafeature.protein[datafeature.protein$Run == runid, ]$group.channel,
                                                               levels <- unique(datafeature.protein[datafeature.protein$Run == runid, ]$group.channel),
                                                               labels <- seq(1, length(unique(datafeature.protein[datafeature.protein$Run == runid, ]$group.channel))))
    }

    for (k in seq_along(unique(datafeature.ptm$Run))) {

      runid <- unique(datafeature.ptm$Run)[k]
      datafeature.ptm[datafeature.ptm$Run == runid, ]$xorder <- factor(datafeature.ptm[datafeature.ptm$Run == runid, ]$group.channel,
                                                                               levels <- unique(datafeature.ptm[datafeature.ptm$Run == runid, ]$group.channel),
                                                                               labels <- seq(1, length(unique(datafeature.ptm[datafeature.ptm$Run == runid, ]$group.channel))))
    }


    ## check
    ## unique(datafeature[datafeature$Run == '5', c('Channel', 'Condition', 'Run', 'xorder','group.channel')])

    ## need to make data.frame with same variables for condition name
    datafeature.protein$xorder <- as.numeric(datafeature.protein$xorder)
    datafeature.ptm$xorder <- as.numeric(datafeature.ptm$xorder)

    ## keep unique information for x-axis labeling. will be used in plotting
    tempGroupName.protein <- unique(datafeature.protein[, c("Condition", "xorder", "Run", "Channel")])
    tempGroupName.ptm <- unique(datafeature.ptm[, c("Condition", "xorder", "Run", "Channel")])

    groupline.protein <- tempGroupName.protein %>% dplyr::group_by(Condition, Run) %>% dplyr::mutate(groupAxis = n())
    groupline.protein <- groupline.protein %>% dplyr::select(-xorder, -Channel)
    groupline.protein <- groupline.protein[!duplicated(groupline.protein), ]

    groupline.ptm <- tempGroupName.ptm %>% dplyr::group_by(Condition, Run) %>% dplyr::mutate(groupAxis = n())
    groupline.ptm <- groupline.ptm %>% dplyr::select(-xorder, -Channel)
    groupline.ptm <- groupline.ptm[!duplicated(groupline.ptm), ]


    groupline.protein <- groupline.protein %>% dplyr::group_by(Run) %>% dplyr::mutate(cumGroupAxis = cumsum(groupAxis))
    groupline.ptm <- groupline.ptm %>% dplyr::group_by(Run) %>% dplyr::mutate(cumGroupAxis = cumsum(groupAxis))

    groupline.protein$cumGroupAxis <- groupline.protein$cumGroupAxis + 0.5
    groupline.ptm$cumGroupAxis <- groupline.ptm$cumGroupAxis + 0.5

    ## add coordinate for group id
    groupline.protein$xorder <- groupline.protein$cumGroupAxis - groupline.protein$groupAxis / 2
    groupline.protein$abundance <- y.limup - 0.5
    groupline.ptm$xorder <- groupline.ptm$cumGroupAxis - groupline.ptm$groupAxis / 2
    groupline.ptm$abundance <- y.limup - 0.5

    groupline.all.protein <- groupline.protein
    groupline.all.ptm <- groupline.ptm

    ## remove last condition for vertical line between groups
    groupline.protein <- groupline.protein[-which(groupline.protein$Condition %in%
                                                    levels(groupline.protein$Condition)[nlevels(groupline.protein$Condition)]), ]
    groupline.ptm <- groupline.ptm[-which(groupline.ptm$Condition %in% levels(groupline.ptm$Condition)[nlevels(groupline.ptm$Condition)]), ]

    ## need to fill in incomplete rows for Runlevel data
    haverun <- FALSE

    if (sum(is.element(colnames(datarun.protein), "Run")) != 0) {
      datamat <- reshape2::dcast(Protein + Channel ~ Run, data = datarun.protein, value.var = 'Abundance', keep = TRUE)

      datarun.protein <- reshape2::melt(datamat, id.vars=c('Protein', 'Channel'))
      colnames(datarun.protein)[colnames(datarun.protein) %in% c("variable", "value")] <- c('Run', 'Abundance')

      ## match x axis order
      datarun.protein <- merge(datarun.protein, tempGroupName.protein, by = c('Run', 'Channel'))

      haverun <- TRUE
    }

    if (sum(is.element(colnames(datarun.ptm), "Run")) != 0) {
      datamat <- reshape2::dcast(Protein + Channel ~ Run, data = datarun.ptm, value.var = 'Abundance', keep = TRUE)

      datarun.ptm <- reshape2::melt(datamat, id.vars=c('Protein', 'Channel'))
      colnames(datarun.ptm)[colnames(datarun.ptm) %in% c("variable", "value")] <- c('Run', 'Abundance')

      ## match x axis order
      datarun.ptm <- merge(datarun.ptm, tempGroupName.ptm, by = c('Run', 'Channel'))

      haverun <- TRUE
    }


    ## save the plots as pdf or not
    ## If there are the file with the same name, add next numbering at the end of file name

    ## y-axis labeling
    yaxis.name <- 'Log2-intensities'

    ## Only plot proteins that occur in both datasets
    global_proteins <- (datafeature.protein %>% distinct(Protein))[[1]]
    ptm_proteins <- (datafeature.ptm %>% distinct(GlobalProtein))[[1]]
    plot_proteins <- intersect(ptm_proteins, global_proteins)

    datafeature.ptm <- datafeature.ptm %>% filter(GlobalProtein %in% plot_proteins)

    ## TODO: If protein included in which.protein is not available in both datasets, print that here:
    plot_proteins <- (datafeature.ptm %>% distinct(Protein))[[1]]

    if (originalPlot) {
      if (address != FALSE) {
        allfiles <- list.files()

        num <- 0
        filenaming <- paste0(address, "ProfilePlot")
        finalfile <- paste0(address, "ProfilePlot.pdf")

        while (is.element(finalfile, allfiles)) {
          num <- num + 1
          finalfile <- paste0(paste(filenaming, num, sep = "-"), ".pdf")
        }

        pdf(finalfile, width = width, height = height)
      }

      ## factoring for run, channel, condition should be done before loop

      for (i in 1:length(plot_proteins)) {

        sub.ptm <- datafeature.ptm[datafeature.ptm$Protein == as.character(plot_proteins[i]), ]
        sub.protein <- datafeature.protein[datafeature.protein$Protein == as.character((sub.ptm %>% distinct(GlobalProtein))[[1]]), ]

        sub.protein$PeptideSequence <- factor(as.character(sub.protein$PeptideSequence))
        sub.protein$Charge <- factor(as.character(sub.protein$Charge))
        sub.protein$PSM <- factor(as.character(sub.protein$PSM))

        sub.ptm$PeptideSequence <- factor(as.character(sub.ptm$PeptideSequence))
        sub.ptm$Charge <- factor(as.character(sub.ptm$Charge))
        sub.ptm$PSM <- factor(as.character(sub.ptm$PSM))

        # if all measurements are NA,
        if (nrow(sub.protein) == sum(is.na(sub.protein$abundance))) {
          message(paste0("Can't the Profile plot for ", unique(sub.protein$Protein),
                         "(", i, " of ", length(plot_proteins),
                         ") because all measurements are NAs."))
          next()
        }

        if (nrow(sub.protein) == sum(!is.na(sub.protein$abundance) & sub.protein$abundance == 0)) {
          message(paste0("Can't the Profile plot for ", unique(sub.protein$Protein),
                         "(", i, " of ", length(plot_proteins),
                         ") because all measurements are zeros."))
          next()
        }
        # if all measurements are NA,
        if (nrow(sub.ptm) == sum(is.na(sub.ptm$abundance))) {
          message(paste0("Can't the Profile plot for ", unique(sub.ptm$Protein),
                         "(", i, " of ", length(plot_proteins),
                         ") because all measurements are NAs."))
          next()
        }

        if (nrow(sub.ptm) == sum(!is.na(sub.ptm$abundance) & sub.ptm$abundance == 0)) {
          message(paste0("Can't the Profile plot for ", unique(sub.ptm$Protein),
                         "(", i, " of ", length(plot_proteins),
                         ") because all measurements are zeros."))
          next()
        }


        ## seq for peptide and charge
        ## for seting up color and linetype
        b.protein <- unique(sub.protein[, c("PeptideSequence", "PSM")])
        b.protein <- b.protein[with(b.protein, order(PeptideSequence, PSM)), ] ## add because if there are missing value, orders are different.

        b.ptm <- unique(sub.ptm[, c("PeptideSequence", "PSM")])
        b.ptm <- b.ptm[with(b.ptm, order(PeptideSequence, PSM)), ] ## add because if there are missing value, orders are different.

        temp1.protein <- xtabs(~PeptideSequence, b.protein)
        temp1.ptm <- xtabs(~PeptideSequence, b.ptm)

        ## unique charge id within peptide sequence, for line type
        ss.protein <- NULL
        ss.ptm <- NULL
        ## unique peptide sequence id, for color
        s.protein <- NULL
        s.ptm <- NULL

        for (j in seq_along(temp1.protein)) {
          temp3 <- rep(j, temp1.protein[j])
          s.protein <- c(s.protein, temp3)
          temp2 <- seq(1, temp1.protein[j])
          ss.protein <- c(ss.protein, temp2)
        }

        for (j in seq_along(temp1.ptm)) {
          temp3 <- rep(j, temp1.ptm[j])
          s.ptm <- c(s.ptm, temp3)
          temp2 <- seq(1, temp1.ptm[j])
          ss.ptm <- c(ss.ptm, temp2)
        }

        ## for annotation of condition
        groupline.tmp.protein <- data.frame(groupline.protein,
                                            "PSM" = unique(sub.protein$PSM)[1],
                                            "PeptideSequence" = unique(sub.protein$PeptideSequence)[1])

        groupline.all.tmp.protein <- data.frame(groupline.all.protein,
                                                "PSM" = unique(sub.protein$PSM)[1],
                                                "PeptideSequence" = unique(sub.protein$PeptideSequence)[1])

        groupline.tmp.ptm <- data.frame(groupline.ptm,
                                        "PSM" = unique(sub.ptm$PSM)[1],
                                        "PeptideSequence" = unique(sub.ptm$PeptideSequence)[1])

        groupline.all.tmp.ptm <- data.frame(groupline.all.ptm,
                                            "PSM" = unique(sub.ptm$PSM)[1],
                                            "PeptideSequence" = unique(sub.ptm$PeptideSequence)[1])

        ## 2019. 12. 17, MC : for profile plot, define color for dot
        cbp <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        check.length.protein <- length(unique(s.protein)) %/% length(cbp)
        if ( check.length.protein > 0 ){
          cbp.protein <- rep(cbp, times=check.length.protein + 1)
        } else {
          cbp.protein <- cbp
        }
        check.length.ptm <- length(unique(s.ptm)) %/% length(cbp)
        if ( check.length.ptm > 0 ){
          cbp.ptm <- rep(cbp, times=check.length.ptm + 1)
        } else {
          cbp.ptm <- cbp
        }
        ##

        ## 1st plot for Protein plot
        protein_temp <- ggplot(aes_string(x = 'xorder', y = 'abundance',
                                          color = 'PSM', linetype = 'PSM'), data = sub.protein) +
          facet_grid(~Run) +
          geom_point(size=dot.size.profile) +
          geom_line(size = 0.5) +
          scale_colour_manual(values=cbp[s.protein]) +
          scale_linetype_manual(values = ss.protein) +
          scale_shape_manual(values = c(16)) +
          labs(title = paste0('Protein - ', unique(sub.protein$Protein)),
               x = 'MS runs') +
          scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
          scale_x_continuous('MS runs') +
          geom_vline(data = groupline.tmp.protein,
                     aes(xintercept = cumGroupAxis),
                     colour = "grey", linetype = "longdash") +
          geom_text(data = groupline.all.tmp.protein,
                    aes(x = xorder, y = abundance, label = Condition),
                    size = text.size,
                    angle = text.angle, hjust = .9,
                    color = "black") +
          theme(
            panel.background = element_rect(fill = 'white', colour = "black"),
            legend.key = element_rect(fill = 'white', colour = 'white'),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = 'gray95'),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = y.axis.size, colour = "black"),
            axis.ticks = element_line(colour = "black"),
            axis.title.x = element_text(size = x.axis.size + 5, vjust = -0.4),
            axis.title.y = element_text(size = y.axis.size + 5, vjust = 0.3),
            title = element_text(size = x.axis.size + 8, vjust = 1.5),
            legend.position = "top",
            legend.text = element_text(size = legend.size)) +
          guides(color = guide_legend(title = paste("# peptide:", nlevels(sub.protein$PeptideSequence)),
                                      title.theme = element_text(size = 13, angle = 0),
                                      keywidth = 0.4,
                                      keyheight = 0.1,
                                      default.unit = 'inch',
                                      ncol = ncol.guide),
                 linetype = guide_legend(title = paste("# peptide:", nlevels(sub.protein$PeptideSequence)),
                                         title.theme = element_text(size = 13, angle = 0),
                                         keywidth = 0.4,
                                         keyheight = 0.1,
                                         default.unit = 'inch',
                                         ncol = ncol.guide))

        ## 1st plot for PTM plot
        ptm_temp <- ggplot(aes_string(x = 'xorder', y = 'abundance',
                                      color = 'PSM', linetype = 'PSM'), data = sub.ptm) +
          facet_grid(~Run) +
          geom_point(size=dot.size.profile) +
          geom_line(size = 0.5) +
          scale_colour_manual(values=cbp[s.ptm]) +
          scale_linetype_manual(values = ss.ptm) +
          scale_shape_manual(values = c(16)) +
          labs(title = paste0('PTM - ', unique(sub.ptm$Protein)),
               x = 'MS runs') +
          scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
          scale_x_continuous('MS runs') +
          geom_vline(data = groupline.tmp.ptm,
                     aes(xintercept = cumGroupAxis),
                     colour = "grey", linetype = "longdash") +
          geom_text(data = groupline.all.tmp.ptm,
                    aes(x = xorder, y = abundance, label = Condition),
                    size = text.size,
                    angle = text.angle, hjust = .9,
                    color = "black") +
          theme(
            panel.background = element_rect(fill = 'white', colour = "black"),
            legend.key = element_rect(fill = 'white', colour = 'white'),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = 'gray95'),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = y.axis.size, colour = "black"),
            axis.ticks = element_line(colour = "black"),
            axis.title.x = element_text(size = x.axis.size + 5, vjust = -0.4),
            axis.title.y = element_text(size = y.axis.size + 5, vjust = 0.3),
            title = element_text(size = x.axis.size + 8, vjust = 1.5),
            legend.position = "top",
            legend.text = element_text(size = legend.size)) +
          guides(color = guide_legend(title = paste("# peptide:", nlevels(sub.ptm$PeptideSequence)),
                                      title.theme = element_text(size = 13, angle = 0),
                                      keywidth = 0.4,
                                      keyheight = 0.1,
                                      default.unit = 'inch',
                                      ncol = ncol.guide),
                 linetype = guide_legend(title = paste("# peptide:", nlevels(sub.ptm$PeptideSequence)),
                                         title.theme = element_text(size = 13, angle = 0),
                                         keywidth = 0.4,
                                         keyheight = 0.1,
                                         default.unit = 'inch',
                                         ncol = ncol.guide))

        gridExtra::grid.arrange(ptm_temp, protein_temp, ncol=1)

        message(paste("Drew the Profile plot for ", unique(sub.ptm$Protein),
                      "(", i, " of ", length(plot_proteins), ")"))
      }
      # end-loop for each protein

      if (address != FALSE) {
        dev.off()
      }

    } # end original plot

    ############################################
    ## 2st plot for original plot : summary
    ############################################

    if (summaryPlot) {
      if (address != FALSE) {
        allfiles <- list.files()

        num <- 0
        filenaming <- paste0(address, "ProfilePlot_wSummarization")
        finalfile <- paste0(address, "ProfilePlot_wSummarization.pdf")

        while (is.element(finalfile, allfiles)) {
          num <- num + 1
          finalfile <- paste0(paste(filenaming, num, sep = "-"), ".pdf")
        }

        pdf(finalfile, width = width, height = height)
      }

      for (i in 1:length(plot_proteins)) {

        sub.ptm <- datafeature.ptm[datafeature.ptm$Protein == as.character(plot_proteins[i]), ]
        sub.protein <- datafeature.protein[datafeature.protein$Protein == as.character((sub.ptm %>% distinct(GlobalProtein))[[1]]), ]

        sub.protein$PeptideSequence <- factor(as.character(sub.protein$PeptideSequence))
        sub.protein$Charge <- factor(as.character(sub.protein$Charge))
        sub.protein$PSM <- factor(as.character(sub.protein$PSM))

        sub.ptm$PeptideSequence <- factor(as.character(sub.ptm$PeptideSequence))
        sub.ptm$Charge <- factor(as.character(sub.ptm$Charge))
        sub.ptm$PSM <- factor(as.character(sub.ptm$PSM))

        # if all measurements are NA,
        if (nrow(sub.protein) == sum(is.na(sub.protein$abundance))) {
          message(paste0("Can't the Profile plot for ", unique(sub.protein$Protein),
                         "(", i, " of ", length(plot_proteins),
                         ") because all measurements are NAs."))
          next()
        }

        if (nrow(sub.protein) == sum(!is.na(sub.protein$abundance) & sub.protein$abundance == 0)) {
          message(paste0("Can't the Profile plot for ", unique(sub.protein$Protein),
                         "(", i, " of ", length(plot_proteins),
                         ") because all measurements are zeros."))
          next()
        }
        # if all measurements are NA,
        if (nrow(sub.ptm) == sum(is.na(sub.ptm$abundance))) {
          message(paste0("Can't the Profile plot for ", unique(sub.ptm$Protein),
                         "(", i, " of ", length(plot_proteins),
                         ") because all measurements are NAs."))
          next()
        }

        if (nrow(sub.ptm) == sum(!is.na(sub.ptm$abundance) & sub.ptm$abundance == 0)) {
          message(paste0("Can't the Profile plot for ", unique(sub.ptm$Protein),
                         "(", i, " of ", length(plot_proteins),
                         ") because all measurements are zeros."))
          next()
        }

        ## for annotation of condition
        groupline.tmp.protein <- data.frame(groupline.protein,
                                            "PSM" = unique(sub.protein$PSM)[1],
                                            "PeptideSequence" = unique(sub.protein$PeptideSequence)[1])

        groupline.all.tmp.protein <- data.frame(groupline.all.protein,
                                                "PSM" = unique(sub.protein$PSM)[1],
                                                "PeptideSequence" = unique(sub.protein$PeptideSequence)[1])

        groupline.tmp.ptm <- data.frame(groupline.ptm,
                                        "PSM" = unique(sub.ptm$PSM)[1],
                                        "PeptideSequence" = unique(sub.ptm$PeptideSequence)[1])

        groupline.all.tmp.ptm <- data.frame(groupline.all.ptm,
                                            "PSM" = unique(sub.ptm$PSM)[1],
                                            "PeptideSequence" = unique(sub.ptm$PeptideSequence)[1])

        if (haverun) {
          subrun.ptm <- datarun.ptm[datarun.ptm$Protein == as.character(plot_proteins[i]), ]
          subrun.protein <- datarun.protein[datarun.protein$Protein == as.character((sub.ptm %>% distinct(GlobalProtein))[[1]]), ]

          if (nrow(subrun.protein) != 0) {

            quantrun.ptm <- sub.ptm[1, ]
            quantrun.ptm[, 2:ncol(quantrun.ptm)] <- NA
            quantrun.ptm <- quantrun.ptm[rep(seq_len(nrow(subrun.ptm))), ]

            quantrun.ptm$Protein <- subrun.ptm$Protein
            quantrun.ptm$PeptideSequence <- "Run summary"
            quantrun.ptm$Charge <- "Run summary"
            quantrun.ptm$PSM <- "Run summary"
            quantrun.ptm$Channel <- subrun.ptm$Channel
            quantrun.ptm$Run <- subrun.ptm$Run
            quantrun.ptm$abundance <- subrun.ptm$Abundance
            quantrun.ptm$xorder <- subrun.ptm$xorder

          } else {
            # if there is only one Run measured across all runs, no Run information for linear with censored
            quantrun.ptm <- datafeature.ptm[1, ]
            quantrun.ptm[, 2:ncol(quantrun.ptm)] <- NA

            quantrun.ptm$Protein <- levels(datafeature.ptm$Protein)[i]
            quantrun.ptm$PeptideSequence <- "Run summary"
            quantrun.ptm$Charge <- "Run summary"
            quantrun.ptm$PSM <- "Run summary"
            quantrun.ptm$abundance <- NA
            quantrun.ptm$Intensity <- NA
          }

          if (nrow(subrun.protein) != 0) {

            quantrun.protein <- sub.protein[1, ]
            quantrun.protein[, 2:ncol(quantrun.protein)] <- NA
            quantrun.protein <- quantrun.protein[rep(seq_len(nrow(subrun.protein))), ]

            quantrun.protein$Protein <- subrun.protein$Protein
            quantrun.protein$PeptideSequence <- "Run summary"
            quantrun.protein$Charge <- "Run summary"
            quantrun.protein$PSM <- "Run summary"
            quantrun.protein$Channel <- subrun.protein$Channel
            quantrun.protein$Run <- subrun.protein$Run
            quantrun.protein$abundance <- subrun.protein$Abundance
            quantrun.protein$xorder <- subrun.protein$xorder

          } else {
            # if there is only one Run measured across all runs, no Run information for linear with censored
            quantrun.protein <- datafeature.protein[1, ]
            quantrun.protein[, 2:ncol(quantrun.protein)] <- NA

            quantrun.protein$Protein <- levels(datafeature.protein$Protein)[i]
            quantrun.protein$PeptideSequence <- "Run summary"
            quantrun.protein$Charge <- "Run summary"
            quantrun.protein$PSM <- "Run summary"
            quantrun.protein$abundance <- NA
            quantrun.protein$Intensity <- NA
          }

          quantrun.ptm$analysis <- "Run summary"
          quantrun.protein$analysis <- "Run summary"
          sub.ptm$analysis <- "Processed feature-level data"
          sub.protein$analysis <- "Processed feature-level data"

          final.ptm <- rbind(sub.ptm, quantrun.ptm)
          final.ptm$analysis <- factor(final.ptm$analysis)
          final.ptm$PSM <- factor(final.ptm$PSM)

          final.protein <- rbind(sub.protein, quantrun.protein)
          final.protein$analysis <- factor(final.protein$analysis)
          final.protein$PSM <- factor(final.protein$PSM)

          ## Draw summarized ptm plot
          ptempall.ptm <- ggplot(aes_string(x = 'xorder', y = 'abundance',
                                            color = 'analysis', linetype = 'PSM', size = 'analysis'), data = final.ptm) +
            facet_grid(~Run) +
            geom_point(size = dot.size.profile) +
            geom_line(size = 0.5) +
            scale_colour_manual(values = c("lightgray", "darkred")) +
            scale_shape_manual(values = c(16)) +
            scale_size_manual(values = c(1.7, 2), guide = "none") +
            scale_linetype_manual(values = c(rep(1, times = length(unique(final.ptm$PSM))-1), 2), guide = "none") +
            labs(title = paste0('PTM - ', unique(sub.ptm$Protein)),
                 x = 'MS runs') +
            scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
            geom_vline(data = groupline.tmp.ptm,
                       aes(xintercept = cumGroupAxis),
                       colour = "grey", linetype = "longdash") +
            geom_text(data = groupline.all.tmp.ptm,
                      aes(x = xorder, y = abundance, label = Condition),
                      size = text.size,
                      angle = text.angle, hjust = .9,
                      color = "black") +
            theme(
              panel.background = element_rect(fill = 'white', colour = "black"),
              legend.key = element_rect(fill = 'white', colour = 'white'),
              panel.grid.minor = element_blank(),
              strip.background = element_rect(fill = 'gray95'),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size = y.axis.size, colour = "black"),
              axis.ticks = element_line(colour = "black"),
              axis.title.x = element_text(size = x.axis.size + 5, vjust = -0.4),
              axis.title.y = element_text(size = y.axis.size + 5, vjust = 0.3),
              title = element_text(size = x.axis.size + 8, vjust = 1.5),
              legend.position = "top",
              legend.text = element_text(size = legend.size),
              legend.title = element_blank()) +
            guides(color = guide_legend(order = 1,
                                        title = NULL,
                                        label.theme = element_text(size = 10, angle = 0)))

          ## draw point again because some red summary dots could be hiden
          ptempall.ptm <- ptempall.ptm + geom_point(data = final.ptm, aes(x = xorder, y = abundance, size = analysis, color = analysis))

          ## Draw summarized protein plot
          ptempall.protein <- ggplot(aes_string(x = 'xorder', y = 'abundance',
                                                color = 'analysis', linetype = 'PSM', size = 'analysis'), data = final.protein) +
            facet_grid(~Run) +
            geom_point(size = dot.size.profile) +
            geom_line(size = 0.5) +
            scale_colour_manual(values = c("lightgray", "darkred")) +
            scale_shape_manual(values = c(16)) +
            scale_size_manual(values = c(1.7, 2), guide = "none") +
            scale_linetype_manual(values = c(rep(1, times = length(unique(final.protein$PSM))-1), 2), guide = "none") +
            labs(title = paste0('Protein - ', unique(sub.protein$Protein)),
                 x = 'MS runs') +
            scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
            geom_vline(data = groupline.tmp.protein,
                       aes(xintercept = cumGroupAxis),
                       colour = "grey", linetype = "longdash") +
            geom_text(data = groupline.all.tmp.protein,
                      aes(x = xorder, y = abundance, label = Condition),
                      size = text.size,
                      angle = text.angle, hjust = .9,
                      color = "black") +
            theme(
              panel.background = element_rect(fill = 'white', colour = "black"),
              legend.key = element_rect(fill = 'white', colour = 'white'),
              panel.grid.minor = element_blank(),
              strip.background = element_rect(fill = 'gray95'),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size = y.axis.size, colour = "black"),
              axis.ticks = element_line(colour = "black"),
              axis.title.x = element_text(size = x.axis.size + 5, vjust = -0.4),
              axis.title.y = element_text(size = y.axis.size + 5, vjust = 0.3),
              title = element_text(size = x.axis.size + 8, vjust = 1.5),
              legend.position = "top",
              legend.text = element_text(size = legend.size),
              legend.title = element_blank()) +
            guides(color = guide_legend(order = 1,
                                        title = NULL,
                                        label.theme = element_text(size = 10, angle = 0)))

          ## draw point again because some red summary dots could be hiden
          ptempall.protein <- ptempall.protein + geom_point(data = final.protein, aes(x = xorder,
                                                                                      y = abundance, size = analysis, color = analysis))

          gridExtra::grid.arrange(ptempall.ptm, ptempall.protein, ncol=1)

          message(paste("Drew the Profile plot with summarization for ", unique(sub.ptm$Protein),
                        "(", i, " of ", length(unique(datafeature.ptm$Protein)), ")"))

        }

      } # end-loop for each protein

      if (address!=FALSE) {
        dev.off()
      }
    } # end summarization plot
  } # end Profile plot


  ## QC plot (Quality control plot) ##
  ## ---------------------------------
  if (toupper(type) == "QCPLOT") {

    ## y-axis labeling
    yaxis.name <- 'Log2-intensities'

    ## save the plots as pdf or not
    ## If there are the file with the same name, add next numbering at the end of file name
    if (address != FALSE) {
      allfiles <- list.files()

      num <- 0
      filenaming <- paste0(address,"QCPlot")
      finalfile <- paste0(address,"QCPlot.pdf")

      while (is.element(finalfile, allfiles)) {
        num <- num + 1
        finalfile <- paste0(paste(filenaming, num, sep = "-"), ".pdf")
      }

      pdf(finalfile, width = width, height = height)
    }

    ## assign upper or lower limit
    y.limup <- ceiling(max(datafeature.protein$abundance, datafeature.ptm$abundance, na.rm = TRUE) + 3)

    if (is.numeric(ylimUp)) {
      y.limup <- ylimUp
    }

    y.limdown <- 0
    if (is.numeric(ylimDown)) {
      y.limdown <- ylimDown
    }

    datafeature.protein <- datafeature.protein[with(datafeature.protein, order(Run, Condition, Channel)), ]
    datafeature.ptm <- datafeature.ptm[with(datafeature.ptm, order(Run, Condition, Channel)), ]
    datafeature.protein$Run <- factor(datafeature.protein$Run)
    datafeature.ptm$Run <- factor(datafeature.ptm$Run)
    datarun.protein$Run <- factor(datarun.protein$Run)
    datarun.ptm$Run <- factor(datarun.ptm$Run)

    ## !! important: order of x-axis
    ## can be reorder by group and then channel, WITHIN Run
    ## first make new column for x-axis
    datafeature.protein$group.channel <- paste(datafeature.protein$Condition, datafeature.protein$Channel, sep = "_")
    datafeature.ptm$group.channel <- paste(datafeature.ptm$Condition, datafeature.ptm$Channel, sep = "_")

    ## not sure better way for coding
    ## potentially change it.
    datafeature.protein$xorder <- NA
    datafeature.ptm$xorder <- NA

    for (k in seq_along(unique(datafeature.protein$Run))) {

      runid <- unique(datafeature.protein$Run)[k]
      datafeature.protein[datafeature.protein$Run == runid, ]$xorder <- factor(datafeature.protein[datafeature.protein$Run == runid, ]$group.channel,
                                                                               levels <- unique(datafeature.protein[datafeature.protein$Run == runid, ]$group.channel),
                                                                               labels <- seq(1, length(unique(datafeature.protein[datafeature.protein$Run == runid, ]$group.channel))))
    }

    for (k in seq_along(unique(datafeature.ptm$Run))) {

      runid <- unique(datafeature.ptm$Run)[k]
      datafeature.ptm[datafeature.ptm$Run == runid, ]$xorder <- factor(datafeature.ptm[datafeature.ptm$Run == runid, ]$group.channel,
                                                                       levels <- unique(datafeature.ptm[datafeature.ptm$Run == runid, ]$group.channel),
                                                                       labels <- seq(1, length(unique(datafeature.ptm[datafeature.ptm$Run == runid, ]$group.channel))))
    }


    ## check
    ## unique(datafeature[datafeature$Run == 'PAMI-176_Mouse_K-T', c('Channel', 'Condition', 'Run', 'xorder','group.channel')])

    ## need to make data.frame with same variables for condition name
    datafeature.protein$xorder <- as.numeric(datafeature.protein$xorder)
    datafeature.ptm$xorder <- as.numeric(datafeature.ptm$xorder)

    ## keep unique information for x-axis labeling. will be used in plotting
    tempGroupName.protein <- unique(datafeature.protein[, c("Condition", "xorder", "Run", "Channel")])
    tempGroupName.ptm <- unique(datafeature.ptm[, c("Condition", "xorder", "Run", "Channel")])

    ## count # per condition per Run
    #groupline <- unique(datafeature[, c('Condition', 'Run')])
    #groupline$groupAxis <- as.numeric(xtabs(~Condition+Run, tempGroupName))
    groupline.protein <- tempGroupName.protein %>% dplyr::group_by(Condition, Run) %>% dplyr::mutate(groupAxis = n())
    groupline.protein <- groupline.protein %>% dplyr::select(-xorder, -Channel)
    groupline.protein <- groupline.protein[!duplicated(groupline.protein), ]

    groupline.ptm <- tempGroupName.ptm %>% dplyr::group_by(Condition, Run) %>% dplyr::mutate(groupAxis = n())
    groupline.ptm <- groupline.ptm %>% dplyr::select(-xorder, -Channel)
    groupline.ptm <- groupline.ptm[!duplicated(groupline.ptm), ]

    ## make accumurated # as condition increase
    groupline.protein <- groupline.protein %>% dplyr::group_by(Run) %>% dplyr::mutate(cumGroupAxis = cumsum(groupAxis))
    groupline.ptm <- groupline.ptm %>% dplyr::group_by(Run) %>% dplyr::mutate(cumGroupAxis = cumsum(groupAxis))

    groupline.protein$cumGroupAxis <- groupline.protein$cumGroupAxis + 0.5
    groupline.ptm$cumGroupAxis <- groupline.ptm$cumGroupAxis + 0.5

    ## add coordinate for group id
    groupline.protein$xorder <- groupline.protein$cumGroupAxis - groupline.protein$groupAxis / 2
    groupline.protein$abundance <- y.limup - 0.5

    groupline.ptm$xorder <- groupline.ptm$cumGroupAxis - groupline.ptm$groupAxis / 2
    groupline.ptm$abundance <- y.limup - 0.5

    ## save all information, for labeling group in plot
    groupline.all.protein <- groupline.protein
    groupline.all.ptm <- groupline.ptm

    ## remove last condition for vertical line between groups
    groupline.protein <- groupline.protein[-which(groupline.protein$Condition %in% levels(groupline.protein$Condition)[nlevels(groupline.protein$Condition)]), ]
    groupline.ptm <- groupline.ptm[-which(groupline.ptm$Condition %in% levels(groupline.ptm$Condition)[nlevels(groupline.ptm$Condition)]), ]

    ## all protein
    if (which.Protein == 'all' | which.Protein == 'allonly') {

      ## for annotation of condition
      groupline.tmp.protein <- data.frame(groupline.protein,
                                  "PSM" = unique(datafeature.protein$PSM)[1],
                                  "PeptideSequence" = unique(datafeature.protein$PeptideSequence)[1])

      groupline.all.tmp.protein <- data.frame(groupline.all.protein,
                                      "PSM" = unique(datafeature.protein$PSM)[1],
                                      "PeptideSequence" = unique(datafeature.protein$PeptideSequence)[1])

      ## for annotation of condition
      groupline.tmp.ptm <- data.frame(groupline.ptm,
                                          "PSM" = unique(datafeature.ptm$PSM)[1],
                                          "PeptideSequence" = unique(datafeature.ptm$PeptideSequence)[1])

      groupline.all.tmp.ptm <- data.frame(groupline.all.ptm,
                                              "PSM" = unique(datafeature.ptm$PSM)[1],
                                              "PeptideSequence" = unique(datafeature.ptm$PeptideSequence)[1])

      ## 1st plot for original plot
      ## for boxplot, x-axis, xorder should be factor
      datafeature.protein$xorder <- factor(datafeature.protein$xorder)
      datafeature.ptm$xorder <- factor(datafeature.ptm$xorder)

      ptemp.ptm <- ggplot(aes_string(x = 'xorder', y = 'abundance'), data = datafeature.ptm) +
        facet_grid(~Run) +
        geom_boxplot(aes_string(fill = 'Condition'), outlier.shape = 1, outlier.size = 1.5) +
        labs(title = 'All PTMs',
             x = 'MS runs') +
        scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
        geom_vline(data = groupline.tmp.ptm,
                   aes(xintercept = cumGroupAxis),
                   colour = "grey", linetype = "longdash") +
        geom_text(data = groupline.all.tmp.ptm,
                  aes(x = xorder, y = abundance, label = Condition),
                  size = text.size,
                  angle = text.angle, hjust = .9,
                  color = "black") +
        theme(
          panel.background = element_rect(fill = 'white', colour = "black"),
          legend.key = element_rect(fill = 'white', colour = 'white'),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = 'gray95'),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = y.axis.size, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.title.x = element_text(size = x.axis.size + 5, vjust = -0.4),
          axis.title.y = element_text(size = y.axis.size + 5, vjust = 0.3),
          title = element_text(size = x.axis.size + 8, vjust = 1.5),
          legend.position = "none")

      ptemp.protein <- ggplot(aes_string(x = 'xorder', y = 'abundance'), data = datafeature.protein) +
        facet_grid(~Run) +
        geom_boxplot(aes_string(fill = 'Condition'), outlier.shape = 1, outlier.size = 1.5) +
        labs(title = 'All proteins',
             x = 'MS runs') +
        scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
        geom_vline(data = groupline.tmp.protein,
                   aes(xintercept = cumGroupAxis),
                   colour = "grey", linetype = "longdash") +
        geom_text(data = groupline.all.tmp.protein,
                  aes(x = xorder, y = abundance, label = Condition),
                  size = text.size,
                  angle = text.angle, hjust = .9,
                  color = "black") +
        theme(
          panel.background = element_rect(fill = 'white', colour = "black"),
          legend.key = element_rect(fill = 'white', colour = 'white'),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = 'gray95'),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = y.axis.size, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.title.x = element_text(size = x.axis.size + 5, vjust = -0.4),
          axis.title.y = element_text(size = y.axis.size + 5, vjust = 0.3),
          title = element_text(size = x.axis.size + 8, vjust = 1.5),
          legend.position = "none")

      gridExtra::grid.arrange(ptemp.ptm, ptemp.protein, ncol=1)

      message("Drew the Quality Contol plot(boxplot) for all ptms/proteins.")
    }

    ## each protein
    ## choose Proteins or not
    if (which.Protein != 'allonly') {
      if (which.Protein != "all") {
        ## check which.Protein is name of Protein
        if (is.character(which.Protein)) {
          temp.name <- which.Protein

          ## message if name of Protein is wrong.
          if (length(setdiff(temp.name,unique(datafeature.ptm$Protein))) > 0) {
            stop(paste0("Please check protein name. Data set does not have this protein. - ",
                        toString(temp.name)))
          }
        }

        ## check which.Protein is order number of Protein
        if (is.numeric(which.Protein)) {
          temp.name <- levels(datafeature.ptm$Protein)[which.Protein]

          ## message if name of Protein is wrong.
          if (length(levels(datafeature.ptm$Protein)) < max(which.Protein)) {
            stop(paste0("Please check your ion of proteins. There are ",
                        length(levels(datafeature.ptm$Protein))," proteins in this dataset."))
          }
        }

        ## use only assigned proteins
        datafeature.ptm <- datafeature.ptm[which(datafeature.ptm$Protein %in% temp.name), ]
        temp_proteins <- as.character((datafeature.ptm %>% distinct(GlobalProtein))[[1]])
        datafeature.ptm$Protein <- factor(datafeature.ptm$Protein)

        datafeature.protein <- datafeature.protein[which(datafeature.protein$Protein %in% temp_proteins), ]
        datafeature.protein$Protein <- factor(datafeature.protein$Protein)

        datarun.protein <- datarun.protein[which(datarun.protein$Protein %in% temp_proteins), ]
        datarun.ptm <- datarun.ptm[which(datarun.ptm$Protein %in% temp.name), ]
        datarun.protein$Protein <- factor(datarun.protein$Protein)
        datarun.ptm$Protein <- factor(datarun.ptm$Protein)
      }

      ## Only plot proteins that occur in both datasets
      global_proteins <- (datafeature.protein %>% distinct(Protein))[[1]]
      ptm_proteins <- (datafeature.ptm %>% distinct(GlobalProtein))[[1]]
      plot_proteins <- intersect(ptm_proteins, global_proteins)

      datafeature.ptm <- datafeature.ptm %>% filter(GlobalProtein %in% plot_proteins)

      ## TODO: If protein included in which.protein is not available in both datasets, print that here:
      plot_proteins <- (datafeature.ptm %>% distinct(Protein))[[1]]
      ## factoring for run, channel, condition should be done before loop

      for (i in 1:length(plot_proteins)) {

        sub.ptm <- datafeature.ptm[datafeature.ptm$Protein == as.character(plot_proteins[i]), ]
        sub.protein <- datafeature.protein[datafeature.protein$Protein == as.character((sub.ptm %>% distinct(GlobalProtein))[[1]]), ]

        sub.ptm <- sub.ptm[!is.na(sub.ptm$abundance), ]
        sub.protein <- sub.protein[!is.na(sub.protein$abundance), ]

        ## if all protein measurements are NA,
        if (nrow(sub.ptm) == sub.ptm[!is.na(sub.ptm$abundance), ]) {
          message(paste("Can't the Quality Control plot for ", unique(sub.ptm$Protein),
                        "(", i, " of ", length(plot_proteins),
                        ") because all measurements are NAs."))
          next()
        }
        ## if all ptm measurements are NA,
        if (nrow(sub.protein) == sub.protein[!is.na(sub.protein$abundance), ]) {
          message(paste("Can't the Quality Control plot for ", unique(sub.protein$Protein),
                        "(", i, " of ", length(plot_proteins),
                        ") because all measurements are NAs."))
          next()
        }

        ## for annotation of condition
        groupline.tmp.ptm <- data.frame(groupline.ptm,
                                    "PSM" = unique(sub.ptm$PSM)[1],
                                    "PeptideSequence" = unique(sub.ptm$PeptideSequence)[1])

        groupline.all.tmp.ptm <- data.frame(groupline.all.ptm,
                                        "PSM" = unique(sub.ptm$PSM)[1],
                                        "PeptideSequence" = unique(sub.ptm$PeptideSequence)[1])

        groupline.tmp.protein <- data.frame(groupline.protein,
                                    "PSM" = unique(sub.protein$PSM)[1],
                                    "PeptideSequence" = unique(sub.protein$PeptideSequence)[1])

        groupline.all.tmp.protein <- data.frame(groupline.all.protein,
                                        "PSM" = unique(sub.protein$PSM)[1],
                                        "PeptideSequence" = unique(sub.protein$PeptideSequence)[1])

        ## 1st plot for original plot
        ## for boxplot, x-axis, xorder should be factor
        sub.ptm$xorder <- factor(sub.ptm$xorder)
        sub.protein$xorder <- factor(sub.protein$xorder)

        ptemp.ptm <- ggplot(aes_string(x = 'xorder', y = 'abundance'), data = sub.ptm) +
          facet_grid(~Run) +
          geom_boxplot(aes_string(fill = 'Condition'), outlier.shape = 1, outlier.size = 1.5) +
          labs(title = paste0('PTM - ', unique(sub.ptm$Protein)),
               x = 'MS runs') +
          scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
          geom_vline(data = groupline.tmp.ptm,
                     aes(xintercept = cumGroupAxis),
                     colour = "grey", linetype = "longdash") +
          geom_text(data = groupline.all.tmp.ptm,
                    aes(x = xorder, y = abundance, label = Condition),
                    size = text.size,
                    angle = text.angle, hjust = .9,
                    color = "black") +
          theme(
            panel.background = element_rect(fill = 'white', colour = "black"),
            legend.key = element_rect(fill = 'white', colour = 'white'),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = 'gray95'),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = y.axis.size, colour = "black"),
            axis.ticks = element_line(colour = "black"),
            axis.title.x = element_text(size = x.axis.size + 5, vjust = -0.4),
            axis.title.y = element_text(size = y.axis.size + 5, vjust = 0.3),
            title = element_text(size = x.axis.size + 8, vjust = 1.5),
            legend.position = "none")

        ptemp.protein <- ggplot(aes_string(x = 'xorder', y = 'abundance'), data = sub.protein) +
          facet_grid(~Run) +
          geom_boxplot(aes_string(fill = 'Condition'), outlier.shape = 1, outlier.size = 1.5) +
          labs(title = paste0('Protein - ', unique(sub.protein$Protein)),
               x = 'MS runs') +
          scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
          geom_vline(data = groupline.tmp.protein,
                     aes(xintercept = cumGroupAxis),
                     colour = "grey", linetype = "longdash") +
          geom_text(data = groupline.all.tmp.protein,
                    aes(x = xorder, y = abundance, label = Condition),
                    size = text.size,
                    angle = text.angle, hjust = .9,
                    color = "black") +
          theme(
            panel.background = element_rect(fill = 'white', colour = "black"),
            legend.key = element_rect(fill = 'white', colour = 'white'),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = 'gray95'),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = y.axis.size, colour = "black"),
            axis.ticks = element_line(colour = "black"),
            axis.title.x = element_text(size = x.axis.size + 5, vjust = -0.4),
            axis.title.y = element_text(size = y.axis.size + 5, vjust = 0.3),
            title = element_text(size = x.axis.size + 8, vjust = 1.5),
            legend.position = "none")

        gridExtra::grid.arrange(ptemp.ptm, ptemp.protein, ncol=1)

        message(paste("Drew the Quality Contol plot(boxplot) for ", unique(sub.ptm$Protein),
                      "(", i, " of ", length(unique(datafeature.ptm$Protein)), ")"))

      } # end-loop
    }

    if (address != FALSE) {
      dev.off()
    }
  } # end QC plot

}
