
# ----- Define a function for plotting a matrix ----- #
myImagePlot <- function(x, ...){
     min <- 0
     max <- 1
     yLabels <- rownames(x)
     xLabels <- colnames(x)
     title <-c()
    # check for additional function arguments
    if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
       min <- Lst$zlim[1]
       max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
       yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
       xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
       title <- Lst$title
    }
    }
    # check for null values
    if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
    }
    if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
    }

    layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))

    # Red and green range from 0 to 1 while Blue ranges from 1 to 0
    ColorRamp <- rgb( seq(1,1,length=256),  # Red
                    seq(1,0,length=256),  # Green
                    seq(0,0,length=256))  # Blue
    ColorLevels <- seq(min, max, length=length(ColorRamp))

    # Data Map
    image(1:length(xLabels), 1:length(yLabels), t(x), col = ColorRamp, xlab="",
    ylab="", axes=FALSE, zlim=c(min,max))
    
    axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
    axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
    cex.axis=0.7)

    # Color Scale
    par(mar = c(3,2.5,2.5,2))
    image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")

    layout(1)
}
# ----- END plot function ----- #
# create an empty list that will serve as a container to receive the incoming files

list.data<-list()

setwd(paste("/home/freek/PopGenDensityDependence/", "27-03-2018-18-16-33" , sep = "/"))
    
pars <- read.csv("Parameters.csv", nrows = 11, header = FALSE)
colnames(pars) <- c("parameter", "value")

NMETA = pars$value[9]

for(i in 0: (NMETA-1)){
    datasub <- paste(getwd(), "SubPopulations", sep = "/")
    datasub <- paste(datasub, "Sub", sep = "/")
    datasub <- paste(datasub, i, sep = "")

    datasub <- paste(datasub, "Locus0.csv", sep = "/")
    list.data[[i+1]]<-read.csv(datasub, header = FALSE)
    colnames(list.data[[i+1]]) <- c("Allele0", "Allele1")
}

B = matrix(
    nrow = nrow(list.data[[1]]),
    ncol = NMETA
)

########### PLOT WITH POPULATION SIZE ############

for(i in 1:NMETA){
    B[,i] <- list.data[[i]]$Allele0 / (list.data[[i]]$Allele0 + list.data[[i]]$Allele1)
}

dev.off()
par(mfrow=c(1,1), mar = c(0,0,0,0))
image(B)
par(mfrow=c(NMETA,1), mar = c(0,0,0,0), new = TRUE)
for(i in 1:NMETA){
    par(mfg=c(i+1,1))
    plot(list.data[[NMETA-i]]$Total, type = 'l',xaxt='n', yaxt='n', xaxs = "i",yaxs="i" , lwd = 2,
    ylim = c(0,1.15*max(pars$value[3], pars$value[4])*pars$value[10]),
    bty="n")
}
