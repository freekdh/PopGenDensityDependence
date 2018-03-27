
# ----- Define a function for plotting a matrix ----- #
myImagePlot <- function(x, ...){
     min <- min(x)
     max <- max(x)
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
    ColorRamp <- rgb( seq(0,1,length=256),  # Red
                    seq(0,1,length=256),  # Green
                    seq(1,0,length=256))  # Blue
    ColorLevels <- seq(min, max, length=length(ColorRamp))

    # Reverse Y axis
    reverse <- nrow(x) : 1
    yLabels <- yLabels[reverse]
    x <- x[reverse,]

    # Data Map
    par(mar = c(3,5,2.5,2))
    image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
    ylab="", axes=FALSE, zlim=c(min,max))
    if( !is.null(title) ){
        title(main=title)
    }
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
for(i in 0: 49){
    file <- paste(getwd(), "R_test_data" , sep = "/")
    file <- paste(file, "SubPopulations", sep = "/")
    file <- paste(file, "Sub", sep = "/")
    file <- paste(file, i, sep = "")
    file <- paste(file, "Locus0.csv", sep = "/")

    list.data[[i+1]]<-read.csv(file)
    colnames(list.data[[i+1]]) <- c("Locus1", "Locus2")
    }

B = matrix(
    nrow = 199,
    ncol = 50
)

for(i in 1:50){
    B[,i] <- list.data[[i]]$Locus1 / (list.data[[i]]$Locus1 + list.data[[i]]$Locus2)
}
B[is.na(B)] <- 0
myImagePlot(t(B))
