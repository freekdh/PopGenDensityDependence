# Data Analysis plot

setwd("/home/freek/PopGenDensityDependence/28-03-2018-12-13-08")
library(grid)
library(gridGraphics)
library(gtable)
library(gridExtra)
library(lattice)

list.data<-list()

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

for(i in 1:NMETA){
    B[,i] <- list.data[[i]]$Allele0 / (list.data[[i]]$Allele0 + list.data[[i]]$Allele1)
}

########### PLOT WITH POPULATION SIZE ############

for(i in 1:NMETA){
    list.data[[i]]$Total <- list.data[[i]]$Allele0 + list.data[[i]]$Allele1
    list.data[[i]]$Total[list.data[[i]]$Total == 0] <- NA
}

ColorRamp <- rgb( seq(1,1,length=256),  # Red
                    seq(1,0,length=256),  # Green
                    seq(0,0,length=256))  # Blue

ColorLevels <- seq(0, 1, length=length(ColorRamp))


par(mfrow=c(1,1), cex=1.5 , mar = c(0,0,0,0))
image(1:nrow(B),1:ncol(B),B,col = ColorRamp)
legend("topleft", legend = c(paste("r0 = ", pars$value[1]),
                            paste("r1 = ",pars$value[2]),
                            paste("k0 = ", pars$value[3]),
                            paste("k1 = ",pars$value[4]),
                            paste("n01 = ",pars$value[5]),
                            paste("n02 = ",pars$value[6]),
                            paste("rep = ", pars$value[10])), 
                  bty = "n")
par(mfrow=c(NMETA,1), mar = c(0,0,0,0))
for(i in 1:NMETA-1){
    par(mfg=c(i+1,1))
    plot(list.data[[NMETA-i]]$Total, type = 'l',xaxt='n', yaxt='n', xaxs = "i",yaxs="i" , lwd = 3,
    ylim = c(0,1.05*max(pars$value[3], pars$value[4])*pars$value[10]),
    bty="n")
}

grid.echo()
a <- grid.grab()
grid.newpage()
pushViewport(viewport(y=0.08, height = .9, width = .825, just = "bottom"))
grid.draw(a)
a <- grid.grab()
dev.off()

par(cex=1.5)
image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n")

grid.echo()
b <- grid.grab()
dev.off()
pdf("GraphPDF.pdf",width = 15, height = 10)
grid.arrange(a, b, ncol=2, widths = c(10,1),
               top=textGrob("Title", gp=gpar(fontsize=20,font=8)))
grid.rect(gp=gpar(fill=NA))
dev.off()


# figure out the legend: