args <- commandArgs(trailingOnly = TRUE)
Foldername <- args[1]

setwd(Foldername)

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

pdf("PlotData.pdf")
par(mfrow=c(1,1), mar = c(0,0,0,0))
image(B,col = ColorRamp)
text(0.1, 0.95, paste("r0 = ", pars$value[1]))
text(0.1, 0.9, paste("r1 = ", pars$value[2]))
text(0.1, 0.85, paste("k0 = ", pars$value[3]))
text(0.1, 0.8, paste("k1 = ", pars$value[4]))
par(mfrow=c(NMETA,1), mar = c(0,0,0,0), new = TRUE)
for(i in 1:NMETA-1){
    par(mfg=c(i+1,1))
    plot(list.data[[NMETA-i]]$Total, type = 'l',xaxt='n', yaxt='n', xaxs = "i",yaxs="i" , lwd = 2,
    ylim = c(0,1.15*max(pars$value[3], pars$value[4])*pars$value[10]),
    bty="n")
}
dev.off()

