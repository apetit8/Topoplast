f <- function(x, k=0.15) 1/(1+(1/k-1)*exp(-x/k/(1-k)))

pdf("../../figures/Model_SupfigA.pdf", width=7/2.54, height=5/2.54, pointsize=8)
    par(mar=c(5,4,1,1))

    curve(f(x), xlim=c(-1,1), xlab="Sum of regulations towards gene i", ylab="Expression of gene i", yaxt="n")
    axis(2, at=c(0,0.15,0.5,1), las=2)
    arrows(x0=0, y0=0, y1=0.15, lty=2, col="darkgray", length=0)
    arrows(x0=-1, x1=0, y0=0.15, lty=2, col="darkgray", length=0)
    text(x=-0.5, y=0.15, expression("Constitutive\nexpression "*kappa), pos=3)
    abline(a=0.15, b=1, lty=3, col="blue")
    arrows(x0=0.3-0.15, x1=0.5, y0=0.3, lty=3, col="blue", length=0)
    text(x=0.5, y=0.30, "slope = 1", col="blue", pos=3)

dev.off()

source("../functions/functions.R")

pp <- read.param("../../templates/Full_netw/param1.txt")
mm <- matrix(pp$TYPE_ALLELES == "normal", ncol=pp$GENET_NBLOC)

ll <- c(1, 10, 10, 10, 5)

pdf("../../figures/Model_SupfigB.pdf", width=7/2.54, height=5/2.54, pointsize=8)
    par(mar=c(1, 5, 1, 1))
    image(x=1:pp$GENET_NBLOC, y=1:pp$GENET_NBLOC, z=mm[,pp$GENET_NBLOC:1], xlab="", ylab="", xaxt="n", yaxt="n", asp=1, bty="n")
    abline(v=0.5+c(0, cumsum(ll)), col="darkgray", lty=2)
    abline(h=sum(ll)+0.5-c(0, cumsum(ll)), col="darkgray", lty=2)
    text(x=0, y=sum(ll)-cumsum(ll)+(ll+1)/2, c("sensor", "TF", "target\n(fluct)", "target\n(const)", "neutral"), pos=2, xpd=NA)
    
    text(16, 33, "mutable", col="white")
    text(16, 23, "non-mutable", col="darkred")
dev.off()
