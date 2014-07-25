library(photobiology)
library(ggtern)

oldwd <- setwd("raw.data/bees")
Maxwell.data <- read.table(file="Maxwell.data", header=TRUE)

leaf.conv.spct <- with(Maxwell.data, data.table(w.length = w.length,
                                                U = UV * D65 * Leaf,
                                                B = Blue * D65 * Leaf,
                                                G = Green * D65 * Leaf))

flower.conv.spct <- with(Maxwell.data, data.table(w.length = w.length,
                                                U = UV * D65 * Flower,
                                                B = Blue * D65 * Flower,
                                                G = Green * D65 * Flower))

R.flower <- 1 / integrate_spct(flower.conv.spct)
R.leaf <- 1 / integrate_spct(leaf.conv.spct)

S.leaf <- data.table(U = Maxwell.data$UV * R.leaf[1],
                     B = Maxwell.data$Blue * R.leaf[2],
                     G = Maxwell.data$Green * R.leaf[3])

P.flower <- R.leaf / R.flower

E.flower <- P.flower / (P.flower + 1)

rel.q.abs.flower <- P.flower / sum(P.flower)

S.flower <- data.table(U=rel.q.abs.flower[1], B=rel.q.abs.flower[2], G=rel.q.abs.flower[3])
S.flower <- rbind(data.table(U=0.5, B=0.5, G=0.5), S.flower)

S.leaf[ , spec.loc.U := U / (U + B + G)]
S.leaf[ , spec.loc.B := B / (U + B + G)]
S.leaf[ , spec.loc.G := G / (U + B + G)]


fig.spec0 <- ggtern(data=S.leaf, aes(x=spec.loc.U, y=spec.loc.B, z=spec.loc.G)) + geom_point() + geom_line()
fig.spec0 + geom_point(data=S.flower, aes(x=U, y=B, z=G)) + labs(x="UV", y="Blue", z="Green")


# Paula's data ------------------------------------------------------------

aitovirna <- read.table(file="aitovirna1sivu.spc", col.names=c("w.length", "Rpc"))
aitovirna.spct <- setReflectorSpct(aitovirna)
class(aitovirna.spct)
aitovirna.spct * D65.spct
interpolate_spct(D65.spct, aitovirna$w.length)$data.col * aitovirna[,Rfr]

setwd(oldwd)

cedta()
