
edgetypeslegend<-function(x){
par(mai=c(1.02,0,0,0))
plot(1, ylim=c(2,7), xlim=c(1,5), type="n", xlab="", ylab="", frame=FALSE, axes=FALSE)


text( 0.8, 6.5, "Edge types", pos=4, font=2)
lines(c(1,2), c(4,4), lty="dashed")
text( 2.2, 4, "dissociation", pos=4)
lines(c(1,2), c(4.5,4.5), lty="dashed")
lines(c(1.8,2,1.8), c(4.45,4.5,4.55), lty="solid")
text( 2.2, 4.5, "indirect effect", pos=4)

lines(c(1,2), c(5,5), lty="solid")
text( 2.2, 5, "bindig, association", pos=4)

lines(c(1,2), c(5.5,5.5), lty="solid")
lines(c(2,2), c(5.4,5.6), lty="solid")
text( 2.2, 5.5, "inhibition", pos=4)

lines(c(1,2), c(6,6), lty="solid")
lines(c(1.8,2,1.8), c(5.95,6,6.05), lty="solid")
text( 2.2, 6, "activation", pos=4)

text(0.8, 3.5, "Edge labels", pos=4, font=2)
text( 1.5, 3, "+p")
text( 2.2, 3, "phosphorylation", pos=4)
text( 1.5, 2.5, "-p")
text( 2.2, 2.5, "dephosphorylation", pos=4)

text( 1.5, 2, "+u")
text( 2.2, 2, "ubiquination", pos=4)
}


