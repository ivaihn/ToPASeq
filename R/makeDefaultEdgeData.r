makeDefaultEdgeData<-function(){
graphite2SPIA<-data.frame(
 type=c("binding", "control(In(ACTIVATION))", "control(In(INHIBITION))",
   "control(Out(ACTIVATION))", "control(Out(INHIBITION))", 
   "control(Out(INHIBITION-COMPETITIVE))", "control(Out(ACTIVATION_UNKMECH))",
   "control(Out(unknown))", "control(indirect)", "process", "process(BiochemicalReaction)",
   "process(activation)",  "process(binding/association)", "process(dephosphorylation)",
   "process(dissociation)", "process(expression)", "process(indirect effect)", 
   "process(indirect)", "process(inhibition)", "process(missing interaction)",
   "process(missing)", "process(phosphorylation)", "process(repression)",
   "process(ubiquitination)", "process(methylation)", "process(state change)")
 ,
 spiaType=c("binding/association", "activation", "inhibition", "activation", "inhibition", "inhibition", "activation", 
"indirect effect", "indirect effect", "activation", "activation", "activation", "binding/association", 
"dephosphorylation", "dissociation", "expression", "indirect effect", "indirect effect", "inhibition", 
"indirect effect", "indirect effect", "phosphorylation", "inhibition", "ubiquination", "inhibition", "ubiquination") 
 , stringsAsFactors=FALSE)
 graphite2SPIA<-as.matrix(graphite2SPIA)
 beta<-data.frame(
 rel=c("activation", "compound", "binding/association", "expression", "inhibition", "activation_phosphorylation", 
 "phosphorylation", "inhibition_phosphorylation", "inhibition_dephosphorylation", "dissociation", "dephosphorylation", 
 "activation_dephosphorylation", "state change", "activation_indirect effect", "inhibition_ubiquination", "ubiquination", 
 "expression_indirect effect", "inhibition_indirect effect", "repression", "dissociation_phosphorylation", "indirect effect_phosphorylation", 
 "activation_binding/association", "indirect effect", "activation_compound", "activation_ubiquination")
 ,
 beta=c(1, 0, 0, 1, -1, 1, 0, -1, -1, 0, 0, 1, 0, 1, -1, 0, 1, -1, -1, 0, 0, 1, 0, 1, 1)
 , stringsAsFactors=FALSE
 )
 defaultEdgeAttrs<-list(graphite2SPIA=graphite2SPIA, beta=beta)
 return(defaultEdgeAttrs)
}


 
