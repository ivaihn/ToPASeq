makeDefaultEdgeData<-function(){
graphite2SPIA<-data.frame(
type=c("Binding", "Control(In: ACTIVATION of BiochemicalReaction)", "Control(In: ACTIVATION of ComplexAssembly)",
  "Control(In: ACTIVATION of Degradation)", "Control(In: ACTIVATION of Transport)", "Control(indirect)",           
  "Control(In: INHIBITION-COMPETITIVE of BiochemicalReaction)", "Control(In: INHIBITION of BiochemicalReaction)",
  "Control(In: INHIBITION of ComplexAssembly)", "Control(In: INHIBITION of Transport)", "Control(Out: ACTIVATION of ACTIVATION)",
  "Control(Out: ACTIVATION of BiochemicalReaction)", "Control(Out: ACTIVATION of ComplexAssembly)", "Control(Out: ACTIVATION of TemplateReaction)", "Control(Out: ACTIVATION of Transport)", "Control(Out: INHIBITION-COMPETITIVE of BiochemicalReaction)", "Control(Out: INHIBITION of ACTIVATION)", "Control(Out: INHIBITION of BiochemicalReaction)", "Control(Out: INHIBITION of ComplexAssembly)", 
  "Control(Out: INHIBITION of TemplateReaction)", "Control(Out: INHIBITION of Transport)", "Process", "Process(activation)", "Process(binding/association)", "Process(BiochemicalReaction)", "Process(dephosphorylation)", "Process(dissociation)", "Process(expression)", 
  "Process(glycosylation)", "Process(indirect)", "Process(indirect effect)", "Process(inhibition)", "Process(methylation)", "Process(missing)", "Process(missing interaction)", "Process(phosphorylation)", "Process(repression)", "Process(state change)", "Process(ubiquitination)"),
spiaType=c("binding/association", "activation", "activation", "activation", "activation", "indirect effect", "inhibition", "inhibition", "inhibition", "inhibition", "activation", "activation", "activation", "activation", "activation", "inhibition", "inhibition",
"inhibition", "inhibition", "inhibition", "inhibition", "activation", "activation", "binding/association", "activation",
"dephosphorylation", "dissociation", "expression", "glycosylation", "indirect effect", "indirect effect", "inhibition", "methylation",
"indirect effect", "indirect effect", "phosphorylation", "inhibition", "indirect effect",  "ubiquitination")
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


 
