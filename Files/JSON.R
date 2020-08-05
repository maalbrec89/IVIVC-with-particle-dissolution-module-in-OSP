rm(list = setdiff(ls(), lsf.str()))

require(jsonlite)
require(xlsx)

# load PSV file with particle radius and rel_amount:
PSVfilename <- "PSV_CompoundA_BatchX_lognormal.xlsx"
PSV <- read.xlsx(paste0(file.path(getwd(),PSVfilename)),1)

##################################################
########## PARSE FORMULATION BB ##################
##################################################

# load json formulation file:
filename <- "BB_Formulation_ParticleDissolution_10Bins_input.json"
BB <- fromJSON(file.path(getwd(),filename))

PSV[,3] <- PSV[,3]/1000   # convert from µm to mm

for(i in 1:10){
  radius_i <- PSV[i,3]
  BB$Parameters[[i]]$Value[3] <- radius_i
}

write(toJSON(BB, digits=I(5), pretty=T, auto_unbox=T), file.path(getwd(),"/BB_Formulation_ParticleDissolution_10Bins_output.json"))

##################################################
########## PARSE ADMINISTRATION BB ###############
##################################################

# load json administration file:
filename <- "BB_Administration_ParticleDissolution_10Bins_input.json"
BB <- fromJSON(file.path(getwd(),filename))

Dose <- 0.1      # Dose in g
MW <- 200        # Molecular weight in g/mol
Sum_rel_amountFactor <- sum(PSV[11:20,3])
amountCorrectionFactorForBin = Dose/MW/Sum_rel_amountFactor

for(i in 1:10){
  rel_amountFactor_i <- PSV[10+i,3]
  
  startAmount_i <- rel_amountFactor_i*amountCorrectionFactorForBin
  startMass_i <- startAmount_i*MW
  
  BB$Schemas$SchemaItems[[1]]$Parameters[[i]][2,2] <- startMass_i
}

write(toJSON(BB, digits=I(5), pretty=T, auto_unbox=T), paste0(getwd(),"/BB_Administration_ParticleDissolution_10Bins_output.json"))
