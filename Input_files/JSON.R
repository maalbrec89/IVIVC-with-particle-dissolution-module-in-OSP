rm(list = setdiff(ls(), lsf.str()))

require(jsonlite)
require(xlsx)

# set working directory to "<your local path>/IVIVE-with-particle-dissolution-module-in-OSP"
workingDir <- file.path("C:/OSP","IVIVE-with-particle-dissolution-module-in-OSP") 
setwd(workingDir)

# Define dose and molecular weight
Dose <- 50       # Dose in mg
MW <- 200        # Molecular weight in g/mol

# load PSV file with particle radius and rel_amount:
PSVfilename <- "PSV_CompoundA_BatchX_lognormal.xlsx"
PSV <- read.xlsx(paste0(file.path("Output_files",PSVfilename)),1)

##################################################
########## PARSE FORMULATION BB ##################
##################################################

# load json formulation file:
filename <- "BB_Formulation_ParticleDissolution_10Bins_input.json"
BB <- fromJSON(file.path("Input_files/",filename))
nBin <- length(unique(PSV$Container.Path))

PSV[PSV$Parameter.Name=="radius (at t=0)","Value"] <- PSV[PSV$Parameter.Name=="radius (at t=0)","Value"]/1000   # convert radii from µm to mm
for(i in 1:nBin){
  radius_i <- PSV[PSV$Parameter.Name=="radius (at t=0)",3][i]
  BB$Parameters[[i]]$Value[which(BB$Parameters[[i]]$Name=="Particle radius (mean)")] <- radius_i
}

write(toJSON(BB, digits=I(5), pretty=T, auto_unbox=T), file.path("Output_files","BB_Formulation_ParticleDissolution_10Bins_output.json"))

##################################################
########## PARSE ADMINISTRATION BB ###############
##################################################

# load json administration file:
filename <- "BB_Administration_ParticleDissolution_10Bins_input.json"
BB <- fromJSON(file.path("Input_files",filename))

Sum_rel_amountFactor <- sum(PSV[PSV$Parameter.Name=="rel_amountFactor","Value"])
amountCorrectionFactorForBin = Dose/MW/Sum_rel_amountFactor

for(i in 1:nBin){
  rel_amountFactor_i <- PSV[PSV$Parameter.Name=="rel_amountFactor","Value"][i]
  
  startAmount_i <- rel_amountFactor_i*amountCorrectionFactorForBin
  startMass_i <- startAmount_i*MW
  
  BB$Schemas$SchemaItems[[1]]$Parameters[[i]][BB$Schemas$SchemaItems[[1]]$Parameters[[i]]$Name=="InputDose","Value"] <- startMass_i
}

write(toJSON(BB, digits=I(5), pretty=T, auto_unbox=T), file.path("Output_files","BB_Administration_ParticleDissolution_10Bins_output.json"))
