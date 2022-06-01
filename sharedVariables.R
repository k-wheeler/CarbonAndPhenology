##Shared Variables
options(stringsAsFactors=FALSE)
nchain=5
lat <- 42.5351
long <- -72.1744
leafNames <- c("B1A","B1B","B1C","B1D","B1E","B1F",
               "B2A","B2B","B2C","B2D","B2E","B2F",
               "B3A","B3B","B3C","B3D","B3E","B3F",
               "O1A","O1B","O1C","O1D","O1E","O1F",
               "O2A","O2B","O2C","O2D","O2E","O2F",
               "O3A","O3B","O3C","O3D","O3E","O3F",
               "O4A","O4B","O4C","O4D","O4E","O4F")
trees <- c("B1","B2","B3","O1","O2","O3","O4")
letterSequence <- c("a","b","c","d","e","f")
CCIdataDirectory <- "Data/CCI_Measurements/"
LicorDataDirectory <- "Data/Licor_Measurements/"

licorLeafNames <- c("B1A","B1B","B1E","B2A","B2B","B3A","B3B","B3E",
                    "O1A","O1B","O1E","O2A","O2B","O2E","O3A","O3B","O4A","O4B")
finalDataDirectory <- "Data/finalData/"
oakLicorAugLeafNames <- c("O1A","O1B","O2A","O3A","O3B","O4A","O4B","O1E","O2E")
beechLicorAugLeafNames <- c("B1A","B1B","B1E","B2A","B3A","B2B","B2E","B3B","B3E")

combinations <- as.data.frame(rbind(
  c('Tair_D',TRUE,FALSE,FALSE,FALSE),
  c('Tair_D',FALSE,TRUE,FALSE,FALSE),
  c('Tair_D',FALSE,FALSE,FALSE,FALSE),
  c('TairOnly',TRUE,FALSE,FALSE,FALSE),
  c('TairOnly',FALSE,TRUE,FALSE,FALSE),
  c('TairOnly',FALSE,FALSE,FALSE,FALSE),
  c('NPP',TRUE,FALSE,FALSE,FALSE),
  c('NPP',FALSE,TRUE,FALSE,FALSE),
  c('NPP',FALSE,FALSE,FALSE,FALSE),
  c('NPP',TRUE,FALSE,TRUE,FALSE),
  c('NPP',FALSE,TRUE,TRUE,FALSE),
  c('NPP',FALSE,FALSE,TRUE,FALSE),
  c('GPP',TRUE,FALSE,FALSE,FALSE),
  c('GPP',FALSE,TRUE,FALSE,FALSE),
  c('GPP',FALSE,FALSE,FALSE,FALSE),
  c('GPP',TRUE,FALSE,TRUE,FALSE),
  c('GPP',FALSE,TRUE,TRUE,FALSE),
  c('GPP',FALSE,FALSE,TRUE,FALSE),
  c('NPP',TRUE,FALSE,TRUE,TRUE),
  c('NPP',FALSE,TRUE,TRUE,TRUE),
  c('NPP',FALSE,FALSE,TRUE,TRUE),
  c('GPP',TRUE,FALSE,TRUE,TRUE),
  c('GPP',FALSE,TRUE,TRUE,TRUE),
  c('GPP',FALSE,FALSE,TRUE,TRUE),
  c('NPP_aging',TRUE,FALSE,FALSE,FALSE),
  c('NPP_aging',FALSE,TRUE,FALSE,FALSE),
  c('NPP_aging',FALSE,FALSE,FALSE,FALSE),
  c('GPP_aging',TRUE,FALSE,FALSE,FALSE),
  c('GPP_aging',FALSE,TRUE,FALSE,FALSE),
  c('GPP_aging',FALSE,FALSE,FALSE,FALSE),
  c('NPP_feedback',TRUE,FALSE,FALSE,FALSE),
  c('NPP_feedback',FALSE,TRUE,FALSE,FALSE),
  c('NPP_feedback',FALSE,FALSE,FALSE,FALSE),
  c('GPP_feedback',TRUE,FALSE,FALSE,FALSE),
  c('GPP_feedback',FALSE,TRUE,FALSE,FALSE),
  c('GPP_feedback',FALSE,FALSE,FALSE,FALSE),
  c('NPP_feedback',TRUE,FALSE,TRUE,FALSE),
  c('NPP_feedback',FALSE,FALSE,TRUE,FALSE),
  c('GPP_feedback',TRUE,FALSE,TRUE,FALSE),
  c('GPP_feedback',FALSE,FALSE,TRUE,FALSE)
))
colnames(combinations) <- c("Cov","missingYear","excludePostSOS","Addition","interaction")


harvardSOStransDOYs <- 279.03+c(-1.21,-0.71,-3.61,0.09,0.89,1.39,0.71,1.95,0.81,-0.92,0.03)
