## Run the multispecies SMS model ## 
library(smsR)


wd <- "C:/Users/nsja/Documents/Github/wg_WGSAM/NorthSeaKeyRun_2020"

# load multispecies data # 

years <- 1974:2019
nyears <- length(years)
seasons <- 1:4
nspecies <- 27
maxage <- 10

spp <- c("Fulmar",
         "Guillemot",
         "Her.Gull",
         "Kittiwake",
         "GBB.Gull",
         "Gannet",
         "Puffin",
         "Razorbill",
         "R.radiata",
         "G.gurnards",
         "W.horsemac",
         "N.horsemac",
         "Greyseal",
         "H.porpoise",
         "Hake",
         "Cod",
         "Whiting",
         "Haddock",
         "Saithe",
         "Mackerel",
         "Herring",
         "N.sandeel",
         "S.sandeel",
         "Nor.pout",
         "Sprat",
         "Plaice",
         "Sole")


is.dynamic <- rep(0, length(spp))
is.dynamic[16:27] <- 1



