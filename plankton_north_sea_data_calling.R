d_allsp_data_plankton_north_sea<-function(d0){
  
  name.plankton<-c('Calanus_I-IV','Para-Pseudocalanus_spp','Temora_longicornis','Acartia_spp_(unidentified)','Centropages_typicus',
                   'Oithona_spp','Echinoderm_larvae','Calanus_finmarchicus','Calanus_helgolandicus','Metridia_lucens','Decapoda_larvae_(Total)',
                   'Euphausiacea_Total','Thalassiosira_spp','Rhizosolenia_styliformis','Ceratium_fusus','Ceratium_furca','Ceratium_tripos',
                   'Ceratium_macroceros','Rhizosolenia_alata_alata','Pseudocalanus_elongatus_Adult','Nitzschia_delicatissima','Nitzschia_seriata')
  
  longs<-read.csv("Data/Plankton_North_Sea_data/boxcornerlongs201117.csv",header=F)
  lats<-read.csv("Data/Plankton_North_Sea_data/boxcornerlats201117.csv",header=F)
  longs<-longs[[1]]
  lats<-lats[[1]]
  d0[d0==-999]<-NA
  Years<-1958:2013
  
  d_allsp<-list() #this is a list with 22 elements, all the species
  for (sp_c in 1:22)
  {
    d_thissp<-list() #a list with 26 locations for the current species
    d1<-d0[((sp_c-1)*26+1):(sp_c*26),]
    for (st_c in 1:26)
    {
      d_thissp[[st_c]]<-data.frame(Year=Years,Dat=d1[st_c,])
    }
    d_allsp[[sp_c]]<-d_thissp
  }
  names(d_allsp)<-name.plankton
  
  
  return(d_allsp=d_allsp)
  
}