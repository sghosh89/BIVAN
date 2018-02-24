
d_allsp_data<-function(d0){
  
  name.aphid<-c("Apple-grass aphid","Bird cherry-oat aphid","Black bean aphid", "Blackberry-cereal aphid","Blackcurrant-sowthistle aphid",
                "Corn leaf aphid","Currant-lettuce aphid","Damson-hop aphid","Grain aphid","Green spruce aphid",
                "Leaf-curling plum aphid","Mealy cabbage aphid","Mealy plum aphid","Pea aphid","Peach-potato aphid",
                "Potato aphid","Rose-grain aphid","Shallot aphid","Sycamore aphid","Willow-carrot aphid")
  longs<-c(-4.567,0.57,-3.069,-3.312,-2.637,-1.682,-2.763,-0.356,-3.454,0.427,0.939)
  lats<-c(55.477,52.26,56.457,55.949,52.125,55.213,53.854,51.807,50.628,51.733,51.185)
  loc_name<-c("Ayr","Broom's Barn","Dundee","Edinburgh","Hereford","Newcastle","Preston","Rothamsted",
              "Starcross","Writtle","Wye")
  d0[d0==-999]<-NA
  Years<-1976:2010
  d_allsp<-list() #this is a list with 20 elements, all the species
  for (sp_c in 1:20)
  {
    d_thissp<-list() #a list with 11 locations for the current species
    for (st_c in 1:11)
    {
      d_thissp[[st_c]]<-data.frame(Year=Years,Dat=d0[st_c,((sp_c-1)*35+1):(sp_c*35)])
    }
    d_allsp[[sp_c]]<-d_thissp
  }
  names(d_allsp)<-name.aphid
  
  return(d_allsp=d_allsp)
  
}
