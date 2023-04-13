addGG <- function(d) {
  
  dOut <- 
    d |> 
    mutate(
      riverGG = factor(
        river,
        levels=c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),
        labels = c("West Brook","Open Large","Open Small","Isolated Small"),
        ordered = T
      ),
      seasonGG = factor(
        season, 
        labels = c("Spring","Summer","Autumn","Winter"), 
        ordered = T
      ),
      speciesGG = factor(
        species, 
        levels = c('bkt','bnt','ats'), 
        labels = c("Brook trout", "Brown trout", "Atlantic salmon"), 
        ordered = T
      )
    )
  return(dOut)
}

addGG_noSpecies <- function(d) {
  
  dOut <- 
    d |> 
    mutate(
      riverGG = factor(
        river,
        levels=c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),
        labels = c("West Brook","Open Large","Open Small","Isolated Small"),
        ordered = T
      ),
      seasonGG = factor(
        season, 
        labels = c("Spring","Summer","Autumn","Winter"), 
        ordered = T
      )
    )
  return(dOut)
}


