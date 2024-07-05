
geom_mean <- function(x) {
  y  <-  exp(mean(log(x)))
  return(y)
}

split_data <- function(p_train = NULL){
  if (is.null(p_train)) {
    p_train = 0.5
  }
  p_test = 1 - p_train
  
  model_data <- plates %>% 
    filter(is.na(evaporation_beads)) %>% 
    filter(control == "no") %>% 
    mutate(dataset = sample(c("train", "test"), prob = c(p_train, p_test),
                            size = nrow(.), replace = TRUE))
  
  res <- list(train = filter(model_data, dataset == "train"), 
              test = filter(model_data, dataset == "test"))
  return(res)
  
}

calculate_correction <- function(nd, qb, n){
  ind <- order(nd)[ceiling(seq(1, to = length(nd), length.out = n))]
  cf <- geom_mean(nd[ind]/qb[ind])
  return(cf)
}

remove_epic_suffix <- function(methylSet){
  rownames(methylSet) <- str_extract(string = rownames(methylSet), pattern = "cg[0-9]*")
  return(methylSet)
}