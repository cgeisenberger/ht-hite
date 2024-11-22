
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

# Define a function to extract insert sizes from a BAM file
extract_insert_sizes <- function(bam_file) {
  # Open the BAM file
  bam <- BamFile(bam_file, yieldSize = 100000)  # Adjust yieldSize for memory efficiency
  open(bam)
  
  # Initialize an empty vector to store insert sizes
  insert_sizes <- c()
  
  # Iterate through the BAM file in chunks
  while (length(chunk <- scanBam(bam, param = ScanBamParam(what = c("isize")))[[1]]$isize) > 0) {
    insert_sizes <- c(insert_sizes, chunk)  # Append insert sizes
  }
  
  # Close the BAM file
  close(bam)
  
  # Return the insert sizes as a numeric vector
  return(as.numeric(insert_sizes))
}

read_bed <- function(path){
  colnames_bed <- c("chr", "start", "end", "feature", "coverage", "non_zero_bases", "length", "non_zero_prop")
  id = str_extract(path, "E[0-9]{2}_[0-9]*_(Beads|Control|Columns|Promega)")
  data <- read_table(file = path, col_names = colnames_bed)
  data <- data %>% 
    drop_na() %>% 
    add_column(id, .before = 1)
  return(data)
}

cv <- function(x){
  y = mean(x)/sd(x)
}

read_vcf_file <- function(path){
  raw <- read.vcfR(path)
  sample_name <- str_extract(path, pattern = "E[0-9]{2}_[0-9]+")
  
  variants <- tibble(
    sample = rep(sample_name, nrow(raw@fix)),
    chr = raw@fix[, "CHROM"],
    pos = as.numeric(raw@fix[, "POS"]),
    ref = raw@fix[, "REF"],
    alt = raw@fix[, "ALT"], 
    filter = raw@fix[, "FILTER"],
    depth = extract.info(raw, element = "DP")
  )
  
  variants <- variants %>% 
    mutate(depth = as.numeric(depth))
  
  return(variants)
}

read_flagstat <- function(path){
  
  raw <- read_tsv(file = path, col_names = c("stat", "accessory_stat", "variable"))
  
  sample = str_extract(path, pattern = "E[0-9]{2}_[0-9]*")
  
  method = str_extract(path, pattern = "(Beads|Columns|Promega|Control)")
  method = ifelse(method == "Promega", "Reference", method)
  
  variable <- str_extract(string = raw$variable, pattern = "^[^(]*") %>% 
    trimws %>% 
    str_replace_all(pattern = " ", replacement = "_") %>% 
    str_replace_all(pattern = "%", replacement = "percent")
  print(variable)
  
  reformated <- tibble(sample = sample,
                       method = method, 
                       measure = variable, 
                       data = as.numeric(str_remove(string = raw$stat, pattern = "%")))
  return(reformated[1:18, ])
}


calculate_variants_overlap <- function(vcf_list, min_depth = NULL, sample_names = NULL) {
  
  # iniate matrix
  variant_accordance <- matrix(data = NaN, nrow = length(vcf_list), ncol = length(vcf_list))
  
  # set depth for filtering
  min_depth <- ifelse(is.null(min_depth), 0, min_depth)
  
  for (i in 1:length(vcf_list)) {
    print(i)
    
    vars_i <- vcf_list[[i]] %>% 
      filter(depth > min_depth) %>% 
      pull(var_id)
    
    for (j in 1:length(vcf_list)){
      
      vars_j <- vcf_list[[j]] %>% 
        filter(depth > min_depth) %>% 
        pull(var_id)
      
      variant_accordance[i, j] <- length(intersect(vars_i, vars_j))/
        length(union(vars_i, vars_j))
    }
  }
  
  # add sample names if supplied
  if(!is.null(sample_names)){
    colnames(variant_accordance) <- sample_names
    rownames(variant_accordance) <- sample_names
  }
  
  return(variant_accordance)
}
