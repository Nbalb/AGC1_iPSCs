pretty_path_label <- function(pathways, wrap = NULL){
  
  split_path <- str_replace_all(pathways, "_", " ") %>% 
    str_split(" ")
  
  split_NA <- function(x){
    
    term <- str_c(x, collapse = " ") %>% 
      str_extract("[A-Z]*[D,R]+NA")
    term1 <- str_remove(term, "[D,R,T]+NA") %>% str_to_lower()
    term2 <- str_extract(term, "[D,R,T]+NA")
    str_c(term1, term2)
  }
  
  map_chr(split_path, function(x){
    
  pret <- case_when(
      str_detect(x, "[A-Z]*[D,R]+NA") ~ split_NA(x),
      str_detect(x, "^GO[A-Z]{2}|WP") ~ str_to_upper(x),
      TRUE ~ str_to_title(x)
    ) %>% 
      str_c(collapse = " ")
    
  if(!is.null(wrap)){
    
   pret <- str_wrap(pret, width = wrap) 
   
  }
  
  return(pret)
  
  })
  
}

