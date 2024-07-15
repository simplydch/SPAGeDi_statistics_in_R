### Additional Functions

create_temp_file <- function(file_ext=".txt"){
    tmp_file <- paste0("tmp_",paste(sample(c(LETTERS, letters, rep(1:9, each=3)), 10), 
                           collapse=""), file_ext)
    file.create(tmp_file)
    cat(paste0("Created file: ", tmp_file ))
    tmp_file
}


