setwd("/protestify_data/images/addlexpts")

jpg <- list.files(pattern="*.jpg")
png <- list.files(pattern="*.png")

jpg_num <- length(jpg)
png_num <- length(png)

my.file.rename <- function(from, to) {
  todir <- dirname(to)
  if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
  file.rename(from = from,  to = to)
}

for(i in 1:jpg_num){
  #exclude out all the duplicate
  file_name = jpg[i]
  
  if(!any(grepl("duplicate",file_name))){
    if(any(grepl("poster",file_name))){
      from_folder = "/protestify_data/images/addlexpts/"
      to_folder = "/protestify_data/images/poster/"
     
      my.file.rename(from = paste(from_folder,file_name,sep=""),to =paste(to_folder,file_name,sep=""))
    }
  }
}

for(i in 1:png_num){
  #exclude out all the duplicate
  file_name = png[i]
  
  if(!any(grepl("duplicate",file_name))){
    if(any(grepl("poster",file_name))){
      from_folder = "/protestify_data/images/addlexpts/"
      to_folder = "/protestify_data/images/poster/"
      
      my.file.rename(from = paste(from_folder,file_name,sep=""),to =paste(to_folder,file_name,sep=""))
    }
  }
}

