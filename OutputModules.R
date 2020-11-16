#### Output Modules #####



#### XML export ####
dend_to_text_wrapper <- function(
  R,
  i,
  filename
){
  fileConn<-file(filename)
  writeLines(paste("<root>",dend_to_text(R, i),"</root>"), fileConn)
  close(fileConn)
}

dend_to_text <- function(
  R,
  i
){
  me <- R$HCL$merge
  #Leaf case
  if(i <0){
    index <- abs(i) + dim(me)[1]
    lines <- paste("<Name>", R$annos$name[abs(i)], "</Name>", 
                   "<Type> Leaf </Type>", 
                   "<pval>", R$pvals[i], "</pval>")
  }
  
  else{
    child_add <- paste("<All_children>",paste(R$annos$name[abs(get_members(R, i))], collapse= ", "),'</All_children>')
    
    left_text <- dend_to_text(R, me[i,][1], list_children)
    right_text <-dend_to_text(R, me[i,][2], list_children)
    lines <- paste("<Name>", i, "</Name>", #paste(get_anno_name(R,i),  collapse =","), "</Name>", 
                   "<Type> Internal Node </Type>", 
                   "<pval>", R$pvals[i], "</pval>",
                   child_add,
                   "<Child>", left_text, "</Child>",
                   "<Child>", right_text, "</Child>")
  }
  lines
}

export_xml <- function(
  R
){
  i = dim(R$HCL$merge)[1]
  dend_to_text_wrapper(R, i, filename = "output.xml")
}

#### R structure export ####


#### csv export ####





