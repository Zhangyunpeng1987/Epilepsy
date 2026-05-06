####海马体与颞叶皮层细胞类型间差异基因功能富集点图+热图-----------------
temp <- read.table("./temporal_lobe\\enrich\\celltype_EX_IN_Astro_Micro_Oli_Opcs_Endo.go.txt",sep = '\t',header = T)
colnames(temp) <- c(colnames(temp)[1],paste0(colnames(temp),"_Temp","")[2:15])

hippo <- read.table("./Hippo\\enrich\\celltype_EX_IN_Astro_Micro_Oli_Opcs_Endo.go.txt",sep = '\t',header = T)
colnames(hippo) <- c(colnames(hippo)[1],paste0(colnames(hippo),"_HIPPO","")[2:15])

exp <- merge(hippo, temp, by='Description', all=TRUE)


rownames(exp) <- exp$Description
exp <- exp[,c(2:ncol(exp))]
spot.size.function<-function(x){
  x= (-log(x))/25
}

spot.colour.function <- function(x) {
  colours <- rep("white", length(x));
  colours[x==exp$Excitatory_p_HIPPO] <- "#CC6699"; 
  colours[x==exp$Inhbitory_p_HIPPO] <- "#FF9966"; 
  colours[x==exp$Astro_p_HIPPO] <- "#99CC99";
  colours[x==exp$Microglial_p_HIPPO] <- "#0099CC";
  colours[x==exp$Oli_p_HIPPO] <- "#CCCCFF";
  colours[x==exp$OPC_p_HIPPO] <- "#FFCC00";
  colours[x==exp$Endo_p_HIPPO] <- "#CC9999";
  
  colours[x==exp$Excitatory_p] <- "#CC99CC96"; 
  colours[x==exp$Inhbitory_p] <- "#99CCCC"; 
  colours[x==exp$Astro_p] <- "#6699CC";
  colours[x==exp$Microglial_p] <- "#E58601";
  colours[x==exp$Oli_p] <- "#669999";
  colours[x==exp$OPC_p] <- "#B40F20";
  colours[x==exp$Endo_p] <- "#6CD398";
  return(colours);
}
color_cols <- c("#CC6699", "#FF9966", "#99CC99","#0099CC","#CCCCFF","#FFCC00","#CC9999","#CC99CC96","#99CCCC96","#6699CC96","#E5860196","#66999996","#B40F2096","#6CD39896")
#加图例
key.trans <- list(title = "cell type",space = "right",columns = 1,
                  points = list(pch = c(20,20,20,20,20,20,20,20,20,20,20,20,20),
                                col = color_cols,
                                cex=c(1,1,1,1)),
                  text = list(c("Hippo_Ex","Hippo_In","Hippo_Astro",
                                "Hippo_Micro","Hippo_Oli","Hippo_OPCs","Hippo_Endo",
                                "Temp_Ex","Temp_In","Temp_Astro",
                                "Temp_Micro","Temp_Oli","Temp_ODCs","Temp_Endo")),
                  
                  #lines = list(col = colors,lty = lines),
                  cex.title = 1,cex = .9)
create.dotmap(
  exp[,c(1,3,5,7,9,11,13,15,17,19,21,23,25,27)],bg.data = exp[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28)],
  pch = 20,na.spot.size=2,spot.size.function=spot.size.function,
  spot.colour.function=spot.colour.function,
  ylab.cex=1,yaxis.cex=0.5,xaxis.cex=0.5,
  colour.scheme=colorRampPalette(c("#80B1D3","#FDB462","#E59CC4", "#BC80BD"))(100),
  key = key.trans,total.colours = 10,colourkey = T,xaxis.rot = 45
)
