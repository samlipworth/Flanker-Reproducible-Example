#this should be your clustering csv
flanker_output<-read_csv('~/oxa48/all') 
flanker_output<-flanker_output %>% mutate(window = map_chr(isolate, function(s) rev(strsplit(s, "_")[[1]])[3])) 
flanker_output$window<-as.numeric(flanker_output$window) 

x<-unique(flanker_output$window) %>% sort() 
flanker_output<-flanker_output %>% mutate(guuid = map_chr(isolate, function(s) rev(strsplit(s, "_")[[1]])[6]))

flanker_output<-flanker_output %>% mutate(gene = map_chr(isolate, function(s) rev(strsplit(s, "_")[[1]])[5])) 

gene<-select(flanker_output,guuid,gene)
flanker_output<-select(flanker_output,guuid,window,group) 
flanker_output<-flanker_output %>% pivot_wider(id_cols = guuid,names_from=window,values_from=group,names_sort=TRUE) 

flanker_output$ID<-flanker_output %>% group_indices(flanker_output[,2:50]) 
test<-select(flanker_output,guuid,ID) %>% distinct() 

flanker_output<-select(flanker_output,-ID) 
flanker_output<-flanker_output %>% pivot_longer(-guuid,names_to = "window",values_to="group")

flanker_output<-left_join(flanker_output,test,by=c("guuid"="guuid")) 
flanker_output$window<-as.numeric(flanker_output$window) 
flanker_output<-left_join(flanker_output,gene,by=c("guuid"="guuid")) 
#just check this looks sensible
p0<-ggplot(flanker_output) + aes(x=window,y=interaction(guuid,ID),fill=group) +geom_tile() + theme_minimal() + facet_wrap(~gene) 


#here should be your mashtree of plasmids
tree<-read.tree('~/oxa48/tree.tree')
tree<-midpoint.root(tree)
x<-data.frame(tree$tip.label)

x<-x %>% 
  mutate(id = map_chr(tree.tip.label, function(s) rev(strsplit(s, "_")[[1]])[1]))

tree$tip.label<-x$id
x<-left_join(x,gene,by=c("id"="guuid"))
x<-select(x,id,gene) %>% distinct()

g<-ggtree(tree)
p1<- g %<+% x + geom_tippoint(aes(color=gene))

names(flanker_output)<-c('id','window','group','cluster')

mlst<-read_csv('~/oxa48/mlst.csv')

mlst$mlst<-as.factor(mlst$mlst)


#make sure all gbk files in your working directory
filenames<-list.files("~/oxa48/",pattern="*.gbk",full.names = TRUE)
filenames<-str_replace_all(filenames,'/home/sam/oxa48//','/home/sam/oxa48/')
files<-as.list(filenames)

library(genbankr)
files2<-lapply(files, readGenBank)


for (i in 1:length(files2)){
  assign(paste(paste("df", i, sep="_"), "summary", sep="."), data.frame(genes(files2[[i]])))
  t<-paste(paste("df", i, sep="_"), "summary", sep=".")
  df<-get(t)
  df$seqnames<-paste('g',i,sep='')
  assign(paste(paste("df", i, sep="_"), "summary", sep="."), df)
}

files2<-lapply(ls(pattern='df_*'), get)

gggene_configure<-function(seqname){
  seqname<-select(seqname,seqnames,gene,start,end,strand)
  seqname$gene<-ifelse(is.na(seqname$gene),'hypothetical',seqname$gene)
  seqname$strand<-ifelse(seqname$strand=='+','forward','reverse')
  seqname$direction<-ifelse(seqname$strand =='forward',1,-1)
  
  return(seqname)
  
}

dfs<-lapply(ls(pattern = 'df_'),get)


all<-bind_rows(dfs)
all2<-gggene_configure(all)
all<-all2
all<-select(all,seqnames,gene,start,end,strand,direction)
names(all)<-c('molecule','gene','start','end','strand','direction')
filenames<-data.frame(filenames)
filenames<-filenames %>% 
  mutate(isolate = map_chr(filenames, function(s) rev(strsplit(s, "_")[[1]])[6])) 
filenames$isolate<-str_replace_all(filenames$isolate,'/home/sam/oxa48/','')
#filenames$isolate<-str_replace_all(filenames$isolate,'.gbk','')
filenames<-select(filenames,-filenames)
filenames$molecule<-paste('g',1:nrow(filenames),sep = '')
all<-left_join(all,filenames,by=c("molecule"="molecule"))
all<-select(all,isolate,gene,start,end,strand,direction)
names(all)<-c('id','gene','start','end','strand','direction')


p<-ggtree(tree) +geom_facet(panel = 'Alignment',mapping= aes(xmin=start,xmax=end,fill=gene),
                            data = all, geom = geom_motif, on = 'bla')


library(ggnewscale)
p1<-p +new_scale_fill()
p2<-facet_plot(p1,panel="MLST",data=mlst,geom=geom_tile,aes(x=1,fill=MIC))

z<-p2 + new_scale_fill()
p3<-facet_plot(z,panel="Flank group",data=flanker_output,geom=geom_tile,aes(x=1,fill=cluster))
flanker_output$window=-1*flanker_output$window
p4<-facet_plot(p3, panel="Flankergram",data=flanker_output,geom=geom_tile, aes(x=window,fill=group))

p5<- p4 + theme_bw() +
  xlim_tree(0.08) +
  xlim_expand(c(-4900,0),'Alignment') +
  xlim_expand(c(-4900,0),'Flankergram') + 
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())

y<-select(x,id)

plasmids<-read_tsv('~/oxa48/plasmids') %>% filter(id >=80 & cov >= 80) %>% select(guuid,gene)
pOXA<-filter(plasmids, gene =='IncL/M(pOXA-48)_1_pOXA-48')

y$`IncL/M(pOXA-48)`<-ifelse(y$id %in% pOXA$guuid,'Present','Absent' )

p<-p4 + new_scale_fill()

y$`IncL/M(pOXA-48)`<-as.factor(y$`IncL/M(pOXA-48)`)

p5<-facet_plot(p,panel='IncL/M(pOXA-48)',data=y,geom=geom_tile,aes(x=1,fill=`IncL/M(pOXA-48)`))

library(gggenes)
g<-ggplot(all,aes(xmin=start,xmax=end,y=id,fill=gene)) +
  geom_gene_arrow() + theme_genes() 

p5<- p5 + theme_bw() +
  xlim_tree(0.08) +
  xlim_expand(c(-4900,0),'Alignment') +
  xlim_expand(c(-4900,0),'Flankergram') + 
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())


#here we adjust the size of the columns
library(grid)
library(gtable)


#have to make mlst/flank panels smaller
gt = ggplot_gtable(ggplot_build(p5))
gtable_show_layout(gt) # will show you the layout - very handy function
gt # see plot layout in table format
gt$layout$l[grep('panel-1-2', gt$layout$name)] # you want to find the column specific to panel-2
gt$widths[9] = 0.3*gt$widths[9] # in this case it was colmun 7 - reduce the width by a half
gt$widths[11] = 0.3*gt$widths[11]
gt$widths[15] = 0.3*gt$widths[15]
grid.draw(gt) # plot with grid draw