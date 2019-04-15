library(tidyverse)
library(pipeR)
library(gridExtra)
loadNamespace('cowplot')
setwd("/Volumes/kazcancer/T-DNA_search/HN00100341_hdd1/")
write_df= function(x, path, delim='\t', na='NA', append=FALSE, col_names=!append, ...) {
  file = if (grepl('gz$', path)) {
    gzfile(path, ...)
  } else if (grepl('bz2$', path)) {
    bzfile(path, ...)
  } else if (grepl('xz$', path)) {
    xzfile(path, ...)
  } else {path}
  utils::write.table(x, file,
                     append=append, quote=FALSE, sep=delim, na=na,
                     row.names=FALSE, col.names=col_names)
}


all_peak_tbl=read_tsv("all_peak_table_tidy.tsv",
                      col_names = c("sample","chr","start","end","coverage")) %>>%
  mutate(length=end-start) %>>%
  mutate(average_depth = coverage/length, start=start+1) %>>%
  dplyr::select(-coverage,-length)

test_max_value = function(.tbl){
  list=.tbl$average_depth
  MAX_p=pnorm(max(list),mean=mean(list),sd=sd(list),lower.tail = F)
  SEC_p=pnorm(sort(list)[23],mean=mean(list),sd=sd(list),lower.tail = F)
  spread(.tbl, key = sample, value=average_depth) %>>%
    mutate(max_p=MAX_p,second_p=SEC_p)
}

peak_tbl_with_p=all_peak_tbl %>>%
  group_by(chr,start,end)%>>%
  arrange(average_depth) %>>%
  mutate(max_sample=sample[24],second_sample=sample[23])%>>%
  ungroup()%>>%group_by(chr,start,end,max_sample,second_sample)%>>%
  nest() %>>%
  mutate(out = purrr::map(data, ~test_max_value(.))) %>>%
  select(-data)%>>%unnest()

peak_tbl_with_p %>>%arrange(max_p)%>>%
  write_df("all_peak_table_with_p.tsv")
######################################################################
#sampleごとの偏り
depth_by_sample = all_peak_tbl %>>%
  mutate(width=end-start+1) %>>%
  mutate(all_coverage = width * average_depth) %>>%
  group_by(sample)%>>%
  summarise(width=sum(width),all_coverage=sum(all_coverage)) %>>%
  mutate(average_depth = all_coverage/width) %>>%
  ggplot()+
  geom_bar(aes(x=sample,y=average_depth),stat = "identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -45,hjust = 0))
ggsave("peak_coverage_by_sample.pdf",depth_by_sample,width = 10,height = 5)

plot_by_chr = function(.tbl){
  chr = first(.tbl$chr_name)
  ggplot(.tbl)+
    geom_bar(aes(x=sample,y=average_depth),stat = "identity")+
    ggtitle(label = chr)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = -45,hjust = 0))
}
depth_by_chr = all_peak_tbl %>>%
  mutate(width=end-start+1) %>>%
  mutate(all_coverage = width * average_depth) %>>%
  group_by(sample,chr)%>>%
  summarise(width=sum(width),all_coverage=sum(all_coverage)) %>>%
  mutate(average_depth = all_coverage/width,chr_name=chr) %>>%
  ungroup()%>>%nest(-chr) %>>%
  mutate(plot = purrr::map(data, ~plot_by_chr(.)))
ggsave("peak_coverage_by_chr.pdf",gridExtra::marrangeGrob(depth_by_chr$plot,nrow = 4,ncol=1),width =8,height=12)









