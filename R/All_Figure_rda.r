####FIG 2e#### 

library(cowplot)
library(ggplot2)
library(ggpubr)
my_theme <- function() {
  theme_half_open(font_size = 8, line_size = 0.2, rel_small = 1,rel_tiny = 1,rel_large = 1) 
  
}


load("data/Fig_2e.rda")

upf1_1<-ggplot() + 
  geom_segment(data=PK1, aes(x=Re_peak,y=upf,xend=Re_peak,yend=0),colour=ifelse(PK1$upf>5 ,"#eda83f","gray50"),size=0.5)+
  scale_y_continuous(name = "clone", breaks  =seq(0,8,4),limits = c(0,8),expand = c(0, 0))+
  scale_x_continuous(name = "distance to TC",  breaks =seq(0,50,10),limits = c(0,45),expand =  expand_scale(mult = c(0, 0.02), add = c(0, 0)))+
  ggtitle("upf1-5(35S:AtUPF1G851E, n=12)")+
  my_theme()



DE656<-ggplot() + 
  geom_segment(data=PK1, aes(x=Re_peak,y=DE,xend=Re_peak,yend=0),colour=ifelse(PK1$upf>5 ,"#eda83f","gray50"),size=0.5)+
  scale_y_continuous(name = "clone", breaks  =seq(0,8,4),limits = c(0,8),expand = c(0, 0))+
  scale_x_continuous(name = "distance to TC",  breaks =seq(0,50,10),limits = c(0,45),expand =  expand_scale(mult = c(0, 0.02), add = c(0, 0)))+
  ggtitle("upf1-5(35S:AtUPF1DE656AA, n=14)")+
  my_theme()




load("data/Fig_2e_1.rda")
PK1_B2upf <- ggplot() + 
  geom_segment(data=PK1_proto, aes(x=site, y=B2upf, xend=site, yend=0), colour=ifelse(PK1_proto$B2upf>1500, "#eda83f", "gray50"), size=0.5) +
  scale_y_continuous(name = "5'P end", breaks = seq(0,5000,1500), limits = c(0,5000), expand = c(0, 0)) +
  scale_x_continuous(name = "distance to TC", breaks = seq(0,50,10), limits = c(0,45), expand = expansion(mult = c(0, 0.02), add = c(0, 0))) +
  ggtitle("PK1: Protoplast upf1-1") +
  my_theme()




p2<-ggarrange( PK1_B2upf,upf1_1, DE656 , align = "v",
               ncol = 1, nrow = 3,common.legend = TRUE)
p2

ggsave("output/Fig_2e.pdf",p2,width =4,height =8,dpi = 300, bg='transparent',paper="A4",units = "cm")







####FIG 3a####
library(dplyr)

load("data/Fig_3a.rda")

summary(data.frame(rel$F_UTR_size, rel$CDS_size, rel$T_UTR_size))
utr5.SF <- median(rel$F_UTR_size, na.rm = T)/median(rel$CDS_size, na.rm = T)
utr3.SF <- median(rel$T_UTR_size, na.rm = T)/median(rel$CDS_size, na.rm = T)

library(scales)
library(ggplot2)
library(cowplot)
library(dplyr)

utr5.peak.dist <- rel[rel$position < 1, ]
cds.peak.dist <- rel[rel$position < 2 & rel$position >= 1, ]
utr3.peak.dist <- rel[rel$position >= 2, ]

utr5.peak.dist$rel_location <- rescale(utr5.peak.dist$position, to = c(1-utr5.SF, 1), from = c(0,1))
utr3.peak.dist$rel_location <- rescale(utr3.peak.dist$position, to = c(2, 2+utr3.SF), from = c(2,3))
peak.metagene.coord <- c(utr5.peak.dist$rel_location, cds.peak.dist$position, utr3.peak.dist$rel_location)

qplot1 <- ggplot() +
  geom_area(data = data.frame(x = peak.metagene.coord), aes(x = x),
            stat = "bin", bins = 100, color="#E8A680", fill = "#e8a680", alpha = 0.5,
            position = "identity", size=0.8)

density_data <- ggplot_build(qplot1)$data[[1]]

model.df <- rbind(
  data.frame(
    group=rep("UTR"),
    xmin=min(density_data$xmin),
    xmax=1,
    ymin=1.39*max(density_data$y)*0.98,
    ymax=1.39*max(density_data$y)*1.02),
  data.frame(
    group=rep("CDS"),
    xmin=1,
    xmax=2,
    ymin=1.39*max(density_data$y)*0.97,
    ymax=1.39*max(density_data$y)*1.03),
  data.frame(
    group=rep("UTR"),
    xmin=2,
    xmax=max(density_data$xmax),
    ymin=1.39*max(density_data$y)*0.98,
    ymax=1.39*max(density_data$y)*1.02))

my_theme <- function() {
  theme_half_open(font_size = 8, line_size = 0.25, rel_small = 1, rel_tiny = 1, rel_large = 1) 
}

output <- qplot1 + geom_vline(xintercept = 1:2, col = "gray50", linetype = "dashed") +
  list(
    geom_rect(data = model.df, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=group)),
    geom_text(data = data.frame(
      x=c((min(density_data$xmin)+1)/2, 1.5, (max(density_data$xmax)+2)/2),
      y=rep(1.39*max(density_data$y*1.03), 3),
      labels=c("5'UTR", "CDS", "3'UTR"), size=8),
      aes(x=x, y=y, label=labels), vjust=-0.5, size=8/ggplot2::.pt, check_overlap = TRUE))+
  scale_y_continuous(breaks = seq(0, 4000, 1000),
                     expand = expansion(mult=c(0, 0.05)), name = "Position count") +
  scale_fill_manual(values = c("grey70", "grey50"), breaks = c("UTR", "CDS")) +
  guides(fill="none") +
  labs(x=NULL) +
  coord_cartesian(clip = 'off') +
  my_theme() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.position = 'top',
        legend.justification = "center"
  )

raw_counts <- c(
  "5'UTR" = nrow(utr5.peak.dist),
  "CDS"   = nrow(cds.peak.dist),
  "3'UTR" = nrow(utr3.peak.dist)
)
total_count <- sum(raw_counts)

region_summary <- data.frame(
  region = names(raw_counts),
  count = as.numeric(raw_counts)
) %>%
  dplyr::mutate(
    percent = round(100 * count / total_count, 1),
    label = paste0("n=", count, "\n(", percent, "%)"),
    x = dplyr::case_when(
      region == "5'UTR" ~ (min(density_data$xmin) + 1) / 2,
      region == "CDS" ~ 1.5,
      region == "3'UTR" ~ (max(density_data$xmax) + 2) / 2
    ),
    y = 1.25 * max(density_data$y * 0.8)
  )

output1 <- output +
  geom_text(data = region_summary,
            aes(x = x, y = y, label = label),
            vjust = 0, size = 8 / ggplot2::.pt)

output1

ggsave("output/Fig_3a.pdf",output1,width =12,height = 5,dpi = 300,paper="A4",units = "cm")




####FIG 3b ####


load("data/Fig_3b.rda")
library(eulerr)
v <- venn(gene_sets)


pdf("output/Fig_3b.pdf", width = 7, height = 7)

plot(v,
     quantities = list(cex = 1.5), 
     labels = list(labels = c("Peak Regions", "upf1-1 protein"), cex = 1.2), 
     fills = list(
       fill = c(
         "Peaks" = "gray92",             
         "UPF1" = "gray92",            
         "Peaks&UPF1" = "#7a0a40"       
       ),
       alpha = 1 
     ),
     edges = list(col = "black", lwd = 1) 
)

dev.off() 






####FIG 3c####


library(ggplot2)
library(ggpattern)

my_theme <- function() {
  theme(panel.background =  element_blank(),
        panel.border =  element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background= element_blank(),
        legend.box.background= element_blank(),
        plot.background= element_blank(),
        strip.background = element_blank(),
        legend.text = element_text(face = c("italic", "plain"),size = 17),
        legend.key= element_rect(colour = NA, fill=NA),
        legend.title = element_text(size=17),
        plot.title = element_text(size=12, hjust = 0.5), 
        axis.title.x = element_text(size = 17),
        axis.title.y = element_blank(), 
        axis.text.x  = element_text(size = 17),
        axis.text.y  = element_blank(), 
        axis.line.x  = element_line(size = 1, colour = "black"),
        axis.line.y  = element_blank(), 
        axis.ticks.y = element_blank()
  )
}

categories <- c("Gray Group", "Purple Group") 
values <- c(36, 162) 
total <- sum(values)
percentages <- round(100 * values / total, 1)
labels <- paste0(values, " (", percentages, "%)")

df <- data.frame(
  Bar = "1",
  Category = factor(categories, levels = categories), # The levels dictate the stack order
  Value = values,
  Label = labels
)

pdf("output/Fig_3c_horizontal_ratio.pdf", width = 7, height = 2.5)


p <- ggplot(df, aes(x = Value, y = Bar, fill = Category)) +
  geom_col_pattern(
    aes(pattern = Category),       
    color = "black", 
    width = 1,
    pattern_color = "white",      
    pattern_fill = "white",        
    pattern_angle = 45,            
    pattern_density = 0.15,         
    pattern_spacing = 0.07         
  ) + 
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 5) + 
  
  scale_fill_manual(values = c("Purple Group" = "#7a0a40", "Gray Group" = "gray92")) +
  
  scale_pattern_manual(values = c("Purple Group" = "stripe", "Gray Group" = "none")) +
  
  labs(title = paste("Ratio of Categories (n =", total, ")"),
       x = "Count") +
  my_theme() +
  theme(legend.position = "none")

print(p)

dev.off()


library(dplyr)
library(tidyr)
library(ggplot2)


load("data/Fig_3c.rda")

gene_order <- input %>%
  distinct(gene, region_clean) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = region_clean, values_from = present, values_fill = 0) %>%
  arrange(desc(`5' UTR`), desc(CDS), desc(`3' UTR`), gene) %>%
  pull(gene)

all_genes <- unique(input$gene)
gene_order <- c(gene_order, setdiff(all_genes, gene_order))
total_genes <- length(gene_order)

cols <- 9
rows <- ceiling(total_genes / cols)

gene_map <- data.frame(gene = factor(gene_order, levels = gene_order)) %>%
  arrange(gene) %>% 
  mutate(
    id = row_number(),
    y = ceiling(id / cols),
    x = ifelse(id %% cols == 0, cols, id %% cols)
  )


build_feature_layer <- function(feature_label, target_region) {
  genes_in_region <- input %>%
    filter(region_clean == target_region) %>%
    pull(gene) %>%
    unique()
  
  count <- length(genes_in_region)
  pct <- round((count / total_genes) * 100, 1)
  display_label <- paste0(feature_label, "\n(", pct, "%)")
  
  gene_map %>%
    mutate(
      Feature = display_label,
      Fill_Group = ifelse(gene %in% genes_in_region, display_label, "Absent")
    )
}

plot_data <- bind_rows(
  build_feature_layer("5' UTR", "5' UTR"),
  build_feature_layer("CDS",    "CDS"),
  build_feature_layer("3' UTR", "3' UTR")
)


unique_labels <- unique(plot_data$Feature)

find_label <- function(pattern) grep(pattern, unique_labels, value = TRUE)
label_5 <- find_label("5' UTR")
label_C <- find_label("CDS")
label_3 <- find_label("3' UTR")

plot_data$Feature <- factor(plot_data$Feature, levels = c(label_5, label_C, label_3))

custom_colors <- c(
  "Absent" = "#e0e0e0" 
)
custom_colors[label_5] <- "#84ae92"  
custom_colors[label_C] <- "#4682b4"  
custom_colors[label_3] <- "#fba518"  



my_plot <- ggplot(plot_data, aes(x = x, y = y, fill = Fill_Group)) +
  geom_tile(color = "white", linewidth = 0.2) + 
  facet_wrap(~Feature) +
  scale_fill_manual(values = custom_colors) +
  coord_fixed() +
  theme_void() + 
  labs(title = paste("Gene Presence Map (n =", total_genes, ")")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "#333333", margin = margin(b = 15)),
    strip.text = element_text(size = 14, face = "bold", color = "#333333", margin = margin(b = 10)),
    legend.position = "none",
    plot.margin = margin(20, 20, 20, 20)
  )

print(my_plot)


ggsave("output/Fig_3c_binary.pdf", my_plot, width = 10, height = 7, dpi = 300, units = "cm", paper="A4")





####Fig 4a####
library(ggplot2)
library(cowplot)

my_theme <- function() {
  theme_half_open(font_size = 8, line_size = 0.2, rel_small = 1,rel_tiny = 1,rel_large = 1) 
  
}




load("data/Fig_4a.rda")

library(dplyr)
library(ggplot2)
library(ggpattern)

dens <- density(df$V3)

dens_df <- data.frame(x = dens$x, y = dens$y)

perc_data <- data.frame(
  p_50_130  = mean(df$V3 >= 50 & df$V3 <= 130) * 100,
  p_over350 = mean(df$V3 > 350) * 100
)

perc_data


perc_data <- data.frame(
  p_50_130  = mean( df$V3 < 133) * 100,
  p_over350 = mean(df$V3 > 350) * 100
)

perc_data <- data.frame(
  p_115  = mean( df$V3 <= 115) * 100,
  p_130  = mean( df$V3  >115 & df$V3 <= 140) * 100,
  p_over350 = mean(df$V3 > 350) * 100
)


perc_data




p <- ggplot(dens_df, aes(x, y)) +
  geom_line(color = "black", size = 0.2) +
  
  geom_area(fill = "#fffdf2", alpha = 0.6) +
  
  geom_area(
    data = subset(dens_df, x  <= 115),
    aes(x, y),
    fill = "#91c057"
    
  )+  geom_area(
    data = subset(dens_df, x >=115 & x <= 140),
    aes(x, y),
    fill = "#78a779"
    
  )+
  geom_area(
    data = subset(dens_df, x > 350),
    aes(x, y),
    fill =  "#4fb096"
  )+
  geom_vline(xintercept = c( 115,140, 350),
             color = "#c30d23", linetype = "dashed", size = 0.2) +
  # annotate("text", x = 50,  y = 0.007, label = "50",  size = 2.8) +
  annotate("text", x = 133, y = 0.007, label = "133", size = 2.8) +
  annotate("text", x = 350, y = 0.007, label = "350", size = 2.8) +
  labs(x = "Length", y = "Density") +
  my_theme() +
  theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.007)) +
  coord_cartesian(xlim = c(0, 600))


ggsave("output/Fig_4a.pdf", p, width = 10, height = 7, dpi = 300, units = "cm", paper="A4")





####Fig 4b####


library(gridExtra)
library(grid)
library(ggplot2)
library(ggpubr)
library(cowplot)


my_theme<- function() {
  
  theme_half_open(font_size = 8, line_size = 0.2, rel_small = 1,rel_tiny = 1,rel_large = 1) 
}


load("data/Fig_4b_1.rda")
load("data/Fig_4b_2.rda")
load("data/Fig_4b_3.rda")
load("data/Fig_4b_4.rda")




B2WT_long<-subset(DATA_B2_long,DATA_B2_long$V1=="WT")
B2upf_long<-subset(DATA_B2_long,DATA_B2_long$V1=="upf1-1")
B2WT_short1<-subset(DATA_B2_short1,DATA_B2_short1$V1=="WT")
B2upf_short1<-subset(DATA_B2_short1,DATA_B2_short1$V1=="upf1-1")
B2WT_short2<-subset(DATA_B2_short2,DATA_B2_short2$V1=="WT")
B2upf_short2<-subset(DATA_B2_short2,DATA_B2_short2$V1=="upf1-1")
B2WT_intron<-subset(DATA_B2_intron,DATA_B2_intron$V1=="WT")
B2upf_intron<-subset(DATA_B2_intron,DATA_B2_intron$V1=="upf1-1")

B2WT_data_long<-data.frame(B2WT_long[,2],B2WT_long[,3])
colnames(B2WT_data_long)<-c("position","sum")
B2upf1_data_long<-data.frame(B2upf_long[,2],B2upf_long[,3])
colnames(B2upf1_data_long)<-c("position","sum")

B2WT_data_short2<-data.frame(B2WT_short2[,2],B2WT_short2[,3])
colnames(B2WT_data_short2)<-c("position","sum")
B2upf1_data_short2<-data.frame(B2upf_short2[,2],B2upf_short2[,3])
colnames(B2upf1_data_short2)<-c("position","sum")
B2WT_data_short1<-data.frame(B2WT_short1[,2],B2WT_short1[,3])
colnames(B2WT_data_short1)<-c("position","sum")
B2upf1_data_short1<-data.frame(B2upf_short1[,2],B2upf_short1[,3])
colnames(B2upf1_data_short1)<-c("position","sum")

B2WT_data_intron<-data.frame(B2WT_intron[,2],B2WT_intron[,3])
colnames(B2WT_data_intron)<-c("position","sum")
B2upf1_data_intron<-data.frame(B2upf_intron[,2],B2upf_intron[,3])
colnames(B2upf1_data_intron)<-c("position","sum")


plot_long.model<-ggplot()+
  geom_text(aes(x=-20,y=1.9),label="n=285",size=8/ggplot2::.pt)+
  scale_x_continuous(limits = c(-50,100))+
  my_theme()+ 
  guides(color = guide_legend(override.aes = list(size = 2)))+
  theme_void()
plot_short1.model<-ggplot()+
  geom_text(aes(x=-20,y=1.9),label="n=1254",size=8/ggplot2::.pt)+
  scale_x_continuous(limits = c(-50,100))+
  my_theme()+ 
  guides(color = guide_legend(override.aes = list(size = 2)))+
  theme_void()

plot_short2.model<-ggplot()+
  geom_text(aes(x=-20,y=1.9),label="n=1258",size=8/ggplot2::.pt)+
  scale_x_continuous(limits = c(-50,100))+
  my_theme()+ 
  guides(color = guide_legend(override.aes = list(size = 2)))+
  theme_void()

plot_intron.model<-ggplot()+
  geom_text(aes(x=-20,y=1.9),label="n=1278",size=8/ggplot2::.pt)+
  scale_x_continuous(limits = c(-50,100))+
  my_theme()+ 
  guides(color = guide_legend(override.aes = list(size = 2)))+
  theme_void()

plot_long<-ggplot() + 
  geom_line(data = B2WT_data_long, mapping = aes(x = position, y = sum, group = 1,color="B2WT"),size=0.4)+
  geom_line(data = B2upf1_data_long, mapping = aes(x = position, y = sum, group = 1,color="B2upf1-1"),size=0.4)+
  scale_color_manual(values = c(#light yellow :#ffd60a
    'B2WT' = '#0b318f','B2upf1-1' = '#f5bc1d'
    
  ))+
  labs(x="distance from TC (nt)",y="5'P relative density")+
  scale_y_continuous(expand = c(0, 0),limits = c(0.5,2.5),breaks = seq(0.5,2.5,1)) +
  scale_x_continuous(
    breaks = seq(0, 100, 10), 
    labels = function(x) ifelse(x %% 20 == 0, x, "") 
  ) +
  my_theme()



plot_short1<-ggplot() + 
  geom_line(data = B2WT_data_short1, mapping = aes(x = position, y = sum, group = 1,color="B2WT"),size=0.4)+
  geom_line(data = B2upf1_data_short1, mapping = aes(x = position, y = sum, group = 1,color="B2upf1-1"),size=0.4)+
  scale_color_manual(values = c(
    'B2WT' = '#0b318f','B2upf1-1' = '#f5bc1d'
    
  ))+
  labs(color=NULL,x="distance from TC (nt)",y="5'P relative density")+
  scale_y_continuous(expand = c(0, 0), limits = c(0.5,2.5),breaks = seq(0.5,2.5,1)) +
  scale_x_continuous(
    breaks = seq(0, 100, 10), 
    labels = function(x) ifelse(x %% 20 == 0, x, "") 
  ) +
  my_theme()
plot_short2<-ggplot() + 
  geom_line(data = B2WT_data_short2, mapping = aes(x = position, y = sum, group = 1,color="B2WT"),size=0.4)+
  geom_line(data = B2upf1_data_short2, mapping = aes(x = position, y = sum, group = 1,color="B2upf1-1"),size=0.4)+
  scale_color_manual(values = c(
    'B2WT' = '#0b318f','B2upf1-1' = '#f5bc1d'
    
  ))+
  labs(color=NULL,x="distance from TC (nt)",y="5'P relative density")+
  scale_y_continuous(expand = c(0, 0), limits = c(0.5,2.5),breaks = seq(0.5,2.5,1)) +
  scale_x_continuous(
    breaks = seq(0, 100, 10), 
    labels = function(x) ifelse(x %% 20 == 0, x, "") 
  ) +
  my_theme()


plot_intron<-ggplot() + 
  
  geom_line(data = B2WT_data_intron, mapping = aes(x = position, y = sum, group = 1,color="B2WT"),size=0.4)+
  geom_line(data = B2upf1_data_intron, mapping = aes(x = position, y = sum, group = 1,color="B2upf1-1"),size=0.4)+
  scale_color_manual(values = c(
    'B2WT' = '#0b318f','B2upf1-1' = '#f5bc1d'
  ))+
  labs(x="distance from TC (nt)",y="5'P relative density")+
  scale_y_continuous(expand = c(0, 0), limits = c(0.5,2.5),breaks = seq(0.5,2.5,1)) +
  scale_x_continuous(
    breaks = seq(0, 100, 10), 
    labels = function(x) ifelse(x %% 20 == 0, x, "") 
  ) +
  my_theme()


plot_legend<-ggplot() + 
  
  geom_line(data = B2WT_data_intron, mapping = aes(x = position, y = sum, group = 1,color="WT"),size=0.8)+
  geom_line(data = B2upf1_data_intron, mapping = aes(x = position, y = sum, group = 1,color="upf1-1"),size=0.8)+
  labs(color=NULL)+
  guides(color=guide_legend(nrow = 1,reverse = TRUE))+
  scale_color_manual(values = c(
    'WT' = '#0b318f','upf1-1' = '#f5bc1d'
  ))+theme_bw(base_size = 12)+
  theme(legend.background = element_blank())+
  my_theme()

plot.legend<-get_legend(plot_legend)


final.plot<-plot_grid(plot.legend,NA,NA,NA,
                      plot_intron.model,plot_long.model,plot_short1.model,plot_short2.model,
                      plot_intron+labs(x="")+guides(color = FALSE),
                      plot_long+labs(y="")+guides(color = FALSE),
                      plot_short1+labs(x="")+guides(color = FALSE),
                      plot_short2+labs(y="",x="")+guides(color = FALSE),nrow=3,ncol = 4,rel_heights = c(0.1,0.1,0.3),align = 'vh')
final.plot<-plot_grid(plot.legend,plot_intron.model,plot_intron+guides(color = FALSE), plot_long.model,plot_long+guides(color = FALSE),plot_short2.model,plot_short2+guides(color = FALSE),plot_short1.model,plot_short1+guides(color = FALSE),nrow=9,ncol = 1,rel_heights = c(0.05,0.125,0.28,0.125,0.28,0.125,0.28,0.125,0.28),align = 'lr')




ggsave(filename = "output/Fig_4b.pdf",plot = final.plot,width = 5.5,height = 10.5,units = 'cm',paper="A4")




####Fig 4c####

my_theme<- function() {
  
  theme_half_open(font_size = 8, line_size = 0.2, rel_small = 1,rel_tiny = 1,rel_large = 1) 
}

load("data/Fig_4c.rda")

colnames(DATA_ratio)<-c("class","peak","all","ratio")



DATA_ratio$class <- factor(DATA_ratio$class, levels = c("Intron", "Long", "Short2", "Short1"))

plot<-ggplot(DATA_ratio, aes(x = class, y = ratio, fill = class)) +
  geom_bar(stat = "identity", width = 0.82) +
  scale_y_continuous(breaks = seq(0, 80, 20), limits = c(0, 80),
                     expand = expansion(mult = c(0, 0.05)), name = "Percentage") +
  scale_fill_manual(values = c("#254494", "#4FB096", "#78a779", "#91c057")) +
  geom_text(aes(label = paste0(sprintf("%.1f", ratio), "")),
            position = position_stack(vjust = 0.5), size = 8 / ggplot2::.pt, color = "white") +
  my_theme() +
  theme(legend.position = "none")

ggsave("output/Fig_4c_order.pdf",plot,width =7,height = 3.7,dpi = 300,paper="A4",units = "cm")



####Fig 4d####
library(ggplot2) 
library(drc)
library(cowplot)
library(ggrastr)



abline_to_segment_safe <- function(slope, intercept,
                                   x_range = c(-Inf, Inf),
                                   y_range = c(-Inf, Inf),
                                   color = "black", size = 0.5, ...) {

  x1 <- x_range[1]
  x2 <- x_range[2]
  y1 <- slope * x1 + intercept
  y2 <- slope * x2 + intercept
  

  points <- data.frame(
    x = c(x1, x2),
    y = c(y1, y2)
  )
  

  points$y[points$y < y_range[1]] <- y_range[1]
  points$y[points$y > y_range[2]] <- y_range[2]
  

  points$x <- (points$y - intercept) / slope
  
  geom_segment(
    aes(x = points$x[1], xend = points$x[2],
        y = points$y[1], yend = points$y[2]),
    color = color, size = size, ...
  )
}



load("data/Fig_4d.rda")



colnames(total) <- c("AVG", "degradome")
total_cut1 <- subset(total, total$AVG < 350)
total_cut2 <- subset(total, total$AVG >= 350)

S1 <- total_cut1$AVG
S2 <- total_cut2$AVG
v1 <- total_cut1$degradome
v2 <- total_cut2$degradome


linear_model <- lm(degradome ~ AVG, data = total_cut1)


lm_eqn <- function(model) {
  a  <- format(as.numeric(coef(model)[1]), digits =4) 
  b  <- format(as.numeric(coef(model)[2]), digits = 4)
  r2 <- format(summary(model)$r.squared, digits = 4)
  
  eq <- substitute(italic(y) == a + b * italic(x) * "," ~~ italic(R)^2 == r2, 
                   list(a = a, b = b, r2 = r2))
  
  return(as.character(as.expression(eq)))
}


my_formula_string <- lm_eqn(linear_model)
print(my_formula_string)

my_theme<- function() {
  
  theme_half_open(font_size = 8, line_size = 0.2, rel_small = 1,rel_tiny = 1,rel_large = 1) 
}





plotB1 <- ggplot() +
  
  geom_vline(xintercept = 350, lwd = 0.3, linetype = "dashed", color = "gray65") +
  

  abline_to_segment_safe(
    slope = 0.1694,
    intercept = 6.557,
    x_range = c(0, 350), 
    y_range = c(0, 80),
    color = "#c30d23",
    size = 0.3
  ) +
  

  scale_x_continuous(
    breaks = seq(0, 1000, 175),
    limits = c(0, 550),
    name = "Average 3'UTR Length (nt)",
    expand = expansion(mult = c(0, 0.05))
  ) +
  

  scale_y_continuous(
    breaks = seq(0, 80, 20),
    limits = c(0, 85), 
    labels = function(x) paste0(x, ""), 
    name = "UIF detection rate \n(count/100 transcripts)", 
    expand = expansion(mult = c(0, 0.05))
  ) +
  

  geom_point_rast(aes(x = S1, y = v1), 
                  color = "#208091", 
                  alpha = 0.8, 
                  size = 2.4,
                  shape = 16) +
  geom_point_rast(aes(x = S2, y = v2), 
                  color = "#208091", 
                  alpha = 0.8, 
                  size = 2.4,
                  shape = 16) +
  my_theme()


ggsave("output/Fig_4d_in_350.pdf", plotB1, width = 10, height = 7.6, dpi = 300, units = "cm", paper = "A4")



####Fig 4e####

library(tidyverse)
library(ggplot2)
library(dplyr)
library(purrr)
library(stringr)
library(tidyverse)
library(patchwork)
library(stringr)
library(cowplot)


load("data/Fig_4e.rda")


my_plot_theme <- function() {
  theme_minimal() +
    theme_half_open(font_size = 8, line_size = 0.2, rel_small = 1, rel_tiny = 1, rel_large = 1) +
    theme(strip.text.y = element_blank())
}

redundant_go_terms <- c(
  "mRNA metabolic process", "mRNA binding", "RNA biosynthetic process", 
  "spliceosomal snRNP complex", "mRNA processing", "RNA splicing", 
  "respiratory chain complex", "oxidoreductase complex",
  "kinase activity", "signal transduction", "Cul4-RING E3 ubiquitin ligase complex",
  "plastid", "plastid membrane", "chloroplast thylakoid", 
  "chloroplast thylakoid membrane", "organelle envelope", "chloroplast envelope",
  "vacuole", "plant-type vacuole"
)


final_uif_data <- raw_uif_data %>%
  mutate(
    neg_log10p = as.numeric(MinusLog10AdjP),
    fold = as.numeric(OddsRatio),
    GO_term_clean = str_remove_all(GO_term, "\\s*\\(?GO:\\d+\\)?[:~\\s-]*"),
    GO_term_clean = str_wrap(str_trim(GO_term_clean), width = 55),
    GO_class = str_to_title(GO_class)
  ) %>%
  mutate(GO_term = GO_term_clean) %>%
  filter(!is.na(neg_log10p)) %>%
  filter(neg_log10p > -log10(0.05)) %>%
  filter(GO_term_clean %in% redundant_go_terms) %>%
  mutate(OddsRatio = ifelse(is.infinite(OddsRatio), 15, OddsRatio))



make_uif_plot <- function(df, group_filter, x_max, show_x = FALSE) {
  plot_df <- df %>% 
    filter(Group == group_filter) %>%
    arrange(GO_class, neg_log10p) %>%
    mutate(GO_term = factor(GO_term, levels = unique(GO_term)))
  
  if(nrow(plot_df) == 0) return(NULL)
  
  ggplot(plot_df, aes(x = neg_log10p, y = GO_term)) +
    geom_point(aes(size = OddsRatio, color = GO_class)) +
    facet_grid(GO_class ~ ., scales = "free_y", space = "free_y") +
    scale_x_continuous(limits = c(0.5,4), breaks = c(1,2.5,4)) +
    scale_size_continuous(range = c(0, 4.5), limits = c(0, 20),breaks = c(3, 6, 9, 15),
                          labels = c("3", "6", "9", "Inf")) +
     scale_color_manual(values = c(
      "Biological Process" = "#FF6B6B", 
      "Molecular Function" = "#FBA518", 
      "Cellular Component" = "steelblue"
    )
    )+
    labs(
      title = paste("Group:", group_filter),
      x = if(show_x) expression(-log[10] * "(Adjusted P-Value)") else NULL, 
      y = NULL,
      size = "Odds Ratio"
    ) +
    my_plot_theme()
}


global_x_max <- max(final_uif_data$neg_log10p) + 1

p1 <- make_uif_plot(final_uif_data, "Intron", global_x_max)
p2 <- make_uif_plot(final_uif_data, "Long", global_x_max)
p3 <- make_uif_plot(final_uif_data, "Short I", global_x_max)
p4 <- make_uif_plot(final_uif_data, "Short II", global_x_max, show_x = TRUE)


final_combined <- 
  p1 / p2 / p4 / p3 +
  plot_layout(
    guides = "collect",
    heights = c(
      nrow(filter(final_uif_data, Group == "Intron")),
      nrow(filter(final_uif_data, Group == "Long")),
      nrow(filter(final_uif_data, Group == "Short II")),
      nrow(filter(final_uif_data, Group == "Short I"))
    )
  ) +
  plot_annotation(
    theme = theme(legend.position = "bottom")
  )



ggsave("output/Fig_4e.pdf", 
       plot = final_combined, 
       width = 4, 
       height = 5.75, 
       limitsize = FALSE)





####Fig 5b 5d####



library(dplyr)
library(broom)
library(car)  
library(ggplot2)

my_theme <- function() {
  theme_half_open(font_size = 8, line_size = 0.2, rel_small = 1,rel_tiny = 1,rel_large = 1) 
}


load("data/Fig_5b5d6h_new.rda")

All$sample <- factor(All$sample, levels = unique(All$sample[order(All$test)]))

All$gene<-factor(All$gene, levels = c("NOS", "abc","SAP5","3143", "NPF","MPK6","18500","PRA1","PRA1(T)","CKL2","CKL2(T)"))


fold_changes <- All %>%
  group_by(gene) %>%
  summarize(
    FC = mean(value[sample != "WT"]) / mean(value[sample == "WT"])
  )



All_gc <- All %>% filter(test == "GC")

All_other <- All %>% filter(test != "GC")


All_other_L1 <- All_other %>% filter(test != "Length_R")
All_other_L2 <- All_other %>% filter(test == "Length_R")
All_other_L2_18500<-All_other_L2 %>% filter(gene=="18500")
All_other_L2_MPK6<-All_other_L2 %>% filter(gene=="MPK6")



run_tests <- function(data) {
  
  lev <- leveneTest(value ~ sample, data = data)
  var_equal <- lev$`Pr(>F)`[1] > 0.05
  
  t_res <- t.test(value ~ sample, data = data, var.equal = var_equal)
  
  data.frame(
    levene_p   = lev$`Pr(>F)`[1],
    var_equal  = var_equal,
    t_stat     = t_res$statistic,
    df         = t_res$parameter,
    p_value    = t_res$p.value,
    mean_diff  = diff(rev(t_res$estimate)) 
  )
}


run_tests <- function(data) {
  lev <- leveneTest(value ~ sample, data = data)
  var_equal <- lev$`Pr(>F)`[1] > 0.05
  t_res <- t.test(value ~ sample, data = data, var.equal = var_equal)
  
  data.frame(
    gene      = unique(data$gene),
    p_value   = t_res$p.value
  )
}

pvals <- All_other %>%
  group_by(gene) %>%
  do(run_tests(.))




results <- All_other %>%
  group_by(gene) %>%
  do(run_tests(.))


pvals <- pvals %>%
  mutate(
    signif_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )



label_pos <- All_other %>%
  group_by(gene) %>%
  summarise(y_pos = max(value) + 0.3)

pvals <- left_join(pvals, label_pos, by = "gene")
pvals_fc <- pvals %>%
  left_join(fold_changes, by = c("gene"))

pvals_fc_cut<-pvals_fc[1:5,]


A_fig<-ggplot(All_other_L1, aes(x = sample, y = value, color = sample, fill = sample)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.6, size = 0.2) +
  geom_jitter(width = 0, size = 0.5, color = "black",shape = 20) +
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.2, size = 0.1) +
  facet_wrap(~gene, nrow = 1) +
  geom_text(
    data = pvals_fc_cut,
    aes(x = 1.5, y = 3.6, label = paste0(signif_label, "\n", round(FC, 3))),
    inherit.aes = FALSE,
    size = 2.5, lineheight = 0.65  
  ) +
  labs(
    y = "Expression level",
    x = ""
  ) +
  my_theme() +
  scale_y_continuous(expand = c(0,0), limits = c(0,4.2)) +
  scale_fill_manual(values = c("WT" = "#5EABD6", "upf1-5" = "#E14434")) +
  scale_color_manual(values = c("WT" = "black", "upf1-5" = "black")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1,size = 8)
  )



ggsave("output/Fig_5b_new.pdf",A_fig,width =8,height =3.6,dpi = 300,paper="A4",units = "cm")



pvals_fc_cut<-pvals_fc[7,]
fig18500<-ggplot(All_other_L2_18500, aes(x = sample, y = value, color = sample, fill = sample)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.6, size = 0.2) +
  geom_jitter(width = 0, size = 0.5, color = "black",shape = 20) +
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.2, size = 0.1) +
  facet_wrap(~gene, nrow = 1) +
  geom_text(
    data = pvals_fc_cut,
    aes(x = 1.5, y = 1.8, label = paste0(signif_label, "\n", round(FC, 3))),
    inherit.aes = FALSE,
    size = 2.5, lineheight = 1
  ) +
  labs(
    y = "Expression level",
    x = ""
  ) +
  my_theme() +
  scale_y_continuous(expand = c(0,0), limits = c(0,4.2)) +
  scale_fill_manual(values = c("WT" = "#5EABD6", "upf1-5" = "#E14434")) +
  scale_color_manual(values = c("WT" = "black", "upf1-5" = "black")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1,size = 8)
  )

pvals_fc_cut<-pvals_fc[6,]
figMPK6<-ggplot(All_other_L2_MPK6, aes(x = sample, y = value, color = sample, fill = sample)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.6, size = 0.2) +
  geom_jitter(width = 0, size = 0.5, color = "black",shape = 20) +
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.2, size = 0.1) +
  facet_wrap(~gene, nrow = 1) +
  geom_text(
    data = pvals_fc_cut,
    aes(x = 1.5, y = 1.8, label = paste0(signif_label, "\n", round(FC, 3))),
    inherit.aes = FALSE,
    size = 2.5, lineheight = 1
  ) +
  labs(
    y = "Expression level",
    x = ""
  ) +
  my_theme() +
  scale_y_continuous(expand = c(0,0), limits = c(0,4.2)) +
  scale_fill_manual(values = c("WT" = "#5EABD6", "upf1-5" = "#E14434")) +
  scale_color_manual(values = c("WT" = "black", "upf1-5" = "black")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1,size = 8)
  )


fig18500 <- fig18500 
figMPK6  <- figMPK6  

figB <- plot_grid(
  fig18500, figMPK6,
  nrow = 2, ncol = 1,
  align = "v"
)


ggsave("output/Fig_5d_18500_new.pdf",fig18500,width =2.9,height =3.4,dpi = 300,paper="A4",units = "cm")
ggsave("output/Fig_5d_MPK6_new.pdf",figMPK6,width =2.9,height =3.4,dpi = 300,paper="A4",units = "cm")





####Fig 6b####



my_theme <- function() {
  theme_half_open(font_size = 8, line_size = 0.2, rel_small = 1,rel_tiny = 1,rel_large = 1)}


load("data/Fig_6b_1.rda")
load("data/Fig_6b_2.rda")


library(tidyverse)
library(recipes)
library(caret)
library(glmnet)
library(ROSE)
library(ggplot2)

set.seed(3750)
n_repeat <-100
coef_list <- list()

column_name_all <- colnames(all)[!colnames(all) %in% "peak"]


rec <- recipe(peak ~ ., data = all) %>%
  update_role(peak, new_role = "outcome") %>%
  update_role(all_of(column_name_all), new_role = "predictor") %>%
  step_YeoJohnson(all_numeric_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  prep()

scaled.all <- bake(rec, new_data = NULL)


for (i in 1:n_repeat) {
  balanced_data <- ovun.sample(peak ~ ., data = scaled.all, method = "under", seed = i)$data
  train_idx <- createDataPartition(balanced_data$peak, p = 0.7, list = FALSE)
  train_data <- balanced_data[train_idx, ]
  x_train <- as.matrix(train_data[, -ncol(train_data)])
  y_train <- as.factor(train_data[, ncol(train_data)])
  cv_model <- cv.glmnet(
    x = x_train,
    y = y_train,
    alpha = 1,
    family = "binomial",
    grouped = FALSE
  )
  
  best_lambda <- cv_model$lambda.min
  coef_i <- coef(cv_model, s = best_lambda)
  coef_df <- as.data.frame(as.matrix(coef_i))
  colnames(coef_df) <- paste0("rep", i)
  coef_list[[i]] <- coef_df
}

coef_all <- do.call(cbind, coef_list)
coef_stats <- data.frame(
  feature = rownames(coef_all),
  mean = rowMeans(coef_all),
  sd = apply(coef_all, 1, sd)
) %>%
  filter(feature != "(Intercept)") %>%
  mutate(direction = mean < 0)

coef_stats_top <- coef_stats %>%
  slice_max(order_by = abs(mean), n =10)


p<-ggplot(coef_stats_top, aes(x = abs(mean), y = reorder(feature, abs(mean)), fill = direction)) +
  geom_bar(stat = "identity", color = "black", width = 0.75, size = 0.25) +   # thinner bars
  geom_errorbarh(aes(xmin = abs(mean) - sd, xmax = abs(mean) + sd), height = 0.2) +
  scale_fill_manual(values = c("FALSE" = "#f7c2ad", "TRUE" = "#698f9e")) +
  scale_x_continuous(
    expand = expansion(mult = 0) 
  ) +
  labs(
    x = "Absolute Coefficient (± SD)",
    y = "Feature"
  ) +
  theme_minimal() +
  my_theme() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

ggsave("output/Fig_6b.pdf",p,width =3,height = 2,dpi = 300, bg='transparent',paper="A4")




####Fig 6c#####
library(ggplot2)
library(cowplot)

my_theme <- function() {
  theme_half_open(font_size = 8, line_size = 0.2, rel_small = 1,rel_tiny = 1,rel_large = 1) 
  
}

load("data/Fig_6c.rda")
density(df$V2)


q10 <- quantile(df$V2, 0.1)
q90 <- quantile(df$V2, 0.9)


dens <- density(df$V2)


dens_df <- data.frame(x = dens$x, y = dens$y)

perc_data <- data.frame(
  p_50_130  = mean(df$V3 >= 50 & df$V3 <= 130) * 100,
  p_over350 = mean(df$V3 > 350) * 100
)

perc_data



p <- ggplot(df, aes(x = V2)) +
  
  geom_area(
    data = subset(dens_df, x < q10),
    aes(x, y),
    fill = "#ccb6d6",   
    alpha = 1,
    size = 0.1
  ) +
  
  geom_area(
    data = subset(dens_df, x > q90),
    aes(x, y),
    fill = "#844b96",  
    alpha = 1,
    size = 0.1
  ) +
  geom_density(fill = "white", alpha = 0.1) +
  geom_vline(xintercept = q10, color = "#c30d23", linetype = "dashed", size = 0.3) +
  geom_vline(xintercept = q90, color = "#c30d23", linetype = "dashed", size = 0.3) +
  annotate("text", x = 0.4, y = 6, label = "10%", hjust = 0, size = 2.8) +
  annotate("text", x = 0.16, y = 6, label = "10%", hjust = 0, size = 2.8) +
  annotate("text", x = q10, y = 12, label = round(q10, 3),
           hjust = 0.5, angle = 0, color = "black", size = 2.8) +
  annotate("text", x = q90, y = 12, label = round(q90, 3),
           hjust = 0.5, angle = 0, color = "black", size = 2.8) +
  labs(x = "GC", y = "Density") +
  scale_y_continuous(expand = c(0,0), limits = c(0,12)) +
  my_theme() +
  theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10))



ggsave("output/Fig_6c.pdf",p,width =1.8,height = 2,dpi = 300,paper="A4")





####Fig 6d####

library(gridExtra)
library(grid)
library(ggplot2)
library(ggpubr)
library(cowplot)


my_theme<- function() {
  
  theme_half_open(font_size = 8, line_size = 0.2, rel_small = 1,rel_tiny = 1,rel_large = 1) 
}

load("data/Fig_6d_1.rda")
load("data/Fig_6d_2.rda")

B2WT_high<-subset(DATA_B2_high,DATA_B2_high$V1=="WT")
B2upf_high<-subset(DATA_B2_high,DATA_B2_high$V1=="upf1-1")
B2WT_low<-subset(DATA_B2_low,DATA_B2_low$V1=="WT")
B2upf_low<-subset(DATA_B2_low,DATA_B2_low$V1=="upf1-1")


colnames(B2WT_high)<-c("sample","position","sum")
colnames(B2WT_low)<-c("sample","position","sum")
colnames(B2upf_high)<-c("sample","position","sum")
colnames(B2upf_low)<-c("sample","position","sum")


plot_high<-ggplot() + 
  geom_line(data = B2WT_high, mapping = aes(x = position, y = sum, group = 1,color="B2WT"),size=0.4)+
  geom_line(data = B2upf_high, mapping = aes(x = position, y = sum, group = 1,color="B2upf1-1"),size=0.4)+
  scale_color_manual(values = c(#light yellow :#ffd60a
    'B2WT' = '#0b318f','B2upf1-1' = '#f5bc1d'
  ))+
  labs(x="distance from TC(nt)",y="5'P relative density")+
  scale_y_continuous(expand = c(0, 0), limits = c(0.5,2.6)) +
  scale_x_continuous(breaks=seq(0,100,20)) + guides(color = guide_legend(override.aes = list(size = 2)))+
  annotate("text", x = 80, y = 2.1, label = "N = 858", hjust = 1, size = 3.5) +
  my_theme()


plot_low<-ggplot() + 
  geom_line(data = B2WT_low, mapping = aes(x = position, y = sum, group = 1,color="B2WT"),size=0.4)+
  geom_line(data = B2upf_low, mapping = aes(x = position, y = sum, group = 1,color="B2upf1-1"),size=0.4)+
  scale_color_manual(values = c(
    'B2WT' = '#0b318f','B2upf1-1' = '#f5bc1d'
  ))+
  labs(color=NULL,x="distance from TC(nt)",y="5'P relative density")+
  scale_y_continuous(expand = c(0, 0), limits = c(0.5,2.6)) +
  scale_x_continuous(breaks=seq(0,100,20)) +guides(color = guide_legend(override.aes = list(size = 2)))+
  annotate("text", x = 80, y = 2.1, label = "N = 792", hjust = 1, size = 3.5) +
  my_theme()



GC_save<- ggarrange( plot_low,plot_high,  align = "h",
                     font.label = list(size = 18),
                     ncol =2, nrow = 1,common.legend = TRUE)



ggsave(filename = "output/Fig_6d.pdf",plot =GC_save,width = 10,height = 5.5,units = 'cm',paper="A4")



####Fig 6e####

load("data/Fig_6e.rda")

colnames(DATA_ratio)<-c("class","peak","all","ratio")


my_theme<- function() {
  
  theme_half_open(font_size = 8, line_size = 0.2, rel_small = 1,rel_tiny = 1,rel_large = 1) 
}

DATA_ratio$class <- factor(DATA_ratio$class, levels = c("Low", "High"))

plot <- ggplot(DATA_ratio, aes(x = class, y = ratio, fill = class)) +
  geom_bar(stat = "identity", width = 0.85) +
  scale_y_continuous(breaks = seq(0, 80, 20), limits = c(0, 80),
                     expand = expansion(mult = c(0, 0.05)), name = "Percentage") +
  scale_fill_manual(values = c("High" = "#844B96", "Low" = "#B189BC")) +
  geom_text(aes(label = paste0(sprintf("%.1f", ratio), "%")),
            position = position_stack(vjust = 0.5), size=8/ggplot2::.pt, color = "white") +
  theme(legend.position = "none") +
  my_theme()



ggsave("output/Fig_6e.pdf",plot,width =5.4,height =5.3,dpi = 300,paper="A4",units = "cm")








####Fig 6f####

abline_to_segment_safe <- function(slope, intercept,
                                   x_range = c(-Inf, Inf),
                                   y_range = c(-Inf, Inf),
                                   color = "black", size = 0.5, ...) {
  
  x1 <- x_range[1]
  x2 <- x_range[2]
  y1 <- slope * x1 + intercept
  y2 <- slope * x2 + intercept
  
  
  points <- data.frame(
    x = c(x1, x2),
    y = c(y1, y2)
  )
  
  
  points$y[points$y < y_range[1]] <- y_range[1]
  points$y[points$y > y_range[2]] <- y_range[2]
  
  
  points$x <- (points$y - intercept) / slope
  
  geom_segment(
    aes(x = points$x[1], xend = points$x[2],
        y = points$y[1], yend = points$y[2]),
    color = color, size = size, ...
  )
}


load("data/Fig_6f.rda")
colnames(total)<-c("AVG","degradome")

S1 <- c(total$AVG)
v1 <- c(total$degradome)

lm_eqn <- function(df){
  m <- lm(degradome ~ AVG, df);
  eq <- substitute(italic(y) == ~a + ~b %.% italic(x)*","~~italic(R)^2~"="~r2, 
                   list(a = sprintf(format(coef(m)[1], digits = 2)), 
                        b = sprintf(format(coef(m)[2], digits = 2)), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

lm_eqn(total)
my_theme<- function() {
  
  theme_half_open(font_size = 8, line_size = 0.2, rel_small = 1,rel_tiny = 1,rel_large = 1) 
}



plotB1 <- ggplot() +
  abline_to_segment_safe(
    slope = 352,
    intercept = -73,
    x_range = c(0.1, 0.55),  
    y_range = c(0, 80),       
    color = "#e50012",
    size = 0.3
  ) +
  geom_point_rast(aes(x = S1, y = v1), color = "#208091", alpha = 0.8, size = 0.7, shape = 16, raster.dpi = getOption("ggrastr.default.dpi", 900)) +
  scale_x_continuous(breaks = seq(0.1, 0.55, 0.1), limits = c(0.1, 0.55),
                     name = "Average GC content", expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(breaks = seq(0, 80, 20), limits = c(0, 80),
                     name =  "UIF detection rate \n(count/100 transcripts)", 
                     expand = expansion(mult = c(0, 0.05))) +
  my_theme()

plotB1

ggsave("output/Fig_6f_d100.pdf",plotB1,width = 4.8,height = 5.3,dpi = 300,units = "cm",paper="A4")




####Fig 6h####


library(dplyr)
library(broom)
library(car)  
library(ggplot2)
library(cowplot)

my_theme <- function() {
  theme_half_open(font_size = 8, line_size = 0.2, rel_small = 1,rel_tiny = 1,rel_large = 1) 
}

load("data/Fig_5b5d6h_new.rda")
All$sample <- factor(All$sample, levels = unique(All$sample[order(All$test)]))
All$gene<-factor(All$gene, levels = c("NOS", "abc","SAP5","3143", "NPF","MPK6","18500","PRA1","PRA1(T)","CKL2","CKL2(T)"))


fold_changes <- All %>%
  group_by(gene) %>%
  summarize(
    FC = mean(value[sample != "WT"]) / mean(value[sample == "WT"])
  )


All_gc <- All %>% filter(test == "GC")

All_other <- All %>% filter(test != "GC")

run_tests <- function(data) {
  
  lev <- leveneTest(value ~ sample, data = data)
  var_equal <- lev$`Pr(>F)`[1] > 0.05
  
  t_res <- t.test(value ~ sample, data = data, var.equal = var_equal)
  
  data.frame(
    levene_p   = lev$`Pr(>F)`[1],
    var_equal  = var_equal,
    t_stat     = t_res$statistic,
    df         = t_res$parameter,
    p_value    = t_res$p.value,
    mean_diff  = diff(rev(t_res$estimate)) 
  )
}


run_tests <- function(data) {
  lev <- leveneTest(value ~ sample, data = data)
  var_equal <- lev$`Pr(>F)`[1] > 0.05
  t_res <- t.test(value ~ sample, data = data, var.equal = var_equal)
  
  data.frame(
    gene      = unique(data$gene),
    p_value   = t_res$p.value
  )
}

pvals <- All_other %>%
  group_by(gene) %>%
  do(run_tests(.))



results <- All_other %>%
  group_by(gene) %>%
  do(run_tests(.))


pvals <- pvals %>%
  mutate(
    signif_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )



label_pos <- All_other %>%
  group_by(gene) %>%
  summarise(y_pos = max(value) + 0.3)

pvals <- left_join(pvals, label_pos, by = "gene")
pvals_fc <- pvals %>%
  left_join(fold_changes, by = c("gene"))


pvals <- All_gc %>%
  group_by(gene) %>%
  do(run_tests(.))



results <- All_gc %>%
  group_by(gene) %>%
  do(run_tests(.))


pvals <- left_join(pvals, label_pos, by = "gene")
pvals_fc <- pvals %>%
  left_join(fold_changes, by = c("gene"))



pvals <- left_join(pvals, label_pos, by = "gene")
pvals_fc <- pvals %>%
  left_join(fold_changes, by = c("gene"))


pvals_fc <- pvals_fc %>%
  mutate(
    signif_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )



fig_pra1 <- ggplot(
  All_gc,
  aes(x = sample, y = value, color = sample, fill = sample)
) +
  geom_bar(
    stat = "summary",
    fun = "mean",
    width = 0.6
  ) +
  
  geom_jitter(
    width = 0,
    size = 0.5,
    color = "black",
    shape = 20
  ) +
  
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    size = 0.1
  ) +
  
  facet_wrap(~gene, nrow = 1) +
  
  geom_text(
    data = pvals_fc,
    aes(
      x = 1.5,
      y = 2.5,
      label = paste0(signif_label, "\n", round(FC, 3))
    ),
    inherit.aes = FALSE,
    size = 2.5,
    lineheight = 0.65
  ) +
  labs(
    y = "Expression level",
    x = ""
  ) +
  
  my_theme() +
  
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 3)) +
  
  scale_fill_manual(
    values = c("WT" = "#5EABD6", "upf1-5" = "#E14434")
  ) +
  scale_color_manual(
    values = c("WT" = "black", "upf1-5" = "black")
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1, size = 8)
  )

fig_pra1

ggsave("output/Fig_6h_new.pdf",fig_pra1,width =6.4,height =3.4,dpi = 300,paper="A4",units = "cm")




####Fig 7b####

library(dplyr)
library(ggplot2)

load("data/Fig_7b.rda")


conf_data <- all_predictions %>%
  mutate(true_label = factor(true_label, levels = c(0, 1)),
         predicted_label = factor(predicted_label, levels = c(0, 1))) %>%
  group_by(run, true_label, predicted_label) %>%
  summarise(n = n(), .groups = "drop")


overall_acc_per_run <- conf_data %>%
  group_by(run) %>%
  summarise(
    total_correct = sum(n[true_label == predicted_label]),
    total_n = sum(n),
    accuracy = 100 * total_correct / total_n,
    .groups = "drop"
  )


overall_accuracy_summary <- overall_acc_per_run %>%
  summarise(
    mean_acc = mean(accuracy),
    sd_acc = sd(accuracy)
  )

cat(sprintf("Overall Accuracy: %.1f%% ± %.1f%%\n", 
            overall_accuracy_summary$mean_acc, 
            overall_accuracy_summary$sd_acc))


conf_summary <- conf_data %>%
  group_by(run, true_label) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup() %>%
  group_by(true_label, predicted_label) %>%
  summarise(
    mean_percent = mean(percent),
    sd_percent = sd(percent),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("%.1f\n± %.1f%%", mean_percent, sd_percent))


con <- ggplot(conf_summary, aes(x = predicted_label, y = true_label, fill = mean_percent)) +
  geom_tile(color = "black") +
  geom_text(aes(label = label), size = 8 / ggplot2::.pt) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    x = "Predicted Label",
    y = "True Label",
    fill = "Mean %"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 8 ),
    axis.title.y = element_text(size = 8 ),
    axis.text.x = element_text(size = 8 ),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


print(con)
ggsave("output/Fig_7b.pdf",con,width =6,height =6,dpi = 300,paper="A4",units = "cm")

####Fig 7c#### 


install.packages("PRROC") 
library(PRROC)
library(pROC)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(readr)

my_theme <- function() {
  theme_half_open(font_size = 8, line_size = 0.2, rel_small = 1,rel_tiny = 1,rel_large = 1) 
  
}


run_colors <- colorRampPalette(c("#b8dcf5", "#041747"))(7)

load("data/Fig_7c.rda")


pr_df_all_filtered <- pr_df_all %>% filter(Recall > 0.01)

pr_df_all_filtered$Run <- factor(pr_df_all_filtered$Run, levels = paste0("Run_", 1:length(unique(pr_df_all_filtered$Run))))
all_run_labels <- unique(pr_df_all_filtered$Run)

common_recalls <- seq(0.01, 1, length.out = 1000)

interpolated_list <- lapply(split(pr_df_all_filtered, pr_df_all_filtered$Run), function(df) {
  approx(df$Recall, df$Precision, xout = common_recalls, rule = 2)$y
})

mean_precision <- rowMeans(do.call(cbind, interpolated_list), na.rm = TRUE)
mean_pr_df <- data.frame(Recall = common_recalls, Precision = mean_precision, Run = "Mean")
mean_pr_df <- mean_pr_df[order(mean_pr_df$Recall), ]

mean_pr_auc <- sum(diff(mean_pr_df$Recall) * 
                     (head(mean_pr_df$Precision, -1) + tail(mean_pr_df$Precision, -1)) / 2)



color_vector <- c(setNames(run_colors, all_run_labels), "Mean" = "#E69F00")

combined_pr_df <- bind_rows(pr_df_all_filtered, mean_pr_df)

Prc <- ggplot(combined_pr_df, aes(x = Recall, y = Precision, color = Run, group = Run)) +
  geom_line(size = 0.25, data = subset(combined_pr_df, Run != "Mean")) +
  geom_line(size = 0.4, data = subset(combined_pr_df, Run == "Mean")) +
  scale_y_continuous(limit = c(0.5, 1), breaks = seq(0.5, 1, 0.10), labels = scales::number_format(accuracy = 0.01)) +
  scale_color_manual(values = color_vector, name = "Run") +
  annotate("text", x = 0.05, y = 1,
           label = paste0("Mean AUC = ", round(mean_pr_auc, 3)), 
           size = 2.8, hjust = 0) +
  labs(
    title = "PR Curves",
    x = "Recall",
    y = "Precision"
  ) +
  my_theme() +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 8)
  )


roc_df_all$Run <- factor(roc_df_all$Run, levels = paste0("Run_", 1:length(unique(roc_df_all$Run))))
common_fprs <- seq(0, 1, length.out = 1000)

interpolated_list_roc <- lapply(split(roc_df_all, roc_df_all$Run), function(df) {
  approx(df$FPR, df$TPR, xout = common_fprs, rule = 2)$y
})

mean_tpr <- rowMeans(do.call(cbind, interpolated_list_roc), na.rm = TRUE)
mean_roc_df <- data.frame(FPR = common_fprs, TPR = mean_tpr)
mean_roc_df <- mean_roc_df[order(mean_roc_df$FPR), ]
mean_auc <- sum(diff(mean_roc_df$FPR) * 
                  (head(mean_roc_df$TPR, -1) + tail(mean_roc_df$TPR, -1)) / 2)
mean_roc_df$Run <- "Mean"


roc_df_all$Run <- as.character(roc_df_all$Run)
all_run_labels <- unique(roc_df_all$Run)
color_vector <- c(setNames(run_colors, all_run_labels), "Mean" = "#E69F00")


Roc <- ggplot() +
  geom_line(data = roc_df_all, mapping = aes(x = FPR, y = TPR, color = Run, group = Run), size = 0.25) +
  geom_line(data = mean_roc_df, mapping = aes(x = FPR, y = TPR, color = Run, group = Run), size = 0.4) +
  annotate("text", x = 0.05, y = 1, label = paste0("Mean AUC = ", round(mean_auc, 3)), size = 2.8, hjust = 0) +
  scale_color_manual(values = color_vector, name = "Run") +
  labs(title = "ROC Curves", x = "False Positive Rate", y = "True Positive Rate") +
  my_theme() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8)
  )



Prc <- Prc + theme(plot.margin = margin(5.5, 10, 5.5, 5.5))  
Roc <- Roc + theme(plot.margin = margin(5.5, 10, 5.5, 10))  


p <- ggarrange(Prc, Roc, nrow = 1, ncol =2, widths = c(0.8, 1))
p

ggsave("output/Fig_7c.pdf", p, width = 17, height = 5, dpi = 300, units = "cm", paper = "A4")





####Fig 7d####

library(dplyr)
library(reshape2)   
library(pheatmap)
library(magick)


library(ComplexHeatmap)
library(circlize)   
library(data.table) 

load("data/Fig_7d.rda")



col_fun <- colorRamp2(c(0, 1), c("white", "#FBA518"))


heatmap_legend_param = list(
  title_gp = gpar(fontsize = 8),
  labels_gp = gpar(fontsize = 8),
  direction = "horizontal"
)



ht_list <- lapply(1:7, function(i) {
  Heatmap(
    mat_list[[i]],
    name = "Importance",
    col = col_fun,
    show_row_names = FALSE,
    show_column_names = (i == 7),   
    column_names_rot = 0, 
    column_title_gp = gpar(fontsize=8),
    column_names_gp = gpar(fontsize=8),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_title_rot = 0,
    show_heatmap_legend = (i == 7),  
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 8),
      labels_gp = gpar(fontsize = 8),
      direction = "horizontal"
      
    ),
    row_title = paste("Run", i),
    use_raster = T,
    raster_quality =10,#8
    row_title_gp = gpar(fontsize = 8),
    row_names_gp = gpar(fontsize = 8),
  )
})


final_ht <- Reduce(`%v%`, ht_list)


pdf("output/Fig_7d.pdf", width =3.5, height =2.8)
draw(final_ht, heatmap_legend_side = "bottom")
dev.off()







####Fig 7e####

load("data/Fig_7e.rda")

library(dplyr)
library(ggplot2)
library(cowplot)

my_theme <- function() {
  theme_half_open(font_size = 8, line_size = 0.2, rel_small = 1, rel_tiny = 1, rel_large = 1) 
}

my_colors <- c(
  "abc" = "#f3f1ff", "MPK6-abc" = "#f39f16",
  "upf3" =  "#f3f1ff", "MPK6-upf3" ="#f39f16",
  "SAP5" = "#f3f1ff", "MPK6-SAP5" = "#f39f16"
)

data_plot1 <- All %>% filter(gene %in% c("SAP5", "abc"))
pvals_plot1 <- pvals %>% filter(gene %in% c("SAP5", "abc"))

plot1 <- ggplot(data_plot1, aes(x = sample, y = value, fill = sample)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.6, linewidth = 0.2, color = "black") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 0.2, color = "black") +


scale_y_continuous(expand = c(0,0), limits = c(0,5)) +
  geom_jitter(width = 0, size = 0.65, color = "black") +
  facet_wrap(~gene, nrow = 1, scales = "free_x") + 
  geom_text(data = pvals_plot1, aes(x = 1.5, y = y_pos, label = signif_label), inherit.aes = FALSE, size = 4) +
  labs(y = "Fold Change", x = "") +
  scale_fill_manual(values = my_colors) +
  theme_minimal(base_size = 12) +
  my_theme() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(plot1)

ggsave("output/Fig_7e.pdf",plot1,width =9,height =4.5,dpi = 300,paper="A4",units = "cm")




####Fig 7f####



my_theme <- function() {
  theme_half_open(font_size = 8, line_size = 0.2, rel_small = 1,rel_tiny = 1,rel_large = 1) 
  
}

load("data/Fig_7f.rda")

library(ggplot2)
library(dplyr)
custom_order <- c("SMG5", "GABARAPL1", "MFN2", "AP2A2", 
                  "PSMD5", "PRKAB1", "TRAM1", "VAMP3", 
                  "CRIPT", "TMED2")


summary_df <- combined_df %>%
  group_by(Gene) %>%
  summarise(
    Mean_Prob = mean(Probability_Class1),
    SD_Prob = sd(Probability_Class1)
  ) %>%
  ungroup()

summary_df$Gene <- factor(summary_df$Gene, levels = custom_order)
summary_df$Color <- ifelse(summary_df$Gene %in% custom_order[1:4], "orange", "steelblue")

p  <- ggplot(summary_df, aes(x = Gene, y = Mean_Prob, fill = Color)) +
  geom_hline(yintercept = 0.5, lty = "dashed", color = "gray80") +
  geom_bar(stat = "identity", width = 0.85) +
  geom_jitter(data = combined_df, aes(x = Gene, y = Probability_Class1), 
              width = 0, height = 0, size =0.18, inherit.aes = FALSE) +
  # geom_errorbar(aes(ymin = Mean_Prob - SD_Prob, ymax = Mean_Prob + SD_Prob), width = 0.2) +
  scale_fill_identity() +
  theme_minimal() +
  my_theme() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.8), breaks = seq(0, 0.8, 0.25)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Gene",
    y = "Mean Probability (± SD)"
  )

p


ggsave("output/Fig_7f.pdf",p,width =5,height =6.1,dpi = 300,units = "cm",paper="A4")




####Fig 7g####

library(dplyr)
library(ggplot2)
library(cowplot)

my_theme <- function() {
  theme_half_open(font_size = 8, line_size = 0.2, rel_small = 1, rel_tiny = 1, rel_large = 1) 
}

load("data/Fig_7g.rda")


df_all <- df_all %>% 
  filter(count >= 12)

summary_df <- df_all %>%
  group_by(gene, site) %>%
  summarise(
    pred_mean = mean(avg_pred),
    se = sd(avg_pred) / sqrt(n()),
    .groups = "drop"
  )

desired_order <- c("TRAM1")  

p <- ggplot() +
  geom_line(
    data = df_all,
    aes(x = site, y = avg_pred, group = run),
    color = "grey80",
    alpha = 0.6,
    size = 0.25
  ) +
  geom_line(
    data = summary_df,
    aes(x = site, y = pred_mean, color = gene),
    size = 0.25
  ) +
  scale_color_manual(values = c(
    "TRAM1" = "#1f77b4"
  )) +
  geom_hline(yintercept = 0.5, lty = "dashed") +
  labs(
    x = "Site",
    y = "Prediction"
  ) +
  scale_y_continuous(limits = c(0, 0.75), breaks = c(0, 0.25, 0.5, 0.75)) +
  theme(
    plot.title = element_text(size = 20, face = "bold", color = "blue"),
    plot.subtitle = element_text(size = 14, color = "darkgreen"),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  facet_wrap(~ factor(gene, levels = desired_order), scales = "free_x", ncol = 10) +
  my_theme()

p

ggsave(p,filename = "output/Fig_7g.pdf",width =6.4,height = 6,units = "cm",paper="A4")






####Fig S2b#### 
library(ggpubr)
my_theme <- function() {
  theme_half_open(font_size = 8, line_size = 0.2, rel_small = 1,rel_tiny = 1,rel_large = 1) 
  
}

load("data/Fig_S2b.rda")

GRP8$Re_peak<-GRP8$site+180
GRP8<-subset(GRP8,GRP8$Re_peak<70)
sum(GRP8$DE656)
sum(GRP8$upf1.1)
upf1_1<-ggplot() + 
  geom_segment(data=GRP8, aes(x=Re_peak,y=upf1.1,xend=Re_peak,yend=0),colour=ifelse(GRP8$upf1.1>10 ,"#eda83f","gray50"),size=0.7)+
  scale_y_continuous(name = "clone", breaks  =seq(0,20,10),limits = c(0,20),expand = c(0, 0))+
  scale_x_continuous(name = "distance to TC",  breaks =seq(0,90,30),limits = c(0,70),expand =  expand_scale(mult = c(0, 0.02), add = c(0, 0)))+
  ggtitle("upf1-5(35S:AtUPF1G851E, n=42)")+
  my_theme()



DE656<-ggplot() + 
  geom_segment(data=GRP8, aes(x=Re_peak,y=DE656,xend=Re_peak,yend=0),colour=ifelse(GRP8$DE656>10 ,"#eda83f","gray50"),size=0.7)+
  scale_y_continuous(name = "clone", breaks  =seq(0,20,10),limits = c(0,20),expand = c(0, 0))+
  scale_x_continuous(name = "distance to TC",  breaks =seq(0,90,30),limits = c(0,70),expand =  expand_scale(mult = c(0, 0.02), add = c(0, 0)))+
  ggtitle("upf1-5(35S:AtUPF1DE656AA, n=43)")+
  my_theme()


p2<-ggarrange( upf1_1, DE656 , align = "v",
               ncol = 1, nrow = 2,common.legend = TRUE)
p2


load("data/Fig_S2b_1.rda")
B2upf<-ggplot() + 
  geom_segment(data=GRP8_pare, aes(x=site,y=B2upf,xend=site,yend=0),colour=ifelse(GRP8_pare$B2upf>2000 ,"#eda83f","gray50"),size=0.7)+
  scale_y_continuous(name = "5'P end", breaks  =seq(0,3000,1500),limits = c(0,3000),expand = c(0, 0))+
  scale_x_continuous(name = "distance to TC",  breaks =seq(0,90,30),limits = c(0,70),expand =  expand_scale(mult = c(0, 0.02), add = c(0, 0)))+
  my_theme()
B2wt<-ggplot() + 
  #geom_bar(stat='identity',colour="#1F1762",fill=ifelse(GRP8$B1upf>10000 ,"#EECB27","#BEDAE5"))+
  geom_segment(data=GRP8_pare, aes(x=site,y=B2WT,xend=site,yend=0),colour=ifelse(GRP8$B2WT>2000 ,"#eda83f","gray50"),size=0.7)+
  scale_y_continuous(name = "5'P end", breaks  =seq(0,3000,1500),limits = c(0,3000),expand = c(0, 0))+
  scale_x_continuous(name = "distance to TC",  breaks =seq(0,90,30),limits = c(0,70),expand =  expand_scale(mult = c(0, 0.02), add = c(0, 0)))+
  my_theme()



p<-ggarrange(B2wt,B2upf,upf1_1, DE656 , align = "v",
             ncol = 1, nrow = 4,common.legend = TRUE)
p


ggsave("output/Fig_S2b.pdf",p,width =7.2,height = 10,dpi = 300, bg='transparent',paper="A4",units = "cm")





####Fig S5####

library(ComplexHeatmap)
library(grid)        
library(circlize)    
library(gridExtra)   


load("data/Fig_S5.rda")


colnames(data_ptc)    <- c("sample", "gene", -3:-1, 1:150)
colnames(data_normal) <- c("sample", "gene", -3:-1, 1:150)


cbind_ptc <- do.call("cbind", split(data_ptc, data_ptc$sample))
cbind_ptc <- cbind_ptc[, -156:-157]

cbind_normal <- do.call("cbind", split(data_normal, data_normal$sample))
cbind_normal <- cbind_normal[, -156:-157]


gene_total_ptc    <- rowSums(cbind_ptc[,    3:308])
gene_total_normal <- rowSums(cbind_normal[, 3:308])
combined_total    <- gene_total_ptc + gene_total_normal

threshold  <- 100
keep_genes <- which(combined_total >= threshold)

cbind_ptc    <- cbind_ptc[keep_genes, ]
cbind_normal <- cbind_normal[keep_genes, ]


WT_upf_scale_ptc <- t(apply(cbind_ptc[, 3:308], 1, scale))
Zscore_ptc       <- WT_upf_scale_ptc[, 1:153] - WT_upf_scale_ptc[, 154:306]
rownames(Zscore_ptc) <- cbind_ptc$`upf1-1.gene`

WT_upf_scale_normal <- t(apply(cbind_normal[, 3:308], 1, scale))
Zscore_normal       <- WT_upf_scale_normal[, 1:153] - WT_upf_scale_normal[, 154:306]
rownames(Zscore_normal) <- cbind_normal$`upf1-1.gene`


wt_cols  <- 3:155
upf_cols <- 156:308

wt_mat  <- as.matrix(cbind_ptc[, wt_cols])    # n_genes × 153
upf_mat <- as.matrix(cbind_ptc[, upf_cols])   # n_genes × 153

n_genes <- nrow(wt_mat)

wt_norm  <- colSums(wt_mat,  na.rm = TRUE) / n_genes
upf_norm <- colSums(upf_mat, na.rm = TRUE) / n_genes

profile_mat <- cbind(WT = wt_norm, UPF = upf_norm)   # 153 × 2


labvec_ptc <- rep('', 153)
labvec_ptc[c(3, 53, 103, 153)] <- c("0", "50", "100", "150")


col_fun <- colorRamp2(
  c(-9.5, -6, 0, 6, 9.5),
  c("#4096a6", "#93bbc2", "white", "#fabf8c", "#ed9187")
)


max_positions <- apply(Zscore_ptc, 1, which.max)
row_order_idx <- order(max_positions)


y_max <- max(profile_mat, na.rm = TRUE) * 1.1  

top_anno <- HeatmapAnnotation(
  Profile = anno_lines(
    profile_mat,                        
    ylim       = c(0, 10),
    smooth     = TRUE,
    add_points = FALSE,
    height     = unit(2, "cm"),
    gp         = gpar(
      col = c('#f5bc1d', '#0b318f'),  
      lwd = 2
    ),
    axis_param = list(side = "left")
  ),
  Labels_Top = anno_text(
    labvec_ptc,
    rot  = 0,
    just = "center",
    gp   = gpar(fontsize = 8)
  ),
  show_annotation_name = FALSE,
  border = FALSE
)

bottom_anno <- HeatmapAnnotation(
  Ticks = anno_text(
    labvec_ptc,
    rot  = 0,
    just = "center",
    gp   = gpar(fontsize = 8)
  ),
  show_annotation_name = FALSE,
  border = FALSE
)


ht <- Heatmap(
  Zscore_ptc,
  name             = "Δ Z-score",
  col              = col_fun,
  
  cluster_columns  = FALSE,
  cluster_rows     = FALSE,
  row_order        = row_order_idx,
  
  show_column_names = FALSE,
  row_names_gp      = gpar(fontsize = 10),
  row_title_gp      = gpar(fontsize = 12),
  
  column_title      = "Distance downstream of PTC (nt)",
  column_title_side = "bottom",
  column_title_gp   = gpar(fontsize = 12, fontface = "plain"),
  
  row_labels        = gsub("_0$", "", rownames(Zscore_ptc)),
  row_names_side    = "left",
  row_names_max_width = unit(10, "cm"),
  
  top_annotation    = top_anno,
  bottom_annotation = bottom_anno,
  
  column_split = NULL,
  
  heatmap_legend_param = list(
    title            = "ΔZ-score",
    legend_direction = "vertical",
    title_gp         = gpar(fontsize = 10, fontface = "bold"),
    labels_gp        = gpar(fontsize = 10)
  )
)


draw(ht, heatmap_legend_side = "right")


pdf("output/Fig_S5.pdf", width = 7, height = 10)
draw(ht, heatmap_legend_side = "right")
dev.off()







####Fig S6####

library(ggplot2)
library(ggrastr)
library(cowplot)
library(ggpubr)
library(dplyr)

my_theme <- function() {
  theme_half_open(font_size = 8, line_size = 0.2, rel_small = 1, rel_tiny = 1, rel_large = 1) 
}

load("data/Fig_S6.rda")

filtered_data_upf <- all_data[!is.na(all_data$upf_logFC) & !is.na(all_data$upf_padj), ]
filtered_data_upf$logP <- -log10(filtered_data_upf$upf_padj)

filtered_data_smg <- all_data[!is.na(all_data$smg_logFC) & !is.na(all_data$smg_padj), ]
filtered_data_smg$logP <- -log10(filtered_data_smg$smg_padj)
class_colors <- c(
  
  "all"             = "gray85",
  
  "long_no_peak"    = "#844B96",
  
  "long_peak"       = "#E69F00",
  
  "intron_no_peak"  ="#7570B3",
  
  "intron_peak"     ="#E7298A",
  
  "shortI_no_peak"  = "#1B9E77",
  
  "shortI_peak"     = "#D95F02",
  
  "shortII_no_peak" =  "#0072B2",
  
  "shortII_peak"    =  "#D55E00"
  
)

calc_stats <- function(df, prefix) {
  fc_col <- paste0(prefix, "_logFC")
  padj_col <- paste0(prefix, "_padj")
  classes <- c("all", "long_no_peak", "long_peak", 
               "intron_no_peak", "intron_peak",
               "shortI_no_peak", "shortI_peak",
               "shortII_no_peak", "shortII_peak")
  
  stats_list <- list()
  for (cl in classes) {
    n_tot <- sum(df$class == cl)
    n_up  <- sum(df$class == cl & df[[padj_col]] < 0.01 & df[[fc_col]] > 1, na.rm = TRUE)
    n_dn  <- sum(df$class == cl & df[[padj_col]] < 0.01 & df[[fc_col]] < -1, na.rm = TRUE)
    r_up  <- ifelse(n_tot > 0, round((n_up / n_tot) * 100, 2), 0)
    r_dn  <- ifelse(n_tot > 0, round((n_dn / n_tot) * 100, 2), 0)
    stats_list[[cl]] <- list(total = n_tot, up = n_up, dn = n_dn, r_up = r_up, r_dn = r_dn)
  }
  return(stats_list)
}

stats_upf <- calc_stats(filtered_data_upf, "upf")
stats_smg <- calc_stats(filtered_data_smg, "smg")

thresholds <- seq(0.1, 2.0, by = 0.1) 
classes <- c("long_no_peak", "long_peak", "intron_no_peak", "intron_peak", "shortI_no_peak", "shortI_peak", "shortII_no_peak", "shortII_peak")

legend_labels <- sapply(classes, function(cl) {
  n_upf <- stats_upf[[cl]]$total
  n_smg <- stats_smg[[cl]]$total
  paste0(cl, " (UPF: N=", n_upf, " | SMG: N=", n_smg, ")")
})
names(legend_labels) <- classes

line_data_list <- list()
for (cl in classes) {
  n_tot_upf <- sum(filtered_data_upf$class == cl)
  n_tot_smg <- sum(filtered_data_smg$class == cl)
  
  for (th in thresholds) {
    n_up_upf <- sum(filtered_data_upf$class == cl & filtered_data_upf$upf_padj < 0.01 & filtered_data_upf$upf_logFC > th, na.rm = TRUE)
    n_dn_upf <- sum(filtered_data_upf$class == cl & filtered_data_upf$upf_padj < 0.01 & filtered_data_upf$upf_logFC < -th, na.rm = TRUE)
    n_up_smg <- sum(filtered_data_smg$class == cl & filtered_data_smg$smg_padj < 0.01 & filtered_data_smg$smg_logFC > th, na.rm = TRUE)
    n_dn_smg <- sum(filtered_data_smg$class == cl & filtered_data_smg$smg_padj < 0.01 & filtered_data_smg$smg_logFC < -th, na.rm = TRUE)
    
    line_data_list[[length(line_data_list) + 1]] <- data.frame(dataset = "UPF", class = cl, direction = "Up-regulated (logFC > Threshold)", threshold = th, ratio = ifelse(n_tot_upf > 0, (n_up_upf / n_tot_upf) * 100, 0))
    line_data_list[[length(line_data_list) + 1]] <- data.frame(dataset = "UPF", class = cl, direction = "Down-regulated (logFC < -Threshold)", threshold = th, ratio = ifelse(n_tot_upf > 0, (n_dn_upf / n_tot_upf) * 100, 0))
    line_data_list[[length(line_data_list) + 1]] <- data.frame(dataset = "SMG", class = cl, direction = "Up-regulated (logFC > Threshold)", threshold = th, ratio = ifelse(n_tot_smg > 0, (n_up_smg / n_tot_smg) * 100, 0))
    line_data_list[[length(line_data_list) + 1]] <- data.frame(dataset = "SMG", class = cl, direction = "Down-regulated (logFC < -Threshold)", threshold = th, ratio = ifelse(n_tot_smg > 0, (n_dn_smg / n_tot_smg) * 100, 0))
  }
}

line_data <- do.call(rbind, line_data_list)
line_data$class <- factor(line_data$class, levels = classes)

line_data$dataset <- factor(line_data$dataset, levels = c("SMG", "UPF"))

line_data$category <- NA
line_data$category[grepl("long", line_data$class)]    <- "Long"
line_data$category[grepl("intron", line_data$class)]  <- "Intron"
line_data$category[grepl("shortI_", line_data$class)] <- "Short I"  
line_data$category[grepl("shortII", line_data$class)] <- "Short II"

max_ratio <- max(line_data$ratio, na.rm = TRUE)
global_y_limit <- ceiling(max_ratio / 5) * 5 

plot_category_lines <- function(cat_name) {
  df_sub <- line_data[line_data$category == cat_name, ]
  p <- ggplot(df_sub, aes(x = threshold, y = ratio, color = class, group = class)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = class_colors, labels = legend_labels) +
    facet_grid(dataset ~ direction) +  
    scale_x_continuous(breaks = seq(0, 2.0, by = 0.5), limits = c(0, 2.0)) + 
    scale_y_continuous(limits = c(0, global_y_limit)) + 
    labs(title = paste0(cat_name), x = "Absolute log2 Fold Change Threshold", y = "Ratio of genes pass Threshold (%)", color = NULL) +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "gray90", color = "black"), 
      strip.text = element_text(face = "bold"), 
      legend.position = "bottom", 
      legend.direction = "vertical", 
      legend.text = element_text(size = 8),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  return(p)
}

p_long    <- plot_category_lines("Long")
p_intron  <- plot_category_lines("Intron")
p_shortI  <- plot_category_lines("Short I")
p_shortII <- plot_category_lines("Short II")

line_grid <- ggarrange(p_long, p_intron, p_shortI, p_shortII, ncol = 2, nrow = 2, align = "hv", labels = c("a", "b", "c", "d"))

ggsave("output/Fig_S6.pdf", line_grid, width = 14, height = 10, dpi = 300)






####Fig S7####

load("data/Fig_S7.rda")


library(tidyverse)
library(cowplot)

my_theme <- function() {
  theme_half_open(font_size = 8, line_size = 0.2, rel_small = 1, rel_tiny = 1, rel_large = 1)
}

final_plot <- ggplot(all_data, aes(x = Condition, y = GO_term)) +
  geom_tile(aes(fill = neg_log10p), alpha = 0.1) + 
  geom_point(aes(size = fold, color = neg_log10p)) +
  scale_size_continuous(
    range = c(0.4, 2.2), 
    name = "Fold Enrichment",
    breaks = c(3, 6,9) 
  ) + 
  scale_color_viridis_c(option = "magma", name = "-log10(p)") +
  scale_fill_viridis_c(option = "magma", guide = "none") +
  facet_grid(GO_class ~ Group, scales = "free_y", space = "free_y") +
  theme_bw(base_size = 9) +
  theme(
    axis.text.y = element_text(
      size = 6.5,
      color = "black",
      lineheight = 0.8
    ),
    plot.margin = margin(t = 1, r = 2, b = 1, l = 1.5, unit = "cm"),
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(), 
    axis.text.x = element_text(size = 7),
    strip.text.y = element_text(angle = -90, size = 7, face = "bold"),
    panel.spacing.y = unit(0.05, "lines"),
    legend.position = "bottom"
  )

ggsave(final_plot, file = "output/Fig_S7.pdf", height = 28, width = 20, units = "cm", paper = "A4")









####Fig S8a####


load("data/Fig_S8a.rda")

colnames(total1) <- c("length","GC")
total1$GC <- total1$GC * 0.01

S1 <- c(total1$GC)
v1 <- c(total1$length)

linear_model <- lm(v1 ~ S1, data=total1)
summary(linear_model)

library(ggplot2)
library(cowplot)
library(ggrastr)

my_theme <- function() {
  theme_half_open(
    font_size = 8,
    line_size = 0.2,
    rel_small = 1,
    rel_tiny = 1,
    rel_large = 1
  )
}

eq <- paste0(
  "y = ", round(coef(linear_model)[2], 4),
  "x + ", round(coef(linear_model)[1], 2),
  "\nR² = ", round(summary(linear_model)$r.squared, 4)
)

plotB3 <- ggplot() +
  

  geom_abline(
    slope = linear_model$coefficients[2],
    intercept = linear_model$coefficients[1],
    color = "red", 
    linewidth = 0.5
  ) +
  
  scale_y_continuous(
    breaks = seq(0, 550, 250),
    limits = c(0, 550),
    name = "3'UTR Length (nt)",
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  scale_x_continuous(
    breaks = seq(0.25, 0.4, 0.05),
    limits = c(0.25, 0.4),
    name = "GC Content",
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  geom_point_rast(
    aes(x = S1, y = v1),
    shape = 21,          
    fill = "#9895bd",  
    color = "#6a66a0",  
    alpha = 0.95,
    size = 1.8,
    stroke = 0.25
  ) +
  
  my_theme()

ggsave(
  "output/Fig_S8a.pdf",
  plotB3,
  width =6.5,
  height = 7.5,
  dpi = 300,
  units = "cm",
  paper = "A4"
)





####Fig S8b####

library(dplyr)
library(ggplot2)
library(cowplot)



load("data/Fig_S8b.rda")


coef_data <- binned_data %>%
  group_by(length_bin) %>%
  summarize(
    model = list(lm(Detection_Rate_Pct ~ Avg_GC)),
    intercept = coef(model[[1]])[1],
    slope = coef(model[[1]])[2],
    r_squared = summary(model[[1]])$r.squared,
    eq_label = sprintf(
      "y = %.2fx %s %.2f\nR² = %.3f", 
      slope, 
      ifelse(intercept >= 0, "+", "-"), 
      abs(intercept), 
      r_squared
    ),
    .groups = "drop"
  )
for(i in 1:3) {
  cat("\n--- Results for", levels(binned_data$length_bin)[i], "---\n")
  print(summary(lm(Detection_Rate_Pct ~ Avg_GC, data = filter(binned_data, length_bin == levels(binned_data$length_bin)[i])))$coefficients)
}


p_stratified <- ggplot(binned_data, aes(x = Avg_GC, y = Detection_Rate_Pct)) +
  geom_point(alpha = 0.8, size = 0.75, shape = 16, color = "#208091") +
  geom_abline(data = coef_data, aes(intercept = intercept, slope = slope), 
              color = "#D65A47", linewidth = 0.5) +
  geom_text(data = coef_data, aes(x = -Inf, y = Inf, label = eq_label), 
            hjust = -0.1, vjust = 1.3, size = 3, color = "#208091") +
  facet_wrap(~length_bin, ncol = 1, axes = "all_x", axis.labels = "all_x") +  
  scale_x_continuous(name = "Average GC Ratio (in bins of 50 transcripts)") +
  scale_y_continuous(name = "UIF detection rate (count/100 transcripts)") +
  theme_half_open(font_size = 9) +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p_stratified)

ggsave("output/Fig_S8b.pdf", p_stratified, width = 6, height = 8, units = "cm")









####Fig S9####

load("data/Fig_S9.rda")

library(ComplexHeatmap)
library(circlize)  
library(data.table) 


col_fun <- colorRamp2(c(0, 1), c("white", "#FBA518"))

heatmap_legend_param = list(
  title_gp = gpar(fontsize = 8),
  labels_gp = gpar(fontsize = 8),
  direction = "vertical"
)


ht_list <- lapply(1:5, function(i) {
  Heatmap(
    mat_list[[i]],
    name = "Importance",
    col = col_fun,
    show_row_names = FALSE,
    show_column_names = (i == 5),   
    column_names_rot = 0, 
    column_title_gp = gpar(fontsize=8),
    column_names_gp = gpar(fontsize=8),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_title_rot = 0,
    show_heatmap_legend = (i == 5),  
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 8),
      labels_gp = gpar(fontsize = 8),
      direction = "vertical"
      
    ),
    row_title = paste("Run", i),
    use_raster = T,
    raster_quality =10,#8
    row_title_gp = gpar(fontsize = 8),
    row_names_gp = gpar(fontsize = 8),
  )
})


final_ht <- Reduce(`%v%`, ht_list)

pdf("output/Fig_S9.pdf", width = 6.5, height =1.8)
draw(final_ht, heatmap_legend_side = "right")
dev.off()


