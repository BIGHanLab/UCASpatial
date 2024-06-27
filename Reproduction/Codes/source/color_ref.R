library(RColorBrewer)
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)

# from Yawei
mycolors <- brewer.pal(11, "Paired")
pt8_colors <- c(mycolors[3],'#38b000','#008000',mycolors[c(1,7,8,11,9,10,5,6)])

pt2_colors <- c(mycolors[4],mycolors[1:2],mycolors[5],'#ff4d6d',mycolors[6])
pt2_colors <- pt2_colors[c(6,1,2,3,4,5)]

wtpt3_colors <- c(pt8_colors[1],pt8_colors[8:9],pt8_colors[5],pt8_colors[7],'#247ba0',mycolors[5],'#ff4d6d',mycolors[6])

# from Guanghao
color2 <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", 
            "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")

# from Zurui
# SpatialPlot(DSS_14h_spa_result,features = c('Smooth.muscle.cell'),pt.size.factor = 5,image.alpha = 0,
#             crop = T,min.cutoff = 0,stroke = NA)+
#   scale_fill_gradientn(colours = c("black","#010107","#3A0F6F","#631980","#882781","#B3367A","#DC4868","#F66C5C","#FE9F6D","#FED093","#FCFABC"))+
#   scale_color_manual(values = "white")+theme(panel.background = element_rect(fill = "black"),legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
#   guides(fill = guide_colourbar(title = " ",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 23,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 1,
#                                 frame.colour = "black",frame.linewidth = 0.5,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=10)))+labs(title = "Smooth Muscle Cell")


# Xuyin
cols_p8 <- c("#BEBADA","#D9D9D9","#CCEBC5","#B15928","#FF7F00","#1F78B4","#BC80BD","#33A02C","#B2DF8A","#E31A1C","#FB8072")
sc_cols_p8 <- c("#081D58","#A6CEE3" ,"#1F78B4","#33A02C","#FC4E2A","#E31A1C","#BD0026","#B15928","#1D91C0","#CAB2D6","#6A3D9A","#B2DF8A")

cols_p3 <- c(pt8_colors[1],pt8_colors[8:9],pt8_colors[5],pt8_colors[7],'#247ba0',mycolors[5],'#ff4d6d',mycolors[6])
sc_cols_p3 <- c("#081D58","#A6CEE3" ,"#1F78B4","#33A02C","#FC4E2A","#BD0026","#B15928","#1D91C0","#CAB2D6","#6A3D9A","#B2DF8A")

cols_p2 <- c(mycolors[4],mycolors[1:2],mycolors[5],'#ff4d6d',mycolors[6])
cols_p2 <- cols_p3[c(6,1,2,3,4,5)]
sc_cols_p2 <-  c("#081D58","#A6CEE3" ,"#1F78B4","#33A02C","#BD0026","#B15928","#1D91C0","#CAB2D6","#6A3D9A","#B2DF8A")

cols_p8_kmeans <- c(cols_p8[1],cols_p8[8],sc_cols_p8[12],sc_cols_p8[2],sc_cols_p8[3],
                    sc_cols_p8[11],sc_cols_p8[1],cols_p8[5])


cols_p8_layer <- c("#00AFBB", "#FC4E07", "#E7B800")

cols_p8_stromal_old <- brewer.pal(10, "Paired")
cols_p8_stromal <- c("#A6CEE3","#1F78B4","#E31A1C","#FF7F00","#6A3D9A")
cols_p8_stromal_newident <- c("#A6CEE3","#FF7F00","#E31A1C","#1F78B4","#6A3D9A")

adapt_cols <- brewer.pal(12,"Paired")

# install.packages("ggsci")
library("ggsci")

library(scales)
#show_col(brewer.pal(12,"Paired"))




