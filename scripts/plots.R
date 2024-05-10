suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(QCEWAS)
})

# Define command line arguments
parser <- argparse::ArgumentParser(description="Script for plotting ewas results")
parser$add_argument('--input-file', '-i',
                    required=TRUE,
                    help="Path to annotated ewas results data")
parser$add_argument('--out-dir',
                    required=TRUE,
                    help="Path to output directory")
parser$add_argument('--stratified',
                    choices=c("yes", "no"), 
                    default="no",
                    help="Results from a stratified analysis: yes or no")
parser$add_argument("--assoc", 
                    required=TRUE,
                    type="character", 
                    nargs=1, 
                    help="Association variable EWAS was performed with.")


# parse arguments
args <- parser$parse_args()
results <- args$input_file
out_dir <- args$out_dir
stratified <- args$stratified
assoc <- args$assoc

# Read in annotated EWAS results
ewas <- fread(results)

# Wrangle data for plotting
if (stratified == "yes"){
  ewas <- ewas %>% 
    dplyr::select(MarkerName, CpG_chrm, CpG_beg, "P-value") %>% 
    rename("Pvalue" = "P-value") %>% 
    mutate("CHR" = as.numeric(sub("chr", "", CpG_chrm)),
           "MAPINFO" = CpG_beg)
} else{
  ewas <- ewas %>% 
    dplyr::select(cpgid, CpG_chrm, CpG_beg, bacon.pval)%>% 
    rename(Pvalue = bacon.pval) %>% 
    mutate(CHR = as.numeric(sub("chr", "", CpG_chrm)),
           MAPINFO = CpG_beg)
}

#################################################################
#                  MANHATTAN & QQ PLOT                          #
#################################################################


# Calculate the cumulative position of each chromosome
chr.pos <- ewas %>% 
  group_by(CHR) %>% 
  summarize(chr_len=max(as.numeric(MAPINFO))) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>% 
  dplyr::select(-chr_len) 


# Add cumulative position data to results
ewas <- left_join(ewas, chr.pos, by = "CHR") %>% 
  arrange(CHR, as.numeric(MAPINFO)) %>% 
  mutate(POS = as.numeric(MAPINFO)+tot)


x_axis <- ewas %>% group_by(CHR) %>% 
  summarise(center = (max(POS) + min(POS))/2) 

ewas %>% 
  filter(-log10(Pvalue)>1) %>% 
  ggplot(aes(x=POS, y=-log10(Pvalue))) +
  geom_point(aes(color=as.factor(CHR)), alpha = 0.8, size = 2) +
  scale_color_manual(values= rep(c("steelblue1", "steelblue4"), 
                                 length.out = 23)) + 
  scale_x_continuous(label = x_axis$CHR, breaks = x_axis$center,
                     guide = guide_axis(check.overlap = T), expand=c(0,0)) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(1,28)) +
  geom_hline(yintercept = -log10(0.05/573388),
             linetype = 'solid',
             color = "red",
             linewidth = 0.5) +
  theme_bw(base_size = 16) +
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(y= expression(-log[10]("P-value")),
       x= "Chromosome") -> manh.plot


StatQQplot <- ggproto("StatQQplot", Stat,
                      default_aes = aes(y = stat(observed), x = stat(expected)),
                      
                      required_aes = c("observed"),
                      
                      compute_group = function(data, scales, dparams = list(),
                                               na.rm = FALSE) {
                        
                        observed <- data$observed#[!is.na(data$x)]
                        N <- length(observed)
                        
                        ## expected
                        expected <- sort(-log10((1:N)/N-1/(2*N)))
                        observed <- sort(-log10(observed))
                        data.frame(observed, expected)
                        
                      }
)

stat_qqplot <- function(mapping = NULL, data = NULL, geom = "point",
                        position = "identity", na.rm = FALSE, show.legend = NA, 
                        inherit.aes = TRUE, ...) {
  layer(
    stat = StatQQplot, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

lambda <- QCEWAS::P_lambda(ewas$Pvalue)
lambda_label <- paste0("lambda==", round(lambda, digits = 2))

ggplot(ewas, aes(observed= Pvalue)) +
  stat_qqplot() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  theme_bw(base_size = 16) +
  annotate("text", x=1, y=25, label = lambda_label, parse=T, size = 6)+
  labs(y= expression(Observed ~ -log[10]("P-value")),
      x = expression(Expected ~ -log[10]("P-value"))) -> qq.plot

left.panel <- plot_grid(NULL, qq.plot, labels= c("B", ""), label_size = 22, 
                        ncol=1, rel_heights = c(1,5))
full.plot <- plot_grid(manh.plot, left.panel, labels = c("A", ""), ncol = 2,
          rel_widths = c(2.5,1), label_size = 22)

filename <- paste0(out_dir, assoc, "_ewas_manhattan_qq_plots.jpg")
ggsave(filename,
      plot = full.plot,
      width = 30,
      height = 12,
      units = "cm")
