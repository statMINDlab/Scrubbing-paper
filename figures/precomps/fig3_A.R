source("code/0_SharedCode.R")

library(magick)
library(cowplot)

bases <- c("P36", "FIX_GSR")
p <- setNames(vector("list", length(bases)), bases)
for (base in bases) {
  pdir <- file.path(dir_github, "plots/QCFC", base)
  p[[base]] <- cowplot::plot_grid(
   ggdraw()+draw_image(image_read(file.path(pdir, "QCFC_Naive_None.png"))),
   ggdraw()+draw_image(image_read(file.path(pdir, "QCFC_Naive_FD5.png"))),
   ggdraw()+draw_image(image_read(file.path(pdir, "QCFC_Naive_FD2.png"))),
   ggdraw()+draw_image(image_read(file.path(pdir, "QCFC_Naive_FD2+.png"))),

   NULL, NULL, NULL, NULL,

   ggdraw()+draw_image(image_read(file.path(pdir, "QCFC_Between_None.png"))),
   ggdraw()+draw_image(image_read(file.path(pdir, "QCFC_Between_FD5.png"))),
   ggdraw()+draw_image(image_read(file.path(pdir, "QCFC_Between_FD2.png"))),
   ggdraw()+draw_image(image_read(file.path(pdir, "QCFC_Naive_FD2+.png"))),

   ggdraw()+draw_image(image_read(file.path(pdir, "QCFC_Within_None.png"))),
   ggdraw()+draw_image(image_read(file.path(pdir, "QCFC_Within_FD5.png"))),
   ggdraw()+draw_image(image_read(file.path(pdir, "QCFC_Within_FD2.png"))),
   ggdraw()+draw_image(image_read(file.path(pdir, "QCFC_Naive_FD2+.png"))),

    nrow=4, rel_heights = c(6,1,6,6)
  )
}

png(file.path(dir_github, "figures/precomps/Fig3A.png"), width=2400, height=1900, res=200)
plot(p$P36)
dev.off()

png(file.path(dir_github, "figures/precomps/Fig3A_App_FIX_GSR.png"), width=2400, height=1900, res=200)
plot(p$FIX_GSR)
dev.off()
