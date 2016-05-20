

create.plot = function(title, file, pdf = T, corr = T ){
  res = read.table(file, header = T)
  
  if(pdf){
    pdf(paste("loglog_plots_", algo, "_nnods_", nnodes, ".pdf", sep = ""))
  }
  
  cols = c(rep("green4", floor(nrow(res)/2)), rep("red", (nrow(res) - floor(nrow(res)/2))))
  plot(res$est.shd, res$stars.subs,
       pch = 20, cex = 0.7, col = cols,
       main = title,
       xlab = "shd to true DAG",
       ylab = "SSCI")
  if(corr){
    fit = lm(res$stars.sub~res$est.shd)
    abline(fit, lty = 2)
    legend("topleft", bty="n", legend=paste("R² =", format(summary(fit)$adj.r.squared, digits=3)))
    
  }
  
  plot(res$est.eqclass.shd, res$stars.subs,
       pch = 20, cex = 0.7, col = cols,
       main = title,
       xlab = "shd to true CPDAG",
       ylab = "SSCI")
  if(corr){
    fit = lm(res$stars.subs~res$est.eqclass.shd)
    abline(fit, lty = 2)
    legend("topleft", bty="n", legend=paste("R² =", format(summary(fit)$adj.r.squared, digits=3)))
    
  }
  
#   plot(res$est.shd, res$stars.resamp,
#        pch = 20, cex = 0.7, col = cols,
#        main = title,
#        xlab = "shd to true DAG",
#        ylab = "SSCI")
#   if(corr){
#     fit = lm(res$stars.resamp~res$est.shd)
#     abline(fit, lty = 2)
#     legend("topleft", bty="n", legend=paste("R² =", format(summary(fit)$adj.r.squared, digits=3)))
#     
#   }
#   
#   
#   plot(res$est.eqclass.shd, res$stars.resamp,
#        pch = 20, cex = 0.7, col = cols,
#        main = title,
#        xlab = "shd to true CPDAG",
#        ylab = "SSCI")
#   if(corr){
#     fit = lm(res$stars.resamp~res$est.eqclass.shd)
#     abline(fit, lty = 2)
#     legend("topleft", bty="n", legend=paste("R² =", format(summary(fit)$adj.r.squared, digits=3)))
#   }
#   
  
  if(pdf){
    dev.off()
  }
  
}

