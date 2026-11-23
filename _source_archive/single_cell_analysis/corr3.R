(base) [uqsnegi@bunya1 corr_log_results]$ vim corr3.R 

  }

  # Prepare correlation matrix
  x2 = dcast(rg, p1 ~ p2, value.var = "rg",fun.aggregate=mean)
  mat2 = as.matrix(x2[, 2:ncol(x2)])
  rownames(mat2) = x2$p1
  mat2[mat2 > 1] = 1
  mat2[mat2 < -1] = -1
  mat2[is.na(mat2)] = 0

  # Prepare FDR matrix
  x2 = dcast(rg, p1 ~ p2, value.var = "q", fun.aggregate=mean)
  qmat2 = as.matrix(x2[, 2:ncol(x2)])
  rownames(qmat2) = x2$p1
  qmat2[is.na(qmat2)] = 1

  if (nrow(mat2) == ncol(mat2)) {
    diag(mat2) = 1
    diag(qmat2) = -1
  }

  trait1_idx = match(trait1, rownames(mat2))
  trait2_idx = match(trait2, colnames(mat2))
  mat2 = mat2[trait1_idx, trait2_idx]
  qmat2 = qmat2[trait1_idx, trait2_idx]

  # Calculate the number of significant correlations per trait
  sig_counts = apply(qmat2, 1, function(x) sum(x < 0.05))

  # Set up the layout to include the barplot
  layout(matrix(c(1, 2), nrow = 2), heights = c(3, 1))

  # Plot correlation matrix with significance
  corrplot(mat2, method = "square", order = order,
           p.mat = qmat2, sig.level = 0.05, insig = "label_sig",
           pch = "*", pch.cex = 1.5, full_col=FALSE,
           na.label = "square", na.label.col = "grey30")

  # Add legend for square size
  par(new = TRUE)
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(-1, 1), ylim = c(-1, 1))
  legend("bottomleft", legend = c("FDR < 0.05", "FDR < 0.1", "FDR >= 0.1"),
         pch = 22, pt.cex = c(2, 1.5, 1), col = "black",
         title = "Square Size", cex = 1.2, bty = "n")

  # Plot the barplot
  par(mar = c(5, 4, 0, 2))  # Adjust margins to fit the barplot
  barplot(sig_counts, main = "Significant Correlations per Trait",
          xlab = "Trait", ylab = "Number of Significant Correlations",
          las = 2, col = "skyblue")
}

# Example usage
x11()
pdf("corrplot2_rg_disease_with_barplot.pdf", width = 16, height = 12, family = "Helvetica")
  corrplot_nsquare_with_barplot(rg, traitlist1$TRAIT, traitlist2$TRAIT, landscape=TRUE)
dev.off()

png("corrplot2_rg_disease_with_barplot.png", width = 16, height = 12, units = 'in', res = 300, family = "Helvetica")
  corrplot_nsquare_with_barplot(rg, traitlist1$TRAIT, traitlist2$TRAIT, landscape=TRUE)
dev.off()
