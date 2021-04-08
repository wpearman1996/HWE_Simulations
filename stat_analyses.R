
plot_grid(low_struc_nucdist +
          stat_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x),
                                jittered_points = TRUE,position="points_sina",
                                alpha = 0.7, scale = 0.9,bandwidth = 0.00023)+xlim(0,0.004),
          low_rand_struc_nucdist +
          stat_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x),
                                jittered_points = TRUE,position="points_sina",
                                alpha = 0.7, scale = 0.9,bandwidth = 0.00023)+xlim(0,0.004),
          ncol = 1,labels=c("A","B"),rel_widths = c(1,1))

t.test(low_rand_struc_mats$struc_mats[low_rand_struc_mats$filter=="No Filter"],
       low_struc_mats$struc_mats[low_struc_mats$filter=="No Filter"])
t.test(low_rand_struc_mats$struc_mats[low_rand_struc_mats$filter=="HWE Out Any"],
       low_struc_mats$struc_mats[low_struc_mats$filter=="HWE Out Any"])
t.test(low_rand_struc_mats$struc_mats[low_rand_struc_mats$filter=="HWE Out All"],
       low_struc_mats$struc_mats[low_struc_mats$filter=="HWE Out All"])
t.test(low_rand_struc_mats$struc_mats[low_rand_struc_mats$filter=="HWE Out Across"],
       low_struc_mats$struc_mats[low_struc_mats$filter=="HWE Out Across"])

plot_grid(extreme_struc_nucdist + 
         stat_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x),
                                jittered_points = TRUE,position="points_sina",
                                alpha = 0.7, scale = 0.9,bandwidth = 0.00156) + xlim(0,0.05),
          extreme_rand_struc_nucdist +
          stat_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x),
                                jittered_points = TRUE,position="points_sina",
                                alpha = 0.7, scale = 0.9,bandwidth = 0.00156)+ xlim(0,0.05),
          ncol = 1,labels=c("A","B"),rel_widths = c(1,1))
t.test(extreme_rand_struc_mats$struc_mats[extreme_rand_struc_mats$filter=="No Filter"],
       extreme_struc_mats$struc_mats[extreme_struc_mats$filter=="No Filter"])
t.test(extreme_rand_struc_mats$struc_mats[extreme_rand_struc_mats$filter=="HWE Out Any"],
       extreme_struc_mats$struc_mats[extreme_struc_mats$filter=="HWE Out Any"])
t.test(extreme_rand_struc_mats$struc_mats[extreme_rand_struc_mats$filter=="HWE Out All"],
       extreme_struc_mats$struc_mats[extreme_struc_mats$filter=="HWE Out All"])
t.test(extreme_rand_struc_mats$struc_mats[extreme_rand_struc_mats$filter=="HWE Out Across"],
       extreme_struc_mats$struc_mats[extreme_struc_mats$filter=="HWE Out Across"])
