library(WindR)
w.start()
universe <- c('600507.SH', '600103.SH', '600401.SH', '601088.SH', '601006.SH', '600395.SH', '601369.SH', '601398.SH', '000568.SH', '601228.SH', '600397.SH', '601939.SH', '601988.SH', '601566.SH', '000157.SZ', '002507.SZ', '600403.SH', '600377.SH', '600348.SH', '601328.SH', '601717.SH', '600028.SH', '603366.SH', '600139.SH')
ptd <- data.frame()
n <- length(universe)
for(i in 1 : n){
  ptd[universe[i]] <- w.wsd(universe[i], 'close', '2016-01-01', '2016-12-31')
}
for(i in 1 : n){
  for(j in 1 : i){
    jpeg(file = paste(universe[i], '.jpeg', sep = paste0('&', universe[j])))
    plot(ptd[,i], ptd[,j], type = 'b', xlab = universe[i], ylab = universe[j])
    dev.off()
  }
}
