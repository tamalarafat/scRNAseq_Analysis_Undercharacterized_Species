df = data.frame("peak" = c(9336, 20), "nopeak" = c(18401, 16), row.names = c("chip", "markers")) # 20 out of 23

fisher.test(df) # p-value = 0.007578

df = data.frame("peak" = c(9290, 66), "nopeak" = c(18373, 44), row.names = c("chip", "markers")) # 66 out of 87

fisher.test(df) # 2.108e-08

df = data.frame("peak" = c(9333, 23), "nopeak" = c(18410, 7), row.names = c("chip", "markers")) # 23 out of 26

fisher.test(df) # 1.793e-06
