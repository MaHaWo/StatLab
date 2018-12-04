# coding: utf-8
import pandas as pd 
import numpy as np
data = np.fromfile('data.bin', np.float64)
data = pd.DataFrame({"value":data})
dds = pd.concat([data.describe(), pd.DataFrame({'skewness':data.skew()}).T, pd.DataFrame({'excess_kurtosis':data.kurtosis()}).T])
ddsv = dds.values
ddsv = [d[0] for d in ddsv]
dds.to_csv("description.csv")
np.array(ddsv).tofile("description.bin")
