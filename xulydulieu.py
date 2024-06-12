import pandas as pd
import numpy as np
df = pd.read_csv("D:/Na/ReGeneL/gene_diff_data.csv")
sid =list(df["SID"])
gse = pd.read_csv("D:/Na/ReGeneL/GSE66099.csv")
list1= ["SID"] + sid 
gen = gse[list1]
print(gen)
#gen.to_csv("D:/Na/ReGeneL/dulieudaxuly_nolabel.csv")




