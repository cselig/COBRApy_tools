# creates a dataframe with both experimental and model (FBA) fitness data

import pandas
import cobra.io.sbml
from cobra.flux_analysis import single_deletion
import matplotlib.pyplot as plt
import math

# Experimental data
s2 = pandas.read_excel('../experimental_data/dunham_S2.xlsx', 'Sheet1')

# Create model, record wt fitness
model = cobra.io.sbml.create_cobra_model_from_sbml_file('../SBML/Yeast_YPD_Medium.xml', 
	print_time=True)
biomassRxn = model.reactions[1574] # biomass reaction
model.optimize()
wtFitness = model.solution.x_dict[biomassRxn.id]

# perform deletions
growth_rates, statuses = single_deletion(model, element_list = model.genes[1:905])
sdFrame = pandas.DataFrame.from_dict(growth_rates, 'index')
sdFrame.columns = ['Growth_Rate']

sdFrame['wt_Fitness'] = pandas.Series()
sdFrame['% change'] = pandas.Series()

# calculate % changes for FBA data
for index, row in sdFrame.iterrows():
	row['wt_Fitness'] = wtFitness
	row['% change'] = (row['Growth_Rate'] - wtFitness) / wtFitness

sdFrame['M_Fitness'] = pandas.Series()

# MM1N-Phosphate
s2 = s2.set_index(['Name'])

# populate experimental fitness data
for index, row in s2.iterrows():
	sdFrame.ix[index, 'M_Fitness'] = s2['MM1N-Sulfate'][str(index)]

# get rid of genes where we don't have data from both FBA and experiment
for index, row in sdFrame.iterrows():
	if math.isnan(sdFrame['M_Fitness'][index]) or math.isnan(sdFrame['Growth_Rate'][index]):
		sdFrame = sdFrame.drop(index)

# export data to txt file
sdFrame.to_csv("data")

# plot model vs experimental data
plt.plot(sdFrame['% change'], sdFrame['M_Fitness'], 'o')
plt.xlabel('FBA Fitness')
plt.ylabel('Experimental Fitness')
plt.suptitle('S-lim')
plt.show()