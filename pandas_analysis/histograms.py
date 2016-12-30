# plots histogram for experimental growth rates

import pandas
import cobra.io.sbml
from cobra.flux_analysis import single_deletion
import matplotlib.pyplot as plt
import math

# Create model, record wt fitness
model = cobra.io.sbml.create_cobra_model_from_sbml_file('../SBML/S-lim.xml', 
	print_time=True)
biomassRxn = model.reactions[1574] # biomass reaction

# store wild type fitness
model.optimize()
wtFitness = model.solution.x_dict[biomassRxn.id]

# perform deletions
growth_rates, statuses = single_deletion(model, element_list = model.genes[1:905])
sdFrame = pandas.DataFrame.from_dict(growth_rates, 'index')
sdFrame.columns = ['Growth_Rate']

sdFrame['wt_Fitness'] = pandas.Series()
sdFrame['% change'] = pandas.Series()

# calculate % changes for experimental data
for index, row in sdFrame.iterrows():
	row['wt_Fitness'] = wtFitness
	row['% change'] = (row['Growth_Rate'] - wtFitness) / wtFitness

# plot histograms
plt.hist(sdFrame['% change'])
plt.title("S-lim")
plt.xlabel("Growth Rate")
plt.ylabel("Frequency")
plt.show()
