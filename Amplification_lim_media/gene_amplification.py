import cobra.io.sbml
from cobra.flux_analysis import single_deletion

# Takes gene, COBRApy model, and file. Prints to console and file the resulting
# growth rate after increasing the bound of all reactions associated with gene
# by 1%
def amp(gene, model, f, original, biomassRxn):
	rxnList = []
	for rxn in model.reactions:
		if gene in rxn.genes and not rxn in rxnList:
			rxn.upper_bound = rxn.upper_bound * 1.01 # 1% increase
			rxn.lower_bound = rxn.lower_bound * 1.01 # for uptake reactions
			rxnList.append(rxn)
	model.optimize()
	growthRate = model.solution.x_dict[biomassRxn.id]
	for rxn in rxnList: # reset bounds
		rxn.upper_bound = rxn.upper_bound / 1.01
		rxn.lower_bound = rxn.lower_bound / 1.01
	print(gene.id + ": " + str(abs(growthRate - original) / original))
	f.write(gene.id + ": " + str(abs(growthRate - original) / original) + "\n")

# Write all amplification data for a model to file
def makeData(model, name):
	fileName = "amp_data_" + name + "lim"
	lim_file = open(fileName, "w")
	model.optimize()
	biomassRxn = model.reactions[1574] # biomass reaction
	original = model.solution.x_dict[biomassRxn.id]
	lim_file.write(name + 
		"-lim: Percent changes in growth rate for 1% amplification of genes\n")
	for gene in model.genes:
		amp(gene, model, lim_file, original, biomassRxn)
	lim_file.close()

# Perform all gene amplifications for phosphate, sulfate, and 
# glucose limited media
def main():
	P_model = cobra.io.sbml.create_cobra_model_from_sbml_file('../SBML/P-lim.xml', 
		print_time=True)
	S_model = cobra.io.sbml.create_cobra_model_from_sbml_file('../SBML/S-lim.xml', 
		print_time=True)
	D_model = cobra.io.sbml.create_cobra_model_from_sbml_file('../SBML/D-lim.xml', 
		print_time=True)
	makeData(P_model, "P")
	makeData(S_model, "S")
	makeData(D_model, "D")

if __name__ == "__main__":
	main()