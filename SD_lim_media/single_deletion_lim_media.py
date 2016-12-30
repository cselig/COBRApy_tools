import cobra.io.sbml as IO
from cobra.flux_analysis import single_deletion
import pandas

# Params:
#	model: a cobrapy model built from SBML file
#	name: string of limited nutrient eg "P"
# Writes growth rates of all single deletion mutants to file
def make_SD_data(model, name):
	outFile = "SD_data_" + name + "lim"
	single_deletion_data = open(outFile, 'w')
	for x in range(1,905):
		growth_rates, statuses = single_deletion(model, element_list = model.genes[x-1:x])
		table = pandas.DataFrame.from_dict({"growth_rates": growth_rates,
		"statuses": statuses})
		single_deletion_data.write(str(x) + ": " + str(growth_rates) + '\n')
		print table
	single_deletion_data.close()

# Performs single deletions on all genes for each of the three 
# limited media, S, P, and D-lim and writes fitnesses to a file
def main():
	P_model = IO.create_cobra_model_from_sbml_file('../SBML/P-lim.xml', print_time=True)
	S_model = IO.create_cobra_model_from_sbml_file('../SBML/S-lim.xml', print_time=True)
	D_model = IO.create_cobra_model_from_sbml_file('../SBML/D-lim.xml', print_time=True)
	make_SD_data(P_model, "P")
	make_SD_data(S_model, "S")
	make_SD_data(D_model, "D")

if __name__ == "__main__":
	main()





