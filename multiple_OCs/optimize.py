import cobra.io.sbml as IO

# Optimizes BM/S ratio for model
# Parameters: takes cobrapy model, uptake reaction of lim substrate, and 
# biomass reaction
# Returns: dictionary of reaction fluxes
def optimize(model, uptakeRxn, biomassRxn): 
	# Base case, OC(BM) = 1
	model.optimize()
	# Base case for OCs
	uptakeRxn.objective_coefficient = 1
	biomassRxn.objective_coefficient = 1
	# Change OC (biomass) until objective function reaches desired tolerance
	# Assumes that solving for OC(BM) at every iteration converges 
	count = 0
	while(abs(model.solution.f) > 0.01):
		biomassRxn.objective_coefficient = -model.solution.x_dict[uptakeRxn.id] * \
			uptakeRxn.objective_coefficient / model.solution.x_dict[biomassRxn.id]
		model.optimize()
		count = count + 1
	# print(count)
	return model.solution.x_dict

def main():
	yeast = IO.create_cobra_model_from_sbml_file('../SBML/Yeast_YPGly_Medium.xml', 
		print_time=True)
	reactionNum = 553 # reaction number of glycerol (substrate) uptake
	uptakeRxn = yeast.reactions[reactionNum] # uptake of substrate
	biomassRxn = yeast.reactions[1574] # biomass reactions
	xdict = optimize(yeast, uptakeRxn, biomassRxn)
	# Output ratio BM/S
	print("Flux ratio BM/S: " + str(-xdict[biomassRxn.id] / xdict[uptakeRxn.id]))
	print("BM flux: " + str(xdict[biomassRxn.id]))

if __name__ == "__main__":
	main()