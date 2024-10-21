sp = openmc.StatePoint('statepoint.100.h5')
tally = sp.get_tally()
print(tally)

#Getting the heat
heat_scores = sp.get_tally(scores=['kappa-fission'])
flux_scores = sp.get_tally(scores=['flux'])

print(heat_scores)
print(flux_scores)

#Reshaping
heat_scores.std_dev.shape = mesh_dimensions
heat_scores.mean.shape = mesh_dimensions

flux_scores.std_dev.shape = mesh_dimensions
flux_scores.mean.shape = mesh_dimensions

#Plotting
#lim = 500
fig1 = plt.subplot()
#fig1.set_xlim(500-lim,500+lim)
#fig1.set_ylim(500-lim,500+lim)
fig1.imshow(heat_scores.mean)
plt.savefig("heatscores.png", dpi=600)
plt.show()

fig2 = plt.subplot()
fig2.imshow(flux_scores.mean)
plt.savefig("fluxscores.png", dpi=600)
plt.show()