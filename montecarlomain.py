import montecarlo

a = montecarlo.montecarlo()
a.diffusion_tree()
a.structure_geometry()
a.rebar_calc()
a.crack_calc()
a.derating_calc()
a.parameter_fill()
a.monte_zeroes()
a.monte_fill()
a.histogram()
a.outputs()
a.normalization()
a.plot_damage()

print a.monte_bincount
print a.x

