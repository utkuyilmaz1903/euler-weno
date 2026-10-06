# Design decisions

Notes on why the code looks the way it does, mostly for future me.

## Fluxes get their own type

I first used `Conserved` for fluxes too, since both have three numbers. I split them because they mean different things: `Conserved` is what sits in a cell (per unit volume), `Flux` is what goes through a face (per unit area, per second). Adding one to the other makes no sense, and now the compiler stops me if I try. The only place a flux turns into a change of state is the update in `advance`, and there it is written out field by field.

The face flux (`rusanov_flux`) and the flux of a single state (`physical_flux`) stay the same type, because they measure the same thing in the same units. They are told apart by their names.

## The scheme picks the number of ghost cells

At first the ghost cell count was a number you passed in. I didn't like that: whoever runs the solver would have to know that first order needs 1 ghost cell per side and WENO5 needs 3. Now you only pick the scheme (`Reconstruction`) and `ghost_cells_per_side` works out the count, so nobody has to think about it or can get it wrong.
