=== regions{n}.png ===
The darkest region is region 1.
The lightest region is region n, with n the amount of islands.
This includes fixed islands, as they are also their own separate region.

=== table{n}.txt ===
The output table for the corresponding geometry from regions{n}.png.
This contains the columns
- (t): Time, but there is no time stepping
- (mx, my, mz): Average magnetization of the whole simulation
- E_total: Total energy, which is of course very important for the balancedness.
- E_Zeeman: This only concerns B_ext, which is only present in the 'fixed' islands to keep them fixed. Hence, the energy we are interested in is E_total-E_Zeeman.
- a{m}: Initial magnetization angle of island with region index {m}. Fixed islands are NOT listed here.
- m.region{m}x, m.region{m}y, m.region{m}z: Relaxed magnetization angle of island with reigon index {m}, including fixed islands.
- roundness: Roundness of every island, because for now every island is exactly the same.
- size: Major axis of the constituent ellipses of all islands, because for now every island is exactly the same.
- Cell_size: Cell size of the simulation grid.

=== log{n}.txt ===
The mumax3 log.txt file that was generated when calculating table{n}.txt.