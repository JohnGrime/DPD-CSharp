# DPD C#

_A simple dissipative particle dyanmics code written in C#_

**Note: this document is obviously incomplete at this stage!**

## Requirements

* [.NET Core](https://www.microsoft.com/net/learn/get-started/macos)

## Platform compatibility

In principle, anything that supports .NET Core!

## Example usage

For a quick and simple test, run the provided `example.dpd` simulation:

	dotnet run -c Release example.dpd

On my machine, this produced the following (partial) output:

	DPD simulation parameters:
		step_no is 1
		max_steps is 100000
		save_every is 1000
		print_every is 1000
		delta_t is 0.01
		ran1 initialised with -1
		rcut is 1
		lambda is 0.65
		sigma is 3
		kBT is 1
		cell is 10, 10, 10
	3 site types:
		w; (internal 0)
		h; (internal 1)
		t; (internal 2)
	2 molecule types:
		lipid; 115 counts in system:
			13 sites
				h, (id 1)
				h, (id 1)
				h, (id 1)
				t, (id 2)
				t, (id 2)
				t, (id 2)
				t, (id 2)
				t, (id 2)
				t, (id 2)
				t, (id 2)
				t, (id 2)
				t, (id 2)
				t, (id 2)
			12 bonds
				1 - 2 : eq length 0.700, k 100.000
				2 - 3 : eq length 0.700, k 100.000
				4 - 5 : eq length 0.700, k 100.000
				5 - 6 : eq length 0.700, k 100.000
				6 - 7 : eq length 0.700, k 100.000
				7 - 8 : eq length 0.700, k 100.000
				9 - 10 : eq length 0.700, k 100.000
				10 - 11 : eq length 0.700, k 100.000
				11 - 12 : eq length 0.700, k 100.000
				12 - 13 : eq length 0.700, k 100.000
				3 - 4 : eq length 0.700, k 100.000
				3 - 9 : eq length 0.700, k 100.000
			8 angles
				1 - 2 - 2 : eq length 3.142, k 20.000
				3 - 4 - 4 : eq length 3.142, k 20.000
				5 - 6 - 6 : eq length 3.142, k 20.000
				5 - 6 - 6 : eq length 3.142, k 20.000
				7 - 6 - 6 : eq length 3.142, k 20.000
				7 - 8 - 8 : eq length 3.142, k 20.000
				9 - 10 - 10 : eq length 3.142, k 20.000
				11 - 10 - 10 : eq length 1.571, k 3.000
		water; 1505 counts in system:
			1 sites
				w, (id 0)
			0 bonds
			0 angles
	3 interactions (defined in kBT):
	              	           w	           h	           t
		           w	   25.000	   15.000	   80.000
		           h	   15.000	   35.000	   80.000
		           t	   80.000	   80.000	   25.000
	Sites (3000):
		( printing truncated to first 10 sites )

		h (internal 1), position = -0.054, -0.054, -0.054
		h (internal 1), position = 0.404, 0.404, 0.404
		h (internal 1), position = 0.066, 0.066, 0.066
		t (internal 2), position = 0.216, 0.216, 0.216
		t (internal 2), position = -0.434, -0.434, -0.434
		t (internal 2), position = -0.136, -0.136, -0.136
		t (internal 2), position = -0.275, -0.275, -0.275
		t (internal 2), position = 0.015, 0.015, 0.015
		t (internal 2), position = 0.072, 0.072, 0.072
		t (internal 2), position = 0.352, 0.352, 0.352
	Friction coefficient is 4.5 ( as sigma = 3 )
	Bead density is 3 per cubic Rc
	Step 1/100000, sim time         0.01, CPU time 0s:
		Total energy     : 19651557.916507
		Kinetic energy   :     0.000000 ( Average     0.000000, target kBT     1.000000, current kBT     0.000000 )
		Nonbonded energy : 19608421.785025 ( Average     4.796279 from 4088257 collisions )
		Bond energy      :  4527.837056 ( Average     1.640521 )
		Angle energy     : 38608.294426 ( Average    13.988512 )
		Pressure         : 12976292.716443 -71477.737153 75707.862834
		                   -62895.830638 12887770.397620 55693.329895
		                   83465.750873 55693.329895 12881352.649045
		System centre of mass =    -0.007375     0.001635    -0.005913
		Net system momentum   =     0.000000     0.000000     0.000000


	Step 1000/100000, sim time        10.00, CPU time 13s:
		Total energy     : 19001.265820
		Kinetic energy   :  5208.739275 ( Average     1.736246, target kBT     1.000000, current kBT     1.157884 )
		Nonbonded energy : 11933.669336 ( Average     0.715662 from 16675 collisions )
		Bond energy      :   989.162290 ( Average     0.358392 )
		Angle energy     :   869.694919 ( Average     0.315107 )
		Pressure         :    26.545026    -0.053005     0.202331
		                       0.124360    22.807272    -0.020332
		                       0.668986    -0.020332    23.542248
		System centre of mass =     0.009292    -0.001699    -0.029246
		Net system momentum   =     0.000000     0.000000     0.000000


## Notes
