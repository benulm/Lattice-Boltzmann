#ifndef STANDARD_COLLIDE_
#define STANDARD_COLLIDE_

namespace lb
{

	void simulation::collide()
	{
		for (int j=0; j<static_cast<int>(l.ny); ++j)
		{
			for (int i=0; i<static_cast<int>(l.nx); ++i)
			{
				float_type f_old[9];
				for (unsigned k =0; k<9; ++k) {
					f_old[k] = l.get_node(i,j).f(k);
				}

				velocity_set().equilibrate(l.get_node(i,j));

				for (unsigned k =0; k<9; ++k) {
					l.get_node(i,j).f(k) = f_old[k]+2.*beta*(l.get_node(i,j).f(k)-f_old[k]);
				}
			}
		}
	}

}

#endif // STANDARD_COLLIDE_