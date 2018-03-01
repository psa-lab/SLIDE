#ifndef _GENERATE_POINTS_H
#define _GENERATE_POINTS_H

extern void  generate_points ( point_pt surface_points,
			       int      number_of_surface_points,
			       point_pt center_points,
			       int      number_of_center_points,
			       FILE     *fp );

extern point_pt  generate_grid_points ( int      *number_of_surface_points,
					point_pt center_points,
					int      number_of_center_points,
					double   grid_spacing,
					FILE     *fp );

#endif
