import numpy as np

__all__ = ["PolygonS"]

class PolygonS(object):

    """
    Bi-dimensional spherical polygon generator.

    PolygonS generates a polygon in bi-dimensional spherical surface 
    with vertices regularly or irregularly placed.  

    
    Attributes:
        n_vertices (int): 
            The number of vertices.

        anchor_point (tuple): 
            The tuple with the latitude and longitude of the anchor point, in degrees.

        is_anchor_point_inside (bool): 
            It is True if the anchor point is inside the polygon.
        
        direction (str, default = 'anticlockwise'): 
            The direction of the circuit formed by the sides of the polygon.
            It can assume the values 'clockwise' or 'anticlockwise'.
        
        vertices_lat (np.ndarray): array of shape (n_vertices + 1,).
            The vertices latitude, in degrees.

        vertices_lon (np.ndarray): array of shape (n_vertices + 1,).
            The vertices longitude, in degrees.

            
    Methods:
        create(anchor_point, n_vertices, angular_dist_min, angular_dist_max, angular_separation, direction):
            Creates the polygon by defining its vertices latitudes and longitudes, in degrees.

        internal_point(p_lat, p_lon):
            Checks weather the point p = (p_lat, p_lon) is inside or outside the polygon, 
            with the point's latitude and longitude given in degrees.

    """

    def __init__(self):
        """
        Initialises the PolygonS class.

        """
        
        self.n_vertices = None
        self.direction = None
        self.vertices_lat = None
        self.vertices_lon = None
        self.anchor_point = None
        self.is_anchor_point_inside = None


    def _3d_vector(self, lat, lon):
        """
        Creates a 3-d radius vector in the embeding Eclidean space with norm one
        normal to the 2-d spherical surface.

        Returns the array of with the Cartesian components of the 
        radius vector.

        Keyword arguments:
            lat (float):
                The latitude of the point on the spherical surface in degrees.

            lon (float):
                The longitude of the point on the spherical surface in degrees.

        Returns:
            array (np.ndarray): array of shape (3,).
        """

        lat_rad = np.radians(lat)
        lon_rad = np.radians(lon)

        x = np.sin(np.pi/2. - lat_rad)*np.cos(lon_rad)
        y = np.sin(np.pi/2. - lat_rad)*np.sin(lon_rad)
        z = np.cos(np.pi/2. - lat_rad)

        return x, y, z
    

    def _tangent_space_basis(self, lat, lon):
        """
        Creates the orthonormal basis of the tangent space at a given point
        on the spherical surface.

        Returns the 2-d orthonormal basis of the tangent space as vectors 
        in the embeding 3-d Euclidean space.

        Keyword arguments:
            lat (float):
                The latitude of the point on the spherical surface in degrees.

            lon (float):
                The longitude of the point on the spherical surface in degrees.

        Returns:
            tuple of arrays (np.ndarray, np.ndarray):
                A tuple of arrays of shape (3,) and (3,) with the 
                xyz-components of the vectors making the orthonomal basis.
        """

        lat_rad = np.radians(lat)
        lon_rad = np.radians(lon)

        e1_x = np.cos(np.pi/2. - lat_rad)*np.cos(lon_rad)
        e1_y = np.cos(np.pi/2. - lat_rad)*np.sin(lon_rad)
        e1_z = -np.sin(np.pi/2. - lat_rad)

        e2_x = -np.sin(np.pi/2. - lat_rad)*np.sin(lon_rad)
        e2_y = np.sin(np.pi/2. - lat_rad)*np.cos(lon_rad)
        e2_z = 0.

        return np.array([e1_x, e1_y, e1_z]), np.array([e2_x, e2_y, e2_z])


    def _intersect(
            self, 
            vertices_lat, 
            vertices_lon, 
            fixed_point_lat, 
            fixed_point_lon, 
            point_lat, 
            point_lon
            ):
        
        """
        Checks whether the shorter segment of the great circle r, connecting the fixed point
        x = (fixed_point_lat, fixed_point_lon) to the point of interest p = (point_lat, point_lon),
        crosses the shorter segment of the great circle s defined by the points 
        v1 = (vertex1_x, vertex1_y) and v2 = (vertex2_x, vertex2_y). The check is performed 
        in the embeding 3-d Euclidean space. It compares the components os radius vectors of 
        the two great circles at the intersection point using the basis of the plane that 
        contains the great circles formed by the radius and tangent vectors of the great circle.
        
        Returns True if the shorter segment of the great circle r connecting the fixed point x 
        and the point of interest p intersects the shorter segment of the great circle segment s 
        defined by the vertices v1 and v2.


        Keyword arguments:
            vertex1_{lat/lon} (float):
                The {latitude/longitude} of the vertex 1 in degrees.

            vertex2_{lat/lon} (float):
                The {latitude/longitude} of the vertex 2 in degrees.

            fixed_point_{lat/lon} (float):
                The {latitude/longitude} of the fixed point in degrees.

            p_{lat/lon} (float):
                The {latitude/longitude} of the point of interest p in degrees.


        Returns: 
            bool: 
                Whether the great circle r crosses the great circle s.

        """

        vertices_x, vertices_y, vertices_z = self._3d_vector(vertices_lat, vertices_lon)
        q_x, q_y, q_z = self._3d_vector(np.repeat(fixed_point_lat, self.n_vertices), np.repeat(fixed_point_lon, self.n_vertices))
        p_x, p_y, p_z = self._3d_vector(np.repeat(point_lat, self.n_vertices), np.repeat(point_lon, self.n_vertices))

        cos_angle_v1v2 = vertices_x[:self.n_vertices]*vertices_x[1:] + vertices_y[:self.n_vertices]*vertices_y[1:] + vertices_z[:self.n_vertices]*vertices_z[1:]
        cos_angle_qp = q_x*p_x + q_y*p_y + q_z*p_z

        norm_v1v2 = np.sqrt(1. - cos_angle_v1v2**2)
        t_12_x = (vertices_x[1:] - cos_angle_v1v2*vertices_x[:self.n_vertices])/norm_v1v2
        t_12_y = (vertices_y[1:] - cos_angle_v1v2*vertices_y[:self.n_vertices])/norm_v1v2
        t_12_z = (vertices_z[1:] - cos_angle_v1v2*vertices_z[:self.n_vertices])/norm_v1v2
        
        norm_qp = np.sqrt(1. - cos_angle_qp**2)
        t_qp_x = (p_x - cos_angle_qp*q_x)/norm_qp
        t_qp_y = (p_y - cos_angle_qp*q_y)/norm_qp
        t_qp_z = (p_z - cos_angle_qp*q_z)/norm_qp

        angle_v1v2 = np.arccos(cos_angle_v1v2)
        angle_qp = np.arccos(cos_angle_qp)

        cos_angle_qv1 = q_x*vertices_x[:self.n_vertices] + q_y*vertices_y[:self.n_vertices] + q_z*vertices_z[:self.n_vertices]
        cos_angle_tqpv1 = t_qp_x*vertices_x[:self.n_vertices] + t_qp_y*vertices_y[:self.n_vertices] + t_qp_z*vertices_z[:self.n_vertices]
        cos_angle_qt12 = q_x*t_12_x + q_y*t_12_y + q_z*t_12_z
        cos_angle_tqpt12 = t_qp_x*t_12_x + t_qp_y*t_12_y + t_qp_z*t_12_z

        alpha = np.arctan2(1. - cos_angle_qv1**2 - cos_angle_tqpv1**2,cos_angle_qt12*cos_angle_qv1 + cos_angle_tqpt12*cos_angle_tqpv1) 
        alpha_prime = np.arctan2(1. - cos_angle_qv1**2 - cos_angle_qt12**2,cos_angle_tqpv1*cos_angle_qv1 + cos_angle_tqpt12*cos_angle_qt12)

        inside = (alpha >= 0.)&(alpha <= angle_v1v2)&(alpha_prime >= 0)&(alpha_prime <= angle_qp)
        
        return inside       
    

    def create(self, anchor_point, n_vertices, angular_dist_min, angular_dist_max, angular_separation = 'regular', direction='clockwise'):

        """
        Creates a polygon in the 2-d spherical surgace with radius 1. The polygon is build
        by first choosing an anchor point on that surface and then "launching" the vertices
        from it. The vertices are launched along a great circle defined by the chosend angular 
        distance from the anchor point and the tangent vector at the anchor point. The angular
        distance between the anchor point and the vertex can be random within the interval 
        [angular_dist_min, angular_dist_max] and angular separation between the vertices
        can be random by generating a tangent vector at a random direction or with fixed 
        angular distance when angular_dist_min = angular_dist_max and angular separation 
        2*Pi/n_vertices.

        Keyword arguments:
            anchor_point (array-like): array of shape (2,).
                The latitude and longitude of the anchor point, in degrees.

            n_vertices (int, default = 3):
                The number of vertices.

            angular_dist_min (float): 
                The minimum angular distance between the anchor point and vertices, in degrees.

            angular_dist_max (float): 
                The maximum angular distance between the anchor point and vertices, in degrees.

            angular_separation (str, default = 'regular'): the possible values 'regular' or 'random'.
                Whether the angular separation between two consecutive vertices is regular, i.e.
                is 2*Pi/n_vertices, or random in the range [0, 2*Pi[.

            direction (str, default = 'anticlockwise'): the possible values are 'clockwise' or 'anticlockwise'.
                Whether the sequence of vertices runs clockwise or anticlockwise.

        Returns:
            self (object):
                The polygon.
        
        """

        self.n_vertices = n_vertices
        self.direction = direction
        self.anchor_point = anchor_point

        lat_0, lon_0 = self.anchor_point

        r_0_x, r_0_y, r_0_z = self._3d_vector(lat_0, lon_0)
        r_0 = np.array([r_0_x, r_0_y, r_0_z])
        e1_hat, e2_hat = self._tangent_space_basis(lat_0, lon_0)

        if angular_separation == 'regular':
            alpha = np.linspace(0., 2.*np.pi, n_vertices, endpoint=False)

        elif angular_separation == 'random':
            alpha = np.sort(2.*np.pi*np.random.rand(n_vertices))

        t_0 = np.array([np.cos(alpha_i)*e1_hat + np.sin(alpha_i)*e2_hat for alpha_i in alpha])

        beta = (angular_dist_max - angular_dist_min)*np.random.rand(n_vertices) + angular_dist_min
        beta = np.radians(beta)

        if (beta < 0.).sum():
            self.is_anchor_point_inside = False
        
        else:
            self.is_anchor_point_inside = True

        r = np.array([np.cos(beta_i)*r_0 + np.sin(beta_i)*t_0_i for beta_i, t_0_i in zip(beta, t_0)])

        self.vertices_lat = np.array([np.pi/2. - np.arccos(z_i) for z_i in r.T[2]])
        self.vertices_lon = np.array([np.arctan2(y_i,x_i) for x_i, y_i in zip(r.T[0], r.T[1])])

        self.vertices_lat = np.append(self.vertices_lat, self.vertices_lat[0])
        self.vertices_lon = np.append(self.vertices_lon, self.vertices_lon[0])

        self.vertices_lat = np.rad2deg(self.vertices_lat)
        self.vertices_lon = np.rad2deg(self.vertices_lon)

        if direction == 'clockwise':
            self.vertices_lat = np.flip(self.vertices_lat)
            self.vertices_lon = np.flip(self.vertices_lon)
        elif direction == 'anticlockwise':
            pass
        
        
        self.vertices = np.array([elem for elem in zip(self.vertices_lat, self.vertices_lon)])

        return

    
    def internal_point(self,p_lat, p_lon):

        """
        Checks whether the point p = (p_lat, p_lon) is inside the polygon.

        For a given point p, it checks how many times the shorter segment of 
        the great circle connecting p to the anchor point. If the anchor 
        point is inside (outside) the polygon, then an odd number of 
        crossings indicates that p is inside (outside) the polygon and 
        an even number of crossings indicates p is outside (inside).

        Returns Irue if the point is inside the polygon and False otherwise.

        Keyword arguments:
            p_x (float):
                The x-coordinate of the point one whishes to check.

            p_y (float):
                The y-coordinate of the point one whishes to check.


        Returns:
            bool
            
        """

        count = self._intersect(
            self.vertices_lat, 
            self.vertices_lon, 
            self.anchor_point[0], 
            self.anchor_point[1], 
            p_lat, 
            p_lon
            ).sum()
        
        if self.is_anchor_point_inside:
            if count % 2:
                return False
            else:
                return True
        
        else:
            if count % 2:
                return True
            else:
                return False
