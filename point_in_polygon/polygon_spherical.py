import numpy as np
from prettytable import PrettyTable

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
        
        direction (str, default = 'anticlockwise'): 
            The direction of the circuit formed by the sides of the polygon.
            It can assume the values 'clockwise' or 'anticlockwise'.
        
        vertices_lat (np.ndarray): array of shape (n_vertices + 1,).
            The vertices latitude, in degrees.

        vertices_lon (np.ndarray): array of shape (n_vertices + 1,).
            The vertices longitude, in degrees.

        inside (list): a list of bools.
            The check of whether points are inside or outside 
            the polygon.


            
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
        self.inside = None


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
            point_lat, 
            point_lon
            ):
        
        """
        Checks whether the west-east great circle starting at the point p = (point_lat, point_lon)
        crosses the sides of the polygon (i.e. any great circle segment between two consecutive 
        polygon vertices). The number of crossings is determined by checking 
        point_lon with respect to the longitude of the crossing point, if it exists.

        Returns an array-like of bools, where True means a crossing.


        Keyword arguments:
            vertices_{lat/lon} (float):
                The {latitude/longitude} of the polygon vertices in degrees.

            p_{lat/lon} (float):
                The {latitude/longitude} of the point of interest p in degrees.


        Returns: 
            crossed (np.ndarray): ndarray of bools with shape (n_vertices,) 
                Whether the great circle starting at p crossed the sides of the polygon.

        """

        vertices_x, vertices_y, vertices_z = self._3d_vector(vertices_lat, vertices_lon)
        p_x, p_y, p_z = self._3d_vector(np.repeat(point_lat, self.n_vertices), np.repeat(point_lon, self.n_vertices))
        cos_angle_v1v2 = vertices_x[:self.n_vertices]*vertices_x[1:] + vertices_y[:self.n_vertices]*vertices_y[1:] + vertices_z[:self.n_vertices]*vertices_z[1:]

        epsilon = 1e-15
        
        norm_v1v2 = np.sqrt(1. - cos_angle_v1v2**2)
        norm_v1v2 += (norm_v1v2 < np.finfo(float).eps)*epsilon
        t_12_x = (vertices_x[1:] - cos_angle_v1v2*vertices_x[:self.n_vertices])/norm_v1v2
        t_12_y = (vertices_y[1:] - cos_angle_v1v2*vertices_y[:self.n_vertices])/norm_v1v2
        t_12_z = (vertices_z[1:] - cos_angle_v1v2*vertices_z[:self.n_vertices])/norm_v1v2

        t_p_x = -np.sin(np.pi/2. - np.radians(point_lat))*np.sin(np.radians(point_lon))
        t_p_y = np.sin(np.pi/2. - np.radians(point_lat))*np.cos(np.radians(point_lon))
        t_p_z = 0

        cos_angle_pv1 = p_x*vertices_x[:self.n_vertices] + p_y*vertices_y[:self.n_vertices] + p_z*vertices_z[:self.n_vertices]
        cos_angle_tpv1 = t_p_x*vertices_x[:self.n_vertices] + t_p_y*vertices_y[:self.n_vertices] + t_p_z*vertices_z[:self.n_vertices]
        cos_angle_pt12 = p_x*t_12_x + p_y*t_12_y + p_z*t_12_z
        cos_angle_tpt12 = t_p_x*t_12_x + t_p_y*t_12_y + t_p_z*t_12_z

        alpha_prime = np.arctan2(1. - cos_angle_pv1**2 - cos_angle_pt12**2,cos_angle_tpv1*cos_angle_pv1 + cos_angle_tpt12*cos_angle_pt12)

        cond1 = (vertices_lat[:self.n_vertices] > point_lat) != (vertices_lat[1:] > point_lat)
        cond2 = np.radians(point_lon) < alpha_prime

        crossed = cond1 & cond2
        
        return crossed       
    

    def create(self, anchor_point, n_vertices, angular_dist_min, angular_dist_max, angular_separation = 'regular', direction='clockwise'):

        """
        Creates a polygon in the 2-d spherical surface with radius 1. The polygon is build
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

    
    def internal_point(self, p , labels = None):

        """
        Checks whether points p are inside the polygon.

        For a given point, it counts how many times the great circle 
        starting at it crosses the sides of the polygon. An odd number of crosses 
        indicates that the point is inside the polygon and an even number of crosses 
        indicates it is outside.

        Returns True if the point is inside the polygon and False otherwise.

        Keyword arguments:
            p (array): array-like of shape (2,) in the case of one point or
            of shape (n_points, 2) in the case of n_points.
                The latitude/longitude of the point or a list of pairs for a 
                number of points

            labels(array, default = None): array-like of shape (n_points,)
                The labels of the points.

        Returns:
            Table of the point labels, their coordinates and whether they 
            are inside or ouside the polygon.
            
        """

        self.inside = []

        if len(p) == 2:
            points = np.array([p])
        elif (len(p) > 2):
            points = np.array([i for i in p])

        n_points = len(points)

        for pt in points:

            crossings = self._intersect(
                self.vertices_lat, 
                self.vertices_lon,
                pt[0], 
                pt[1]
                ).sum()
            
            if crossings % 2:
                self.inside.append(True)
            else:
                self.inside.append(False)

        table = PrettyTable(["point",  "(lat, lon)", "is inside"])
        if labels:
            for i in range(n_points):
                table.add_row([labels[i], f"{(points[i][0], points[i][1])}", self.inside[i]])

        else:
            for i in range(n_points):
                table.add_row([i + 1, f"{(points[i][0], points[i][1])}", self.inside[i]])

        return table
