import numpy as np

__all__ = ["PolygonE"]

class PolygonE(object):

    """
    Bi-dimensional Euclidean polygon generator.

    PolygonE generates a polygon in bi-dimensional Euclidean space 
    with vertices regularly or irregularly placed.  

    
    Attributes:
        n_vertices (int): 
            The number of vertices.
        
        direction (str, default = 'anticlockwise'): 
            The direction of the circuit formed by the sides of the polygon.
            It can assume the values 'clockwise' or 'anticlockwise'.
        
        vertices_x (np.ndarray): array of shape (n_vertices + 1,).
            The vertices x-coordinate.

        vertices_y (np.ndarray): array of shape (n_vertices + 1,).
            The vertices y-coordinate.

            
    Methods:

        create(anchor_point, n_vertices, dist_min, dist_max, angular_separation, direction):
            Creates the polygon by defining its vertices xy-coordinates.

        internal_point(p_x, p_y):
            Checks weather the point p = (p_x, p_y) is inside or outside the polygon.

    """

    def __init__(self):

        """
        Initialises the PolygonE class.

        """
        
        self.n_vertices = None
        self.direction = None
        self.vertices_x = None
        self.vertices_y = None


    def _intersect(
            self, 
            vertices_x, 
            vertices_y,  
            p_x, 
            p_y
            ):

        """
        Checks whether the horizontal half line starting at the point p = (p_x, py) and positive
        oriented with respect to the y-axis crosses the sides of the polygon (i.e. the segments 
        formed by two consecutives vertices). The number of crossings is determined by checking 
        the p_x with respect to the x-coordinate of the crossing point, if it exists.

        Returns an array-like of bools, where True means a crossing.


        Keyword arguments:
            vertices_{x/y} (float):
                The {x/y}-coordinate of the polygon vertices.

            p_{x/y} (float):
                The {x/y}-coordinate of the half line starting point p.


        Returns: 
            crossed (np.ndarray): ndarray of bools with shape (n_vertices,) 
                Whether the horizontal half-line starting at p crossed the sides of the polygon.
        """

        epsilon = 1e-15

        cond1 = (vertices_y[:self.n_vertices] > p_y) != (vertices_y[1:] > p_y)

        x_intersect = (vertices_x[1:] - vertices_x[:self.n_vertices]
                       )*(p_y - vertices_y[:self.n_vertices]
                          )/(vertices_y[1:] - vertices_y[:self.n_vertices] + epsilon
                             ) + vertices_x[:self.n_vertices]
        cond2 = p_x < x_intersect
        crossed = cond1 & cond2

        return crossed
        

    def create(self, anchor_point = (0.,0.), n_vertices=3, dist_min=0.5, dist_max=1., angular_separation = 'regular', direction='anticlockwise'):

        """
        Creates a polygon in bi-dimensional Euclidean space. The polygon is build
        by first choosing an anchor point and then placing the vertices by 
        constructing the radius vectors connecting the anchor point to the vertices 
        with random length within the interval [dist_min, dist_max] and/or random 
        angular separation in the interval [0, 2*Pi[ or with fixed length when 
        dist_min = dist_max and angular distance 2*Pi/n_vertices.

        Keyword arguments:
            anchor_point (array-like, default = (0., 0.)): array of shape (2,).
                The xy-coordinates of the point from which the polygon will be built.

            n_vertices (int, default = 3):
                The number of vertices.

            dist_min (float, default = 0.5): 
                The minimum length of the radius vector connecting the anchor point to the vertex.

            dist_max (float, default = 1.): 
                The maximum length of the radius vector connecting the anchor point to the vertex.

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

        if angular_separation == 'regular':
            theta = np.linspace(0,2.*np.pi,self.n_vertices, endpoint = False)
        
        elif angular_separation == 'random':
            theta = np.sort(2.*np.pi*np.random.rand(self.n_vertices))
        
        radius = (dist_max - dist_min)*np.random.rand(self.n_vertices) + dist_min
        x_0, y_0 = anchor_point

        if direction=='clockwise':
            self.vertices_x = radius*np.cos(theta) + x_0
            self.vertices_x = np.append(self.vertices_x, self.vertices_x[0])
            self.vertices_y = radius*np.sin(theta) + y_0
            self.vertices_y = np.append(self.vertices_y, self.vertices_y[0])

            self.vertices_x = np.flip(self.vertices_x)
            self.vertices_y = np.flip(self.vertices_y)

        elif direction in ['anticlockwise', 'counterclockwise']:
            self.vertices_x = radius*np.cos(theta) + x_0
            self.vertices_x = np.append(self.vertices_x, self.vertices_x[0])
            self.vertices_y = radius*np.sin(theta) + y_0
            self.vertices_y = np.append(self.vertices_y, self.vertices_y[0])
    
    
        return
    

    def internal_point(self, p_x, p_y):
        
        """
        Checks whether the point p = (p_x, p_y) is inside the polygon.

        For a given point p, it counts how many times the horizontal half-line 
        starting at p crosses the sides of the polygon. An odd number of crosses 
        indicates that p is inside the polygon and an even number of crosses 
        indicates p is outside.

        Returns True if the point is inside the polygon and False otherwise.

        Keyword arguments:
            p_x (float):
                The x-coordinate of the point one whishes to check.

            p_y (float):
                The y-coordinate of the point one whishes to check.


        Returns:
            bool

        """

        crossings = self._intersect(
            self.vertices_x, 
            self.vertices_y,
            p_x, 
            p_y
            ).sum()
            
        if crossings % 2:
            return True
        else:
            return False
