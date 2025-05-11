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


    def _intersect(self, vertex1_x, vertex1_y, vertex2_x, vertex2_y, p_x, p_y, tan_vec_x, tan_vec_y):

        """
        Checks whether the half line r starting at the point p = (p_x, py) and with 
        tangent vector t = (tan_vec_x, tan_vec_y) crosses the segment s defined by the 
        points v1 = (vertex1_x, vertex1_y) and v2 = (vertex2_x, vertex2_y). It solves
        the linear system defined by defined by the crossing of r and s.

        Returns True if the half line r intersects the line segment s and False otherwise.


        Keyword arguments:
            vertex1_{x/y} (float):
                The {x/y}-coordinate of the vertex 1.

            vertex2_{x/y} (float):
                The {x/y}-coordinate of the vertex 2.

            p_{x/y} (float):
                The {x/y}-coordinate of the half line starting point p.

            tan_vec_{x/y} (float):
                The {x/y}-component of the tangent vector of the half 
                line starting point p.


        Returns: 
            bool: 
                Whether the half line r crosses the segment s.
        """

        # segment tangent vector components
        l_x = vertex2_x - vertex1_x
        l_y = vertex2_y - vertex1_y

        # determinant of the coefficient matrix
        det = -tan_vec_x*l_y + tan_vec_y*l_x
        if det == 0.:
            return False
        else: 
            alpha = (tan_vec_y*(p_x - vertex1_x) - tan_vec_x*(p_y - vertex1_y))/det
            beta = (l_y*(p_x - vertex1_x) - l_x*(p_y - vertex1_y))/det

            return (np.abs(det) > np.finfo(float).eps)&(alpha >= 0.)&(alpha <= 1.)&(beta >= 0.)
        

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

        For a given point p, it checks how many times the half line starting 
        at p with a random tangent vector crosses the sides of the polygon. An
        odd number of crosses indicates that p is inside the polygon and an even 
        number of crosses indicates p is outside.

        Returns Irue if the point is inside the polygon and False otherwise.

        Keyword arguments:
            p_x (float):
                The x-coordinate of the point one whishes to check.

            p_y (float):
                The y-coordinate of the point one whishes to check.


        Returns:
            bool

        """

        vector_x, vector_y = 2.*np.random.rand(2) - 1.
        count = 0

        for i in range(self.n_vertices):
            count += self._intersect(self.vertices_x[i], 
                                     self.vertices_y[i], 
                                     self.vertices_x[i+1], 
                                     self.vertices_y[i+1], 
                                     p_x, p_y, 
                                     vector_x, vector_y
                                    )
            
        if count % 2:
            return True
        else:
            return False
