from zernike_like import (base_input, ZernikeParent,
                          ZBasis, KBasis)


class ApproxParent():
    """Parent Class for the ZApprox, KApprox, and HApprox
    classes"""

    def __init__(self, n):
        """
        Parameters
        ----------
        n : int
            Zernike Order
        """
        self.n = n

class KApprox(ApproxParent):

    def __init__(n, num_sides, polygon_radius):
        """
        Parameters
        ----------
        num_sides : int
            number of sides on the regular polygon
        polygon_radius : float
            radius of the polygon
        """
        # inheriting from ApproxParent
        super().__init__(n)

        self.num_sides = num_sides
        self.polygon_radius = polygon_radius
    
