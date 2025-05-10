"""Class for plotting arrows to represent genes.

License
-------
This file is part of HomologyViz
BSD 3-Clause License
Copyright (c) 2024, Ivan Munoz Gutierrez
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy import ndarray


class Arrow:
    """
    Generate coordinates for plotting horizontal arrow shapes.

    .. image:: _static/arrow_description.png
       :width: 60%
       :alt: Arrow attribute diagram
       :align: center

    This class computes the (x, y) coordinates needed to plot an arrow
    pointing left or right, representing genes or features along DNA sequences.

    The arrow consists of a rectangular tail and a triangular head. If the
    arrow is shorter than the specified head height, only the triangular head
    is drawn. The head and tail dimensions are customizable.

    Attributes
    ----------
    x1 : int or float
        Start position along the x-axis.
    x2 : int or float
        End position along the x-axis.
    y : int or float
        Vertical position along the y-axis.
    ratio_tail_head_width : float
        Ratio of tail width to head width. Default is 0.5.
    head_width : int or float
        Width of the arrowhead in the y-axis.
    head_height : int or float
        Length of the arrowhead in the x-axis.
    tail_width : float
        Computed width of the tail based on `head_width` and `ratio_tail_head_width`.
    head_shoulder : float
        Distance from the top/bottom of the tail to the top/bottom of the head.
    """

    def __init__(
        self,
        x1: float,
        x2: float,
        y: float,
        ratio_tail_head_width: float = 0.5,
        head_width: float = 2,
        head_height: float = 200,
    ):
        self.x1 = x1
        self.x2 = x2
        self.y = y
        self.ratio_tail_head_width = ratio_tail_head_width
        self.tail_width = head_width * ratio_tail_head_width
        self.head_width = head_width
        self.head_height = head_height
        self.head_shoulder = (head_width - self.tail_width) / 2

    def coordinates_arrow_forward(self) -> tuple[ndarray, ndarray]:
        """
        Compute coordinates for a right-pointing arrow.

        Returns
        -------
        tuple of np.ndarray
            Arrays representing x and y coordinates of the arrow polygon.
        """
        # TODO: Add __repr__ method for better debugging and logging.

        height = self.x2 - self.x1
        # If total height is smaller or equal to the arrow's head hight, then plot
        # only head
        if height <= self.head_height:
            x_1 = self.x1
            x_2 = self.x1
            x_3 = self.x2
            x_4 = self.x1
            x_5 = self.x1
            y_1 = self.y
            y_2 = self.y + (self.head_width * 0.5)
            y_3 = self.y
            y_4 = self.y - (self.head_width * 0.5)
            y_5 = self.y
            x_values = np.array([x_1, x_2, x_3, x_4, x_5])
            y_values = np.array([y_1, y_2, y_3, y_4, y_5])
            return (x_values, y_values)
        # Tail x-values
        x_1 = self.x1
        x_2 = self.x1
        # Head x-values
        x_3 = self.x2 - self.head_height
        x_4 = self.x2 - self.head_height
        x_5 = self.x2
        x_6 = self.x2 - self.head_height
        x_7 = self.x2 - self.head_height
        # Tail x-values
        x_8 = self.x1
        x_9 = self.x1
        # Tail y-values
        y_1 = self.y
        y_2 = self.y + (self.tail_width * 0.5)
        # Head y-values
        y_3 = self.y + (self.tail_width * 0.5)
        y_4 = self.y + (self.tail_width * 0.5) + self.head_shoulder
        y_5 = self.y
        y_6 = self.y - (self.tail_width * 0.5) - self.head_shoulder
        y_7 = self.y - (self.tail_width * 0.5)
        # Tail y-values
        y_8 = self.y - (self.tail_width * 0.5)
        y_9 = self.y
        # make datapoints for plotting
        x_values = np.array([x_1, x_2, x_3, x_4, x_5, x_6, x_7, x_8, x_9])
        y_values = np.array([y_1, y_2, y_3, y_4, y_5, y_6, y_7, y_8, y_9])
        return (x_values, y_values)

    def coordinates_arrow_reverse(self) -> tuple[ndarray, ndarray]:
        """
        Compute coordinates for a left-pointing arrow.

        Returns
        -------
        tuple of np.ndarray
            Arrays representing x and y coordinates of the arrow polygon.
        """
        height = self.x1 - self.x2
        # If total height is smaller or equal to the arrow's head hight plot
        # only head
        if height <= self.head_height:
            x_1 = self.x1
            x_2 = self.x1
            x_3 = self.x2
            x_4 = self.x1
            x_5 = self.x1
            y_1 = self.y
            y_2 = self.y - (self.head_width * 0.5)
            y_3 = self.y
            y_4 = self.y + (self.head_width * 0.5)
            y_5 = self.y
            x_values = np.array([x_1, x_2, x_3, x_4, x_5])
            y_values = np.array([y_1, y_2, y_3, y_4, y_5])
            return (x_values, y_values)
        # Tail x-values
        x_1 = self.x1
        x_2 = self.x1
        # Head x-values
        x_3 = self.x2 + self.head_height
        x_4 = self.x2 + self.head_height
        x_5 = self.x2
        x_6 = self.x2 + self.head_height
        x_7 = self.x2 + self.head_height
        # Tail x-values
        x_8 = self.x1
        x_9 = self.x1
        # Tail y-values
        y_1 = self.y
        y_2 = self.y - (self.tail_width * 0.5)
        # Head y-values
        y_3 = self.y - (self.tail_width * 0.5)
        y_4 = self.y - (self.tail_width * 0.5) - self.head_shoulder
        y_5 = self.y
        y_6 = self.y + (self.tail_width * 0.5) + self.head_shoulder
        y_7 = self.y + (self.tail_width * 0.5)
        # Tail y-values
        y_8 = self.y + (self.tail_width * 0.5)
        y_9 = self.y
        # make datapoints for plotting
        x_values = np.array([x_1, x_2, x_3, x_4, x_5, x_6, x_7, x_8, x_9])
        y_values = np.array([y_1, y_2, y_3, y_4, y_5, y_6, y_7, y_8, y_9])
        return (x_values, y_values)

    def get_coordinates(self) -> tuple[ndarray, ndarray]:
        """
        Get coordinates for either a forward or reverse arrow depending on direction.

        Returns
        -------
        tuple of np.ndarray
            Arrays representing x and y coordinates of the arrow polygon.
        """
        if self.x1 < self.x2:
            return self.coordinates_arrow_forward()
        else:
            return self.coordinates_arrow_reverse()


if __name__ == "__main__":
    ####################
    # Example of usage #
    ####################

    # Make arrow
    arrow1 = Arrow(x1=750, x2=2000, y=10)
    x_values, y_values = arrow1.get_coordinates()
    print(x_values)
    print(y_values)
    plt.fill(x_values, y_values, "b")

    arrow2 = Arrow(x1=1250, x2=250, y=20)
    x_values, y_values = arrow2.get_coordinates()
    plt.fill(x_values, y_values, "r")

    arrow3 = Arrow(x1=550, x2=500, y=30)
    x_values, y_values = arrow3.get_coordinates()
    plt.fill(x_values, y_values, "r")

    plt.show()
