"""Class for plotting arrows to represent genes.

Notes
-----
- This file is part of HomologyViz
- BSD 3-Clause License
- Copyright (c) 2024, Iván Muñoz Gutiérrez
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
        x1: int | float,
        x2: int | float,
        y: int | float,
        ratio_tail_head_width: float = 0.5,
        head_width: int | float = 2,
        head_height: int | float = 200,
    ):
        # Validate arrow dimensions.
        if not 0 < ratio_tail_head_width <= 1:
            raise ValueError(
                "ratio_tail_head_width must be greater than 0 and less or equal to 1."
            )

        if head_width <= 0:
            raise ValueError("head_width must be greater than 0.")

        if head_height <= 0:
            raise ValueError("head_height must be greater than 0.")

        if x1 == x2:
            raise ValueError("x1 and x2 must be different.")

        # Initialize arrow attributes.
        self.x1 = x1
        self.x2 = x2
        self.y = y
        self.ratio_tail_head_width = ratio_tail_head_width
        self.head_width = head_width
        self.head_height = head_height

        # Calculate arrow geometry
        self.tail_width = head_width * ratio_tail_head_width

    def __repr__(self) -> str:
        """Return a developer-friendly representation of the arrow."""
        return (
            f"Arrow(x1={self.x1!r}, x2={self.x2!r}, "
            f"y={self.y!r}, "
            f"ratio_tail_head_width={self.ratio_tail_head_width!r}, "
            f"head_width={self.head_width!r}, "
            f"head_height={self.head_height!r})"
        )

    def get_coordinates(self) -> tuple[ndarray, ndarray]:
        """
        Generate polygon coordinates for a forward or reverse arrow.

        Returns
        -------
        tuple of np.ndarray
            Arrays containing the x and y coordinates of the arrow polygon.
        """
        direction = 1 if self.x2 > self.x1 else -1

        length = abs(self.x2 - self.x1)

        # Short arrows consist only of a triangular head.
        if length <= self.head_height:
            x_values = np.array(
                [
                    self.x1,
                    self.x1,
                    self.x2,
                    self.x1,
                    self.x1,
                ]
            )

            y_values = np.array(
                [
                    self.y,
                    self.y + direction * self.head_width / 2,
                    self.y,
                    self.y - direction * self.head_width / 2,
                    self.y,
                ]
            )

            return x_values, y_values

        # Calculate the beginning of the arrowhead.
        head_start = self.x2 - direction * self.head_height

        # Generate x-coordinates.
        x_values = np.array(
            [
                self.x1,
                self.x1,
                head_start,
                head_start,
                self.x2,
                head_start,
                head_start,
                self.x1,
                self.x1,
            ]
        )

        # Generate y-coordinates.
        y_values = np.array(
            [
                self.y,
                self.y + direction * self.tail_width / 2,
                self.y + direction * self.tail_width / 2,
                self.y + direction * self.head_width / 2,
                self.y,
                self.y - direction * self.head_width / 2,
                self.y - direction * self.tail_width / 2,
                self.y - direction * self.tail_width / 2,
                self.y,
            ]
        )

        return x_values, y_values


if __name__ == "__main__":
    ####################
    # Example of usage #
    ####################

    # Make arrow
    arrow1 = Arrow(x1=750, x2=2000, y=10, ratio_tail_head_width=1.0)
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
