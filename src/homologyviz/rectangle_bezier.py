"""Class for plotting rectangles with curve heights to represent sequence homologies.

License
-------
This file is part of HomologyViz
BSD 3-Clause License
Copyright (c) 2024, Ivan Munoz Gutierrez
"""

import bezier
import numpy as np
import plotly.graph_objects as go


class RectangleCurveHeight:
    """Make coordinates to represent a rectangle with curve side heights.

    Attributes
    ----------
    x_coordinates : list[float]
    y_coordintaes : list[float]
    proportions : list[float]
        The proportions list must have values between 0 and 1. The list must start at 0
        and end at 1. The values must grow from 0 to 1. The values determine the Bezier
        curve.
    num_points : int
        Number of points to plot the height curves.
    """

    def __init__(
        self,
        x_coordinates: list[float],
        y_coordinates: list[float],
        proportions: list[float] = [0, 0.1, 0.5, 0.9, 1],
        num_points: int = 100,
    ):
        self.x_coordinates = x_coordinates
        self.y_coordinates = y_coordinates
        self.proportions = proportions
        self.degree = len(proportions) - 1
        self.num_points = num_points

    def coordinates_rectangle_height_bezier(self) -> tuple[np.ndarray, np.ndarray]:
        """Get coordinates to plot a rectangle with curve side heights."""
        x_right, y_right = self.get_bezier_nodes_vertical(
            x1=self.x_coordinates[1],
            x2=self.x_coordinates[2],
            y1=self.y_coordinates[1],
            y2=self.y_coordinates[2],
            proportions=self.proportions,
        )

        x_left, y_left = self.get_bezier_nodes_vertical(
            x1=self.x_coordinates[3],
            x2=self.x_coordinates[0],
            y1=self.y_coordinates[3],
            y2=self.y_coordinates[0],
            proportions=self.proportions,
        )
        x_points = np.concatenate((x_left, x_right))
        y_points = np.concatenate((y_left, y_right))

        x_points = np.append(x_points, self.x_coordinates[3])
        y_points = np.append(y_points, self.y_coordinates[3])

        return (x_points, y_points)

    def get_bezier_curve(self, curve, num_points=100):
        s_vals = np.linspace(0.0, 1.0, num_points)
        curve_points = curve.evaluate_multi(s_vals)
        return curve_points[0, :], curve_points[1, :]

    def get_bezier_nodes_vertical(
        self,
        x1: float,
        x2: float,
        y1: float,
        y2: float,
        proportions: list[float] = [0, 0.1, 0.5, 0.9, 1],
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        The proportions list values must be in the 0 to 1 range. The list must start with
        0 and end with 1.
        """
        degree = len(proportions) - 1
        x_coordinates = self.x_points_bezier_vertical(x1, x2, proportions)
        y_coordinates = self.y_points_bezier_vertical(y1, y2, degree)
        nodes = list(map(list, zip(x_coordinates, y_coordinates)))
        nodes = np.asfortranarray(nodes).T
        curve = bezier.Curve(nodes, degree)
        curve_x, curve_y = self.get_bezier_curve(curve)
        return (curve_x, curve_y)

    def y_points_bezier_vertical(
        self, y1: float, y2: float, degree: int
    ) -> list[float]:
        """Make y values to create a bezier vertical line."""
        delta = y2 - y1
        proportion = delta / degree
        values = [y1]
        for i in range(degree - 1):
            y = values[i] + proportion
            values.append(y)
        values.append(y2)
        return values

    def x_points_bezier_vertical(
        self, x1: float, x2: float, proportions: list[float] = [0, 0.1, 0.5, 0.9, 1]
    ) -> list[float]:
        """Make x values to create a bezier vertical line following the list proportions.

        The proportions list values must be in the 0 to 1 range. The list must start with
        0 and end with 1.
        """
        delta = x2 - x1
        return [x1 + proportion * delta for proportion in proportions]


if __name__ == "__main__":
    proportions = [0, 0.1, 0.5, 0.9, 1]
    x_coordinates = [0, 1, 1.5, 0.5]
    y_coordinates = [10, 10, 0, 0]
    rectangle = RectangleCurveHeight(
        x_coordinates=x_coordinates,
        y_coordinates=y_coordinates,
        proportions=proportions,
    )
    # Get coordinates for rectange with a bezier height
    x_points, y_points = rectangle.coordinates_rectangle_height_bezier()
    # Create figure
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=x_points, y=y_points, mode="lines", fill="toself", line=dict(color="blue")
        )
    )

    # Update layout
    fig.update_layout(
        title="Rectangle with Bezier Curved Sides and Straight Top/Bottom",
        xaxis_title="X-axis",
        yaxis_title="Y-axis",
        showlegend=False,
    )

    # Show figure
    fig.show()
