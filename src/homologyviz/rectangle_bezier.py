"""
Generate coordinates for drawing rectangles with curved sides using Bézier curves.

This module defines the `RectangleCurveHeight` class, which calculates the x and y
coordinates needed to plot homology regions as rectangles with smoothly curved vertical
sides. This feature is intended for enhancing the visual distinction of homologous
sequences in the HomologyViz app, offering an alternative to traditional straight-edge
renderings.

Dependencies
------------
- `bezier`: Used to compute Bézier curve coordinates.
- `numpy`: For efficient numerical operations and array handling.
- `plotly.graph_objects` (optional): Used for plotting, though not directly invoked.

Usage
-----
This module is not yet integrated into the HomologyViz GUI, but future versions may allow
users to toggle between straight and curved homology representations.

License
-------
This file is part of HomologyViz
BSD 3-Clause License
Copyright (c) 2024, Iván Muñoz Gutiérrez
"""

import bezier
import numpy as np


class RectangleCurveHeight:
    """
    Generate coordinates for a rectangle with curved vertical sides using Bézier curves.

    This class is designed to construct the shape of a homology region where the left and
    right vertical edges are represented as Bézier curves, giving the region a smoother
    and more dynamic appearance in graphical sequence alignment plots.

    Parameters
    ----------
    x_coordinates : list[float]
        A list of 4 x-coordinates defining the corners of the rectangular region, ordered
        clockwise or counter-clockwise.
    y_coordinates : list[float]
        A list of 4 y-coordinates corresponding to `x_coordinates`.
    proportions : list[float], optional
        Proportional values between 0 and 1 defining the control points for the Bézier
        curve. Defaults to [0, 0.1, 0.5, 0.9, 1]. The curve shape depends on these.
    num_points : int, optional
        Number of points used to render each Bézier curve. More points yield smoother
        curves. Default is 100.

    Attributes
    ----------
    x_coordinates : list[float]
        Stores the x-values of the rectangle corners.
    y_coordinates : list[float]
        Stores the y-values of the rectangle corners.
    proportions : list[float]
        Used to shape the Bézier curves on the sides of the rectangle.
    degree : int
        Degree of the Bézier curve, inferred from the number of proportions.
    num_points : int
        Number of points used to evaluate and render the Bézier curves.

    Notes
    -----
    The Bézier curve rendering is powered by the `bezier` Python library. Ensure it is
    installed in your environment (e.g., via `pip install bezier`).

    This class is intended for internal use by the HomologyViz plotting system.
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
        """
        Get coordinates to plot a polygon resembling a rectangle with curved vertical
        sides.

        This method constructs the full x and y coordinate arrays needed to draw a
        homology region shaped like a rectangle, but with both vertical edges replaced by
        Bézier curves. The top and bottom edges are straight.

        Returns
        -------
        tuple : [numpy.ndarray, numpy.ndarray]
            A tuple containing:
            - x_points: The x-coordinates of the polygon.
            - y_points: The y-coordinates of the polygon.

        Notes
        -----
        The resulting polygon starts at the top-left, curves down the left edge,
        then follows the bottom edge to the right, curves up the right edge,
        and finally closes the shape by returning to the start.
        """
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

    def get_bezier_curve(
        self,
        curve: bezier.Curve,
        num_points: int = 100,
    ) -> tuple[np.array, np.array]:
        """
        Evaluate a Bézier curve at evenly spaced intervals.

        Parameters
        ----------
        curve : bezier.Curve
            A Bézier curve object created from control points using the `bezier` library.
        num_points : int, optional
            Number of points to evaluate along the curve (default is 100).

        Returns
        -------
        tuple : [numpy.ndarray, numpy.ndarray]
            A tuple containing the x and y coordinates of the evaluated Bézier curve.
        """
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
        Generate the x and y coordinates of a vertical Bézier curve between two points.

        This function computes the Bézier curve using x-coordinates interpolated from
        `x1` to `x2` based on the given `proportions`, and y-coordinates spaced evenly
        from `y1` to `y2` for the curve degree determined by the proportions list.

        Parameters
        ----------
        x1 : float
            Starting x-coordinate of the curve.
        x2 : float
            Ending x-coordinate of the curve.
        y1 : float
            Starting y-coordinate of the curve.
        y2 : float
            Ending y-coordinate of the curve.
        proportions : list of float, optional
            List of float values between 0 and 1 representing how control points are
            spaced along the x-axis. Must start at 0 and end at 1.

        Returns
        -------
        tuple : [numpy.ndarray, numpy.ndarray]
            The x and y coordinates of the Bézier curve evaluated at evenly spaced
            intervals.
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
        """
        Generate y-coordinates for a vertical Bézier curve of given degree.

        This function computes evenly spaced y-values between `y1` and `y2` for use as
        control points in a vertical Bézier curve. The number of output points equals
        `degree + 1`.

        Parameters
        ----------
        y1 : float
            Starting y-coordinate of the curve.
        y2 : float
            Ending y-coordinate of the curve.
        degree : int
            Degree of the Bézier curve (determines the number of control points as
            `degree + 1`).

        Returns
        -------
        list : [float]
            A list of y-coordinates evenly spaced between `y1` and `y2`.
        """
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
        """
        Generate x-coordinates for control points of a vertical Bézier curve.

        This function calculates a list of x-values spaced according to the specified
        `proportions` between `x1` and `x2`. These values are used to shape the curve
        horizontally, while the corresponding y-values are distributed vertically.

        Parameters
        ----------
        x1 : float
            Starting x-coordinate of the curve.
        x2 : float
            Ending x-coordinate of the curve.
        proportions : list of float, default=[0, 0.1, 0.5, 0.9, 1]
            List of normalized positions (between 0 and 1) to interpolate between `x1`
            and `x2`. Must start with 0 and end with 1. These determine the curvature
            profile.

        Returns
        -------
        list : [float]
            A list of x-coordinates for Bézier control points, matching the provided
            proportions.
        """
        delta = x2 - x1
        return [x1 + proportion * delta for proportion in proportions]
