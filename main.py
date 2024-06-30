# ruff: noqa: F403, F405

from functools import partial
from manim import *  # type: ignore
import numpy as np
from numpy.typing import NDArray
from manim.typing import Point3D


class SphericalPoint:
    def __init__(self, latitude: float, longitude: float):
        self.phi: float = np.radians(90 - latitude)
        self.theta: float = (
            np.radians(longitude) if latitude >= 0 else np.radians(180 - longitude)
        )


# Convert spherical coordinates to Cartesian coordinates
def spherical_to_cartesian(point: SphericalPoint) -> Point3D:
    phi, theta = point.phi, point.theta
    x: int = np.cos(theta) * np.sin(phi)
    y: int = np.sin(theta) * np.sin(phi)
    z: int = np.cos(phi)
    return (x, y, z)


def geodesic_path(
    t: float,
    v1: NDArray[np.float64],
    v2: NDArray[np.float64],
) -> Point3D:
    w: NDArray[np.float64] = v2 - (v1 @ v2) * v1
    u: NDArray[np.float64] = w / np.linalg.norm(w)

    v = np.cos(t) * v1 + np.sin(t) * u
    return tuple(v)


class SphereWithGeodesicScene(ThreeDScene):
    CONFIG = {
        "x_axis_label": "$x$",
        "y_axis_label": "$y$",
        "z_axis_label": "$z$",
    }

    def construct(self):
        # Set up the 3D axes
        axes = ThreeDAxes()
        axes.add(axes.get_axis_labels())

        # Create a sphere
        sphere = Sphere(radius=1, color=BLUE, resolution=(50, 50))

        # Latitude and Longitude for the two points
        point1 = SphericalPoint(45, 90)
        point2 = SphericalPoint(-30, 60)

        cartesian_point1 = spherical_to_cartesian(point1)
        cartesian_point2 = spherical_to_cartesian(point2)

        v1: NDArray[np.float64] = np.array(cartesian_point1)
        v2: NDArray[np.float64] = np.array(cartesian_point2)
        c = np.arccos(v1 @ v2)

        # Create dots at the given latitude and longitude
        dot1 = Dot3D(point=list(cartesian_point1), color=RED)
        dot2 = Dot3D(point=list(cartesian_point2), color=GREEN)

        # Create the geodesic (great circle) line between the two points
        geodesic = ParametricFunction(
            partial(geodesic_path, v1=v1, v2=v2),
            t_range=np.array([0, c]),
            color=ORANGE,
        )

        ## Add the sphere, dots, and geodesic to the scene
        self.add(axes, sphere, dot1, dot2, geodesic)

        # Rotate the camera to give a better view of the sphere
        self.set_camera_orientation(phi=75 * DEGREES, theta=30 * DEGREES)
        self.wait()
        # self.set_camera_orientation(
        #     phi=(point1.phi + point2.phi) / 2,
        #     theta=(point1.theta + point2.theta) / 2,
        # )
