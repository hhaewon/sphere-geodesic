# ruff: noqa: F403, F405

from functools import partial
from manim import *  # type: ignore
import numpy as np
from numpy.typing import NDArray
from manim.typing import Point3D


class SphericalPoint:
    def __init__(self, latitude: float, longitude: float):
        self.latitude = np.radians(latitude)
        self.longitude = np.radians(longitude)
        self.phi: float = np.radians(90 - latitude)
        self.theta: float = (
            np.radians(longitude) if latitude >= 0 else np.radians(180 - longitude)
        )


# Convert spherical coordinates to Cartesian coordinates
def convert_spherical_to_cartesian(point: SphericalPoint) -> Point3D:
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


def get_parametric_geodesic_function(
    v1: NDArray[np.float64],
    v2: NDArray[np.float64],
) -> ParametricFunction:
    c = np.arccos(v1 @ v2)
    geodesic = ParametricFunction(
        partial(geodesic_path, v1=v1, v2=v2),
        t_range=np.array([0, c]),
        color=ORANGE,
    )
    return geodesic


def convert_sphere_to_Mercator(point: SphericalPoint, R: float) -> tuple[float, float]:
    latitude, longitude = point.latitude, point.longitude
    x = R * longitude
    y = R * np.log(np.tan(np.pi / 4 + latitude / 2))

    return (x, y)


def convert_Mercator_to_sphere(point: tuple[float, float], R: float) -> SphericalPoint:
    x, y = point
    latitude = x / R
    longitude = 2 * np.arctan(np.exp(y / R)) - np.pi / 2
    return SphericalPoint(latitude=latitude, longitude=longitude)


class SphereWithGeodesicScene(ThreeDScene):
    CONFIG = {
        "x_axis_label": "$x$",
        "y_axis_label": "$y$",
        "z_axis_label": "$z$",
    }

    def construct(self):
        global_map = ImageMobject("world_map.jpg", z_index=0).scale(1.5)
        R = global_map.get_right()[0] / np.pi

        # Latitude and Longitude for the two points
        point1 = SphericalPoint(45, 90)
        point2 = SphericalPoint(-30, 60)

        mercator_point1 = convert_sphere_to_Mercator(point=point1, R=R)
        mercator_point2 = convert_sphere_to_Mercator(point=point2, R=R)

        origin = (0, 0)
        dot_mercator1 = Dot3D(point=list(mercator_point1) + [0], color=RED, z_index=1)
        dot_mercator2 = Dot3D(point=list(mercator_point2) + [0], color=GREEN, z_index=1)
        dot_origin = Dot3D(point=list(origin) + [0], color=BLACK, z_index=1)

        # Set up the 3D axess
        axes = ThreeDAxes()
        # Create a sphere
        sphere = Sphere(radius=1, color=BLUE, resolution=(50, 50))

        cartesian_point1 = convert_spherical_to_cartesian(point1)
        cartesian_point2 = convert_spherical_to_cartesian(point2)

        v1: NDArray[np.float64] = np.array(cartesian_point1)
        v2: NDArray[np.float64] = np.array(cartesian_point2)

        # Create dots at the given latitude and longitude
        dot1 = Dot3D(point=list(cartesian_point1), color=RED)
        dot2 = Dot3D(point=list(cartesian_point2), color=GREEN)

        # Create the geodesic (great circle) line between the two points
        geodesic = get_parametric_geodesic_function(v1=v1, v2=v2)

        # animation
        self.add(global_map)
        self.play(FadeIn(global_map, run_time=1))
        vgroup1 = VGroup(dot_mercator1, dot_mercator2, dot_origin)
        self.add(vgroup1)
        self.play(Create(vgroup1))
        self.wait(duration=2)
        self.play(FadeOut(global_map, vgroup1, run_time=1))

        axes.add(axes.get_axis_labels())

        ## Add the sphere, dots, and geodesic to the scene
        vgroup2 = VGroup(axes, sphere, dot1, dot2, geodesic)
        self.add(vgroup2)

        # Rotate the camera to give a better view of the sphere
        self.set_camera_orientation(phi=75 * DEGREES, theta=30 * DEGREES)
        self.play(FadeIn(vgroup2))
        self.wait()
        # self.set_camera_orientation(
        #     phi=(point1.phi + point2.phi) / 2,
        #     theta=(point1.theta + point2.theta) / 2,
        # )
