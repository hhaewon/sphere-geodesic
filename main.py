from dataclasses import dataclass, field
from functools import partial
from manim import *
import numpy as np


@dataclass()
class Point:
    """
    latitude: degrees
    longitude: degrees
    """

    latitude: int
    longitude: int
    radian_latitude: int = field(init=False)
    radian_logitude: int = field(init=False)

    def __post_init__(self):
        self.radian_latitude, self.radian_logitude = np.radians(
            [self.latitude, self.longitude]
        )

    @property
    def radian_angles(self):
        return (self.radian_latitude, self.radian_logitude)


# Convert spherical coordinates to Cartesian coordinates
def spherical_to_cartesian(point: Point) -> list[int]:
    lat, lon = point.radian_angles
    x = np.cos(lat) * np.cos(lon)
    y = np.cos(lat) * np.sin(lon)
    z = np.sin(lat)
    return [x, y, z]


def geodesic_path(point1: Point, point2: Point, t: int):
    lat1, lon1 = point1.radian_angles
    lat2, lon2 = point2.radian_angles
    # Interpolate between the two points on the sphere
    lat = np.arctan2(
        np.sin(lat1) * (1 - t) + np.sin(lat2) * t,
        np.sqrt(
            (np.cos(lat1) * np.cos(lon1) * (1 - t) + np.cos(lat2) * np.cos(lon2) * t)
            ** 2
            + (np.cos(lat1) * np.sin(lon1) * (1 - t) + np.cos(lat2) * np.sin(lon2) * t)
            ** 2
        ),
    )
    lon = np.arctan2(
        np.cos(lat1) * np.sin(lon1) * (1 - t) + np.cos(lat2) * np.sin(lon2) * t,
        np.cos(lat1) * np.cos(lon1) * (1 - t) + np.cos(lat2) * np.cos(lon2) * t,
    )
    return spherical_to_cartesian(lat, lon)


class SphereWithGeodesicScene(ThreeDScene):
    def construct(self):
        # Set up the 3D axes
        axes = ThreeDAxes()

        # Create a sphere
        sphere = Sphere(radius=1, color=BLUE, resolution=(50, 50))

        # Latitude and Longitude for the two points
        point1 = Point(45, 90)
        point2 = Point(-30, 60)

        point1 = spherical_to_cartesian(point1)
        point2 = spherical_to_cartesian(point2)

        # Create dots at the given latitude and longitude
        dot1 = Dot3D(point=point1, color=RED)
        dot2 = Dot3D(point=point2, color=GREEN)

        # Create a simple straight line between the two points
        line = Line3D(start=point1, end=point2, color=YELLOW)

        # Create the geodesic (great circle) line between the two points
        geodesic = ParametricFunction(
            partial(geodesic_path, point1=point1, point2=point2),
            t_range=np.array([0, 1]),
            color=ORANGE,
        )

        ## Add the sphere, dots, and geodesic to the scene
        self.add(axes, sphere, dot1, dot2, geodesic)

        # Rotate the camera to give a better view of the sphere
        self.set_camera_orientation(phi=75 * DEGREES, theta=30 * DEGREES)

        # Animate the rotation of the sphere
        self.play(Rotate(sphere, angle=PI, axis=RIGHT))
        self.wait()
